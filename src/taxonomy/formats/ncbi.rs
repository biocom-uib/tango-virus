use std::{
    collections::{hash_map::Entry, HashMap, HashSet},
    fmt,
    fs::File,
    io::{self, BufRead, BufReader},
    num::ParseIntError,
    path::Path,
};

use crate::taxonomy::NodeId;
use newick_rs::SimpleTree;
use serde::{Deserialize, Serialize, ser::SerializeMap, Serializer, Deserializer};
use string_interner::{
    backend::{Backend, StringBackend},
    StringInterner, Symbol,
};
use thiserror::Error;

pub type TaxId = NodeId;
pub type RankSymbol = <StringBackend as Backend>::Symbol;

#[derive(Deserialize, Serialize)]
pub struct NcbiTaxonomy<Names> {
    pub root: TaxId,

    pub parent_ids: HashMap<TaxId, TaxId>,
    children_lookup: HashMap<TaxId, Vec<TaxId>>,

    #[serde(deserialize_with = "deserialize_ranks", serialize_with = "serialize_ranks")]
    pub ranks: HashMap<TaxId, RankSymbol>,
    pub ranks_interner: StringInterner,

    // old to new
    pub merged_taxids: Option<HashMap<TaxId, TaxId>>,

    pub names: Names,
}

pub fn serialize_ranks<S: Serializer>(ranks: &HashMap<TaxId, RankSymbol>, serializer: S) -> Result<S::Ok, S::Error> {
    let mut ranks_s = serializer.serialize_map(Some(ranks.len()))?;

    for (node, rank_sym) in ranks.iter() {
        ranks_s.serialize_entry(node, &rank_sym.to_usize())?;
    }

    ranks_s.end()
}

pub fn deserialize_ranks<'de, D: Deserializer<'de>>(
    deserializer: D,
) -> Result<HashMap<TaxId, RankSymbol>, D::Error> {
    use serde::de::{MapAccess, Visitor};

    struct RanksVisitor {}

    impl<'de> Visitor<'de> for RanksVisitor {
        type Value = HashMap<TaxId, RankSymbol>;

        fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
            formatter.write_str("ranks dictionary")
        }

        fn visit_map<M: MapAccess<'de>>(self, mut access: M) -> Result<Self::Value, M::Error> {
            let mut res = HashMap::with_capacity(access.size_hint().unwrap_or(0));

            while let Some((node, rank_sym_usize)) = access.next_entry()? {
                let rank_sym = RankSymbol::try_from_usize(rank_sym_usize).ok_or_else(|| {
                    <M::Error as serde::de::Error>::invalid_value(
                        serde::de::Unexpected::Unsigned(rank_sym_usize as u64),
                        &"a valid StringInterner symbol (as usize)",
                    )
                })?;

                res.insert(node, rank_sym);
            }

            Ok(res)
        }
    }

    deserializer.deserialize_map(RanksVisitor {})
}

#[derive(Debug, Error)]
pub enum DmpError {
    #[error("Root not found")]
    MissingRoot,

    #[error("Error parsing taxid: {}", .0)]
    TaxIdParseError(#[from] ParseIntError),

    #[error("Missing .dmp field: {}", .0)]
    MissingField(&'static str),

    #[error("error reading .dmp file")]
    ReadError(#[from] io::Error),

    #[error("multiple parents found for {}: {} and {}", .0, .1, .2)]
    MultipleParents(TaxId, TaxId, TaxId),
}

macro_rules! parse_fields {
    ($line:expr => { $(let $fields:pat = next as $field_names:expr;)+ }; $($body:tt)+) => {{
        let mut line_iter = $line.split("\t|\t");

        $(
            let $fields = line_iter.next().ok_or(DmpError::MissingField($field_names))?;
        )+

        $($body)+
    }};
}

impl<Names> NcbiTaxonomy<Names> {
    pub fn with_merged_taxids(self, merged_dmp: &Path) -> Result<Self, DmpError> {
        let mut merged = HashMap::new();

        for line in BufReader::new(File::open(merged_dmp)?).lines() {
            let line = line?;
            let line = line.strip_suffix("\t|").unwrap_or(&line);

            parse_fields! {
                &line => {
                    let old_taxid  = next as "old_tax_id";
                    let new_taxid  = next as "new_tax_id";
                };

                let old_taxid = NodeId(old_taxid.parse()?);
                let new_taxid = NodeId(new_taxid.parse()?);

                merged.insert(old_taxid, new_taxid);
            }
        }

        Ok(Self {
            merged_taxids: Some(merged),
            ..self
        })
    }

    pub fn with_names<NewNames>(self, new_names: NewNames) -> NcbiTaxonomy<NewNames> {
        NcbiTaxonomy {
            root: self.root,
            parent_ids: self.parent_ids,
            children_lookup: self.children_lookup,
            ranks: self.ranks,
            ranks_interner: self.ranks_interner,
            merged_taxids: self.merged_taxids,
            names: new_names,
        }
    }

    pub fn node_count(&self) -> usize {
        1 /* root */ + self.parent_ids.len()
    }

    pub fn is_leaf(&self, node: TaxId) -> bool {
        self.children(node).is_empty()
    }

    pub fn children(&self, node: TaxId) -> &[TaxId] {
        self.children_lookup
            .get(&node)
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    pub fn parents(&self, start: TaxId) -> impl Iterator<Item = TaxId> + '_ {
        super::super::util::ParentsIter::new_with(start, |node| {
            if node == self.root {
                None
            } else {
                self.parent_ids.get(&node).copied()
            }
        })
    }

    pub fn preorder_descendants(&self, start: TaxId) -> impl Iterator<Item = TaxId> + '_ {
        super::super::util::PreOrderIter::new_with(start, |node| {
            self.children(node).iter().copied()
        })
    }

    pub fn postorder_descendants(&self, start: TaxId) -> impl Iterator<Item = TaxId> + '_ {
        super::super::util::PostOrderIter::new_with(start, |node| {
            self.children(node).iter().copied()
        })
    }

    pub fn contract<S: AsRef<str>>(&mut self, ranks: &[S]) {
        let mut ranks_syms = HashSet::new();

        for rank in ranks {
            if let Some(rank_sym) = self.ranks_interner.get(rank.as_ref()) {
                ranks_syms.insert(rank_sym);
            } else {
                eprintln!(
                    "NcbiTaxonomy::contract: warning: rank {:?} does not appear in the taxonomy",
                    rank.as_ref()
                );
            }
        }

        let mut parent_ids = HashMap::new();
        let mut children_lookup = HashMap::new();

        let mut add_node = |parent, child| {
            parent_ids.insert(child, parent);
            children_lookup
                .entry(parent)
                .or_insert_with(|| Vec::new())
                .push(child);
        };

        let mut stack = vec![(self.root, self.root, self.children(self.root).iter())];

        while let Some(&mut (cur_parent, cur, ref mut cur_children)) = stack.last_mut() {
            if let Some(&cur_child) = cur_children.next() {
                let cur_grandchildren = self.children(cur);

                // leaf
                if cur_grandchildren.is_empty() {
                    add_node(cur_parent, cur_child);

                // internal
                } else {
                    match self.ranks.get(&cur_child) {
                        Some(rank) if ranks_syms.contains(rank) => {
                            add_node(cur_parent, cur_child);

                            stack.push((cur_child, cur_child, cur_grandchildren.iter()));
                        }
                        _ => {
                            // cur is not a valid parent
                            stack.push((cur_parent, cur_child, cur_grandchildren.iter()));
                        }
                    }
                }
            } else {
                stack.pop();
            }
        }

        self.parent_ids = parent_ids;
        self.children_lookup = children_lookup;
    }
}

impl NcbiTaxonomy<NoNames> {
    pub fn load_nodes(nodes_dmp: &Path) -> Result<Self, DmpError> {
        let mut root = None;

        let mut parent_ids = HashMap::new();
        let mut children_lookup = HashMap::new();

        let mut ranks = HashMap::new();
        let mut ranks_interner = StringInterner::new();

        for line in BufReader::new(File::open(nodes_dmp)?).lines() {
            let line = line?;

            parse_fields! {
                &line => {
                    let taxid  = next as "tax id";
                    let ptaxid = next as "parent tax_id";
                    let rank   = next as "rank";
                };

                let taxid = NodeId(taxid.parse()?);
                let ptaxid = NodeId(ptaxid.parse()?);

                ranks.insert(taxid, ranks_interner.get_or_intern(rank));

                if taxid == ptaxid {
                    root = Some(taxid);

                } else {
                    match parent_ids.entry(taxid) {
                        Entry::Vacant(vac) => {
                            vac.insert(ptaxid);
                        },
                        Entry::Occupied(occ) => {
                            return Err(DmpError::MultipleParents(taxid, ptaxid, *occ.get()));
                        },
                    }

                    children_lookup.entry(ptaxid).or_insert_with(|| Vec::new()).push(taxid);
                }
            }
        }

        Ok(NcbiTaxonomy {
            root: root.ok_or(DmpError::MissingRoot)?,
            parent_ids,
            children_lookup,
            ranks,
            ranks_interner,
            merged_taxids: None,
            names: NoNames::new(),
        })
    }
}

pub type NameClassSymbol = <StringBackend as Backend>::Symbol;

pub struct AllNames {
    pub name_classes: StringInterner,
    pub unique_names: HashMap<TaxId, HashMap<NameClassSymbol, String>>,
    pub names_lookup: HashMap<String, Vec<(NameClassSymbol, TaxId)>>,
}

impl AllNames {
    pub fn load_names_dmp<P: AsRef<Path>>(names_dmp: P) -> Result<Self, DmpError> {
        let mut this = AllNames {
            name_classes: StringInterner::new(),
            unique_names: HashMap::new(),
            names_lookup: HashMap::new(),
        };

        this.fill_from(names_dmp)?;
        Ok(this)
    }

    pub fn only_of_class(self, name_class: &str) -> Option<SingleClassNames> {
        let single_class = self.name_classes.get(name_class)?;

        let unique_names = self
            .unique_names
            .into_iter()
            .filter_map(|(taxid, mut names)| Some((taxid, names.remove(&single_class)?)))
            .collect();

        let names_lookup = self
            .names_lookup
            .into_iter()
            .filter_map(|(name, taxids)| {
                let new_taxids: Vec<TaxId> = taxids
                    .into_iter()
                    .filter(|(class, _)| *class == single_class)
                    .map(|(_, taxid)| taxid)
                    .collect();

                if new_taxids.is_empty() {
                    None
                } else {
                    Some((name, new_taxids))
                }
            })
            .collect();

        Some(SingleClassNames {
            name_class: name_class.to_owned(),
            unique_names,
            names_lookup,
        })
    }
}

#[derive(Deserialize, Serialize)]
pub struct SingleClassNames {
    pub name_class: String,
    pub unique_names: HashMap<TaxId, String>,
    pub names_lookup: HashMap<String, Vec<TaxId>>,
}

impl SingleClassNames {
    pub fn load_names_dmp<P: AsRef<Path>>(
        name_class: String,
        names_dmp: P,
    ) -> Result<Self, DmpError> {
        let mut this = SingleClassNames {
            name_class,
            unique_names: HashMap::new(),
            names_lookup: HashMap::new(),
        };

        this.fill_from(names_dmp)?;
        Ok(this)
    }
}

#[derive(Copy, Clone, Deserialize, Serialize)]
pub struct NoNames {}

impl NoNames {
    pub fn new() -> Self {
        NoNames {}
    }
}

pub trait NamesAssoc {
    fn insert(&mut self, taxid: TaxId, name: &str, unique_name: &str, unique_name: &str);

    type NameLookup;
    fn lookup_names<'a>(&'a self, taxid: TaxId) -> Option<&'a Self::NameLookup>;

    type NamesLookupIter<'a>: Iterator<Item = &'a str>;
    fn iter_lookup_names<'a>(lookup: &'a Self::NameLookup) -> Self::NamesLookupIter<'a>;

    type TaxIdsLookup;
    fn lookup_taxids<'a>(&'a self, name: &str) -> Option<&'a Self::TaxIdsLookup>;

    type TaxIdsLookupIter<'a>: Iterator<Item = TaxId>;
    fn iter_lookup_taxids<'a>(lookup: &'a Self::TaxIdsLookup) -> Self::TaxIdsLookupIter<'a>;

    fn fill_from<P: AsRef<Path>>(&mut self, names_dmp: P) -> Result<(), DmpError> {
        for line in BufReader::new(File::open(names_dmp)?).lines() {
            let line = line?;

            parse_fields! {
                &line => {
                    let taxid = next as "tax_id";
                    let name = next as "name_txt";
                    let unique = next as "unique name";
                    let name_class = next as "name class";
                };

                let taxid = NodeId(taxid.parse()?);

                self.insert(taxid, name_class, name, unique);
            }
        }

        Ok(())
    }
}

impl NamesAssoc for AllNames {
    fn insert(&mut self, taxid: TaxId, name: &str, unique_name: &str, name_class: &str) {
        let name_class = self.name_classes.get_or_intern(name_class);

        self.unique_names
            .entry(taxid)
            .or_default()
            .entry(name_class)
            .or_insert_with(|| unique_name.to_owned());

        self.names_lookup
            .entry(name.to_owned())
            .or_default()
            .push((name_class, taxid));
    }

    type NameLookup = HashMap<NameClassSymbol, String>;
    fn lookup_names<'a>(&'a self, taxid: TaxId) -> Option<&'a Self::NameLookup> {
        self.unique_names.get(&taxid)
    }

    type NamesLookupIter<'a> = impl Iterator<Item = &'a str>;
    fn iter_lookup_names<'a>(lookup: &'a Self::NameLookup) -> Self::NamesLookupIter<'a> {
        lookup.values().map(|name| name.as_str())
    }

    type TaxIdsLookup = Vec<(NameClassSymbol, TaxId)>;
    fn lookup_taxids<'a>(&'a self, name: &str) -> Option<&'a Self::TaxIdsLookup> {
        self.names_lookup.get(name)
    }

    type TaxIdsLookupIter<'a> = impl Iterator<Item = TaxId>;
    fn iter_lookup_taxids<'a>(lookup: &'a Self::TaxIdsLookup) -> Self::TaxIdsLookupIter<'a> {
        lookup.iter().map(|(_, taxid)| *taxid)
    }
}

impl NamesAssoc for SingleClassNames {
    fn insert(&mut self, taxid: TaxId, name: &str, unique_name: &str, name_class: &str) {
        if name_class == self.name_class {
            self.unique_names
                .entry(taxid)
                .or_insert_with(|| unique_name.to_owned());

            self.names_lookup
                .entry(name.to_owned())
                .or_default()
                .push(taxid);
        }
    }

    type NameLookup = String;
    fn lookup_names<'a>(&'a self, taxid: TaxId) -> Option<&'a Self::NameLookup> {
        self.unique_names.get(&taxid)
    }

    type NamesLookupIter<'a> = std::iter::Once<&'a str>;
    fn iter_lookup_names<'a>(lookup: &'a Self::NameLookup) -> Self::NamesLookupIter<'a> {
        std::iter::once(&*lookup)
    }

    type TaxIdsLookup = Vec<TaxId>;
    fn lookup_taxids<'a>(&'a self, name: &str) -> Option<&'a Self::TaxIdsLookup> {
        self.names_lookup.get(name)
    }

    type TaxIdsLookupIter<'a> = std::iter::Cloned<std::slice::Iter<'a, TaxId>>;
    fn iter_lookup_taxids<'a>(lookup: &'a Self::TaxIdsLookup) -> Self::TaxIdsLookupIter<'a> {
        lookup.iter().cloned()
    }
}

impl NamesAssoc for NoNames {
    fn insert(&mut self, _taxid: TaxId, _name: &str, _unique_name: &str, _name_class: &str) {}

    type NameLookup = !;
    fn lookup_names<'a>(&'a self, _taxid: TaxId) -> Option<&'a Self::NameLookup> {
        None
    }

    type NamesLookupIter<'a> = std::iter::Empty<&'a str>;
    fn iter_lookup_names<'a>(_lookup: &'a Self::NameLookup) -> Self::NamesLookupIter<'a> {
        std::iter::empty()
    }

    type TaxIdsLookup = !;
    fn lookup_taxids<'a>(&'a self, _name: &str) -> Option<&'a Self::TaxIdsLookup> {
        None
    }

    type TaxIdsLookupIter<'a> = std::iter::Empty<TaxId>;
    fn iter_lookup_taxids<'a>(_lookup: &'a Self::TaxIdsLookup) -> Self::TaxIdsLookupIter<'a> {
        std::iter::empty()
    }

    fn fill_from<P: AsRef<Path>>(&mut self, _names_dmp: P) -> Result<(), DmpError> {
        Ok(())
    }
}

fn make_simple_tree(
    node_id: TaxId,
    children_lookup: &mut HashMap<TaxId, Vec<TaxId>>,
    names: &mut HashMap<TaxId, String>,
) -> SimpleTree {
    let children = children_lookup
        .get_mut(&node_id)
        .map(std::mem::take)
        .unwrap_or_default()
        .into_iter()
        .map(|child_id| make_simple_tree(child_id, children_lookup, names))
        .collect();

    SimpleTree {
        name: names.remove(&node_id).unwrap_or_default(),
        children,
        length: None,
    }
}

impl From<NcbiTaxonomy<SingleClassNames>> for SimpleTree {
    fn from(taxonomy: NcbiTaxonomy<SingleClassNames>) -> Self {
        let root = taxonomy.root;
        let mut children_lookup = taxonomy.children_lookup;
        let mut names = taxonomy.names.unique_names;

        make_simple_tree(root, &mut children_lookup, &mut names)
    }
}
