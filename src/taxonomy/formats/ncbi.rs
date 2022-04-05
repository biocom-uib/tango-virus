use std::{
    collections::HashMap,
    fs::File,
    io::{self, BufRead, BufReader},
    num::ParseIntError,
    path::Path, slice,
};

use crate::taxonomy::{NodeId, Taxonomy, TaxonomyMut, generic::{GenericTaxonomy, TaxonomyBuildError}, TopologyReplacer, LabelledTaxonomy};
use itertools::Either;
use newick_rs::SimpleTree;
use serde::{Deserialize, Serialize};
use string_interner::{
    backend::{Backend, StringBackend},
    StringInterner, symbol::SymbolU32, Symbol,
};
use thiserror::Error;

pub type TaxId = NodeId;

#[derive(Deserialize, Serialize)]
pub struct NcbiTaxonomy<Names> {
    pub(crate) tree: GenericTaxonomy,

    // old to new
    pub(crate) merged_taxids: Option<HashMap<TaxId, TaxId>>,

    pub(crate) names: Names,
}

#[derive(Debug, Error)]
pub enum DmpError {
    #[error("Tree consistency error")]
    TreeBuildError(#[from] TaxonomyBuildError),

    #[error("Error parsing taxid: {}", .0)]
    TaxIdParseError(#[from] ParseIntError),

    #[error("Missing .dmp field: {}", .0)]
    MissingField(&'static str),

    #[error("error reading .dmp file")]
    ReadError(#[from] io::Error),
}

macro_rules! parse_fields {
    ($line:expr => { $(let $fields:pat = next as $field_names:expr;)+ }; $($body:tt)+) => {{
        let mut line_iter = $line.strip_suffix("\t|").unwrap_or($line).split("\t|\t");

        $(
            let $fields = line_iter.next().ok_or(DmpError::MissingField($field_names))?;
        )+

        $($body)+
    }};
}

impl<Names> NcbiTaxonomy<Names> {
    pub fn node_count(&self) -> usize {
        1 /* root */ + self.tree.parent_ids.len()
    }

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

    pub fn fixup_node(&self, node: NodeId) -> Option<NodeId> {
        if self.tree.has_node(node) {
            Some(node)
        } else if let Some(merged) = &self.merged_taxids {
            merged.get(&node).copied()
        } else {
            None
        }
    }

    pub fn with_names<NewNames>(self, new_names: NewNames) -> NcbiTaxonomy<NewNames> {
       NcbiTaxonomy {
          tree: self.tree,
          merged_taxids: self.merged_taxids,
          names: new_names,
       }
    }

    pub fn map_names<F, NewNames>(self, f: F) -> NcbiTaxonomy<NewNames>
    where
        F: FnOnce(Names) -> NewNames
    {
       NcbiTaxonomy {
          tree: self.tree,
          merged_taxids: self.merged_taxids,
          names: f(self.names),
       }
    }
}

impl NcbiTaxonomy<NoNames> {
    pub fn load_nodes<P: AsRef<Path>>(nodes_dmp: P) -> Result<Self, DmpError> {
        let mut builder = GenericTaxonomy::builder();

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

                builder.set_rank(taxid, rank);

                if taxid == ptaxid {
                    builder.set_root(taxid)?;
                } else {
                    builder.insert_edge(ptaxid, taxid)?;
                }
            }
        }

        let tree = builder.build()?;

        Ok(NcbiTaxonomy {
            tree,
            merged_taxids: None,
            names: NoNames::new(),
        })
    }
}

impl<Names: 'static> Taxonomy for NcbiTaxonomy<Names> {
    fn get_root(&self) -> TaxId {
        self.tree.get_root()
    }

    fn is_leaf(&self, node: TaxId) -> bool {
        self.tree.is_leaf(node)
    }

    fn has_uniform_depths(&self) -> Option<usize> {
        self.tree.has_uniform_depths()
    }

    type Children<'a> = <GenericTaxonomy as Taxonomy>::Children<'a>;

    fn iter_children<'a>(&'a self, node: TaxId) -> Self::Children<'a> {
        self.tree.iter_children(node)
    }

    fn find_parent(&self, node: TaxId) -> Option<TaxId> {
        self.tree.find_parent(node)
    }

    type RankSym = <StringBackend as Backend>::Symbol;

    fn rank_sym_str(&self, rank_sym: Self::RankSym) -> Option<&str> {
        self.tree.rank_sym_str(rank_sym)
    }

    fn lookup_rank_sym(&self, rank: &str) -> Option<Self::RankSym> {
        self.tree.lookup_rank_sym(rank)
    }

    fn find_rank(&self, node: TaxId) -> Option<Self::RankSym> {
        self.tree.find_rank(node)
    }

    fn get_rank(&self, node: TaxId) -> Self::RankSym {
        self.tree.get_rank(node)
    }

    type NodeRanks<'a> = <GenericTaxonomy as Taxonomy>::NodeRanks<'a>;

    fn node_ranks<'a>(&'a self) -> Self::NodeRanks<'a> {
        self.tree.node_ranks()
    }
}

impl<Names: NamesAssoc + Send + 'static> LabelledTaxonomy for NcbiTaxonomy<Names> {
    type Labels<'a> = Either<
        std::iter::Empty<&'a str>,
        Names::NamesLookupIter<'a>
    >;

    fn labels_of<'a>(&'a self, node: NodeId) -> Self::Labels<'a> {
        if let Some(result) = self.names.lookup_names(node) {
            Either::Right(Names::iter_lookup_names(result))
        } else {
            Either::Left(std::iter::empty())
        }
    }

    type NodesWithLabel<'a> = Either<
        std::iter::Empty<NodeId>,
        Names::TaxIdsLookupIter<'a>
    >;

    fn nodes_with_label<'a>(&'a self, label: &'a str) -> Self::NodesWithLabel<'a> {
        if let Some(result) = self.names.lookup_taxids(label) {
            Either::Right(Names::iter_lookup_taxids(result))
        } else {
            Either::Left(std::iter::empty())
        }
    }
}

impl<Names: NamesAssoc + Send + 'static> TaxonomyMut for NcbiTaxonomy<Names> {
    type UnderlyingTopology = <GenericTaxonomy as TaxonomyMut>::UnderlyingTopology;

    fn replace_topology_with<'a, Replacer>(&'a mut self, replacer: Replacer)
    where
        Replacer: TopologyReplacer<Self::UnderlyingTopology>
    {
        self.tree.replace_topology_with(replacer);

        let known_taxid = |taxid| self.tree.has_node(taxid);

        rayon::join(
            || self.names.forget_taxids(&known_taxid),
            || {
                if let Some(merged) = &mut self.merged_taxids {
                    merged.retain(|_old_taxid, new_taxid| known_taxid(*new_taxid));
                }
            },
        );
    }
}


#[derive(Copy, Clone, Debug, Hash, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct NameClassSymbol(u32);

impl From<SymbolU32> for NameClassSymbol {
    fn from(sym: SymbolU32) -> Self {
       NameClassSymbol(sym.to_usize() as u32)
    }
}

impl From<NameClassSymbol> for SymbolU32 {
   fn from(sym: NameClassSymbol) -> Self {
      <SymbolU32 as Symbol>::try_from_usize(sym.0 as usize).unwrap()
   }
}

#[derive(Default, Serialize, Deserialize)]
pub struct AllNames {
    pub name_classes: StringInterner,
    pub taxid_to_names: HashMap<TaxId, Vec<(NameClassSymbol, String)>>,
    pub name_to_taxids: HashMap<String, Vec<(NameClassSymbol, TaxId)>>,
}

impl AllNames {
    pub fn load_filtered_names_dmp<S, P>(classes: &[S], names_dmp: P) -> Result<Self, DmpError>
    where
        S: AsRef<str>,
        P: AsRef<Path>,
    {
        let mut this = AllNames {
            name_classes: StringInterner::new(),
            taxid_to_names: HashMap::new(),
            name_to_taxids: HashMap::new(),
        };

        this.fill_from(names_dmp, |class| classes.iter().any(|c| c.as_ref() == class))?;

        let order: HashMap<NameClassSymbol, usize>  = classes.into_iter()
            .enumerate()
            .map(|(i, s)| (this.name_classes.get_or_intern(s.as_ref()).into(), i))
            .collect();

        for names in this.taxid_to_names.values_mut() {
            names.sort_by_key(|(class, _)| order.get(class).copied().unwrap_or(classes.len()));
        }

        for taxids in this.name_to_taxids.values_mut() {
            taxids.sort_by_key(|(class, _)| order.get(class).copied().unwrap_or(classes.len()));
        }

        Ok(this)
    }

    pub fn only_of_class(self, name_class: &str) -> Result<SingleClassNames, Self> {
        let single_class = if let Some(sym) = self.name_classes.get(name_class) {
           sym.into()
        } else {
           return Err(self)
        };

        let taxid_to_name = self
            .taxid_to_names
            .into_iter()
            .filter_map(|(taxid, names)| {
                let single_name = names.into_iter()
                    .filter(|(class, _)| *class == single_class)
                    .map(|(_, name)| name)
                    .next()?;

                Some((taxid, single_name))
            })
            .collect();

        let name_to_taxids = self
            .name_to_taxids
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

        Ok(SingleClassNames {
            name_class: name_class.to_owned(),
            taxid_to_name,
            name_to_taxids,
        })
    }
}

#[derive(Debug, Deserialize, Serialize)]
pub struct SingleClassNames {
    pub name_class: String,
    pub taxid_to_name: HashMap<TaxId, String>,
    pub name_to_taxids: HashMap<String, Vec<TaxId>>,
}

impl SingleClassNames {
    pub fn load_names_dmp<P: AsRef<Path>>(
        name_class: String,
        names_dmp: P,
    ) -> Result<Self, DmpError> {
        let mut this = SingleClassNames {
            name_class,
            taxid_to_name: HashMap::new(),
            name_to_taxids: HashMap::new(),
        };

        this.fill_from(names_dmp, |_| true)?;
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
    fn insert(&mut self, taxid: TaxId, name: &str, unique_name: &str, name_class: &str);

    fn forget_taxids<Keep>(&mut self, keep: Keep) where Keep: FnMut(TaxId) -> bool;

    type NameLookup;
    fn lookup_names<'a>(&'a self, taxid: TaxId) -> Option<&'a Self::NameLookup>;

    type NamesLookupIter<'a>: Iterator<Item = &'a str> + 'a;
    fn iter_lookup_names<'a>(lookup: &'a Self::NameLookup) -> Self::NamesLookupIter<'a>;

    type TaxIdsLookup;
    fn lookup_taxids<'a>(&'a self, name: &str) -> Option<&'a Self::TaxIdsLookup>;

    type TaxIdsLookupIter<'a>: Iterator<Item = TaxId> + 'a;
    fn iter_lookup_taxids<'a>(lookup: &'a Self::TaxIdsLookup) -> Self::TaxIdsLookupIter<'a>;

    fn fill_from<P, F>(&mut self, names_dmp: P, class_filter: F) -> Result<(), DmpError>
    where
        P: AsRef<Path>,
        F: Fn(&str) -> bool,
    {
        for line in BufReader::new(File::open(names_dmp)?).lines() {
            let line = line?;

            parse_fields! {
                &line => {
                    let taxid = next as "tax_id";
                    let name = next as "name_txt";
                    let unique_name = next as "unique name";
                    let name_class = next as "name class";
                };

                let taxid = NodeId(taxid.parse()?);

                if class_filter(name_class) {
                    self.insert(taxid, name, unique_name, name_class);
                }
            }
        }

        Ok(())
    }
}

impl NamesAssoc for AllNames {
    fn insert(&mut self, taxid: TaxId, name: &str, unique_name: &str, name_class: &str) {
        let chosen_name = if !name.is_empty() {
            name
        } else if !unique_name.is_empty() {
            unique_name
        } else {
            return
        };

        let name_class = self.name_classes.get_or_intern(name_class).into();

        self.taxid_to_names
            .entry(taxid)
            .or_default()
            .push((name_class, chosen_name.to_owned()));

        self.name_to_taxids
            .entry(chosen_name.to_owned())
            .or_default()
            .push((name_class, taxid));
    }

    fn forget_taxids<Keep: FnMut(TaxId) -> bool>(&mut self, mut keep: Keep) {
        self.taxid_to_names.retain(|taxid, _| keep(*taxid));

        self.name_to_taxids.retain(|_, taxids| {
            taxids.retain(|(_, taxid)| keep(*taxid));
            !taxids.is_empty()
        });
    }

    type NameLookup = Vec<(NameClassSymbol, String)>;
    fn lookup_names<'a>(&'a self, taxid: TaxId) -> Option<&'a Self::NameLookup> {
        self.taxid_to_names.get(&taxid)
    }

    type NamesLookupIter<'a> = IterStrSnd<slice::Iter<'a, (NameClassSymbol, String)>>;
    fn iter_lookup_names<'a>(lookup: &'a Self::NameLookup) -> Self::NamesLookupIter<'a> {
        IterStrSnd(lookup.iter())
    }

    type TaxIdsLookup = Vec<(NameClassSymbol, TaxId)>;
    fn lookup_taxids<'a>(&'a self, name: &str) -> Option<&'a Self::TaxIdsLookup> {
        self.name_to_taxids.get(name)
    }

    type TaxIdsLookupIter<'a> = IterCopySnd<slice::Iter<'a, (NameClassSymbol, NodeId)>>;
    fn iter_lookup_taxids<'a>(lookup: &'a Self::TaxIdsLookup) -> Self::TaxIdsLookupIter<'a> {
        IterCopySnd(lookup.iter())
    }
}

pub struct IterStrSnd<I>(I);

impl<'a, T: 'a, I: Iterator<Item = &'a (T, String)>> Iterator for IterStrSnd<I> {
    type Item = &'a str;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(|s| &*s.1)
    }
}

pub struct IterCopySnd<I>(I);

impl<'a, T: 'a, U: 'a + Copy, I: Iterator<Item = &'a (T, U)>> Iterator for IterCopySnd<I> {
    type Item = U;

    fn next(&mut self) -> Option<Self::Item> {
        self.0.next().map(|(_, b)| *b)
    }
}

impl NamesAssoc for SingleClassNames {
    fn insert(&mut self, taxid: TaxId, name: &str, unique_name: &str, name_class: &str) {
        if name_class == self.name_class {
            self.taxid_to_name
                .entry(taxid)
                .or_insert_with(|| if unique_name.is_empty() { name } else { unique_name }.to_owned());

            self.name_to_taxids
                .entry(name.to_owned())
                .or_default()
                .push(taxid);
        }
    }

    fn forget_taxids<Keep>(&mut self, mut keep: Keep) where Keep: FnMut(TaxId) -> bool {
        self.taxid_to_name.retain(|taxid, _| keep(*taxid));

        self.name_to_taxids.retain(|_, taxids| {
            taxids.retain(|taxid| keep(*taxid));
            !taxids.is_empty()
        })
    }

    type NameLookup = String;
    fn lookup_names<'a>(&'a self, taxid: TaxId) -> Option<&'a Self::NameLookup> {
        self.taxid_to_name.get(&taxid)
    }

    type NamesLookupIter<'a> = std::iter::Once<&'a str>;
    fn iter_lookup_names<'a>(lookup: &'a Self::NameLookup) -> Self::NamesLookupIter<'a> {
        std::iter::once(&*lookup)
    }

    type TaxIdsLookup = Vec<TaxId>;
    fn lookup_taxids<'a>(&'a self, name: &str) -> Option<&'a Self::TaxIdsLookup> {
        self.name_to_taxids.get(name)
    }

    type TaxIdsLookupIter<'a> = std::iter::Cloned<std::slice::Iter<'a, TaxId>>;
    fn iter_lookup_taxids<'a>(lookup: &'a Self::TaxIdsLookup) -> Self::TaxIdsLookupIter<'a> {
        lookup.iter().cloned()
    }
}

impl NamesAssoc for NoNames {
    fn insert(&mut self, _taxid: TaxId, _name: &str, _unique_name: &str, _name_class: &str) {}

    fn forget_taxids<Keep>(&mut self, _keep: Keep) where Keep: FnMut(TaxId) -> bool {}

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

    fn fill_from<P, F>(&mut self, _names_dmp: P, _class_filter: F) -> Result<(), DmpError>
    where
        F: Fn(&str) -> bool,
        P: AsRef<Path>,
    {
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

impl<Names> From<NcbiTaxonomy<Names>> for GenericTaxonomy {
    fn from(taxonomy: NcbiTaxonomy<Names>) -> Self {
        taxonomy.tree
    }
}

impl From<NcbiTaxonomy<SingleClassNames>> for SimpleTree {
    fn from(taxonomy: NcbiTaxonomy<SingleClassNames>) -> Self {
        let root = taxonomy.tree.root;
        let mut children_lookup = taxonomy.tree.children_lookup;
        let mut names = taxonomy.names.taxid_to_name;

        make_simple_tree(root, &mut children_lookup, &mut names)
    }
}
