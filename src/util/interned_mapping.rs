use std::{collections::{HashMap, HashSet}, hash::Hash, io};

use csv::StringRecord;
use lending_iterator::prelude::*;
use string_interner::{StringInterner, DefaultSymbol};

use super::csv_stream::CsvReaderIter;


pub trait Interner: Default {
    type Value<'a>: Copy;
    type Symbol: Copy + Eq + Ord + Hash;

    fn get_or_intern(&mut self, value: Self::Value<'_>) -> Self::Symbol;
    fn resolve(&self, sym: Self::Symbol) -> Option<Self::Value<'_>>;
}

impl Interner for StringInterner {
    type Value<'a> = &'a str;
    type Symbol = DefaultSymbol;

    fn get_or_intern(&mut self, value: &str) -> Self::Symbol {
        StringInterner::get_or_intern(self, value)
    }

    fn resolve(&self, symbol: Self::Symbol) -> Option<&str> {
        StringInterner::resolve(self, symbol)
    }
}


pub struct InternedMultiMapping<ValueInterner = StringInterner>
where
    ValueInterner: Interner,
{
    key_interner: StringInterner,
    value_interner: ValueInterner,

    mapping: HashMap<DefaultSymbol, HashSet<ValueInterner::Symbol>>,
}

impl<ValueInterner: Interner> Default for InternedMultiMapping<ValueInterner> {
    fn default() -> Self {
        Self {
            key_interner: Default::default(),
            value_interner: Default::default(),
            mapping: Default::default(),
        }
    }
}


impl<ValueInterner> InternedMultiMapping<ValueInterner>
where
    ValueInterner: Interner,
{
    pub fn add(&mut self, key: &str, value: ValueInterner::Value<'_>) {
        let value_sym = self.value_interner.get_or_intern(value);
        let key_sym = self.key_interner.get_or_intern(key);

        self.mapping.entry(key_sym).or_default().insert(value_sym);
    }

    #[apply(Gat!)]
    pub fn try_add_many<'a, Iter, E>(&'a mut self, mut iter: Iter) -> Result<(), E>
    where
        Iter: for<'n> LendingIterator<Item<'n> = Result<(&'n str, ValueInterner::Value<'n>), E>>,
    {
        let mut entry_cache_key = String::new();
        let mut entry_cache_values: Option<&mut HashSet<ValueInterner::Symbol>> = None;

        while let Some(result) = iter.next() {
            let (key, value) = result?;
            let value_sym = self.value_interner.get_or_intern(value);

            entry_cache_values = match entry_cache_values {
                Some(values) if entry_cache_key == key => {
                    values.insert(value_sym);
                    Some(values)
                },
                _ => {
                    let key_sym = self.key_interner.get_or_intern(key);

                    let values = self.mapping.entry(key_sym).or_default();
                    values.insert(value_sym);

                    entry_cache_key.replace_range(.., key);
                    Some(values)
                }
            };
        }

        Ok(())
    }

    pub fn lookup_syms(&self, key: &str) -> impl Iterator<Item = ValueInterner::Symbol> + '_ {
        self
            .key_interner
            .get(key)
            .and_then(|key_sym| self.mapping.get(&key_sym))
            .into_iter()
            .flatten()
            .copied()
    }

    pub fn lookup(&self, key: &str) -> impl Iterator<Item = ValueInterner::Value<'_>> + '_ {
        self.lookup_syms(key)
            .filter_map(|value_sym| self.value_interner.resolve(value_sym))
    }

    pub fn read_tsv_with<R, F>(reader: R, header: bool, record_parser: F) -> anyhow::Result<Self>
    where
        R: io::Read,
        F: for<'a> Fn(&'a StringRecord) -> anyhow::Result<(&'a str, ValueInterner::Value<'a>)>,
    {
        let csv_reader = csv::ReaderBuilder::new()
            .delimiter(b'\t')
            .has_headers(header)
            .from_reader(reader);

        let csv_iter =
            CsvReaderIter::new(csv_reader)
                .map::<HKT!(anyhow::Result<(&str, ValueInterner::Value<'_>)>), _>(|[], result| {
                    result.map_err(|e| e.into()).and_then(&record_parser)
                });

        let mut mapping = Self::default();

        mapping.try_add_many(csv_iter)?;

        Ok(mapping)
    }

    pub fn write_tsv_with<W, F>(
        &self,
        writer: W,
        header: &[&str],
        serializer_fn: F,
    ) -> io::Result<()>
    where
        W: io::Write,
        F: for<'a> Fn(&mut csv::Writer<W>, &'a str, ValueInterner::Value<'a>) -> io::Result<()>,
    {
        let mut csv_writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .from_writer(writer);

        csv_writer.write_record(header)?;

        for (&key_sym, value_syms) in &self.mapping {
            let key = self.key_interner.resolve(key_sym).unwrap();

            for &value_sym in value_syms {
                let value = self.value_interner.resolve(value_sym).unwrap();

                serializer_fn(&mut csv_writer, key, value)?;
            }
        }

        csv_writer.flush()?;
        Ok(())
    }
}


impl InternedMultiMapping<StringInterner> {
    pub fn read_tsv<R: io::Read>(reader: R, header: bool) -> anyhow::Result<Self> {
        Self::read_tsv_with(reader, header, |record| Ok((&record[0], &record[1])))
    }

    pub fn write_tsv<W: io::Write>(&self, writer: W, header: &[&str]) -> io::Result<()> {
        self.write_tsv_with(writer, header, |csv_writer, key, value| {
            csv_writer.write_record(&[key, value])?;
            Ok(())
        })
    }
}
