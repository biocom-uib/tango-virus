use std::{collections::{HashMap, HashSet}, hash::Hash, io, marker::PhantomData};

use csv::StringRecord;
use lending_iterator::prelude::*;
use string_interner::{backend::Backend, DefaultStringInterner, DefaultSymbol, StringInterner};

use super::csv_stream::CsvReaderIter;


pub trait Interner: Default {
    type Value<'a>: Copy;
    type Symbol: Copy + Eq + Ord + Hash;

    fn get(&self, value: Self::Value<'_>) -> Option<Self::Symbol>;
    fn get_or_intern(&mut self, value: Self::Value<'_>) -> Self::Symbol;
    fn resolve(&self, sym: Self::Symbol) -> Option<Self::Value<'_>>;
}

impl<B: Backend> Interner for StringInterner<B>
where
    StringInterner<B>: Default,
    B::Symbol: Ord + Hash,
{
    type Value<'a> = &'a str;
    type Symbol = B::Symbol;

    fn get(&self, value: &str) -> Option<Self::Symbol> {
        StringInterner::get(self, value)
    }

    fn get_or_intern(&mut self, value: &str) -> Self::Symbol {
        StringInterner::get_or_intern(self, value)
    }

    fn resolve(&self, symbol: Self::Symbol) -> Option<&str> {
        StringInterner::resolve(self, symbol)
    }
}

pub struct TrivialInterner<T: Copy + Eq + Ord + Hash>(PhantomData<T>);

impl<T: Copy + Eq + Ord + Hash> Default for TrivialInterner<T> {
    fn default() -> Self {
        Self(PhantomData)
    }
}

impl<T: Copy + Eq + Ord + Hash> Interner for TrivialInterner<T> {
    type Value<'a> = T;

    type Symbol = T;

    fn get(&self, value: Self::Value<'_>) -> Option<Self::Symbol> {
        Some(value)
    }

    fn get_or_intern(&mut self, value: Self::Value<'_>) -> Self::Symbol {
        value
    }

    fn resolve(&self, sym: Self::Symbol) -> Option<Self::Value<'_>> {
        Some(sym)
    }
}

#[derive(Default)]
pub struct PairInterner<I1: Interner = DefaultStringInterner, I2: Interner = I1>(pub I1, pub I2);

impl<I1: Interner, I2: Interner> Interner for PairInterner<I1, I2> {
    type Value<'a> = (I1::Value<'a>, I2::Value<'a>);

    type Symbol = (I1::Symbol, I2::Symbol);

    fn get(&self, value: Self::Value<'_>) -> Option<Self::Symbol> {
        Some((self.0.get(value.0)?, self.1.get(value.1)?))
    }

    fn get_or_intern(&mut self, value: Self::Value<'_>) -> Self::Symbol {
        (self.0.get_or_intern(value.0), self.1.get_or_intern(value.1))
    }

    fn resolve(&self, sym: Self::Symbol) -> Option<Self::Value<'_>> {
        Some((self.0.resolve(sym.0)?, self.1.resolve(sym.1)?))
    }
}

pub struct InternedMultiMapping<ValueInterner = DefaultStringInterner>
where
    ValueInterner: Interner,
{
    key_interner: DefaultStringInterner,
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

pub trait IntoValues<'a, I: Interner>: IntoIterator<Item = I::Value<'a>> {}

impl<'a, I, Iter> IntoValues<'a, I> for Iter
where
    I: Interner,
    Iter: IntoIterator<Item = I::Value<'a>>,
{
}

impl<ValueInterner> InternedMultiMapping<ValueInterner>
where
    ValueInterner: Interner,
{
    pub fn key_interner(&self) -> &DefaultStringInterner {
        &self.key_interner
    }

    pub fn value_interner(&self) -> &ValueInterner {
        &self.value_interner
    }

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

    #[apply(Gat!)]
    pub fn try_add_many_grouped<VS: HKT, E, Iter>(&mut self, mut iter: Iter) -> Result<(), E>
    where
        Iter: for<'n> LendingIterator<Item<'n> = Result<(&'n str, Feed<'n, VS>), E>>,
        for<'n> Feed<'n, VS>: IntoValues<'n, ValueInterner>,
    {
        let mut entry_cache_key = String::new();
        let mut entry_cache_values: Option<&mut HashSet<ValueInterner::Symbol>> = None;

        while let Some(result) = iter.next() {
            let (key, values) = result?;

            let value_syms = values
                .into_iter()
                .map(|value| self.value_interner.get_or_intern(value));

            entry_cache_values = match entry_cache_values {
                Some(values) if entry_cache_key == key => {
                    values.extend(value_syms);
                    Some(values)
                },
                _ => {
                    let key_sym = self.key_interner.get_or_intern(key);

                    let values = self.mapping.entry(key_sym).or_default();
                    values.extend(value_syms);

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
            .comment(Some(b'#'))
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

    pub fn read_grouped_tsv_with<VS: HKT, R, F>(reader: R, header: bool, record_parser: F) -> anyhow::Result<Self>
    where
        R: io::Read,
        F: Fn(&'_ StringRecord) -> anyhow::Result<(&'_ str, Feed<'_, VS>)>,
        for<'a> Feed<'a, VS>: IntoValues<'a, ValueInterner>,
    {
        let csv_reader = csv::ReaderBuilder::new()
            .comment(Some(b'#'))
            .delimiter(b'\t')
            .has_headers(header)
            .from_reader(reader);

        let csv_iter =
            CsvReaderIter::new(csv_reader)
                .map::<HKT!(anyhow::Result<(&str, Feed<'_, VS>)>), _>(|[], result| {
                    result.map_err(|e| e.into()).and_then(&record_parser)
                });

        let mut mapping = Self::default();

        mapping.try_add_many_grouped(csv_iter)?;

        Ok(mapping)
    }

    pub fn write_tsv_with<W, F>(
        &self,
        writer: W,
        header: &[&str],
        serializer_fn: F,
    ) -> csv::Result<()>
    where
        W: io::Write,
        F: for<'a> Fn(&mut csv::Writer<W>, &'a str, ValueInterner::Value<'a>) -> csv::Result<()>,
    {
        let mut csv_writer = csv::WriterBuilder::new()
            .delimiter(b'\t')
            .has_headers(false)
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


impl InternedMultiMapping<DefaultStringInterner> {
    pub fn read_tsv<R: io::Read>(reader: R, header: bool) -> anyhow::Result<Self> {
        Self::read_tsv_with(reader, header, |record| Ok((&record[0], &record[1])))
    }

    pub fn write_tsv<W: io::Write>(&self, writer: W, header: &[&str]) -> csv::Result<()> {
        self.write_tsv_with(writer, header, |csv_writer, key, value| {
            csv_writer.write_record(&[key, value])?;
            Ok(())
        })
    }
}
