use std::io;

use csv::StringRecord;
use lending_iterator::prelude::*;

pub struct CsvReaderIter<R> {
    buf_record: StringRecord,
    reader: csv::Reader<R>,
}

impl<R: io::Read> CsvReaderIter<R> {
    pub fn new(reader: csv::Reader<R>) -> Self {
        Self {
            buf_record: StringRecord::new(),
            reader,
        }
    }
}

#[gat]
impl<R: io::Read> LendingIterator for CsvReaderIter<R> {
    type Item<'next> = io::Result<&'next StringRecord>;

    fn next(&mut self) -> Option<Self::Item<'_>> {
        match self.reader.read_record(&mut self.buf_record) {
            Err(e) => {
                Some(Err(e.into()))
            }
            Ok(true) => {
                Some(Ok(&self.buf_record))
            }
            Ok(false) => {
                None
            }
        }
    }
}
