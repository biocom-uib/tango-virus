use std::{io::BufRead, marker::PhantomData, borrow::Cow};

use quick_xml::{Reader, events::{Event, BytesStart, BytesEnd}, reader::Span};


type XmlResult<T> = quick_xml::Result<T>;

pub trait EntryBuilder: Default {
    type Entry;

    // accession
    fn accession(&mut self, acc: &str);

    // name
    fn name(&mut self, name: &str);

    // sequence
    fn sequence(&mut self, sequence: &str);

    // dbreference
    fn begin_dbreference(&mut self, ref_type: &str, ref_id: &str);
    fn dbreference_property(&mut self, prop_type: &str, prop_value: &str);
    fn finish_dbreference(&mut self);

    // organism
    fn begin_organism(&mut self);

    fn organism_scientific_name(&mut self, name: &str);
    fn organism_common_name(&mut self, name: &str);

    //  organism/dbReference
    fn begin_organism_dbreference(&mut self, ref_type: &str, ref_id: &str);
    fn organism_dbreference_property(&mut self, prop_type: &str, prop_value: &str);
    fn finish_organism_dbreference(&mut self);

    fn finish_organism(&mut self);

    // finish building
    fn finish(self) -> Option<Self::Entry>;
}

pub struct UniProtXmlReader<R: BufRead> {
    xml_reader: Reader<R>,
    buffer: Vec<u8>,
}

pub mod dbref_types {
    pub const NCBI_TAXONOMY: &str = "NCBI Taxonomy";
}

#[allow(dead_code)]
mod attrs {
    pub const TYPE: &str  = "type";
    pub const ID: &str    = "id";
    pub const VALUE: &str = "value";

    pub mod u8 {
        pub const TYPE:  &[u8] = super::TYPE.as_bytes();
        pub const ID:    &[u8] = super::ID.as_bytes();
        pub const VALUE: &[u8] = super::VALUE.as_bytes();
    }
}

#[allow(dead_code)]
mod tags {
    pub const ACCESSION: &str   = "accession";
    pub const COMMON: &str      = "common";
    pub const DBREFERENCE: &str = "dbReference";
    pub const ENTRY: &str       = "entry";
    pub const NAME: &str        = "name";
    pub const ORGANISM: &str    = "organism";
    pub const PROPERTY: &str    = "property";
    pub const SCIENTIFIC: &str  = "scientific";
    pub const SEQUENCE: &str    = "sequence";
    pub const UNIPROT: &str     = "uniprot";

    pub mod u8 {
        pub const ACCESSION: &[u8]   = super::ACCESSION.as_bytes();
        pub const COMMON: &[u8]      = super::COMMON.as_bytes();
        pub const DBREFERENCE: &[u8] = super::DBREFERENCE.as_bytes();
        pub const ENTRY: &[u8]       = super::ENTRY.as_bytes();
        pub const NAME: &[u8]        = super::NAME.as_bytes();
        pub const ORGANISM: &[u8]    = super::ORGANISM.as_bytes();
        pub const PROPERTY: &[u8]    = super::PROPERTY.as_bytes();
        pub const SCIENTIFIC: &[u8]  = super::SCIENTIFIC.as_bytes();
        pub const SEQUENCE: &[u8]    = super::SEQUENCE.as_bytes();
        pub const UNIPROT: &[u8]     = super::UNIPROT.as_bytes();
    }
}

macro_rules! read_dbreference {
    ($self:ident, $builder:ident, $e:expr, $begin:ident, $add_prop:ident, $finish:ident) => {{
        UniProtXmlReader::begin_read_dbreference(&$self.xml_reader, $e, |ref_type, ref_id| {
            $builder.$begin(ref_type, ref_id)
        })?;

        $self.read_dbreference_properties(|prop_type, prop_value| {
            $builder.$add_prop(prop_type, prop_value)
        })?;

        $builder.$finish();
    }}
}

macro_rules! read_text {
    ($self:ident, $builder:ident, $tag:ident, $setter:ident) => {{
        let end = BytesEnd::new(tags::$tag);
        let span = $self
            .xml_reader
            .read_to_end_into(end.name(), &mut $self.buffer)?;
        $builder.accession(&$self.decode_span(span)?);
    }}
}

impl<R: BufRead> UniProtXmlReader<R> {
    pub fn new(reader: R) -> Self {
        let mut xml_reader = Reader::from_reader(reader);

        xml_reader
            .expand_empty_elements(true)
            .trim_text(true);

        Self {
            xml_reader,
            buffer: Vec::new(),
        }
    }

    pub fn decode_bytes<'a>(&self, bytes: &'a [u8]) -> XmlResult<Cow<'a, str>> {
        self.xml_reader.decoder().decode(bytes)
    }

    pub fn decode_span(&self, span: Span) -> XmlResult<Cow<'_, str>> {
        self.decode_bytes(&self.buffer[span])
    }

    pub fn into_entries_iter<EB>(self) -> EntriesIter<R, EB> {
        EntriesIter {
            reader: self,
            error: false,
            marker: PhantomData,
        }
    }

    fn find_next_entry_start(&mut self) -> XmlResult<bool> {
        loop {
            match self.xml_reader.read_event_into(&mut self.buffer)? {
                Event::DocType(_) => {},
                Event::Decl(_) => {},
                Event::Start(e) => {
                    match e.name().as_ref() {
                        tags::u8::UNIPROT => {},
                        tags::u8::ENTRY => return Ok(true),
                        name => {
                            eprintln!("Unknown tag found while scanning entries: {name:?}");
                        },
                    }
                },
                Event::End(e) => {
                    match e.name().as_ref() {
                        tags::u8::UNIPROT => {},
                        name => {
                            eprintln!("Unknown tag close found while scanning entries: {name:?}");
                        },
                    }
                }
                Event::Eof => return Ok(false),
                evt => {
                    eprintln!("Unexpected event while scanning entries: {evt:?}");
                }
            };
            self.buffer.clear();
        }
    }

    fn read_next_entry<EB: EntryBuilder>(&mut self) -> XmlResult<Option<EB::Entry>> {
        if self.find_next_entry_start()? {
            let mut builder = EB::default();
            self.read_entry(&mut builder)?;
            Ok(builder.finish())
        } else {
            Ok(None)
        }
    }

    fn read_entry<EB: EntryBuilder>(&mut self, builder: &mut EB) -> XmlResult<()> {
        loop {
            match self.xml_reader.read_event_into(&mut self.buffer)? {
                Event::Start(e) => {
                    match e.name().as_ref() {
                        tags::u8::ACCESSION => {
                            read_text!(self, builder, ACCESSION, accession);
                        }
                        tags::u8::NAME => {
                            read_text!(self, builder, NAME, name);
                        }
                        tags::u8::DBREFERENCE => {
                            read_dbreference!(self, builder, e,
                                begin_dbreference,
                                dbreference_property,
                                finish_dbreference
                            );
                        }
                        tags::u8::SEQUENCE => {
                            let end = BytesEnd::new(tags::SEQUENCE);
                            let span = self
                                .xml_reader
                                .read_to_end_into(end.name(), &mut self.buffer)?;
                            builder.sequence(&self.decode_span(span)?);
                        }
                        tags::u8::ORGANISM => self.read_organism(builder)?,
                        _ => {}
                    };
                }
                Event::Eof => break,
                Event::End(e) if e.name().as_ref() == tags::u8::ENTRY => break,
                evt => {
                    eprintln!("Unexpected event while reading <entry>: {evt:?}");
                }
            }

            self.buffer.clear();
        }

        Ok(())
    }

    fn read_organism<EB: EntryBuilder>(&mut self, builder: &mut EB) -> XmlResult<()> {
        builder.begin_organism();

        loop {
            match self.xml_reader.read_event_into(&mut self.buffer)? {
                Event::Start(e) => match e.name().as_ref() {
                    tags::u8::SCIENTIFIC => {
                        read_text!(self, builder, SCIENTIFIC, organism_scientific_name);
                    }
                    tags::u8::COMMON => {
                        read_text!(self, builder, COMMON, organism_common_name);
                    }
                    tags::u8::DBREFERENCE => {
                        read_dbreference!(self, builder, e,
                            begin_organism_dbreference,
                            organism_dbreference_property,
                            finish_organism_dbreference
                        );
                    }
                    _ => {}
                },
                Event::End(e) if e.name().as_ref() == b"organism" => {
                    builder.finish_organism();
                    break;
                }
                evt => {
                    eprintln!("Unexpected event while reading <organism>: {evt:?}",);
                }
            }

            self.buffer.clear();
        }

        Ok(())
    }

    fn begin_read_dbreference(
        xml_reader: &Reader<R>,
        start: BytesStart<'_>,
        begin: impl FnOnce(&str, &str),
    ) -> XmlResult<()> {

        let mut ref_type = None;
        let mut id = None;

        for attr in start.attributes() {
            let attr = attr?;

            match attr.key.as_ref() {
                b"type" => ref_type = Some(attr.value),
                b"id" => id = Some(attr.value),
                _ => {}
            };
        }

        if let (Some(ref_type), Some(id)) = (ref_type, id) {
            let ref_type = xml_reader.decoder().decode(&ref_type)?;
            let id = xml_reader.decoder().decode(&id)?;
            begin(&ref_type, &id);
        }

        Ok(())
    }

    fn read_dbreference_properties(
        &mut self,
        mut add_prop: impl FnMut(&str, &str),
    ) -> XmlResult<()> {
        let decoder = self.xml_reader.decoder();
        let property_end = BytesEnd::new(tags::PROPERTY);

        loop {
            match self.xml_reader.read_event_into(&mut self.buffer)? {
                Event::Start(e) => {
                    if e.name().as_ref() == b"property" {
                        let mut prop_type = None;
                        let mut prop_value = None;

                        for attr in e.attributes() {
                            let attr = attr?;

                            match attr.key.as_ref() {
                                attrs::u8::TYPE => prop_type = Some(attr.value),
                                attrs::u8::VALUE => prop_value = Some(attr.value),
                                _ => {}
                            };
                        }

                        if let (Some(prop_type), Some(prop_value)) = (prop_type, prop_value) {
                            let prop_type = decoder.decode(&prop_type)?;
                            let prop_value = decoder.decode(&prop_value)?;
                            add_prop(&prop_type, &prop_value);
                        }
                    }

                    self.xml_reader
                        .read_to_end_into(property_end.name(), &mut self.buffer)?;
                }
                Event::End(e) if e.name().as_ref() == tags::u8::DBREFERENCE => break,
                evt => {
                    eprintln!("Unexpected event while reading <dbReference>: {evt:?}");
                }
            }

            self.buffer.clear();
        }

        Ok(())
    }
}

pub struct EntriesIter<R: BufRead, EB> {
    reader: UniProtXmlReader<R>,
    error: bool,
    marker: PhantomData<fn() -> EB>,
}

impl<R, EB> Iterator for EntriesIter<R, EB>
where
    R: BufRead,
    EB: EntryBuilder,
{
    type Item = XmlResult<EB::Entry>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.error {
            return None;
        }

        let r = self.reader.read_next_entry::<EB>();

        if r.is_err() {
            self.error = true;
        }

        r.transpose()
    }
}
