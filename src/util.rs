pub mod blastout;
pub mod cli_tools;
pub mod csv_flatten_fix;
pub mod csv_stream;
pub mod filter;
pub mod interned_mapping;
pub mod progress_monitor;
pub mod vpf_class_record;

macro_rules! writing_new_file_or_stdout {
    ($path:expr, $writer:pat => $body:expr $(,)?) => {{
        let path = $path;

        if path == "-" {
            let $writer = Ok::<_, std::io::Error>(std::io::stdout());
            $body
        } else {
            let $writer = std::fs::File::create(path);
            $body
        }
    }};
}

use flate2::read::GzDecoder;
pub(crate) use writing_new_file_or_stdout;


#[allow(clippy::large_enum_variant)]
pub enum MaybeGzDecoder<R> {
    GzDecoder(GzDecoder<R>),
    Reader(R),
}

macro_rules! maybe_gzdecoder {
    ($reader_expr:expr, $reader_pat:pat => $body:expr $(,)?) => {{
        match $reader_expr {
            crate::util::MaybeGzDecoder::GzDecoder($reader_pat) => $body,
            crate::util::MaybeGzDecoder::Reader($reader_pat) => $body,
        }
    }}
}

pub(crate) use maybe_gzdecoder;
