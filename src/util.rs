#[allow(dead_code)]
pub mod blastout;
#[allow(dead_code)]
pub mod cli_tools;
pub mod csv_flatten_fix;
#[allow(dead_code)]
pub mod csv_stream;
pub mod filter;
#[allow(dead_code)]
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

#[cfg(target_family = "unix")]
fn is_broken_pipe_impl(mut err: &(dyn std::error::Error + 'static)) -> bool {
    loop {
        if let Some(io_err) = err.downcast_ref::<std::io::Error>() {
            if io_err.kind() == std::io::ErrorKind::BrokenPipe {
                return true;
            }
        } else if let Some(cause) = err.source() {
            err = cause;
        } else {
            break;
        }
    }

    false
}

#[cfg(not(target_family = "unix"))]
fn is_broken_pipe_impl(result: &(dyn std::error::Error + 'static)) -> bool {
    false
}

pub fn is_broken_pipe<E>(err: &E) -> bool
where
    E: std::error::Error + 'static
{
    is_broken_pipe_impl(err)
}

pub fn ignore_broken_pipe<E>(result: Result<(), E>) -> Result<(), E>
where
    E: std::error::Error + 'static
{
    if let Err(err) = &result {
        if is_broken_pipe(err) {
            return Ok(())
        }
    }

    result
}

pub fn is_broken_pipe_anyhow<E>(err: &E) -> bool
where
    E: AsRef<dyn std::error::Error + 'static>,
{
    is_broken_pipe_impl(err.as_ref())
}

pub fn ignore_broken_pipe_anyhow<E>(result: Result<(), E>) -> Result<(), E>
where
    E: AsRef<dyn std::error::Error + 'static>,
{
    if let Err(err) = &result {
        if is_broken_pipe_anyhow(err) {
            return Ok(())
        }
    }

    result
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
