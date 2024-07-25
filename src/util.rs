use flate2::read::GzDecoder;

#[allow(dead_code)]
pub mod cli_tools;
pub mod csv_flatten_fix;
#[allow(dead_code)]
pub mod csv_stream;
pub mod filter;
#[allow(dead_code)]
pub mod interned_mapping;
#[allow(dead_code)]
pub mod progress_monitor;

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

pub(crate) use writing_new_file_or_stdout;

#[cfg(target_family = "unix")]
fn is_broken_pipe_impl(mut err: &(dyn std::error::Error + 'static)) -> bool {
    loop {
        if let Some(io_err) = err.downcast_ref::<std::io::Error>() {
            if io_err.kind() == std::io::ErrorKind::BrokenPipe {
                return true;
            }
        } else if let Some(cause) = err.source() {
            err = cause;
        } else if let Some(polars_err) = err.downcast_ref::<PolarsError>() {
            // currently (0.41) polars does not report IO { error } as .cause()
            if let PolarsError::IO { error: io_err, msg: _ } = polars_err {
                err = io_err.as_ref();
            }
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

#[allow(clippy::large_enum_variant, dead_code)]
pub enum MaybeGzDecoder<R> {
    GzDecoder(GzDecoder<R>),
    Reader(R),
}

#[allow(unused_macros)]
macro_rules! maybe_gzdecoder {
    ($reader_expr:expr, $reader_pat:pat => $body:expr $(,)?) => {{
        match $reader_expr {
            crate::util::MaybeGzDecoder::GzDecoder($reader_pat) => $body,
            crate::util::MaybeGzDecoder::Reader($reader_pat) => $body,
        }
    }}
}

#[allow(unused_imports)]
pub(crate) use maybe_gzdecoder;
use polars::prelude::PolarsError;

