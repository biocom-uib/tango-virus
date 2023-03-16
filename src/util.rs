pub mod blastout;
pub mod cli_tools;
pub mod csv_stream;
pub mod filter;
pub mod interned_mapping;
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
