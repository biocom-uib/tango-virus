macro_rules! writing_new_file_or_stdout {
    ($path:expr, $writer:pat => $body:expr $(,)?) => {{
        let path = $path;

        if path == "-" {
            let $writer = std::io::stdout();
            $body
        } else {
            let $writer = std::fs::File::create(path)?;
            $body
        }
    }};
}

pub(crate) use writing_new_file_or_stdout;
