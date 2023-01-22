use std::{sync::{atomic::{AtomicU64, AtomicBool, self}, Arc}, io::{self, Write}, time::Duration, thread, panic};

use progress_streams::{ProgressWriter, ProgressReader};


#[derive(Default)]
pub struct TransferState {
    transferred: AtomicU64,
    complete: AtomicBool,
}

impl TransferState {
    pub fn is_complete(&self) -> bool {
        self.complete.load(atomic::Ordering::Acquire)
    }

    pub fn transferred(&self) -> u64 {
        self.transferred.load(atomic::Ordering::Acquire)
    }
}

pub fn default_progress_callback(file_size: u64) -> impl Fn(&TransferState) {
    move |transfer_state| {
        let transferred = transfer_state.transferred() as f64;
        let total = file_size as f64;
        let progress = transferred / total;

        print!(
            "\r  {:>6.02}%  {:>8.02}/{:.02} MiB",
            progress*100., transferred / (1024. * 1024.), total / (1024. * 1024.)
        );

        let _ = io::stdout().flush();
    }
}

pub fn default_finish_callback() -> impl Fn(&TransferState) {
    |_| { println!(); }
}

pub trait ReaderConsumer {
    type Output = io::Result<()>;

    fn read_from<R: io::Read>(self, reader: R) -> Self::Output;
}

pub fn monitor_reader<R, ProgressCallback, FinishCallback, Op>(
    reader: R,
    update_interval: Duration,
    mut progress_callback: ProgressCallback,
    finish_callback: FinishCallback,
    op: Op
) -> Op::Output
where
    R: io::Read,
    ProgressCallback: FnMut(&TransferState) + Send + 'static,
    FinishCallback: FnOnce(&TransferState) + Send + 'static,
    Op: ReaderConsumer,
{
    let state = Arc::new(TransferState::default());
    progress_callback(&state);

    let state_clone = Arc::clone(&state);

    let handle = thread::spawn(move || {
        while !state_clone.is_complete() {
            thread::sleep(update_interval);
            progress_callback(&state_clone);
        }
        finish_callback(&state_clone);
    });

    let progress_reader = ProgressReader::new(reader, |progress| {
        state
            .transferred
            .fetch_add(progress as u64, atomic::Ordering::Release);
    });

    let op_result = op.read_from(progress_reader);
    state.complete.store(true, atomic::Ordering::Release);

    if let Err(panic_msg) = handle.join() {
        panic::panic_any(panic_msg);
    }

    op_result
}


pub trait WriterConsumer {
    type Output = io::Result<()>;

    fn write_into<W: Write>(self, writer: W) -> Self::Output;
}

pub fn monitor_writer<W, ProgressCallback, FinishCallback, Op>(
    writer: W,
    update_interval: Duration,
    mut progress_callback: ProgressCallback,
    finish_callback: FinishCallback,
    op: Op
) -> Op::Output
where
    W: Write,
    ProgressCallback: FnMut(&TransferState) + Send + 'static,
    FinishCallback: FnOnce(&TransferState) + Send + 'static,
    Op: WriterConsumer,
{
    let state = Arc::new(TransferState::default());
    progress_callback(&state);

    let state_clone = Arc::clone(&state);

    let handle = thread::spawn(move || {
        while !state_clone.is_complete() {
            thread::sleep(update_interval);
            progress_callback(&state_clone);
        }
        finish_callback(&state_clone);
    });

    let progress_writer = ProgressWriter::new(writer, |progress| {
        state
            .transferred
            .fetch_add(progress as u64, atomic::Ordering::Release);
    });

    let op_result = op.write_into(progress_writer);
    state.complete.store(true, atomic::Ordering::Release);

    if let Err(panic_msg) = handle.join() {
        panic::panic_any(panic_msg);
    }

    op_result
}

struct CopyFrom<R: io::Read>(R);

impl<R: io::Read> WriterConsumer for CopyFrom<R> {
    type Output = io::Result<u64>;

    fn write_into<W: Write>(mut self, mut writer: W) -> Self::Output {
        io::copy(&mut self.0, &mut writer)
    }
}

pub fn monitor_copy<R, W, ProgressCallback, FinishCallback>(
    reader: R,
    writer: W,
    update_interval: Duration,
    progress_callback: ProgressCallback,
    finish_callback: FinishCallback,
) -> io::Result<u64>
where
    R: io::Read,
    W: io::Write,
    ProgressCallback: FnMut(&TransferState) + Send + 'static,
    FinishCallback: FnOnce(&TransferState) + Send + 'static,
{
    monitor_writer(writer, update_interval, progress_callback, finish_callback, CopyFrom(reader))
}
