use indicatif::{ProgressBar, ProgressStyle};

pub(crate) fn progress_bar(len: usize, steady: bool) -> ProgressBar {
    let p_bar = ProgressBar::new(len as u64);
    p_bar.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:80} {pos:>7}/{len:7}")
            .progress_chars("##-"),
    );
    if steady {
        p_bar.enable_steady_tick(1000);
    }
    p_bar
}
