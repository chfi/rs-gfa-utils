use structopt::StructOpt;

use std::path::PathBuf;

use gfautil::{
    commands,
    commands::{
        convert_names::GfaIdConvertOptions, gaf2paf::GAF2PAFArgs,
        gfa2vcf::GFA2VCFArgs, subgraph::SubgraphArgs, Result,
    },
};

#[derive(StructOpt, Debug)]
enum Command {
    Subgraph(SubgraphArgs),
    EdgeCount,
    #[structopt(name = "gaf2paf")]
    Gaf2Paf(GAF2PAFArgs),
    GfaSegmentIdConversion(GfaIdConvertOptions),
    #[structopt(name = "gfa2vcf")]
    Gfa2Vcf(GFA2VCFArgs),
}

#[derive(StructOpt, Debug)]
struct LogOpt {
    /// Show no messages.
    #[structopt(long)]
    quiet: bool,
    /// Show info messages.
    #[structopt(long)]
    info: bool,
    /// Show debug messages.
    #[structopt(long)]
    debug: bool,
}

#[derive(StructOpt, Debug)]
struct Opt {
    #[structopt(
        name = "input GFA file",
        short,
        long = "gfa",
        parse(from_os_str)
    )]
    in_gfa: PathBuf,
    #[structopt(subcommand)]
    command: Command,
    #[structopt(flatten)]
    log_opts: LogOpt,
    /// The number of threads to use when applicable. If omitted,
    /// Rayon's default will be used, based on the RAYON_NUM_THREADS
    /// environment variable, or the number of logical CPUs.
    #[structopt(short, long)]
    threads: Option<usize>,
}

fn init_logger(opt: &LogOpt) {
    let mut builder = pretty_env_logger::formatted_builder();
    if !opt.quiet {
        let mut log_level = log::LevelFilter::Error;
        if opt.info {
            log_level = log::LevelFilter::Info;
        }
        if opt.debug {
            log_level = log::LevelFilter::Debug;
        }
        builder.filter_level(log_level);
    }

    builder.init();
}

fn main() -> Result<()> {
    let opt = Opt::from_args();

    init_logger(&opt.log_opts);

    if let Some(threads) = &opt.threads {
        log::info!("Initializing threadpool to use {} threads", threads);
        rayon::ThreadPoolBuilder::new()
            .num_threads(*threads)
            .build_global()?;
    }

    match opt.command {
        Command::Gfa2Vcf(args) => {
            commands::gfa2vcf::gfa2vcf(&opt.in_gfa, &args)?;
        }

        Command::Subgraph(args) => {
            commands::subgraph::subgraph(&opt.in_gfa, &args)?;
        }
        Command::Gaf2Paf(args) => {
            commands::gaf2paf::gaf2paf(&opt.in_gfa, &args)?;
        }
        Command::EdgeCount => {
            commands::stats::edge_count(&opt.in_gfa)?;
        }
        Command::GfaSegmentIdConversion(args) => {
            commands::convert_names::convert_segment_ids(&opt.in_gfa, &args)?;
        }
    }
    Ok(())
}
