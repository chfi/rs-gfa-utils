use clap::arg_enum;
use structopt::{clap::ArgGroup, StructOpt};

use bstr::{io::*, BString, ByteSlice, ByteVec};
use std::{
    fs::File,
    io::{BufReader, Read, Write},
    path::PathBuf,
};

use gfa::{
    gfa::{name_conversion::NameMap, SegmentId, GFA},
    optfields::{OptFields, OptionalFields},
    parser::GFAParser,
    writer::{gfa_string, write_gfa},
};

use handlegraph::hashgraph::HashGraph;

use gfautil::{edges, subgraph};

use gfautil::commands;
use gfautil::commands::*;
use gfautil::commands::{
    convert_names::GfaIdConvertOptions, gaf2paf::GAF2PAFArgs,
    gfa2vcf::GFA2VCFArgs, subgraph::SubgraphArgs,
};

use log::{debug, info, warn};

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

    match opt.command {
        Command::Gfa2Vcf(args) => {
            let gfa: GFA<usize, ()> = commands::load_gfa(&opt.in_gfa)?;
            gfa2vcf::gfa2vcf(&opt.in_gfa, &gfa, &args)?;
        }

        Command::Subgraph(args) => {
            let gfa: GFA<BString, OptionalFields> =
                commands::load_gfa(&opt.in_gfa)?;
            commands::subgraph::subgraph(&gfa, &args)?;
        }
        Command::Gaf2Paf(args) => {
            let gfa: GFA<BString, OptionalFields> =
                commands::load_gfa(&opt.in_gfa)?;
            commands::gaf2paf::gaf2paf(gfa, &args)?;
        }
        Command::EdgeCount => {
            let gfa: GFA<usize, ()> = commands::load_gfa(&opt.in_gfa)?;

            let hashgraph = HashGraph::from_gfa(&gfa);
            let edge_counts = edges::graph_edge_count(&hashgraph);
            println!("nodeid,inbound,outbound,total");
            edge_counts
                .iter()
                .for_each(|(id, i, o, t)| println!("{},{},{},{}", id, i, o, t));
        }
        Command::GfaSegmentIdConversion(args) => {
            commands::convert_names::convert_segment_ids(opt.in_gfa, &args)?;
        }
    }
    Ok(())
}
