use clap::arg_enum;
use structopt::{clap::ArgGroup, StructOpt};

use bstr::{io::*, BString};
use std::{
    fs::File,
    io::{BufReader, Read, Write},
    path::PathBuf,
};

use gfa::gfa::GFA;
use gfa::optfields::OptionalFields;
use gfa::parser::GFAParser;
use gfa::writer::gfa_string;
use handlegraph::hashgraph::HashGraph;

use gfautil::{edges, gaf_convert, subgraph};

arg_enum! {
    #[derive(StructOpt, Debug)]
    enum SubgraphBy {
        Paths,
        Segments,
    }
}

#[derive(StructOpt, Debug)]
#[structopt(group = ArgGroup::with_name("names").required(true))]
struct SubgraphArgs {
    subgraph_by: SubgraphBy,
    #[structopt(long, group = "names")]
    file: Option<PathBuf>,
    #[structopt(long, group = "names")]
    list: Option<Vec<String>>,
}

#[derive(StructOpt, Debug)]
struct GAF2PAFArgs {
    #[structopt(name = "path to GAF file", long = "gaf", parse(from_os_str))]
    gaf: PathBuf,
    #[structopt(name = "PAF output paf", short = "o", long = "paf")]
    out: Option<PathBuf>,
}

#[derive(StructOpt, Debug)]
enum Command {
    Subgraph(SubgraphArgs),
    EdgeCount,
    #[structopt(name = "gaf2paf")]
    Gaf2Paf(GAF2PAFArgs),
}

#[derive(StructOpt, Debug)]
struct Opt {
    #[structopt(
        name = "path to GFA input file",
        short,
        long = "gfa",
        parse(from_os_str)
    )]
    in_gfa: PathBuf,
    #[structopt(subcommand)]
    command: Command,
}

fn read_byte_lines<R: Read>(reader: R) -> Vec<Vec<u8>> {
    BufReader::new(reader)
        .byte_lines()
        .map(|l| l.unwrap())
        .collect()
}

fn main() {
    let opt = Opt::from_args();
    match opt.command {
        Command::Subgraph(subgraph_args) => {
            let parser = GFAParser::new();
            let gfa: GFA<BString, OptionalFields> =
                parser.parse_file(&opt.in_gfa).unwrap();
            let names: Vec<Vec<u8>> = if let Some(list) = subgraph_args.list {
                list.into_iter().map(|s| s.bytes().collect()).collect()
            } else if let Some(path) = subgraph_args.file {
                read_byte_lines(File::open(path).unwrap())
            } else {
                read_byte_lines(std::io::stdin())
            };

            let new_gfa = match subgraph_args.subgraph_by {
                SubgraphBy::Paths => subgraph::paths_new_subgraph(&gfa, &names),
                SubgraphBy::Segments => {
                    subgraph::segments_subgraph(&gfa, &names)
                }
            };
            println!("{}", gfa_string(&new_gfa));
        }
        Command::Gaf2Paf(args) => {
            let parser = GFAParser::new();
            let gfa: GFA<BString, OptionalFields> =
                parser.parse_file(&opt.in_gfa).unwrap();
            let paf_lines = gaf_convert::gaf_to_paf(&gfa, &args.gaf);

            if let Some(out_path) = args.out {
                let mut out_file = File::create(&out_path)
                    .expect("Error creating PAF output file");

                paf_lines.iter().for_each(|p| {
                    writeln!(out_file, "{}", p).unwrap();
                });
            } else {
                paf_lines.iter().for_each(|p| println!("{}", p));
            }
        }
        Command::EdgeCount => {
            let parser = GFAParser::new();
            let gfa: GFA<BString, ()> = parser.parse_file(&opt.in_gfa).unwrap();

            let hashgraph = HashGraph::from_gfa(&gfa);
            let edge_counts = edges::graph_edge_count(&hashgraph);
            println!("nodeid,inbound,outbound,total");
            edge_counts
                .iter()
                .for_each(|(id, i, o, t)| println!("{},{},{},{}", id, i, o, t));
        }
    }
}
