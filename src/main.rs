use clap::arg_enum;
use structopt::{clap::ArgGroup, StructOpt};

use bstr::{io::*, BString, ByteSlice, ByteVec};
use std::{
    fs::File,
    io::{BufRead, BufReader, Read, Write},
    path::PathBuf,
};

use gfa::gfa::GFA;
use gfa::optfields::OptionalFields;
use gfa::parser::GFAParser;
use gfa::writer::gfa_string;
use handlegraph::hashgraph::HashGraph;

use gfautil::{edges, gaf_convert, subgraph};

arg_enum! {
    #[derive(Debug, PartialEq)]
    enum SubgraphBy {
        Paths,
        Segments,
    }
}

/// Generate a subgraph of the input GFA.
///
/// The output will be the lines of the input GFA that include the
/// provided segment or path names.
#[derive(StructOpt, Debug)]
#[structopt(group = ArgGroup::with_name("names").required(true))]
struct SubgraphArgs {
    /// Choose between providing a list of path names, or a list of
    /// components of segment names
    #[structopt(name = "paths|segments", possible_values = &["paths", "segments"], case_insensitive = true)]
    subgraph_by: SubgraphBy,
    /// File containing a list of names
    #[structopt(
        name = "File containing names",
        long = "file",
        group = "names"
    )]
    file: Option<PathBuf>,
    /// Provide a list of names on the command line
    #[structopt(name = "List of names", long = "names", group = "names")]
    list: Option<Vec<String>>,
}

/*
#[derive(StructOpt, Debug)]
#[structopt(group = ArgGroup::with_name("names").required(true))]
struct InputOptions {
    #[structopt(name = "path to file", long, group = "names")]
    file: Option<PathBuf>,
    #[structopt(name = "list of path or segment names", long, group = "names")]
    list: Option<Vec<String>>,
}
*/

/// Convert a file of GAF records into PAF records.
///
/// The provided GFA file should be the same as the one used to create the GAF.
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
        name = "input GFA file",
        short,
        long = "gfa",
        parse(from_os_str)
    )]
    in_gfa: PathBuf,
    #[structopt(subcommand)]
    command: Command,
}

fn byte_lines_iter<'a, R: Read + 'a>(
    reader: R,
) -> Box<dyn Iterator<Item = Vec<u8>> + 'a> {
    Box::new(BufReader::new(reader).byte_lines().map(|l| l.unwrap()))
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
            } else {
                let in_lines = if let Some(path) = subgraph_args.file {
                    byte_lines_iter(File::open(path).unwrap())
                } else {
                    byte_lines_iter(std::io::stdin())
                };

                if subgraph_args.subgraph_by == SubgraphBy::Segments {
                    in_lines
                        .flat_map(|line| {
                            line.split_str("\t")
                                .map(|n| Vec::from_slice(n))
                                .collect::<Vec<_>>()
                        })
                        .collect()
                } else {
                    in_lines.collect()
                }
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
            let paf_lines = gaf_convert::gaf_to_paf(gfa, &args.gaf);

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
