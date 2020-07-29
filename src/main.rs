use structopt::{clap::ArgGroup, StructOpt};

use std::fs::File;
use std::path::PathBuf;

use bstr::io::*;
use bstr::BString;
use clap::arg_enum;
use std::io;
use std::io::BufReader;

use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use gfa::writer::gfa_string;

use gfautil::edges;
use gfautil::subgraph;

use handlegraph::hashgraph::HashGraph;

arg_enum! {
    #[derive(StructOpt, Debug)]
    enum NameList {
        File,
        List,
        Stdin,
    }
}

arg_enum! {
    #[derive(StructOpt, Debug)]
    enum SubgraphBy {
        Paths,
        Segments,
    }
}

#[derive(StructOpt, Debug)]
#[structopt(group = ArgGroup::with_name("names").required(true))]
struct SubgraphCmd {
    subgraph_by: SubgraphBy,
    #[structopt(long, group = "names")]
    file: Option<PathBuf>,
    #[structopt(long, group = "names")]
    list: Option<Vec<String>>,
    // #[structopt(long, group = "names")]
    // stdin: bool,
    // input: NameList,
    // list_file: Option<PathBuf>,
    // #[structopt(required_unless("))
    // arg_list: Option<Vec<String>>,
}

#[derive(StructOpt, Debug)]
enum Command {
    // Subgraph { cmd: SubgraphBy, input: NameList },
    Subgraph(SubgraphCmd),
    // Subgraph { paths: Vec<String> },
    EdgeCount,
    // IterTest,
}

#[derive(StructOpt, Debug)]
struct Opt {
    #[structopt(short, long, parse(from_os_str))]
    in_file: PathBuf,
    #[structopt(flatten)]
    command: Command,
}

fn main() {
    let opt = Opt::from_args();
    println!("parsing GFA");
    let parser = GFAParser::new();
    let gfa = parser.parse_file(&opt.in_file).unwrap();
    println!("parsed GFA");
    match opt.command {
        // Command::IterTest => {
        //     gfautil::iter_test::iter_test(&gfa);
        // }
        // Command::Subgraph { paths } => {
        // Command::Subgraph { cmd, input } => {
        Command::Subgraph(sub_cmd) => {
            let names: Vec<Vec<u8>> = if let Some(path) = sub_cmd.file {
                let file = File::open(path).unwrap();
                let lines = BufReader::new(file)
                    .byte_lines()
                    .map(|l| l.unwrap())
                    .collect();
                lines
            } else if let Some(list) = sub_cmd.list {
                list.into_iter().map(|s| s.bytes().collect()).collect()
            } else {
                // else read from stdin
                let reader = BufReader::new(std::io::stdin());
                let lines = reader.byte_lines().map(|l| l.unwrap()).collect();
                lines
            };

            match sub_cmd.subgraph_by {
                SubgraphBy::Paths => {
                    let new_gfa = subgraph::paths_new_subgraph(&gfa, &names);
                    println!("{}", gfa_string(&new_gfa));
                }
                SubgraphBy::Segments => {
                    let new_gfa = subgraph::segments_subgraph(&gfa, &names);
                    println!("{}", gfa_string(&new_gfa));
                }
            }
        }
        Command::EdgeCount => {
            let hashgraph = HashGraph::from_gfa(&gfa);
            let edge_counts = edges::graph_edge_count(&hashgraph);
            println!("nodeid,inbound,outbound,total");
            edge_counts
                .iter()
                .for_each(|(id, i, o, t)| println!("{},{},{},{}", id, i, o, t));
        }
    }
}
