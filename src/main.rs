use structopt::StructOpt;

use std::path::PathBuf;

use gfa::gfa::GFA;
use gfa::parser::parse_gfa;
use gfa::writer::gfa_string;

use gfautil::edges;
use gfautil::subgraph;

use handlegraph::hashgraph::HashGraph;

#[derive(StructOpt, Debug)]
enum Command {
    Subgraph { paths: Vec<String> },
    EdgeCount,
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
    let gfa = parse_gfa(&opt.in_file).unwrap();
    match opt.command {
        Command::Subgraph { paths } => {
            let new_gfa = subgraph::paths_new_subgraph(&gfa, &paths);
            println!("{}", gfa_string(&new_gfa));
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
