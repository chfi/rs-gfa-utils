use std::path::PathBuf;

use gfa::gfa::GFA;

use handlegraph::hashgraph::HashGraph;

use crate::edges;

use super::{load_gfa, Result};

pub fn edge_count(gfa_path: &PathBuf) -> Result<()> {
    let gfa: GFA<usize, ()> = load_gfa(gfa_path)?;

    let hashgraph = HashGraph::from_gfa(&gfa);
    let edge_counts = edges::graph_edge_count(&hashgraph);
    println!("nodeid,inbound,outbound,total");
    edge_counts
        .iter()
        .for_each(|(id, i, o, t)| println!("{},{},{},{}", id, i, o, t));

    Ok(())
}
