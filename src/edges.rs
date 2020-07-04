use gfa::gfa::GFA;

use handlegraph::handle::Direction;
use handlegraph::handlegraph::{handles_iter, HandleGraph};

pub fn graph_edge_count<T: HandleGraph>(graph: &T) -> Vec<(u64, usize, usize, usize)> {
    handles_iter(graph)
        .map(|h| {
            let inbound = graph.degree(h, Direction::Left);
            let outbound = graph.degree(h, Direction::Right);
            let total = inbound + outbound;
            let id: u64 = h.unpack_number();

            (id, inbound, outbound, total)
        })
        .collect()
}
