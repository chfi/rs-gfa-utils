use gfa::gfa::GFA;

use handlegraph::handle::Direction;
use handlegraph::handlegraph::{handle_iter, HandleGraph};

pub fn graph_edge_count<T: HandleGraph>(graph: &T) -> Vec<(u64, usize, usize, usize)> {
    handle_iter(graph)
        .map(|h| {
            let inbound = graph.get_degree(&h, Direction::Left);
            let outbound = graph.get_degree(&h, Direction::Right);
            let total = inbound + outbound;
            let id: u64 = h.unpack_number();

            (id, inbound, outbound, total)
        })
        .collect()
}
