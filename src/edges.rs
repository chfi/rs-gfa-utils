use handlegraph::handlegraph::*;

/// Return the inbound and outbound edge counts for each node in the
/// graph
pub fn graph_edge_count<T: HandleGraphRef>(
    graph: T,
) -> Vec<(u64, usize, usize, usize)> {
    use handlegraph::handle::Direction;
    graph
        .handles()
        .map(|h| {
            let inbound = graph.degree(h, Direction::Left);
            let outbound = graph.degree(h, Direction::Right);
            let total = inbound + outbound;
            let id: u64 = h.unpack_number();

            (id, inbound, outbound, total)
        })
        .collect()
}
