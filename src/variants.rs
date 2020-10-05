use handlegraph::{handle::*, handlegraph::*};

use fnv::{FnvHashMap, FnvHashSet};

// Finds all the nodes between two given nodes
pub fn extract_bubble_nodes<T>(
    graph: &T,
    from: NodeId,
    to: NodeId,
) -> FnvHashSet<NodeId>
where
    T: HandleGraph,
{
    let mut visited: FnvHashSet<NodeId> = FnvHashSet::default();
    let mut parents: FnvHashMap<NodeId, NodeId> = FnvHashMap::default();
    let mut stack: Vec<NodeId> = Vec::new();

    stack.push(from);

    while let Some(current) = stack.pop() {
        if !visited.contains(&current) {
            visited.insert(current);

            let handle = Handle::pack(current, false);

            if current != to {
                let neighbors =
                    graph.handle_edges_iter(handle, Direction::Right);

                for h in neighbors {
                    let node = h.id();
                    if !visited.contains(&node) {
                        stack.push(node);
                        parents.insert(node, current);
                    }
                }
            }
        }
    }

    visited
}

pub fn extract_paths_between_nodes<T>(
    graph: &T,
    from: NodeId,
    to: NodeId,
) -> FnvHashSet<Vec<NodeId>>
where
    T: HandleGraph,
{
    let mut visited: FnvHashSet<NodeId> = FnvHashSet::default();
    let mut parents: FnvHashMap<NodeId, NodeId> = FnvHashMap::default();
    let mut stack: Vec<NodeId> = Vec::new();

    let mut paths = FnvHashSet::default();

    stack.push(from);

    while let Some(current) = stack.pop() {
        if !visited.contains(&current) {
            visited.insert(current);

            let handle = Handle::pack(current, false);

            if current != to {
                let neighbors =
                    graph.handle_edges_iter(handle, Direction::Right);

                for h in neighbors {
                    let node = h.id();
                    if !visited.contains(&node) || node == to {
                        stack.push(node);
                        parents.insert(node, current);
                    }
                }
            } else {
                let mut cur_step = to;
                let mut cur_path = Vec::new();
                while cur_step != from {
                    cur_path.push(cur_step);
                    cur_step = *parents.get(&cur_step).unwrap();
                }
                cur_path.push(from);
                paths.insert(cur_path);
            }
        } else {
            let mut cur_step = current;
            let mut cur_path = Vec::new();
            while cur_step != from {
                cur_path.push(cur_step);
                cur_step = *parents.get(&cur_step).unwrap();
            }
            cur_path.push(from);
            paths.insert(cur_path);
        }
    }

    paths
}
