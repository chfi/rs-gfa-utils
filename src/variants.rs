use handlegraph::{handle::*, handlegraph::*, hashgraph::HashGraph};

use fnv::{FnvHashMap, FnvHashSet};

use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};

use bstr::{BStr, BString, ByteSlice, ByteVec};

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
                    stack.push(node);
                    parents.insert(node, current);
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

/*
  Variant identification code from https://github.com/HopedWall/rs-gfatovcf
*/

pub fn find_all_paths_between(
    g: &HashGraph,
    start_node_id: &NodeId,
    end_node_id: &NodeId,
    max_edges: i32,
) -> Vec<Vec<NodeId>> {
    let mut all_paths_list: Vec<Vec<NodeId>> = Vec::new();

    // Put a limit on the maximum amount of edges that can be traversed
    // this should prevent eccessive memory usage
    // info!("Max edges is {:#?}", max_edges);
    let mut curr_edges = 0;
    let mut edges_limit_reached = false;

    // Keep a set of visited nodes so that loops are avoided
    let mut visited_node_id_set: FnvHashSet<NodeId> = FnvHashSet::default();

    // Create queue
    // NOTE: this is a Queue based implementation, this was done
    // in order not to get a stack overflow (the previous recursion-based
    // version was often experiencing this kind of issue)
    let mut q: VecDeque<NodeId> = VecDeque::new();

    // Insert first value
    q.push_back(*start_node_id);
    all_paths_list.push(vec![*start_node_id]);

    while !q.is_empty() {
        // info!("All paths is {:#?}", all_paths_list);
        // info!("Q is: {:#?}", q);

        let curr_node = q.pop_front().unwrap();
        // info!("Curr node is {:#?}", curr_node);

        if curr_node == *end_node_id {
            continue;
        }

        visited_node_id_set.insert(curr_node);
        let current_handle = Handle::pack(curr_node, false);

        // Get all paths that end in curr_node
        let mut curr_paths_list: Vec<_> = all_paths_list.clone();
        curr_paths_list.retain(|x| x.ends_with(&[curr_node]));

        // Only keep those which don't
        all_paths_list.retain(|x| !x.ends_with(&[curr_node]));

        // info!("Curr_paths_list: {:#?}", curr_paths_list);
        //io::stdin().read_line(&mut String::new());

        for neighbor in g.handle_edges_iter(current_handle, Direction::Right) {
            // info!("Neighbor: {:#?}", neighbor.id());
            // Append, for each current_path, this neighbor
            let mut temp = curr_paths_list.clone();
            temp.iter_mut().for_each(|x| x.push(neighbor.id()));
            all_paths_list.append(&mut temp);

            // Add new node to queue
            if !visited_node_id_set.contains(&neighbor.id())
                && !q.contains(&neighbor.id())
            {
                q.push_back(neighbor.id());
            }

            // Break if too many edges have been visited
            curr_edges += 1;
            if curr_edges > max_edges {
                edges_limit_reached = true;
                break;
            }
        }

        if edges_limit_reached {
            break;
        }

        // info!("All_paths_list: {:#?}", all_paths_list);
        //io::stdin().read_line(&mut String::new());
    }

    // Only keep paths that end in end_node_id
    // start_node_id does not have to be checked
    // TODO: maybe not needed?
    all_paths_list.retain(|x| x.ends_with(&[*end_node_id]));

    // info!(
    //     "All paths between {} and {} are: {:#?}",
    //     start_node_id, end_node_id, all_paths_list
    // );

    //io::stdin().read_line(&mut String::new());

    all_paths_list
}

/// A struct that holds Variants, as defined in the VCF format
#[derive(Debug, PartialEq)]
pub struct Variant {
    chromosome: BString,
    position: i32,
    id: Option<BString>,
    reference: BString,
    alternate: Option<BString>,
    quality: Option<i32>,
    filter: Option<BString>,
    info: Option<BString>,
    format: Option<BString>,
    sample_name: Option<BString>,
}

/// Detects variants from a list of bubbles
pub fn detect_all_variants(
    path_to_steps_map: &HashMap<BString, Vec<BString>>,
    possible_bubbles_list: &[(NodeId, NodeId)],
    graph: &HashGraph,
    node_id_to_path_and_pos_map: &BTreeMap<NodeId, HashMap<BString, usize>>,
    verbose: bool,
    max_edges: i32,
    reference_paths: &[BString],
) -> Vec<Variant> {
    let mut stuff_to_alts_map: HashMap<BString, HashSet<BString>> =
        HashMap::new();

    // For each reference path, explore all bubbles in order to find variants;
    // these will be stored in stuff_to_alts_map
    for current_ref in reference_paths {
        // Obtain all steps for current_ref
        let ref_path: Vec<u64> = path_to_steps_map[current_ref]
            .iter()
            .map(|x| {
                let s = x.to_str().unwrap();
                s.parse::<u64>().unwrap()
            })
            .collect();

        if verbose {
            println!("path_to_steps_map: {:?}", path_to_steps_map);
        }

        // info!("BEFORE DETECT");

        // Loop through all bubbles in order to find all variants
        // for the given reference
        detect_variants_per_reference(
            &current_ref,
            &ref_path,
            possible_bubbles_list,
            graph,
            node_id_to_path_and_pos_map,
            &mut stuff_to_alts_map,
            verbose,
            max_edges,
        );

        // info!("AFTER DETECT");
    }

    // Convert stuff_to_alts_map to a more readable format
    let mut vcf_list: Vec<Variant> = Vec::new();
    for (chrom_pos_ref, alt_type_set) in &stuff_to_alts_map {
        let vec: Vec<_> = chrom_pos_ref.split_str("_").collect();
        // let vec: Vec<&[u8]> = chrom_pos_ref.split('_').collect();
        let chrom = vec[0];
        let pos = vec[1];
        let refr = vec[2];

        let (alt_list, type_set): (Vec<_>, Vec<_>) = alt_type_set
            .iter()
            .map(|x| {
                let split: Vec<_> = x.split_str("_").collect();
                (split[0], split[1])
            })
            .unzip();

        let alts = alt_list.join(&b","[..]);
        let mut types: BString = "TYPE=".into();
        types.extend_from_slice(&type_set.join(&b";TYPE="[..]));
        // types.push_str(&type_set.join(&b";TYPE="[..]));

        let pos = pos.to_str().unwrap();
        let pos = pos.parse().unwrap();

        let v = Variant {
            chromosome: chrom.into(),
            position: pos,
            id: None,
            reference: refr.into(),
            alternate: Some(alts.into()),
            quality: None,
            filter: None,
            info: Some(types),
            format: Some("GT".into()),
            sample_name: Some("0|1".into()),
        };

        vcf_list.push(v);
    }

    // Sort vcf_list for printing variants in the correct order
    vcf_list.sort_by(|a, b| match a.chromosome.cmp(&b.chromosome) {
        std::cmp::Ordering::Equal => a.position.cmp(&b.position),
        other => other,
    });

    vcf_list
}
/// Detect variants for a specific reference
fn detect_variants_per_reference(
    current_ref: &[u8],
    ref_path: &[u64],
    possible_bubbles_list: &[(NodeId, NodeId)],
    graph: &HashGraph,
    node_id_to_path_and_pos_map: &BTreeMap<NodeId, HashMap<BString, usize>>,
    stuff_to_alts_map: &mut HashMap<BString, HashSet<BString>>,
    verbose: bool,
    max_edges: i32,
) {
    // info!("BEFORE GET LAST");

    // Create closure that will be used later
    let get_last = |prec_node_seq_ref: &[u8], node_seq_ref: &[u8]| {
        let last: &[u8] = &prec_node_seq_ref[prec_node_seq_ref.len() - 1..];
        let mut last: Vec<u8> = last.into();
        last.extend(node_seq_ref);
        last
    };

    // Check all bubbles
    for &(start, end) in possible_bubbles_list {
        if verbose {
            println!("ref_path: {:?}", ref_path);
            println!("Bubble [{},{}]", start, end);
        }

        // info!("BEFORE FIND START");

        let start_node_index_in_ref_path: usize;
        match ref_path.iter().position(|&r| NodeId::from(r) == start) {
            None => continue, //ignore, start not found in ref path
            Some(r) => start_node_index_in_ref_path = r,
        };

        // info!("BEFORE FIND ALL PATHS BETWEEN");

        let all_path_list: Vec<Vec<NodeId>> =
            find_all_paths_between(&graph, &start, &end, max_edges);

        // info!("AFTER FIND ALL PATHS BETWEEN");

        // info!("All paths list: {:?}", all_path_list);
        for path in &all_path_list {
            if verbose {
                println!("\tPath: {:?}", path);
            }

            //println!("INSIDE FOR LOOP");

            let x: &[u8] = current_ref;
            let path_map: &HashMap<BString, usize> =
                node_id_to_path_and_pos_map.get(&start).unwrap();

            // let ucurr
            let x: &BStr = x.as_ref();
            let pos_ref = path_map.get(x).unwrap();

            let mut pos_ref = pos_ref + 1;

            // let mut pos_ref =
            //     node_id_to_path_and_pos_map[&start][current_ref.as_ref()] + 1;
            let mut pos_path = pos_ref;

            let max_index = std::cmp::min(path.len(), ref_path.len());

            let mut current_index_step_path = 0;
            let mut current_index_step_ref = 0;

            for _i in 0..max_index {
                //Check if ref_path goes out of bounds
                //TODO: check how paths_to_steps is created, there may be some problems there
                // since ref_path is obtained from paths_to_steps
                if current_index_step_ref + start_node_index_in_ref_path
                    >= ref_path.len()
                {
                    continue;
                }

                let mut current_node_id_ref = NodeId::from(
                    ref_path
                        [current_index_step_ref + start_node_index_in_ref_path],
                );
                let mut current_node_id_path = path[current_index_step_path];

                if verbose {
                    println!(
                        "{} {} ---> {} {}",
                        pos_ref,
                        pos_path,
                        current_node_id_ref,
                        current_node_id_path
                    );
                }

                if current_node_id_ref == current_node_id_path {
                    if verbose {
                        println!("REFERENCE");
                    }

                    let node_seq = graph
                        .sequence(Handle::pack(current_node_id_ref, false));
                    pos_ref += node_seq.len();
                    pos_path = pos_ref;

                    current_index_step_ref += 1;
                    current_index_step_path += 1;
                } else {
                    // Shouldn't be happening anymore
                    if current_index_step_path + 1 >= path.len() {
                        break;
                    }
                    if current_index_step_ref + start_node_index_in_ref_path + 1
                        >= ref_path.len()
                    {
                        break;
                    }

                    let succ_node_id_path = path[current_index_step_path + 1];
                    let succ_node_id_ref = NodeId::from(
                        ref_path[current_index_step_ref
                            + start_node_index_in_ref_path
                            + 1],
                    );
                    if succ_node_id_ref == current_node_id_path {
                        if verbose {
                            println!("DEL");
                        }

                        let node_seq_ref = graph
                            .sequence(Handle::pack(current_node_id_ref, false));

                        let prec_node_id_ref = NodeId::from(
                            ref_path[current_index_step_ref
                                + start_node_index_in_ref_path
                                - 1],
                        );
                        let prec_nod_seq_ref = graph
                            .sequence(Handle::pack(prec_node_id_ref, false));

                        let last = get_last(&prec_nod_seq_ref, &node_seq_ref);

                        let pos_path_bs = (pos_path - 1).to_string();
                        let pos_path_bs = Vec::from(pos_path_bs.as_bytes());
                        let key = [current_ref.into(), pos_path_bs, last]
                            .join(&b"_"[..]);

                        let key = key.into();

                        stuff_to_alts_map
                            .entry(key)
                            .or_insert_with(HashSet::new);
                        //TODO: find a better way to do this
                        let last = get_last(&prec_nod_seq_ref, &node_seq_ref);

                        let pos_path_bs = (pos_path - 1).to_string();
                        let pos_path_bs = Vec::from(pos_path_bs.as_bytes());
                        let key = [current_ref.into(), pos_path_bs, last]
                            .join(&b"_"[..]);

                        let key: BString = key.into();

                        let last =
                            &prec_nod_seq_ref[prec_nod_seq_ref.len() - 1..];
                        // .to_string();
                        let mut string_to_insert = Vec::from(last);
                        string_to_insert.extend(b"_del");
                        stuff_to_alts_map
                            .get_mut(&key)
                            .unwrap()
                            .insert(string_to_insert.into());

                        pos_ref += node_seq_ref.len();

                        current_index_step_ref += 1;
                        current_node_id_ref = NodeId::from(
                            ref_path[current_index_step_ref
                                + start_node_index_in_ref_path
                                - 1],
                        );
                        if verbose {
                            println!("\t {}", current_node_id_ref);
                        }

                        continue;
                    } else if succ_node_id_path == current_node_id_ref {
                        if verbose {
                            println!("INS");
                        }

                        let node_seq_path = graph.sequence(Handle::pack(
                            current_node_id_path,
                            false,
                        ));

                        let prec_node_id_ref = NodeId::from(
                            ref_path[current_index_step_ref
                                + start_node_index_in_ref_path
                                - 1],
                        );
                        let prec_nod_seq_ref = graph
                            .sequence(Handle::pack(prec_node_id_ref, false));

                        let last = Vec::from(
                            &prec_nod_seq_ref[prec_nod_seq_ref.len() - 1..],
                        );
                        //let key = [current_ref.to_string(), (pos_ref-1).to_string(), String::from(prec_nod_seq_ref)].join("_");
                        let key = [
                            current_ref.into(),
                            Vec::from((pos_ref - 1).to_string().as_bytes()),
                            last,
                        ]
                        .join(&b"_"[..]);

                        stuff_to_alts_map
                            .entry(key.into())
                            .or_insert_with(HashSet::new);

                        //Re-create key since it goes out of scope
                        let last = Vec::from(
                            &prec_nod_seq_ref[prec_nod_seq_ref.len() - 1..],
                        );

                        let key = [
                            current_ref.into(),
                            Vec::from((pos_ref - 1).to_string().as_bytes()),
                            last,
                        ]
                        .join(&b"_"[..]);

                        let last = prec_nod_seq_ref
                            [prec_nod_seq_ref.len() - 1..]
                            .into();

                        let mut string_to_insert: Vec<u8> = last;
                        string_to_insert.push_str(&node_seq_path);
                        // string_to_insert.extend_from_slice(&node_seq_path);
                        string_to_insert.push_str("_ins");
                        // string_to_insert.push_

                        let key: BString = key.into();
                        stuff_to_alts_map
                            .get_mut(&key)
                            .unwrap()
                            .insert(string_to_insert.into());

                        pos_path += node_seq_path.len();

                        current_index_step_path += 1;
                        current_node_id_path = path[current_index_step_path];
                        if verbose {
                            println!("\t{}", current_node_id_path);
                        }

                        continue;
                    } else {
                        let node_seq_ref = graph
                            .sequence(Handle::pack(current_node_id_ref, false));
                        let node_seq_path = graph.sequence(Handle::pack(
                            current_node_id_path,
                            false,
                        ));

                        if node_seq_ref == node_seq_path {
                            if verbose {
                                println!("REFERENCE");
                            }
                        } else {
                            if verbose {
                                println!("SNV");
                            }
                        }

                        let key: BString = [
                            current_ref,
                            pos_path.to_string().as_bytes(),
                            node_seq_ref.as_bytes(),
                        ]
                        .join(&b"_"[..])
                        .into();

                        stuff_to_alts_map
                            .entry(key.into())
                            .or_insert_with(HashSet::new);

                        //TODO: find a better way to do this
                        let key: BString = [
                            current_ref,
                            pos_path.to_string().as_bytes(),
                            node_seq_ref.as_bytes(),
                        ]
                        .join(&b"_"[..])
                        .into();

                        let mut string_to_insert =
                            node_seq_path.chars().last().unwrap().to_string();
                        string_to_insert.push_str("_snv");
                        stuff_to_alts_map
                            .get_mut(&key)
                            .unwrap()
                            .insert(string_to_insert.into());

                        pos_ref += node_seq_ref.len();
                        pos_path += node_seq_path.len();
                        current_index_step_ref += 1;
                        current_index_step_path += 1;
                    }
                }
            }
            if verbose {
                println!("---");
            }
        }
    }
    if verbose {
        println!("==========================================");
    }
}
