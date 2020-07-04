#[allow(dead_code)]
use gfa::gfa::GFA;

use std::collections::HashSet;

macro_rules! filtered {
    ($coll:expr, $pred:expr) => {
        $coll.iter().filter($pred).cloned().collect();
    };
}

/// Build a new GFA consisting of subgraphs of the given GFA
pub fn paths_new_subgraph(gfa: &GFA, paths: &[String]) -> GFA {
    let path_names: HashSet<&String> = paths.iter().collect();

    // Filter out the paths in the GFA we don't want
    let paths: Vec<_> =
        filtered!(gfa.paths, |p| path_names.contains(&p.path_name));

    // Set of the segments in the paths we're keeping
    let segment_names: HashSet<&str> = paths
        .iter()
        .flat_map(|path| path.segment_names.iter().map(|(seg, _)| seg.as_str()))
        .collect();

    // Filter out the segments in the GFA we don't want
    let segments =
        filtered!(gfa.segments, |s| segment_names.contains(s.name.as_str()));

    // Filter out the links in the GFA we don't want
    let links = filtered!(&gfa.links, |l| {
        segment_names.contains(l.from_segment.as_str())
            && segment_names.contains(l.to_segment.as_str())
    });

    let containments = filtered!(&gfa.containments, |l| {
        segment_names.contains(l.container_name.as_str())
            && segment_names.contains(l.contained_name.as_str())
    });

    GFA {
        segments,
        links,
        paths,
        containments,
    }
}

macro_rules! filtered_into {
    ($coll:expr, $pred:expr) => {
        $coll.into_iter().filter($pred).collect()
    };
}

// Consume the given GFA to create the subgraph GFA
pub fn paths_subgraph(gfa: GFA, paths: &[String]) -> GFA {
    let path_names: HashSet<&String> = paths.iter().collect();
    // Filter out the paths in the GFA we don't want
    let paths: Vec<_> =
        filtered_into!(gfa.paths, |p| path_names.contains(&p.path_name));

    // Set of the segments in the paths we're keeping
    let segment_names: HashSet<&str> = paths
        .iter()
        .flat_map(|path| path.segment_names.iter().map(|(seg, _)| seg.as_str()))
        .collect();

    // Filter out the segments in the GFA we don't want
    let segments = filtered_into!(gfa.segments, |s| segment_names
        .contains(s.name.as_str()));

    // Filter out the links in the GFA we don't want
    let links = filtered_into!(gfa.links, |l| {
        segment_names.contains(l.from_segment.as_str())
            && segment_names.contains(l.to_segment.as_str())
    });

    let containments = filtered_into!(gfa.containments, |l| {
        segment_names.contains(l.container_name.as_str())
            && segment_names.contains(l.contained_name.as_str())
    });

    GFA {
        segments,
        links,
        paths,
        containments,
    }
}
