#[allow(dead_code)]
use gfa::gfa::GFA;
use gfa::optfields::OptFields;

use std::collections::HashSet;

use bstr::{BStr, BString};

macro_rules! filtered {
    ($coll:expr, $pred:expr) => {
        $coll.iter().filter($pred).cloned().collect();
    };
}

/// Build a new GFA consisting of subgraphs of the given GFA that only
/// include segments that are in the provided paths. Paths are to be
/// provided as a slice of path names
pub fn paths_new_subgraph<T: OptFields + Clone>(
    gfa: &GFA<BString, T>,
    paths: &[Vec<u8>],
) -> GFA<BString, T> {
    let path_names: HashSet<&[u8]> =
        paths.into_iter().map(|p| p.as_ref()).collect();

    // Filter out the paths in the GFA we don't want
    let paths: Vec<_> =
        filtered!(gfa.paths, |p| path_names.contains(p.path_name.as_slice()));

    // Set of the segments in the paths we're keeping
    let segment_names: HashSet<&[u8]> = paths
        .iter()
        .flat_map(|path| path.iter().map(|(seg, _)| seg.as_ref()))
        // .flat_map(|path| path.segment_names.iter().map(|(seg, _)| seg.as_str()))
        .collect();

    // Filter out the segments in the GFA we don't want
    let segments =
        filtered!(gfa.segments, |s| segment_names.contains(s.name.as_slice()));

    // Filter out the links in the GFA we don't want
    let links = filtered!(&gfa.links, |l| {
        segment_names.contains(l.from_segment.as_slice())
            && segment_names.contains(l.to_segment.as_slice())
    });

    let containments = filtered!(&gfa.containments, |l| {
        segment_names.contains(l.container_name.as_slice())
            && segment_names.contains(l.contained_name.as_slice())
    });

    GFA {
        header: gfa.header.clone(),
        segments,
        links,
        paths,
        containments,
    }
}

/// Returns a subgraph GFA that only contains elements with the
/// provided segment names
pub fn segments_subgraph<T: OptFields + Clone>(
    gfa: &GFA<BString, T>,
    segment_names: &[Vec<u8>],
    // path_names: &[&[u8]]
) -> GFA<BString, T> {
    let segment_names: HashSet<&[u8]> =
        segment_names.into_iter().map(|s| s.as_ref()).collect();

    let segments =
        filtered!(gfa.segments, |s| segment_names.contains(s.name.as_slice()));

    let links = filtered!(&gfa.links, |l| {
        segment_names.contains(l.from_segment.as_slice())
            && segment_names.contains(l.to_segment.as_slice())
    });

    let containments = filtered!(&gfa.containments, |l| {
        segment_names.contains(l.container_name.as_slice())
            && segment_names.contains(l.contained_name.as_slice())
    });

    let paths: Vec<_> = filtered!(gfa.paths, |p| p
        .iter()
        .any(|(s, _)| segment_names.contains(s.as_ref())));

    GFA {
        header: gfa.header.clone(),
        segments,
        links,
        paths,
        containments,
    }
}

/*
macro_rules! filtered_into {
    ($coll:expr, $pred:expr) => {
        $coll.into_iter().filter($pred).collect()
    };
}

// Consume the given GFA to create the subgraph GFA
pub fn paths_subgraph<T: OptFields>(
    gfa: GFA<BString, T>,
    paths: &[BString],
) -> GFA<BString, T> {
    let path_names: HashSet<&[u8]> = paths.iter().collect();
    // Filter out the paths in the GFA we don't want
    let paths: Vec<_> =
        filtered_into!(gfa.paths, |p| path_names.contains(&p.path_name));

    // Set of the segments in the paths we're keeping
    let segment_names: HashSet<&[u8]> = paths
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
        version: gfa.version,
        segments,
        links,
        paths,
        containments,
    }
}

*/