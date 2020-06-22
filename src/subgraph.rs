use gfa::gfa::GFA;

use std::collections::HashSet;

// Build a new GFA consisting of subgraphs of the given GFA
pub fn paths_new_subgraph(gfa: &GFA, paths: &[String]) -> GFA {
    let path_names: HashSet<&String> = paths.iter().collect();
    // Filter out the paths in the GFA we don't want
    let paths: Vec<_> = gfa
        .paths
        .iter()
        .filter(|p| path_names.contains(&p.path_name))
        .cloned()
        .collect();

    // Set of the segments in the paths we're keeping
    let mut segment_names: HashSet<&str> = HashSet::new();

    paths.iter().for_each(|path| {
        path.segment_names.iter().for_each(|(seg, _)| {
            segment_names.insert(seg);
        });
    });

    // Filter out the segments in the GFA we don't want
    let segments: Vec<_> = gfa
        .segments
        .iter()
        .filter(|s| segment_names.contains(s.name.as_str()))
        .cloned()
        .collect();

    // Filter out the links in the GFA we don't want
    let links: Vec<_> = gfa
        .links
        .iter()
        .filter(|l| {
            segment_names.contains(l.from_segment.as_str())
                && segment_names.contains(l.to_segment.as_str())
        })
        .cloned()
        .collect();

    let containments: Vec<_> = gfa
        .containments
        .iter()
        .filter(|l| {
            segment_names.contains(l.container_name.as_str())
                && segment_names.contains(l.contained_name.as_str())
        })
        .cloned()
        .collect();

    GFA {
        segments,
        links,
        paths,
        containments,
    }
}

// Consume the given GFA to create the subgraph GFA
pub fn paths_subgraph(gfa: GFA, paths: &[String]) -> GFA {
    let path_names: HashSet<&String> = paths.iter().collect();
    // Filter out the paths in the GFA we don't want
    let paths: Vec<_> = gfa
        .paths
        .into_iter()
        .filter(|p| path_names.contains(&p.path_name))
        .collect();

    // Set of the segments in the paths we're keeping
    let mut segment_names: HashSet<&str> = HashSet::new();

    paths.iter().for_each(|path| {
        path.segment_names.iter().for_each(|(seg, _)| {
            segment_names.insert(seg);
        });
    });

    // Filter out the segments in the GFA we don't want
    let segments: Vec<_> = gfa
        .segments
        .into_iter()
        .filter(|s| segment_names.contains(s.name.as_str()))
        .collect();

    // Filter out the links in the GFA we don't want
    let links: Vec<_> = gfa
        .links
        .into_iter()
        .filter(|l| {
            segment_names.contains(l.from_segment.as_str())
                && segment_names.contains(l.to_segment.as_str())
        })
        .collect();

    let containments: Vec<_> = gfa
        .containments
        .into_iter()
        .filter(|l| {
            segment_names.contains(l.container_name.as_str())
                && segment_names.contains(l.contained_name.as_str())
        })
        .collect();

    GFA {
        segments,
        links,
        paths,
        containments,
    }
}
