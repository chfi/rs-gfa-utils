use handlegraph::{handle::*, handlegraph::*, hashgraph::HashGraph};

use fnv::{FnvHashMap, FnvHashSet};

use std::collections::{BTreeMap, HashMap, HashSet, VecDeque};

use bstr::{BStr, BString, ByteSlice, ByteVec};

use bio::alphabets::dna;

use gfa::{
    cigar::CIGAR,
    gfa::{Orientation, Path, GFA},
    optfields::OptFields,
};

#[derive(Debug, Default, Clone, PartialEq)]
pub struct SubPath<'a> {
    pub path_name: BString,
    pub steps: Vec<(usize, Orientation, Option<&'a CIGAR>)>,
}

impl<'a> SubPath<'a> {
    pub fn segment_ids(&self) -> impl Iterator<Item = usize> + '_ {
        self.steps.iter().map(|x| x.0)
    }
}

pub fn oriented_sequence<T: AsRef<[u8]>>(
    seq: T,
    orient: Orientation,
) -> BString {
    let seq: &[u8] = seq.as_ref();
    if orient.is_reverse() {
        dna::revcomp(seq).into()
    } else {
        seq.into()
    }
}

pub fn path_segments_sequences<'a, T, I>(
    gfa: &'a GFA<usize, T>,
    subpaths: I,
) -> FnvHashMap<usize, BString>
where
    T: OptFields,
    I: IntoIterator<Item = &'a SubPath<'a>> + 'a,
{
    let all_segments: FnvHashSet<usize> = subpaths
        .into_iter()
        .flat_map(|sub| sub.steps.iter().map(|step| step.0))
        .collect();

    gfa.segments
        .iter()
        .filter(|&seg| all_segments.contains(&seg.name))
        .map(|seg| (seg.name, seg.sequence.clone()))
        .collect()
}

pub fn bubble_sub_paths<T: OptFields>(
    gfa: &GFA<usize, T>,
    from: usize,
    to: usize,
) -> Vec<SubPath<'_>> {
    gfa.paths
        .iter()
        .filter_map(|path| {
            let mut steps = path
                .iter()
                .zip(path.overlaps.iter())
                .skip_while(|&((x, _o), _cg)| x != from && x != to)
                .peekable();

            let &((first, _), _) = steps.peek()?;
            let end = if first == from { to } else { from };

            let steps: Vec<_> = steps
                .scan(first, |previous, ((step, orient), overlap)| {
                    if *previous == end {
                        None
                    } else {
                        *previous = step;
                        Some((step, orient, overlap.as_ref()))
                    }
                })
                .collect();

            Some(SubPath {
                path_name: path.path_name.clone(),
                steps,
            })
        })
        .collect()
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VariantKey {
    pub ref_name: BString,
    pub sequence: BString,
    pub pos: usize,
}

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Variant {
    Del(BString),
    Ins(BString),
    Snv(u8),
    Mnv(BString),
}

impl std::fmt::Display for Variant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Variant::Del(b) => write!(f, "Del({})", b),
            Variant::Ins(b) => write!(f, "Ins({})", b),
            Variant::Snv(b) => write!(f, "Snv({})", char::from(*b)),
            Variant::Mnv(b) => write!(f, "Mnv({})", b),
        }
    }
}

pub fn detect_variants_against_ref(
    segment_sequences: &FnvHashMap<usize, BString>,
    ref_name: &[u8],
    ref_path: &[usize],
    query_path: &[usize],
) -> FnvHashMap<VariantKey, Variant> {
    let mut variants = FnvHashMap::default();

    let mut ref_ix = 0;
    let mut query_ix = 0;

    let mut ref_seq_ix = 0;
    let mut query_seq_ix = 0;

    loop {
        if ref_ix >= ref_path.len() || query_ix >= query_path.len() {
            break;
        }

        let ref_node = ref_path[ref_ix];
        let ref_seq = segment_sequences.get(&ref_node).unwrap();

        let query_node = query_path[query_ix];
        let query_seq = segment_sequences.get(&query_node).unwrap();

        if ref_node == query_node {
            ref_ix += 1;
            ref_seq_ix += ref_seq.len();

            query_ix += 1;
            query_seq_ix += query_seq.len();
        } else {
            if ref_ix + 1 >= ref_path.len() || query_ix + 1 >= query_path.len()
            {
                break;
            }
            let next_ref_node = ref_path[ref_ix + 1];
            let next_query_node = query_path[query_ix + 1];

            if next_ref_node == query_node {
                // Deletion
                let prev_ref_node = ref_path[ref_ix - 1];
                let prev_ref_seq =
                    segment_sequences.get(&prev_ref_node).unwrap();

                let last_prev_seq: u8 = *prev_ref_seq.last().unwrap();

                let key_ref_seq: BString = std::iter::once(last_prev_seq)
                    .chain(ref_seq.iter().copied())
                    .collect();

                let var_key = VariantKey {
                    ref_name: ref_name.into(),
                    pos: ref_seq_ix,
                    sequence: key_ref_seq,
                };

                let variant = Variant::Del(BString::from(&[last_prev_seq][..]));
                variants.insert(var_key, variant);

                ref_ix += 1;
                ref_seq_ix += ref_seq.len();
            } else if next_query_node == ref_node {
                // Insertion
                let prev_ref_node = ref_path[ref_ix - 1];
                let prev_ref_seq =
                    segment_sequences.get(&prev_ref_node).unwrap();

                let last_prev_seq: u8 = *prev_ref_seq.last().unwrap();

                let key_ref_seq: BString =
                    std::iter::once(last_prev_seq).collect();

                let var_key = VariantKey {
                    ref_name: ref_name.into(),
                    pos: ref_seq_ix,
                    sequence: key_ref_seq,
                };

                let var_seq: BString = std::iter::once(last_prev_seq)
                    .chain(query_seq.iter().copied())
                    .collect();
                let variant = Variant::Ins(var_seq);

                variants.insert(var_key, variant);

                query_ix += 1;
                query_seq_ix += query_seq.len();
            } else {
                let var_key = VariantKey {
                    ref_name: ref_name.into(),
                    pos: ref_seq_ix + 1,
                    sequence: ref_seq.clone(),
                };

                let variant = if ref_seq.len() == 1 {
                    let last_query_seq: u8 = *query_seq.last().unwrap();
                    Variant::Snv(last_query_seq)
                } else {
                    Variant::Mnv(query_seq.clone())
                };

                variants.insert(var_key, variant);

                ref_ix += 1;
                ref_seq_ix += ref_seq.len();

                query_ix += 1;
                query_seq_ix += query_seq.len();
            }
        }
    }

    variants
}

pub fn detect_variants_in_sub_paths(
    segment_sequences: &FnvHashMap<usize, BString>,
    sub_paths: &[SubPath<'_>],
) -> FnvHashMap<BString, FnvHashMap<VariantKey, Variant>> {
    let mut variants = FnvHashMap::default();

    for ref_path in sub_paths.iter() {
        let ref_name = ref_path.path_name.clone();
        let ref_steps = ref_path.segment_ids().collect::<Vec<_>>();
        for query in sub_paths.iter() {
            if ref_path.path_name != query.path_name {
                let query_path = query.segment_ids().collect::<Vec<_>>();
                let vars = detect_variants_against_ref(
                    segment_sequences,
                    &ref_name,
                    &ref_steps,
                    &query_path,
                );

                let var_map: &mut FnvHashMap<_, _> =
                    variants.entry(ref_name.clone()).or_default();
                var_map.extend(vars.into_iter());
            }
        }
    }

    variants
}

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

pub fn extract_nodes_in_bubble<T>(
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

/// A struct that holds Variants, as defined in the VCF format
#[derive(Debug, PartialEq)]
pub struct VCFRecord {
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

impl std::fmt::Display for VCFRecord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        fn display_field<T: std::fmt::Display>(field: Option<T>) -> String {
            if let Some(x) = field {
                x.to_string()
            } else {
                ".".to_string()
            }
        }

        write!(f, "{}\t", self.chromosome)?;
        write!(f, "{}\t", self.position)?;
        write!(f, "{}\t", display_field(self.id.as_ref()))?;
        write!(f, "{}\t", self.reference)?;
        write!(f, "{}\t", display_field(self.alternate.as_ref()))?;
        write!(f, "{}\t", display_field(self.quality.as_ref()))?;
        write!(f, "{}\t", display_field(self.filter.as_ref()))?;
        write!(f, "{}\t", display_field(self.info.as_ref()))?;
        write!(f, "{}\t", display_field(self.format.as_ref()))?;
        writeln!(f, "{}", display_field(self.sample_name.as_ref()))
    }
}
