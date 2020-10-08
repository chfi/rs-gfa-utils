use handlegraph::{handle::*, handlegraph::*};

use fnv::{FnvHashMap, FnvHashSet};

use bstr::{BStr, BString, ByteSlice};

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

pub fn bubble_path_indices(
    paths: &FnvHashMap<BString, Vec<(usize, usize, Orientation)>>,
    // paths: &FnvHashMap<BString, Vec<(usize, usize)>>,
    vertices: &FnvHashSet<u64>,
) -> FnvHashMap<u64, FnvHashMap<BString, usize>> {
    let mut path_map: FnvHashMap<u64, FnvHashMap<BString, usize>> =
        FnvHashMap::default();

    for &node in vertices.iter() {
        for (path_name, path) in paths.iter() {
            let node_ix = path.iter().position(|&(x, _, _)| x == node as usize);
            if let Some(ix) = node_ix {
                let path_name = path_name.clone().to_owned();
                let entry = path_map.entry(node).or_default();
                entry.insert(path_name, ix);
            }
        }
    }

    path_map
}

pub fn gfa_paths_with_offsets(
    gfa: &GFA<usize, ()>,
    seg_seq_map: &FnvHashMap<usize, &[u8]>,
) -> FnvHashMap<BString, Vec<(usize, usize, Orientation)>> {
    gfa.paths
        .iter()
        .map(|path| {
            let name = path.path_name.clone();
            let steps: Vec<(usize, usize, Orientation)> = path
                .iter()
                .scan(1, |offset, (step, orient)| {
                    let step_offset = *offset;
                    let step_len = seg_seq_map.get(&step).unwrap().len();
                    *offset += step_len;
                    Some((step, step_offset, orient))
                })
                .collect();

            (name, steps)
        })
        .collect()
}

pub fn gfa_paths<T: OptFields>(
    gfa: &GFA<usize, T>,
) -> FnvHashMap<BString, Vec<(usize, Orientation)>> {
    gfa.paths
        .iter()
        .map(|path| {
            let name = path.path_name.clone();
            let steps = path.iter().collect::<Vec<_>>();
            (name, steps)
        })
        .collect()
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

pub fn sub_path_base_offset<T: OptFields>(
    segment_map: &FnvHashMap<usize, &[u8]>,
    path: &Path<usize, T>,
    step_ix: usize,
) -> Option<usize> {
    let result =
        path.iter()
            .enumerate()
            .try_fold(0, |offset, (ix, (seg, _))| {
                if ix == step_ix {
                    Err(Some(offset))
                } else {
                    let seq = segment_map.get(&seg).ok_or(None)?;
                    Ok(offset + seq.len())
                }
            });
    result.err().unwrap()
}

pub fn many_sub_paths<T>(
    gfa: &GFA<usize, T>,
    indices: &FnvHashSet<u64>,
) -> FnvHashMap<BString, FnvHashMap<u64, Vec<usize>>>
where
    T: OptFields,
{
    gfa.paths
        .iter()
        .map(|path| {
            let mut index_map: FnvHashMap<u64, Vec<usize>> =
                FnvHashMap::default();

            for (ix, (step, _orient)) in path.iter().enumerate() {
                let step_u64 = step as u64;
                if indices.contains(&step_u64) {
                    let entry = index_map.entry(step_u64).or_default();
                    entry.push(ix);
                }
            }

            (path.path_name.clone(), index_map)
        })
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

            let mut found_end = false;

            let steps: Vec<_> = steps
                .scan(first, |previous, ((step, orient), overlap)| {
                    if *previous == end {
                        found_end = true;
                        None
                    } else {
                        *previous = step;
                        Some((step, orient, overlap.as_ref()))
                    }
                })
                .collect();

            if !found_end {
                return None;
            }

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
    segment_sequences: &FnvHashMap<usize, &[u8]>,
    ref_name: &[u8],
    ref_path: &[(usize, usize, Orientation)],
    query_path: &[(usize, usize, Orientation)],
) -> FnvHashMap<VariantKey, FnvHashSet<Variant>> {
    let mut variants: FnvHashMap<_, FnvHashSet<_>> = FnvHashMap::default();

    let mut ref_ix = 0;
    let mut query_ix = 0;

    let mut ref_seq_ix = 0;

    loop {
        if ref_ix >= ref_path.len() || query_ix >= query_path.len() {
            break;
        }

        let (ref_node, ref_offset, _) = ref_path[ref_ix];
        let ref_seq = segment_sequences.get(&ref_node).unwrap();

        ref_seq_ix = ref_offset;

        let (query_node, query_offset, _) = query_path[query_ix];
        let query_seq = segment_sequences.get(&query_node).unwrap();

        if ref_node == query_node {
            ref_ix += 1;
            query_ix += 1;
        } else {
            if ref_ix + 1 >= ref_path.len() || query_ix + 1 >= query_path.len()
            {
                break;
            }
            let (next_ref_node, next_ref_offset, _) = ref_path[ref_ix + 1];
            let (next_query_node, next_query_offset, _) =
                query_path[query_ix + 1];

            if next_ref_node == query_node {
                // Deletion
                let (prev_ref_node, prev_ref_offset, _) = if ref_ix == 0 {
                    ref_path[ref_ix]
                } else {
                    ref_path[ref_ix - 1]
                };

                let prev_ref_seq =
                    segment_sequences.get(&prev_ref_node).unwrap();

                let last_prev_seq: u8 = *prev_ref_seq.last().unwrap();

                let key_ref_seq: BString = std::iter::once(last_prev_seq)
                    .chain(ref_seq.iter().copied())
                    .collect();

                let var_key = VariantKey {
                    ref_name: ref_name.into(),
                    pos: ref_seq_ix - 1,
                    sequence: key_ref_seq,
                };

                let variant = Variant::Del(BString::from(&[last_prev_seq][..]));

                let entry = variants.entry(var_key).or_default();
                entry.insert(variant);

                ref_ix += 1;
            } else if next_query_node == ref_node {
                // Insertion

                let (prev_ref_node, prev_ref_offset, _) = if ref_ix == 0 {
                    ref_path[ref_ix]
                } else {
                    ref_path[ref_ix - 1]
                };
                let prev_ref_seq =
                    segment_sequences.get(&prev_ref_node).unwrap();

                let last_prev_seq: u8 = *prev_ref_seq.last().unwrap();

                let key_ref_seq: BString =
                    std::iter::once(last_prev_seq).collect();

                let var_key = VariantKey {
                    ref_name: ref_name.into(),
                    pos: ref_seq_ix - 1,
                    sequence: key_ref_seq,
                };

                let var_seq: BString = std::iter::once(last_prev_seq)
                    .chain(query_seq.iter().copied())
                    .collect();
                let variant = Variant::Ins(var_seq);

                let entry = variants.entry(var_key).or_default();
                entry.insert(variant);

                query_ix += 1;
            } else {
                if ref_seq != query_seq {
                    let var_key = VariantKey {
                        ref_name: ref_name.into(),
                        pos: ref_seq_ix,
                        sequence: ref_seq.as_bstr().to_owned(),
                    };

                    let variant = if ref_seq.len() == 1 {
                        let last_query_seq: u8 = *query_seq.last().unwrap();
                        Variant::Snv(last_query_seq)
                    } else {
                        Variant::Mnv(query_seq.as_bstr().to_owned())
                    };

                    let entry = variants.entry(var_key).or_default();
                    entry.insert(variant);
                }

                ref_ix += 1;

                query_ix += 1;
            }
        }
    }

    variants
}

fn sub_path_edge_orient(
    path: &[(usize, usize, Orientation)],
) -> (Orientation, Orientation) {
    let from = path.first().unwrap().2;
    let to = path.last().unwrap().2;
    (from, to)
}

pub fn detect_variants_in_sub_paths(
    segment_sequences: &FnvHashMap<usize, &[u8]>,
    paths: &FnvHashMap<BString, Vec<(usize, usize, Orientation)>>,
    path_indices: &FnvHashMap<u64, FnvHashMap<BString, usize>>,
    from: u64,
    to: u64,
) -> FnvHashMap<BString, FnvHashMap<VariantKey, FnvHashSet<Variant>>> {
    let mut variants: FnvHashMap<BString, FnvHashMap<_, FnvHashSet<_>>> =
        FnvHashMap::default();

    let from_indices = path_indices.get(&from).unwrap();
    let to_indices = path_indices.get(&to).unwrap();

    let sub_paths: FnvHashMap<&BStr, &[(usize, usize, Orientation)]> = paths
        .iter()
        .filter_map(|(path_name, path)| {
            let from_ix = *from_indices.get(path_name)?;
            let to_ix = *to_indices.get(path_name)?;
            let from = from_ix.min(to_ix);
            let to = from_ix.max(to_ix);
            let sub_path = &path[from..=to];
            Some((path_name.as_bstr(), sub_path))
        })
        .collect();

    for (ref_name, ref_path) in sub_paths.iter() {
        let ref_orient = sub_path_edge_orient(ref_path);
        for (query_name, query_path) in sub_paths.iter() {
            let query_orient = sub_path_edge_orient(query_path);

            if ref_name != query_name && ref_orient == query_orient {
                let vars = detect_variants_against_ref(
                    segment_sequences,
                    ref_name,
                    ref_path,
                    query_path,
                );

                let ref_name: BString = ref_name.clone().to_owned();
                let var_map: &mut FnvHashMap<_, _> =
                    variants.entry(ref_name).or_default();
                var_map.extend(vars.into_iter());
            }
        }
    }

    variants
}

pub fn variant_vcf_record(
    variants: &FnvHashMap<BString, FnvHashMap<VariantKey, FnvHashSet<Variant>>>,
) -> Vec<VCFRecord> {
    let mut vcf_records = Vec::new();

    for (_, variant_map) in variants.iter() {
        for (key, var_set) in variant_map.iter() {
            let (alt_list, type_set): (Vec<BString>, Vec<BString>) = var_set
                .iter()
                .map(|var| match var {
                    Variant::Del(seq) => (seq.clone(), "del".into()),
                    Variant::Ins(seq) => (seq.clone(), "ins".into()),
                    Variant::Snv(base) => {
                        let base_seq =
                            std::iter::once(*base).collect::<BString>();
                        (base_seq, "snv".into())
                    }
                    Variant::Mnv(seq) => (seq.clone(), "mnv".into()),
                })
                .unzip();

            let alts = bstr::join(",", alt_list);
            let mut types: BString = "TYPE=".into();
            let types_temp = bstr::join(";TYPE=", type_set);
            types.extend(types_temp);

            let vcf = VCFRecord {
                chromosome: key.ref_name.clone(),
                position: key.pos as i64,
                id: None,
                reference: key.sequence.clone(),
                alternate: Some(alts.into()),
                quality: None,
                filter: None,
                info: Some(types),
                format: Some("GT".into()),
                sample_name: Some("0|1".into()),
            };

            vcf_records.push(vcf);
        }
    }

    vcf_records
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
    position: i64,
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
        write!(f, "{}", display_field(self.sample_name.as_ref()))
    }
}
