pub mod vcf;

use vcf::VCFRecord;

use bio::alphabets::dna;
use bstr::{BStr, BString, ByteSlice};
use fnv::{FnvHashMap, FnvHashSet};
use gfa::{
    cigar::CIGAR,
    gfa::{Orientation, Path, GFA},
    optfields::OptFields,
};
use handlegraph::{handle::*, handlegraph::*};

use indicatif::{
    ParallelProgressIterator, ProgressBar, ProgressIterator, ProgressStyle,
};
use rayon::prelude::*;

fn progress_bar(len: usize, steady: bool) -> ProgressBar {
    let p_bar = ProgressBar::new(len as u64);
    p_bar.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:80} {pos:>7}/{len:7}")
            .progress_chars("##-"),
    );
    if steady {
        p_bar.enable_steady_tick(1000);
    }
    p_bar
}

#[allow(unused_imports)]
use log::{debug, info, trace, warn};

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
    paths: &Vec<Vec<(usize, usize, Orientation)>>,
    vertices: &FnvHashSet<u64>,
) -> FnvHashMap<u64, FnvHashMap<usize, usize>> {
    let mut transposed: FnvHashMap<usize, FnvHashMap<u64, usize>> =
        FnvHashMap::default();

    {
        debug!("Finding ultrabubble node indices for {} paths", paths.len());
        let p_bar = progress_bar(paths.len(), false);
        transposed.par_extend(
            paths.par_iter().enumerate().progress_with(p_bar).map(
                |(path_ix, path)| {
                    let node_indices: FnvHashMap<u64, usize> = path
                        .iter()
                        .enumerate()
                        .filter_map(|(ix, &(step, _, _))| {
                            let step = step as u64;
                            if vertices.contains(&step) {
                                Some((step, ix))
                            } else {
                                None
                            }
                        })
                        .collect();

                    (path_ix, node_indices)
                },
            ),
        );
    }

    debug!("Transposing path/ultrabubble node index map");
    let p_bar = progress_bar(vertices.len(), true);

    let path_map: FnvHashMap<u64, FnvHashMap<usize, usize>> = vertices
        .par_iter()
        .progress_with(p_bar)
        .map(|&node| {
            let inner = transposed
                .iter()
                .filter_map(|(path_ix, step_map)| {
                    let ix = step_map.get(&node)?;
                    Some((*path_ix, *ix))
                })
                .collect();
            (node, inner)
        })
        .collect();

    path_map
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
    Mnp(BString),
}

impl std::fmt::Display for Variant {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Variant::Del(b) => write!(f, "Del({})", b),
            Variant::Ins(b) => write!(f, "Ins({})", b),
            Variant::Snv(b) => write!(f, "Snv({})", char::from(*b)),
            Variant::Mnp(b) => write!(f, "Mnp({})", b),
        }
    }
}

pub fn detect_variants_against_ref(
    segment_sequences: &FnvHashMap<usize, BString>,
    ref_name: &[u8],
    ref_path: &[(usize, usize, Orientation)],
    query_path: &[(usize, usize, Orientation)],
) -> FnvHashMap<VariantKey, FnvHashSet<Variant>> {
    let mut variants: FnvHashMap<_, FnvHashSet<_>> = FnvHashMap::default();

    let mut ref_ix = 0;
    let mut query_ix = 0;

    let mut ref_seq_ix;

    loop {
        if ref_ix >= ref_path.len() || query_ix >= query_path.len() {
            break;
        }

        let (ref_node, ref_offset, _) = ref_path[ref_ix];
        let ref_seq = segment_sequences.get(&ref_node).unwrap();

        ref_seq_ix = ref_offset;

        let (query_node, _query_offset, _) = query_path[query_ix];
        let query_seq = segment_sequences.get(&query_node).unwrap();

        if ref_node == query_node {
            ref_ix += 1;
            query_ix += 1;
        } else {
            if ref_ix + 1 >= ref_path.len() || query_ix + 1 >= query_path.len()
            {
                trace!("At end of ref or query");
                break;
            }
            let (next_ref_node, _next_ref_offset, _) = ref_path[ref_ix + 1];
            let (next_query_node, _next_query_offset, _) =
                query_path[query_ix + 1];

            if next_ref_node == query_node {
                trace!("Deletion at ref {}\t query {}", ref_ix, query_ix);
                // Deletion
                let (prev_ref_node, _prev_ref_offset, _) = if ref_ix == 0 {
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
                trace!("Insertion at ref {}\t query {}", ref_ix, query_ix);
                // Insertion

                let (prev_ref_node, _prev_ref_offset, _) = if ref_ix == 0 {
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
                        trace!("SNV at ref {}\t query {}", ref_ix, query_ix);
                        let last_query_seq: u8 = *query_seq.last().unwrap();
                        Variant::Snv(last_query_seq)
                    } else {
                        trace!("MNP at ref {}\t query {}", ref_ix, query_ix);
                        Variant::Mnp(query_seq.as_bstr().to_owned())
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct VariantConfig {
    pub ignore_inverted_paths: bool,
}

impl VariantConfig {
    pub fn ignore_path(
        &self,
        ref_orient: (Orientation, Orientation),
        query_orient: (Orientation, Orientation),
    ) -> bool {
        if self.ignore_inverted_paths && ref_orient != query_orient {
            trace!("Ignoring inverted path");
            true
        } else {
            false
        }
    }
}

impl Default for VariantConfig {
    fn default() -> Self {
        Self {
            ignore_inverted_paths: true,
        }
    }
}

pub fn detect_variants_in_sub_paths(
    variant_config: &VariantConfig,
    segment_sequences: &FnvHashMap<usize, BString>,
    path_names: &[BString],
    paths: &[Vec<(usize, usize, Orientation)>],
    ref_path_names: Option<&FnvHashSet<BString>>,
    path_indices: &FnvHashMap<u64, FnvHashMap<usize, usize>>,
    from: u64,
    to: u64,
) -> Option<FnvHashMap<BString, FnvHashMap<VariantKey, FnvHashSet<Variant>>>> {
    let from_indices = path_indices.get(&from)?;
    let to_indices = path_indices.get(&to)?;

    let mut variants: FnvHashMap<BString, FnvHashMap<_, FnvHashSet<_>>> =
        FnvHashMap::default();

    let sub_paths: Vec<(usize, &[(usize, usize, Orientation)])> = paths
        .iter()
        .enumerate()
        .filter_map(|(path_ix, path)| {
            let from_ix = *from_indices.get(&path_ix)?;
            let to_ix = *to_indices.get(&path_ix)?;
            let from = from_ix.min(to_ix);
            let to = from_ix.max(to_ix);
            let sub_path = &path[from..=to];
            if sub_path.len() > 2 {
                Some((path_ix, sub_path))
            } else {
                None
            }
        })
        .collect();

    let mut query_paths = sub_paths.clone();

    query_paths.sort_by(|(_, v), (_, w)| {
        let v_iter = v.iter().map(|(a, b, c)| (a, b, bool::from(*c)));
        let w_iter = w.iter().map(|(a, b, c)| (a, b, bool::from(*c)));
        v_iter.cmp(w_iter)
    });

    query_paths.dedup_by(|(_, v), (_, w)| v == w);

    let is_ref_path = |p: &BStr| {
        if let Some(ref_path_names) = ref_path_names {
            ref_path_names.contains(p)
        } else {
            true
        }
    };

    variants.extend(sub_paths.iter().filter_map(|(ref_ix, ref_path)| {
        let ref_name = path_names.get(*ref_ix)?;
        if !is_ref_path(ref_name.as_ref()) {
            return None;
        }
        let ref_orient = sub_path_edge_orient(ref_path);

        let mut ref_map: FnvHashMap<VariantKey, FnvHashSet<_>> =
            FnvHashMap::default();

        for (query_ix, query_path) in query_paths.iter() {
            let query_name = path_names.get(*query_ix)?;
            let query_orient = sub_path_edge_orient(query_path);

            if ref_name != query_name
                && !variant_config.ignore_path(ref_orient, query_orient)
            {
                let vars = detect_variants_against_ref(
                    segment_sequences,
                    ref_name,
                    ref_path,
                    query_path,
                );

                ref_map.extend(vars)
            }
        }

        let ref_name: BString = ref_name.clone().to_owned();
        Some((ref_name, ref_map))
    }));

    Some(variants)
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
                    Variant::Mnp(seq) => (seq.clone(), "mnp".into()),
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
                format: None,
                sample_name: None,
            };

            vcf_records.push(vcf);
        }
    }

    vcf_records
}
