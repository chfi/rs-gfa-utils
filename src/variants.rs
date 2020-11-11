pub mod vcf;

use vcf::VCFRecord;

use bio::alphabets::dna;
use bstr::{BStr, BString, ByteSlice};
use fnv::{FnvHashMap, FnvHashSet};
use indicatif::ParallelProgressIterator;
use rayon::prelude::*;

use gfa::gfa::{Orientation, GFA};

use crate::util::progress_bar;

#[allow(unused_imports)]
use log::{debug, info, trace, warn};

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

pub type PathStep = (usize, usize, Orientation);

pub struct PathData {
    pub segment_map: FnvHashMap<usize, BString>,
    pub path_names: Vec<BString>,
    pub paths: Vec<Vec<PathStep>>,
}

pub fn gfa_path_data(mut gfa: GFA<usize, ()>) -> PathData {
    let segments = std::mem::take(&mut gfa.segments);

    info!("Building map from segment IDs to sequences");
    let segment_map: FnvHashMap<usize, BString> = segments
        .into_iter()
        .map(|seg| (seg.name, seg.sequence))
        .collect();

    let gfa_paths = std::mem::take(&mut gfa.paths);

    let p_bar = progress_bar(gfa_paths.len(), false);

    info!("Extracting paths and offsets from GFA");
    let (path_names, paths): (Vec<_>, Vec<_>) = gfa_paths
        .into_par_iter()
        .progress_with(p_bar)
        .map(|mut path| {
            let steps: Vec<(usize, usize, Orientation)> = path
                .iter()
                .scan(1, |offset, (step, orient)| {
                    let step_offset = *offset;
                    let step_len = segment_map.get(&step).unwrap().len();
                    *offset += step_len;
                    Some((step, step_offset, orient))
                })
                .collect();

            let path_name = std::mem::take(&mut path.path_name);

            (path_name, steps)
        })
        .unzip();

    PathData {
        segment_map,
        path_names,
        paths,
    }
}

pub fn bubble_path_indices(
    paths: &[Vec<(usize, usize, Orientation)>],
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

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct VariantKey {
    pub ref_name: BString,
    pub sequence: BString,
    pub pos: usize,
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
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

/// Abstraction to handle the different cases in
/// `detect_variants_against_ref_with`
trait VariantHandler {
    fn deletion(
        &mut self,
        ref_ix: usize,
        query_ix: usize,
        ref_seq_ix: usize,
        query_seq_ix: usize,
    );

    fn insertion(
        &mut self,
        ref_ix: usize,
        query_ix: usize,
        ref_seq_ix: usize,
        query_seq_ix: usize,
    );

    fn mismatch(
        &mut self,
        ref_ix: usize,
        query_ix: usize,
        ref_seq_ix: usize,
        query_seq_ix: usize,
    );

    fn match_(
        &mut self,
        ref_ix: usize,
        query_ix: usize,
        ref_seq_ix: usize,
        query_seq_ix: usize,
    );
}

fn detect_variants_against_ref_with<H: VariantHandler>(
    segment_sequences: &FnvHashMap<usize, BString>,
    ref_path: &[(usize, usize, Orientation)],
    query_path: &[(usize, usize, Orientation)],
    handler: &mut H,
) {
    let mut ref_ix = 0;
    let mut query_ix = 0;

    let mut ref_seq_ix;
    let mut query_seq_ix;

    loop {
        if ref_ix >= ref_path.len() || query_ix >= query_path.len() {
            break;
        }

        let (ref_node, ref_offset, _) = ref_path[ref_ix];
        let ref_seq = segment_sequences.get(&ref_node).unwrap();

        ref_seq_ix = ref_offset;

        let (query_node, query_offset, _) = query_path[query_ix];
        let query_seq = segment_sequences.get(&query_node).unwrap();

        query_seq_ix = query_offset;

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
                handler.deletion(ref_ix, query_ix, ref_seq_ix, query_seq_ix);

                ref_ix += 1;
            } else if next_query_node == ref_node {
                trace!("Insertion at ref {}\t query {}", ref_ix, query_ix);
                // Insertion
                handler.insertion(ref_ix, query_ix, ref_seq_ix, query_seq_ix);

                query_ix += 1;
            } else {
                if ref_seq != query_seq {
                    handler.mismatch(
                        ref_ix,
                        query_ix,
                        ref_seq_ix,
                        query_seq_ix,
                    );
                } else {
                    handler.match_(ref_ix, query_ix, ref_seq_ix, query_seq_ix);
                }

                ref_ix += 1;
                query_ix += 1;
            }
        }
    }
}

/// Implementation of `VariantHandler` that fills a hashmap of
/// variants, same as the original `detect_variants_against_ref`
#[derive(Debug, Clone)]
struct VCFVariantHandler<'a> {
    segment_sequences: &'a FnvHashMap<usize, BString>,
    ref_name: &'a [u8],
    ref_path: &'a [(usize, usize, Orientation)],
    query_path: &'a [(usize, usize, Orientation)],
    variants: FnvHashMap<VariantKey, FnvHashSet<Variant>>,
}

impl<'a> VCFVariantHandler<'a> {
    fn new(
        segment_sequences: &'a FnvHashMap<usize, BString>,
        ref_name: &'a [u8],
        ref_path: &'a [(usize, usize, Orientation)],
        query_path: &'a [(usize, usize, Orientation)],
    ) -> Self {
        Self {
            segment_sequences,
            ref_name,
            ref_path,
            query_path,
            variants: FnvHashMap::default(),
        }
    }
}

impl<'a> VariantHandler for VCFVariantHandler<'a> {
    fn deletion(
        &mut self,
        ref_ix: usize,
        _query_ix: usize,
        ref_seq_ix: usize,
        _query_seq_ix: usize,
    ) {
        let (ref_node, _ref_offset, _) = self.ref_path[ref_ix];
        let ref_seq = self.segment_sequences.get(&ref_node).unwrap();

        // Deletion
        let (prev_ref_node, _prev_ref_offset, _) = if ref_ix == 0 {
            self.ref_path[ref_ix]
        } else {
            self.ref_path[ref_ix - 1]
        };

        let prev_ref_seq = self.segment_sequences.get(&prev_ref_node).unwrap();

        let last_prev_seq: u8 = *prev_ref_seq.last().unwrap();

        let key_ref_seq: BString = std::iter::once(last_prev_seq)
            .chain(ref_seq.iter().copied())
            .collect();

        let var_key = VariantKey {
            ref_name: self.ref_name.into(),
            pos: ref_seq_ix - 1,
            sequence: key_ref_seq,
        };

        let variant = Variant::Del(BString::from(&[last_prev_seq][..]));

        let entry = self.variants.entry(var_key).or_default();
        entry.insert(variant);
    }

    fn insertion(
        &mut self,
        ref_ix: usize,
        query_ix: usize,
        ref_seq_ix: usize,
        _query_seq_ix: usize,
    ) {
        let (query_node, _query_offset, _) = self.query_path[query_ix];
        let query_seq = self.segment_sequences.get(&query_node).unwrap();

        let (prev_ref_node, _prev_ref_offset, _) = if ref_ix == 0 {
            self.ref_path[ref_ix]
        } else {
            self.ref_path[ref_ix - 1]
        };
        let prev_ref_seq = self.segment_sequences.get(&prev_ref_node).unwrap();

        let last_prev_seq: u8 = *prev_ref_seq.last().unwrap();

        let key_ref_seq: BString = std::iter::once(last_prev_seq).collect();

        let var_key = VariantKey {
            ref_name: self.ref_name.into(),
            pos: ref_seq_ix - 1,
            sequence: key_ref_seq,
        };

        let var_seq: BString = std::iter::once(last_prev_seq)
            .chain(query_seq.iter().copied())
            .collect();
        let variant = Variant::Ins(var_seq);

        let entry = self.variants.entry(var_key).or_default();
        entry.insert(variant);
    }

    fn mismatch(
        &mut self,
        ref_ix: usize,
        query_ix: usize,
        ref_seq_ix: usize,
        _query_seq_ix: usize,
    ) {
        let (ref_node, _ref_offset, _) = self.ref_path[ref_ix];
        let ref_seq = self.segment_sequences.get(&ref_node).unwrap();

        let (query_node, _query_offset, _) = self.query_path[query_ix];
        let query_seq = self.segment_sequences.get(&query_node).unwrap();

        let var_key = VariantKey {
            ref_name: self.ref_name.into(),
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

        let entry = self.variants.entry(var_key).or_default();
        entry.insert(variant);
    }

    fn match_(
        &mut self,
        _ref_ix: usize,
        _query_ix: usize,
        _ref_seq_ix: usize,
        _query_seq_ix: usize,
    ) {
    }
}

/// NB Should be equivalent to `detect_variants_against_ref` in
/// function, but I need to test it, and I'm not sure about how
/// performance may differ
pub fn detect_variants_against_ref_(
    segment_sequences: &FnvHashMap<usize, BString>,
    ref_name: &[u8],
    ref_path: &[(usize, usize, Orientation)],
    query_path: &[(usize, usize, Orientation)],
) -> FnvHashMap<VariantKey, FnvHashSet<Variant>> {
    let mut handler = VCFVariantHandler::new(
        segment_sequences,
        ref_name,
        ref_path,
        query_path,
    );

    detect_variants_against_ref_with(
        segment_sequences,
        ref_path,
        query_path,
        &mut handler,
    );

    handler.variants
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct SNPRow {
    pub ref_pos: usize,
    pub query_pos: usize,
    pub ref_base: u8,
    pub query_base: u8,
}

#[derive(Debug, Clone)]
struct SNPVariantHandler<'a> {
    segment_sequences: &'a FnvHashMap<usize, BString>,
    ref_path: &'a [(usize, usize, Orientation)],
    query_path: &'a [(usize, usize, Orientation)],
    snp_rows: Vec<SNPRow>,
}

impl<'a> SNPVariantHandler<'a> {
    fn new(
        segment_sequences: &'a FnvHashMap<usize, BString>,
        ref_path: &'a [(usize, usize, Orientation)],
        query_path: &'a [(usize, usize, Orientation)],
    ) -> Self {
        Self {
            segment_sequences,
            ref_path,
            query_path,
            snp_rows: Vec::new(),
        }
    }
}

impl<'a> VariantHandler for SNPVariantHandler<'a> {
    fn deletion(&mut self, _: usize, _: usize, _: usize, _: usize) {}
    fn insertion(&mut self, _: usize, _: usize, _: usize, _: usize) {}

    fn mismatch(
        &mut self,
        ref_ix: usize,
        query_ix: usize,
        ref_seq_ix: usize,
        query_seq_ix: usize,
    ) {
        let (ref_node, _ref_offset, _) = self.ref_path[ref_ix];
        let ref_seq = self.segment_sequences.get(&ref_node).unwrap();

        let (query_node, _query_offset, _) = self.query_path[query_ix];
        let query_seq = self.segment_sequences.get(&query_node).unwrap();

        if ref_seq.len() == 1 && query_seq.len() == 1 {
            let ref_base = ref_seq[0];
            let query_base = query_seq[0];

            let snp_row = SNPRow {
                ref_pos: ref_seq_ix,
                query_pos: query_seq_ix,
                ref_base,
                query_base,
            };
            self.snp_rows.push(snp_row);
        } else {
            debug!("TODO: SNPVariantHandler ignoring mismatch with ref and/or query nodes not being length 1");

            /*
            let ref_base = ref_seq[0];
            let query_base = query_seq[0];

            let snp_row = SNPRow {
                ref_pos: ref_seq_ix,
                query_pos: query_seq_ix,
                ref_base,
                query_base,
            };
            self.snp_rows.push(snp_row);
            */
        }
    }

    fn match_(&mut self, _: usize, _: usize, _: usize, _: usize) {}
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

pub type PathIndices = FnvHashMap<u64, FnvHashMap<usize, usize>>;

fn path_data_sub_paths<'a, 'b>(
    path_data: &'a PathData,
    path_indices: &'b PathIndices,
    from: u64,
    to: u64,
) -> Option<Vec<(usize, &'a [PathStep])>> {
    let from_indices = path_indices.get(&from)?;
    let to_indices = path_indices.get(&to)?;

    let sub_paths = path_data
        .paths
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

    Some(sub_paths)
}

pub fn detect_variants_in_sub_paths(
    variant_config: &VariantConfig,
    path_data: &PathData,
    ref_path_names: Option<&FnvHashSet<BString>>,
    path_indices: &FnvHashMap<u64, FnvHashMap<usize, usize>>,
    from: u64,
    to: u64,
) -> Option<FnvHashMap<BString, FnvHashMap<VariantKey, FnvHashSet<Variant>>>> {
    let mut variants: FnvHashMap<BString, FnvHashMap<_, FnvHashSet<_>>> =
        FnvHashMap::default();

    let sub_paths = path_data_sub_paths(path_data, path_indices, from, to)?;

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
        let ref_name = path_data.path_names.get(*ref_ix)?;
        if !is_ref_path(ref_name.as_ref()) {
            return None;
        }
        let ref_orient = sub_path_edge_orient(ref_path);

        let mut ref_map: FnvHashMap<VariantKey, FnvHashSet<_>> =
            FnvHashMap::default();

        for (query_ix, query_path) in query_paths.iter() {
            let query_name = path_data.path_names.get(*query_ix)?;
            let query_orient = sub_path_edge_orient(query_path);

            if ref_name != query_name
                && !variant_config.ignore_path(ref_orient, query_orient)
            {
                let vars = detect_variants_against_ref(
                    &path_data.segment_map,
                    ref_name,
                    ref_path,
                    query_path,
                );

                ref_map.extend(vars)
            }
        }

        let ref_name: BString = ref_name.clone();
        Some((ref_name, ref_map))
    }));

    Some(variants)
}

pub fn find_snps_in_sub_paths(
    path_data: &PathData,
    ref_path_ix: usize,
    path_indices: &PathIndices,
    from: u64,
    to: u64,
) -> Option<FnvHashMap<BString, Vec<SNPRow>>> {
    let mut query_snp_map: FnvHashMap<BString, Vec<SNPRow>> =
        FnvHashMap::default();

    let sub_paths = path_data_sub_paths(path_data, path_indices, from, to)?;

    let ref_sub_path = sub_paths.iter().find(|&(ix, _)| ix == &ref_path_ix)?;
    let ref_sub_path = ref_sub_path.1;

    for (path_ix, query_path) in sub_paths.iter() {
        if let Some(query_name) = path_data.path_names.get(*path_ix) {
            let mut snp_handler = SNPVariantHandler::new(
                &path_data.segment_map,
                ref_sub_path,
                query_path,
            );

            detect_variants_against_ref_with(
                &path_data.segment_map,
                ref_sub_path,
                query_path,
                &mut snp_handler,
            );

            let snp_rows = snp_handler.snp_rows;

            let query_name: BString = query_name.clone();
            let entry = query_snp_map.entry(query_name).or_default();
            entry.extend(snp_rows);
        }
    }

    Some(query_snp_map)
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
