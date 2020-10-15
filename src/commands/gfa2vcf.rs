use bstr::{io::*, BStr, BString, ByteSlice, ByteVec};
use fnv::{FnvHashMap, FnvHashSet};
use std::{
    fs::File,
    io::Write,
    io::{BufReader, Read},
    path::PathBuf,
};
use structopt::StructOpt;

use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use rayon::prelude::*;

use gfa::gfa::{Orientation, GFA};

#[allow(unused_imports)]
use log::{debug, info, log_enabled, warn};

use crate::variants;

use super::{load_gfa, Result};

/// Output a VCF for the given GFA, using the graph's ultrabubbles to
/// identify areas of variation. (experimental!)
#[derive(StructOpt, Debug)]
pub struct GFA2VCFArgs {
    /// Load ultrabubbles from a file instead of calculating them.
    #[structopt(
        name = "ultrabubbles file",
        long = "ultrabubbles",
        short = "ub"
    )]
    ultrabubbles_file: Option<PathBuf>,
    /// Don't compare two paths if their start and end orientations
    /// don't match each other
    #[structopt(name = "ignore inverted paths", long = "no-inv")]
    ignore_inverted_paths: bool,
    #[structopt(
        name = "file containing paths to use as references",
        long = "paths-file"
    )]
    ref_paths_file: Option<PathBuf>,
    #[structopt(name = "list of paths to use as references", long = "refs")]
    ref_paths_vec: Option<Vec<String>>,
}

fn load_paths_file(file_path: PathBuf) -> Result<Vec<BString>> {
    let file = File::open(file_path)?;
    let reader = BufReader::new(file);
    let lines = reader.byte_lines();

    let mut paths = Vec::new();
    for line in lines {
        let line = line?;
        paths.push(line.into());
    }

    Ok(paths)
}

fn paths_list(paths: Vec<String>) -> Vec<BString> {
    paths.into_iter().map(BString::from).collect()
}

pub fn gfa2vcf(gfa_path: &PathBuf, args: GFA2VCFArgs) -> Result<()> {
    let ref_paths_list = args.ref_paths_vec.map(paths_list).unwrap_or_default();

    let ref_paths_file = args
        .ref_paths_file
        .map(load_paths_file)
        .transpose()?
        .unwrap_or_default();

    let ref_path_names: Option<FnvHashSet<BString>> = {
        let ref_paths: FnvHashSet<BString> = ref_paths_list
            .into_iter()
            .chain(ref_paths_file.into_iter())
            .collect();
        if ref_paths.is_empty() {
            None
        } else {
            if log_enabled!(log::Level::Debug) {
                debug!("Using reference paths:");
                for p in ref_paths.iter() {
                    debug!("\t{}", p);
                }
            }
            Some(ref_paths)
        }
    };

    let gfa: GFA<usize, ()> = load_gfa(&gfa_path)?;

    if gfa.paths.len() < 2 {
        panic!("GFA must contain at least two paths");
    }

    info!("GFA has {} paths", gfa.paths.len());

    info!("Building map from segment IDs to sequences");
    let segment_map: FnvHashMap<usize, &[u8]> = gfa
        .segments
        .iter()
        .map(|seg| (seg.name, seg.sequence.as_ref()))
        .collect();

    info!("Extracting paths and offsets from GFA");
    let all_paths = variants::gfa_paths_with_offsets(&gfa, &segment_map);

    info!("Finding graph ultrabubbles");
    let mut ultrabubbles = if let Some(path) = &args.ultrabubbles_file {
        super::saboten::load_ultrabubbles(path)
    } else {
        super::saboten::find_ultrabubbles(gfa_path)
    }?;

    info!("Using {} ultrabubbles", ultrabubbles.len());

    ultrabubbles.sort();

    let mut representative_paths = Vec::new();

    let mut remaining_ultrabubbles: FnvHashMap<u64, u64> =
        ultrabubbles.iter().copied().collect();


    for (count, (path, steps)) in all_paths.iter().enumerate() {
        if remaining_ultrabubbles.is_empty() {
            break;
        }
        let path_name = path.clone().to_owned();

        let maybe_contained: Vec<(u64, u64)> = steps
            .iter()
            .filter_map(|&(step, _, _)| {
                let x = step as u64;
                let y = remaining_ultrabubbles.get(&x)?;
                Some((*y, x))
            })
            .collect();

        if maybe_contained.is_empty() {
            continue;
        }

        let contained = steps
            .iter()
            .flat_map(|&(step, _, _)| {
                let y = step as u64;
                maybe_contained.iter().filter_map(
                    move |&(a, b)| {
                        if a == y {
                            Some((b, a))
                        } else {
                            None
                        }
                    },
                )
            })
            .collect::<Vec<_>>();

        for &(x, y) in contained.iter() {
            remaining_ultrabubbles.remove(&x);
        }

        let bubbles = contained.into_iter().collect::<Vec<_>>();

        if !bubbles.is_empty() {
            representative_paths.push((path_name, bubbles));
        }
    }

    info!("Chose {} reference paths", representative_paths.len());
    info!(
        "{} ultrabubbles did not exist in any paths",
        remaining_ultrabubbles.len()
    );

    let ultrabubble_nodes = ultrabubbles
        .iter()
        .flat_map(|&(a, b)| {
            use std::iter::once;
            once(a).chain(once(b))
        })
        .collect::<FnvHashSet<_>>();

    info!("Finding ultrabubble path indices");
    let path_indices = variants::bubble_path_indices(&all_paths, &ultrabubble_nodes);

    info!(
        "Identifying variants in {} ultrabubbles",
        ultrabubbles.len() - remaining_ultrabubbles.len()
    );

    let mut all_vcf_records = Vec::new();

    let var_config = variants::VariantConfig {
        ignore_inverted_paths: args.ignore_inverted_paths,
    };

    let path_count = representative_paths.len();
    for (ix, (path_name, bubbles)) in representative_paths.into_iter().enumerate() {
        let ref_path_name: FnvHashSet<BString> =
            [path_name].iter().cloned().collect::<FnvHashSet<_>>();
        all_vcf_records.par_extend(
            bubbles
                .par_iter()
                .filter_map(|&(from, to)| {
                    let vars = variants::detect_variants_in_sub_paths(
                        &var_config,
                        &segment_map,
                        Some(&ref_path_name),
                        &all_paths,
                        &path_indices,
                        from,
                        to,
                    )?;

                    let vcf_records = variants::variant_vcf_record(&vars);
                    Some(vcf_records)
                })
                .flatten(),
        )
    }

    // all_vcf_records.sort_by(|v0, v1| v0.vcf_cmp(v1));

    let vcf_header = variants::vcf::VCFHeader::new(gfa_path);

    println!("{}", vcf_header);

    for vcf in all_vcf_records {
        println!("{}", vcf);
    }

    Ok(())
}
