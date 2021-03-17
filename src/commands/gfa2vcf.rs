use bstr::{io::*, BString, ByteSlice};
use fnv::{FnvHashMap, FnvHashSet};
use indicatif::{
    ParallelProgressIterator, ProgressBar, ProgressIterator, ProgressStyle,
};
use rayon::prelude::*;
use std::{fs::File, io::BufReader, path::PathBuf};
use structopt::StructOpt;

#[allow(unused_imports)]
use log::{debug, info, log_enabled, warn};

use gfa::gfa::GFA;

use crate::{util::progress_bar, variants, variants::PathStep};

use super::{load_gfa, Result};

/// Output a VCF for the given GFA, using the graph's ultrabubbles to
/// identify areas of variation.
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

    let path_data = {
        let gfa: GFA<usize, ()> = load_gfa(&gfa_path)?;

        if gfa.paths.len() < 2 {
            panic!("GFA must contain at least two paths");
        }

        if let Some(ref_paths) = ref_path_names.as_ref() {
            let gfa_paths = gfa
                .paths
                .iter()
                .map(|path| path.path_name.as_bstr())
                .collect::<FnvHashSet<_>>();

            for path in ref_paths.iter() {
                if !gfa_paths.contains(path.as_bstr()) {
                    eprintln!(
                        "Reference path does not exist in graph: {}",
                        path.as_bstr()
                    );
                    std::process::exit(1);
                }
            }
        }

        info!("GFA has {} paths", gfa.paths.len());

        variants::gfa_path_data(gfa)
    };

    let mut ultrabubbles = if let Some(path) = &args.ultrabubbles_file {
        super::saboten::load_ultrabubbles(path)
    } else {
        super::saboten::find_ultrabubbles(gfa_path)
    }?;

    info!("Using {} ultrabubbles", ultrabubbles.len());

    ultrabubbles.sort();

    let ultrabubble_nodes = ultrabubbles
        .iter()
        .flat_map(|&(a, b)| {
            use std::iter::once;
            once(a).chain(once(b))
        })
        .collect::<FnvHashSet<_>>();

    let path_indices =
        variants::bubble_path_indices(&path_data.paths, &ultrabubble_nodes);

    let mut all_vcf_records = Vec::new();

    let var_config = variants::VariantConfig {
        ignore_inverted_paths: args.ignore_inverted_paths,
    };

    info!(
        "Identifying variants in {} ultrabubbles",
        ultrabubbles.len()
    );

    let p_bar = progress_bar(ultrabubbles.len(), false);

    all_vcf_records.par_extend(
        ultrabubbles
            .par_iter()
            .progress_with(p_bar)
            .filter_map(|&(from, to)| {
                let vars = variants::detect_variants_in_sub_paths(
                    &var_config,
                    &path_data,
                    ref_path_names.as_ref(),
                    &path_indices,
                    from,
                    to,
                )?;

                let vcf_records = variants::variant_vcf_record(&vars);
                Some(vcf_records)
            })
            .flatten(),
    );
    info!("Variant identification complete");

    all_vcf_records.sort_by(|v0, v1| v0.vcf_cmp(v1));
    all_vcf_records.dedup();

    info!("Writing {} unique VCF records", all_vcf_records.len());

    let vcf_header = variants::vcf::VCFHeader::new(gfa_path);

    println!("{}", vcf_header);

    for vcf in all_vcf_records {
        println!("{}", vcf);
    }

    Ok(())

    /*
    for (path_name, bubbles) in representative_paths.into_iter().progress_with(p_bar) {
        let ref_path_name: FnvHashSet<BString> =
            [path_name].iter().cloned().collect::<FnvHashSet<_>>();
        let before_len = all_vcf_records.len();
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
        );
        let after_len = all_vcf_records.len();
        // info!("\tfound {} variants", after_len - before_len);
    }
    */
}

#[allow(dead_code)]
fn find_representative_paths(
    ultrabubbles: &[(u64, u64)],
    all_paths: &FnvHashMap<BString, Vec<PathStep>>,
) -> Vec<(BString, Vec<(u64, u64)>)> {
    let mut representative_paths = Vec::new();

    let mut remaining_ultrabubbles: FnvHashMap<u64, u64> =
        ultrabubbles.iter().copied().collect();

    info!("Building set of reference paths");
    let p_bar = ProgressBar::new(all_paths.len() as u64);
    p_bar.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:80} {pos:>7}/{len:7}")
            .progress_chars("##-"),
    );

    for (path, steps) in all_paths.iter().progress_with(p_bar) {
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
                maybe_contained.iter().filter_map(move |&(a, b)| {
                    if a == y {
                        Some((b, a))
                    } else {
                        None
                    }
                })
            })
            .collect::<Vec<_>>();

        for &(x, _y) in contained.iter() {
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

    representative_paths
}
