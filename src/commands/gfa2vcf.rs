use bstr::io::*;
use bstr::BString;
use fnv::{FnvHashMap, FnvHashSet};
use std::io::{BufReader, Read};
use std::{fs::File, io::Write, path::PathBuf};
use structopt::StructOpt;

use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use rayon::prelude::*;

use gfa::gfa::GFA;

#[allow(unused_imports)]
use log::{debug, info, warn};

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
        name = "file containing paths to include",
        long = "paths-file"
    )]
    included_paths_file: Option<PathBuf>,
    #[structopt(name = "list of paths to include", long = "paths")]
    included_paths_vec: Option<Vec<String>>,
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
    let included_paths_list =
        args.included_paths_vec.map(paths_list).unwrap_or_default();

    let included_paths_file = args
        .included_paths_file
        .map(load_paths_file)
        .transpose()?
        .unwrap_or_default();

    let included_paths: FnvHashSet<BString> = included_paths_list
        .into_iter()
        .chain(included_paths_file.into_iter())
        .collect();

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
    let ultrabubbles = if let Some(path) = &args.ultrabubbles_file {
        let ub = crate::ultrabubbles::load_ultrabubbles(path)?;
        ub
    } else {
        crate::ultrabubbles::gfa_ultrabubbles(&gfa)
    };

    info!("Using {} ultrabubbles", ultrabubbles.len());

    let ultrabubble_nodes = ultrabubbles
        .iter()
        .flat_map(|&(a, b)| {
            use std::iter::once;
            once(a).chain(once(b))
        })
        .collect::<FnvHashSet<_>>();

    info!("Finding ultrabubble path indices");
    let path_indices =
        variants::bubble_path_indices(&all_paths, &ultrabubble_nodes);

    let mut all_vcf_records = Vec::new();

    let var_config = variants::VariantConfig {
        ignore_inverted_paths: args.ignore_inverted_paths,
    };

    info!(
        "Identifying variants in {} ultrabubbles",
        ultrabubbles.len()
    );

    let p_bar = ProgressBar::new(ultrabubbles.len() as u64);
    p_bar.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise}] {bar:80} {pos:>7}/{len:7}")
            .progress_chars("##-"),
    );

    all_vcf_records.par_extend(
        ultrabubbles
            .par_iter()
            .progress_with(p_bar)
            .filter_map(|&(from, to)| {
                let vars = variants::detect_variants_in_sub_paths(
                    &var_config,
                    &segment_map,
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

    all_vcf_records.sort_by(|v0, v1| v0.vcf_cmp(v1));

    let vcf_header = variants::vcf::VCFHeader::new(gfa_path);

    println!("{}", vcf_header);

    for vcf in all_vcf_records {
        println!("{}", vcf);
    }

    Ok(())
}
