use bstr::BString;
use fnv::{FnvHashMap, FnvHashSet};
use std::path::PathBuf;
use structopt::StructOpt;

use indicatif::ProgressIterator;

use gfa::gfa::GFA;

#[allow(unused_imports)]
use log::{debug, info, log_enabled, warn};

use crate::{
    util::progress_bar,
    variants,
    variants::{PathStep, SNPRow},
};

use super::{load_gfa, Result};

/// Given a reference path from the GFA, by name, find and report the
/// SNPs for all other paths compared to the reference. Uses the
/// graph's ultrabubbles to identify areas of variation.
#[derive(StructOpt, Debug)]
pub struct SNPArgs {
    #[structopt(name = "name of reference path", long = "ref", short = "r")]
    /// The name of the path to be used as reference.
    ref_path: String,
    /// A list of SNP positions to use.
    #[structopt(
        name = "SNP positions",
        long = "snps",
        required_unless_one(&["SNP positions file", "ultrabubbles file"])
    )]
    snp_positions: Option<Vec<usize>>,
    /// Path to a file containing SNP positions to use, one position
    /// per line.
    #[structopt(
        name = "SNP positions file",
        long = "snps-file",
        required_unless_one(&["SNP positions", "ultrabubbles file"])
    )]
    snp_positions_file: Option<PathBuf>,
    /// Path to a file containing bubbles to use, if not providing SNP
    /// positions.
    #[structopt(
        name = "ultrabubbles file",
        long = "ultrabubbles",
        short = "u",
        required_unless_one(&["SNP positions", "SNP positions file"])
    )]
    ultrabubbles_file: Option<PathBuf>,
}

fn snp_positions(args: &SNPArgs) -> Result<Vec<usize>> {
    let mut res = Vec::new();

    if let Some(positions) = args.snp_positions.as_ref() {
        res.extend(positions.iter().copied());
    }

    if let Some(file_path) = args.snp_positions_file.as_ref() {
        let positions = load_snp_positions_file(file_path)?;
        res.extend(positions);
    }

    if res.is_empty() {
        panic!("No SNPs were provided");
    }

    Ok(res)
}

fn load_snp_positions_file(file_path: &PathBuf) -> Result<Vec<usize>> {
    use bstr::{io::*, ByteSlice};
    use std::{fs::File, io::BufReader};

    let mut res = Vec::new();

    let file = File::open(file_path)?;
    let reader = BufReader::new(file);

    for line in reader.byte_lines() {
        let line = line?;
        let line = line.trim().to_str()?;
        let pos = line.parse::<usize>()?;
        res.push(pos);
    }

    Ok(res)
}

fn build_snp_reference_bubbles(
    path: &[PathStep],
    positions: &mut [usize],
) -> Vec<(u64, u64)> {
    let mut res = Vec::with_capacity(positions.len());

    positions.sort();
    let mut steps_iter = path.iter().enumerate();

    for &snp_pos in positions.iter() {
        if let Some((ix, _)) =
            steps_iter.find(|(_, (_, pos, _))| *pos == snp_pos)
        {
            if ix > 0 && ix < path.len() {
                let (prev, _, _) = path[ix - 1];
                let (next, _, _) = path[ix + 1];
                res.push((prev as u64, next as u64));
            }
        }
    }

    res.shrink_to_fit();
    res
}

pub fn gfa2snps(gfa_path: &PathBuf, args: SNPArgs) -> Result<()> {
    let ref_path_name: BString = BString::from(args.ref_path.as_str());

    let path_data = {
        let gfa: GFA<usize, ()> = load_gfa(&gfa_path)?;

        if gfa.paths.len() < 2 {
            panic!("GFA must contain at least two paths");
        }

        info!("GFA has {} paths", gfa.paths.len());

        variants::gfa_path_data(gfa)
    };

    info!("Using reference path {}", ref_path_name);

    let ref_path_ix = path_data
        .path_names
        .iter()
        .position(|name| name == &ref_path_name)
        .expect("Reference path does not exist in graph");

    let ref_path = &path_data.paths[ref_path_ix];

    let ultrabubbles = if let Ok(mut positions) = snp_positions(&args) {
        Ok(build_snp_reference_bubbles(&ref_path, &mut positions))
    } else if let Some(path) = &args.ultrabubbles_file {
        super::saboten::load_ultrabubbles(path)
    } else {
        unreachable!()
    }?;

    info!("Found ultrabubbles for {} SNPs", ultrabubbles.len());

    let ultrabubble_nodes = ultrabubbles
        .iter()
        .flat_map(|&(a, b)| {
            use std::iter::once;
            once(a).chain(once(b))
        })
        .collect::<FnvHashSet<_>>();

    let path_indices =
        variants::bubble_path_indices(&path_data.paths, &ultrabubble_nodes);

    let p_bar = progress_bar(ultrabubbles.len(), false);

    let mut path_snp_rows: FnvHashMap<BString, Vec<SNPRow>> =
        FnvHashMap::default();

    for &(from, to) in ultrabubbles.iter().progress_with(p_bar) {
        let results = variants::find_snps_in_sub_paths(
            &path_data,
            ref_path_ix,
            &path_indices,
            from,
            to,
        );

        if let Some(snp_results) = results {
            for (name, snp_rows) in snp_results.into_iter() {
                let entry = path_snp_rows.entry(name).or_default();
                entry.extend(snp_rows);
            }
        }
    }

    println!("path\treference base\treference pos\tquery base\tquery pos");
    for (name, snp_rows) in path_snp_rows.into_iter() {
        for snp in snp_rows.into_iter() {
            let ref_base = char::from(snp.ref_base);
            let query_base = char::from(snp.query_base);
            println!(
                "{}\t{}\t{}\t{}\t{}",
                &name, ref_base, snp.ref_pos, query_base, snp.query_pos
            );
        }
    }

    Ok(())
}
