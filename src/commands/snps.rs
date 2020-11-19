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
    variants::{PathData, PathStep, SNPRow},
};

use super::{load_gfa, Result};

/// Given a reference path from the GFA, by name, find and report the
/// SNPs for all other paths compared to the reference. Uses the
/// graph's ultrabubbles to identify areas of variation.
#[derive(StructOpt, Debug)]
pub struct SNPArgs {
    /// Load ultrabubbles from a file instead of calculating them.
    #[structopt(
        name = "ultrabubbles file",
        long = "ultrabubbles",
        short = "u"
    )]
    ultrabubbles_file: Option<PathBuf>,
    #[structopt(name = "name of reference path", long = "ref", short = "r")]
    /// The name of the path to be used as reference.
    ref_path: String,
}

fn build_snp_reference_bubbles(
    path: &[PathStep],
    positions: &mut [usize],
) -> Vec<(usize, usize)> {
    positions.sort();

    let mut res = Vec::with_capacity(positions.len());

    let mut steps_iter = path.iter();

    for &snp_pos in positions.iter() {
        if let Some(&(node, pos, _)) =
            steps_iter.find(|(_, pos, _)| *pos == snp_pos)
        {
            res.push((node, pos));
        }
    }

    res.shrink_to_fit();
    res
}

pub fn gfa2snps(gfa_path: &PathBuf, args: SNPArgs) -> Result<()> {
    let ref_path_name = BString::from(args.ref_path);

    let mut snps = vec![3037usize, 14408, 19724, 23403, 26258, 29573];

    let path_data = {
        let gfa: GFA<usize, ()> = load_gfa(&gfa_path)?;

        if gfa.paths.len() < 2 {
            panic!("GFA must contain at least two paths");
        }

        info!("GFA has {} paths", gfa.paths.len());

        variants::gfa_path_data(gfa)
    };

    let path_ix = path_data
        .path_names
        .iter()
        .position(|name| name == &ref_path_name)
        .unwrap();

    let path = &path_data.paths[path_ix];

    let snp_nodes = build_snp_reference_bubbles(path, &mut snps);

    println!("{:5}  {:5}", "Pos", "Node");
    for (node, pos) in snp_nodes {
        println!("{:5}  {:5}", pos, node);
    }

    Ok(())
}

pub fn _gfa2snps(gfa_path: &PathBuf, args: SNPArgs) -> Result<()> {
    let ref_path_name = BString::from(args.ref_path);

    let path_data = {
        let gfa: GFA<usize, ()> = load_gfa(&gfa_path)?;

        if gfa.paths.len() < 2 {
            panic!("GFA must contain at least two paths");
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

    let p_bar = progress_bar(ultrabubbles.len(), false);

    let mut path_snp_rows: FnvHashMap<BString, Vec<SNPRow>> =
        FnvHashMap::default();

    info!("Using reference path {}", ref_path_name);

    let ref_path_ix = path_data
        .path_names
        .iter()
        .position(|name| name == &ref_path_name)
        .expect("Reference path does not exist in graph");

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
