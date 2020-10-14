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

use gfa::gfa::GFA;

#[allow(unused_imports)]
use log::{debug, info, log_enabled, warn};

use crate::variants;

use super::{load_gfa, Result};

/// Output a VCF for the given GFA, using the graph's ultrabubbles to
/// identify areas of variation. (experimental!)
#[derive(StructOpt, Debug)]
pub struct GFA2VCFArgs {
    /// Load ultrabubbles from a file instead of calculating them.
    #[structopt(name = "ultrabubbles file", long = "ultrabubbles", short = "ub")]
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

    let mut remaining_ultrabubbles = ultrabubbles.clone();

    for (path, steps) in all_paths.iter() {
        let indices = remaining_ultrabubbles.iter().position(|&(x, y)| {
            let a = x as usize;
            let b = y as usize;
            let mut has_a = false;
            let mut has_b = false;
            for (u, _, _) in steps.iter() {
                if u == &a {
                    has_a = true;
                } else if u == &b {
                    has_b = true;
                }
                if has_a && has_b {
                    break;
                }
            }
            has_a && has_b
        });

        let path_name = path.clone().to_owned();
        let mut bubbles = Vec::new();
        for ix in indices {
            let ub = remaining_ultrabubbles.remove(ix);
            bubbles.push(ub);
        }
        let path_name_str: String = path_name.to_string();
        println!(
            "{:<40}\t{} bubbles\t{} remaining",
            path_name,
            bubbles.len(),
            remaining_ultrabubbles.len()
        );
        representative_paths.push((path_name, bubbles));
    }

    println!("{} paths", representative_paths.len());

    /*

    let ultrabubble_nodes = ultrabubbles
        .iter()
        .flat_map(|&(a, b)| {
            use std::iter::once;
            once(a).chain(once(b))
        })
        .collect::<FnvHashSet<_>>();

    let path_index = args.path_index.unwrap_or(0);

    let mut paths = gfa
        .paths
        .iter()
        .map(|p| p.iter().map(|(x, _)| x).collect::<Vec<_>>())
        .collect::<Vec<_>>();

    let mut first_path = paths.remove(path_index);

    // let mut first_path = gfa.paths[0].iter().map(|(x, _)| x).collect::<Vec<_>>();

    // let mut other_paths =

    let end_index = args.index_arg.unwrap_or(1000);

    println!("{} paths", gfa.paths.len());
    let mut same_vec: Vec<FnvHashSet<usize>> = Vec::new();
    for i in 0..end_index {
        let same = paths
            .iter()
            .enumerate()
            .filter_map(|(ix, steps)| {
                if steps[i] == first_path[i] {
                    Some(ix)
                } else {
                    None
                }
            })
            .collect::<FnvHashSet<usize>>();
        if i % 100 == 0 {
            println!("{}\t{}", i, same.len());
        }
        same_vec.push(same);
    }

    // let first_vec = same_vec[0].clone();
    let mut total_set = same_vec.pop().unwrap();
    println!("total set before: {}", total_set.len());

    for other in same_vec.into_iter() {
        let intersection = total_set
            .intersection(&other)
            .copied()
            .collect::<FnvHashSet<usize>>();
        total_set = intersection;
    }

    println!("intersection: {}", total_set.len());
    // println!("minimum: {:?}", total_set.iter().min());
    println!("maximum: {:?}", total_set.iter().max());

    */

    /*

    let mut path_vecs = Vec::new();
    for i in 0..=4 {
        let path = &gfa.paths[i];
        let path_vec = path.iter().collect::<Vec<_>>();
        println!("{}\t{}", i, path_vec.len());
        path_vecs.push(path_vec);
    }


    println!();
    // let mut step = 0;
    // let mut step = 3576;
    // let mut step = 6027;
    // let mut step = 8783;
    // let mut step = 12474;
    // let mut step = 14434;
    // let mut step = 15103;
    // let mut step = 18061;
    // let mut step = 21995;
    // let mut step = 24325;
    // let mut step = 25433;
    // let mut step = 27351;
    let mut step = 29000;
    loop {
        // let mut bases = path_vecs.iter().map(|v| v[step]).collect::<Vec<_>>();
        let mut bases = path_vecs
            .iter()
            .take(4)
            .map(|v| v[step])
            .collect::<Vec<_>>();
        bases.push(path_vecs[4][step + 1]);
        // bases[4] = path_vecs[
        let base_clone = bases.clone();
        bases.sort_by(|x, y| x.0.cmp(&y.0));
        bases.dedup();
        if bases.len() != 1 {
            for (i, (b, o)) in base_clone.iter().enumerate() {
                println!("{}\t{:4}{}", i, b, o);
            }
            break;
        }
        step += 1;
    }

    println!("first {} steps of first 5 paths were equal", step);

    for i in 0..5 {
        let path = &gfa.paths[i];
        let ten_bases = path
            .iter()
            .map(|(x, _)| x)
            .skip(step - 5)
            .take(25)
            .collect::<Vec<_>>();
        // let ten_bases = path.iter().take(25).collect::<Vec<_>>();

        println!(" ~~ {} ~~ ", i);
        for b in ten_bases.iter() {
            print!("{:>5} ", b);
        }
        println!();
        for b in ten_bases.iter() {
            let s = segment_map.get(&b).unwrap();
            let s: BString = s.iter().copied().collect();
            let s = s.to_string();
            print!("{:>5} ", s);
        }
        println!();
    }
    */

    /*
    println!("x\ty\toccurrences");
    for ((x, y), count) in bub_path_vec {
        println!("{}\t{}\t{}", x, y, count);
    }
    */

    /*
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
        p_bar.enable_steady_tick(1000);

        all_vcf_records.par_extend(
            ultrabubbles
                .par_iter()
                .progress_with(p_bar)
                .filter_map(|&(from, to)| {
                    let vars = variants::detect_variants_in_sub_paths(
                        &var_config,
                        &segment_map,
                        ref_path_names.as_ref(),
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
    */

    Ok(())
}
