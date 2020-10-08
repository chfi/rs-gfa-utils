use clap::arg_enum;
use structopt::{clap::ArgGroup, StructOpt};

use bstr::{io::*, BString, ByteSlice, ByteVec};
use std::{
    fs::File,
    io::{BufReader, Read, Write},
    path::PathBuf,
};

use gfa::{
    gfa::{name_conversion::NameMap, GFA},
    optfields::OptionalFields,
    parser::GFAParser,
    writer::{gfa_string, write_gfa},
};

use fnv::{FnvHashMap, FnvHashSet};

use handlegraph::hashgraph::HashGraph;

use gfautil::{edges, gaf_convert, subgraph, variants};

arg_enum! {
    #[derive(Debug, PartialEq)]
    enum SubgraphBy {
        Paths,
        Segments,
    }
}

/// Generate a subgraph of the input GFA.
///
/// The output will be the lines of the input GFA that include the
/// provided segment or path names.
#[derive(StructOpt, Debug)]
#[structopt(group = ArgGroup::with_name("names").required(true))]
struct SubgraphArgs {
    /// Choose between providing a list of path names, or a list of
    /// components of segment names
    #[structopt(name = "paths|segments", possible_values = &["paths", "segments"], case_insensitive = true)]
    subgraph_by: SubgraphBy,
    /// File containing a list of names
    #[structopt(
        name = "File containing names",
        long = "file",
        group = "names"
    )]
    file: Option<PathBuf>,
    /// Provide a list of names on the command line
    #[structopt(name = "List of names", long = "names", group = "names")]
    list: Option<Vec<String>>,
}

/// Convert a file of GAF records into PAF records.
///
/// The provided GFA file should be the same as the one used to create the GAF.
#[derive(StructOpt, Debug)]
struct GAF2PAFArgs {
    #[structopt(name = "path to GAF file", long = "gaf", parse(from_os_str))]
    gaf: PathBuf,
    #[structopt(name = "PAF output paf", short = "o", long = "paf")]
    out: Option<PathBuf>,
}

#[derive(StructOpt, Debug)]
/// Convert a GFA with string names to one with integer names, and
/// back.
struct GfaIdConvertOptions {
    /// Path to a name map that was previously generated for the given GFA.
    /// Required if transforming to the original segment names. If not
    /// provided, a new map is generated and saved to disk.
    #[structopt(
        name = "path to name map",
        long = "namemap",
        parse(from_os_str),
        required_unless("to_usize")
    )]
    name_map_path: Option<PathBuf>,

    #[structopt(name = "convert to integer names", long = "to-int")]
    to_usize: bool,

    #[structopt(name = "check result hash", long = "hash")]
    check_hash: bool,
}

/// Output a VCF for the given GFA, using the graph's ultrabubbles to
/// identify areas of variation. (experimental!)
#[derive(StructOpt, Debug)]
struct GFAToVCFArgs {
    /// Don't compare two paths if their start and end orientations
    /// don't match each other
    #[structopt(name = "ignore inverted paths", long = "no-inv")]
    ignore_inverted_paths: bool,
}

fn gfa_to_name_map_path(path: &PathBuf) -> PathBuf {
    let mut new_path: PathBuf = path.clone();
    let old_name = new_path.file_stem().and_then(|p| p.to_str()).unwrap();
    let new_name = format!("{}.name_map.json", old_name);
    new_path.set_file_name(&new_name);
    new_path
}

fn converted_gfa_path(path: &PathBuf) -> PathBuf {
    let mut new_path: PathBuf = path.clone();
    let old_name = new_path.file_stem().and_then(|p| p.to_str()).unwrap();
    let new_name = format!("{}.uint_ids.gfa", old_name);
    new_path.set_file_name(&new_name);
    new_path
}

fn restored_gfa_path(path: &PathBuf) -> PathBuf {
    let mut new_path: PathBuf = path.clone();
    let old_name = new_path.file_stem().and_then(|p| p.to_str()).unwrap();
    let new_name = format!("{}.str_ids.gfa", old_name);
    new_path.set_file_name(&new_name);
    new_path
}

#[derive(StructOpt, Debug)]
enum Command {
    Subgraph(SubgraphArgs),
    EdgeCount,
    #[structopt(name = "gaf2paf")]
    Gaf2Paf(GAF2PAFArgs),
    GfaSegmentIdConversion(GfaIdConvertOptions),
    #[structopt(name = "gaf2vcf")]
    GfaToVcf(GFAToVCFArgs),
}

#[derive(StructOpt, Debug)]
struct Opt {
    #[structopt(
        name = "input GFA file",
        short,
        long = "gfa",
        parse(from_os_str)
    )]
    in_gfa: PathBuf,
    #[structopt(subcommand)]
    command: Command,
}

fn byte_lines_iter<'a, R: Read + 'a>(
    reader: R,
) -> Box<dyn Iterator<Item = Vec<u8>> + 'a> {
    Box::new(BufReader::new(reader).byte_lines().map(|l| l.unwrap()))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let opt = Opt::from_args();

    match opt.command {
        Command::GfaToVcf(var_args) => {
            let parser = GFAParser::new();
            let gfa: GFA<usize, ()> = parser.parse_file(&opt.in_gfa).unwrap();

            let segment_map: FnvHashMap<usize, &[u8]> = gfa
                .segments
                .iter()
                .map(|seg| (seg.name, seg.sequence.as_ref()))
                .collect();

            let all_paths =
                variants::gfa_paths_with_offsets(&gfa, &segment_map);

            let ultrabubbles = gfautil::ultrabubbles::gfa_ultrabubbles(&gfa);
            let ultrabubble_nodes = ultrabubbles
                .iter()
                .flat_map(|&(a, b)| {
                    use std::iter::once;
                    once(a).chain(once(b))
                })
                .collect::<FnvHashSet<_>>();

            let path_indices =
                variants::bubble_path_indices(&all_paths, &ultrabubble_nodes);

            let mut all_vcf_records = Vec::new();

            let var_config = variants::VariantConfig {
                ignore_inverted_paths: var_args.ignore_inverted_paths,
            };

            for &(from, to) in ultrabubbles.iter() {
                let vars = variants::detect_variants_in_sub_paths(
                    &var_config,
                    &segment_map,
                    &all_paths,
                    &path_indices,
                    from,
                    to,
                );

                let vcf_records = variants::variant_vcf_record(&vars);
                all_vcf_records.extend(vcf_records);

                /*
                let from_indices = path_indices.get(&from).unwrap();
                let to_indices = path_indices.get(&to).unwrap();

                let sub_paths: FnvHashMap<
                    &BStr,
                    &[(usize, usize, Orientation)],
                > = all_paths
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
                */
            }

            all_vcf_records.sort_by(|v0, v1| v0.vcf_cmp(v1));

            let vcf_header = variants::vcf::VCFHeader::new(&opt.in_gfa);

            println!("{}", vcf_header);

            for vcf in all_vcf_records {
                println!("{}", vcf);
            }
        }

        Command::Subgraph(subgraph_args) => {
            let parser = GFAParser::new();
            let gfa: GFA<BString, OptionalFields> =
                parser.parse_file(&opt.in_gfa).unwrap();

            let names: Vec<Vec<u8>> = if let Some(list) = subgraph_args.list {
                list.into_iter().map(|s| s.bytes().collect()).collect()
            } else {
                let in_lines = if let Some(path) = subgraph_args.file {
                    byte_lines_iter(File::open(path).unwrap())
                } else {
                    byte_lines_iter(std::io::stdin())
                };

                if subgraph_args.subgraph_by == SubgraphBy::Segments {
                    in_lines
                        .flat_map(|line| {
                            line.split_str("\t")
                                .map(Vec::from_slice)
                                .collect::<Vec<_>>()
                        })
                        .collect()
                } else {
                    in_lines.collect()
                }
            };

            let new_gfa = match subgraph_args.subgraph_by {
                SubgraphBy::Paths => subgraph::paths_new_subgraph(&gfa, &names),
                SubgraphBy::Segments => {
                    subgraph::segments_subgraph(&gfa, &names)
                }
            };
            println!("{}", gfa_string(&new_gfa));
        }
        Command::Gaf2Paf(args) => {
            let parser = GFAParser::new();
            let gfa: GFA<BString, OptionalFields> =
                parser.parse_file(&opt.in_gfa).unwrap();
            let paf_lines = gaf_convert::gaf_to_paf(gfa, &args.gaf);

            if let Some(out_path) = args.out {
                let mut out_file = File::create(&out_path)
                    .expect("Error creating PAF output file");

                paf_lines.iter().for_each(|p| {
                    writeln!(out_file, "{}", p).unwrap();
                });
            } else {
                paf_lines.iter().for_each(|p| println!("{}", p));
            }
        }
        Command::EdgeCount => {
            let parser = GFAParser::new();
            let gfa: GFA<usize, ()> = parser.parse_file(&opt.in_gfa).unwrap();

            let hashgraph = HashGraph::from_gfa(&gfa);
            let edge_counts = edges::graph_edge_count(&hashgraph);
            println!("nodeid,inbound,outbound,total");
            edge_counts
                .iter()
                .for_each(|(id, i, o, t)| println!("{},{},{},{}", id, i, o, t));
        }
        Command::GfaSegmentIdConversion(conv_opt) => {
            // Converting from string to integer names
            if !conv_opt.to_usize && conv_opt.name_map_path.is_none() {
                eprintln!("this shouldn't happen");
            }

            if conv_opt.to_usize {
                let parser = GFAParser::new();
                let gfa: GFA<BString, OptionalFields> =
                    parser.parse_file(&opt.in_gfa).unwrap();

                let name_map = if let Some(ref path) = conv_opt.name_map_path {
                    let map = NameMap::load_json(&path)?;
                    map
                } else {
                    NameMap::build_from_gfa(&gfa)
                };

                if let Some(new_gfa) =
                    name_map.gfa_bstring_to_usize(&gfa, conv_opt.check_hash)
                {
                    let new_gfa_path = converted_gfa_path(&opt.in_gfa);
                    let mut new_gfa_file = File::create(new_gfa_path.clone())?;
                    let mut gfa_str = String::new();
                    write_gfa(&new_gfa, &mut gfa_str);
                    writeln!(new_gfa_file, "{}", gfa_str)?;
                    println!(
                        "Saved converted GFA to {}",
                        new_gfa_path.display()
                    );

                    if conv_opt.name_map_path.is_none() {
                        let name_map_path = gfa_to_name_map_path(&opt.in_gfa);
                        name_map.save_json(&name_map_path)?;
                        println!(
                            "Saved new name map to {}",
                            name_map_path.display()
                        );
                    }
                } else {
                    println!("Could not convert the GFA segment IDs");
                }
            } else {
                // Converting from integer to string names
                let name_map_path = conv_opt
                    .name_map_path
                    .expect("Need name map to convert back");
                let name_map = NameMap::load_json(&name_map_path)?;

                let parser = GFAParser::new();
                let gfa: GFA<usize, OptionalFields> =
                    parser.parse_file(&opt.in_gfa).unwrap();

                let new_gfa: GFA<BString, OptionalFields> =
                    name_map.gfa_usize_to_bstring(&gfa).expect(
                        "Error during conversion -- is it the right name map?",
                    );

                let new_gfa_path = restored_gfa_path(&opt.in_gfa);
                let mut new_gfa_file = File::create(new_gfa_path.clone())?;
                let mut gfa_str = String::new();
                write_gfa(&new_gfa, &mut gfa_str);
                writeln!(new_gfa_file, "{}", gfa_str)?;
                println!("Saved restored GFA to {}", new_gfa_path.display());
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn name_map_path_correct() {
        let gfa_path =
            PathBuf::from("../../some_folder/another/some_gfa_file.gfa");
        let new_path = gfa_to_name_map_path(&gfa_path);
        assert_eq!(
            Some("../../some_folder/another/some_gfa_file.name_map.json"),
            new_path.to_str()
        );
    }

    #[test]
    fn converted_gfa_path_correct() {
        let gfa_path = PathBuf::from("some_gfa_file.gfa");
        let new_path = converted_gfa_path(&gfa_path);
        assert_eq!(Some("some_gfa_file.uint_ids.gfa"), new_path.to_str());
    }
}
