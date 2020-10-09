use bstr::BString;
use std::{fs::File, io::Write, path::PathBuf};
use structopt::StructOpt;

use gfa::{gfa::GFA, optfields::OptionalFields};

use crate::gaf_convert;

/// Convert a file of GAF records into PAF records.
///
/// The provided GFA file should be the same as the one used to create the GAF.
#[derive(StructOpt, Debug)]
pub struct GAF2PAFArgs {
    #[structopt(name = "path to GAF file", long = "gaf", parse(from_os_str))]
    gaf: PathBuf,
    #[structopt(name = "PAF output paf", short = "o", long = "paf")]
    out: Option<PathBuf>,
}

pub fn gaf2paf(
    gfa: GFA<BString, OptionalFields>,
    args: &GAF2PAFArgs,
) -> Result<(), Box<dyn std::error::Error>> {
    let paf_lines = gaf_convert::gaf_to_paf(gfa, &args.gaf);

    if let Some(out_path) = &args.out {
        let mut out_file =
            File::create(&out_path).expect("Error creating PAF output file");

        paf_lines.iter().for_each(|p| {
            writeln!(out_file, "{}", p).unwrap();
        });
    } else {
        paf_lines.iter().for_each(|p| println!("{}", p));
    }

    Ok(())
}
