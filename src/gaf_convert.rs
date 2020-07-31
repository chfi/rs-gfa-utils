use std::{
    fs::File,
    io::{prelude::*, BufReader},
    path::Path,
};

use bstr::{io::*, BStr, BString, ByteSlice, ByteVec};

use gfa::{
    gafpaf::{parse_gaf, CIGAROp, GAFPath, GAFStep, CIGAR, GAF, PAF},
    gfa::GFA,
    optfields::OptFields,
};

pub fn gaf_to_paf<T: OptFields>(
    gfa: &GFA<BString, T>,
    gaf_path: &Path,
) -> Vec<PAF<()>> {
    let file = File::open(gaf_path).unwrap();
    let lines = BufReader::new(file).byte_lines().map(|l| l.unwrap());
    let mut gafs: Vec<GAF<()>> = Vec::new();

    for (i, line) in lines.enumerate() {
        let fields = line.split_str(b"\t");
        if let Some(gaf) = parse_gaf(fields) {
            gafs.push(gaf);
        } else {
            eprintln!("Error parsing GAF line {}", i);
        }
    }

    gafs.iter().for_each(|gaf| {
        println!("{}", gaf);
    });

    Vec::new()
}
