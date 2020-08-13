use std::{fs::File, io, io::BufReader, path::PathBuf};

use bstr::{io::*, BString, ByteSlice, ByteVec};

use gfa::{
    gafpaf::{parse_gaf, GAFPath, GAFStep, CIGAR},
    gfa::GFA,
    optfields::{OptFieldVal, OptFields, OptionalFields},
    parser::GFAParser,
};

use gfautil::gaf_convert::gaf_to_paf;

type GAF = gfa::gafpaf::GAF<OptionalFields>;
type PAF = gfa::gafpaf::PAF<OptionalFields>;

// fn load_gfa(path: &str) -> GFA {
//     let path = PathBuf::from(path);
//     let parser = GFAParser::new();
//     parser.parse_file(path).unwrap()
// }

fn load_pafs(gfa_path: &str, gaf_path: &str) -> Vec<PAF> {
    let gfa_path = PathBuf::from(gfa_path);
    let parser = GFAParser::new();
    let gfa: GFA<BString, OptionalFields> =
        parser.parse_file(gfa_path).unwrap();

    let gaf_path = PathBuf::from(gaf_path);
    let pafs = gaf_to_paf(gfa, &gaf_path);

    pafs
}

// fn load_gaf(path: &PathBuf) -> io::Result<GAF> {
//     let
// }

#[test]
fn gafpaf_no_overlaps() {
    let pafs = load_pafs("./tests/data/ov1.gfa", "./tests/data/ov1.gaf");

    for paf in pafs {
        println!("{}", paf);
    }
}

// fn gafpaf_overlaps() {}
