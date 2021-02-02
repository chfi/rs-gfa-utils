use std::path::PathBuf;

use gfa::{
    cigar::CIGAR,
    gfa::GFA,
    optfields::{OptFieldVal, OptFields, OptionalFields},
    parser::GFAParser,
};

use gfautil::gaf_convert::gaf_to_paf;

type PAF = gfa::gafpaf::PAF<OptionalFields>;

fn load_pafs(gfa_path: &str, gaf_path: &str) -> Vec<PAF> {
    let gfa_path = PathBuf::from(gfa_path);
    let parser = GFAParser::new();
    let gfa: GFA<Vec<u8>, OptionalFields> =
        parser.parse_file(gfa_path).unwrap();

    let gaf_path = PathBuf::from(gaf_path);
    let pafs = gaf_to_paf(gfa, &gaf_path);

    pafs
}

fn get_cigar(opts: &OptionalFields) -> Option<CIGAR> {
    let cg = opts.get_field(b"cg")?;
    if let OptFieldVal::Z(cg) = &cg.value {
        CIGAR::from_bytestring(&cg)
    } else {
        None
    }
}

fn compare_paf_query(
    paf: &PAF,
    q_name: &str,
    q_len: usize,
    q_range: (usize, usize),
) {
    assert_eq!(paf.query_seq_name, q_name.as_bytes());
    assert_eq!(paf.query_seq_len, q_len);
    assert_eq!(paf.query_seq_range, q_range);
}

fn compare_paf_target(
    paf: &PAF,
    t_name: &str,
    t_len: usize,
    t_range: (usize, usize),
) {
    assert_eq!(paf.target_seq_name, t_name.as_bytes());
    assert_eq!(paf.target_seq_len, t_len);
    assert_eq!(paf.target_seq_range, t_range);
}

fn compare_paf_rest(
    paf: &PAF,
    res_matches: usize,
    block_length: usize,
    cigar: &str,
) {
    let paf_cigar = get_cigar(&paf.optional).unwrap();
    assert_eq!(paf.residue_matches, res_matches);
    assert_eq!(paf.block_length, block_length);
    assert_eq!(paf_cigar.to_string(), cigar);
}

#[test]
fn gafpaf_no_overlaps() {
    let pafs = load_pafs("./tests/data/ov1.gfa", "./tests/data/ov1.gaf");
    let mut iter = pafs.iter();

    // read1
    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read1", 6, (0, 1));
    compare_paf_target(&paf, "2", 3, (2, 3));
    compare_paf_rest(&paf, 1, 1, "1M");

    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read1", 6, (1, 5));
    compare_paf_target(&paf, "3", 4, (0, 4));
    compare_paf_rest(&paf, 4, 4, "4M");

    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read1", 6, (5, 6));
    compare_paf_target(&paf, "4", 5, (0, 1));
    compare_paf_rest(&paf, 1, 1, "1M");

    // read2
    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read2", 7, (0, 2));
    compare_paf_target(&paf, "2", 3, (1, 3));
    compare_paf_rest(&paf, 2, 2, "2M");

    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read2", 7, (2, 6));
    compare_paf_target(&paf, "5", 4, (0, 4));
    compare_paf_rest(&paf, 4, 4, "4M");

    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read2", 7, (6, 7));
    compare_paf_target(&paf, "6", 4, (0, 1));
    compare_paf_rest(&paf, 1, 1, "1M");

    assert!(iter.next().is_none());
}

#[test]
fn gafpaf_overlaps() {
    let pafs = load_pafs("./tests/data/ov2.gfa", "./tests/data/ov2.gaf");
    let mut iter = pafs.iter();

    // read1
    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read1", 6, (0, 1));
    compare_paf_target(&paf, "2", 3, (2, 3));
    compare_paf_rest(&paf, 1, 1, "1M");

    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read1", 6, (1, 5));
    compare_paf_target(&paf, "3", 4, (0, 4));
    compare_paf_rest(&paf, 3, 4, "1I3M");

    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read1", 6, (5, 6));
    compare_paf_target(&paf, "4", 5, (0, 1));
    compare_paf_rest(&paf, 1, 1, "1M");

    // read2
    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read2", 7, (0, 2));
    compare_paf_target(&paf, "2", 3, (1, 3));
    compare_paf_rest(&paf, 2, 2, "2M");

    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read2", 7, (2, 6));
    compare_paf_target(&paf, "5", 4, (0, 4));
    compare_paf_rest(&paf, 3, 4, "2M1I1M");

    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read2", 7, (6, 7));
    compare_paf_target(&paf, "6", 4, (0, 1));
    compare_paf_rest(&paf, 1, 1, "1M");

    assert!(iter.next().is_none());
}

#[test]
fn gafpaf_dels() {
    let pafs = load_pafs("./tests/data/ov1.gfa", "./tests/data/dels.gaf");
    let mut iter = pafs.iter();

    // read1
    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read1", 5, (0, 1));
    compare_paf_target(&paf, "2", 3, (2, 3));
    compare_paf_rest(&paf, 1, 1, "1M");

    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read1", 5, (1, 5));
    compare_paf_target(&paf, "3", 4, (0, 4));
    compare_paf_rest(&paf, 3, 4, "1M1D2M");

    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read1", 5, (5, 6));
    compare_paf_target(&paf, "4", 5, (0, 1));
    compare_paf_rest(&paf, 1, 1, "1M");

    // read2
    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read2", 6, (0, 2));
    compare_paf_target(&paf, "2", 3, (1, 3));
    compare_paf_rest(&paf, 2, 2, "2M");

    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read2", 6, (2, 6));
    compare_paf_target(&paf, "5", 4, (0, 4));
    compare_paf_rest(&paf, 3, 4, "1D3M");

    let paf = iter.next().unwrap();
    compare_paf_query(&paf, "read2", 6, (6, 7));
    compare_paf_target(&paf, "6", 4, (0, 1));
    compare_paf_rest(&paf, 1, 1, "1M");
}
