pub mod convert_names;
pub mod gaf2paf;
pub mod gfa2vcf;
pub mod saboten;
pub mod snps;
pub mod stats;
pub mod subgraph;

use std::io::{BufReader, Read};

use bstr::io::*;
use gfa::{
    gfa::{SegmentId, GFA},
    optfields::OptFields,
    parser::GFAParser,
};

#[allow(unused_imports)]
use log::{debug, info, warn};

pub type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

pub fn byte_lines_iter<'a, R: Read + 'a>(
    reader: R,
) -> Box<dyn Iterator<Item = Vec<u8>> + 'a> {
    Box::new(BufReader::new(reader).byte_lines().map(|l| l.unwrap()))
}

pub fn load_gfa<N, T, P>(path: P) -> Result<GFA<N, T>>
where
    N: SegmentId,
    T: OptFields,
    P: AsRef<std::path::Path>,
{
    let parser = GFAParser::new();
    info!("Parsing GFA from {}", path.as_ref().display());
    let gfa = parser.parse_file(path.as_ref())?;
    Ok(gfa)
}
