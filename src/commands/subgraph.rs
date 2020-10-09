use clap::arg_enum;
use structopt::{clap::ArgGroup, StructOpt};

use bstr::{io::*, BString, ByteSlice, ByteVec};
use std::{
    fs::File,
    io::{BufReader, Read, Write},
    path::PathBuf,
};

use gfa::{
    gfa::{name_conversion::NameMap, SegmentId, GFA},
    optfields::{OptFields, OptionalFields},
    parser::GFAParser,
    writer::{gfa_string, write_gfa},
};

use handlegraph::hashgraph::HashGraph;

use crate::subgraph;

use super::{byte_lines_iter, Result};

use log::{debug, info, warn};

arg_enum! {
    #[derive(Debug, PartialEq)]
    pub enum SubgraphBy {
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
pub struct SubgraphArgs {
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

pub fn subgraph(
    gfa: &GFA<BString, OptionalFields>,
    args: &SubgraphArgs,
) -> Result<()> {
    let names: Vec<Vec<u8>> = if let Some(list) = &args.list {
        list.into_iter().map(|s| s.bytes().collect()).collect()
    } else {
        let in_lines = if let Some(path) = &args.file {
            byte_lines_iter(File::open(path).unwrap())
        } else {
            byte_lines_iter(std::io::stdin())
        };

        if args.subgraph_by == SubgraphBy::Segments {
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

    let new_gfa = match args.subgraph_by {
        SubgraphBy::Paths => subgraph::paths_new_subgraph(&gfa, &names),
        SubgraphBy::Segments => subgraph::segments_subgraph(&gfa, &names),
    };
    println!("{}", gfa_string(&new_gfa));

    Ok(())
}
