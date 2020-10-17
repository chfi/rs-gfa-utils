use saboten::{
    biedgedgraph::*,
    cactusgraph,
    cactusgraph::{BridgeForest, CactusGraph, CactusTree},
};

use bstr::{io::*, ByteSlice};
use std::{
    fs::File,
    io::BufReader,
    path::{Path, PathBuf},
};

use gfa::{
    gfa::GFA,
    parser::{GFAParser, GFAParserBuilder},
};

#[allow(unused_imports)]
use log::{debug, info, log_enabled, warn};

use super::Result;

pub fn run_saboten(gfa_path: &PathBuf) -> Result<()> {
    let ultrabubbles = find_ultrabubbles(gfa_path)?;
    print_ultrabubbles(ultrabubbles.iter())
}

pub fn print_ultrabubbles<'a, I>(ultrabubbles: I) -> Result<()>
where
    I: Iterator<Item = &'a (u64, u64)> + 'a,
{
    for (x, y) in ultrabubbles {
        println!("{}\t{}", x, y);
    }

    Ok(())
}

pub fn find_ultrabubbles(gfa_path: &PathBuf) -> Result<Vec<(u64, u64)>> {
    let mut parser_builder = GFAParserBuilder::all();
    parser_builder.paths = false;
    parser_builder.containments = false;
    let parser: GFAParser<usize, ()> = parser_builder.build();
    let gfa: GFA<usize, ()> = parser.parse_file(gfa_path)?;

    info!("Computing ultrabubbles");
    debug!("Building biedged graph");

    let be_graph = BiedgedGraph::from_gfa(&gfa);

    debug!("Building cactus graph");
    let cactus_graph = CactusGraph::from_biedged_graph(&be_graph);

    debug!("Building cactus tree");
    let cactus_tree = CactusTree::from_cactus_graph(&cactus_graph);

    debug!("Building bridge forest");
    let bridge_forest = BridgeForest::from_cactus_graph(&cactus_graph);

    debug!("Finding ultrabubbles");
    let ultrabubbles =
        cactusgraph::find_ultrabubbles(&cactus_tree, &bridge_forest);

    let ultrabubbles = cactusgraph::inverse_map_ultrabubbles(ultrabubbles);

    debug!("Done computing ultrabubbles");
    Ok(ultrabubbles.into_iter().map(|(x_y, _cont)| x_y).collect())
}

static LINE_ERROR: &str = "Ultrabubble record was missing fields";

pub fn load_ultrabubbles<P: AsRef<Path>>(path: P) -> Result<Vec<(u64, u64)>> {
    info!("Loading ultrabubbles from file {}", path.as_ref().display());
    let file = File::open(path.as_ref())?;
    let reader = BufReader::new(file);
    let lines = reader.byte_lines();

    let mut ultrabubbles = Vec::new();

    for line in lines {
        let line = line?;
        let mut fields = line.split_str("\t");
        let start = fields.next().ok_or(LINE_ERROR)?.to_str()?;
        let start = start.parse::<u64>()?;

        let end = fields.next().ok_or(LINE_ERROR)?.to_str()?;
        let end = end.parse::<u64>()?;

        ultrabubbles.push((start, end));
    }

    Ok(ultrabubbles)
}
