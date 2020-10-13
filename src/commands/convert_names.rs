use structopt::StructOpt;

use bstr::BString;
use std::{fs::File, io::Write, path::PathBuf};

use gfa::{
    gfa::{name_conversion::NameMap, GFA},
    optfields::OptionalFields,
    writer::write_gfa,
};

use super::{load_gfa, Result};

#[derive(StructOpt, Debug)]
/// Convert a GFA with string names to one with integer names, and
/// back.
pub struct GfaIdConvertArgs {
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

fn segment_id_to_usize(
    gfa_path: &PathBuf,
    gfa: &GFA<BString, OptionalFields>,
    args: &GfaIdConvertArgs,
) -> Result<()> {
    let name_map = if let Some(ref path) = &args.name_map_path {
        let map = NameMap::load_json(&path)?;
        map
    } else {
        NameMap::build_from_gfa(&gfa)
    };

    if let Some(new_gfa) = name_map.gfa_bstring_to_usize(&gfa, args.check_hash)
    {
        let new_gfa_path = converted_gfa_path(&gfa_path);
        let mut new_gfa_file = File::create(new_gfa_path.clone())?;
        let mut gfa_str = String::new();
        write_gfa(&new_gfa, &mut gfa_str);
        writeln!(new_gfa_file, "{}", gfa_str)?;
        println!("Saved converted GFA to {}", new_gfa_path.display());

        if args.name_map_path.is_none() {
            let name_map_path = gfa_to_name_map_path(&gfa_path);
            name_map.save_json(&name_map_path)?;
            println!("Saved new name map to {}", name_map_path.display());
        }
    } else {
        println!("Could not convert the GFA segment IDs");
    }

    Ok(())
}

fn segment_id_to_bstring(
    gfa_path: &PathBuf,
    gfa: &GFA<usize, OptionalFields>,
    args: &GfaIdConvertArgs,
) -> Result<()> {
    let name_map_path = args
        .name_map_path
        .as_ref()
        .expect("Need name map to convert back");
    let name_map = NameMap::load_json(&name_map_path)?;

    let new_gfa: GFA<BString, OptionalFields> = name_map
        .gfa_usize_to_bstring(&gfa)
        .expect("Error during conversion -- is it the right name map?");

    let new_gfa_path = restored_gfa_path(gfa_path);
    let mut new_gfa_file = File::create(new_gfa_path.clone())?;
    let mut gfa_str = String::new();
    write_gfa(&new_gfa, &mut gfa_str);
    writeln!(new_gfa_file, "{}", gfa_str)?;
    println!("Saved restored GFA to {}", new_gfa_path.display());

    Ok(())
}

pub fn convert_segment_ids(
    gfa_path: &PathBuf,
    args: &GfaIdConvertArgs,
) -> Result<()> {
    if !args.to_usize && args.name_map_path.is_none() {
        eprintln!("this shouldn't happen");
    }

    if args.to_usize {
        let gfa: GFA<BString, OptionalFields> = load_gfa(&gfa_path)?;
        segment_id_to_usize(&gfa_path, &gfa, args)
    } else {
        // Converting from integer to string names
        let gfa: GFA<usize, OptionalFields> = load_gfa(&gfa_path)?;
        segment_id_to_bstring(&gfa_path, &gfa, args)
    }
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
