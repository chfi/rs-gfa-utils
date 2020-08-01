use std::{
    cmp::Ordering,
    fs::File,
    io::{prelude::*, BufReader},
    path::Path,
};

use bstr::{io::*, BStr, BString, ByteSlice, ByteVec};

use gfa::{
    gafpaf::{parse_gaf, CIGAROp, GAFPath, GAFStep, CIGAR},
    gfa::{Link, Orientation, Segment, GFA},
    optfields::{OptFieldVal, OptFields, OptionalFields},
};

type GAF = gfa::gafpaf::GAF<OptionalFields>;
type PAF = gfa::gafpaf::PAF<OptionalFields>;

fn get_cigar(opts: &OptionalFields) -> Option<CIGAR> {
    let cg = opts.get_field(b"cg")?;
    if let OptFieldVal::Z(cg) = &cg.value {
        let (_, cigar) = CIGAR::parse(&cg).ok()?;
        Some(cigar)
    } else {
        None
    }
}

fn get_gaf_cigar(gaf: &GAF) -> Option<CIGAR> {
    get_cigar(&gaf.optional)
}

fn get_paf_cigar(paf: &PAF) -> Option<CIGAR> {
    get_cigar(&paf.optional)
}

fn gaf_to_paf_clone(gaf: &GAF) -> PAF {
    PAF {
        query_seq_name: gaf.seq_name.clone(),
        query_seq_len: gaf.seq_len.clone(),
        query_seq_range: gaf.seq_range,
        strand: gaf.strand,
        target_seq_name: Default::default(),
        target_seq_len: gaf.path_len,
        target_seq_range: gaf.path_range,
        residue_matches: gaf.residue_matches,
        block_length: gaf.block_length,
        quality: gaf.quality,
        optional: gaf.optional.clone(),
    }
}

fn find_segment<'a, T: OptFields>(
    segs: &'a [Segment<BString, T>],
    name: &[u8],
) -> Option<&'a Segment<BString, T>> {
    let ix = segs
        .binary_search_by(|s| {
            let seg: &[u8] = s.name.as_ref();
            seg.cmp(name)
        })
        .ok()?;
    segs.get(ix)
}

fn cmp_links_find<T: OptFields, B: AsRef<[u8]>>(
    link: &Link<BString, T>,
    from: B,
    to: B,
) -> Ordering {
    let link_from: &[u8] = link.from_segment.as_ref();
    let link_to: &[u8] = link.to_segment.as_ref();
    let from_cmp = link_from.cmp(from.as_ref());
    if from_cmp == Ordering::Equal {
        link_to.cmp(to.as_ref())
    } else {
        from_cmp
    }
}

fn cmp_links<T: OptFields>(
    l1: &Link<BString, T>,
    l2: &Link<BString, T>,
) -> Ordering {
    cmp_links_find(l1, &l2.from_segment, &l2.to_segment)
}

fn find_link<'a, T: OptFields>(
    links: &'a [Link<BString, T>],
    from: &[u8],
    to: &[u8],
) -> Option<&'a Link<BString, T>> {
    let ix = links
        .binary_search_by(|l| cmp_links_find(l, from, to))
        .ok()?;
    links.get(ix)
}

fn unwrap_step(step: &GAFStep) -> (Orientation, &[u8]) {
    match step {
        GAFStep::SegId(o, id) => (*o, id.as_ref()),
        GAFStep::StableIntv(o, id, _from, _to) => (*o, id.as_ref()),
    }
}

// must take sorted segment and link slices
fn gaf_line_to_pafs<T: OptFields>(
    segments: &[Segment<BString, T>],
    links: &[Link<BString, T>],
    gaf: &GAF,
) -> Vec<PAF> {
    match &gaf.path {
        GAFPath::StableId(id) => {
            // TODO this will likely be a bit more complex, not sure
            let paf = PAF {
                target_seq_name: id.clone(),
                ..gaf_to_paf_clone(gaf)
            };
            vec![paf]
        }
        GAFPath::OrientIntv(steps) => {
            let mut pafs = Vec::new();
            let mut query_start = gaf.seq_range.0;
            let mut query_end = gaf.seq_range.1;
            // for (i, step) in steps.iter().enumerate() {
            for (i, step) in steps[1..].iter().enumerate() {
                let (o, id) = unwrap_step(step);

                let prev = steps.get(i - 1).unwrap();
                let (_, prev_id) = unwrap_step(prev);

                let seg_prev = find_segment(segments, prev_id).unwrap();
                let seg = find_segment(segments, id).unwrap();

                let link = find_link(links, id, prev_id).unwrap();

                let seg_len = seg_prev.sequence.len();

                let mut paf = gaf_to_paf_clone(gaf);
                paf.target_seq_name = id.into();
                paf.target_seq_len = seg_len;
                // TODO handle overlaps and offsets to figure out
                // target seq range
                paf.query_seq_range.0 = query_start;
                query_start = query_start + seg_len;
                paf.query_seq_range.1 = query_start; // handle case of offset

                let next_step = steps.get(i + 1);

                // query_start = gaf.
                /*

                if let Some(next) = next_step {
                    let (_, next_id) = unwrap_step(next);
                    let seg_next = find_segment(segments, next_id).unwrap();
                    let link = find_link(links, id, next_id).unwrap();
                } else {
                }
                */
            }
            pafs
        }
    }
}

pub fn gaf_to_paf<T: OptFields>(
    gfa: GFA<BString, T>,
    gaf_path: &Path,
) -> Vec<PAF> {
    let mut segments = gfa.segments;
    segments.sort_by(|s1, s2| s1.name.cmp(&s2.name));
    let mut links = gfa.links;
    links.sort_by(cmp_links);

    let file = File::open(gaf_path).unwrap();
    let lines = BufReader::new(file).byte_lines().map(|l| l.unwrap());
    let mut gafs: Vec<GAF> = Vec::new();

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
