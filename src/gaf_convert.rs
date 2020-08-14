use std::{cmp::Ordering, fs::File, io::BufReader, path::Path};

use bstr::{io::*, BString, ByteSlice, ByteVec};

use gfa::{
    gafpaf::{parse_gaf, GAFPath, GAFStep, CIGAR},
    gfa::{Link, Orientation, Segment, GFA},
    optfields::{OptFieldVal, OptFields, OptionalFields},
};

type GAF = gfa::gafpaf::GAF<OptionalFields>;
type PAF = gfa::gafpaf::PAF<OptionalFields>;

fn set_cigar(opts: &mut OptionalFields, cg: CIGAR) {
    let cg_tag = opts.iter_mut().find(|o| &o.tag == b"cg").unwrap();
    cg_tag.value = OptFieldVal::Z(cg.to_string().into());
}

fn get_cigar(opts: &OptionalFields) -> Option<CIGAR> {
    let cg = opts.get_field(b"cg")?;
    if let OptFieldVal::Z(cg) = &cg.value {
        CIGAR::from_bytes(&cg)
    } else {
        None
    }
}

fn get_gaf_cigar(gaf: &GAF) -> Option<CIGAR> {
    get_cigar(&gaf.optional)
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

fn unwrap_step(step: &GAFStep) -> (Orientation, &[u8]) {
    match step {
        GAFStep::SegId(o, id) => (*o, id.as_ref()),
        GAFStep::StableIntv(o, id, _from, _to) => (*o, id.as_ref()),
    }
}

// must take sorted segment and link slices
fn gaf_line_to_pafs<T: OptFields>(
    segments: &[Segment<BString, T>],
    gaf: &GAF,
) -> Vec<PAF> {
    match &gaf.path {
        GAFPath::StableId(id) => {
            let paf = PAF {
                target_seq_name: id.clone(),
                ..gaf_to_paf_clone(gaf)
            };
            vec![paf]
        }
        GAFPath::OrientIntv(steps) => {
            let seg_steps: Vec<(Orientation, &Segment<_, _>)> = steps
                .iter()
                .map(|s| {
                    let (o, id) = unwrap_step(s);
                    let segment = find_segment(segments, id).unwrap();
                    (o, segment)
                })
                .collect();

            let mut query_index = gaf.seq_range.0;
            let mut tgt_offset = gaf.path_range.0;
            let mut query_remaining = gaf.seq_len;

            let mut seqs: Vec<BString> = Vec::new();

            let mut pafs = Vec::new();

            let mut gaf_cigar =
                get_gaf_cigar(gaf).expect("missing cigar in GAF record");

            for (orient, target) in seg_steps {
                let seg_len = target.sequence.len();

                let step_len = query_remaining.min(seg_len - tgt_offset);
                query_remaining -= step_len;

                let query_start = query_index;
                let query_end = query_start + step_len;

                let target_seq_name = target.name.clone();
                let target_seq_len = seg_len;

                let target_seq_range = (tgt_offset, tgt_offset + step_len);

                let sequence =
                    target.sequence[tgt_offset..tgt_offset + step_len].into();

                let cg_ix = gaf_cigar.query_index(step_len);
                let split_cg = gaf_cigar.split_with_index(cg_ix);

                seqs.push(sequence);

                query_index = query_end;

                let mut optional = gaf.optional.clone();

                if split_cg.0.is_empty() {
                    set_cigar(&mut optional, split_cg.1);
                } else {
                    set_cigar(&mut optional, split_cg.0);
                    gaf_cigar = split_cg.1;
                }

                use Orientation::*;

                let strand = match (gaf.strand, orient) {
                    (Forward, Forward) => Forward,
                    (Forward, Backward) => Backward,
                    (Backward, Forward) => Backward,
                    (Backward, Backward) => Forward,
                };

                // TODO several of these fields need to be changed,
                // including strand and everything after the target
                // sequence fields
                let paf = PAF {
                    query_seq_name: gaf.seq_name.clone(),
                    query_seq_len: gaf.seq_len,
                    query_seq_range: (query_start, query_end),
                    strand,
                    target_seq_name,
                    target_seq_len,
                    target_seq_range,
                    residue_matches: gaf.residue_matches,
                    block_length: gaf.block_length,
                    quality: gaf.quality,
                    optional,
                };

                pafs.push(paf);
                tgt_offset = 0;
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

    let mut pafs: Vec<PAF> = Vec::new();

    gafs.iter().for_each(|gaf| {
        let cur_pafs = gaf_line_to_pafs(&segments, &gaf);
        pafs.extend(cur_pafs);
    });

    pafs
}
