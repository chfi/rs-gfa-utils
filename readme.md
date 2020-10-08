Rust gfautil
================

Command line tool for various operations on GFA and related files.

## Usage

Install with cargo:

```bash
cargo install gfautil
```

Or clone and build it:

```bash
git clone https://github.com/chfi/rs-gfa-utils.git
cd rs-gfa-utils
cargo build --release
```

The compiled binary will be located at `target/release/gfautil`.

```bash
$ gfautil
gfautil 0.2.0

USAGE:
    gfautil --gfa <input GFA file> <SUBCOMMAND>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information

OPTIONS:
    -i, --gfa <input GFA file>

SUBCOMMANDS:
    edge-count
    gaf2paf                      Convert a file of GAF records into PAF records
    gfa-segment-id-conversion    Convert a GFA with string names to one with integer names, and back
    gfa2vcf                      Output a VCF for the given GFA, using the graph's ultrabubbles to identify areas of
                                 variation. (experimental!)
    help                         Prints this message or the help of the given subcommand(s)
    subgraph                     Generate a subgraph of the input GFA
```


## GAF -> PAF

Given a GAF file, and the GFA used to create it, output a PAF file
derived from the GAF records. For every path segment in each GAF
record, a corresponding PAF record is produced.

Convert `example.gaf`, via `example.gfa`, with output on stdout:

```bash
gfautil --gfa ./example.gfa gaf2paf --gaf ./example.gaf
```

Save output to `out.paf`:

```bash
gfautil --gfa ./example.gfa gaf2paf --gaf ./example.gaf -o out.paf
```


## GFA -> VCF

Find the ultrabubbles in the input GFA, then use those to identify
variants. For each ultrabubble, the section covered by the bubble is
extracted from each embedded path. Those sub-paths are then compared
pairwise.

Currently the variant identification is mostly based on the nodes that
make up each path, and only barely takes the sequences into account.

Outputs is in the VCF format, on stdout.

```bash
gfautil --gfa ./example.gfa variant
```

There's a setting to skip comparing a pair of paths if their
orientations at the start and end of the bubble don't match:

```bash
gfautil --gfa ./example.gfa variant --no-inv
```


## Subgraph

Return a subgraph of the given GFA. Provide either a list of segment
names, or a list of path names. If segment names are provided, the
resulting subgraph will include the lines that contain at least one
of those segments. If path names are provided, the segments in the
given paths are used instead.

```bash
gfautil --gfa example.gfa subgraph segments --names s1 s2 s3
```

```bash
cat names.txt
s1
s2
s3
gfautil --gfa example.gfa subgraph segments --file names.txt
```

```bash
gfautil --gfa example.gfa subgraph paths --names p1 p2
```
