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
gfautil 0.3.2

USAGE:
    gfautil [FLAGS] [OPTIONS] -i <input GFA file> <SUBCOMMAND>

FLAGS:
        --debug      Show debug messages
    -h, --help       Prints help information
        --info       Show info messages
        --quiet      Show no messages
    -V, --version    Prints version information

OPTIONS:
    -i <input GFA file>
    -t, --threads <threads>    The number of threads to use when applicable. If omitted, Rayon's default will be used,
                               based on the RAYON_NUM_THREADS environment variable, or the number of logical CPUs

SUBCOMMANDS:
    edge-count
    gaf2paf         Convert a file of GAF records into PAF records
    gfa2vcf         Output a VCF for the given GFA, using the graph's ultrabubbles to identify areas of variation
    help            Prints this message or the help of the given subcommand(s)
    id-convert      Convert a GFA with string names to one with integer names, and back
    snps            Given a reference path from the GFA, by name, find and report the SNPs for all other paths
                    compared to the reference.
    subgraph        Generate a subgraph of the input GFA
    ultrabubbles
```


## GAF -> PAF

Given a GAF file, and the GFA used to create it, output a PAF file
derived from the GAF records. For every path segment in each GAF
record, a corresponding PAF record is produced.

Convert `example.gaf`, via `example.gfa`, with output on stdout:

```bash
gfautil -i ./example.gfa gaf2paf --gaf ./example.gaf
```

Save output to `out.paf`:

```bash
gfautil -i ./example.gfa gaf2paf --gaf ./example.gaf -o out.paf
```


## GFA -> VCF

Find the ultrabubbles in the input GFA, then use those to identify
variants. For each ultrabubble, the section covered by the bubble is
extracted from each embedded path. Those sub-paths are then compared
pairwise.

The `-u` option can be used to load the ultrabubbles from a file (output
by the `ultrabubbles` command) instead of computing them.

Currently the variant identification is mostly based on the nodes that
make up each path, and only barely takes the sequences into account.

Outputs is in the VCF format, on stdout.

```bash
gfautil -i ./example.gfa gfa2vcf
```

There's a setting to skip comparing a pair of paths if their
orientations at the start and end of the bubble don't match:

```bash
gfautil -i ./example.gfa gfa2vcf --no-inv
```

Loading the list of ultrabubbles from a file:
```bash
gfautil -i ./example.gfa gfa2vcf -u example.ultrabubbles
```

## Identify SNPs in GFA against reference path

Given the name of a path in the input GFA to use as reference,
identify SNPs among all other paths, using either a list of
ultrabubbles constructed using the `gfautil ultrabubbles` command, or
a list of SNP positions.

Outputs a tab-delimited list in the format:

```
<query-path-name>\t<reference base>\t<reference pos>\t<query base>\t<query pos>
```

SNP positions can be provided as a list in the arguments to `gfautil`:

```bash
gfautil --debug -t 8 -i ./example.gfa snps --ref "reference path name" --snps 1234 5677 1> example.gfa.snps
```

SNP positions can also be provided as a file, with one position per line:

```bash
gfautil --debug -t 8 -i ./example.gfa snps --ref "reference path name" --snps-file ./positions.txt 1> example.gfa.snps
```

Using ultrabubbles from a file:
```bash
gfautil -i ./example.gfa snps --ref the_path -u example.bubbles
```


## Subgraph

Return a subgraph of the given GFA. Provide either a list of segment
names, or a list of path names. If segment names are provided, the
resulting subgraph will include the lines that contain at least one
of those segments. If path names are provided, the segments in the
given paths are used instead.

```bash
gfautil -i example.gfa subgraph segments --names s1 s2 s3
```

```bash
cat names.txt
s1
s2
s3
gfautil -i example.gfa subgraph segments --file names.txt
```

```bash
gfautil -i example.gfa subgraph paths --names p1 p2
```
