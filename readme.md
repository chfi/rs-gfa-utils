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

The compiled binary will be located at `target/release/gfautils`.

```bash
$ gfautil --help
gfautil 0.1.0

USAGE:
    gfautil --gfa <input GFA file> <SUBCOMMAND>

FLAGS:
    -h, --help
            Prints help information

    -V, --version
            Prints version information


OPTIONS:
    -i, --gfa <input GFA file>



SUBCOMMANDS:
    edge-count
    gaf2paf       Convert a file of GAF records into PAF records
    help          Prints this message or the help of the given subcommand(s)
    subgraph      Generate a subgraph of the input GFA
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
