# eQTLplotter

A Rust command-line toolkit for producing plots from eQTL summary statistics. It
supports Manhattan, violin, proportional bar, and QQ plots through dedicated
subcommands.

## Installation

This project uses [Cargo](https://doc.rust-lang.org/cargo/) and requires a Rust
compiler. Build the binary with:

```bash
cargo build --release
```

The compiled executable will be at `target/release/eQTLplotter`.

## General usage

Run `eQTLplotter <SUBCOMMAND> [OPTIONS]` to generate plots. Each subcommand
produces a PNG file and accepts common options such as `--width`, `--height`
and optional `--threads` to set the number of Rayon worker threads.

### Manhattan plot (`man`)

Generates Manhattan plots from a QTL results file or multiple files. Points are
shown as −log10 *p*-values along genomic coordinates.

```
eQTLplotter man --qtl results.tsv
```

Key options include:

- `--qtl <FILE>`: single QTL TSV input.
- `--multi <CSV>`: comma-separated list of QTL files. Supports three modes:
  - `sep` (default): stack plots for each file.
  - `inter`: only plot points that are significant in all files.
  - `compare`: overlay points from different files with color coding.
- `--output <PNG>`: output path; defaults to the input stem with `.png`.
- `--thresh <FLOAT>`: plot points with *p* ≤ threshold (0 disables filtering).
- `--over <FLOAT>`: only emit figures containing points with *p* ≤ value.
- `--roof <U32>`: cap the y‑axis at this −log10(*p*) value.
- `--no-tag`: hide chromosome labels.

### Violin plot (`vio`)

Displays the distribution of variants per phenotype across multiple QTL files.
Each file contributes one violin and is labeled by its stem name.

```
eQTLplotter vio --multi file1.tsv,file2.tsv
```

Options:

- `--multi <CSV>`: comma-separated list of 2+ QTL files (required).
- `--output <PNG>`: output path; defaults to a `violin_*.png` next to the first
  file.
- `--roof <U32>`: maximum y-axis value for `#variants per phenotype`.

### Proportional bar plot (`pro`)

Compares two QTL files by the proportion of significant variants per phenotype.
A histogram shows the percentage leaning toward the second file.

```
eQTLplotter pro --multi left.tsv,right.tsv
```

Options:

- `--multi <FILE_A,FILE_B>`: exactly two QTL files (required).
- `--thresh <FLOAT>`: significance cutoff for counting variants (default
  `5e-8`).
- `--bin <PERCENT>`: bin width for the histogram of percentages.
- `--diff <RATIO>`: color tails for values below `(1-RATIO)*100%` or above`RATIO*100%`.
- `--output <FILE>`: output file name (default `out.png`).
- `--venn`: plot a Venn diagram in the middle of the Histogram.(Not Fully Implemented)


### QQ plot (`qq`)

Creates a QQ scatter comparing −log10 *p*-values from two QTL files.

```
eQTLplotter qq --multi a.tsv,b.tsv
```

Options:

- `--multi <FILE_A,FILE_B>`: exactly two QTL files (first is y‑axis).
- `--group`: group by phenotype (column 0) using minimum *p* per phenotype;
  otherwise pairs rows by variant ID (column 1).
- `--qcut <PERCENT>`: drop top percentage of points by max(−log10 *p*).
- `--diff`: color points above vs. below the diagonal differently.
- `--roof <FLOAT>`: cap −log10 *p* before plotting.
- `--output <PNG>`: output path; defaults to a `qq_*.png` next to the first

## Examples

```bash
# Manhattan plot for a single result set
cargo run -- man --qtl example.tsv

# Overlay comparison of two result files
cargo run -- man --multi a.tsv,b.tsv --mode compare --thresh 5e-8

# Violin plot for multiple cohorts
cargo run -- vio --multi cohort1.tsv,cohort2.tsv,cohort3.tsv

# Proportional bar plot comparing datasets A and B
cargo run -- pro --multi A.tsv,B.tsv --thresh 1e-6

# QQ plot by phenotype grouping
cargo run -- qq --multi expr.tsv,prot.tsv --group
```

## License

This project is distributed under the terms of the MIT license.
