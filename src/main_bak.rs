use anyhow::{Context, Result};
use plotters::prelude::*;
use rayon::prelude::*;
use regex::Regex;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};

const REPORT_EVERY: usize = 100_000;

#[derive(Debug)]
struct Point {
    chr_num: u32,
    bp: u32,
    ind: u32,
    y: f64,   // -log10(p)
}

use plotters::style::RGBColor;

/// Generate a soft, pastel-like palette of size `n`
fn pastel_palette(n: usize) -> Vec<RGBColor> {
    (0..n).map(|i| {
        let hue = (i as f64 / n as f64) * 360.0;
        hsl_to_rgb(hue, 0.55, 0.65)
    }).collect()
}

/// Minimal HSL -> RGB conversion that returns a Plotters RGBColor
fn hsl_to_rgb(h: f64, s: f64, l: f64) -> RGBColor {
    let c = (1.0 - (2.0 * l - 1.0).abs()) * s;
    let h_prime = h / 60.0;
    let x = c * (1.0 - ((h_prime % 2.0) - 1.0).abs());
    let (r1, g1, b1) = if (0.0..1.0).contains(&h_prime) {
        (c, x, 0.0)
    } else if (1.0..2.0).contains(&h_prime) {
        (x, c, 0.0)
    } else if (2.0..3.0).contains(&h_prime) {
        (0.0, c, x)
    } else if (3.0..4.0).contains(&h_prime) {
        (0.0, x, c)
    } else if (4.0..5.0).contains(&h_prime) {
        (x, 0.0, c)
    } else {
        (c, 0.0, x)
    };
    let m = l - c / 2.0;
    let (r, g, b) = (r1 + m, g1 + m, b1 + m);
    RGBColor(
        (r.clamp(0.0, 1.0) * 255.0) as u8,
        (g.clamp(0.0, 1.0) * 255.0) as u8,
        (b.clamp(0.0, 1.0) * 255.0) as u8,
    )
}

/// Map CHR string (possibly with "chr" prefix) to numeric ordering.
/// 1..22 => 1..22, X=>23, Y=>24, M/MT=>25. Others => None (skip).
fn chr_to_num(chr_raw: &str) -> Option<u32> {
    let c = chr_raw.trim().to_ascii_lowercase();
    let c = c.strip_prefix("chr").unwrap_or(&c);
    match c {
        "x" => Some(23),
        "y" => Some(24),
        "m" | "mt" => Some(25),
        _ => c.parse::<u32>().ok().filter(|v| (1..=22).contains(v)),
    }
}

/// Parse a TSV line like your Python code:
/// names=["gene","variant","dist","MA_samples","MA_count","MAF","p","beta","se"]
/// We only need `variant` and `p`. The `variant` may look like:
///  - 1_12345, chr1_12345, 1:12345, chr1:12345, X_999, MT:123, etc.
fn parse_line(
    line: &str,
    re_variant: &Regex,
    eps: f64,
) -> Option<(u32, u32, f64)> {
    // Split into (up to) 9 fields; skip malformed lines
    let mut it = line.split('\t');
    let _gene = it.next()?;
    let variant = it.next()?.trim();
    // skip: dist, MA_samples, MA_count, MAF
    it.next()?; it.next()?; it.next()?; it.next()?;
    let p_str = it.next()?; // p
    // skip: beta, se (if present)
    // (Don’t rely on them existing)

    // Extract CHR/BP from variant via regex:
    // (?i) => case-insensitive; optional "chr" prefix; separator ":" or "_"
    // Groups: chr, bp
    if let Some(caps) = re_variant.captures(variant) {
        let chr_raw = caps.name("chr")?.as_str();
        let bp_raw = caps.name("bp")?.as_str();
        let chr_num = chr_to_num(chr_raw)?;
        let bp = bp_raw.parse::<u32>().ok()?;

        // p -> -log10(p), clamp p to [eps, +inf)
        let p = p_str.parse::<f64>().ok()?;
        let p_adj = if p.is_finite() && p > 0.0 { p } else { eps };
        let y = -p_adj.log10();
        Some((chr_num, bp, y))
    } else {
        None
    }
}

fn main() -> Result<()> {
    // Usage:
    //   cargo run --release -- [input.tsv] [output.png]
    // Defaults:
    //   input  = GeQTLResults_all.tsv.renamed.tsv
    //   output = manhattan.png
    let args: Vec<String> = env::args().collect();
    let input_path = args.get(1)
        .map(|s| s.as_str())
        .unwrap_or("GeQTLResults_all.tsv.renamed.tsv");
    let output_path = args.get(2)
        .map(|s| s.as_str())
        .unwrap_or("manhattan.png");

    println!("[INFO] Input: {input_path}");
    println!("[INFO] Output: {output_path}");

    // Prepare regex for variant extraction (case-insensitive):
    // Accepts optional "chr" prefix, then (1..22|X|Y|M|MT), then ":" or "_", then digits
    let re_variant = Regex::new(r"(?i)^(?:chr)?(?P<chr>[0-9]{1,2}|x|y|m|mt)[:_](?P<bp>\d+)")
        .expect("invalid regex");

    // We’ll parallelize the *parsing* step. Fastest way: read all lines, then par_iter.
    // For huge files, this uses memory; if you need lower memory, we can switch to a
    // streaming + per-chromosome bucketing approach.
    let file = File::open(input_path)
        .with_context(|| format!("Failed to open input file: {input_path}"))?;
    let reader = BufReader::new(file);

    // Collect all lines (skip header if present by checking if it starts with '#')
    let lines: Vec<String> = reader
        .lines()
        .filter_map(|l| l.ok())
        .filter(|l| !l.starts_with('#')) // ignore commented headers
        .collect();

    // Parallel parse & compute -log10(p)
    let eps = f64::MIN_POSITIVE; // similar to np.finfo(float).tiny
    let mut points: Vec<(u32, u32, f64)> = lines
        .par_iter()
        .filter_map(|line| parse_line(line, &re_variant, eps))
        .collect();

    if points.is_empty() {
        anyhow::bail!("No valid rows parsed. Check the input format or regex.");
    }

    // Parallel sort by (chr_num, bp)
    points.par_sort_unstable_by(|a, b| {
        match a.0.cmp(&b.0) {
            std::cmp::Ordering::Equal => a.1.cmp(&b.1),
            other => other,
        }
    });

    // Assign running index "ind" (genomic position across chromosomes)
    let mut max_y = 0.0_f64;
    let mut out = Vec::with_capacity(points.len());
    for (i, (chr_num, bp, y)) in points.into_iter().enumerate() {
        if y.is_finite() && y > max_y { max_y = y; }
        out.push(Point {
            chr_num,
            bp,
            ind: i as u32,
            y,
        });
        if i % REPORT_EVERY == 0 && i > 0 {
            println!("[INFO] Processed {i} points...");
        }
    }

    println!("[INFO] Total points parsed: {}", out.len());

    // ==== Plot with Plotters ====
    let width = 1800;
    let height = 900;
    let root = BitMapBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let x_max = out.len().saturating_sub(1) as u32;
    let y_max = (max_y.ceil() as i32 + 1).max(10) as u32; // a little headroom

    // Pastel-like colors for chromosomes 1..=25 (including X=23, Y=24, M/MT=25)
    let palette = pastel_palette(25);
    println!("[INFO] Using {} pastel colors (per chromosome)", palette.len());

    let mut chart = ChartBuilder::on(&root)
        .caption("Manhattan (multi-threaded preprocessing)", ("sans-serif", 28))
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(0u32..x_max, 0f64..y_max as f64)?;

    chart.configure_mesh()
        .x_desc("Genomic position (concatenated chromosomes)")
        .y_desc("-log10(p)")
        .light_line_style(&TRANSPARENT)
        .draw()?;

    // Draw in chunks; each point colored by its chromosome (pastel palette)
    let mut drawn: usize = 0;
    for chunk in out.chunks(REPORT_EVERY) {
        chart.draw_series(
            chunk.iter().map(|p| {
                // chr_num: 1..=25 -> index 0..=24
                let idx = (p.chr_num.saturating_sub(1) as usize) % palette.len();
                Circle::new((p.ind, p.y), 1, palette[idx].filled())
            })
        )?;
        drawn += chunk.len();
        println!("[INFO] Plotted {} / {} points", drawn, out.len());
    }

    // Optional: light genome-wide threshold line (e.g., 7.3 ~ 5e-8)
    let gw_line = 7.30103_f64; // -log10(5e-8)
    if gw_line < y_max as f64 {
        chart.draw_series(std::iter::once(PathElement::new(
            vec![(0u32, gw_line), (x_max, gw_line)],
            &BLACK.mix(0.4),
        )))?
        .label("5e-8")
        .legend(|(x, y)| PathElement::new(vec![(x - 8, y), (x + 8, y)], &BLACK));
    }

    chart.configure_series_labels().border_style(&BLACK).draw()?;
    root.present()?;

    println!("[INFO] Plot complete. Total {} points drawn.", out.len());
    println!("[INFO] Saved plot to: {output_path}");
    Ok(())
}