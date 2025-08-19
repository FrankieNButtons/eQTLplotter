use anyhow::{Context, Result};
use plotters::coord::types::{RangedCoordf64, RangedCoordu32};
use plotters::prelude::*;
use plotters::style::{RGBAColor, RGBColor};
use rayon::prelude::*;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Create a soft, pastel-like color palette (RGB)
fn pastel_palette(n: usize) -> Vec<RGBColor> {
    (0..n)
        .map(|i| {
            let hue = (i as f64 / n as f64) * 360.0;
            hsl_to_rgb(hue, 0.55, 0.65)
        })
        .collect()
}

/// Minimal HSL -> RGB conversion returning a Plotters RGBColor
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

/// Return a nicer legend label for a path (file stem without suffix)
fn file_stem_label<P: AsRef<Path>>(p: P) -> String {
    p.as_ref()
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| p.as_ref().to_string_lossy().to_string())
}

/// Load one QTL tsv and count how many variants each phenotype (first column) has.
/// Returns a vector of per-phenotype counts.
fn load_counts_for_file(path: &str) -> Result<Vec<u32>> {
    let file = File::open(path).with_context(|| format!("Failed to open QTL file: {}", path))?;
    let reader = BufReader::new(file);
    let mut map: HashMap<String, u32> = HashMap::new();
    for line in reader.lines().filter_map(|l| l.ok()) {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut it = line.split('\t');
        if let Some(phen) = it.next() {
            *map.entry(phen.trim().to_string()).or_insert(0) += 1;
        }
    }
    if map.is_empty() {
        anyhow::bail!("No valid phenotype rows parsed from {}", path);
    }
    Ok(map.into_values().collect())
}

#[allow(dead_code)]
/// Convert a vector of per-phenotype counts into a frequency map: f[y] = number of phenotypes with y variants.
fn to_frequency(counts: &[u32]) -> (Vec<u32>, u32) {
    let &max_y = counts.iter().max().unwrap_or(&0);
    let mut freq = vec![0u32; (max_y as usize) + 1];
    for &c in counts {
        if let Some(slot) = freq.get_mut(c as usize) {
            *slot += 1;
        }
    }
    (freq, max_y)
}

/// Gaussian KDE over discrete integer support 0..=y_max for the raw counts.
/// Returns density values normalized to [0,1] so that the maximum equals 1.
fn kde_density(counts: &[u32], y_max: u32) -> Vec<f64> {
    if counts.is_empty() || y_max == 0 {
        return vec![0.0; (y_max as usize) + 1];
    }
    // compute mean and std
    let n = counts.len() as f64;
    let mean = counts.iter().map(|&c| c as f64).sum::<f64>() / n;
    let var = counts
        .iter()
        .map(|&c| {
            let d = c as f64 - mean;
            d * d
        })
        .sum::<f64>()
        / n;
    let std = var.sqrt();
    // heuristic bandwidth (avoid 0):
    let bw = std.max(1.0) * 0.6; // a bit smoother for discrete counts
    let inv_two_bw2 = 1.0 / (2.0 * bw * bw);
    // precompute Gaussian contributions per sample if y_max small; otherwise do naive sum
    let mut dens = vec![0.0f64; (y_max as usize) + 1];
    for &c in counts {
        let center = c as i64;
        // limit to 4*bw window for speed
        let w = (4.0 * bw).ceil() as i64;
        let start = (0i64).max(center - w) as usize;
        let end = (y_max as i64).min(center + w) as usize;
        for y in start..=end {
            let dy = y as f64 - c as f64;
            dens[y] += (-dy * dy * inv_two_bw2).exp();
        }
    }
    // normalize to max 1.0
    let maxv = f64::max(dens.iter().cloned().fold(0.0f64, f64::max), 1e-12);
    for v in &mut dens {
        *v /= maxv;
    }
    dens
}

/// Draw a symmetric KDE-based "violin" at x = `x_center`.
/// - `density[y]` is the KDE density at y (normalized to max 1).
/// - The half-width per y is: half_max * density[y].
fn draw_violin<DB: DrawingBackend>(
    chart: &mut ChartContext<DB, Cartesian2d<RangedCoordf64, RangedCoordu32>>,
    x_center: f64,
    density: &[f64],
    half_max: f64,
    fill: RGBAColor,
    outline: RGBColor,
) -> Result<()>
where
    DB::ErrorType: 'static,
{
    // Build left side (y ascending) and right side (y descending) for a closed polygon
    let mut left: Vec<(f64, u32)> = Vec::with_capacity(density.len());
    let mut right: Vec<(f64, u32)> = Vec::with_capacity(density.len());
    for (y, &d) in density.iter().enumerate() {
        let y_u = y as u32;
        let w = half_max * d;
        left.push((x_center - w, y_u));
        right.push((x_center + w, y_u));
    }
    right.reverse();

    let mut poly: Vec<(f64, u32)> = Vec::with_capacity(left.len() + right.len());
    poly.extend(left.into_iter());
    poly.extend(right.into_iter());

    // fill
    chart.draw_series(std::iter::once(Polygon::new(poly.clone(), fill)))?;
    // outline: draw as a polyline around
    let mut border: Vec<(f64, u32)> = poly;
    if let Some(first) = border.first().cloned() {
        border.push(first);
    }
    chart.draw_series(std::iter::once(PathElement::new(
        border,
        outline.stroke_width(1),
    )))?;
    Ok(())
}

/// Violin plot for multiple QTL files.
/// X-axis: file names (stem without suffix).
/// Y-axis: number of variants per phenotype.
/// Violin shape: smooth KDE-based density with half-width = 0.4 x-unit at max density.
pub fn run_vio_multi(
    qtl_paths_csv: &str,
    output_path: &str,
    width: u32,
    height: u32,
    threads: Option<usize>,
    y_roof: Option<u32>,
) -> Result<()> {
    if let Some(n) = threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .ok();
    }

    let paths: Vec<String> = qtl_paths_csv
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    if paths.is_empty() {
        anyhow::bail!("--multi was provided but no valid paths were found");
    }

    println!("[INFO] VIO: {} file(s)", paths.len());

    // Load per-file phenotype counts in parallel
    let per_file_counts: Vec<Vec<u32>> = paths
        .par_iter()
        .map(|p| {
            println!("[INFO] Loading {}", p);
            load_counts_for_file(p)
        })
        .collect::<Result<Vec<_>>>()?;

    // Build KDE densities per file and determine global y_max
    let mut per_file_density: Vec<Vec<f64>> = Vec::with_capacity(per_file_counts.len());

    let max_y_global = per_file_counts
        .iter()
        .filter_map(|counts| counts.iter().max())
        .max()
        .cloned()
        .unwrap_or(0);
    if max_y_global == 0 {
        anyhow::bail!("All files have zero max variants per phenotype? Check inputs.");
    }

    let y_max = y_roof.unwrap_or(max_y_global);

    for counts in &per_file_counts {
        let dens = kde_density(counts, y_max);
        per_file_density.push(dens);
    }

    let x_min = -0.6f64; // leave margin for half-width
    let x_max = (paths.len() as f64) - 0.4;

    let root = BitMapBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(60)
        .y_label_area_size(120)
        .build_cartesian_2d(x_min..x_max, 0u32..y_max)?;

    let stems: Vec<String> = paths.iter().map(|p| file_stem_label(p)).collect();
    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        // .y_desc("#variants per phenotype")
        .x_labels(stems.len())
        .x_label_formatter(&|x: &f64| {
            let i = (x + 0.5).floor() as isize;
            if i >= 0 && (i as usize) < stems.len() {
                stems[i as usize].clone()
            } else {
                String::new()
            }
        })
        .y_label_formatter(&|y: &u32| format!("{}", y))
        .label_style(("sans-serif", 50))
        .axis_desc_style(("sans-serif", 30))
        .draw()?;

    // Colors & alpha fills
    let palette = pastel_palette(stems.len().max(2));

    // Each violin centered at integer x = idx
    let half_max = 0.4f64; // half width in x-units at max density

    for (idx, dens) in per_file_density.iter().enumerate() {
        let x_center = idx as f64;
        let base = palette[idx % palette.len()];
        // Light alpha fill (0.35) and darker outline
        let fill = RGBAColor(base.0, base.1, base.2, 0.35);
        let outline = RGBColor(base.0 / 2, base.1 / 2, base.2 / 2);
        draw_violin(&mut chart, x_center, dens, half_max, fill, outline)?;
    }

    // Add a simple legend-like labels near the bottom (optional)
    // We'll skip formal legend API to keep it clean; x tick labels already show stems.

    root.present()?;

    println!("[OK] VIO: wrote {}", output_path);
    println!("[INFO] Files (left -> right):");
    for (i, p) in paths.iter().enumerate() {
        println!("  [{}] {}", i + 1, p);
    }

    Ok(())
}
