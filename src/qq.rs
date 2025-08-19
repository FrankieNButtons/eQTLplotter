use anyhow::{Context, Result};
use plotters::prelude::*;
use plotters::style::{RGBAColor, RGBColor};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

/// Get a file stem (without extension) to use in labels
fn file_stem_label<P: AsRef<Path>>(p: P) -> String {
    p.as_ref()
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| p.as_ref().to_string_lossy().to_string())
}

/// Parse a CSV string that MUST contain exactly two paths separated by a single comma.
fn parse_two_paths(csv: &str) -> Result<(String, String)> {
    let comma_cnt = csv.matches(',').count();
    if comma_cnt != 1 {
        anyhow::bail!(
            "--multi must contain exactly two paths separated by a single comma (found {} comma[s])",
            comma_cnt
        );
    }
    let parts: Vec<String> = csv
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    if parts.len() != 2 {
        anyhow::bail!("--multi must contain exactly two non-empty paths");
    }
    Ok((parts[0].clone(), parts[1].clone()))
}

/// Try to detect p-value column index by header tokens; fallback to index 6 (7th col, 0-based) if no header.
fn detect_p_index_from_header(tokens: &[&str]) -> Option<usize> {
    let keys = ["p", "P", "pval", "pvalue", "p_val", "p_value"]; // be generous
    for (i, t) in tokens.iter().enumerate() {
        let tt = t.trim();
        if keys.iter().any(|k| tt.eq_ignore_ascii_case(k)) {
            return Some(i);
        }
    }
    None
}

/// Load p-values keyed by either phenotype (column 0) for grouped mode, or variant id (column 1) for genotype-level.
/// For each key we keep the MIN p-value within that file.
fn load_p_map(path: &str, grouped: bool) -> Result<HashMap<String, f64>> {
    let file = File::open(path).with_context(|| format!("Failed to open QTL file: {}", path))?;
    let mut reader = BufReader::new(file);

    let mut first_line = String::new();
    let mut has_header = false;
    let mut p_idx: usize = 6; // default to 7th column if no header
    let mut key_idx: usize = if grouped { 0 } else { 1 };

    // Peek first non-empty, non-comment line
    loop {
        first_line.clear();
        let bytes = reader.read_line(&mut first_line)?;
        if bytes == 0 {
            anyhow::bail!("{} appears empty", path);
        }
        if first_line.trim().is_empty() || first_line.starts_with('#') {
            continue;
        }
        let toks: Vec<&str> = first_line.trim_end().split('\t').collect();
        if let Some(idx) = detect_p_index_from_header(&toks) {
            has_header = true;
            p_idx = idx;
            if !grouped {
                // Detect a likely variant column by common header names
                // Fallback remains column 1 if not found
                let vnames = ["variant", "snp", "var", "snp_id", "variant_id", "rsid"];
                if let Some(vidx) = toks.iter().position(|t| {
                    let tt = t.trim();
                    vnames.iter().any(|vn| tt.eq_ignore_ascii_case(vn))
                }) {
                    key_idx = vidx;
                }
            }
        }
        break;
    }

    let mut map: HashMap<String, f64> = HashMap::new();

    // helper to process a single logical line
    let mut process_line = |line: &str| {
        if line.is_empty() || line.starts_with('#') {
            return;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < p_idx + 1 {
            return;
        }
        if cols.len() <= key_idx {
            return;
        }
        let key = cols[key_idx].trim();
        if key.is_empty() {
            return;
        }
        if let Ok(p) = cols[p_idx].trim().parse::<f64>() {
            let e = map.entry(key.to_string()).or_insert(f64::INFINITY);
            if p < *e {
                *e = p;
            }
        }
    };

    // If we had a header, skip it and process subsequent lines; otherwise the first_line is data.
    if !has_header {
        process_line(first_line.trim_end());
    }

    let mut buf = String::new();
    loop {
        buf.clear();
        let n = reader.read_line(&mut buf)?;
        if n == 0 {
            break;
        }
        process_line(buf.trim_end());
    }

    if map.is_empty() {
        anyhow::bail!("No valid p-values parsed from {}", path);
    }
    Ok(map)
}

/// Load per-row p-values by **row index** (no header). Uses column 6 (0-based) as p-value.
fn load_p_list(path: &str) -> Result<Vec<f64>> {
    const P_IDX: usize = 6; // 7th column is p-value
    let file = File::open(path).with_context(|| format!("Failed to open QTL file: {}", path))?;
    let reader = BufReader::new(file);

    let mut out: Vec<f64> = Vec::new();
    let mut raw_lines: usize = 0;
    let mut skipped_comments: usize = 0;
    let mut bad_lines: usize = 0;

    for line in reader.lines().filter_map(|l| l.ok()) {
        raw_lines += 1;
        if line.is_empty() || line.starts_with('#') {
            skipped_comments += 1;
            continue;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() <= P_IDX {
            bad_lines += 1;
            continue;
        }
        match cols[P_IDX].trim().parse::<f64>() {
            Ok(p) if p.is_finite() && p >= 0.0 => out.push(p),
            _ => {
                bad_lines += 1;
            }
        }
    }

    if out.is_empty() {
        anyhow::bail!(
            "No valid p-values parsed from {} (raw={}, comments/empty={}, bad={})",
            path,
            raw_lines,
            skipped_comments,
            bad_lines
        );
    }

    let (mut pmin, mut pmax) = (f64::INFINITY, f64::NEG_INFINITY);
    let mut zeros = 0usize;
    for &p in &out {
        if p < pmin {
            pmin = p;
        }
        if p > pmax {
            pmax = p;
        }
        if p == 0.0 {
            zeros += 1;
        }
    }
    println!(
        "[INFO] Loaded {} rows from {} (raw={}, comments/empty={}, bad={}; min_p={:.3e}, max_p={:.3e}, zeros={})",
        out.len(),
        path,
        raw_lines,
        skipped_comments,
        bad_lines,
        pmin,
        pmax,
        zeros
    );

    Ok(out)
}

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

/// Run QQ-like scatter comparing p-values from two QTL files.
/// Y-axis: p from first file; X-axis: p from second file.
/// With `grouped=true`, we collapse to phenotype-level by taking min p per phenotype in each file and join on phenotype.
/// Otherwise, we join on the variant id in column 1 (chr:pos:ref:alt as per your convention).
pub fn run_qq(
    qtl_paths_csv: &str,
    output_path: &str,
    width: u32,
    height: u32,
    threads: Option<usize>,
    grouped: bool,
    qcut_percent: Option<f64>, // drop top percent by max(x,y)
    diff: bool,                // color by above/below diagonal
    roof: Option<f64>,         // cap -log10(p) at this ceiling
) -> Result<()> {
    if let Some(n) = threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .ok();
    }

    let (path_a, path_b) = parse_two_paths(qtl_paths_csv)?; // a = Y, b = X
    println!("[INFO] QQ: first (Y) = {}", path_a);
    println!("[INFO] QQ: second (X) = {}", path_b);
    if let Some(p) = qcut_percent {
        println!("[INFO] Qcut percent: {:.3}%", p);
    }
    if diff {
        println!("[INFO] Diff-color: enabled (y>x vs y<=x)");
    }
    if let Some(r) = roof {
        println!("[INFO] Roof: cap -log10(p) at {:.3}", r);
    }

    // Build pairs (px, py) depending on grouping mode
    let mut points: Vec<(f64, f64)> = Vec::new(); // stored as (-log10(px), -log10(py))

    if grouped {
        // Grouped: take min p per phenotype in each file and join on phenotype keys
        let (map_a_res, map_b_res) =
            rayon::join(|| load_p_map(&path_a, true), || load_p_map(&path_b, true));
        let map_a = map_a_res?;
        let map_b = map_b_res?;
        println!(
            "[INFO] Loaded phenotypes: A={} B={}",
            map_a.len(),
            map_b.len()
        );

        let mut overlap = 0usize;
        for (k, &py) in &map_a {
            if let Some(&px) = map_b.get(k) {
                if py.is_finite() && px.is_finite() {
                    overlap += 1;
                }
            }
        }
        if overlap == 0 {
            anyhow::bail!("No overlapping phenotypes between files.");
        }
        println!("[INFO] Overlap phenotypes: {}", overlap);

        points.reserve(overlap);
        for (k, &py) in &map_a {
            if let Some(&px) = map_b.get(k) {
                if py.is_finite() && px.is_finite() {
                    points.push((-px.log10(), -py.log10())); // X from second, Y from first
                }
            }
        }
    } else {
        // Non-grouped: pair strictly by **row index** using 7th column as p-value
        let (list_a_res, list_b_res) =
            rayon::join(|| load_p_list(&path_a), || load_p_list(&path_b));
        let list_a = list_a_res?;
        let list_b = list_b_res?;
        let n_a = list_a.len();
        let n_b = list_b.len();
        let n = n_a.min(n_b);
        if n == 0 {
            anyhow::bail!("No overlapping rows after parsing");
        }
        if n_a != n_b {
            println!(
                "[INFO] Row count mismatch: A={} vs B={}. Using first {} rows in each for QQ.",
                n_a, n_b, n
            );
        } else {
            println!("[INFO] Row counts match: {} rows", n);
        }

        let mut nonfinite = 0usize;
        points.reserve(n);
        for i in 0..n {
            let py = list_a[i];
            let px = list_b[i];
            if py.is_finite() && px.is_finite() {
                points.push((-px.log10(), -py.log10()));
            } else {
                nonfinite += 1;
            }
        }
        if nonfinite > 0 {
            println!("[INFO] Skipped {} non-finite pairs", nonfinite);
        }
        if points.is_empty() {
            anyhow::bail!("All pairs were non-finite; nothing to plot");
        }
        println!("[INFO] QQ: {} matched pairs by row index", points.len());
    }

    // Optional: cap values at `roof`
    if let Some(r) = roof {
        for p in &mut points {
            if p.0 > r { p.0 = r; }
            if p.1 > r { p.1 = r; }
        }
    }

    // Optional: drop the top qcut% of points by max(x,y) to make the plot denser
    if let Some(pct) = qcut_percent {
        if pct.is_finite() && pct > 0.0 {
            let n = points.len();
            let k = ((pct / 100.0) * n as f64).floor() as usize;
            if k > 0 && k < n {
                let mut mags: Vec<f64> = points.iter().map(|(x, y)| x.max(*y)).collect();
                mags.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let cutoff = mags[n - k];
                let before = points.len();
                points = points
                    .into_iter()
                    .filter(|(x, y)| x.max(*y) <= cutoff)
                    .collect();
                println!(
                    "[INFO] Qcut: removed top {:.3}% ({} of {}), threshold={:.3}",
                    pct,
                    before - points.len(),
                    before,
                    cutoff
                );
            } else {
                println!("[INFO] Qcut: {:.3}% ignored (k={} not in 1..n-1).", pct, k);
            }
        }
    }

    // We already pushed (-log10(px), -log10(py)) into `points` above
    let max_x = points.iter().map(|&(x, _)| x).fold(0. / 0., f64::max);
    let max_y = points.iter().map(|&(_, y)| y).fold(0. / 0., f64::max);

    // Determine axis upper bound according to roof and/or qcut
    let mut axis_upper: f64 = max_x.max(max_y);
    // Find qcut threshold if qcut_percent is provided and resulted in any points being dropped
    let qcut = if let Some(pct) = qcut_percent {
        if pct.is_finite() && pct > 0.0 {
            // Find the cutoff used (after qcut filtering)
            let n = points.len();
            let k = ((pct / 100.0) * n as f64).floor() as usize;
            if k > 0 && k < n {
                let mut mags: Vec<f64> = points.iter().map(|(x, y)| x.max(*y)).collect();
                mags.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                Some(mags[n - k])
            } else {
                None
            }
        } else {
            None
        }
    } else {
        None
    };
    axis_upper = match (roof, qcut) {
        (Some(r), Some(q)) => r.min(q),
        (Some(r), None) => r,
        (None, Some(q)) => q,
        (None, None) => axis_upper,
    };

    let stem_a = file_stem_label(&path_a);
    let stem_b = file_stem_label(&path_b);

    let root = BitMapBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(60)
        .y_label_area_size(80)
        .build_cartesian_2d(0f64..axis_upper, 0f64..axis_upper)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .x_desc(&format!("-log10(P) ({})", stem_b)) // 修改后
        .y_desc(&format!("-log10(P) ({})", stem_a)) // 修改后
        .label_style(("sans-serif", 28))
        .axis_desc_style(("sans-serif", 36))
        .draw()?;

    if diff {
        let pal = pastel_palette(2);
        let c_over = RGBAColor(pal[0].0, pal[0].1, pal[0].2, 0.65); // y > x
        let c_under = RGBAColor(pal[1].0, pal[1].1, pal[1].2, 0.65); // y <= x

        let under: Vec<(f64, f64)> = points.iter().cloned().filter(|(x, y)| y <= x).collect();
        let over: Vec<(f64, f64)> = points.iter().cloned().filter(|(x, y)| y > x).collect();

        chart.draw_series(under.iter().map(|(x, y)| Circle::new((*x, *y), 3, c_under.filled())))?;
        chart.draw_series(over.iter().map(|(x, y)| Circle::new((*x, *y), 3, c_over.filled())))?;
    } else {
        let base = pastel_palette(3)[0];
        let fill = RGBAColor(base.0, base.1, base.2, 0.65);
        chart.draw_series(points.iter().map(|(x, y)| Circle::new((*x, *y), 3, fill.filled())))?;
    }

    // Draw a thin red dashed diagonal y = x using short segments (Plotters doesn't expose a direct dashed line builder here)
    let diag_end = axis_upper;
    let seg = (diag_end / 80.0).max(0.05); // segment length in axis units
    let mut t = 0f64;
    while t < diag_end {
        let t2 = (t + seg).min(diag_end);
        chart.draw_series(std::iter::once(PathElement::new(
            vec![(t, t), (t2, t2)],
            RED.stroke_width(2),
        )))?;
        t += seg * 2.0; // skip one segment to create a dash pattern
    }

    root.present()?;
    println!("[OK] QQ: wrote {} (points: {})", output_path, points.len());
    println!(
        "[INFO] Mode: {} | Key: {}",
        if grouped {
            "grouped (phenotype)"
        } else {
            "genotype-level"
        },
        if grouped {
            "phenotype (min p per phenotype)"
        } else {
            "variant (second column)"
        }
    );

    Ok(())
}
