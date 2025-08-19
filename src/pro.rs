use anyhow::{Context, Result};
use plotters::prelude::*;
use plotters::style::text_anchor::{HPos, Pos, VPos};
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
    let keys = ["p", "P", "pval", "pvalue", "p_val", "p_value"];
    for (i, t) in tokens.iter().enumerate() {
        let tt = t.trim();
        if keys.iter().any(|k| tt.eq_ignore_ascii_case(k)) {
            return Some(i);
        }
    }
    None
}

/// Collect all phenotype names (first column) without p-value thresholding.
fn collect_all_phenos(path: &str) -> Result<std::collections::HashSet<String>> {
    use std::collections::HashSet;
    let file = File::open(path).with_context(|| format!("Failed to open QTL file: {}", path))?;
    let mut reader = BufReader::new(file);
    let mut first_line = String::new();
    let mut has_header = false;

    // find first non-empty, non-comment line and decide if it's a header using p-column detection
    loop {
        first_line.clear();
        let bytes = reader.read_line(&mut first_line)?;
        if bytes == 0 { anyhow::bail!("{} appears empty", path); }
        if first_line.trim().is_empty() || first_line.starts_with('#') { continue; }
        let toks: Vec<&str> = first_line.trim_end().split('\t').collect();
        if detect_p_index_from_header(&toks).is_some() { has_header = true; }
        break;
    }

    let mut set: HashSet<String> = HashSet::new();
    let mut process_line = |line: &str| {
        if line.is_empty() || line.starts_with('#') { return; }
        let cols: Vec<&str> = line.split('\t').collect();
        if let Some(pheno) = cols.get(0).map(|s| s.trim()) {
            if !pheno.is_empty() { set.insert(pheno.to_string()); }
        }
    };

    if !has_header { process_line(first_line.trim_end()); }

    let mut buf = String::new();
    loop {
        buf.clear();
        let n = reader.read_line(&mut buf)?;
        if n == 0 { break; }
        process_line(buf.trim_end());
    }
    Ok(set)
}

/// Count significant variants PER PHENOTYPE in a QTL tsv file.
/// - phenotype column is assumed to be column 0 (first column)
/// - p-value column is auto-detected by header; if no header, defaults to column 6 (0-based)
fn count_sig_per_pheno(path: &str, p_thresh: f64) -> Result<HashMap<String, usize>> {
    let file = File::open(path).with_context(|| format!("Failed to open QTL file: {}", path))?;
    let mut reader = BufReader::new(file);

    let mut first_line = String::new();
    let mut has_header = false;
    let mut p_idx: usize = 6; // default to 7th column if no header

    // find first non-empty, non-comment line and determine p_idx
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
        }
        break;
    }

    let mut map: HashMap<String, usize> = HashMap::new();

    let mut process_line = |line: &str| {
        if line.is_empty() || line.starts_with('#') {
            return;
        }
        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() <= p_idx {
            return;
        }
        let pheno = cols.get(0).map(|s| s.trim()).unwrap_or("");
        if pheno.is_empty() {
            return;
        }
        if let Ok(p) = cols[p_idx].trim().parse::<f64>() {
            if p.is_finite() && p >= 0.0 && p <= p_thresh {
                *map.entry(pheno.to_string()).or_insert(0) += 1;
            }
        }
    };

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
        anyhow::bail!(
            "No significant variants (p <= {}) parsed from {}",
            p_thresh,
            path
        );
    }
    Ok(map)
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

/// Run proportional bar plot.
///
/// X-axis: 0..100 where 0 == "100% from file A", 50 == balanced, 100 == "100% from file B".
/// For each phenotype, we compute counts of significant variants in A and B (p <= thresh).
/// The phenotype contributes one vote to bin = round((B / (A+B)) * 100). Y-axis counts phenotypes per bin.
///
/// Arguments:
/// - `qtl_paths_csv`: exactly two file paths separated by ','
/// - `output_path`: png output
/// - `width`, `height`: canvas size
/// - `threads`: optional rayon worker count
/// - `thresh`: significance threshold in **p-value** (default suggested: 5e-8)
/// - `bin_width_percent`: width of bins in percent (e.g., 5.0 gives 21 bins at 0,5,...,100)
pub fn run_pro(
    qtl_paths_csv: &str,
    output_path: &str,
    width: u32,
    height: u32,
    threads: Option<usize>,
    thresh: f64,
    bin_width_percent: Option<f64>,
    diff: Option<f64>,
) -> Result<()> {
    if let Some(n) = threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .ok();
    }

    let (path_a, path_b) = parse_two_paths(qtl_paths_csv)?; // A on the LEFT, B on the RIGHT
    let stem_a = file_stem_label(&path_a);
    let stem_b = file_stem_label(&path_b);
    println!("[INFO] PRO: left = {} | right = {}", path_a, path_b);
    println!("[INFO] Threshold (p): {:.3e}", thresh);

    // Count significant variants per phenotype in parallel (one file per thread group)
    let (map_a_res, map_b_res) = rayon::join(
        || count_sig_per_pheno(&path_a, thresh),
        || count_sig_per_pheno(&path_b, thresh),
    );
    let map_a = map_a_res?;
    let map_b = map_b_res?;

    // Unthreshed phenotype counts
    let set_a = collect_all_phenos(&path_a)?;
    let set_b = collect_all_phenos(&path_b)?;
    let union_total = set_a.union(&set_b).count();
    let paired_total = set_a.intersection(&set_b).count();
    println!(
        "[INFO] Phenotypes (unthreshed): left={}, right={}, union={}, paired={}",
        set_a.len(), set_b.len(), union_total, paired_total
    );

    // Build proportions per phenotype among the UNION of phenotypes that have any significant hits
    let mut proportions: Vec<f64> = Vec::new(); // values in 0..=100 (percentage leaning to B)
    let mut has_any = 0usize;

    let mut keys: Vec<&String> = map_a.keys().collect();
    for k in map_b.keys() {
        if !map_a.contains_key(k) {
            keys.push(k);
        }
    }

    proportions.reserve(keys.len());
    for k in keys {
        let a = *map_a.get(k).unwrap_or(&0) as f64;
        let b = *map_b.get(k).unwrap_or(&0) as f64;
        let total = a + b;
        if total > 0.0 {
            has_any += 1;
            let toward_b = (b / total) * 100.0; // 0 == all A; 50 == balanced; 100 == all B
            proportions.push(toward_b);
        }
    }

    if has_any == 0 {
        anyhow::bail!(
            "No phenotypes with significant variants at p <= {} in either file.",
            thresh
        );
    }
    println!("[INFO] Phenotypes with any significant hits: {}", has_any);

    // Bin the proportions along 0..100 using histogram-style binning
    let bw = bin_width_percent.unwrap_or(5.0).max(0.5).min(50.0);
    let n_bins = ((100.0 / bw).ceil() as usize).max(1);

    let mut bins: Vec<usize> = vec![0; n_bins];

    for v in proportions {
        let mut idx = (v / bw).floor() as isize; // histogram-style: floor into [i*bw, (i+1)*bw)
        if idx < 0 { idx = 0; }
        if idx as usize >= n_bins { idx = (n_bins - 1) as isize; }
        bins[idx as usize] += 1;
    }

    let plotted_count: usize = bins.iter().sum();

    // Prepare plot
    let root = BitMapBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let x_max = 100.0f64;
    let y_max = bins.iter().copied().max().unwrap_or(0) as f64 * 1.1 + 1.0;

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(150)
        .y_label_area_size(240)
        .build_cartesian_2d(0f64..x_max, 0f64..y_max)?;

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .x_desc(&format!(
            "proportion of Significant varient in phenotype({} vs. {})",
            stem_a, stem_b
        ))
        .y_desc("Phenotype Count")
        .x_labels(0) // we'll draw custom labels
        .axis_desc_style(("sans-serif", 70))
        .label_style(("sans-serif", 70))
        .draw()?;

    // Custom x-axis labels: left-aligned at start, centered at middle, right-aligned at end
    let (xr, yr) = chart.plotting_area().get_pixel_range();
    let left_x = xr.start;
    let right_x = xr.end;
    let bottom_y = yr.end;
    let xlab_style = ("sans-serif", 70)
        .into_text_style(&root)
        .color(&RGBColor(40, 40, 40));

    // small vertical padding below axis line
    let pad = 6;
    // Left corner label: left-aligned
    root.draw(&Text::new(
        format!("100% {}", stem_a),
        (left_x, bottom_y + pad),
        xlab_style.clone().pos(Pos::new(HPos::Left, VPos::Top)),
    ))?;
    // Middle label: centered
    root.draw(&Text::new(
        "50%".to_string(),
        ((left_x + right_x) / 2, bottom_y + pad),
        xlab_style.clone().pos(Pos::new(HPos::Center, VPos::Top)),
    ))?;
    // Right corner label: right-aligned
    root.draw(&Text::new(
        format!("100% {}", stem_b),
        (right_x, bottom_y + pad),
        xlab_style.clone().pos(Pos::new(HPos::Right, VPos::Top)),
    ))?;

    // Colors (match vio palette: two sides + a neutral grey for the exact middle bin)
    let pal = pastel_palette(2);
    let left_tail = pal[0];
    let right_tail = pal[1];
    let left_color = RGBAColor(left_tail.0, left_tail.1, left_tail.2, 0.85);
    let right_color = RGBAColor(right_tail.0, right_tail.1, right_tail.2, 0.85);
    let mid_gray = RGBAColor(120, 120, 120, 0.85);

    // Draw bars using histogram-style bin edges
    for i in 0..n_bins {
        let count = bins[i] as f64;
        if count <= 0.0 { continue; }
        let x0 = (i as f64) * bw;
        let x1 = ((i + 1) as f64 * bw).min(100.0);

        // diff band coloring
        let d = diff.unwrap_or(0.9).clamp(0.0, 1.0);
        let lower = (1.0 - d) * 100.0;
        let upper = d * 100.0;
        // use bin midpoint for coloring decision
        let center = (x0 + x1) / 2.0;
        let fill = if center <= lower { left_color } else if center >= upper { right_color } else { mid_gray };

        chart.draw_series(std::iter::once(Rectangle::new(
            [(x0, 0.0), (x1, count)],
            fill.filled(),
        )))?;

        // dark grey numeric label above bar
        let label_y = (count + y_max * 0.02).min(y_max);
        let num_style = ("sans-serif", 70)
            .into_text_style(&root)
            .color(&RGBColor(70, 70, 70))
            .pos(Pos::new(HPos::Center, VPos::Bottom));
        chart.draw_series(std::iter::once(Text::new(
            format!("{}", bins[i]),
            ((x0 + x1) / 2.0, label_y),
            num_style.clone(),
        )))?;
    }

    // Top-right statistics: number of plotted phenotypes and its proportion relative to all with any hits
    let stat_style = ("sans-serif", 70)
        .into_text_style(&root)
        .color(&RGBColor(50, 50, 50));
    let denominator = paired_total.max(1) as f64;
    let proportion = (plotted_count as f64) * 100.0 / denominator;
    let stat_text = format!(
        "{} phenotypes plotted({:.1}% of all {})",
        plotted_count,
        proportion,
        paired_total
    );
    // place near top-right; adjust x so it doesn't go off-canvas
    let tr_x = (width as i32).saturating_sub(135 + 960);
    root.draw(&Text::new(stat_text, (tr_x, 70), stat_style))?;

    let legend_text_style = ("sans-serif", 70).into_text_style(&root)
        .color(&BLACK)
        .pos(Pos::new(HPos::Right, VPos::Center));
    let right_inset = 20; // stick to right edge
    let box_size = 45;
    let gap = 10;
    let box_left_x = width as i32 - right_inset - box_size; // box hugs right edge
    let text_x = box_left_x - gap; // text sits to the left of the box

    // Left file legend row (move down to avoid overlap with stats)
    let y1 = 170;
    root.draw(&Text::new(stem_a.clone(), (text_x, y1), legend_text_style.clone()))?;
    root.draw(&Rectangle::new([(box_left_x, y1 - box_size/2), (box_left_x + box_size, y1 + box_size/2)], left_color.filled()))?;

    // Right file legend row (move further down)
    let y2 = 220;
    root.draw(&Text::new(stem_b.clone(), (text_x, y2), legend_text_style.clone()))?;
    root.draw(&Rectangle::new([(box_left_x, y2 - box_size/2), (box_left_x + box_size, y2 + box_size/2)], right_color.filled()))?;

    root.present()?;
    println!(
        "[OK] PRO: wrote {} (phenos with hits: {}, bin width = {}%)",
        output_path, has_any, bw
    );

    Ok(())
}
