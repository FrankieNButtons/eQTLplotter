use anyhow::{Context, Result};
use plotters::chart::SeriesLabelPosition;
use plotters::prelude::*;
use plotters::style::text_anchor::{HPos, Pos, VPos};
use plotters::style::{FontStyle, RGBColor};
use rayon::prelude::*;
use regex::Regex;
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};

/// Compare multiple QTL files by group (e.g., gene), plotting OVERLAY per group:
/// points from different files are drawn on the same axes, with color encoding the file index.
/// `--thresh` filters points; `--over` filters which groups (figures) are emitted.
pub fn run_compare_overlay(
    qtl_paths_csv: &str,
    out_dir: &str,
    width: u32,
    height: u32,
    threads: Option<usize>,
    roof: Option<u32>,
    thresh_p: Option<f64>,
    over_p: Option<f64>,
    no_tag: bool,
    sv_thresh: Option<usize>,
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
    if paths.len() < 2 {
        anyhow::bail!("Overlay compare requires at least two input files.");
    }

    let re_variant = Regex::new(r"(?i)^(?:chr)?(?P<chr>[0-9]{1,2}|x|y|m|mt)[:_](?P<bp>\d+)")
        .expect("invalid regex");
    let eps = f64::MIN_POSITIVE;
    let y_cut = thresh_p
        .filter(|tp| tp.is_finite() && *tp > 0.0)
        .map(|tp| -tp.log10());
    let over_cut = over_p
        .filter(|ov| ov.is_finite() && *ov > 0.0)
        .map(|ov| -ov.log10());
    let sv_active = sv_thresh.unwrap_or(0) > 0;
    let sv_cut = sv_thresh.unwrap_or(0);

    // Load: for each file -> HashMap<group, Vec<(chr_num, bp, y)>>, applying point-level thresh if provided
    let mut per_file: Vec<HashMap<String, Vec<(u32, u32, f64)>>> = Vec::with_capacity(paths.len());
    let mut per_file_sv: Vec<HashMap<String, Vec<(u32, u32, f64)>>> = Vec::with_capacity(paths.len());
    for (i, path) in paths.iter().enumerate() {
        println!("[INFO] [File {}/{}] Loading {}", i + 1, paths.len(), path);
        let file =
            File::open(path).with_context(|| format!("Failed to open QTL file: {}", path))?;
        let reader = BufReader::new(file);
        let mut group_map: HashMap<String, Vec<(u32, u32, f64)>> = HashMap::new();
        let mut sv_group_map: HashMap<String, Vec<(u32, u32, f64)>> = HashMap::new();
        for line in reader.lines().filter_map(|l| l.ok()) {
            if line.starts_with('#') {
                continue;
            }
            let mut it = line.split('\t');
            let group = match it.next() {
                Some(s) => s.trim().to_string(),
                _none => continue,
            };
            let rest_joined = it.collect::<Vec<&str>>().join("\t");
            let pseudo = format!("{}\t{}", &group, rest_joined);
            if let Some((chr_num, bp, y)) = parse_line(&pseudo, &re_variant, eps) {
                if let Some(cut) = y_cut {
                    if y < cut { continue; }
                }
                let is_sv = sv_active && is_sv_from_line(&pseudo, sv_cut);
                if is_sv {
                    sv_group_map.entry(group.clone()).or_default().push((chr_num, bp, y));
                } else {
                    group_map.entry(group.clone()).or_default().push((chr_num, bp, y));
                }
            }
        }
        println!("[INFO] File {}: {} groups loaded.", path, group_map.len());
        per_file.push(group_map);
        per_file_sv.push(sv_group_map);
    }

    // Gather all groups
    let mut all_groups: BTreeSet<String> = BTreeSet::new();
    for gmap in &per_file {
        for g in gmap.keys() {
            all_groups.insert(g.clone());
        }
    }
    if all_groups.is_empty() {
        anyhow::bail!("No groups found in any input file.");
    }

    let n_files = paths.len();
    let file_palette = pastel_palette(n_files.max(2)); // color per file

    let mut n_written = 0usize;
    'group_loop: for group in &all_groups {
        // Merge keys across files for this group (include normal and SV to ensure mapping exists)
        let mut key_set: BTreeSet<(u32, u32)> = BTreeSet::new();
        for gmap in &per_file {
            if let Some(v) = gmap.get(group) {
                for &(c, b, _) in v {
                    key_set.insert((c, b));
                }
            }
        }
        for gmap_sv in &per_file_sv {
            if let Some(v) = gmap_sv.get(group) {
                for &(c, b, _) in v {
                    key_set.insert((c, b));
                }
            }
        }
        if key_set.is_empty() {
            // Nothing to draw for this group
            continue 'group_loop;
        }
        let keys_sorted: Vec<(u32, u32)> = key_set.iter().cloned().collect();

        // Determine chromosome coverage and bp range
        let mut chr_set: BTreeSet<u32> = BTreeSet::new();
        let mut min_bp: u32 = u32::MAX;
        let mut max_bp: u32 = 0;
        for (c, b) in &key_set {
            chr_set.insert(*c);
            if *b < min_bp {
                min_bp = *b;
            }
            if *b > max_bp {
                max_bp = *b;
            }
        }
        let single_chr = chr_set.len() == 1;

        // Build key -> x coordinate map
        let mut key2x: HashMap<(u32, u32), u32> = HashMap::with_capacity(keys_sorted.len());
        if single_chr {
            // Use true genomic position as x
            for k in keys_sorted.iter() {
                key2x.insert(*k, k.1);
            }
        } else {
            // Fallback to unified index across chromosomes
            for (i, k) in keys_sorted.iter().enumerate() {
                key2x.insert(*k, i as u32);
            }
        }

        // Figure-level OVER filter: does any file have any point with y >= over_cut?
        if let Some(oc) = over_cut {
            let mut has_over = false;
            'scan: for gmap in &per_file {
                if let Some(v) = gmap.get(group) {
                    for &(_, _, y) in v {
                        if y >= oc {
                            has_over = true;
                            break 'scan;
                        }
                    }
                }
            }
            if !has_over {
                println!(
                    "[SKIP] Group '{}' skipped (no points with p ≤ {} across files)",
                    group,
                    over_p.unwrap()
                );
                continue 'group_loop;
            }
        }

        // Build per-file point vectors mapped to x coordinate, including SVs
        let mut per_file_points: Vec<Vec<Point>> = Vec::with_capacity(n_files);
        let mut per_file_sv_points: Vec<Vec<Point>> = Vec::with_capacity(n_files);
        let mut global_chr_bounds: BTreeMap<u32, (u32, u32)> = BTreeMap::new();
        let mut max_y_all = 0.0f64;
        for (fi, gmap) in per_file.iter().enumerate() {
            let mut v: Vec<Point> = Vec::new();
            let mut v_sv: Vec<Point> = Vec::new();
            if let Some(raw) = gmap.get(group) {
                for &(c, b, y) in raw {
                    if let Some(&x) = key2x.get(&(c, b)) {
                        if y.is_finite() && y > max_y_all { max_y_all = y; }
                        v.push(Point { chr_num: c, _bp: b, ind: x, y });
                        if !single_chr {
                            global_chr_bounds.entry(c).and_modify(|e| { if x > e.1 { e.1 = x; } if x < e.0 { e.0 = x; } }).or_insert((x, x));
                        }
                    } else {
                        // Key missing from merged set; skip to avoid panic
                        // (This can occur if filters differ across files.)
                        continue;
                    }
                }
            }
            if let Some(raw_sv) = per_file_sv[fi].get(group) {
                for &(c, b, y) in raw_sv {
                    if let Some(&x) = key2x.get(&(c, b)) {
                        v_sv.push(Point { chr_num: c, _bp: b, ind: x, y });
                    } else {
                        // Key missing from merged set; skip gracefully
                        continue;
                    }
                }
            }
            v.sort_unstable_by_key(|p| p.ind);
            v_sv.sort_unstable_by_key(|p| p.ind);
            per_file_points.push(v);
            per_file_sv_points.push(v_sv);
            println!(
                "[INFO] Group '{}' file {}: {} points ({} SV)",
                group,
                fi + 1,
                per_file_points[fi].len(),
                per_file_sv_points[fi].len()
            );
        }

        let (x_min, x_max) = if single_chr {
            (min_bp, max_bp)
        } else {
            (0u32, keys_sorted.len().saturating_sub(1) as u32)
        };
        let y_max: u32 = roof.unwrap_or(60u32);
        let mut y_min: f64 = y_cut.unwrap_or(0.0);
        if y_min >= y_max as f64 {
            y_min = (y_max as f64 - 1.0).max(0.0);
        }

        // In compare mode, set x-axis label independently
        let chr_label_for_axis = if single_chr {
            let only_chr = *chr_set.iter().next().unwrap();
            match only_chr {
                23 => "chrX".to_string(),
                24 => "chrY".to_string(),
                25 => "chrMT".to_string(),
                n => format!("chr{}", n),
            }
        } else {
            "chr*".to_string()
        };
        let x_label = format!("position({} on {})", group, chr_label_for_axis);

        let fname = format!("{}/{}.png", out_dir, sanitize_file_name(group));
        let root = BitMapBackend::new(&fname, (width, height)).into_drawing_area();
        root.fill(&WHITE)?;
        let mut chart = ChartBuilder::on(&root)
            .margin(10)
            .x_label_area_size(40)
            .y_label_area_size(60)
            .build_cartesian_2d(x_min..x_max, y_min..y_max as f64)?;

        chart
            .configure_mesh()
            .x_desc(x_label)
            .y_desc("-log10(p)")
            .disable_x_mesh()
            .disable_y_mesh()
            .draw()?;

        // Draw GW line first (lowest layer), thin & dashed
        let gw_line = 7.30103_f64;
        let gw_draw = gw_line.max(y_min);
        if gw_draw <= y_max as f64 {
            // dashed by drawing short segments with gaps
            let seg = (x_max.saturating_sub(x_min) / 80).max(2); // segment length
            let mut x = x_min;
            while x < x_max {
                let x2 = x.saturating_add(seg).min(x_max);
                chart.draw_series(std::iter::once(PathElement::new(
                    vec![(x, gw_draw), (x2, gw_draw)],
                    RED.stroke_width(1),
                )))?;
                x = x2 + seg; // skip a segment to create dash
            }
        }

        let label_y = (y_min + 0.8).min(y_max as f64);
        if y_min > 0.0 {
            chart.draw_series(std::iter::once(PathElement::new(
                vec![(x_min, y_min), (x_max, y_min)],
                BLACK.stroke_width(3),
            )))?;
        }

        // Draw points per file with its own color and register legend label for each file
        for (fi, pts) in per_file_points.iter().enumerate() {
            let color = file_palette[fi % file_palette.len()];
            chart
                .draw_series(
                    pts.iter()
                        .map(|p| Circle::new((p.ind, p.y.min(y_max as f64)), 1, color.filled())),
                )?
                .label(
                    std::path::Path::new(&paths[fi])
                        .file_stem()
                        .unwrap_or_default()
                        .to_string_lossy()
                        .to_string(),
                )
                .legend(move |(x, y)| Circle::new((x, y), 5, color.filled()));
        }
        // Overlay SV points and legends per file
        let sv_color = RGBColor(255, 0, 255);
        for (fi, pts_sv) in per_file_sv_points.iter().enumerate() {
            if !pts_sv.is_empty() {
                chart
                    .draw_series(pts_sv.iter().map(|p| Circle::new((p.ind, p.y.min(y_max as f64)), 1, sv_color.filled())))?
                    .label(
                        std::path::Path::new(&paths[fi])
                            .file_stem()
                            .unwrap_or_default()
                            .to_string_lossy()
                            .to_string() + "(SV)"
                    )
                    .legend(move |(x, y)| Circle::new((x, y), 5, sv_color.filled()));
            }
        }

        // GW line already drawn above as dashed, thin, and lowest layer.

        if !single_chr {
            // Chromosome borders and labels based on global bounds
            let mut starts: Vec<u32> = global_chr_bounds.values().map(|(s, _e)| *s).collect();
            starts.sort_unstable();
            for start in starts.iter().skip(1) {
                chart.draw_series(std::iter::once(PathElement::new(
                    vec![(*start, 0f64), (*start, y_max as f64)],
                    RGBColor(210, 210, 210).stroke_width(1),
                )))?;
            }
            for (&chr_num, &(start, end)) in &global_chr_bounds {
                let mid = (start + end) / 2;
                let chr_label = match chr_num {
                    23 => "X".to_string(),
                    24 => "Y".to_string(),
                    25 => "MT".to_string(),
                    n => n.to_string(),
                };

                // Only draw chromosome labels if no_tag is false
                if !no_tag {
                    chart.draw_series(std::iter::once(Text::new(
                        chr_label,
                        (mid, label_y),
                        ("sans-serif", 50)
                            .into_font()
                            .style(FontStyle::Bold)
                            .color(&RGBColor(64, 64, 64))
                            .pos(Pos::new(HPos::Center, VPos::Center)),
                    )))?;
                }
            }
        }

        chart
            .configure_series_labels()
            .border_style(&TRANSPARENT)
            .position(SeriesLabelPosition::UpperRight)
            .label_font(("sans-serif", 25).into_font())
            .draw()?;
        root.present()?;
        println!("[OK] Wrote {}", fname);
        n_written += 1;
    }

    println!(
        "[INFO] Overlay comparison complete: {} groups written to {}",
        n_written, out_dir
    );
    Ok(())
}

// Helper to sanitize group names into safe filenames
fn sanitize_file_name(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '/' | '\\' | ' ' | ':' | '*' | '?' | '"' | '<' | '>' | '|' => out.push('_'),
            _ => out.push(ch),
        }
    }
    if out.is_empty() {
        "group".to_string()
    } else {
        out
    }
}

// const REPORT_EVERY: usize = 100_000;

#[derive(Debug, Clone)]
struct Point {
    chr_num: u32,
    _bp: u32,
    ind: u32,
    y: f64, // -log10(p)
}
// Extract the length of the segment between the 2nd and 3rd ':' in a variant string.
fn sv_len_from_variant(variant: &str) -> Option<usize> {
    let mut it = variant.split(':');
    let _ = it.next()?; // chr
    let _ = it.next()?; // pos
    let seg = it.next()?; // node seq between 2nd and 3rd ':'
    Some(seg.len())
}

// Determine if a line represents an SV (segment length >= sv_thresh).
fn is_sv_from_line(line: &str, sv_thresh: usize) -> bool {
    if sv_thresh == 0 { return false; }
    let mut it = line.split('\t');
    let _gene = it.next();
    if let Some(variant) = it.next() {
        if let Some(len) = sv_len_from_variant(variant.trim()) {
            return len >= sv_thresh;
        }
    }
    false
}

// #[derive(Debug)]
// struct PointWithIndex {
//     chr_num: u32,
//     ind: u32,
//     y: f64,
// }

/// Generate a soft, pastel-like palette of size `n`
fn pastel_palette(n: usize) -> Vec<RGBColor> {
    (0..n)
        .map(|i| {
            let hue = (i as f64 / n as f64) * 360.0;
            hsl_to_rgb(hue, 0.55, 0.65)
        })
        .collect()
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
fn parse_line(line: &str, re_variant: &Regex, eps: f64) -> Option<(u32, u32, f64)> {
    // Split into (up to) 9 fields; skip malformed lines
    let mut it = line.split('\t');
    let _gene = it.next()?;
    let variant = it.next()?.trim();
    // skip: dist, MA_samples, MA_count, MAF
    it.next()?;
    it.next()?;
    it.next()?;
    it.next()?;
    let p_str = it.next()?; // p
    // skip: beta, se (if present)
    // (Don't rely on them existing)

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

// fn parse_line(
//     re_variant: &Regex,
//     line: &str,
//     p_idx: usize,
//     variant_idx: Option<usize>,
//     chr_idx: Option<usize>,
//     pos_idx: Option<usize>,
// ) -> Option<Point> {
//     let fields: Vec<&str> = line.split('\t').collect();
//     if p_idx >= fields.len() {
//         return None;
//     }
//     let p_raw = fields[p_idx].trim();
//     let p: f64 = p_raw.parse().ok()?;
//     if !(p.is_finite()) || p <= 0.0 {
//         return None;
//     }
//     let (chr, pos) = if let Some(vidx) = variant_idx {
//         if vidx >= fields.len() {
//             return None;
//         }
//         let caps = re_variant.captures(fields[vidx].trim())?;
//         let chr = caps.name("chr")?.as_str();
//         let pos: u64 = caps.name("bp")?.as_str().parse().ok()?;
//         (normalize_chr(chr), pos)
//     } else {
//         let cidx = chr_idx?;
//         let pidx = pos_idx?;
//         if cidx >= fields.len() || pidx >= fields.len() {
//             return None;
//         }
//         let chr = normalize_chr(fields[cidx].trim());
//         let pos: u64 = fields[pidx].trim().parse().ok()?;
//         (chr, pos)
//     };
//     Some(Point { chr, pos, p })
// }

// fn normalize_chr(s: &str) -> String {
//     let s = s.trim();
//     let s = s.strip_prefix("chr").unwrap_or(s);
//     match s.to_ascii_uppercase().as_str() {
//         "M" | "MT" => "MT".to_string(),
//         other => other.to_string(),
//     }
// }

// pub fn run(
//     qtl_path: &str,
//     output_path: &str,
//     width: u32,
//     height: u32,
//     threads: Option<usize>,
// ) -> Result<()> {
//     run_with_roof(qtl_path, output_path, width, height, threads, None)
// }

pub fn run_with_roof(
    qtl_path: &str,
    output_path: &str,
    width: u32,
    height: u32,
    threads: Option<usize>,
    roof: Option<u32>,
    thresh_p: Option<f64>,
    no_tag: bool,
    sv_thresh: Option<usize>,
) -> Result<()> {
    if let Some(n) = threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .ok();
    }

    // Prepare regex for variant extraction (case-insensitive):
    // Accepts optional "chr" prefix, then (1..22|X|Y|M|MT), then ":" or "_", then digits
    let re_variant = Regex::new(r"(?i)^(?:chr)?(?P<chr>[0-9]{1,2}|x|y|m|mt)[:_](?P<bp>\d+)")
        .expect("invalid regex");

    let file =
        File::open(qtl_path).with_context(|| format!("Failed to open QTL file: {}", qtl_path))?;
    let reader = BufReader::new(file);

    // Collect all lines (skip header if present by checking if it starts with '#')
    let lines: Vec<String> = reader
        .lines()
        .filter_map(|l| l.ok())
        .filter(|l| !l.starts_with('#')) // ignore commented headers
        .collect();

    let sv_active = sv_thresh.unwrap_or(0) > 0;
    let sv_cut = sv_thresh.unwrap_or(0);

    // Parallel parse & compute -log10(p) and SV flag
    let eps = f64::MIN_POSITIVE; // similar to np.finfo(float).tiny
    let mut raw: Vec<(u32, u32, f64, bool)> = lines
        .par_iter()
        .filter_map(|line| {
            let sv = if sv_active { is_sv_from_line(line, sv_cut) } else { false };
            parse_line(line, &re_variant, eps).map(|(c, b, y)| (c, b, y, sv))
        })
        .collect();

    // Optional p-value threshold filter (p ≤ thresh_p  <=>  y = -log10(p) ≥ -log10(thresh_p))
    if let Some(tp) = thresh_p {
        if tp > 0.0 && tp.is_finite() {
            let y_cut = -tp.log10();
            raw.retain(|&(_, _, y, _)| y >= y_cut);
        }
    }

    if raw.is_empty() {
        anyhow::bail!(
            "No valid rows parsed after applying threshold filter. Check the input format or relax --thresh."
        );
    }

    // Parallel sort by (chr_num, bp)
    raw.par_sort_unstable_by(|a, b| match a.0.cmp(&b.0) {
        std::cmp::Ordering::Equal => a.1.cmp(&b.1),
        other => other,
    });
    let total_points = raw.len();
    // Dynamic report frequency: 1/100 of total (at least 1)
    let report_every: usize = (total_points / 100).max(1);
    // Assign running index "ind" (genomic position across chromosomes)
    let mut max_y = 0.0_f64;
    let mut out: Vec<Point> = Vec::with_capacity(total_points);
    let mut sv_out: Vec<Point> = Vec::new();
    for (i, (chr_num, _bp, y, is_sv)) in raw.into_iter().enumerate() {
        if y.is_finite() && y > max_y {
            max_y = y;
        }
        let pt = Point {
            chr_num,
            _bp,
            ind: i as u32,
            y,
        };
        if is_sv {
            sv_out.push(pt.clone());
        } else {
            out.push(pt);
        }
        if i % report_every == 0 && i > 0 {
            let pct = (i as f64 * 100.0 / total_points as f64).round() as u32;
            println!(
                "[INFO] [Panel 1/1 ({}%)] processed {:>3}/{:>3} points",
                pct, i, total_points
            );
        }
    }

    println!("[INFO] Total points parsed: {}", out.len() + sv_out.len());
    let show_chr_borders: bool = true; // toggle chromosome borders & labels

    // ==== Plot with Plotters ====
    let root = BitMapBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let x_max = out.len().saturating_sub(1) as u32;
    let y_max: u32 = roof.unwrap_or(60u32); // NEW: configurable roof

    // Determine dynamic lower bound: if thresh is active, start Y from that level
    let mut y_min: f64 = 0.0;
    if let Some(tp) = thresh_p {
        if tp > 0.0 && tp.is_finite() {
            let y_cut = -tp.log10();
            if y_cut.is_finite() && y_cut > 0.0 {
                y_min = y_cut;
            }
        }
    }
    // Ensure y_min < y_max
    if y_min >= y_max as f64 {
        y_min = (y_max as f64 - 1.0).max(0.0);
    }

    // Pastel-like colors for chromosomes 1..=25 (including X=23, Y=24, M/MT=25)
    let palette = pastel_palette(25);
    println!(
        "[INFO] Using {} pastel colors (per chromosome)",
        palette.len()
    );

    // Compute axis label according to threshold
    let thresh_active = matches!(thresh_p, Some(tp) if tp.is_finite() && tp > 0.0);
    let x_label = if thresh_active {
        "theshed Variants on Chromosomes"
    } else {
        "Variants on Chromosomes"
    };

    let mut chart = ChartBuilder::on(&root)
        // Removed caption to avoid plot title
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(0u32..x_max, y_min..y_max as f64)?;

    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc("-log10(p)")
        .disable_x_mesh()
        .disable_y_mesh()
        .draw()?;

    // GW line first (lowest layer), thin & dashed
    let gw_line = 7.30103_f64;
    let gw_draw = gw_line.max(y_min);
    if gw_draw <= y_max as f64 {
        let seg = (x_max / 80).max(2);
        let mut x = 0u32;
        while x < x_max {
            let x2 = (x + seg).min(x_max);
            chart.draw_series(std::iter::once(PathElement::new(
                vec![(x, gw_draw), (x2, gw_draw)],
                RED.stroke_width(1),
            )))?;
            x = x2 + seg;
        }
    }

    // Visual "x-axis" at threshold (now also the bottom bound)
    let label_y = (y_min + 0.8).min(y_max as f64);
    if y_min > 0.0 {
        chart.draw_series(std::iter::once(PathElement::new(
            vec![(0u32, y_min), (x_max, y_min)],
            BLACK.stroke_width(3),
        )))?;
    }

    // Draw in chunks; each point colored by its chromosome (pastel palette)
    let mut drawn: usize = 0;
    for chunk in out.chunks(report_every) {
        chart.draw_series(chunk.iter().map(|p| {
            // chr_num: 1..=25 -> index 0..=24
            let idx = (p.chr_num.saturating_sub(1) as usize) % palette.len();
            Circle::new((p.ind, p.y.min(y_max as f64)), 1, palette[idx].filled())
        }))?;
        drawn += chunk.len();
        let pct = (drawn as f64 * 100.0 / (out.len() + sv_out.len()) as f64).round() as u32;
        println!(
            "[INFO] [Panel 1/1 ({}%)] plotted {:>3}/{:>3} points",
            pct,
            drawn,
            out.len() + sv_out.len()
        );
    }

    // Draw SV points on top if present
    let sv_color = RGBColor(255, 0, 255);
    if !sv_out.is_empty() {
        chart
            .draw_series(sv_out.iter().map(|p| Circle::new((p.ind, p.y.min(y_max as f64)), 1, sv_color.filled())))?
            .label(
                std::path::Path::new(qtl_path)
                    .file_stem()
                    .unwrap_or_default()
                    .to_string_lossy()
                    .to_string() + "(SV)"
            )
            .legend(move |(x, y)| Circle::new((x, y), 5, sv_color.filled()));
    }

    // GW line already drawn above as dashed, thin, and lowest layer.

    if show_chr_borders {
        // Compute chromosome boundaries (start and end indices)
        let mut chr_boundaries: BTreeMap<u32, (u32, u32)> = BTreeMap::new();
        for p in &out {
            chr_boundaries
                .entry(p.chr_num)
                .and_modify(|e| e.1 = p.ind)
                .or_insert((p.ind, p.ind));
        }

        let mut starts: Vec<u32> = chr_boundaries
            .values()
            .map(|(start, _end)| *start)
            .collect();
        starts.sort_unstable();

        // Light grey chromosome boundary guides (toggle by `show_chr_borders`)
        for start in starts.iter().skip(1) {
            chart.draw_series(std::iter::once(PathElement::new(
                vec![(*start, 0f64), (*start, y_max as f64)],
                RGBColor(210, 210, 210).stroke_width(1),
            )))?;
        }

        // Only draw chromosome labels if no_tag is false
        if !no_tag {
            for (&chr_num, &(start, end)) in &chr_boundaries {
                let mid = (start + end) / 2;
                let chr_label = match chr_num {
                    23 => "X".to_string(),
                    24 => "Y".to_string(),
                    25 => "MT".to_string(),
                    n => n.to_string(),
                };
                chart.draw_series(std::iter::once(Text::new(
                    chr_label,
                    (mid, label_y),
                    ("sans-serif", 50)
                        .into_font()
                        .style(FontStyle::Bold)
                        .color(&RGBColor(64, 64, 64))
                        .pos(Pos::new(HPos::Center, VPos::Center)),
                )))?;
            }
        }
    }

    chart
        .configure_series_labels()
        .border_style(&TRANSPARENT) // 移除图例边框
        .draw()?;
    root.present()?;

    println!("[INFO] Plot complete. Total {} points drawn.", out.len());
    println!("[OK] Wrote {}", output_path);
    Ok(())
}

/// Draw multiple Manhattan plots stacked vertically in one image.
/// `qtl_paths_csv`: comma-separated list of QTL TSV paths.
/// `roof`: if Some(v), fixes the y-axis upper bound to `v` for each panel; otherwise defaults to 60.
pub fn run_multi(
    qtl_paths_csv: &str,
    output_path: &str,
    width: u32,
    height: u32,
    threads: Option<usize>,
    roof: Option<u32>,
    thresh_p: Option<f64>,
    no_tag: bool,
    sv_thresh: Option<usize>,
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

    // Compile regex once for all panels
    let re_variant = Regex::new(r"(?i)^(?:chr)?(?P<chr>[0-9]{1,2}|x|y|m|mt)[:_](?P<bp>\d+)")
        .expect("invalid regex");

    let root = BitMapBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;
    let n_panels = paths.len();
    let areas = root.split_evenly((n_panels, 1));

    for (panel_idx, (area, qpath)) in areas.into_iter().zip(paths.iter()).enumerate() {
        println!(
            "[INFO] [Panel {}/{}] Loading {}",
            panel_idx + 1,
            n_panels,
            qpath
        );

        // Read and preprocess this QTL
        let file =
            File::open(qpath).with_context(|| format!("Failed to open QTL file: {}", qpath))?;
        let reader = BufReader::new(file);
        let lines: Vec<String> = reader
            .lines()
            .filter_map(|l| l.ok())
            .filter(|l| !l.starts_with('#'))
            .collect();

        let sv_active = sv_thresh.unwrap_or(0) > 0;
        let sv_cut = sv_thresh.unwrap_or(0);
        let eps = f64::MIN_POSITIVE;
        let mut raw: Vec<(u32, u32, f64, bool)> = lines
            .par_iter()
            .filter_map(|line| {
                let sv = if sv_active { is_sv_from_line(line, sv_cut) } else { false };
                parse_line(line, &re_variant, eps).map(|(c, b, y)| (c, b, y, sv))
            })
            .collect();

        // Optional p-value threshold filter
        if let Some(tp) = thresh_p {
            if tp > 0.0 && tp.is_finite() {
                let y_cut = -tp.log10();
                raw.retain(|&(_, _, y, _)| y >= y_cut);
            }
        }

        if raw.is_empty() {
            anyhow::bail!(
                "No valid rows parsed for panel {} ({}) after threshold filter. Check the input format or relax --thresh.",
                panel_idx + 1,
                qpath
            );
        }

        let total_points = raw.len();
        // Dynamic report frequency: 1/100 of total (at least 1)
        let report_every: usize = (total_points / 100).max(1);

        raw.par_sort_unstable_by(|a, b| match a.0.cmp(&b.0) {
            std::cmp::Ordering::Equal => a.1.cmp(&b.1),
            other => other,
        });

        let mut out: Vec<Point> = Vec::with_capacity(total_points);
        let mut sv_out: Vec<Point> = Vec::new();
        let mut max_y = 0.0_f64;
        for (i, (chr_num, _bp, y, is_sv)) in raw.into_iter().enumerate() {
            if y.is_finite() && y > max_y {
                max_y = y;
            }
            let pt = Point {
                chr_num,
                _bp,
                ind: i as u32,
                y,
            };
            if is_sv {
                sv_out.push(pt.clone());
            } else {
                out.push(pt);
            }
            if i % report_every == 0 && i > 0 {
                let pct = (i as f64 * 100.0 / total_points as f64).round() as u32;
                println!(
                    "[INFO] [Panel {}/{} ({}%)] processed {:>3}/{:>3} points",
                    panel_idx + 1,
                    n_panels,
                    pct,
                    i,
                    total_points
                );
            }
        }

        let x_max = out.len().saturating_sub(1) as u32;
        let y_max: u32 = roof.unwrap_or(60u32);

        // Determine dynamic lower bound: if thresh is active, start Y from that level
        let mut y_min: f64 = 0.0;
        if let Some(tp) = thresh_p {
            if tp > 0.0 && tp.is_finite() {
                let y_cut = -tp.log10();
                if y_cut.is_finite() && y_cut > 0.0 {
                    y_min = y_cut;
                }
            }
        }
        // Ensure y_min < y_max
        if y_min >= y_max as f64 {
            y_min = (y_max as f64 - 1.0).max(0.0);
        }

        // Compute axis label according to threshold
        let thresh_active = matches!(thresh_p, Some(tp) if tp.is_finite() && tp > 0.0);
        let x_label = if thresh_active {
            "theshed Variants on Chromosomes"
        } else {
            "Variants on Chromosomes"
        };

        // Draw one Manhattan in this sub-area
        let sub = &area;
        let mut chart = ChartBuilder::on(sub)
            // no caption to avoid crowding; keep defaults
            .margin(10)
            .x_label_area_size(40)
            .y_label_area_size(60)
            .build_cartesian_2d(0u32..x_max, y_min..y_max as f64)?;

        chart
            .configure_mesh()
            .x_desc(x_label)
            .y_desc("-log10(p)")
            .disable_x_mesh()
            .disable_y_mesh()
            .draw()?;

        // GW line first (lowest layer), thin & dashed
        let gw_line = 7.30103_f64;
        let gw_draw = gw_line.max(y_min);
        if gw_draw <= y_max as f64 {
            let seg = (x_max / 80).max(2);
            let mut x = 0u32;
            while x < x_max {
                let x2 = (x + seg).min(x_max);
                chart.draw_series(std::iter::once(PathElement::new(
                    vec![(x, gw_draw), (x2, gw_draw)],
                    RED.stroke_width(1),
                )))?;
                x = x2 + seg;
            }
        }

        let label_y = (y_min + 0.8).min(y_max as f64);
        if y_min > 0.0 {
            chart.draw_series(std::iter::once(PathElement::new(
                vec![(0u32, y_min), (x_max, y_min)],
                BLACK.stroke_width(3),
            )))?;
        }

        let palette = pastel_palette(25);

        let mut drawn: usize = 0;
        for chunk in out.chunks(report_every) {
            chart.draw_series(chunk.iter().map(|p| {
                let idx = (p.chr_num.saturating_sub(1) as usize) % palette.len();
                Circle::new((p.ind, p.y.min(y_max as f64)), 1, palette[idx].filled())
            }))?;
            drawn += chunk.len();
            let pct = (drawn as f64 * 100.0 / (out.len() + sv_out.len()) as f64).round() as u32;
            println!(
                "[INFO] [Panel {}/{} ({}%)] plotted {:>3}/{:>3} points",
                panel_idx + 1,
                n_panels,
                pct,
                drawn,
                out.len() + sv_out.len()
            );
        }

        // Draw SV points on top if present
        let sv_color = RGBColor(255, 0, 255);
        if !sv_out.is_empty() {
            chart
                .draw_series(sv_out.iter().map(|p| Circle::new((p.ind, p.y.min(y_max as f64)), 1, sv_color.filled())))?
                .label(
                    std::path::Path::new(qpath)
                        .file_stem()
                        .unwrap_or_default()
                        .to_string_lossy()
                        .to_string() + "(SV)"
                )
                .legend(move |(x, y)| Circle::new((x, y), 5, sv_color.filled()));
        }

        // GW line already drawn above as dashed, thin, and lowest layer.

        // Chromosome borders and labels
        let mut chr_boundaries: BTreeMap<u32, (u32, u32)> = BTreeMap::new();
        for p in &out {
            chr_boundaries
                .entry(p.chr_num)
                .and_modify(|e| e.1 = p.ind)
                .or_insert((p.ind, p.ind));
        }
        let mut starts: Vec<u32> = chr_boundaries
            .values()
            .map(|(start, _end)| *start)
            .collect();
        starts.sort_unstable();
        for start in starts.iter().skip(1) {
            chart.draw_series(std::iter::once(PathElement::new(
                vec![(*start, 0f64), (*start, y_max as f64)],
                RGBColor(210, 210, 210).stroke_width(1),
            )))?;
        }
        if !no_tag {
            for (&chr_num, &(start, end)) in &chr_boundaries {
                let mid = (start + end) / 2;
                let chr_label = match chr_num {
                    23 => "X".to_string(),
                    24 => "Y".to_string(),
                    25 => "MT".to_string(),
                    n => n.to_string(),
                };
                chart.draw_series(std::iter::once(Text::new(
                    chr_label,
                    (mid, label_y),
                    ("sans-serif", 50)
                        .into_font()
                        .style(FontStyle::Bold)
                        .color(&RGBColor(64, 64, 64))
                        .pos(Pos::new(HPos::Center, VPos::Center)),
                )))?;
            }
        }

        chart
            .configure_series_labels()
            .border_style(&TRANSPARENT)
            .draw()?;
    }

    // Present once for the combined image
    root.present()?;

    // Report source order (top -> bottom)
    println!("[INFO] Source order (top -> bottom):");
    for (idx, pth) in paths.iter().enumerate() {
        println!("  [{}] {}", idx + 1, pth);
    }

    println!("[OK] Wrote {}", output_path);
    Ok(())
}

// fn find_first(cols: &[String], names: &[&str]) -> Option<usize> {
//     for (i, c) in cols.iter().enumerate() {
//         let lc = c.to_lowercase();
//         if names.iter().any(|n| lc == n.to_lowercase()) {
//             return Some(i);
//         }
//     }
//     None
// }

/// Plot the intersection of multiple QTL files: only variants present (with y >= -log10(thresh)) in all files are shown.
pub fn run_intersection(
    qtl_paths_csv: &str,
    output_path: &str,
    width: u32,
    height: u32,
    threads: Option<usize>,
    roof: Option<u32>,
    thresh_p: Option<f64>,
    no_tag: bool,
    sv_thresh: Option<usize>,
) -> Result<()> {
    if let Some(n) = threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(n)
            .build_global()
            .ok();
    }

    let thresh = match thresh_p {
        Some(tp) if tp.is_finite() && tp > 0.0 => tp,
        _ => anyhow::bail!("--thresh must be provided and > 0 for intersection mode."),
    };
    let y_cut = -thresh.log10();
    let paths: Vec<String> = qtl_paths_csv
        .split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect();
    if paths.len() < 2 {
        anyhow::bail!("Intersection requires at least two input files.");
    }

    let re_variant = Regex::new(r"(?i)^(?:chr)?(?P<chr>[0-9]{1,2}|x|y|m|mt)[:_](?P<bp>\d+)")
        .expect("invalid regex");
    let eps = f64::MIN_POSITIVE;

    // Helper: parse line, filter by y >= y_cut
    fn parse_intersection_line(
        line: &str,
        re_variant: &Regex,
        eps: f64,
        y_cut: f64,
    ) -> Option<(u32, u32, f64)> {
        parse_line(line, re_variant, eps).and_then(|(chr_num, bp, y)| {
            if y >= y_cut {
                Some((chr_num, bp, y))
            } else {
                None
            }
        })
    }

    // For each file, build HashMap<(chr_num, bp), max y>
    let sv_active = sv_thresh.unwrap_or(0) > 0;
    let sv_cut = sv_thresh.unwrap_or(0);
    let mut maps: Vec<HashMap<(u32, u32), f64>> = Vec::with_capacity(paths.len());
    let mut sv_sets: Vec<BTreeSet<(u32, u32)>> = Vec::with_capacity(paths.len());
    for (i, path) in paths.iter().enumerate() {
        println!("[INFO] [File {}/{}] Loading {}", i + 1, paths.len(), path);
        let file =
            File::open(path).with_context(|| format!("Failed to open QTL file: {}", path))?;
        let reader = BufReader::new(file);
        let mut map: HashMap<(u32, u32), f64> = HashMap::new();
        let mut sv_set: BTreeSet<(u32, u32)> = BTreeSet::new();
        for line in reader.lines().filter_map(|l| l.ok()) {
            if line.starts_with('#') { continue; }
            let is_sv = sv_active && is_sv_from_line(&line, sv_cut);
            if let Some((chr_num, bp, y)) = parse_intersection_line(&line, &re_variant, eps, y_cut) {
                let key = (chr_num, bp);
                map.entry(key).and_modify(|val| { if y > *val { *val = y; } }).or_insert(y);
                if is_sv { sv_set.insert(key); }
            }
        }
        if map.is_empty() {
            anyhow::bail!(
                "No valid points found in file {} after threshold filter.",
                path
            );
        }
        println!(
            "[INFO] File {}: {} points >= -log10({})",
            path,
            map.len(),
            thresh
        );
        maps.push(map);
        sv_sets.push(sv_set);
    }

    // Compute intersection of keys across all maps
    let mut intersect_keys: BTreeSet<(u32, u32)> = maps[0].keys().cloned().collect();
    for map in &maps[1..] {
        let keys: BTreeSet<(u32, u32)> = map.keys().cloned().collect();
        intersect_keys = intersect_keys.intersection(&keys).cloned().collect();
    }
    if intersect_keys.is_empty() {
        anyhow::bail!(
            "No variants found in common (with y >= -log10({})) across all files.",
            thresh
        );
    }
    println!(
        "[INFO] Intersection: {} variants present in all {} files.",
        intersect_keys.len(),
        maps.len()
    );
    // Compute union of SV keys within the intersection
    let mut sv_keys_any: BTreeSet<(u32, u32)> = BTreeSet::new();
    if sv_active {
        for s in &sv_sets {
            for k in s {
                if intersect_keys.contains(k) { sv_keys_any.insert(*k); }
            }
        }
    }

    // For each key, take the maximum y across all maps
    let mut keys_sorted: Vec<(u32, u32)> = intersect_keys.iter().cloned().collect();
    keys_sorted.sort_unstable();
    let mut out: Vec<Point> = Vec::with_capacity(keys_sorted.len());
    let mut max_y = 0.0_f64;
    for (i, &(chr_num, bp)) in keys_sorted.iter().enumerate() {
        let y = maps
            .iter()
            .map(|m| m[&(chr_num, bp)])
            .fold(f64::MIN, f64::max);
        if y > max_y {
            max_y = y;
        }
        out.push(Point {
            chr_num,
            _bp: bp,
            ind: i as u32,
            y,
        });
    }
    println!(
        "[INFO] Ready to plot {} intersection points. Max y = {:.2}",
        out.len(),
        max_y
    );

    let root = BitMapBackend::new(output_path, (width, height)).into_drawing_area();
    root.fill(&WHITE)?;

    let x_max = out.len().saturating_sub(1) as u32;
    let y_max: u32 = roof.unwrap_or(60u32);
    let mut y_min: f64 = y_cut;
    if y_min >= y_max as f64 {
        y_min = (y_max as f64 - 1.0).max(0.0);
    }
    let palette = pastel_palette(25);
    let x_label = "theshed Variants on Chromosomes";

    let mut chart = ChartBuilder::on(&root)
        .margin(10)
        .x_label_area_size(40)
        .y_label_area_size(60)
        .build_cartesian_2d(0u32..x_max, y_min..y_max as f64)?;
    chart
        .configure_mesh()
        .x_desc(x_label)
        .y_desc("-log10(p)")
        .disable_x_mesh()
        .disable_y_mesh()
        .draw()?;

    // GW line first (lowest layer), thin & dashed
    let gw_line = 7.30103_f64;
    let gw_draw = gw_line.max(y_min);
    if gw_draw <= y_max as f64 {
        let seg = (x_max / 80).max(2);
        let mut x = 0u32;
        while x < x_max {
            let x2 = (x + seg).min(x_max);
            chart.draw_series(std::iter::once(PathElement::new(
                vec![(x, gw_draw), (x2, gw_draw)],
                RED.stroke_width(1),
            )))?;
            x = x2 + seg;
        }
    }

    let label_y = (y_min + 0.8).min(y_max as f64);
    if y_min > 0.0 {
        chart.draw_series(std::iter::once(PathElement::new(
            vec![(0u32, y_min), (x_max, y_min)],
            BLACK.stroke_width(3),
        )))?;
    }

    let report_every = (out.len() / 100).max(1);
    let mut drawn = 0;
    for chunk in out.chunks(report_every) {
        chart.draw_series(chunk.iter().map(|p| {
            let idx = (p.chr_num.saturating_sub(1) as usize) % palette.len();
            Circle::new((p.ind, p.y.min(y_max as f64)), 1, palette[idx].filled())
        }))?;
        drawn += chunk.len();
        let pct = (drawn as f64 * 100.0 / out.len() as f64).round() as u32;
        println!(
            "[INFO] [Intersection ({}%)] plotted {:>3}/{:>3} points",
            pct,
            drawn,
            out.len()
        );
    }
    // Overlay SV points and add a legend if any
    if !sv_keys_any.is_empty() {
        let sv_color = RGBColor(255, 0, 255);
        chart
            .draw_series(
                out.iter()
                    .filter(|p| sv_keys_any.contains(&(p.chr_num, p._bp)))
                    .map(|p| Circle::new((p.ind, p.y.min(y_max as f64)), 3, sv_color.filled()))
            )?
            .label("SV")
            .legend(move |(x, y)| Circle::new((x, y), 5, sv_color.filled()));
    }

    // GW line already drawn above as dashed, thin, and lowest layer.

    // Chromosome borders and labels
    let mut chr_boundaries: BTreeMap<u32, (u32, u32)> = BTreeMap::new();
    for p in &out {
        chr_boundaries
            .entry(p.chr_num)
            .and_modify(|e| e.1 = p.ind)
            .or_insert((p.ind, p.ind));
    }
    let mut starts: Vec<u32> = chr_boundaries
        .values()
        .map(|(start, _end)| *start)
        .collect();
    starts.sort_unstable();
    for start in starts.iter().skip(1) {
        chart.draw_series(std::iter::once(PathElement::new(
            vec![(*start, 0f64), (*start, y_max as f64)],
            RGBColor(210, 210, 210).stroke_width(1),
        )))?;
    }
    if !no_tag {
        for (&chr_num, &(start, end)) in &chr_boundaries {
            let mid = (start + end) / 2;
            let chr_label = match chr_num {
                23 => "X".to_string(),
                24 => "Y".to_string(),
                25 => "MT".to_string(),
                n => n.to_string(),
            };
            chart.draw_series(std::iter::once(Text::new(
                chr_label,
                (mid, label_y),
                ("sans-serif", 50)
                    .into_font()
                    .style(FontStyle::Bold)
                    .color(&RGBColor(64, 64, 64))
                    .pos(Pos::new(HPos::Center, VPos::Center)),
            )))?;
        }
    }
    chart
        .configure_series_labels()
        .border_style(&TRANSPARENT)
        .draw()?;
    root.present()?;

    println!(
        "[INFO] Plot complete. Total {} intersection points drawn.",
        out.len()
    );
    println!("[OK] Wrote {}", output_path);
    Ok(())
}



// /// Compare multiple QTL files by group (e.g., gene), plotting one stacked image per group.
// pub fn run_compare_grouped(
//     qtl_paths_csv: &str,
//     out_dir: &str,
//     width: u32,
//     height: u32,
//     threads: Option<usize>,
//     roof: Option<u32>,
//     thresh_p: Option<f64>,
// ) -> Result<()> {
//     if let Some(n) = threads {
//         rayon::ThreadPoolBuilder::new()
//             .num_threads(n)
//             .build_global()
//             .ok();
//     }
//     let paths: Vec<String> = qtl_paths_csv
//         .split(',')
//         .map(|s| s.trim().to_string())
//         .filter(|s| !s.is_empty())
//         .collect();
//     if paths.len() < 2 {
//         anyhow::bail!("Grouped compare requires at least two input files.");
//     }
//     let re_variant = Regex::new(r"(?i)^(?:chr)?(?P<chr>[0-9]{1,2}|x|y|m|mt)[:_](?P<bp>\d+)")
//         .expect("invalid regex");
//     let eps = f64::MIN_POSITIVE;
//     let y_cut = thresh_p
//         .filter(|tp| tp.is_finite() && *tp > 0.0)
//         .map(|tp| -tp.log10());

//     // For each file, group by first column (phenotype/gene)
//     let mut per_file: Vec<HashMap<String, Vec<(u32, u32, f64)>>> = Vec::with_capacity(paths.len());
//     for (i, path) in paths.iter().enumerate() {
//         println!("[INFO] [File {}/{}] Loading {}", i + 1, paths.len(), path);
//         let file = File::open(path)
//             .with_context(|| format!("Failed to open QTL file: {}", path))?;
//         let reader = BufReader::new(file);
//         let mut group_map: HashMap<String, Vec<(u32, u32, f64)>> = HashMap::new();
//         for line in reader.lines().filter_map(|l| l.ok()) {
//             if line.starts_with('#') {
//                 continue;
//             }
//             let mut it = line.split('\t');
//             let group = match it.next() {
//                 Some(s) => s.trim().to_string(),
//                 _none => continue,
//             };
//             // Rebuild the remainder of the line (excluding the first column) for reuse in parse_line
//             let rest_joined = it.collect::<Vec<&str>>().join("\t");
//             let pseudo = format!("{}\t{}", &group, rest_joined);
//             if let Some((chr_num, bp, y)) = parse_line(&pseudo, &re_variant, eps) {
//                 if let Some(cut) = y_cut {
//                     if y < cut {
//                         continue;
//                     }
//                 }
//                 group_map
//                     .entry(group.to_string())
//                     .or_insert_with(Vec::new)
//                     .push((chr_num, bp, y));
//             }
//         }
//         println!(
//             "[INFO] File {}: {} groups loaded.",
//             path,
//             group_map.len()
//         );
//         per_file.push(group_map);
//     }

//     // Collect all group names across all files
//     let mut group_names: BTreeSet<String> = BTreeSet::new();
//     for gmap in &per_file {
//         for k in gmap.keys() {
//             group_names.insert(k.clone());
//         }
//     }
//     if group_names.is_empty() {
//         anyhow::bail!("No groups found in any input file.");
//     }
//     let n_panels = paths.len();
//     let mut n_written = 0usize;
//     for group in &group_names {
//         let fname = format!("{}/{}.png", out_dir, sanitize_file_name(group));
//         let root = BitMapBackend::new(&fname, (width, height)).into_drawing_area();
//         root.fill(&WHITE)?;
//         let areas = root.split_evenly((n_panels, 1));
//         for (panel_idx, (area, gmap)) in areas.into_iter().zip(per_file.iter()).enumerate() {
//             let vals = gmap.get(group).cloned().unwrap_or_default();
//             if vals.is_empty() {
//                 println!(
//                     "[INFO] [Group {}] Panel {}/{}: group missing in file {}.",
//                     group,
//                     panel_idx + 1,
//                     n_panels,
//                     paths[panel_idx]
//                 );
//             }
//             let mut points = vals;
//             points.sort_unstable_by(|a, b| match a.0.cmp(&b.0) {
//                 std::cmp::Ordering::Equal => a.1.cmp(&b.1),
//                 other => other,
//             });
//             let mut out: Vec<Point> = Vec::with_capacity(points.len());
//             let mut max_y = 0.0_f64;
//             for (i, (chr_num, bp, y)) in points.iter().enumerate() {
//                 if *y > max_y {
//                     max_y = *y;
//                 }
//                 out.push(Point {
//                     chr_num: *chr_num,
//                     _bp: *bp,
//                     ind: i as u32,
//                     y: *y,
//                 });
//             }
//             let x_max = out.len().saturating_sub(1) as u32;
//             let y_max: u32 = roof.unwrap_or(60u32);
//             let mut y_min: f64 = y_cut.unwrap_or(0.0);
//             if y_min >= y_max as f64 {
//                 y_min = (y_max as f64 - 1.0).max(0.0);
//             }
//             let thresh_active = y_cut.is_some();
//             let x_label = if thresh_active {
//                 "theshed Variants on Chromosomes"
//             } else {
//                 "Variants on Chromosomes"
//             };
//             let mut chart = ChartBuilder::on(&area)
//                 .margin(10)
//                 .x_label_area_size(40)
//                 .y_label_area_size(60)
//                 .build_cartesian_2d(0u32..x_max, y_min..y_max as f64)?;
//             chart
//                 .configure_mesh()
//                 .x_desc(x_label)
//                 .y_desc("-log10(p)")
//                 .disable_x_mesh()
//                 .disable_y_mesh()
//                 .draw()?;
//             let label_y = (y_min + 0.8).min(y_max as f64);
//             if y_min > 0.0 {
//                 chart.draw_series(std::iter::once(PathElement::new(
//                     vec![(0u32, y_min), (x_max, y_min)],
//                     BLACK.stroke_width(3),
//                 )))?;
//             }
//             let palette = pastel_palette(25);
//             let report_every = (out.len() / 100).max(1);
//             let mut _drawn = 0;
//             for chunk in out.chunks(report_every.max(1)) {
//                 chart.draw_series(chunk.iter().map(|p| {
//                     let idx = (p.chr_num.saturating_sub(1) as usize) % palette.len();
//                     Circle::new((p.ind, p.y.min(y_max as f64)), 1, palette[idx].filled())
//                 }))?;
//                 _drawn += chunk.len();
//             }
//             // GW line (5e-8) if visible
//             let gw_line = 7.30103_f64;
//             let gw_draw = gw_line.max(y_min);
//             if gw_draw <= y_max as f64 {
//                 chart.draw_series(LineSeries::new(
//                     vec![(0u32, gw_draw), (x_max, gw_draw)],
//                     RED.stroke_width(2),
//                 ))?;
//             }
//             // Chromosome borders and labels
//             let mut chr_boundaries: BTreeMap<u32, (u32, u32)> = BTreeMap::new();
//             for p in &out {
//                 chr_boundaries
//                     .entry(p.chr_num)
//                     .and_modify(|e| e.1 = p.ind)
//                     .or_insert((p.ind, p.ind));
//             }
//             let mut starts: Vec<u32> = chr_boundaries
//                 .values()
//                 .map(|(start, _end)| *start)
//                 .collect();
//             starts.sort_unstable();
//             for start in starts.iter().skip(1) {
//                 chart.draw_series(std::iter::once(PathElement::new(
//                     vec![(*start, 0f64), (*start, y_max as f64)],
//                     RGBColor(210, 210, 210).stroke_width(1),
//                 )))?;
//             }
//             for (&chr_num, &(start, end)) in &chr_boundaries {
//                 let mid = (start + end) / 2;
//                 let chr_label = match chr_num {
//                     23 => "X".to_string(),
//                     24 => "Y".to_string(),
//                     25 => "MT".to_string(),
//                     n => n.to_string(),
//                 };
//                 chart.draw_series(std::iter::once(Text::new(
//                     chr_label,
//                     (mid, label_y),
//                     ("sans-serif", 50)
//                         .into_font()
//                         .style(FontStyle::Bold)
//                         .color(&RGBColor(64, 64, 64))
//                         .pos(Pos::new(HPos::Center, VPos::Center)),
//                 )))?;
//             }
//             chart
//                 .configure_series_labels()
//                 .border_style(&TRANSPARENT)
//                 .draw()?;
//         }
//         root.present()?;
//         println!("[OK] Wrote {}", fname);
//         n_written += 1;
//     }
//     println!(
//         "[INFO] Grouped comparison complete: {} groups written to {}",
//         n_written, out_dir
//     );
//     Ok(())
// }
