use std::path::{Path, PathBuf};

use anyhow::Result;
use clap::{Arg, ArgGroup, Command, value_parser};

mod man; // provides: pub fn run(qtl:&str, output:&str, width:u32, height:u32, threads: Option<usize>) -> anyhow::Result<()>
mod pro;
mod qq;
mod vio;

fn main() {
    if let Err(e) = real_main() {
        eprintln!("[ERROR] {e:?}");
        std::process::exit(1);
    }
}

fn real_main() -> Result<()> {
    let app = Command::new("eQTLplotter")
        .version(env!("CARGO_PKG_VERSION"))
        .about("Plot eQTL figures (Manhattan, Violin, QQ Plot, etc.)")
        .subcommand(
            Command::new("man")
                .about("Manhattan plot from QTL results")
                .alias("manhattan")
                .arg(
                    Arg::new("qtl")
                        .help("QTL TSV file (needs p-value and variant or chr+pos)")
                        .short('q')
                        .long("qtl")
                        .required(false)
                        .value_name("FILE"),
                )
                .arg(
                    Arg::new("output")
                        .help("Output PNG path. If omitted, use same dir/name as --qtl with .png")
                        .short('o')
                        .long("output")
                        .required(false)
                        .value_name("PNG"),
                )
                .arg(
                    Arg::new("width")
                        .help("Output width in pixels (default ~4K)")
                        .long("width")
                        .required(false)
                        .value_parser(value_parser!(u32))
                        .default_value("3840"),
                )
                .arg(
                    Arg::new("height")
                        .help("Output height in pixels (default ~4K 16:9)")
                        .long("height")
                        .required(false)
                        .value_parser(value_parser!(u32))
                        .default_value("1920"),
                )
                .arg(
                    Arg::new("threads")
                        .help("Rayon worker threads (optional)")
                        .long("threads")
                        .short('T')
                        .required(false)
                        .value_parser(value_parser!(usize)),
                )
                .arg(
                    Arg::new("roof")
                        .help("Upper y-axis cap for -log10(p) (applies to single and each panel in --multi)")
                        .long("roof")
                        .required(false)
                        .value_parser(value_parser!(u32))
                        .value_name("U32"),
                )
                .arg(
                    Arg::new("multi")
                        .help("Comma-separated QTL TSV paths. In 'sep' mode, behaves like before (stacked panels). In 'inter'/'compare' modes, must contain at least two files.")
                        .long("multi")
                        .required(false)
                        .value_name("CSV"),
                )
                .group(
                    ArgGroup::new("input")
                        .args(["qtl", "multi"])
                        .required(true),
                )
                .arg(
                    Arg::new("thresh")
                        .help("P-value threshold to plot (keep points with p ≤ THRESH). Use 0 to disable filtering.")
                        .long("thresh")
                        .required(false)
                        .value_parser(value_parser!(f64))
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::new("sv")
                        .help("SV node-length threshold (usize). If > 0, highlight variants whose node sequence (between 2nd and 3rd ':') length ≥ SV; plot them on top and add '(SV)' legend.")
                        .long("sv")
                        .required(false)
                        .value_parser(value_parser!(usize))
                        .value_name("USIZE"),
                )
                .arg(
                    Arg::new("over")
                        .help("Only emit figures that contain at least one point with p ≤ OVER (overlay filter). Different from --thresh, which filters points; --over filters which figures are produced.")
                        .long("over")
                        .required(false)
                        .value_parser(value_parser!(f64))
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::new("mode")
                        .help("Plot mode: 'sep' (default, current behavior), 'inter' (intersection of significant genotypes across files), 'compare' (overlay: plot points from different files on the same axes; color encodes file; outputs per-group PNGs to a directory)")
                        .long("mode")
                        .required(false)
                        .value_parser(["sep", "inter", "compare"]) 
                        .default_value("sep")
                        .value_name("MODE"),
                )
                .arg(
                    Arg::new("no-tag")
                        .help("Disable chromosome number labeling in Manhattan plots")
                        .long("no-tag")
                        .required(false)
                        .action(clap::ArgAction::SetTrue),
                )
        )
        .subcommand(
            Command::new("vio")
                .about("Violin plot over phenotype-variant counts across multiple QTL files")
                .arg(
                    Arg::new("multi")
                        .help("Comma-separated QTL TSV paths (2+). Each file contributes one violin; x-axis uses file stem names.")
                        .long("multi")
                        .required(true)
                        .value_name("CSV"),
                )
                .arg(
                    Arg::new("output")
                        .help("Output PNG path. If omitted, default is <dir_of_first_qtl>/violin_{stem1}_{stem2}_... .png")
                        .short('o')
                        .long("output")
                        .required(false)
                        .value_name("PNG"),
                )
                .arg(
                    Arg::new("width")
                        .help("Output width in pixels (default ~4K)")
                        .long("width")
                        .required(false)
                        .value_parser(value_parser!(u32))
                        .default_value("3840"),
                )
                .arg(
                    Arg::new("height")
                        .help("Output height in pixels (default 2160)")
                        .long("height")
                        .required(false)
                        .value_parser(value_parser!(u32))
                        .default_value("2160"),
                )
                .arg(
                    Arg::new("threads")
                        .help("Rayon worker threads (optional)")
                        .long("threads")
                        .short('T')
                        .required(false)
                        .value_parser(value_parser!(usize)),
                )
                .arg(
                    Arg::new("roof")
                        .help("Upper y-axis cap for '#variants per phenotype'")
                        .long("roof")
                        .required(false)
                        .value_parser(value_parser!(u32))
                        .value_name("U32"),
                )
                .arg(
                    Arg::new("box")
                        .help("Overlay a box plot (median/IQR/whiskers) on each violin")
                        .long("box")
                        .required(false)
                        .action(clap::ArgAction::SetTrue),
                )
        )
        .subcommand(
            Command::new("pro")
                .about("Proportional bar plot comparing two QTL files by phenotype-level significant counts")
                .arg(
                    Arg::new("multi")
                        .help("CSV with exactly two QTL TSV paths (left,right)")
                        .long("multi")
                        .required(true)
                        .value_name("FILE_A,FILE_B"),
                )
                .arg(
                    Arg::new("output")
                        .help("Output PNG path. If omitted, default is pro_{stemA}_vs_{stemB}.png next to the first file")
                        .short('o')
                        .long("output")
                        .required(false)
                        .value_name("PNG"),
                )
                .arg(
                    Arg::new("width")
                        .help("Output width in pixels (default ~FHD)")
                        .long("width")
                        .required(false)
                        .value_parser(value_parser!(u32))
                        .default_value("2560"),
                )
                .arg(
                    Arg::new("height")
                        .help("Output height in pixels (default 1080)")
                        .long("height")
                        .required(false)
                        .value_parser(value_parser!(u32))
                        .default_value("1440"),
                )
                .arg(
                    Arg::new("threads")
                        .help("Rayon worker threads (optional)")
                        .long("threads")
                        .short('T')
                        .required(false)
                        .value_parser(value_parser!(usize)),
                )
                .arg(
                    Arg::new("thresh")
                        .help("Significance threshold in p-value (default 5e-8)")
                        .long("thresh")
                        .required(false)
                        .value_parser(value_parser!(f64))
                        .default_value("5e-8")
                        .value_name("FLOAT"),
                )
                .arg(
                    Arg::new("bin")
                        .help("Bin width in percent along x-axis (default 1.0)")
                        .long("bin")
                        .required(false)
                        .value_parser(value_parser!(f64))
                        .default_value("5.0")
                        .value_name("PERCENT"),
                )
                .arg(
                    Arg::new("diff")
                        .help("Bidirectional tail color threshold in [0,1]; e.g., 0.9 colors ≤10% and ≥90% bins differently")
                        .long("diff")
                        .required(false)
                        .value_parser(value_parser!(f64))
                        .value_name("RATIO"),
                )
                .arg(
                    Arg::new("venn")
                        .help("Show a Venn diagram (phenotype-level) at the center of the plot")
                        .long("venn")
                        .required(false)
                        .action(clap::ArgAction::SetTrue),
                )
        )
        .subcommand(
            Command::new("qq")
                .about("QQ scatter comparing p-values from exactly two QTL files")
                .arg(
                    Arg::new("multi")
                        .help("CSV with exactly two QTL TSV paths (first=Y axis, second=X axis)")
                        .long("multi")
                        .required(true)
                        .value_name("FILE_A,FILE_B"),
                )
                .arg(
                    Arg::new("output")
                        .help("Output PNG path. If omitted, default is qq_{stemA}_vs_{stemB}.png next to the first file")
                        .short('o')
                        .long("output")
                        .required(false)
                        .value_name("PNG"),
                )
                .arg(
                    Arg::new("width")
                        .help("Output width in pixels (default ~4K)")
                        .long("width")
                        .required(false)
                        .value_parser(value_parser!(u32))
                        .default_value("1800"),
                )
                .arg(
                    Arg::new("height")
                        .help("Output height in pixels (default 2160)")
                        .long("height")
                        .required(false)
                        .value_parser(value_parser!(u32))
                        .default_value("1800"),
                )
                .arg(
                    Arg::new("threads")
                        .help("Rayon worker threads (optional)")
                        .long("threads")
                        .short('T')
                        .required(false)
                        .value_parser(value_parser!(usize)),
                )
                .arg(
                    Arg::new("group")
                        .help("Group to phenotype-level (column 0) by taking min p per phenotype; otherwise join on variant id (column 1)")
                        .long("group")
                        .required(false)
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("qcut")
                        .help("Drop the top PERCENT of points by max(-log10(px), -log10(py)) to make plot denser (e.g., 1 = top 1%)")
                        .long("qcut")
                        .required(false)
                        .value_parser(value_parser!(f64))
                        .value_name("PERCENT"),
                )
                .arg(
                    Arg::new("diff")
                        .help("Use different colors for points above (y>x) vs below (y<=x) the diagonal")
                        .long("diff")
                        .required(false)
                        .action(clap::ArgAction::SetTrue),
                )
                .arg(
                    Arg::new("roof")
                        .help("Cap -log10(p) values at this ceiling before plotting")
                        .long("roof")
                        .required(false)
                        .value_parser(value_parser!(f64))
                        .value_name("FLOAT"),
                )
        );

    let matches = app.get_matches();

    match matches.subcommand() {
        Some(("man", sub)) => {
            let multi_opt: Option<String> = sub.get_one::<String>("multi").cloned();

            let mode: String = sub
                .get_one::<String>("mode")
                .cloned()
                .unwrap_or_else(|| "sep".to_string());

            // Decide input mode:
            let (use_multi, qtl_str_or_csv): (bool, String) = if let Some(multi) = multi_opt.clone()
            {
                if multi.contains(',') {
                    (true, multi)
                } else {
                    // Treat single path the same as --qtl
                    (false, multi)
                }
            } else {
                let qtl: &String = sub
                    .get_one::<String>("qtl")
                    .expect("either --qtl or --multi is required");
                (false, qtl.to_string())
            };

            let width = *sub
                .get_one::<u32>("width")
                .expect("default provided by clap");
            let height = *sub
                .get_one::<u32>("height")
                .expect("default provided by clap");
            let threads = sub.get_one::<usize>("threads").copied();
            let roof_opt = sub.get_one::<u32>("roof").copied();
            let thresh_opt = sub.get_one::<f64>("thresh").copied();
            let sv_opt = sub.get_one::<usize>("sv").copied();
            let over_opt = sub.get_one::<f64>("over").copied();

            let no_tag = sub.get_flag("no-tag");
            if no_tag {
                println!("[INFO] NoTag : enabled");
            }

            if use_multi {
                println!("[INFO] MULTI : {}", qtl_str_or_csv);
            } else {
                println!("[INFO] QTL   : {}", qtl_str_or_csv);
            }
            println!("[INFO] Mode  : {}", mode);
            println!("[INFO] Size  : {}x{}", width, height);
            if let Some(n) = threads {
                println!("[INFO] rayon threads = {}", n);
            }
            if let Some(roof) = roof_opt {
                println!("[INFO] Roof  : {}", roof);
            }
            if let Some(tp) = thresh_opt {
                println!("[INFO] Thresh: {}", tp);
            }
            if let Some(ov) = over_opt {
                println!("[INFO] Over  : {}", ov);
            }
            if let Some(sv) = sv_opt {
                println!("[INFO] SV    : {}", sv);
            }

            match mode.as_str() {
                "sep" => {
                    // Keep current behavior
                    let output_path: PathBuf = if let Some(out) = sub.get_one::<String>("output") {
                        PathBuf::from(out)
                    } else {
                        let seed = if use_multi {
                            qtl_str_or_csv
                                .split(',')
                                .next()
                                .unwrap_or("output")
                                .trim()
                                .to_string()
                        } else {
                            qtl_str_or_csv.clone()
                        };
                        derive_png_path(&seed)
                    };
                    let output_str = output_path.to_string_lossy().to_string();
                    println!("[INFO] Output: {}", output_str);
                    if use_multi {
                        man::run_multi(
                            &qtl_str_or_csv,
                            &output_str,
                            width,
                            height,
                            threads,
                            roof_opt,
                            thresh_opt,
                            no_tag,
                            sv_opt,
                        )
                    } else {
                        man::run_with_roof(
                            &qtl_str_or_csv,
                            &output_str,
                            width,
                            height,
                            threads,
                            roof_opt,
                            thresh_opt,
                            no_tag,
                            sv_opt,
                        )
                    }
                }
                "inter" => {
                    // Require multi with at least two files and a thresh
                    if !use_multi {
                        anyhow::bail!("--mode inter requires --multi with at least two files");
                    }
                    let files: Vec<String> = qtl_str_or_csv
                        .split(',')
                        .map(|s| s.trim().to_string())
                        .filter(|s| !s.is_empty())
                        .collect();
                    if files.len() < 2 {
                        anyhow::bail!("--mode inter requires --multi with at least two files");
                    }
                    if thresh_opt.is_none() {
                        anyhow::bail!("--mode inter requires --thresh to define significance");
                    }
                    // Single output file, same as sep derivation if not provided
                    let output_path: PathBuf = if let Some(out) = sub.get_one::<String>("output") {
                        PathBuf::from(out)
                    } else {
                        let seed = files
                            .first()
                            .cloned()
                            .unwrap_or_else(|| "output".to_string());
                        derive_png_path(&seed)
                    };
                    let output_str = output_path.to_string_lossy().to_string();
                    println!("[INFO] Output: {}", output_str);
                    // Call into man.rs specialized intersector
                    man::run_intersection(
                        &qtl_str_or_csv,
                        &output_str,
                        width,
                        height,
                        threads,
                        roof_opt,
                        thresh_opt,
                        no_tag,
                        sv_opt,
                    )
                }
                "compare" => {
                    // Require multi with at least two files
                    if !use_multi {
                        anyhow::bail!("--mode compare requires --multi with at least two files");
                    }
                    let files: Vec<String> = qtl_str_or_csv
                        .split(',')
                        .map(|s| s.trim().to_string())
                        .filter(|s| !s.is_empty())
                        .collect();
                    if files.len() < 2 {
                        anyhow::bail!("--mode compare requires --multi with at least two files");
                    }
                    // Determine output directory from -o (dir part), or default <first_qtl_dir>/manout
                    let outdir = derive_compare_outdir_from_o(
                        sub.get_one::<String>("output"),
                        files.first().unwrap(),
                    );
                    // Make sure directory exists
                    if let Err(e) = std::fs::create_dir_all(&outdir) {
                        anyhow::bail!(
                            "Failed to create output directory {}: {}",
                            outdir.to_string_lossy(),
                            e
                        );
                    }
                    println!("[INFO] OutDir: {}", outdir.to_string_lossy());
                    // Call into man.rs compare mode: expects directory path, will emit files per phenotype group
                    man::run_compare_overlay(
                        &qtl_str_or_csv,
                        &outdir.to_string_lossy(),
                        width,
                        height,
                        threads,
                        roof_opt,
                        thresh_opt,
                        over_opt,
                        no_tag,
                        sv_opt,
                    )
                }
                other => {
                    anyhow::bail!(format!("Unknown mode: {}", other));
                }
            }
        }
        Some(("vio", sub)) => {
            let multi_csv: String = sub
                .get_one::<String>("multi")
                .expect("--multi is required")
                .to_string();

            let width = *sub
                .get_one::<u32>("width")
                .expect("default provided by clap");
            let height = *sub
                .get_one::<u32>("height")
                .expect("default provided by clap");
            let threads = sub.get_one::<usize>("threads").copied();
            let roof_opt = sub.get_one::<u32>("roof").copied();
            let with_box = sub.get_flag("box");

            let output_path: PathBuf = if let Some(out) = sub.get_one::<String>("output") {
                PathBuf::from(out)
            } else {
                derive_violin_output_path(&multi_csv)
            };
            let output_str = output_path.to_string_lossy().to_string();

            println!("[INFO] MULTI : {}", multi_csv);
            println!("[INFO] Mode  : vio");
            println!("[INFO] Size  : {}x{}", width, height);
            if let Some(n) = threads {
                println!("[INFO] rayon threads = {}", n);
            }
            if let Some(roof) = roof_opt {
                println!("[INFO] Roof  : {}", roof);
            }
            if with_box {
                println!("[INFO] Box   : enabled");
            }
            println!("[INFO] Output: {}", output_str);

            vio::run_vio_multi(
                &multi_csv,
                &output_str,
                width,
                height,
                threads,
                roof_opt,
                with_box,
            )
        }
        Some(("qq", sub)) => {
            let multi_csv: String = sub
                .get_one::<String>("multi")
                .expect("--multi is required")
                .to_string();

            // parse simple width/height/threads/group
            let width = *sub
                .get_one::<u32>("width")
                .expect("default provided by clap");
            let height = *sub
                .get_one::<u32>("height")
                .expect("default provided by clap");
            let threads = sub.get_one::<usize>("threads").copied();
            let grouped = sub.get_flag("group");
            let qcut = sub.get_one::<f64>("qcut").copied();
            let diff = sub.get_flag("diff");
            let roof = sub.get_one::<f64>("roof").copied();

            // Decide output path
            let output_path: PathBuf = if let Some(out) = sub.get_one::<String>("output") {
                PathBuf::from(out)
            } else {
                derive_qq_output_path(&multi_csv)
            };
            let output_str = output_path.to_string_lossy().to_string();

            println!("[INFO] MULTI : {}", multi_csv);
            println!("[INFO] Mode  : qq");
            println!("[INFO] Size  : {}x{}", width, height);
            if let Some(n) = threads {
                println!("[INFO] rayon threads = {}", n);
            }
            println!(
                "[INFO] Group : {}",
                if grouped { "phenotype" } else { "genotype" }
            );
            if let Some(r) = roof {
                println!("[INFO] Roof  : {}", r);
            }
            println!("[INFO] Output: {}", output_str);

            qq::run_qq(
                &multi_csv,
                &output_str,
                width,
                height,
                threads,
                grouped,
                qcut,
                diff,
                roof,
            )
        }
        Some(("pro", sub)) => {
            let multi_csv: String = sub
                .get_one::<String>("multi")
                .expect("--multi is required")
                .to_string();

            let width = *sub
                .get_one::<u32>("width")
                .expect("default provided by clap");
            let height = *sub
                .get_one::<u32>("height")
                .expect("default provided by clap");
            let threads = sub.get_one::<usize>("threads").copied();
            let thresh = *sub
                .get_one::<f64>("thresh")
                .expect("default provided by clap");
            let binw = *sub.get_one::<f64>("bin").expect("default provided by clap");
            let diff = sub.get_one::<f64>("diff").copied();
            let venn = sub.get_flag("venn");

            let output_path: PathBuf = if let Some(out) = sub.get_one::<String>("output") {
                PathBuf::from(out)
            } else {
                derive_pro_output_path(&multi_csv)
            };
            let output_str = output_path.to_string_lossy().to_string();

            println!("[INFO] MULTI : {}", multi_csv);
            println!("[INFO] Mode  : pro");
            println!("[INFO] Size  : {}x{}", width, height);
            if let Some(n) = threads {
                println!("[INFO] rayon threads = {}", n);
            }
            println!("[INFO] Thresh: {}", thresh);
            println!("[INFO] Bin   : {}%", binw);
            if let Some(d) = diff {
                println!("[INFO] Diff  : {}", d);
            }
            if venn {
                println!("[INFO] Venn  : enabled");
            }
            println!("[INFO] Output: {}", output_str);

            pro::run_pro(
                &multi_csv,
                &output_str,
                width,
                height,
                threads,
                thresh,
                Some(binw),
                diff,
                venn,
            )
        }
        _ => {
            // No subcommand: print help
            let _ = Command::new("eQTLplotter").print_help();
            println!();
            Ok(())
        }
    }
}

fn derive_png_path(qtl_path: &str) -> PathBuf {
    let p = Path::new(qtl_path);
    let parent = p.parent().unwrap_or_else(|| Path::new("."));
    let stem = p
        .file_stem()
        .unwrap_or_else(|| std::ffi::OsStr::new("output"));
    parent.join(format!("{}.png", stem.to_string_lossy()))
}

fn is_path_dir_hint(p: &Path) -> bool {
    // If it exists and is a dir, true; if it ends with a path separator or has no file stem, treat as dir hint
    if p.exists() {
        return p.is_dir();
    }
    let s = p.to_string_lossy();
    s.ends_with('/') || p.file_name().is_none()
}

fn derive_compare_outdir_from_o(o_opt: Option<&String>, seed_first_qtl: &str) -> PathBuf {
    if let Some(o) = o_opt {
        let p = PathBuf::from(o);
        if is_path_dir_hint(&p) {
            return p;
        } else {
            // Use the parent of the provided file path as directory
            return p.parent().unwrap_or_else(|| Path::new(".")).to_path_buf();
        }
    }
    // default: <first_qtl_dir>/manout
    let seed = Path::new(seed_first_qtl);
    let parent = seed.parent().unwrap_or_else(|| Path::new("."));
    parent.join("manout")
}

fn derive_violin_output_path(multi_csv: &str) -> PathBuf {
    let files: Vec<&str> = multi_csv
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect();
    let first = files.first().cloned().unwrap_or("output");
    let first_p = Path::new(first);
    let parent = first_p.parent().unwrap_or_else(|| Path::new("."));

    // Collect stems for all files
    let mut stems: Vec<String> = Vec::with_capacity(files.len());
    for f in &files {
        let stem = Path::new(f)
            .file_stem()
            .and_then(|s| Some(s.to_string_lossy().to_string()))
            .unwrap_or_else(|| f.to_string());
        stems.push(stem);
    }
    let name = format!("violin_{}.png", stems.join("_"));
    parent.join(name)
}

fn derive_qq_output_path(multi_csv: &str) -> PathBuf {
    let files: Vec<&str> = multi_csv
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect();
    let (a, b) = match (files.get(0), files.get(1)) {
        (Some(a), Some(b)) => (*a, *b),
        _ => ("outputA", "outputB"),
    };
    let a_p = Path::new(a);
    let parent = a_p.parent().unwrap_or_else(|| Path::new("."));
    let stem_a = a_p
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| a.to_string());
    let stem_b = Path::new(b)
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| b.to_string());
    parent.join(format!("qq_{}_vs_{}.png", stem_a, stem_b))
}

fn derive_pro_output_path(multi_csv: &str) -> PathBuf {
    let files: Vec<&str> = multi_csv
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .collect();
    let (a, b) = match (files.get(0), files.get(1)) {
        (Some(a), Some(b)) => (*a, *b),
        _ => ("outputA", "outputB"),
    };
    let a_p = Path::new(a);
    let parent = a_p.parent().unwrap_or_else(|| Path::new("."));
    let stem_a = a_p
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| a.to_string());
    let stem_b = Path::new(b)
        .file_stem()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| b.to_string());
    parent.join(format!("pro_{}_vs_{}.png", stem_a, stem_b))
}
