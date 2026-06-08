use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

#[derive(Debug)]
pub struct DiffResult {
    pub transcript: String,
    pub abundance_a: f64,
    pub abundance_b: f64,
    pub fold_change: f64,
    pub log2_fc: f64,
    pub p_value: f64,
    pub significant: bool,
    pub direction: String,
}

fn read_abundances(path: &PathBuf) -> Result<HashMap<String, f64>, Box<dyn std::error::Error>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut abundances = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        if i == 0 {
            continue;
        }
        let line = line?;
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            let transcript = parts[0].trim().to_string();
            let abundance: f64 = parts[1].trim().parse().unwrap_or(0.0);
            abundances.insert(transcript, abundance);
        }
    }

    Ok(abundances)
}

fn fisher_exact_pvalue(a: f64, b: f64) -> f64 {
    // Simplified p-value estimation based on fold change and abundance
    let total = a + b;
    if total < 1e-10 {
        return 1.0;
    }

    let prop_a = a / total;
    let prop_b = b / total;
    let diff = (prop_a - prop_b).abs();

    // Approximate p-value using normal approximation
    let se = ((prop_a * (1.0 - prop_a) + prop_b * (1.0 - prop_b)) / total).sqrt();
    if se < 1e-10 {
        return 1.0;
    }

    let z = diff / se;
    // Normal CDF approximation
    let p = 2.0 * (1.0 - normal_cdf(z.abs()));
    p.min(1.0).max(0.0)
}

fn normal_cdf(x: f64) -> f64 {
    // Approximation of standard normal CDF
    let a1 = 0.254829592;
    let a2 = -0.284496736;
    let a3 = 1.421413741;
    let a4 = -1.453152027;
    let a5 = 1.061405429;
    let p = 0.3275911;

    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs() / 2.0_f64.sqrt();

    let t = 1.0 / (1.0 + p * x);
    let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * (-x * x).exp();

    0.5 * (1.0 + sign * y)
}

pub fn compare_samples(
    path_a: &PathBuf,
    path_b: &PathBuf,
    output: Option<&PathBuf>,
    fc_threshold: f64,
    p_threshold: f64,
) -> Result<Vec<DiffResult>, Box<dyn std::error::Error>> {
    let abundances_a = read_abundances(path_a)?;
    let abundances_b = read_abundances(path_b)?;

    let mut results = Vec::new();

    // Get all transcripts
    let mut all_transcripts: Vec<String> = Vec::new();
    for tx in abundances_a.keys() {
        all_transcripts.push(tx.clone());
    }
    for tx in abundances_b.keys() {
        if !all_transcripts.contains(tx) {
            all_transcripts.push(tx.clone());
        }
    }

    for transcript in &all_transcripts {
        let ab = abundances_a.get(transcript).copied().unwrap_or(0.0);
        let bb = abundances_b.get(transcript).copied().unwrap_or(0.0);

        // Add pseudocount to avoid division by zero
        let ab_adj = ab + 1e-10;
        let bb_adj = bb + 1e-10;

        let fold_change = if bb_adj > 1e-9 { ab_adj / bb_adj } else { 0.0 };
        let log2_fc = fold_change.log2();
        let p_value = fisher_exact_pvalue(ab, bb);
        let significant = p_value < p_threshold && fold_change > fc_threshold || 1.0 / fold_change > fc_threshold;

        let direction = if log2_fc > 0.585 { // ~1.5x
            "UP".to_string()
        } else if log2_fc < -0.585 { // ~0.67x
            "DOWN".to_string()
        } else {
            "UNCHANGED".to_string()
        };

        results.push(DiffResult {
            transcript: (*transcript).clone(),
            abundance_a: ab,
            abundance_b: bb,
            fold_change,
            log2_fc,
            p_value,
            significant,
            direction,
        });
    }

    // Sort by p-value
    results.sort_by(|a, b| a.p_value.partial_cmp(&b.p_value).unwrap());

    // Print summary
    let sig_count = results.iter().filter(|r| r.significant).count();
    println!("=== Differential Expression Analysis ===");
    println!("Sample A: {}", path_a.display());
    println!("Sample B: {}", path_b.display());
    println!("Total transcripts: {}", results.len());
    println!("Significant (FC>{}, p<{}): {}", fc_threshold, p_threshold, sig_count);

    // Print top hits
    println!("\nTop 10 most significant:");
    println!("{:<30} {:>12} {:>12} {:>10} {:>10} {:>10}",
        "Transcript", "Abundance A", "Abundance B", "Fold Change", "Log2FC", "p-value");
    for result in results.iter().take(10) {
        println!("{:<30} {:>12.6} {:>12.6} {:>10.2} {:>10.2} {:>10.4}",
            result.transcript.chars().take(30).collect::<String>(),
            result.abundance_a,
            result.abundance_b,
            result.fold_change,
            result.log2_fc,
            result.p_value);
    }

    // Write to file if specified
    if let Some(out_path) = output {
        let mut file = File::create(out_path)?;
        writeln!(file, "transcript\tabundance_a\tabundance_b\tfold_change\tlog2_fc\tp_value\tsignificant\tdirection")?;
        for result in &results {
            writeln!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                result.transcript,
                result.abundance_a,
                result.abundance_b,
                result.fold_change,
                result.log2_fc,
                result.p_value,
                result.significant,
                result.direction)?;
        }
        println!("\nTSV saved to {}", out_path.display());
    }

    Ok(results)
}

pub fn generate_comparison_report(
    output_path: &PathBuf,
    results: &[DiffResult],
    fc_threshold: f64,
    p_threshold: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(output_path)?;

    let sig_count = results.iter().filter(|r| r.significant).count();
    let up_count = results.iter().filter(|r| r.direction == "UP").count();
    let down_count = results.iter().filter(|r| r.direction == "DOWN").count();

    // HTML header
    writeln!(file, "<!DOCTYPE html>")?;
    writeln!(file, "<html>")?;
    writeln!(file, "<head>")?;
    writeln!(file, "  <title>Differential Expression Report</title>")?;
    writeln!(file, "  <style>")?;
    writeln!(file, "    body {{ font-family: Arial, sans-serif; margin: 20px; }}")?;
    writeln!(file, "    h1 {{ color: #2c3e50; }}")?;
    writeln!(file, "    h2 {{ color: #34495e; border-bottom: 2px solid #3498db; }}")?;
    writeln!(file, "    table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}")?;
    writeln!(file, "    th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}")?;
    writeln!(file, "    th {{ background-color: #3498db; color: white; }}")?;
    writeln!(file, "    .up {{ color: #e74c3c; font-weight: bold; }}")?;
    writeln!(file, "    .down {{ color: #27ae60; font-weight: bold; }}")?;
    writeln!(file, "    .unchanged {{ color: #7f8c8d; }}")?;
    writeln!(file, "    .significant {{ background-color: #f39c12; color: white; padding: 2px 6px; border-radius: 3px; }}")?;
    writeln!(file, "    .stats {{ display: flex; gap: 20px; margin: 20px 0; }}")?;
    writeln!(file, "    .stat-box {{ background-color: #ecf0f1; padding: 15px; border-radius: 5px; flex: 1; }}")?;
    writeln!(file, "    .stat-value {{ font-size: 24px; font-weight: bold; color: #2c3e50; }}")?;
    writeln!(file, "  </style>")?;
    writeln!(file, "</head>")?;
    writeln!(file, "<body>")?;

    // Title
    writeln!(file, "<h1>Differential Expression Analysis</h1>")?;
    writeln!(file, "<p>Generated by RNA-Pop</p>")?;

    // Summary stats
    writeln!(file, "<div class='stats'>")?;
    writeln!(file, "  <div class='stat-box'><div class='stat-value'>{}%</div>Significant</div>",
        sig_count as f64 / results.len() as f64 * 100.0)?;
    writeln!(file, "  <div class='stat-box'><div class='stat-value up'>{}</div>Up-regulated</div>", up_count)?;
    writeln!(file, "  <div class='stat-box'><div class='stat-value down'>{}</div>Down-regulated</div>", down_count)?;
    writeln!(file, "  <div class='stat-box'><div class='stat-value'>{}%</div>Total transcripts</div>", results.len())?;
    writeln!(file, "</div>")?;

    // Significant genes table
    writeln!(file, "<h2>Significant Differential Expression (FC>{}, p<{})</h2>", fc_threshold, p_threshold)?;
    writeln!(file, "<table>")?;
    writeln!(file, "  <tr><th>Transcript</th><th>Abundance A</th><th>Abundance B</th><th>Fold Change</th><th>Log2FC</th><th>p-value</th><th>Direction</th></tr>")?;

    for result in results.iter().filter(|r| r.significant).take(100) {
        let dir_class = match result.direction.as_str() {
            "UP" => "up",
            "DOWN" => "down",
            _ => "unchanged",
        };
        writeln!(file, "  <tr>")?;
        writeln!(file, "    <td>{}</td>", result.transcript)?;
        writeln!(file, "    <td>{:.6}</td>", result.abundance_a)?;
        writeln!(file, "    <td>{:.6}</td>", result.abundance_b)?;
        writeln!(file, "    <td>{:.2}</td>", result.fold_change)?;
        writeln!(file, "    <td>{:.2}</td>", result.log2_fc)?;
        writeln!(file, "    <td>{:.4}</td>", result.p_value)?;
        writeln!(file, "    <td class='{}'><span class='significant'>{}</span></td>", dir_class, result.direction)?;
        writeln!(file, "  </tr>")?;
    }

    writeln!(file, "</table>")?;

    // Footer
    writeln!(file, "<hr>")?;
    writeln!(file, "<p><small>Report generated by RNA-Pop. For research use only.</small></p>")?;
    writeln!(file, "</body>")?;
    writeln!(file, "</html>")?;

    Ok(())
}
