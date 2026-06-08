use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

#[derive(Debug)]
pub struct QCMetrics {
    pub total_reads: usize,
    pub mapped_reads: usize,
    pub mapping_rate: f64,
    pub multi_mapped: usize,
    pub unique_mapped: usize,
    pub avg_mapping_quality: f64,
    pub avg_read_length: f64,
    pub transcripts_with_coverage: usize,
    pub total_transcripts: usize,
    pub coverage_uniformity: f64,
}

pub fn run_qc(sam_path: &PathBuf, output: Option<&PathBuf>) -> Result<QCMetrics, Box<dyn std::error::Error>> {
    let file = File::open(sam_path)?;
    let reader = BufReader::new(file);

    let mut total_reads = 0usize;
    let mut mapped_reads = 0usize;
    let mut multi_mapped = 0usize;
    let mut unique_mapped = 0usize;
    let mut total_mapping_quality = 0f64;
    let mut total_read_length = 0f64;
    let mut transcript_coverage: std::collections::HashMap<String, usize> = std::collections::HashMap::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') {
            if line.starts_with("@SQ") {
                // Count transcripts from SAM header
                let parts: Vec<&str> = line.split('\t').collect();
                for part in parts {
                    if part.starts_with("SN:") {
                        let _ = transcript_coverage.entry(part[3..].to_string()).or_insert(0);
                    }
                }
            }
            continue;
        }

        total_reads += 1;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }

        let flag: u16 = fields[1].parse().unwrap_or(0);
        let is_mapped = (flag & 4) == 0;

        if is_mapped {
            mapped_reads += 1;
            let mapq: i32 = fields[4].parse().unwrap_or(0);
            total_mapping_quality += mapq as f64;

            let seq = fields[9];
            total_read_length += seq.len() as f64;

            // Count multi-mapping
            let nm: i32 = fields.iter()
                .position(|&f| f.starts_with("NM:i:"))
                .map(|i| fields[i][4..].parse().unwrap_or(0))
                .unwrap_or(0);

            if nm > 0 {
                multi_mapped += 1;
            } else {
                unique_mapped += 1;
            }

            // Track transcript coverage
            let transcript = fields[2].to_string();
            if !transcript.eq("*") {
                *transcript_coverage.entry(transcript).or_insert(0) += 1;
            }
        }
    }

    let total_transcripts = transcript_coverage.len();
    let transcripts_with_coverage = transcript_coverage.values().filter(|&&v| v > 0).count();

    // Calculate coverage uniformity (coefficient of variation)
    let coverages: Vec<usize> = transcript_coverage.values().filter(|&&v| v > 0).cloned().collect();
    let coverage_uniformity = if coverages.len() > 1 {
        let mean = coverages.iter().sum::<usize>() as f64 / coverages.len() as f64;
        let variance = coverages.iter().map(|&x| (x as f64 - mean).powi(2)).sum::<f64>() / coverages.len() as f64;
        let std_dev = variance.sqrt();
        if mean > 0.0 {
            1.0 - (std_dev / mean).min(1.0)
        } else {
            0.0
        }
    } else {
        1.0
    };

    let metrics = QCMetrics {
        total_reads,
        mapped_reads,
        mapping_rate: if total_reads > 0 { mapped_reads as f64 / total_reads as f64 * 100.0 } else { 0.0 },
        multi_mapped,
        unique_mapped,
        avg_mapping_quality: if mapped_reads > 0 { total_mapping_quality / mapped_reads as f64 } else { 0.0 },
        avg_read_length: if total_reads > 0 { total_read_length / total_reads as f64 } else { 0.0 },
        transcripts_with_coverage,
        total_transcripts,
        coverage_uniformity,
    };

    // Print metrics
    println!("=== RNA-Pop QC Report ===");
    println!("Total reads:            {}", metrics.total_reads);
    println!("Mapped reads:           {} ({:.1}%)", metrics.mapped_reads, metrics.mapping_rate);
    println!("Unique mapped:          {}", metrics.unique_mapped);
    println!("Multi-mapped:           {}", metrics.multi_mapped);
    println!("Avg mapping quality:    {:.1}", metrics.avg_mapping_quality);
    println!("Avg read length:        {:.0} bp", metrics.avg_read_length);
    println!("Transcripts covered:    {}/{}", metrics.transcripts_with_coverage, metrics.total_transcripts);
    println!("Coverage uniformity:    {:.2}", metrics.coverage_uniformity);

    // Write to file if specified
    if let Some(out_path) = output {
        let mut file = File::create(out_path)?;
        writeln!(file, "=== RNA-Pop QC Report ===")?;
        writeln!(file, "Total reads:            {}", metrics.total_reads)?;
        writeln!(file, "Mapped reads:           {} ({:.1}%)", metrics.mapped_reads, metrics.mapping_rate)?;
        writeln!(file, "Unique mapped:          {}", metrics.unique_mapped)?;
        writeln!(file, "Multi-mapped:           {}", metrics.multi_mapped)?;
        writeln!(file, "Avg mapping quality:    {:.1}", metrics.avg_mapping_quality)?;
        writeln!(file, "Avg read length:        {:.0} bp", metrics.avg_read_length)?;
        writeln!(file, "Transcripts covered:    {}/{}", metrics.transcripts_with_coverage, metrics.total_transcripts)?;
        writeln!(file, "Coverage uniformity:    {:.2}", metrics.coverage_uniformity)?;
    }

    Ok(metrics)
}
