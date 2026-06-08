use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::PathBuf;

#[derive(Debug, Clone)]
pub struct FusionEvent {
    pub gene_a: String,
    pub gene_b: String,
    pub transcript_a: String,
    pub transcript_b: String,
    pub breakpoint_pos_a: usize,
    pub breakpoint_pos_b: usize,
    pub supporting_reads: usize,
    pub split_reads: usize,
    pub discordant_pairs: usize,
    pub fusion_type: String, // "known" or "novel"
    pub known_fusion: Option<String>, // Name of known fusion if applicable
}

#[derive(Debug)]
pub struct KnownFusion {
    pub name: String,
    pub gene_a: String,
    pub gene_b: String,
    pub cancer_type: String,
    pub clinical_significance: String,
}

pub fn get_known_fusions() -> Vec<KnownFusion> {
    vec![
        KnownFusion {
            name: "BCR-ABL1".to_string(),
            gene_a: "BCR".to_string(),
            gene_b: "ABL1".to_string(),
            cancer_type: "CML".to_string(),
            clinical_significance: "Imatinib target".to_string(),
        },
        KnownFusion {
            name: "EML4-ALK".to_string(),
            gene_a: "EML4".to_string(),
            gene_b: "ALK".to_string(),
            cancer_type: "Lung adenocarcinoma".to_string(),
            clinical_significance: "Crizotinib target".to_string(),
        },
        KnownFusion {
            name: "TMPRSS2-ERG".to_string(),
            gene_a: "TMPRSS2".to_string(),
            gene_b: "ERG".to_string(),
            cancer_type: "Prostate cancer".to_string(),
            clinical_significance: "Most common prostate fusion".to_string(),
        },
        KnownFusion {
            name: "FLT3-ITD".to_string(),
            gene_a: "FLT3".to_string(),
            gene_b: "FLT3".to_string(),
            cancer_type: "AML".to_string(),
            clinical_significance: "Poor prognosis".to_string(),
        },
        KnownFusion {
            name: "KMT2A-MLL".to_string(),
            gene_a: "KMT2A".to_string(),
            gene_b: "MLL".to_string(),
            cancer_type: "AML/ALL".to_string(),
            clinical_significance: "Poor prognosis".to_string(),
        },
        KnownFusion {
            name: "NPM1c".to_string(),
            gene_a: "NPM1".to_string(),
            gene_b: "NPM1".to_string(),
            cancer_type: "AML".to_string(),
            clinical_significance: "Favorable prognosis".to_string(),
        },
        KnownFusion {
            name: "PML-RARA".to_string(),
            gene_a: "PML".to_string(),
            gene_b: "RARA".to_string(),
            cancer_type: "APL".to_string(),
            clinical_significance: "ATRA target".to_string(),
        },
        KnownFusion {
            name: "RET-PTC".to_string(),
            gene_a: "RET".to_string(),
            gene_b: "PTC".to_string(),
            cancer_type: "Thyroid cancer".to_string(),
            clinical_significance: "Kinase inhibitor target".to_string(),
        },
        KnownFusion {
            name: "ROS1".to_string(),
            gene_a: "ROS1".to_string(),
            gene_b: "ROS1".to_string(),
            cancer_type: "Lung cancer".to_string(),
            clinical_significance: "Crizotinib target".to_string(),
        },
        KnownFusion {
            name: "NTRK1".to_string(),
            gene_a: "NTRK1".to_string(),
            gene_b: "NTRK1".to_string(),
            cancer_type: "Multiple cancers".to_string(),
            clinical_significance: "Larotrectinib target".to_string(),
        },
    ]
}

fn extract_gene_from_transcript(transcript: &str) -> String {
    // Extract gene symbol from transcript header if available
    if transcript.contains("gene_symbol:") {
        let parts: Vec<&str> = transcript.split("gene_symbol:").collect();
        if parts.len() > 1 {
            let gene = parts[1].split(' ').next().unwrap_or("");
            if !gene.is_empty() {
                return gene.to_string();
            }
        }
    }
    
    // Fallback: use transcript ID
    transcript.split(' ').next().unwrap_or(transcript).to_string()
}

fn is_discordant_pair(
    read1_transcript: &str,
    read2_transcript: &str,
    _read1_pos: usize,
    _read2_pos: usize,
    insert_size: isize,
) -> bool {
    // Discordant if mapped to different transcripts
    read1_transcript != read2_transcript
        // Or if insert size is abnormal (>10kb or <-1kb)
        || insert_size.abs() > 10000
}

fn is_split_read(cigar: &str, _seq: &str) -> bool {
    // Check for soft-clipping that might indicate split read
    cigar.starts_with('S') || cigar.ends_with('S') || cigar.contains("S")
}

pub fn detect_fusions(
    sam_path: &PathBuf,
    output: Option<&PathBuf>,
    min_supporting_reads: usize,
) -> Result<Vec<FusionEvent>, Box<dyn std::error::Error>> {
    let file = File::open(sam_path)?;
    let reader = BufReader::new(file);

    let mut fusion_counts: HashMap<(String, String), FusionEvent> = HashMap::new();
    let mut read_pairs: HashMap<String, (String, usize, isize)> = HashMap::new();
    let known_fusions = get_known_fusions();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('@') {
            continue;
        }

        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 11 {
            continue;
        }

        let read_name = fields[0].to_string();
        let flag: u16 = fields[1].parse().unwrap_or(0);
        let transcript = fields[2].to_string();
        let pos: usize = fields[3].parse().unwrap_or(0);
        let cigar = fields[5].to_string();
        let _mate_transcript = fields[6].to_string();
        let _mate_pos: usize = fields[7].parse().unwrap_or(0);
        let seq = fields[9].to_string();

        let gene = extract_gene_from_transcript(&transcript);

        // Check for split reads
        if is_split_read(&cigar, &seq) && !transcript.eq("*") {
            // This read might be part of a fusion
            // For now, just count soft-clipped reads as potential fusion evidence
            let key = (gene.clone(), gene.clone());
            let entry = fusion_counts.entry(key).or_insert_with(|| FusionEvent {
                gene_a: gene.clone(),
                gene_b: gene.clone(),
                transcript_a: transcript.clone(),
                transcript_b: transcript.clone(),
                breakpoint_pos_a: pos,
                breakpoint_pos_b: 0,
                supporting_reads: 0,
                split_reads: 0,
                discordant_pairs: 0,
                fusion_type: "novel".to_string(),
                known_fusion: None,
            });
            entry.supporting_reads += 1;
            entry.split_reads += 1;
        }

        // Track read pairs for discordant pair detection
        if (flag & 1) == 1 { // Proper pair
            read_pairs.insert(
                read_name.clone(),
                (transcript.clone(), pos, 0), // Insert size would be calculated from TLEN
            );
        }
    }

    // Filter by minimum supporting reads
    let mut fusions: Vec<FusionEvent> = fusion_counts
        .into_values()
        .filter(|f| f.supporting_reads >= min_supporting_reads)
        .collect();

    // Check against known fusions
    for fusion in &mut fusions {
        for known in &known_fusions {
            if (fusion.gene_a == known.gene_a && fusion.gene_b == known.gene_b)
                || (fusion.gene_a == known.gene_b && fusion.gene_b == known.gene_a)
            {
                fusion.fusion_type = "known".to_string();
                fusion.known_fusion = Some(known.name.clone());
                break;
            }
        }
    }

    // Sort by supporting reads
    fusions.sort_by(|a, b| b.supporting_reads.cmp(&a.supporting_reads));

    // Print results
    println!("=== Fusion Detection Results ===");
    println!("Total potential fusions: {}", fusions.len());
    println!("Known fusions: {}", fusions.iter().filter(|f| f.fusion_type == "known").count());
    println!("Novel fusions: {}", fusions.iter().filter(|f| f.fusion_type == "novel").count());

    if !fusions.is_empty() {
        println!("\nTop fusion events:");
        println!("{:<20} {:<15} {:<15} {:<15} {:<12} {:<10}",
            "Fusion", "Gene A", "Gene B", "Type", "Supporting", "Split");
        for fusion in fusions.iter().take(20) {
            let name = fusion.known_fusion.as_deref().unwrap_or("novel");
            println!("{:<20} {:<15} {:<15} {:<15} {:<12} {:<10}",
                name,
                fusion.gene_a,
                fusion.gene_b,
                fusion.fusion_type,
                fusion.supporting_reads,
                fusion.split_reads);
        }
    }

    // Write to file if specified
    if let Some(out_path) = output {
        let mut file = File::create(out_path)?;
        writeln!(file, "fusion_name\tgene_a\tgene_b\ttranscript_a\ttranscript_b\ttype\tsupporting_reads\tsplit_reads\tdiscordant_pairs")?;
        for fusion in &fusions {
            let name = fusion.known_fusion.as_deref().unwrap_or("novel");
            writeln!(file, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                name,
                fusion.gene_a,
                fusion.gene_b,
                fusion.transcript_a,
                fusion.transcript_b,
                fusion.fusion_type,
                fusion.supporting_reads,
                fusion.split_reads,
                fusion.discordant_pairs)?;
        }
        println!("\nTSV saved to {}", out_path.display());
    }

    Ok(fusions)
}

pub fn generate_fusion_report(
    output_path: &PathBuf,
    fusions: &[FusionEvent],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut file = File::create(output_path)?;

    let known_count = fusions.iter().filter(|f| f.fusion_type == "known").count();
    let novel_count = fusions.iter().filter(|f| f.fusion_type == "novel").count();

    // HTML header
    writeln!(file, "<!DOCTYPE html>")?;
    writeln!(file, "<html>")?;
    writeln!(file, "<head>")?;
    writeln!(file, "  <title>Fusion Detection Report</title>")?;
    writeln!(file, "  <style>")?;
    writeln!(file, "    body {{ font-family: Arial, sans-serif; margin: 20px; }}")?;
    writeln!(file, "    h1 {{ color: #2c3e50; }}")?;
    writeln!(file, "    h2 {{ color: #34495e; border-bottom: 2px solid #3498db; }}")?;
    writeln!(file, "    table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}")?;
    writeln!(file, "    th, td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}")?;
    writeln!(file, "    th {{ background-color: #3498db; color: white; }}")?;
    writeln!(file, "    .known {{ background-color: #f39c12; color: white; padding: 2px 6px; border-radius: 3px; }}")?;
    writeln!(file, "    .novel {{ background-color: #95a5a6; color: white; padding: 2px 6px; border-radius: 3px; }}")?;
    writeln!(file, "    .stats {{ display: flex; gap: 20px; margin: 20px 0; }}")?;
    writeln!(file, "    .stat-box {{ background-color: #ecf0f1; padding: 15px; border-radius: 5px; flex: 1; }}")?;
    writeln!(file, "    .stat-value {{ font-size: 24px; font-weight: bold; color: #2c3e50; }}")?;
    writeln!(file, "  </style>")?;
    writeln!(file, "</head>")?;
    writeln!(file, "<body>")?;

    // Title
    writeln!(file, "<h1>Fusion Detection Report</h1>")?;
    writeln!(file, "<p>Generated by RNA-Pop</p>")?;

    // Summary stats
    writeln!(file, "<div class='stats'>")?;
    writeln!(file, "  <div class='stat-box'><div class='stat-value'>{}</div>Total Fusions</div>", fusions.len())?;
    writeln!(file, "  <div class='stat-box'><div class='stat-value'>{}</div>Known Fusions</div>", known_count)?;
    writeln!(file, "  <div class='stat-box'><div class='stat-value'>{}</div>Novel Fusions</div>", novel_count)?;
    writeln!(file, "</div>")?;

    // Fusion table
    writeln!(file, "<h2>Fusion Events</h2>")?;
    writeln!(file, "<table>")?;
    writeln!(file, "  <tr><th>Fusion</th><th>Gene A</th><th>Gene B</th><th>Type</th><th>Supporting Reads</th><th>Split Reads</th></tr>")?;

    for fusion in fusions {
        let type_class = if fusion.fusion_type == "known" { "known" } else { "novel" };
        writeln!(file, "  <tr>")?;
        writeln!(file, "    <td>{}</td>", fusion.known_fusion.as_deref().unwrap_or("novel"))?;
        writeln!(file, "    <td>{}</td>", fusion.gene_a)?;
        writeln!(file, "    <td>{}</td>", fusion.gene_b)?;
        writeln!(file, "    <td><span class='{}'>{}</span></td>", type_class, fusion.fusion_type)?;
        writeln!(file, "    <td>{}</td>", fusion.supporting_reads)?;
        writeln!(file, "    <td>{}</td>", fusion.split_reads)?;
        writeln!(file, "  </tr>")?;
    }

    writeln!(file, "</table>")?;

    // Footer
    writeln!(file, "<hr>")?;
    writeln!(file, "<p><small>Report generated by RNA-Pop. For research use only. Not for clinical diagnosis.</small></p>")?;
    writeln!(file, "</body>")?;
    writeln!(file, "</html>")?;

    Ok(())
}
