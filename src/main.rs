use clap::{Parser, Subcommand};
use std::collections::HashMap;
use std::io::BufRead;
use std::path::PathBuf;
use std::time::Instant;

use rna_pop::chunk_consensus::MultiChunkConsensus;
use rna_pop::consensus::MultiKConsensus;
use rna_pop::fastcon::FastCon;
use rna_pop::fastq::parse_fastq;
use rna_pop::qc::run_qc;
use rna_pop::report::generate_clinical_report;
use rna_pop::compare::{compare_samples, generate_comparison_report};
use rna_pop::{AlignMode, BitPop};

#[derive(Parser)]
#[command(name = "rna-pop", about = "RNA-seq read quantification: transcript mapping, splice-aware alignment, EM-based abundance estimation")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// One-command workflow: build index + map reads
    Run(RunArgs),
    /// Build FM-Index from transcript FASTA file(s)
    Build(BuildArgs),
    /// Map reads to indexed transcripts
    Map(MapArgs),
    /// Show index statistics
    Stats(StatsArgs),
    /// Apply EM algorithm for soft-assignment transcript quantification
    Em(EmArgs),
    /// Multi-k consensus: map reads against multiple k-indexes with voting
    Consensus(ConsensusArgs),
    /// Fast consensus: run `rna-pop map` for each index, then combine
    FastCon(FastConArgs),
    /// Multi chunk-% consensus: same index, different chunk sizes, voting
    ChunkConsensus(ChunkConsensusArgs),
    /// Quality control metrics for SAM file
    Qc(QcArgs),
    /// Generate clinical report with biomarker analysis
    Report(ReportArgs),
    /// Compare two samples for differential expression
    Compare(CompareArgs),
}

// --- Run ---

#[derive(clap::Args)]
struct RunArgs {
    /// Transcriptome source: FASTA file or folder of FASTA files
    #[arg(short, long)]
    genome: Option<String>,

    /// Use existing index file (instead of building from transcripts)
    #[arg(short, long)]
    index: Option<PathBuf>,

    /// Reads file (FASTQ) for single-end mode
    #[arg(short, long)]
    reads: Option<PathBuf>,

    /// R1 FASTQ file for paired-end mapping
    #[arg(short = '1', long)]
    reads_1: Option<PathBuf>,

    /// R2 FASTQ file for paired-end mapping
    #[arg(short = '2', long)]
    reads_2: Option<PathBuf>,

    /// K-mer size (default: auto)
    #[arg(short, long)]
    k: Option<usize>,

    /// Number of parallel threads
    #[arg(short, long)]
    threads: Option<usize>,

    /// Output SAM file path
    #[arg(short, long)]
    sam: Option<PathBuf>,

    /// Alignment mode: xor, sw, hybrid, softclip, chain
    #[arg(long, default_value = "xor")]
    align: String,

    /// Number of top rarest k-mers to try as anchors
    #[arg(long, default_value_t = 1)]
    top_n: usize,

    /// Read type: short or long
    #[arg(long, default_value = "short")]
    read_type: String,
}

// --- Build ---

#[derive(clap::Args)]
struct BuildArgs {
    /// Transcriptome source: FASTA file or folder of FASTA files
    #[arg(short, long)]
    genome: Option<String>,

    /// K-mer size (default: auto)
    #[arg(short, long)]
    k: Option<usize>,

    /// Output index file path
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Number of parallel threads
    #[arg(short, long)]
    threads: Option<usize>,
}

// --- Map ---

#[derive(clap::Args)]
struct MapArgs {
    /// Index file
    #[arg(short, long)]
    index: PathBuf,

    /// Reads file (FASTQ)
    #[arg(short, long)]
    reads: Option<PathBuf>,

    /// R1 FASTQ file for paired-end mapping
    #[arg(short = '1', long)]
    reads_1: Option<PathBuf>,

    /// R2 FASTQ file for paired-end mapping
    #[arg(short = '2', long)]
    reads_2: Option<PathBuf>,

    /// Number of parallel threads
    #[arg(short, long)]
    threads: Option<usize>,

    /// Output SAM file path
    #[arg(short, long)]
    sam: Option<PathBuf>,

    /// Alignment mode: xor, sw, hybrid, softclip, chain
    #[arg(long, default_value = "xor")]
    align: String,

    /// Number of top rarest k-mers to try as anchors
    #[arg(long, default_value_t = 1)]
    top_n: usize,

    /// Read type: short or long
    #[arg(long, default_value = "short")]
    read_type: String,
}

// --- Stats ---

#[derive(clap::Args)]
struct StatsArgs {
    /// Index file
    #[arg(short, long)]
    index: PathBuf,
}

// --- EM ---

#[derive(clap::Args)]
struct EmArgs {
    /// SAM file with mapping results
    #[arg(short, long)]
    sam: PathBuf,

    /// Index file for genome information
    #[arg(short, long)]
    index: Option<PathBuf>,

    /// EM iterations
    #[arg(long, default_value_t = 20)]
    iterations: usize,

    /// EM convergence threshold
    #[arg(long, default_value_t = 0.001)]
    threshold: f64,

    /// Output TSV file with abundances
    #[arg(short, long)]
    tsv: Option<PathBuf>,
}

// --- Consensus ---

#[derive(clap::Args)]
struct ConsensusArgs {
    /// Index file
    #[arg(short, long)]
    index: PathBuf,

    /// Reads file (FASTQ)
    #[arg(short, long)]
    reads: Option<PathBuf>,

    /// R1 FASTQ file for paired-end mapping
    #[arg(short = '1', long)]
    reads_1: Option<PathBuf>,

    /// R2 FASTQ file for paired-end mapping
    #[arg(short = '2', long)]
    reads_2: Option<PathBuf>,

    /// K-mer sizes to try (comma-separated, e.g. "8,10,12")
    #[arg(long, default_value = "8,10,12")]
    kmers: String,

    /// Number of parallel threads
    #[arg(short, long)]
    threads: Option<usize>,

    /// Output SAM file path
    #[arg(short, long)]
    sam: Option<PathBuf>,
}

// --- FastCon ---

#[derive(clap::Args)]
struct FastConArgs {
    /// Index files to combine results from
    #[arg(short, long)]
    indexes: Vec<PathBuf>,

    /// Reads file (FASTQ)
    #[arg(short, long)]
    reads: Option<PathBuf>,

    /// R1 FASTQ file for paired-end mapping
    #[arg(short = '1', long)]
    reads_1: Option<PathBuf>,

    /// R2 FASTQ file for paired-end mapping
    #[arg(short = '2', long)]
    reads_2: Option<PathBuf>,

    /// Number of parallel threads
    #[arg(short, long)]
    threads: Option<usize>,

    /// Output SAM file path
    #[arg(short, long)]
    sam: Option<PathBuf>,

    /// Path to rna-pop executable
    #[arg(long)]
    rna_pop_exe: Option<PathBuf>,
}

// --- ChunkConsensus ---

#[derive(clap::Args)]
struct ChunkConsensusArgs {
    /// Index file
    #[arg(short, long)]
    index: PathBuf,

    /// Reads file (FASTQ)
    #[arg(short, long)]
    reads: Option<PathBuf>,

    /// R1 FASTQ file for paired-end mapping
    #[arg(short = '1', long)]
    reads_1: Option<PathBuf>,

    /// R2 FASTQ file for paired-end mapping
    #[arg(short = '2', long)]
    reads_2: Option<PathBuf>,

    /// Chunk percentages to try (comma-separated, e.g. "10,20,30")
    #[arg(long, default_value = "10,20,30")]
    chunks: String,

    /// Number of parallel threads
    #[arg(short, long)]
    threads: Option<usize>,

    /// Output SAM file path
    #[arg(short, long)]
    sam: Option<PathBuf>,
}

// --- QC ---

#[derive(clap::Args)]
struct QcArgs {
    /// SAM file to analyze
    #[arg(short, long)]
    sam: PathBuf,

    /// Output QC report file (optional)
    #[arg(short, long)]
    output: Option<PathBuf>,
}

// --- Report ---

#[derive(clap::Args)]
struct ReportArgs {
    /// SAM file from mapping
    #[arg(short, long)]
    sam: PathBuf,

    /// TSV file with abundances from EM
    #[arg(short, long)]
    abundances: PathBuf,

    /// Output HTML report path
    #[arg(short, long, default_value = "report.html")]
    output: PathBuf,

    /// Cancer panels to evaluate (comma-separated: breast,lung,prostate,colorectal,pancancer)
    #[arg(long, default_value = "pancancer")]
    panels: String,
}

// --- Compare ---

#[derive(clap::Args)]
struct CompareArgs {
    /// First sample abundance TSV
    #[arg(short = '1', long)]
    abundances_1: PathBuf,

    /// Second sample abundance TSV
    #[arg(short = '2', long)]
    abundances_2: PathBuf,

    /// Output TSV file with comparison results
    #[arg(short, long)]
    output: Option<PathBuf>,

    /// Output HTML report
    #[arg(long)]
    report: Option<PathBuf>,

    /// Fold change threshold for significance
    #[arg(long, default_value_t = 1.5)]
    fc_threshold: f64,

    /// p-value threshold for significance
    #[arg(long, default_value_t = 0.05)]
    p_threshold: f64,
}

fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::Run(args) => cmd_run(args),
        Commands::Build(args) => cmd_build(args),
        Commands::Map(args) => cmd_map(args),
        Commands::Stats(args) => cmd_stats(args),
        Commands::Em(args) => cmd_em(args),
        Commands::Consensus(args) => cmd_consensus(args),
        Commands::FastCon(args) => cmd_fastcon(args),
        Commands::ChunkConsensus(args) => cmd_chunk_consensus(args),
        Commands::Qc(args) => cmd_qc(args),
        Commands::Report(args) => cmd_report(args),
        Commands::Compare(args) => cmd_compare(args),
    }
}

fn get_reads_path(args_reads: &Option<PathBuf>, args_reads_1: &Option<PathBuf>, args_reads_2: &Option<PathBuf>) -> PathBuf {
    match (args_reads, args_reads_1, args_reads_2) {
        (Some(r), _, _) => r.clone(),
        (_, Some(r1), _) => r1.clone(),
        _ => {
            eprintln!("Error: --reads or --reads_1 required");
            std::process::exit(1);
        }
    }
}

fn get_output_path(args_sam: &Option<PathBuf>) -> String {
    args_sam.as_ref().map(|p| p.display().to_string())
        .unwrap_or_else(|| "output.sam".to_string())
}

fn cmd_run(args: RunArgs) {
    let start = Instant::now();

    let align_mode = match args.align.to_lowercase().as_str() {
        "sw" => AlignMode::Sw,
        "hybrid" => AlignMode::Hybrid,
        "softclip" => AlignMode::Softclip,
        "chain" => AlignMode::Chain,
        _ => AlignMode::Xor,
    };

    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap();
    }

    let mut bp = if let Some(ref index_path) = args.index {
        println!("Loading index from {}", index_path.display());
        BitPop::deserialize_from_file(index_path.to_str().unwrap()).unwrap()
    } else if let Some(ref genome) = args.genome {
        println!("Building index from {}", genome);
        let mut bp = BitPop::new(args.k.unwrap_or(10));
        bp.set_auto_k(args.k.is_none());
        bp.set_top_n(args.top_n);
        bp.set_read_type(&args.read_type);
        bp.set_align_mode(align_mode);

        let genome_path = PathBuf::from(genome);
        if genome_path.is_dir() {
            for f in std::fs::read_dir(genome_path).unwrap() {
                let p = f.unwrap().path();
                if p.extension().map(|e| e.to_string_lossy().to_lowercase().contains("fasta")).unwrap_or(false)
                    || p.extension().map(|e| e.to_string_lossy().contains("fa")).unwrap_or(false)
                {
                    bp.load_genome_fasta(p.to_str().unwrap()).unwrap();
                }
            }
        } else {
            bp.load_genome_fasta(genome).unwrap();
        }

        if args.threads.unwrap_or(1) > 1 {
            bp.build_parallel();
        } else {
            bp.build();
        }
        bp
    } else {
        eprintln!("Error: --genome or --index required");
        std::process::exit(1);
    };

    bp.set_top_n(args.top_n);
    bp.set_read_type(&args.read_type);
    bp.set_align_mode(align_mode);

    let reads_path = get_reads_path(&args.reads, &args.reads_1, &args.reads_2);
    let output = get_output_path(&args.sam);

    let parsed = parse_fastq(reads_path.to_str().unwrap()).unwrap();
    let reads: Vec<(String, String)> = parsed.into_iter().map(|(n, s, _)| (n, s)).collect();
    println!("Loaded {} reads", reads.len());

    let start_map = Instant::now();
    let reads_ref: Vec<(&str, &str)> = reads.iter().map(|(n, s)| (n.as_str(), s.as_str())).collect();
    let mapped = bp.map_reads_parallel(&reads_ref, &output, 10).unwrap();
    let map_time = start_map.elapsed();

    println!("Mapped {} reads in {:.2}s ({:.1} reads/s)",
        mapped, map_time.as_secs_f64(), mapped as f64 / map_time.as_secs_f64());

    let total_time = start.elapsed();
    println!("Total time: {:.2}s", total_time.as_secs_f64());
}

fn cmd_build(args: BuildArgs) {
    let start = Instant::now();

    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap();
    }

    let mut bp = BitPop::new(args.k.unwrap_or(10));
    bp.set_auto_k(args.k.is_none());

    let genome = args.genome.unwrap_or_else(|| {
        eprintln!("Error: --genome required");
        std::process::exit(1);
    });

    let genome_path = PathBuf::from(&genome);
    if genome_path.is_dir() {
        for f in std::fs::read_dir(&genome_path).unwrap() {
            let p = f.unwrap().path();
            if p.extension().map(|e| e.to_string_lossy().to_lowercase().contains("fasta")).unwrap_or(false)
                || p.extension().map(|e| e.to_string_lossy().contains("fa")).unwrap_or(false)
            {
                bp.load_genome_fasta(p.to_str().unwrap()).unwrap();
            }
        }
    } else {
        bp.load_genome_fasta(&genome).unwrap();
    }

    println!("Building index...");
    if args.threads.unwrap_or(1) > 1 {
        bp.build_parallel();
    } else {
        bp.build();
    }

    let output = args.output.unwrap_or_else(|| {
        let stem = genome_path.file_stem()
            .unwrap_or_default()
            .to_string_lossy()
            .to_string();
        PathBuf::from(format!("{}.rnapop", stem))
    });

    bp.serialize_to_file(output.to_str().unwrap()).unwrap();
    println!("Index saved to {}", output.display());

    let elapsed = start.elapsed();
    println!("Build time: {:.2}s", elapsed.as_secs_f64());
}

fn cmd_map(args: MapArgs) {
    let start = Instant::now();

    let align_mode = match args.align.to_lowercase().as_str() {
        "sw" => AlignMode::Sw,
        "hybrid" => AlignMode::Hybrid,
        "softclip" => AlignMode::Softclip,
        "chain" => AlignMode::Chain,
        _ => AlignMode::Xor,
    };

    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .unwrap();
    }

    println!("Loading index from {}", args.index.display());
    let mut bp = BitPop::deserialize_from_file(args.index.to_str().unwrap()).unwrap();
    bp.set_top_n(args.top_n);
    bp.set_read_type(&args.read_type);
    bp.set_align_mode(align_mode);

    let reads_path = get_reads_path(&args.reads, &args.reads_1, &args.reads_2);
    let output = get_output_path(&args.sam);

    let parsed = parse_fastq(reads_path.to_str().unwrap()).unwrap();
    let reads: Vec<(String, String)> = parsed.into_iter().map(|(n, s, _)| (n, s)).collect();
    println!("Loaded {} reads", reads.len());

    let start_map = Instant::now();
    let reads_ref: Vec<(&str, &str)> = reads.iter().map(|(n, s)| (n.as_str(), s.as_str())).collect();
    let mapped = bp.map_reads_parallel(&reads_ref, &output, 10).unwrap();
    let map_time = start_map.elapsed();

    println!("Mapped {} reads in {:.2}s ({:.1} reads/s)",
        mapped, map_time.as_secs_f64(), mapped as f64 / map_time.as_secs_f64());

    let total_time = start.elapsed();
    println!("Total time: {:.2}s", total_time.as_secs_f64());
}

fn cmd_stats(args: StatsArgs) {
    println!("Loading index from {}", args.index.display());
    let bp = BitPop::deserialize_from_file(args.index.to_str().unwrap()).unwrap();

    println!("=== Index Statistics ===");
    println!("K-mer size: {}", bp.k());
    println!("Number of genomes: {}", bp.genome_count());

    for gid in 0..bp.genome_count() {
        let gid = gid as u32;
        if let Some(name) = bp.genome_name(gid) {
            let len = bp.genome_seq_len(gid).unwrap_or(0);
            println!("  {} (ID: {}): {} bp", name, gid, len);
        }
    }
}

fn cmd_em(args: EmArgs) {
    use rna_pop::em::{EMClassifier, EMConfig, ReadMappings};

    println!("Loading SAM file from {}", args.sam.display());

    let mut mappings: ReadMappings = Vec::new();
    let mut line = String::new();
    let file = std::fs::File::open(&args.sam).unwrap();
    let mut reader = std::io::BufReader::new(file);

    while reader.read_line(&mut line).unwrap() > 0 {
        if line.starts_with('@') {
            line.clear();
            continue;
        }
        let fields: Vec<&str> = line.trim().split('\t').collect();
        if fields.len() >= 3 {
            let read_name = fields[0].to_string();
            let genome_name = fields[2].to_string();
            if genome_name == "*" {
                line.clear();
                continue;
            }
            let score = fields.get(5).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.0);
            mappings.push((read_name, genome_name, score));
        }
        line.clear();
    }

    println!("Loaded {} mappings", mappings.len());

    // Build genome lengths from index if available
    let mut genome_lengths: std::collections::HashMap<String, usize> = std::collections::HashMap::new();
    if let Some(ref index_path) = args.index {
        println!("Loading index for transcript lengths...");
        let bp = BitPop::deserialize_from_file(index_path.to_str().unwrap()).unwrap();
        for gid in 0..bp.genome_count() {
            let gid = gid as u32;
            if let Some(name) = bp.genome_name(gid) {
                let len = bp.genome_seq_len(gid).unwrap_or(1000);
                genome_lengths.insert(name.to_string(), len);
            }
        }
        println!("Loaded lengths for {} transcripts", genome_lengths.len());
    }

    let config = EMConfig {
        convergence_threshold: args.threshold,
        max_iterations: args.iterations,
        genome_lengths,
        ..EMConfig::default()
    };

    let mut em = EMClassifier::new(config);
    let _results = em.classify(&mappings);

    let report = em.get_abundance_report();
    println!("EM converged in {} iterations (KL={:.6})", em.iterations_run, em.final_kl);
    println!("\n{:40} {}", "Genome", "Abundance");
    println!("{}", "-".repeat(60));
    for (genome, abundance) in &report {
        println!("{:40} {:.6}", genome, abundance);
    }

    if let Some(ref tsv_path) = args.tsv {
        use std::io::Write;
        let mut tsv = std::fs::File::create(tsv_path).unwrap();
        writeln!(tsv, "transcript\tabundance").unwrap();
        for (genome, abundance) in &report {
            writeln!(tsv, "{}\t{:.6}", genome, abundance).unwrap();
        }
        println!("\nTSV saved to {}", tsv_path.display());
    }
}

fn cmd_consensus(args: ConsensusArgs) {
    let start = Instant::now();

    let reads_path = get_reads_path(&args.reads, &args.reads_1, &args.reads_2);
    let output = get_output_path(&args.sam);

    let consensus = MultiKConsensus::from_paths(
        &[args.index.clone()],
        0.0,
    ).unwrap();
    let start_map = Instant::now();
    let (mapped, _) = consensus.map_reads_to_sam(&reads_path, std::path::Path::new(&output), args.threads.unwrap_or(1)).unwrap();
    let map_time = start_map.elapsed();

    println!("Mapped {} reads in {:.2}s ({:.1} reads/s)",
        mapped, map_time.as_secs_f64(), mapped as f64 / map_time.as_secs_f64());

    let total_time = start.elapsed();
    println!("Total time: {:.2}s", total_time.as_secs_f64());
}

fn cmd_fastcon(args: FastConArgs) {
    let start = Instant::now();

    if args.indexes.is_empty() {
        eprintln!("Error: at least one --index required");
        std::process::exit(1);
    }

    let reads_path = get_reads_path(&args.reads, &args.reads_1, &args.reads_2);
    let output = get_output_path(&args.sam);

    let rna_pop_exe = args.rna_pop_exe.unwrap_or_else(|| {
        std::env::current_exe().unwrap().parent().unwrap().to_path_buf()
    });

    let fastcon = FastCon::new(
        args.indexes,
        rna_pop_exe,
        0.0,
        10,
        0.5,
        false,
    ).unwrap();

    let start_map = Instant::now();
    let (mapped, _) = fastcon.run(&reads_path, &PathBuf::from(&output), args.threads.unwrap_or(1)).unwrap();
    let map_time = start_map.elapsed();

    println!("Mapped {} reads in {:.2}s ({:.1} reads/s)",
        mapped, map_time.as_secs_f64(), mapped as f64 / map_time.as_secs_f64());

    let total_time = start.elapsed();
    println!("Total time: {:.2}s", total_time.as_secs_f64());
}

fn cmd_chunk_consensus(args: ChunkConsensusArgs) {
    let start = Instant::now();

    let chunk_pcts: Vec<f64> = args.chunks.split(',').map(|s| s.trim().parse().unwrap()).collect();

    let reads_path = get_reads_path(&args.reads, &args.reads_1, &args.reads_2);
    let output = get_output_path(&args.sam);

    let chunk_consensus = MultiChunkConsensus::from_path(
        args.index.to_str().unwrap(),
        &chunk_pcts,
        50,
        200,
        0.0,
        0.5,
    ).unwrap();

    let start_map = Instant::now();
    let (mapped, _) = chunk_consensus.map_reads_to_sam(&reads_path, std::path::Path::new(&output), args.threads.unwrap_or(1)).unwrap();
    let map_time = start_map.elapsed();

    println!("Mapped {} reads in {:.2}s ({:.1} reads/s)",
        mapped, map_time.as_secs_f64(), mapped as f64 / map_time.as_secs_f64());

    let total_time = start.elapsed();
    println!("Total time: {:.2}s", total_time.as_secs_f64());
}

fn cmd_qc(args: QcArgs) {
    match run_qc(&args.sam, args.output.as_ref()) {
        Ok(_) => {},
        Err(e) => {
            eprintln!("Error running QC: {}", e);
            std::process::exit(1);
        }
    }
}

fn cmd_report(args: ReportArgs) {
    // Run QC first
    let qc = match run_qc(&args.sam, None) {
        Ok(metrics) => metrics,
        Err(e) => {
            eprintln!("Error running QC: {}", e);
            std::process::exit(1);
        }
    };

    // Read abundances
    let abundances_file = match std::fs::File::open(&args.abundances) {
        Ok(f) => f,
        Err(e) => {
            eprintln!("Error opening abundances file: {}", e);
            std::process::exit(1);
        }
    };

    let reader = std::io::BufReader::new(abundances_file);
    let mut abundances = HashMap::new();

    for (i, line) in reader.lines().enumerate() {
        if i == 0 {
            continue; // Skip header
        }
        let line = match line {
            Ok(l) => l,
            Err(_) => continue,
        };
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() >= 2 {
            let transcript = parts[0].trim().to_string();
            let abundance: f64 = parts[1].trim().parse().unwrap_or(0.0);
            abundances.insert(transcript, abundance);
        }
    }

    // Parse panels
    let panels: Vec<&str> = args.panels.split(',').map(|s| s.trim()).collect();

    // Generate report
    match generate_clinical_report(&args.output, &qc, &abundances, &panels) {
        Ok(_) => println!("Report generated: {}", args.output.display()),
        Err(e) => {
            eprintln!("Error generating report: {}", e);
            std::process::exit(1);
        }
    }
}

fn cmd_compare(args: CompareArgs) {
    // Run comparison
    let results = match compare_samples(&args.abundances_1, &args.abundances_2, args.output.as_ref(), args.fc_threshold, args.p_threshold) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Error comparing samples: {}", e);
            std::process::exit(1);
        }
    };

    // Generate HTML report if requested
    if let Some(report_path) = args.report {
        match generate_comparison_report(&report_path, &results, args.fc_threshold, args.p_threshold) {
            Ok(_) => println!("Comparison report generated: {}", report_path.display()),
            Err(e) => {
                eprintln!("Error generating report: {}", e);
                std::process::exit(1);
            }
        }
    }
}
