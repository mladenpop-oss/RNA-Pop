# RNA-Pop

**RNA-seq read quantification and differential expression analysis** using FM-index mapping and EM-based abundance estimation.

RNA-Pop maps RNA-seq reads to transcriptome references using a compact FM-index with bit-level parallelism, quantifies expression at transcript level using an EM algorithm with length normalization, and provides quality control, biomarker analysis, and differential expression tools for cancer research.

## Features

- **FM-index mapping** — SA-IS constructed FM-index with 2-bit DNA encoding for compact memory usage
- **Splice-aware alignment** — XOR, Smith-Waterman, hybrid, softclip, and chain alignment modes
- **EM quantification** — Soft-assignment EM algorithm with transcript length normalization (Spearman ρ = 0.99)
- **Parallel mapping** — Multi-threaded read mapping with Rayon (21.8K reads/s on 16 threads)
- **Consensus modes** — Multi-k, multi-chunk, and fast consensus for improved accuracy
- **Quality control** — Mapping statistics, coverage analysis, and uniformity metrics
- **Cancer biomarkers** — Pan-cancer panels for breast, lung, prostate, colorectal cancers
- **Differential expression** — Multi-sample comparison with fold change and significance testing
- **Clinical reports** — HTML reports with biomarker analysis and QC metrics
- **SAM output** — Standard SAM format for downstream analysis
- **FASTQ/FASTA** — Supports both single-end and paired-end reads

## Installation

```bash
cargo install --path .
```

Or build from source:

```bash
cargo build --release
./target/release/rna-pop
```

## Quick Start

### Single command: build + map

```bash
rna-pop run --genome transcripts.fa --reads reads.fastq --sam output.sam --threads 8
```

### Multi-step workflow

```bash
# 1. Build index
rna-pop build --genome transcripts.fa --output transcripts.rnp --threads 8

# 2. Map reads
rna-pop map --index transcripts.rnp --reads reads.fastq --sam output.sam --threads 8

# 3. Quantify with EM
rna-pop em --sam output.sam --index transcripts.rnp --tsv abundances.tsv
```

## Commands

| Command | Description |
|---------|-------------|
| `run` | Build index + map reads in one step |
| `build` | Build FM-index from transcript FASTA |
| `map` | Map reads to indexed transcripts |
| `stats` | Show index statistics |
| `em` | EM-based transcript quantification |
| `qc` | Quality control metrics for SAM file |
| `report` | Generate clinical report with biomarker analysis |
| `compare` | Differential expression between two samples |
| `consensus` | Multi-k consensus mapping |
| `fastcon` | Fast consensus across multiple indexes |
| `chunkconsensus` | Multi-chunk consensus mapping |

### `build`

```bash
rna-pop build --genome transcripts.fa --output transcripts.rnp --k 22 --threads 8
```

| Option | Description | Default |
|--------|-------------|---------|
| `--genome` | Transcript FASTA file or directory | required |
| `--output` | Output index path | `<genome>.rnp` |
| `--k` | K-mer size | auto |
| `--threads` | Parallel threads | 1 |

### `map`

```bash
rna-pop map --index transcripts.rnp --reads reads.fastq --sam output.sam --threads 8 --align hybrid
```

| Option | Description | Default |
|--------|-------------|---------|
| `--index` | Index file | required |
| `--reads` | FASTQ file (single-end) | - |
| `--reads_1`, `--reads_2` | Paired-end FASTQ | - |
| `--sam` | Output SAM file | `output.sam` |
| `--align` | Alignment mode: `xor`, `sw`, `hybrid`, `softclip`, `chain` | `xor` |
| `--threads` | Parallel threads | 1 |
| `--top_n` | Top rarest k-mers as anchors | 1 |

### `em`

```bash
rna-pop em --sam output.sam --index transcripts.rnp --tsv abundances.tsv --iterations 50
```

| Option | Description | Default |
|--------|-------------|---------|
| `--sam` | SAM file from mapping | required |
| `--index` | Index file (for transcript lengths) | - |
| `--tsv` | Output TSV with abundances | - |
  | `--iterations` | Max EM iterations | 20 |
| `--threshold` | Convergence threshold | 0.001 |

### `qc`

```bash
rna-pop qc --sam output.sam --output qc_report.txt
```

| Option | Description | Default |
|--------|-------------|---------|
| `--sam` | SAM file from mapping | required |
| `--output` | Output QC report file | console |

**Output:** Mapping rate, unique/multi-mapped reads, average mapping quality, transcript coverage, coverage uniformity.

### `report`

```bash
rna-pop report --sam output.sam --abundances abundances.tsv --output report.html --panels pancancer,breast,lung
```

| Option | Description | Default |
|--------|-------------|---------|
| `--sam` | SAM file from mapping | required |
| `--abundances` | TSV file with abundances from EM | required |
| `--output` | Output HTML report path | `report.html` |
| `--panels` | Cancer panels: `breast`,`lung`,`prostate`,`colorectal`,`pancancer` | `pancancer` |

**Output:** HTML report with QC metrics and cancer biomarker expression levels.

### `compare`

```bash
rna-pop compare -1 tumour.tsv -2 normal.tsv --output diff.tsv --report comparison.html
```

| Option | Description | Default |
|--------|-------------|---------|
| `-1, --abundances-1` | First sample abundance TSV | required |
| `-2, --abundances-2` | Second sample abundance TSV | required |
| `--output` | Output TSV with comparison results | - |
| `--report` | Output HTML report | - |
| `--fc-threshold` | Fold change threshold | 1.5 |
| `--p-threshold` | p-value threshold | 0.05 |

**Output:** Differentially expressed transcripts with fold change, log2FC, p-values, and significance.

## Cancer Research Pipeline

Complete workflow for cancer biomarker analysis:

```bash
# 1. Map tumour and normal samples
rna-pop run --genome transcripts.fa --reads tumour.fastq --sam tumour.sam --threads 16
rna-pop run --genome transcripts.fa --reads normal.fastq --sam normal.sam --threads 16

# 2. Quantify transcript abundances
rna-pop em --sam tumour.sam --index transcripts.rnp --tsv tumour.tsv
rna-pop em --sam normal.sam --index transcripts.rnp --tsv normal.tsv

# 3. Quality control
rna-pop qc --sam tumour.sam --output tumour_qc.txt

# 4. Differential expression analysis
rna-pop compare -1 tumour.tsv -2 normal.tsv --report differential.html

# 5. Clinical biomarker report
rna-pop report --sam tumour.sam --abundances tumour.tsv --output clinical.html --panels breast,lung,pancancer
```

## Accuracy

Benchmarked on simulated RNA-seq data (polyester, 11,567 human chr19 transcripts + ERCC spike-ins):

| Metric | 84K reads | 2M reads |
|--------|-----------|----------|
| Map rate | 97.0% | 97.1% |
| Speed | 16.6K reads/s | 21.8K reads/s (16 threads) |
| Spearman ρ (ranking) | 0.85 | **0.99** |
| Pearson r (linear) | 0.88 | **0.98** |
| Top-10 overlap | 9/10 | 7/10 |

### Comparison with Salmon

| Tool | Spearman ρ | Pearson r | Top-10 overlap |
|------|-----------|-----------|----------------|
| **RNA-Pop** | **0.9933** | **0.9787** | **7/10** |
| Salmon (NumReads) | 0.9889 | 0.9731 | 7/10 |

RNA-Pop achieves comparable accuracy to Salmon with a simpler architecture and no external dependencies.

## Output Format

### SAM output

Standard SAM format with mapping scores:

```
read1	0	ENST00000397910.8	100	60	100M	*	0	0	AGCTAGCT...	IIIIIIII...	NM:i:2
```

### TSV output (EM)

Tab-separated transcript abundances:

```
transcript	abundance
ENST00000397910.8	0.002104
ENST00000270460.10	0.001238
...
```

## Architecture

```
┌─────────────┐     ┌──────────────┐     ┌─────────────┐
│  Transcript  │ ──► │  FM-index    │ ──► │   Mapping   │
│   FASTA      │     │  (SA-IS)     │     │  (parallel) │
└─────────────┘     └──────────────┘     └──────┬──────┘
                                                 │
┌─────────────┐     ┌──────────────┐             │
│    EM        │ ◄── │ Length       │ ◄──────────┘
│Quantification│     │Normalization │
└─────────────┘     └──────────────┘
       │
       ▼
┌─────────────┐
│  Abundance  │
│   TSV       │
└─────────────┘
```

## Dependencies

- [Rayon](https://crates.io/crates/rayon) — parallel mapping
- [SA-IS](https://crates.io/crates/libsais) — suffix array construction
- [Clap](https://crates.io/crates/clap) — CLI parsing
- [Zstd](https://crates.io/crates/zstd) — compression support
- [Serde](https://crates.io/crates/serde) — serialization

## License

MIT

## Citation

DOI: https://zenodo.org/doi/10.5281/zenodo.20578611

## Author

Mladen Popović <mladenpop@gmail.com>
