# RNA-Pop

**RNA-seq read quantification** using FM-index mapping and EM-based abundance estimation.

RNA-Pop maps RNA-seq reads to transcriptome references using a compact FM-index with bit-level parallelism, then quantifies expression at transcript level using an EM algorithm with length normalization.

## Features

- **FM-index mapping** — SA-IS constructed FM-index with 2-bit DNA encoding for compact memory usage
- **Splice-aware alignment** — XOR, Smith-Waterman, hybrid, softclip, and chain alignment modes
- **EM quantification** — Soft-assignment EM algorithm with transcript length normalization (Spearman ρ = 0.85)
- **Parallel mapping** — Multi-threaded read mapping with Rayon
- **Consensus modes** — Multi-k, multi-chunk, and fast consensus for improved accuracy
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

## Accuracy

Benchmarked on simulated RNA-seq data (polyester, 84,742 reads, 11,567 human chr19 transcripts):

| Metric | Value |
|--------|-------|
| Map rate | 97% |
| Spearman ρ (ranking) | 0.85 |
| Pearson r (linear) | 0.88 |
| Top-10 overlap | 9/10 |

Length normalization in the EM algorithm is critical for accurate quantification and is enabled automatically when `--index` is provided.

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

## Author

Mladen Popović <mladenpop@gmail.com>
