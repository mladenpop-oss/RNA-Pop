# Changelog

All notable changes to this project will be documented in this file.

## [0.1.0] - 2026-06-07

### Added
- FM-index mapping with SA-IS construction
- EM-based transcript quantification with length normalization
- Multi-threaded parallel mapping with Rayon
- SAM output support
- CLI commands: run, build, map, stats, em
- Support for single-end and paired-end reads
- Multiple alignment modes: xor, sw, hybrid, softclip, chain

### Performance
- 21.8K reads/s on 16 threads (2M reads benchmark)
- 97.1% mapping rate on simulated RNA-seq data
- Spearman ρ = 0.99, Pearson r = 0.98 vs ground truth
