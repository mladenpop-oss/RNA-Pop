library(Biostrings)

# Read the simulated FASTA files
r1 <- readDNAStringSet("/mnt/c/Users/Daddy/Documents/GitHub/test_data/sim_reads/sample_01_1.fasta")
r2 <- readDNAStringSet("/mnt/c/Users/Daddy/Documents/GitHub/test_data/sim_reads/sample_01_2.fasta")

cat("R1 reads:", length(r1), "\n")
cat("R2 reads:", length(r2), "\n")

# Parse transcript IDs from headers
# Format: >read1/ENST00000339924.12 cdna chromosome:...;mate1:253-352;mate2:427-525
parse_tx_id <- function(name) {
    # Extract ENST ID
    parts <- strsplit(as.character(name), "/")[[1]]
    if (length(parts) >= 2) {
        # Get the ENST ID (first part after '/')
        tx <- strsplit(parts[2], " ")[[1]][1]
        return(tx)
    }
    return(NA)
}

tx_ids_r1 <- sapply(names(r1), parse_tx_id)
tx_ids_r2 <- sapply(names(r2), parse_tx_id)

# Count reads per transcript (ground truth)
truth_counts <- table(tx_ids_r1)
cat("\nGround truth: reads per transcript\n")
cat("Total reads:", sum(as.numeric(truth_counts)), "\n")
cat("Transcripts with reads:", length(truth_counts), "\n")
cat("Top 10 transcripts:\n")
print(head(sort(truth_counts, decreasing=TRUE), 10))

# Save ground truth
saveRDS(truth_counts, "/mnt/c/Users/Daddy/Documents/GitHub/test_data/sim_reads/truth_counts.rds")

# Write FASTQ files (simple conversion, quality = all 'I' = Q40)
quality <- DNAString(paste(rep("I", 100), collapse=""))

# Write R1 FASTQ
r1_fq <- paste0("@", seq_along(r1), "\n", 
                 as.character(r1), "\n+\n", 
                 as.character(quality))
writeLines(r1_fq, "/mnt/c/Users/Daddy/Documents/GitHub/test_data/sim_reads/sample_01_R1.fastq")

# Write R2 FASTQ
r2_fq <- paste0("@", seq_along(r2), "\n", 
                 as.character(r2), "\n+\n", 
                 as.character(quality))
writeLines(r2_fq, "/mnt/c/Users/Daddy/Documents/GitHub/test_data/sim_reads/sample_01_R2.fastq")

cat("\nFASTQ files written!\n")
