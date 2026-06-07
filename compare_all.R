# Read truth
truth <- read.table("/mnt/c/Users/Daddy/Documents/GitHub/test_data/sim_reads/truth_counts.tsv", header=TRUE, stringsAsFactors=FALSE, sep="\t")
truth$enst <- sub("\\..*", "", truth$transcript)

# Read RNA-Pop EM
rnapop_lines <- readLines("/mnt/c/Users/Daddy/Documents/GitHub/test_data/sim_reads/rnapop_1M_em.tsv")
rnapop_lines <- rnapop_lines[-1]
enst_vec <- c()
ab_vec <- c()
for (line in rnapop_lines) {
    parts <- strsplit(line, "\t")[[1]]
    if (length(parts) == 2) {
        tx <- parts[1]
        ab <- as.numeric(parts[2])
        enst <- sub("\\..*", "", sub(" .*", "", tx))
        enst_vec <- c(enst_vec, enst)
        ab_vec <- c(ab_vec, ab)
    }
}
rnapop <- data.frame(enst=enst_vec, abundance=ab_vec, stringsAsFactors=FALSE)

# Read Salmon quant — use NumReads, not TPM!
salmon <- read.table("/mnt/c/Users/Daddy/Documents/GitHub/test_data/sim_reads/salmon_out/quant.sf", header=TRUE, stringsAsFactors=FALSE)
salmon$enst <- sub("\\..*", "", salmon$Name)
salmon <- salmon[, c("enst", "NumReads")]
names(salmon)[2] <- "counts"

# Merge all
merged <- merge(truth, rnapop, by="enst", all.x=TRUE)
merged$abundance[is.na(merged$abundance)] <- 0
merged <- merge(merged, salmon, by="enst", all.x=TRUE)
merged$counts[is.na(merged$counts)] <- 0
merged <- merged[merged$true_counts > 0, ]

# Correlations
cat("=== RNA-Pop vs Truth ===\n")
cor_rnapop <- cor.test(merged$true_counts, merged$abundance, method="spearman")
cat("Spearman:", round(cor_rnapop$estimate, 4), "\n")
cor_rnapop_p <- cor.test(merged$true_counts, merged$abundance, method="pearson")
cat("Pearson:", round(cor_rnapop_p$estimate, 4), "\n\n")

cat("=== Salmon (NumReads) vs Truth ===\n")
cor_salmon <- cor.test(merged$true_counts, merged$counts, method="spearman")
cat("Spearman:", round(cor_salmon$estimate, 4), "\n")
cor_salmon_p <- cor.test(merged$true_counts, merged$counts, method="pearson")
cat("Pearson:", round(cor_salmon_p$estimate, 4), "\n\n")

# Top 10 comparison
top_truth <- truth[order(-truth$true_counts), ][1:10, ]
top_rnapop <- rnapop[order(-rnapop$abundance), ][1:10, ]
top_salmon <- salmon[order(-salmon$counts), ][1:10, ]

cat("=== Top 10 Truth ===\n")
print(top_truth[, c("enst", "true_counts")])

cat("\n=== Top 10 RNA-Pop ===\n")
print(top_rnapop[, c("enst", "abundance")])

cat("\n=== Top 10 Salmon (NumReads) ===\n")
print(top_salmon[, c("enst", "counts")])

# Overlap
overlap_rnapop <- intersect(top_truth$enst, top_rnapop$enst)
overlap_salmon <- intersect(top_truth$enst, top_salmon$enst)
cat("\nTop 10 overlap (Truth vs RNA-Pop):", length(overlap_rnapop), "/10\n")
cat("Top 10 overlap (Truth vs Salmon):", length(overlap_salmon), "/10\n")
