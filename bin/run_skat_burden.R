#!/usr/bin/env Rscript
# run_skat_burden.R — Stage 14 SKAT-O gene-burden (CallForge, production/alt engine)
#
# Per-gene SKAT-O on variant-level dosages, adjusting for covariates + ancestry PCs.
# Inputs : genotypes.tsv (variant gene burden_group s1 s2 ...), phenotype.tsv,
#          covariates.tsv, out_results.tsv, out_summary.json
# Output : per-gene p-values + genomic-inflation lambda.
# Requires R + SKAT (env/burden.yml). Validate on a real cohort — SKAT needs an
# adequate sample size; small-N results are not valid association.
suppressMessages({ library(SKAT); library(jsonlite) })
a <- commandArgs(trailingOnly = TRUE)
geno <- read.delim(a[1], check.names = FALSE)
phe  <- read.delim(a[2]); cov <- read.delim(a[3])
res_out <- a[4]; sum_out <- a[5]

samples <- colnames(geno)[-(1:3)]
phe <- phe[match(samples, phe$sample_id), ]
cov <- cov[match(samples, cov$sample_id), ]
y <- as.numeric(phe$phenotype)
covmat <- as.matrix(cov[, setdiff(colnames(cov), "sample_id"), drop = FALSE])
storage.mode(covmat) <- "numeric"

# Null model (logistic: covariates + PCs)
obj <- if (ncol(covmat) > 0) SKAT_Null_Model(y ~ covmat, out_type = "D") else
                             SKAT_Null_Model(y ~ 1, out_type = "D")

genes <- unique(geno$gene)
rows <- lapply(genes, function(g) {
  Z <- t(as.matrix(geno[geno$gene == g, samples, drop = FALSE]))
  storage.mode(Z) <- "numeric"
  p <- tryCatch(SKAT(Z, obj, method = "SKATO")$p.value, error = function(e) NA)
  data.frame(unit = paste0("gene:", g), kind = "gene", engine = "SKAT-O",
             n_variants = sum(geno$gene == g), p = p)
})
df <- do.call(rbind, rows)
df <- df[order(df$p), ]
write.table(df, res_out, sep = "\t", quote = FALSE, row.names = FALSE)

chisq <- qchisq(df$p[!is.na(df$p)], df = 1, lower.tail = FALSE)
lambda <- if (length(chisq)) round(median(chisq) / qchisq(0.5, 1), 4) else NA
writeLines(toJSON(list(engine = "SKAT-O", n_units = nrow(df),
  genomic_inflation_lambda = lambda,
  caveat = "SKAT-O; validate on a real cohort — small-N is not valid association."),
  auto_unbox = TRUE, pretty = TRUE), sum_out)
cat(sprintf("[run_skat_burden] %d genes, lambda=%s\n", nrow(df), lambda))
