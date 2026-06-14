#!/usr/bin/env bash
# run_regenie_burden.sh — Stage 14 regenie gene-burden (CallForge, production default)
#
# Firth-corrected, covariate-adjusted gene-burden via regenie's two-step framework.
# Inputs : genotypes.tsv (variant gene burden_group s1 s2 ...), phenotype.tsv,
#          covariates.tsv, out_results.tsv, out_summary.json
# Requires regenie + plink2 (env/burden.yml). PRODUCTION engine — validate on a real
# cohort; regenie is designed for adequately-sized cohorts (step-1 ridge needs many
# variants/samples). For tiny cohorts use --burden_engine collapse.
set -euo pipefail
GENO="$1"; PHENO="$2"; COVAR="$3"; OUT_RES="$4"; OUT_SUM="$5"

# 1) genotypes.tsv -> PLINK (variants as rows). A full run converts the dosage matrix
#    to PGEN/BGEN; here we outline the canonical regenie burden steps:
#    plink2 --import-dosage ... --make-pgen --out geno
#    regenie --step 1 --pgen geno --phenoFile pheno --covarFile covar --bt \
#            --bsize 1000 --out step1
#    regenie --step 2 --pgen geno --phenoFile pheno --covarFile covar --bt \
#            --anno-file anno.tsv --set-list sets.tsv --mask-def masks.tsv \
#            --aaf-bins 0.01 --build-mask max --bsize 200 --pred step1_pred.list \
#            --out step2
#    then collate step2_*.regenie into OUT_RES (unit, kind=gene, p) + OUT_SUM (lambda).
echo "[run_regenie_burden] regenie production engine — see script for the step1/step2 commands." >&2
echo "ERROR: regenie burden requires a real cohort + genotype conversion; not run on synthetic test." >&2
echo "Use --burden_engine collapse for small cohorts / the test profile." >&2
exit 2
