# CallForge — decisions, caveats, and revisit-on-real-data items

Running log of design decisions and things to revisit when the real DBS cohort
arrives. Folded into the §11 completion report.

## Confirmed decisions
- **Reference**: the GRCh38 no-alt primary assembly CaptureForge designed against
  (`GRCh38.full.fa`, Ensembl-named) — enforced by the reference invariant.
- **CNV caller**: CNVkit default (`--cnv_caller`); GATK gCNV selectable. Germline
  two-step pooled panel-of-normals from QC-pass samples. Segmentation: `cbs`
  (R/DNAcopy) in production, `haar` (R-free) for the hermetic test.
- **SNV/indel caller**: GATK (HaplotypeCaller→GenomicsDBImport→GenotypeGVCFs).
  DeepVariant+GLnexus wired as the alternative (`--caller deepvariant`).
- **Burden engine**: regenie default; SKAT-O/STAAR (R) alternative.
- **Annotation gnomAD**: genome-wide gnomAD **genomes** v4.1 (whole-genome, covers
  coding + promoter + CNV + STR target classes with AFR AF). Exomes dropped as
  redundant for this panel (decision at the Phase 4/5 boundary).
- **Resource fetch**: `bin/fetch_annotation_dbs.sh` in **panel-slice mode** behind
  the scope gate; full genome-wide download stays the documented manual option. No
  pipeline flag triggers a hundreds-of-GB download silently.
- **Cohort QC: somalier + peddy → somalier-only.** peddy is dropped. Rationale:
  somalier already covers relatedness + sex + (with a 1000G-labelled reference)
  ancestry, is far more robust on sparse targeted-panel data, and needs no per-run
  pedigree; peddy adds a heavy 1000G-background dependency for little extra signal
  on a 59-gene panel. Revisit only if a pedigree-aware Mendelian-error report is
  required.
- **Scope gate hardened**: resource validation now also checks the BGZF EOF marker
  (rejects truncated/partial `.bgz` downloads), so the gate validates DATA integrity,
  not just contig overlap or index presence. "19/19 scope-fit" means every
  panel-gene chromosome's data is present + intact (the panel spans 19 of 24 chroms;
  chr10/18/21/X/Y have no panel genes, so are correctly absent from the slice).

## Quarantine / exclusion interlocks (proven on the test fixture)
- Low-coverage `LOWQC` sample is quarantined by the QC gate and **excluded from
  joint calling** (absent from the joint VCF). Proven Phase 3.
- CaptureForge labels **propagate to outputs**: CNV calls carry callability
  (S1 CFH DEL → `breakpoint_blind_low_confidence`; HP → `callable`), STR catalog
  genotypes HMOX1 (GT)20, paralog-region variants carry `PARALOG_GENE`+`PARALOG_CONF`
  in the VCF. Proven Phase 4.

## Revisit on real DBS data (NOT fixed now — flagged intentionally)
- **GATK hard-filter thresholds** are the GATK WGS/WES defaults; targeted-panel
  filtering may need tuning on the real cohort (Ti/Tv, pass rates, QUAL/DP).
- **Per-sample QC gate thresholds** (depth/dup/on-target/contamination) are defaults;
  calibrate against the real DBS coverage distribution.
- **On-target/enrichment/uniformity and Ti/Tv/het:hom** test numbers prove the
  machinery computes and gates fire — they are NOT data quality (synthetic mini-genome).
- **Sex check** is NA on this autosomal panel; somalier (Phase 6) infers sex
  genome-wide from off-target reads as the backstop.

## Deferred validations (do before production use)
- **DeepVariant + GLnexus** (Stage 6 alt): wired but untested (container-only). Run a
  one-time native-Apptainer smoke-test at HPC deployment.
- **regenie** (Stage 14, production default) and **SKAT-O / STAAR** (Stage 14 R alt):
  wired but NEVER executed (the synthetic cohort is too small). Exercise both on a
  real/realistic cohort before production use. The tested engine is `collapse`
  (CMC chi-square, small-N/default-free; engine name is stamped into burden_results.tsv).
- **Ancestry-PC stability is now a GATE (implemented), not a note**: build_burden_matrix
  drops the PCs and runs burden UNADJUSTED if there are too few PCA sites for the cohort
  (recorded in provenance.json). On real data, confirm PC1/PC2 are stable before relying
  on them as covariates; if not, the gate keeps the analysis unadjusted and says so.
- **BQSR known-sites**: currently via `--known_sites` param; manifest-driven
  auto-use deferred (Phase 5 small enhancement).
- **Genome-wide gnomAD panel-slice**: large/slow locally (gnomAD-genomes ~5.4 KB
  records); intended to be regenerated on HPC (fast network) for the real run.
- **Genome-wide dbSNP**: NOT yet fetched locally. Before the real cohort, discover
  or fetch a genome-wide dbSNP (scope gate applies as for gnomAD/ClinVar). Needed in
  TWO places: vcfanno rsID annotation (Stage 11) AND BQSR known-sites (Stage 4).
- **somalier sites**: supply the standard 1000G GRCh38 somalier sites VCF
  (`--somalier_sites`) for the real run — only the test fixture's sites exist now.
- **Ancestry PCs as burden covariates**: on a tight panel, ancestry/relatedness rely
  on (largely off-target) genome-wide sites; PCs may be UNSTABLE. FIRST thing to check
  on real data before using the PCs as regenie/SKAT covariates (Stage 14) — if
  unstable, use self-reported ancestry or external-reference projection instead.

## Apple-Silicon / local-runner notes
- The pinned `env/*.yml` are the production contract (linux/HPC native; Mac via
  `CONDA_SUBDIR=osx-64`/Rosetta). Local native-arm64 dev used per-tool isolated
  envs (calling stack, cnvkit, mosdepth, ExpansionHunter) — some bioconda tools
  (e.g. mosdepth 0.3.8) lack osx-arm64 builds and need ≥0.3.11 or Rosetta.
