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
- **DeepVariant + GLnexus**: wired but untested locally (container-only). Run a
  one-time native-Apptainer smoke-test at HPC deployment; flagged
  "validate before production use".
- **BQSR known-sites**: currently via `--known_sites` param; manifest-driven
  auto-use deferred (Phase 5 small enhancement).
- **Genome-wide gnomAD panel-slice**: large/slow locally (gnomAD-genomes ~5.4 KB
  records); intended to be regenerated on HPC (fast network) for the real run.
- **Genome-wide dbSNP**: NOT yet fetched locally. Before the real cohort, discover
  or fetch a genome-wide dbSNP (scope gate applies as for gnomAD/ClinVar). Needed in
  TWO places: vcfanno rsID annotation (Stage 11) AND BQSR known-sites (Stage 4).

## Apple-Silicon / local-runner notes
- The pinned `env/*.yml` are the production contract (linux/HPC native; Mac via
  `CONDA_SUBDIR=osx-64`/Rosetta). Local native-arm64 dev used per-tool isolated
  envs (calling stack, cnvkit, mosdepth, ExpansionHunter) — some bioconda tools
  (e.g. mosdepth 0.3.8) lack osx-arm64 builds and need ≥0.3.11 or Rosetta.
