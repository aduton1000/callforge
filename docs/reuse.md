# CallForge — reused vs built

Per the engineering stance: stand on validated foundations (nf-core/sarek,
nf-core/raredisease) and add the bespoke layers. This table is updated as each
phase lands.

| Stage | Reused / adapted from | Built bespoke for CallForge |
|-------|-----------------------|------------------------------|
| 0 Inputs & discovery | nf-core sample-sheet validation idiom; nf-core `NO_FILE` optional-input idiom | `discover_resources.py` (organism-generic manifest), `check_reference_invariant.py` (no-alt + BED build assert, from CaptureForge), `ingest_captureforge.py` (per-gene metadata) |
| 1 Raw QC | nf-core FastQC module pattern | aggregation into the cohort dashboard |
| 2 Trimming | nf-core fastp module pattern | DBS-appropriate defaults |
| 3 Alignment | sarek bwa-mem2 + sort subworkflow | — |
| 4 Post-align | sarek MarkDuplicates + BQSR; Picard CollectHsMetrics | per-gene depth tracks (CaptureForge track style), graceful BQSR skip when no known-sites |
| 5 QC gate | — | flag/quarantine scorecard, thresholds, VerifyBamID2 sex/contamination wiring |
| 6 SNV/indel | sarek HaplotypeCaller→GenomicsDBImport→GenotypeGVCFs; DeepVariant+GLnexus alt | — |
| 7 Filtering | GATK hard-filter expressions | panel-tuned thresholds, documented |
| 8 CNV | CNVkit / GATK gCNV | callability labelling from CaptureForge `cnv_bin_spacing` |
| 9 STR | ExpansionHunter | catalog built from CaptureForge STR-class targets |
| 10 Paralog | — | confidence flagging in paralog-ambiguous regions (CaptureForge paralog flags) |
| 11 Annotation | raredisease VEP + vcfanno pattern | manifest-driven, organism-generic, gnomAD-AFR namespacing |
| 12 Cohort QC | somalier + peddy | ancestry PCA / relatedness plots in the dashboard |
| 13 GIAB | hap.py | restriction to the panel BED; PR/F1 plots |
| 14 Burden | regenie / SKAT-O / STAAR | gene-set collapse via CaptureForge burden groups; phenotype gating |
| 15 Reporting | MultiQC | self-contained cohort QC dashboard with captions, provenance.json |

> Stages 1–15 are scaffolded; the table marks intent and lands concretely as each
> phase is implemented at the extension point in `subworkflows/callforge.nf`.
