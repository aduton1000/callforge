#!/usr/bin/env python3
"""test_stage_entries.py — tests for CallForge stage subcommands (Gate 1: anno, burden).

Two layers:
  * Portable (always run): the thin wrapper scripts (validate_stage_inputs / stage_report
    / stage_provenance) on synthesized fixtures, and the CLI command construction + help.
  * Nextflow (auto-skip if no nextflow binary): -preview compiles the baseline pipeline
    and the anno/burden entries, and a real `--stage burden` run on results/ fixtures
    (pure-python + bcftools, no conda) asserting outputs + report + provenance appear and
    that a phenotype-less sheet fails loud.

Run:  python3 test/test_stage_entries.py
Set CALLFORGE_NEXTFLOW=/path/to/nextflow (or have `nextflow` on PATH) to enable layer 2.
"""
import csv, gzip, os, shutil, subprocess, sys, tempfile, unittest
from pathlib import Path

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
BIN = os.path.join(REPO, "bin")

def _nextflow():
    return os.environ.get("CALLFORGE_NEXTFLOW") or shutil.which("nextflow")

# Make a skipped Nextflow layer LOUD (never a silent pass): announce it up front.
if not _nextflow():
    sys.stderr.write(
        "\n[test_stage_entries] NOTE: no nextflow binary found — the Nextflow layer "
        "(baseline/anno preview, real burden run + fail-loud) will be reported as "
        "SKIPPED, not run.\n  Set CALLFORGE_NEXTFLOW=/path/to/nextflow (or put `nextflow` "
        "on PATH) to enable it.\n\n")

def run(cmd, **kw):
    return subprocess.run(cmd, capture_output=True, text=True, **kw)

def make_vcf_gz(path, contigs=("chr1",), with_index=True):
    """Write a minimal gzip'd VCF (gzip magic satisfies the bgzip check) + a .tbi."""
    body = ["##fileformat=VCFv4.2"]
    body += [f"##contig=<ID={c},length=1000>" for c in contigs]
    body += ["#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO",
             f"{contigs[0]}\t10\t.\tA\tT\t50\tPASS\t."]
    with gzip.open(path, "wt") as fh:
        fh.write("\n".join(body) + "\n")
    if with_index:
        open(path + ".tbi", "wb").close()


# ───────────────────────── layer 1: wrapper scripts ─────────────────────────
class ValidateInputs(unittest.TestCase):
    def val(self, *args):
        return run([sys.executable, os.path.join(BIN, "validate_stage_inputs.py"), *args])

    def test_index_missing_fails(self):
        with tempfile.TemporaryDirectory() as d:
            v = os.path.join(d, "x.vcf.gz"); make_vcf_gz(v, with_index=False)
            r = self.val("--stage", "anno", "--vcf", v)
            self.assertEqual(r.returncode, 2)
            self.assertIn("index missing", r.stderr)

    def test_contig_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as d:
            v = os.path.join(d, "x.vcf.gz"); make_vcf_gz(v, contigs=("chr1",))
            fai = os.path.join(d, "ref.fai")
            Path(fai).write_text("chrZ_only\t1000\t6\t60\t61\n")
            r = self.val("--stage", "anno", "--vcf", v, "--fai", fai)
            self.assertEqual(r.returncode, 2)
            self.assertIn("contig", r.stderr.lower())

    def test_contig_chr_reconciled_passes(self):
        # VCF 'chr1' vs reference '1' must reconcile (strip chr) and pass
        with tempfile.TemporaryDirectory() as d:
            v = os.path.join(d, "x.vcf.gz"); make_vcf_gz(v, contigs=("chr1",))
            fai = os.path.join(d, "ref.fai"); Path(fai).write_text("1\t1000\t6\t60\t61\n")
            r = self.val("--stage", "anno", "--vcf", v, "--fai", fai)
            self.assertEqual(r.returncode, 0, r.stderr)

    def test_phenotype_required_fails_when_blank(self):
        with tempfile.TemporaryDirectory() as d:
            s = os.path.join(d, "s.csv")
            Path(s).write_text("sample_id,fastq_1,fastq_2,phenotype\nA,a.fq,b.fq,\n")
            r = self.val("--stage", "burden", "--samplesheet", s, "--require-phenotype")
            self.assertEqual(r.returncode, 2)
            self.assertIn("phenotype", r.stderr.lower())

    def test_phenotype_present_passes(self):
        with tempfile.TemporaryDirectory() as d:
            s = os.path.join(d, "s.csv")
            Path(s).write_text("sample_id,fastq_1,fastq_2,phenotype\nA,a.fq,b.fq,case\n")
            r = self.val("--stage", "burden", "--samplesheet", s, "--require-phenotype")
            self.assertEqual(r.returncode, 0, r.stderr)


class ReportAndProvenance(unittest.TestCase):
    def test_stage_report_embeds_plots_and_metrics(self):
        with tempfile.TemporaryDirectory() as d:
            open(os.path.join(d, "burden_qq.png"), "wb").close()
            Path(os.path.join(d, "captions_burden.tsv")).write_text("burden_qq\tQQ plot of p-values.\n")
            Path(os.path.join(d, "burden_summary.json")).write_text('{"lambda": 1.0, "n_units": 4}')
            r = run([sys.executable, os.path.join(BIN, "stage_report.py"),
                     "--stage", "burden", "--title", "Burden", "--captions",
                     os.path.join(d, "captions_burden.tsv"), "--summary",
                     os.path.join(d, "burden_summary.json"), "--outdir", d,
                     "--out", "burden_report.md"])
            self.assertEqual(r.returncode, 0, r.stderr)
            md = Path(os.path.join(d, "burden_report.md")).read_text()
            self.assertIn("![QQ plot of p-values.](burden_qq.png)", md)
            self.assertIn("burden_summary.json", md)

    def test_stage_provenance_records_inputs_and_params(self):
        import json
        with tempfile.TemporaryDirectory() as d:
            f = os.path.join(d, "in.vcf.gz"); make_vcf_gz(f)
            r = run([sys.executable, os.path.join(BIN, "stage_provenance.py"),
                     "--stage", "burden", "--out", "p.json", "--outdir", d,
                     "--input", f"vcf={f}", "--param", "burden_engine=collapse",
                     "--version", "0.1.0", "--git-commit", "abc1234"])
            self.assertEqual(r.returncode, 0, r.stderr)
            prov = json.loads(Path(os.path.join(d, "p.json")).read_text())
            self.assertEqual(prov["stage"], "burden")
            self.assertEqual(prov["params"]["burden_engine"], "collapse")
            self.assertEqual(len(prov["inputs"]["vcf"]["sha256"]), 64)


# ─────────── layer 1: CaptureForge-label propagation on STANDALONE outputs ───────────
class LabelPropagation(unittest.TestCase):
    """The label-applying module commands, run exactly as the standalone cnv/str entries
    invoke them, on synthesized CaptureForge metadata — proving the labels survive when a
    stage is run alone (not just in-pipeline)."""

    def _meta_and_bed(self, d):
        meta = os.path.join(d, "gene_metadata.tsv")
        Path(meta).write_text(
            "gene\tcnv_callable\tis_str_target\tstr_loci\n"
            "CFH\tno\tno\t\n"          # breakpoint-blind (CaptureForge) -> low-confidence
            "HP\tyes\tno\t\n"          # depth-callable -> callable
            "CD209\tn/a\tyes\ttest_CD209:400-700;test_CD209:813-933\n")  # CaptureForge STR loci
        bed = os.path.join(d, "targets.bed")
        Path(bed).write_text("test_CFH\t400\t4532\tCFH|cnv\ntest_HP\t400\t2395\tHP|cnv\n"
                             "test_CD209\t400\t933\tCD209|str\n")
        return meta, bed

    def test_cnv_annotate_labels_breakpoint_blind(self):
        with tempfile.TemporaryDirectory() as d:
            meta, bed = self._meta_and_bed(d)
            cns = os.path.join(d, "SPIKE.call.cns")   # spike CFH + HP deletions
            Path(cns).write_text("chromosome\tstart\tend\tgene\tlog2\tcn\n"
                                 "test_CFH\t400\t4532\tCFH\t-0.7\t1\n"
                                 "test_HP\t400\t2395\tHP\t-0.6\t1\n")
            r = run([sys.executable, os.path.join(BIN, "cnv_annotate.py"), "--cns", cns,
                     "--bed", bed, "--metadata", meta, "--outdir", d])
            self.assertEqual(r.returncode, 0, r.stderr)
            calls = Path(os.path.join(d, "cnv_calls.tsv")).read_text()
            cfh = [l for l in calls.splitlines() if l.startswith("SPIKE\tCFH\t")][0]
            self.assertIn("breakpoint_blind_low_confidence", cfh)
            hp = [l for l in calls.splitlines() if l.startswith("SPIKE\tHP\t")][0]
            self.assertIn("callable", hp)

    def test_str_catalog_built_from_captureforge_loci(self):
        import json
        with tempfile.TemporaryDirectory() as d:
            meta, _ = self._meta_and_bed(d)
            out = os.path.join(d, "str_catalog.json")
            r = run([sys.executable, os.path.join(BIN, "build_str_catalog.py"),
                     "--metadata", meta, "--motif", "GT", "--out", out])
            self.assertEqual(r.returncode, 0, r.stderr)
            regions = [e.get("ReferenceRegion") for e in json.loads(Path(out).read_text())]
            self.assertTrue(any(str(x).startswith("test_CD209") for x in regions),
                            f"CaptureForge STR locus not in catalog: {regions}")


# ───────────────────────── layer 1: CLI ─────────────────────────
class CliStages(unittest.TestCase):
    def cli(self, *args):
        env = dict(os.environ, PYTHONPATH=REPO)
        return run([sys.executable, "-m", "callforge_cli.cli", *args], cwd=REPO, env=env)

    def test_top_help_lists_stage_subcommands(self):
        r = self.cli("--help")
        self.assertEqual(r.returncode, 0)
        self.assertIn("anno", r.stdout)
        self.assertIn("burden", r.stdout)

    def test_anno_help_lists_input_contract(self):
        r = self.cli("anno", "--help")
        self.assertEqual(r.returncode, 0)
        for tok in ("--input_vcf", "--genome_fasta", "--target_bed"):
            self.assertIn(tok, r.stdout)

    def test_anno_constructs_stage_command(self):
        r = self.cli("anno", "--print-cmd", "--input_vcf", "x.vcf.gz",
                     "--genome_fasta", "ref.fa", "--target_bed", "t.bed",
                     "-profile", "mac_local,conda")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("nextflow run", r.stdout)
        self.assertIn("--stage anno", r.stdout)
        self.assertIn("-profile mac_local,conda", r.stdout)

    def test_burden_passes_engine_through(self):
        r = self.cli("burden", "--print-cmd", "--input_vcf", "a.vcf.gz",
                     "--input", "s.csv", "--target_bed", "t.bed",
                     "--burden_engine", "regenie")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("--stage burden", r.stdout)
        self.assertIn("--burden_engine regenie", r.stdout)

    def test_default_injects_nondestructive_outdir(self):
        r = self.cli("anno", "--print-cmd", "--input_vcf", "x.vcf.gz",
                     "--genome_fasta", "ref.fa", "--target_bed", "t.bed")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("--outdir results/standalone/anno_", r.stdout)

    def test_in_place_opts_into_results(self):
        r = self.cli("burden", "--print-cmd", "--in-place", "--input_vcf", "a.vcf.gz",
                     "--input", "s.csv", "--target_bed", "t.bed")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("--in_place true", r.stdout)
        self.assertNotIn("results/standalone", r.stdout)

    def test_explicit_outdir_respected(self):
        r = self.cli("anno", "--print-cmd", "--input_vcf", "x.vcf.gz", "--genome_fasta",
                     "ref.fa", "--target_bed", "t.bed", "--outdir", "/tmp/myrun")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("--outdir /tmp/myrun", r.stdout)
        self.assertNotIn("results/standalone", r.stdout)

    def test_all_stage_subcommands_listed(self):
        r = self.cli("--help")
        for st in ("align", "coverage", "call", "cnv", "str", "paralog", "anno", "burden"):
            self.assertIn(st, r.stdout)

    def test_cnv_help_mentions_callability(self):
        r = self.cli("cnv", "--help")
        self.assertEqual(r.returncode, 0)
        self.assertIn("breakpoint_blind_low_confidence", r.stdout)

    def test_paralog_help_mentions_info_fields(self):
        r = self.cli("paralog", "--help")
        self.assertEqual(r.returncode, 0)
        self.assertIn("PARALOG_GENE", r.stdout)
        self.assertIn("--input_vcf", r.stdout)

    def test_str_help_mentions_catalog(self):
        r = self.cli("str", "--help")
        self.assertEqual(r.returncode, 0)
        self.assertIn("--str_catalog", r.stdout)

    def test_cohortqc_help_mentions_chaining_and_somalier(self):
        r = self.cli("cohortqc", "--help")
        self.assertEqual(r.returncode, 0)
        self.assertIn("somalier", r.stdout)
        self.assertIn("cohort_qc.json", r.stdout)   # the file burden consumes

    def test_giab_help_mentions_on_target_and_runinfo(self):
        r = self.cli("giab", "--help")
        self.assertEqual(r.returncode, 0)
        self.assertIn("on-target", r.stdout)
        self.assertIn("runinfo", r.stdout)

    def test_call_help_states_passbam_semantics(self):
        r = self.cli("call", "--help")
        self.assertEqual(r.returncode, 0)
        self.assertIn("QC-PASS", r.stdout)
        self.assertIn("does NOT re-run", r.stdout)
        self.assertIn("--input_bams", r.stdout)

    def test_align_help_lists_fastq_or_sheet(self):
        r = self.cli("align", "--help")
        self.assertEqual(r.returncode, 0)
        for tok in ("--input", "--fastq_1", "--genome_fasta", "--target_bed"):
            self.assertIn(tok, r.stdout)

    def test_coverage_constructs_stage_command(self):
        r = self.cli("coverage", "--print-cmd", "--input_bams", "out/*.bam",
                     "--genome_fasta", "ref.fa", "--target_bed", "t.bed")
        self.assertEqual(r.returncode, 0, r.stderr)
        self.assertIn("--stage coverage", r.stdout)
        self.assertIn("--input_bams", r.stdout)
        self.assertIn("out/*.bam", r.stdout)   # shlex-quoted because of the glob


# ───────────────────────── layer 2: Nextflow (auto-skip) ─────────────────────────
@unittest.skipUnless(_nextflow(), "no nextflow binary (set CALLFORGE_NEXTFLOW or PATH)")
class NextflowEntries(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.nf = _nextflow()
        cls.genome = os.path.join(REPO, "test/data/ref/test_genome.fa")
        cls.bed = os.path.join(REPO, "test/data/ref/test_targets.bed")
        cls.sheet = os.path.join(REPO, "test/data/test_samplesheet.csv")
        cls.anno_vcf = os.path.join(REPO, "results/stage10_paralog/paralog.annotated.vcf.gz")
        cls.burden_vcf = os.path.join(REPO, "results/stage11_annotation/annotated.vcf.gz")
        cls.filtered_vcf = os.path.join(REPO, "results/stage7_filter/joint.filtered.vcf.gz")
        cls.cohort_qc_json = os.path.join(REPO, "results/stage12_cohortqc/cohort_qc.json")

    def _preview(self, *extra):
        return run([self.nf, "run", "main.nf", "-profile", "test", "-preview", *extra], cwd=REPO)

    def test_baseline_still_previews_green(self):
        r = self._preview()
        self.assertIn("SUCCESS", r.stdout + r.stderr, r.stderr)

    def test_anno_entry_previews(self):
        if not os.path.exists(self.anno_vcf):
            self.skipTest("results/ fixture missing (run the test pipeline first)")
        # non-destructive guard requires an explicit --outdir for a standalone run
        r = self._preview("--stage", "anno", "--input_vcf", self.anno_vcf,
                          "--genome_fasta", self.genome, "--target_bed", self.bed,
                          "--outdir", "results/standalone/_preview")
        self.assertIn("SUCCESS", r.stdout + r.stderr, r.stderr)

    def test_standalone_guard_blocks_default_results(self):
        # without --outdir/--in_place a standalone run must refuse to overwrite results/
        r = self._preview("--stage", "burden", "--input_vcf", self.burden_vcf,
                          "--input", self.sheet, "--target_bed", self.bed)
        self.assertNotIn("SUCCESS", r.stdout + r.stderr)
        self.assertIn("must not overwrite", r.stdout + r.stderr)

    def test_align_entry_previews(self):
        r = self._preview("--stage", "align", "--input", self.sheet,
                          "--genome_fasta", self.genome, "--target_bed", self.bed,
                          "--outdir", "results/standalone/_preview")
        self.assertIn("SUCCESS", r.stdout + r.stderr, r.stderr)

    def test_coverage_and_call_entries_preview(self):
        # preview only builds channels (no heavy tools needed); dummy BAM+bai satisfy
        # the input-existence checks. Covers gatk + deepvariant calling paths.
        with tempfile.TemporaryDirectory() as d:
            for s in ("S1", "S2"):
                Path(os.path.join(d, f"{s}.analysis.bam")).write_bytes(b"")
                Path(os.path.join(d, f"{s}.analysis.bam.bai")).write_bytes(b"")
            glob = os.path.join(d, "*.analysis.bam")
            for extra in (["--stage", "coverage"],
                          ["--stage", "call"],
                          ["--stage", "call", "--caller", "deepvariant"],
                          ["--stage", "cnv"], ["--stage", "str"]):
                r = self._preview(*extra, "--input_bams", glob, "--genome_fasta", self.genome,
                                  "--target_bed", self.bed, "--outdir", "results/standalone/_preview")
                self.assertIn("SUCCESS", r.stdout + r.stderr, f"{extra}\n{r.stderr[-1500:]}")

    def test_cohortqc_and_giab_preview(self):
        with tempfile.TemporaryDirectory() as d:
            for s in ("CTRL", "S1", "S2"):
                Path(os.path.join(d, f"{s}.analysis.bam")).write_bytes(b"")
                Path(os.path.join(d, f"{s}.analysis.bam.bai")).write_bytes(b"")
            r1 = self._preview("--stage", "cohortqc", "--input_bams", os.path.join(d, "*.analysis.bam"),
                               "--input_vcf", self.filtered_vcf, "--genome_fasta", self.genome,
                               "--outdir", "results/standalone/_preview")
            self.assertIn("SUCCESS", r1.stdout + r1.stderr, r1.stderr[-1500:])
            r2 = self._preview("--stage", "giab", "--input_vcf", self.filtered_vcf,
                               "--target_bed", self.bed, "--genome_fasta", self.genome,
                               "--outdir", "results/standalone/_preview")
            self.assertIn("SUCCESS", r2.stdout + r2.stderr, r2.stderr[-1500:])

    def test_giab_runs_live_and_is_bed_restricted(self):
        # hap.py runs in its container; confirm on-target restriction from runinfo, not assertion.
        if not (os.path.exists(self.filtered_vcf) and shutil.which("docker")):
            self.skipTest("need results VCF + docker (hap.py image)")
        with tempfile.TemporaryDirectory() as d:
            out, work = os.path.join(d, "out"), os.path.join(d, "work")
            r = run([self.nf, "run", "main.nf", "-profile", "test,docker", "--stage", "giab",
                     "--input_vcf", self.filtered_vcf, "--target_bed", self.bed,
                     "--genome_fasta", self.genome, "--outdir", out, "-work-dir", work], cwd=REPO)
            if "SUCCESS" not in (r.stdout + r.stderr):
                self.skipTest(f"giab live run unavailable (docker/hap.py image?): {r.stderr[-400:]}")
            runinfo = Path(os.path.join(out, "stage13_giab", "happy.runinfo.json")).read_text()
            self.assertIn(os.path.basename(self.bed), runinfo)   # panel BED in hap.py's own cmdline
            self.assertRegex(runinfo, r"(-T|--target-regions)")

    def test_cohortqc_to_burden_chaining(self):
        # the cohort_qc.json a cohortqc run emits is the file `burden --cohort_qc_json` consumes
        if not (os.path.exists(self.cohort_qc_json) and os.path.exists(self.burden_vcf)):
            self.skipTest("need results/ cohort_qc.json + annotated VCF fixtures")
        import json
        with tempfile.TemporaryDirectory() as d:
            out, work = os.path.join(d, "out"), os.path.join(d, "work")
            r = run([self.nf, "run", "main.nf", "-profile", "test", "--stage", "burden",
                     "--input_vcf", self.burden_vcf, "--input", self.sheet, "--target_bed", self.bed,
                     "--cohort_qc_json", self.cohort_qc_json, "--outdir", out, "-work-dir", work], cwd=REPO)
            self.assertIn("SUCCESS", r.stdout + r.stderr, r.stderr[-2000:])
            meta = json.loads(Path(os.path.join(out, "stage14_burden", "burden_meta.json")).read_text())
            self.assertIn("n_pca_sites", meta)   # proves the cohortqc PCs were read/consumed

    def test_paralog_runs_live_and_writes_info_fields(self):
        # paralog uses bcftools (no conda needed): run it for real and confirm the
        # CaptureForge PARALOG_GENE/PARALOG_CONF labels land in the standalone VCF.
        if not (os.path.exists(self.filtered_vcf) and shutil.which("bcftools")):
            self.skipTest("need results/stage7_filter/joint.filtered.vcf.gz + bcftools")
        with tempfile.TemporaryDirectory() as d:
            out, work = os.path.join(d, "out"), os.path.join(d, "work")
            r = run([self.nf, "run", "main.nf", "-profile", "test", "--stage", "paralog",
                     "--input_vcf", self.filtered_vcf, "--target_bed", self.bed,
                     "--outdir", out, "-work-dir", work], cwd=REPO)
            self.assertIn("SUCCESS", r.stdout + r.stderr, r.stderr[-2000:])
            ov = os.path.join(out, "stage10_paralog", "paralog.annotated.vcf.gz")
            self.assertTrue(os.path.exists(ov))
            hdr = run(["bcftools", "view", "-h", ov]).stdout
            self.assertIn("ID=PARALOG_GENE", hdr)
            self.assertIn("ID=PARALOG_CONF", hdr)
            q = run(["bcftools", "query", "-f",
                     "%INFO/PARALOG_GENE\t%INFO/PARALOG_CONF\n", ov]).stdout
            self.assertTrue(any(tok not in (".", "") for line in q.splitlines()
                                for tok in line.split("\t")),
                            "no variant carried a PARALOG_GENE/PARALOG_CONF value")
            for f in ("paralog_report.md", "paralog_provenance.json"):
                self.assertTrue(os.path.exists(os.path.join(out, "stage10_paralog", f)), f)

    def test_burden_entry_runs_and_fails_loud(self):
        if not os.path.exists(self.burden_vcf):
            self.skipTest("results/ fixture missing (run the test pipeline first)")
        with tempfile.TemporaryDirectory() as d:
            out, work = os.path.join(d, "out"), os.path.join(d, "work")
            # real run (no conda/docker -> local PATH): pure-python + bcftools
            r = run([self.nf, "run", "main.nf", "-profile", "test", "--stage", "burden",
                     "--input_vcf", self.burden_vcf, "--input", self.sheet,
                     "--target_bed", self.bed, "--outdir", out, "-work-dir", work], cwd=REPO)
            self.assertIn("SUCCESS", r.stdout + r.stderr, r.stderr[-2000:])
            for f in ("burden_results.tsv", "burden_report.md", "burden_provenance.json"):
                self.assertTrue(os.path.exists(os.path.join(out, "stage14_burden", f)), f)
            self.assertTrue(os.path.exists(os.path.join(out, "stage14_burden", "plots", "burden_qq.png")))

            # fail-loud: phenotype-blanked sheet
            nopheno = os.path.join(d, "nopheno.csv")
            rows = list(csv.DictReader(Path(self.sheet).read_text().splitlines()))
            for row in rows:
                row["phenotype"] = ""
            with open(nopheno, "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=rows[0].keys()); w.writeheader(); w.writerows(rows)
            out2, work2 = os.path.join(d, "out2"), os.path.join(d, "work2")
            r2 = run([self.nf, "run", "main.nf", "-profile", "test", "--stage", "burden",
                      "--input_vcf", self.burden_vcf, "--input", nopheno,
                      "--target_bed", self.bed, "--outdir", out2, "-work-dir", work2], cwd=REPO)
            self.assertNotEqual(r2.returncode, 0)
            self.assertFalse(os.path.exists(os.path.join(out2, "stage14_burden", "burden_results.tsv")))


if __name__ == "__main__":
    unittest.main(verbosity=2)
