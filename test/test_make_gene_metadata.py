#!/usr/bin/env python3
"""test_make_gene_metadata.py — tests for bin/make_gene_metadata.py (stdlib-only).

Run:  python3 test/test_make_gene_metadata.py

Covers: the emitted schema is identical to ingest_captureforge.py's (so the stages
accept it); CNV callability is never fabricated (unknown, not yes/no); burden_group is
never auto-filled (blank unless declared); a declared gene absent from the BED is flagged
(soft); duplicates / malformed STR loci fail loud (hard).
"""
import csv, os, re, subprocess, sys, tempfile, unittest
from pathlib import Path

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(REPO, "bin", "make_gene_metadata.py")

def ingest_cols():
    src = Path(os.path.join(REPO, "bin", "ingest_captureforge.py")).read_text()
    return eval(re.search(r"cols = (\[.*?\])", src, re.S).group(1))

def run(args):
    p = subprocess.run([sys.executable, SCRIPT, *args], capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr

def write_bed(path):
    Path(path).write_text(
        "chr1\t100\t200\tCFH|cnv\nchr1\t300\t400\tCFH|coding\n"
        "chr1\t500\t600\tCD209|coding\nchr1\t700\t800\tHP|coding\n")


class MakeGeneMetadata(unittest.TestCase):
    def test_schema_matches_ingest_and_no_fabricated_fields(self):
        with tempfile.TemporaryDirectory() as d:
            bed = os.path.join(d, "panel.bed"); write_bed(bed)
            bg = os.path.join(d, "bg.tsv"); Path(bg).write_text("CFH\tCOMPLEMENT\n")
            out = os.path.join(d, "gm.tsv")
            rc, so, se = run(["--bed", bed, "--cnv-genes", "CFH", "--paralog-genes", "CD209",
                              "--burden-groups", bg, "--out-tsv", out])
            self.assertEqual(rc, 0, se)
            with open(out, newline="") as fh:
                rows = list(csv.DictReader(fh, delimiter="\t"))
            # identical schema to ingest_captureforge.py
            self.assertEqual(list(rows[0].keys()), ingest_cols())
            by = {r["gene"]: r for r in rows}
            # CNV callability NEVER fabricated yes/no for a user-declared target
            self.assertEqual(by["CFH"]["is_cnv_target"], "yes")
            self.assertEqual(by["CFH"]["cnv_callable"], "unknown")
            self.assertEqual(by["CFH"]["cnv_callable_reason"], "user_declared_not_design_derived")
            # paralog flag honoured
            self.assertEqual(by["CD209"]["is_paralog"], "yes")
            # burden_group: declared for CFH, BLANK (never auto) for the rest
            self.assertEqual(by["CFH"]["burden_group"], "COMPLEMENT")
            self.assertEqual(by["CD209"]["burden_group"], "")
            self.assertEqual(by["HP"]["burden_group"], "")
            # design-derived fields not invented
            self.assertEqual(by["CFH"]["off_target_frac"], "")
            self.assertEqual(by["CFH"]["low_coverage"], "")

    def test_declared_gene_not_in_bed_is_flagged_soft(self):
        with tempfile.TemporaryDirectory() as d:
            bed = os.path.join(d, "panel.bed"); write_bed(bed)
            out = os.path.join(d, "gm.tsv"); rep = os.path.join(d, "rep.txt")
            rc, so, se = run(["--bed", bed, "--cnv-genes", "NOTINPANEL",
                              "--out-tsv", out, "--report", rep])
            self.assertEqual(rc, 0, se)                 # soft: still writes
            self.assertIn("DECLARED GENE NOT IN TARGET BED", Path(rep).read_text())
            self.assertIn("NOTINPANEL", Path(rep).read_text())

    def test_duplicate_burden_group_hard_fails(self):
        with tempfile.TemporaryDirectory() as d:
            bed = os.path.join(d, "panel.bed"); write_bed(bed)
            bg = os.path.join(d, "bg.tsv"); Path(bg).write_text("CFH\tA\nCFH\tB\n")
            out = os.path.join(d, "gm.tsv")
            rc, so, se = run(["--bed", bed, "--burden-groups", bg, "--out-tsv", out])
            self.assertNotEqual(rc, 0)
            self.assertIn("duplicate", se.lower())
            self.assertFalse(os.path.exists(out))

    def test_malformed_str_locus_hard_fails(self):
        with tempfile.TemporaryDirectory() as d:
            bed = os.path.join(d, "panel.bed"); write_bed(bed)
            sl = os.path.join(d, "sl.tsv"); Path(sl).write_text("CD209\tnot_a_locus\n")
            out = os.path.join(d, "gm.tsv")
            rc, so, se = run(["--bed", bed, "--str-loci", sl, "--out-tsv", out])
            self.assertNotEqual(rc, 0)
            self.assertIn("malformed locus", se.lower())

    def test_stages_accept_output_cnv_label(self):
        # the generated metadata must drive cnv_annotate.py exactly like ingest's
        with tempfile.TemporaryDirectory() as d:
            bed = os.path.join(d, "panel.bed"); write_bed(bed)
            gm = os.path.join(d, "gm.tsv")
            rc, *_ = run(["--bed", bed, "--cnv-genes", "CFH", "--out-tsv", gm])
            self.assertEqual(rc, 0)
            cns = os.path.join(d, "S.call.cns")
            Path(cns).write_text("chromosome\tstart\tend\tgene\tlog2\tcn\nchr1\t100\t400\tCFH\t-0.7\t1\n")
            r = subprocess.run([sys.executable, os.path.join(REPO, "bin", "cnv_annotate.py"),
                                "--cns", cns, "--bed", bed, "--metadata", gm, "--outdir", d],
                               capture_output=True, text=True)
            self.assertEqual(r.returncode, 0, r.stderr)
            calls = Path(os.path.join(d, "cnv_calls.tsv")).read_text()
            # generated cnv_callable=unknown -> cnv_annotate confidence 'unknown_confidence'
            self.assertIn("unknown_confidence", calls)


if __name__ == "__main__":
    unittest.main(verbosity=2)
