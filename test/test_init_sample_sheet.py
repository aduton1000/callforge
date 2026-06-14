#!/usr/bin/env python3
"""test_init_sample_sheet.py — tests for bin/init_sample_sheet.py (stdlib-only).

Run:  python3 test/test_init_sample_sheet.py   (or: python3 -m unittest -v ...)

Covers:
  - a CLEAN case (Mode B) produces a valid CallForge sheet (exit 0, right columns,
    absolute FASTQ paths, joined phenotype);
  - a case with a deliberate ORPHAN FASTQ and a deliberate ID MISMATCH gets both
    flagged in the reconciliation report, and phenotype is never auto-filled.
"""
import csv, os, subprocess, sys, tempfile, unittest

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT = os.path.join(HERE, "..", "bin", "init_sample_sheet.py")


def touch(path):
    with open(path, "w") as fh:
        fh.write("@read\nACGT\n+\nFFFF\n")


def run(args):
    p = subprocess.run([sys.executable, SCRIPT, *args],
                       capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


class CleanCase(unittest.TestCase):
    def test_mode_b_clean(self):
        with tempfile.TemporaryDirectory() as d:
            fq = os.path.join(d, "fastq"); os.makedirs(fq)
            for s in ("S1", "S2"):
                touch(os.path.join(fq, f"{s}_R1.fastq.gz"))
                touch(os.path.join(fq, f"{s}_R2.fastq.gz"))
            meta = os.path.join(d, "clinical.csv")
            with open(meta, "w") as fh:
                fh.write("sample_id,sex,phenotype,covariate_age,batch\n"
                         "S1,M,case,45,b1\nS2,F,control,52,b2\n")
            out = os.path.join(d, "sheet.csv")
            rc, so, se = run(["--fastq-dir", fq, "--metadata", meta, "--out", out])
            self.assertEqual(rc, 0, f"clean case should pass; stderr=\n{se}")
            self.assertTrue(os.path.isfile(out))
            with open(out, newline="") as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual([r["sample_id"] for r in rows], ["S1", "S2"])
            # required schema present, paths absolute, phenotype joined (not invented)
            for c in ("sample_id", "fastq_1", "fastq_2", "sex", "phenotype",
                      "batch", "covariate_age"):
                self.assertIn(c, rows[0])
            self.assertTrue(rows[0]["fastq_1"].startswith("/"))
            self.assertEqual({r["sample_id"]: r["phenotype"] for r in rows},
                             {"S1": "case", "S2": "control"})


class OrphanAndMismatch(unittest.TestCase):
    def test_flags_orphan_and_id_mismatch(self):
        with tempfile.TemporaryDirectory() as d:
            fq = os.path.join(d, "fastq"); os.makedirs(fq)
            # paired sample
            touch(os.path.join(fq, "S1_R1.fastq.gz"))
            touch(os.path.join(fq, "S1_R2.fastq.gz"))
            # deliberate ORPHAN: R1 with no R2 mate
            touch(os.path.join(fq, "S3_R1.fastq.gz"))
            meta = os.path.join(d, "clinical.csv")
            # S1 matches; S99 is metadata-only (deliberate ID mismatch / no FASTQ);
            # S3 has NO phenotype row -> stays blank (never auto-filled)
            with open(meta, "w") as fh:
                fh.write("sample_id,sex,phenotype\nS1,M,case\nS99,F,control\n")
            out = os.path.join(d, "sheet.csv")
            report = os.path.join(d, "rep.txt")
            rc, so, se = run(["--fastq-dir", fq, "--metadata", meta,
                              "--out", out, "--report", report])
            self.assertEqual(rc, 0, f"soft issues should not hard-fail; stderr=\n{se}")
            self.assertTrue(os.path.isfile(report))
            with open(report) as fh:
                rep = fh.read()
            # orphan flagged
            self.assertIn("ORPHAN FASTQ", rep)
            self.assertIn("S3_R1.fastq.gz", rep)
            # metadata-only / ID mismatch flagged
            self.assertIn("METADATA WITHOUT SEQUENCE", rep)
            self.assertIn("S99", rep)
            # sequence-without-metadata flagged for the orphan-free paired sample S3? no:
            # S3 had no pair, so it is an orphan, not a written row. S1 is fully joined.
            # phenotype never invented: the written sheet has S1=case only, no S3/S99 row
            with open(out, newline="") as fh:
                rows = list(csv.DictReader(fh))
            self.assertEqual([r["sample_id"] for r in rows], ["S1"])
            self.assertEqual(rows[0]["phenotype"], "case")

    def test_duplicate_id_hard_fails(self):
        with tempfile.TemporaryDirectory() as d:
            fq = os.path.join(d, "fastq"); os.makedirs(fq)
            # two lanes of the same sample -> same inferred sample_id -> duplicate
            for lane in ("L001", "L002"):
                touch(os.path.join(fq, f"S1_{lane}_R1_001.fastq.gz"))
                touch(os.path.join(fq, f"S1_{lane}_R2_001.fastq.gz"))
            out = os.path.join(d, "sheet.csv")
            rc, so, se = run(["--fastq-dir", fq, "--out", out])
            self.assertNotEqual(rc, 0, "duplicate sample_id must hard-fail")
            self.assertIn("duplicate", se.lower())
            self.assertFalse(os.path.isfile(out), "no sheet written on hard error")


if __name__ == "__main__":
    unittest.main(verbosity=2)
