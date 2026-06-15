#!/usr/bin/env python3
"""test_callforge_cli.py — tests for the unified `callforge` CLI (stdlib-only).

Run:  python3 test/test_callforge_cli.py

Invokes the CLI as `python3 -m callforge_cli.cli …` from the repo root (no install
needed). Checks that --help / --version and the subcommand help work, and that
`callforge run --print-cmd …` constructs the expected `nextflow run` command without
executing Nextflow.
"""
import os
import subprocess
import sys
import unittest

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def cli(*args):
    env = dict(os.environ, PYTHONPATH=REPO)
    p = subprocess.run([sys.executable, "-m", "callforge_cli.cli", *args],
                       capture_output=True, text=True, cwd=REPO, env=env)
    return p.returncode, p.stdout, p.stderr


class CliHelp(unittest.TestCase):
    def test_top_help_lists_subcommands(self):
        rc, so, se = cli("--help")
        self.assertEqual(rc, 0)
        for cmd in ("run", "samplesheet", "resources"):
            self.assertIn(cmd, so)

    def test_version(self):
        rc, so, se = cli("--version")
        self.assertEqual(rc, 0)
        self.assertEqual(so.strip(), "0.1.0")

    def test_samplesheet_help(self):
        rc, so, se = cli("samplesheet", "--help")
        self.assertEqual(rc, 0)
        self.assertIn("--fastq-dir", so)

    def test_run_help(self):
        rc, so, se = cli("run", "--help")
        self.assertEqual(rc, 0)
        self.assertIn("nextflow", so.lower())


class CliRunCmd(unittest.TestCase):
    def test_run_constructs_nextflow_command(self):
        rc, so, se = cli("run", "--print-cmd",
                         "-profile", "test,conda",
                         "-params-file", "params.full.yaml",
                         "--input", "sample_sheet.csv")
        self.assertEqual(rc, 0, f"stderr=\n{se}")
        # the echoed command must be a real `nextflow run <main.nf> …` with passthrough
        self.assertIn("nextflow run", so)
        self.assertIn("main.nf", so)
        self.assertIn("-profile test,conda", so)
        self.assertIn("-params-file params.full.yaml", so)
        self.assertIn("--input sample_sheet.csv", so)

    def test_run_pipeline_override_passthrough(self):
        rc, so, se = cli("run", "--print-cmd",
                         "--pipeline", "aduton1000/callforge",
                         "-r", "v0.1.0", "-profile", "hpc_slurm,apptainer")
        self.assertEqual(rc, 0, f"stderr=\n{se}")
        self.assertIn("nextflow run aduton1000/callforge", so)
        self.assertIn("-r v0.1.0", so)


if __name__ == "__main__":
    unittest.main(verbosity=2)
