"""callforge_cli — the unified `callforge` command-line front door.

This is a thin wrapper around the existing CallForge entry points (Nextflow
`main.nf`, `bin/init_sample_sheet.py`, `bin/discover_resources.py`). It does NOT
re-architect the pipeline and does NOT collapse the per-stage conda environments:
when the pipeline runs under `-profile conda`/`apptainer`, Nextflow activates each
stage's isolated `env/*.yml` / container automatically. The CLI just launches it.
"""
__version__ = "0.1.0"
