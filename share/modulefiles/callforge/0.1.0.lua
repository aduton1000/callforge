-- Lmod modulefile stub for CallForge (edit INSTALL_ROOT for your cluster).
-- OPTIONAL: only for sites that have Lmod. On a no-module cluster you do NOT need
-- this — put <install_root>/share/bin on PATH and run `callforge-run` directly.
help([[
CallForge 0.1.0 — targeted-capture germline analysis pipeline.
Loads Nextflow + a container engine (Apptainer/Singularity) on PATH and exposes `callforge-run`.
Run:  callforge-run -params-file site.params.yaml --input samplesheet.csv --slurm_partition <queue>
]])

whatis("Name: CallForge")
whatis("Version: 0.1.0")
whatis("Description: Targeted-capture germline analysis pipeline (Nextflow)")

local install_root = "/shared/apps/callforge"      -- EDIT for your cluster
local refs_root    = "/shared/refs/callforge/GRCh38"

-- Dependencies (adjust to your cluster's module names / container engine)
depends_on("nextflow")
depends_on("apptainer")   -- on a SingularityCE site use: depends_on("singularity")

setenv("CALLFORGE_HOME", install_root)
setenv("CALLFORGE_REFS", refs_root)
setenv("CALLFORGE_SIF", pathJoin(install_root, "callforge.sif"))
-- Engine the callforge-run wrapper composes into -profile hpc_slurm,<engine>.
-- The wrapper defaults to singularity; set apptainer here on an Apptainer site.
setenv("CALLFORGE_ENGINE", "apptainer")
setenv("CALLFORGE_CONTAINER_CACHE", pathJoin(os.getenv("HOME") or "/tmp", ".callforge/container_cache"))

prepend_path("PATH", pathJoin(install_root, "share/bin"))
prepend_path("PATH", pathJoin(install_root, "bin"))
