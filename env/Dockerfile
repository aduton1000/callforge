# ═══════════════════════════════════════════════════════════════════════════
#  CallForge — Docker image (Mac/CI development) for the CORE environment.
#  Build:  docker build -t callforge:0.1.0 -f env/Dockerfile .
#  Specialized stages (VEP, somalier, hap.py, regenie, DeepVariant) use their
#  own conda envs / official containers — see env/*.yml and docs/.
#  linux-64 builds are native in the image (no emulation), so callforge.yml
#  resolves cleanly even when authored on Apple Silicon.
# ═══════════════════════════════════════════════════════════════════════════
FROM mambaorg/micromamba:1.5.8

LABEL org.opencontainers.image.title="CallForge" \
      org.opencontainers.image.description="Targeted-capture germline analysis pipeline (core env)" \
      org.opencontainers.image.version="0.1.0" \
      org.opencontainers.image.licenses="MIT"

USER root
COPY env/callforge.yml /tmp/callforge.yml
RUN micromamba install -y -n base -f /tmp/callforge.yml && \
    micromamba clean --all --yes

ENV PATH=/opt/conda/bin:$PATH
ENV MAMBA_DOCKERFILE_ACTIVATE=1

# Sanity-check the core toolchain is on PATH at build time.
RUN fastqc --version && fastp --version 2>&1 | head -1 && bwa-mem2 version 2>&1 | head -1 && \
    samtools --version | head -1 && bcftools --version | head -1 && \
    gatk --version 2>&1 | head -3 && picard MarkDuplicates --version 2>&1 | head -1 || true && \
    mosdepth --version && multiqc --version && \
    python -c "import pandas, numpy, scipy, matplotlib, allel, pysam, pyfaidx, jinja2; print('py libs OK')"

WORKDIR /work
CMD ["/bin/bash"]
