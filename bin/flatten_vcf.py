#!/usr/bin/env python3
"""flatten_vcf.py — flatten the annotated joint VCF to a per-variant TSV (CallForge).

One row per variant with the headline annotations for downstream R/pandas review:
CHROM, POS, REF, ALT, rsID, FILTER, gene + most-severe consequence (from VEP CSQ),
gnomAD AF + AFR AF, ClinVar significance, PhyloP, and the paralog flag.
Missing annotations render as '.'.
"""
import argparse, os, re, subprocess

# VEP severity order (most severe first) for picking a single consequence per variant
SEVERITY = ["transcript_ablation", "splice_acceptor_variant", "splice_donor_variant",
            "stop_gained", "frameshift_variant", "stop_lost", "start_lost",
            "transcript_amplification", "inframe_insertion", "inframe_deletion",
            "missense_variant", "protein_altering_variant", "splice_region_variant",
            "incomplete_terminal_codon_variant", "start_retained_variant",
            "stop_retained_variant", "synonymous_variant", "coding_sequence_variant",
            "mature_miRNA_variant", "5_prime_UTR_variant", "3_prime_UTR_variant",
            "non_coding_transcript_exon_variant", "intron_variant",
            "NMD_transcript_variant", "non_coding_transcript_variant",
            "upstream_gene_variant", "downstream_gene_variant", "regulatory_region_variant",
            "intergenic_variant"]
RANK = {c: i for i, c in enumerate(SEVERITY)}


def header(vcf):
    return subprocess.run(["bcftools", "view", "-h", vcf], capture_output=True, text=True).stdout


def csq_fields(h):
    m = re.search(r"##INFO=<ID=CSQ.*?Format: ([^\">]+)", h)
    return m.group(1).split("|") if m else []


def info_tags(h):
    return set(re.findall(r"##INFO=<ID=([^,]+)", h))


def pick_consequence(csq, fmt):
    """From a CSQ string (comma-sep transcripts), pick gene + most-severe consequence."""
    if not csq or csq == "." or not fmt:
        return ".", "."
    ci = fmt.index("Consequence") if "Consequence" in fmt else 1
    gi = fmt.index("SYMBOL") if "SYMBOL" in fmt else (fmt.index("Gene") if "Gene" in fmt else 0)
    best, best_rank, best_gene = ".", 1e9, "."
    for tx in csq.split(","):
        parts = tx.split("|")
        if len(parts) <= max(ci, gi):
            continue
        for cons in parts[ci].split("&"):
            r = RANK.get(cons, 500)
            if r < best_rank:
                best_rank, best, best_gene = r, cons, parts[gi] or "."
    return best, best_gene


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--gnomad-af-field", default="AF_afr")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    h = header(a.vcf)
    fmt_csq = csq_fields(h)
    present = info_tags(h)
    afr = f"gnomAD_{a.gnomad_af_field}"
    # Optional INFO tags, queried only if defined in the header (bcftools query
    # errors on an undefined tag); absent ones render as '.'.
    opt = ["CSQ", "gnomAD_AF", afr, "ClinVar_CLNSIG", "PhyloP", "PARALOG_GENE", "PARALOG_CONF"]
    use = [t for t in opt if t in present]
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%ID\t%FILTER" + "".join(f"\t%INFO/{t}" for t in use) + "\n"
    rows = subprocess.run(["bcftools", "query", "-f", fmt, a.vcf],
                          capture_output=True, text=True).stdout.splitlines()
    cols = ["chrom", "pos", "ref", "alt", "rsID", "filter", "gene", "consequence",
            "gnomAD_AF", afr, "ClinVar_CLNSIG", "PhyloP", "paralog_gene", "paralog_conf"]
    with open(a.out, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            f = r.split("\t")
            base = f[:6]
            vals = dict(zip(use, f[6:]))
            cons, gene = pick_consequence(vals.get("CSQ", "."), fmt_csq)
            fh.write("\t".join(base[:6] + [gene, cons,
                     vals.get("gnomAD_AF", "."), vals.get(afr, "."),
                     vals.get("ClinVar_CLNSIG", "."), vals.get("PhyloP", "."),
                     vals.get("PARALOG_GENE", "."), vals.get("PARALOG_CONF", ".")]) + "\n")
    print(f"[flatten_vcf] wrote {len(rows)} variants -> {a.out}")


if __name__ == "__main__":
    main()
