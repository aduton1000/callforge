#!/usr/bin/env bash
# build_manual.sh — render the CallForge user manual to PDF + DOCX from Markdown.
#
#   bash bin/build_manual.sh
#
# Requires: pandoc, a LaTeX engine (xelatex; install TinyTeX/texlive if absent),
# and graphviz `dot` (to (re)render the DAG). Outputs land in docs/manual/.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/../docs/manual" && pwd)"
cd "$HERE"

MD="callforge_manual.md"
META="metadata.yaml"
DOT_BIN="${DOT_BIN:-$(command -v dot || true)}"
PANDOC="${PANDOC:-$(command -v pandoc)}"

# 1. (re)render the DAG (PNG + SVG) if dot is available
if [ -n "$DOT_BIN" ]; then
  "$DOT_BIN" -Tpng -Gdpi=200 callforge_dag.dot -o figures/callforge_dag.png
  "$DOT_BIN" -Tsvg            callforge_dag.dot -o figures/callforge_dag.svg
  echo "[build_manual] DAG rendered (figures/callforge_dag.png|svg)"
else
  echo "[build_manual] WARN: graphviz 'dot' not found; using existing DAG image" >&2
fi

# 2. PDF (xelatex; report class -> title page on its own page, TOC on the next)
if "$PANDOC" --version >/dev/null 2>&1 && command -v xelatex >/dev/null 2>&1; then
  "$PANDOC" "$MD" \
    --metadata-file="$META" \
    --pdf-engine=xelatex \
    --toc --toc-depth=3 \
    -o callforge_manual.pdf
  echo "[build_manual] PDF -> callforge_manual.pdf"
else
  echo "[build_manual] WARN: pandoc/xelatex missing; trying weasyprint fallback" >&2
  "$PANDOC" "$MD" --metadata-file="$META" --toc -t html5 -o /tmp/_cf_manual.html
  weasyprint /tmp/_cf_manual.html callforge_manual.pdf
fi

# 3. DOCX (title block + TOC front matter; body page-broken)
"$PANDOC" "$MD" \
  --metadata-file="$META" \
  --toc --toc-depth=3 \
  -o callforge_manual.docx
echo "[build_manual] DOCX -> callforge_manual.docx"

echo "[build_manual] done:"
ls -lh callforge_manual.md callforge_manual.pdf callforge_manual.docx 2>/dev/null | awk '{print "  "$NF, $5}'
