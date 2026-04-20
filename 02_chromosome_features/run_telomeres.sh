#!/usr/bin/env bash
# run_telomeres.sh
# Description: Identifies telomeric repeats using Telociraptor and TIDK.

GENOME_DIR="../genomes"
TELOCIRAPTOR_DIR="../telociraptor_results"
TIDK_DIR="../tidk_results"
TELOMERES_PLOT_DIR="../telomeres_plots"

mkdir -p "$TELOCIRAPTOR_DIR" "$TIDK_DIR" "$TELOMERES_PLOT_DIR"

# Ensure the Telociraptor python script is accessible in your environment
TELOCIRAPTOR_SCRIPT="/path/to/telociraptor/code/telociraptor.py"

cd "$GENOME_DIR" || exit 1

for genome in *.fasta; do
    PREFIX=$(basename "$genome" .fasta)
    echo "[INFO] Processing telomeres for $PREFIX..."

    # 1. Telociraptor Analysis
    python "$TELOCIRAPTOR_SCRIPT" \
        seqin="$genome" \
        basefile="${TELOCIRAPTOR_DIR}/${PREFIX}" \
        i=-1 tweak=F telonull=T

    # 2. TIDK Analysis (2kb window, targeting Haemosporida motif)
    tidk find -c Haemosporida -w 2000 -d "$TIDK_DIR" -o "$PREFIX" "$genome"

    # Rename the TSV output to .csv for downstream ChromSyn compatibility
    cp "${TIDK_DIR}/${PREFIX}_telomeric_repeat_windows.tsv" \
       "${TIDK_DIR}/${PREFIX}.tidk.csv"

    # 3. Plot detected telomeres
    tidk plot \
        --tsv "${TIDK_DIR}/${PREFIX}_telomeric_repeat_windows.tsv" \
        -o "${TELOMERES_PLOT_DIR}/${PREFIX}_telomere"
done

echo "[INFO] Telomere identification complete."
