#!/usr/bin/env bash
# run_mumemto.sh
# Description: Generates a multi-MUM synteny plot for the input assemblies.

INPUT_FASTAS="*.fasta"
OUTPUT_DIR="output/assemblies"
FIGURES_DIR="figures"

mkdir -p "$OUTPUT_DIR" "$FIGURES_DIR"

echo "[INFO] Running MUMemto for sequence alignment..."
mumemto $INPUT_FASTAS -o "$OUTPUT_DIR"

echo "[INFO] Generating synteny visualization..."
mumemto viz \
    -i "$OUTPUT_DIR" \
    --filelist filelist.txt \
    --labels labels.txt \
    --mode gapped \
    --spacer 0.1 \
    -o "$FIGURES_DIR/assemblies.png"

echo "[INFO] Synteny visualization complete."
