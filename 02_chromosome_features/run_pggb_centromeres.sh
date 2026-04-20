#!/usr/bin/env bash
# run_pggb_centromeres.sh
# Description: Builds a pangenome graph using PGGB to align chromosomes.

for i in *1[2-4].fasta; do
    base=$(basename "$i" .fasta)
    echo "[INFO] Processing $base..."
    
    pggb \
        -i "$i" \
        -o "${base}_p85_s5kb" \
        -n 10 \
        -t 30 \
        -p 85 \
        -s 5000 \
        -S

    echo "[INFO] Finished PGGB build for: $base"
done
