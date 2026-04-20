#!/usr/bin/env bash
# extract_centromeres.sh
# Description: Extracts centromere subgraphs based on coordinates and calculates GC content.

set -euo pipefail

OUTPUT_DIR="centromere_extraction"
mkdir -p "$OUTPUT_DIR"

echo "[INFO] Starting centromere extraction..."

# Extended coordinates for each chromosome (Based on GC troughs)
declare -A START END
START=( [01]=675962  [02]=160348  [03]=628663  [04]=562965
        [05]=842198  [06]=333697  [07]=1254418 [08]=1130374
        [09]=975411  [10]=1009338 [11]=1601556 [12]=1046749
        [13]=983315  [14]=2141532 )
END=(   [01]=681789  [02]=166087  [03]=634481  [04]=568927
        [05]=848072  [06]=339569  [07]=1260299 [08]=1136202
        [09]=981222  [10]=1015167 [11]=1607511 [12]=1052261
        [13]=989034  [14]=2147269 )

for chr in $(seq -w 01 14); do
    s=${START[$chr]}
    e=${END[$chr]}
    chr_file="chr${chr}.fasta"
    base_dir="chr${chr}_p90_s5kb"

    echo "[INFO] Processing $chr_file (Range: $s-$e)..."

    # Verify graph existence
    og_file=$(ls "${base_dir}/${chr_file}."*.smooth.final.og 2>/dev/null || true)
    if [[ -z "$og_file" ]]; then
        echo "[WARNING] OG file not found for chr${chr} in $base_dir. Skipping..."
        continue
    fi

    # Find exact path corresponding to the chromosome
    ref_path=$(odgi paths -i "$og_file" -L | grep "PvP01" | grep "#1#${chr}" || true)

    if [[ -z "$ref_path" ]]; then
        echo "[WARNING] No valid path found for PvP01 in chr${chr}. Skipping..."
        continue
    fi

    CHR_OUT="$OUTPUT_DIR/chr${chr}"
    mkdir -p "$CHR_OUT"

    echo "[INFO] Extracting centromere subgraph..."
    odgi extract -i "$og_file" -r "${ref_path}:${s}-${e}" -o "$CHR_OUT/centromere_subgraph.og"

    echo "[INFO] Calculating subgraph statistics..."
    odgi stats -i "$CHR_OUT/centromere_subgraph.og" -S > "$CHR_OUT/centromere_stats.txt"

    echo "[INFO] Generating BED coordinates..."
    odgi view -i "$CHR_OUT/centromere_subgraph.og" -g \
        | grep -e "#${chr}" \
        | awk '{print $2}' \
        | sed 's/:\|-/\t/g' > "$CHR_OUT/centromere_subgraph.coords.bed"

    echo "[INFO] Extracting original sequence..."
    seqkit subseq --bed "$CHR_OUT/centromere_subgraph.coords.bed" "$chr_file" -o "$CHR_OUT/chr${chr}.centromere.fasta"

    echo "[INFO] Calculating GC content and length..."
    seqkit fx2tab -g -l -n "$CHR_OUT/chr${chr}.centromere.fasta" > "$CHR_OUT/chr${chr}.centromere_gc.tsv"

    echo "[INFO] Finished chr${chr} -> $CHR_OUT"
done

echo "[INFO] Extraction complete. Check results in: $OUTPUT_DIR"
