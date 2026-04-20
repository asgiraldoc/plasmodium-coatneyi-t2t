#!/usr/bin/env bash
# prepare_chromsyn_inputs.sh

LINEAGE="/path/to/busco_downloads/lineages/plasmodium_odb12"
cd ../genomes || exit 1

echo "[INFO] Running BUSCO..."
for genome in *.fasta; do
    PREFIX=$(basename "$genome" .fasta)
    busco -i "$genome" -o "run_$PREFIX" -m genome -l "$LINEAGE" --cpu 30
    cp "run_${PREFIX}/run_plasmodium_odb12/full_table.tsv" "../busco_runs/${PREFIX}.busco5.tsv"
done

echo "[INFO] Generating FOFN (File Of File Names) lists..."
ls ../telociraptor_results/*.telomeres.tdt | sed 's#.*/##; s/\.telomeres\.tdt//' | \
  while read id; do echo "$id telociraptor_results/${id}.telomeres.tdt"; done > ../fofn/sequences.fofn

ls ../telociraptor_results/*.gaps.tdt | sed 's#.*/##; s/\.gaps\.tdt//' | \
  while read id; do echo "$id telociraptor_results/${id}.gaps.tdt"; done > ../fofn/gaps.fofn

ls ../busco_runs/*.busco5.tsv | sed 's#.*/##; s/\.busco5\.tsv//' | \
  while read id; do echo "$id busco_runs/${id}.busco5.tsv"; done > ../fofn/busco.fofn

ls ../tidk_results/*.tidk.csv | sed 's#.*/##; s/\.tidk\.csv//' | \
  while read id; do echo "$id tidk_results/${id}.tidk.csv"; done > ../fofn/tidk.fofn

echo "[INFO] Preparation complete."
