#!/usr/bin/env bash
# run_chromsyn.sh

MAX_JOBS=30  # Maximum number of simultaneous processes

# Function to control parallel jobs
function limit_jobs {
  while [ "$(jobs -rp | wc -l)" -ge "$MAX_JOBS" ]; do
    sleep 1
  done
}

SCRIPT="/path/to/chromsyn/chromsyn.R"
BASE_DIR="../fofn"

SEQUENCES="${BASE_DIR}/sequences.fofn"
GAPS="${BASE_DIR}/gaps.fofn"
BUSCO="${BASE_DIR}/busco.fofn"
TIDK="${BASE_DIR}/tidk.fofn"
FEATURES="${BASE_DIR}/features.fofn"

TICKS="0.1e7"
LABELS="F"
NAMESIZE="0"
YPAD="0.25"
PDFWIDTH="15"
PDFHEIGHT="8.3"
CHROMFILL="Clade"
OPACITY="0.15"

echo "[INFO] Launching ChromSyn plotting jobs..."

for YGAP in 2; do
  for FTSIZE in 1.5; do
    limit_jobs

    BASE="chromsyn_ygap${YGAP}_ftsize${FTSIZE}"
    LOG="${BASE}.log"

    echo "[INFO] Launching: $BASE"

    Rscript "$SCRIPT" \
      basefile="$BASE" \
      sequences="$SEQUENCES" \
      ft="$FEATURES" \
      gaps="$GAPS" \
      busco="$BUSCO" \
      tidk="$TIDK" \
      ticks="$TICKS" \
      labels="$LABELS" \
      namesize="$NAMESIZE" \
      ypad="$YPAD" \
      ygap="$YGAP" \
      pdfwidth="$PDFWIDTH" \
      pdfheight="$PDFHEIGHT" \
      chromfill="$CHROMFILL" \
      ftsize="$FTSIZE" \
      opacity="$OPACITY" \
      > "$LOG" 2>&1 &
  done
done

# Wait for all background processes to finish
wait
echo "[INFO] All plotting processes completed successfully."
