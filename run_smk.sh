#!/usr/bin/env bash
# riboseq.smk runner — exec_mode from config/config.yaml: native | conda | container
set -euo pipefail

SNAKEMAKE_CMD="snakemake"

if [[ -f "Snakefile" && -f "config.yaml" ]]; then
    SNAKEFILE="Snakefile"
    CONFIG="config.yaml"
elif [[ -f "workflow/Snakefile" && -f "config/config.yaml" ]]; then
    SNAKEFILE="workflow/Snakefile"
    CONFIG="config/config.yaml"
else
    echo "[ERROR] Snakefile / config.yaml not found" >&2
    exit 1
fi

EXEC_MODE=$(grep -E '^exec_mode:' "$CONFIG" | head -n 1 | sed -E 's/.*:[[:space:]]*"?([^"]*)"?.*/\1/')
EXEC_MODE="${EXEC_MODE:-conda}"
echo "Working directory: $(pwd)"
echo "Snakefile: $SNAKEFILE | Config: $CONFIG | exec_mode=$EXEC_MODE"

SNAKEMAKE_OPTS=""
case "$EXEC_MODE" in
    native)    : ;;
    conda)     SNAKEMAKE_OPTS="$SNAKEMAKE_OPTS --use-conda" ;;
    container) SNAKEMAKE_OPTS="$SNAKEMAKE_OPTS --sdm apptainer" ;;
    *)
        echo "[ERROR] unknown exec_mode: $EXEC_MODE (native|conda|container)" >&2
        exit 1
        ;;
esac

for arg in "$@"; do
    if [[ "$arg" == "--resume" ]]; then
        SNAKEMAKE_OPTS="$SNAKEMAKE_OPTS --rerun-incomplete"
    else
        SNAKEMAKE_OPTS="$SNAKEMAKE_OPTS $arg"
    fi
done

$SNAKEMAKE_CMD -s "$SNAKEFILE" \
    --configfile "$CONFIG" \
    -c all -p \
    --latency-wait 60 \
    $SNAKEMAKE_OPTS
