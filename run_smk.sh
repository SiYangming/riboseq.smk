#!/bin/bash
set -e

# 1. Environment Configuration
echo "Configuring environment..."

# Source bash_profile as requested
if [ -f ~/.bash_profile ]; then
    source ~/.bash_profile
fi

# Activate mamba environment
# Ensuring mamba is initialized for this shell session
if command -v mamba &> /dev/null; then
	eval "$(mamba shell hook --shell bash)"
fi

# Check if environment exists
if mamba info --envs | grep -q "snakemake"; then
    mamba activate snakemake
else
    echo "Warning: 'snakemake' mamba environment not found. Attempting to use current environment."
fi

# Check requirements
echo "Checking requirements..."
if ! command -v snakemake &> /dev/null; then
	echo "Error: snakemake command not found."
	exit 1
fi

echo "Starting Snakemake workflow (dry-run)..."
echo "Input config: config/config.yaml"
echo "Working directory: $(pwd)"

SNAKEMAKE_OPTS=""

for arg in "$@"; do
	if [ "$arg" == "--resume" ]; then
		echo "Resume mode enabled: adding --rerun-incomplete"
		SNAKEMAKE_OPTS="$SNAKEMAKE_OPTS --rerun-incomplete"
	else
		SNAKEMAKE_OPTS="$SNAKEMAKE_OPTS $arg"
	fi
done

snakemake -s workflow/Snakefile \
	--configfile config/config.yaml \
	-c 4 \
	-p \
	--latency-wait 60 \
	-n \
	$SNAKEMAKE_OPTS
