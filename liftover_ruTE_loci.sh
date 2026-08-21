#!/usr/bin/env bash
#SBATCH --job-name=rute_liftover
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

# Lift hg38 ruTE loci to multiple mammalian genomes and build a presence matrix.
#
# Usage:
#   sbatch liftover_ruTE_loci.sh [input_bed] [chain_directory] [output_directory]
#
# Defaults:
#   input_bed       <script_directory>/ruTE_loci_n838_with_milliDiv_and_MYA.bed
#   chain_directory <script_directory>/chain_files
#   output_directory <script_directory>/liftover_results
#
# Add cluster-specific SBATCH options, such as --partition, --qos, --time, and
# --mem, either above or on the sbatch command line.

set -euo pipefail

SCRIPT_DIRECTORY="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
INPUT_BED="${1:-${SCRIPT_DIRECTORY}/ruTE_loci_n838_with_milliDiv_and_MYA.bed}"
CHAIN_DIRECTORY="${2:-${SCRIPT_DIRECTORY}/chain_files}"
OUTPUT_DIRECTORY="${3:-${SCRIPT_DIRECTORY}/liftover_results}"
MAX_PARALLEL="${SLURM_CPUS_PER_TASK:-1}"

if [[ $# -gt 3 ]]; then
    echo "Usage: sbatch $0 [input_bed] [chain_directory] [output_directory]" >&2
    exit 2
fi

if ! command -v liftOver >/dev/null 2>&1; then
    echo "Error: liftOver is not available in PATH." >&2
    exit 1
fi
if [[ ! -f "$INPUT_BED" ]]; then
    echo "Error: input BED file not found: $INPUT_BED" >&2
    exit 1
fi
if [[ ! -d "$CHAIN_DIRECTORY" ]]; then
    echo "Error: chain directory not found: $CHAIN_DIRECTORY" >&2
    exit 1
fi
if [[ ! "$MAX_PARALLEL" =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: SLURM_CPUS_PER_TASK must be a positive integer." >&2
    exit 1
fi

shopt -s nullglob
CHAIN_FILES=("${CHAIN_DIRECTORY}"/hg38To*.over.chain.gz)
shopt -u nullglob
if (( ${#CHAIN_FILES[@]} == 0 )); then
    echo "Error: no hg38To*.over.chain.gz files found in $CHAIN_DIRECTORY" >&2
    exit 1
fi

PER_SPECIES_DIRECTORY="${OUTPUT_DIRECTORY}/per_species"
SUMMARY_DIRECTORY="${OUTPUT_DIRECTORY}/summary"
mkdir -p "$PER_SPECIES_DIRECTORY" "$SUMMARY_DIRECTORY"

TEMP_DIRECTORY="$(mktemp -d "${OUTPUT_DIRECTORY}/.liftover_tmp.XXXXXX")"
cleanup() {
    rm -rf -- "$TEMP_DIRECTORY"
}
trap cleanup EXIT

CLEAN_BED="${TEMP_DIRECTORY}/rute_input_clean.bed"
RUTE_IDS="${TEMP_DIRECTORY}/rute_ids.txt"
MATRIX_TEMP="${TEMP_DIRECTORY}/matrix.tsv"

# Keep the first five BED columns and ignore blank, comment, track, and browser lines.
awk 'BEGIN { OFS = "\t" }
     NF >= 5 && $1 !~ /^(#|track$|browser$)/ { print $1, $2, $3, $4, $5 }' \
    "$INPUT_BED" > "$CLEAN_BED"

if [[ ! -s "$CLEAN_BED" ]]; then
    echo "Error: no valid records were found in the input BED file." >&2
    exit 1
fi

cut -f4 "$CLEAN_BED" > "$RUTE_IDS"
if [[ -n "$(LC_ALL=C sort "$RUTE_IDS" | uniq -d | head -n 1)" ]]; then
    echo "Error: ruTE IDs in BED column 4 must be unique." >&2
    exit 1
fi
cp "$RUTE_IDS" "$MATRIX_TEMP"

species_from_chain() {
    local chain_name
    chain_name="$(basename -- "$1")"
    chain_name="${chain_name#hg38To}"
    printf '%s\n' "${chain_name%.over.chain.gz}"
}

process_chain() {
    local chain_file="$1"
    local species
    local lifted_bed
    local unmapped_bed
    local mapped_ids

    species="$(species_from_chain "$chain_file")"
    lifted_bed="${PER_SPECIES_DIRECTORY}/${species}.lifted.bed"
    unmapped_bed="${PER_SPECIES_DIRECTORY}/${species}.unmapped.bed"
    mapped_ids="${TEMP_DIRECTORY}/${species}.mapped_ids.txt"

    echo "Processing ${species}..."
    liftOver "$CLEAN_BED" "$chain_file" "$lifted_bed" "$unmapped_bed"
    cut -f4 "$lifted_bed" | LC_ALL=C sort -u > "$mapped_ids"
}

# Run independent species conversions concurrently in portable batches.
process_ids=()
for chain_file in "${CHAIN_FILES[@]}"; do
    process_chain "$chain_file" &
    process_ids+=("$!")

    if (( ${#process_ids[@]} >= MAX_PARALLEL )); then
        for process_id in "${process_ids[@]}"; do
            wait "$process_id"
        done
        process_ids=()
    fi
done
for process_id in "${process_ids[@]}"; do
    wait "$process_id"
done

COUNT_FILE="${SUMMARY_DIRECTORY}/liftover_counts.tsv"
FINAL_MATRIX="${SUMMARY_DIRECTORY}/rute_presence_matrix.tsv"
printf 'Species\tMapped_ruTE_count\n' > "$COUNT_FILE"

species_order=()
for chain_file in "${CHAIN_FILES[@]}"; do
    species="$(species_from_chain "$chain_file")"
    mapped_ids="${TEMP_DIRECTORY}/${species}.mapped_ids.txt"
    presence_column="${TEMP_DIRECTORY}/${species}.presence.txt"
    next_matrix="${TEMP_DIRECTORY}/matrix.next.tsv"

    awk 'NR == FNR { mapped[$1] = 1; next }
         { print ($1 in mapped) ? 1 : 0 }' \
        "$mapped_ids" \
        "$RUTE_IDS" > "$presence_column"

    paste "$MATRIX_TEMP" "$presence_column" > "$next_matrix"
    mv "$next_matrix" "$MATRIX_TEMP"

    species_order+=("$species")
    mapped_count="$(wc -l < "$mapped_ids")"
    printf '%s\t%s\n' "$species" "$mapped_count" >> "$COUNT_FILE"
done

{
    printf 'ruTE_ID'
    printf '\t%s' "${species_order[@]}"
    printf '\n'
    cat "$MATRIX_TEMP"
} > "$FINAL_MATRIX"

echo "Completed successfully."
echo "  Species processed: ${#CHAIN_FILES[@]}"
echo "  Presence matrix: $FINAL_MATRIX"
echo "  Mapping counts: $COUNT_FILE"
echo "  Per-species BED files: $PER_SPECIES_DIRECTORY"
