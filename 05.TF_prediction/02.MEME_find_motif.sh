#!/bin/bash
JASPAR_DB="/data2t_2/pathogen_TE_2025_New/08.TF_prediction/combine_JASPAR2024_CORE_vertebrates_ZNF_Barazandeh2017.txt"
INPUT_DIR="TE_sequences"
OUTPUT_DIR="motif_analysis_results"
P_VALUE_THRESHOLD="1e-4"
THREADS=30

mkdir -p $OUTPUT_DIR

echo "Starting FIMO motif analysis with $THREADS threads..."
echo "Using JASPAR database: $JASPAR_DB"
echo "P-value threshold: $P_VALUE_THRESHOLD"
echo "==========================================="

process_single_file() {
    local fasta_file=$1
    local basename=$(basename "$fasta_file" .fa)
    local output_subdir="${OUTPUT_DIR}/${basename}"
    
    echo "Processing: $basename (PID: $$)"
    mkdir -p "$output_subdir"
    
    fimo --oc "$output_subdir" \
         --thresh "$P_VALUE_THRESHOLD" \
         --parse-genomic-coord \
         "$JASPAR_DB" \
         "$fasta_file" 2>/dev/null
    
    if [ $? -eq 0 ]; then
        motif_count=$(tail -n +2 "$output_subdir/fimo.tsv" | wc -l)
        unique_motifs=$(tail -n +2 "$output_subdir/fimo.tsv" | cut -f1 | sort | uniq | wc -l)
        
        echo "✓ Completed: $basename - Found $motif_count motif instances, covering $unique_motifs unique motifs"
        return 0
    else
        echo "✗ Failed: $basename"
        return 1
    fi
}

export -f process_single_file
export JASPAR_DB OUTPUT_DIR P_VALUE_THRESHOLD

find "$INPUT_DIR" -name "*.fa" | parallel -j $THREADS process_single_file {}
# echo -e "$INPUT_DIR/LTR43_loci.fa\n$INPUT_DIR/MIR_loci.fa\n$INPUT_DIR/MLT1A1_loci.fa\n$INPUT_DIR/MLT1H_loci.fa\n$INPUT_DIR/MLT1K_loci.fa" | parallel -j $THREADS process_single_file {}
echo "==========================================="
echo "FIMO analysis completed!"
echo "Results saved in: $OUTPUT_DIR"
