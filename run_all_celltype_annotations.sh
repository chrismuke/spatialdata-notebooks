#!/bin/bash

# Comprehensive cell type annotation script for all liver dataset combinations
# This script runs celltype_annotate_cli_v2.py for all combinations of liver datasets and reference datasets

# Set base directories
LIVER_DIR="/home/chrism/datasets/lv_spatialdat_liu"
REF_DIR="/home/chrism/datasets/scRNAseq_skin"
RESULTS_DIR="/home/chrism/datasets/lv_spatialdat_liu/results"

# Create results directory if it doesn't exist
mkdir -p "$RESULTS_DIR"

# Define liver datasets (zarr files)
LIVER_DATASETS=(
    "lv_0046706_007.zarr"
    "lv_0046706_053.zarr"
    "lv_0046706_117.zarr"
    "lv_0046706_71_control.zarr"
    "lv_0046706_75_control.zarr"
    "lv_0046859_002.zarr"
    "lv_0046859_033.zarr"
    "lv_0046859_039.zarr"
    "lv_0046859_736_control.zarr"
)

# Define reference datasets
REF_DATASETS=(
    "d0c12af4-c0e4-4c7b-873a-70752b449689/original_adata.h5ad"
    "f512b8b6-369d-4a85-a695-116e0806857f/original_adata.h5ad"
)

# Log file for tracking progress
LOG_FILE="$RESULTS_DIR/celltype_annotation_batch.log"
echo "Starting batch cell type annotation at $(date)" > "$LOG_FILE"

# Function to run cell type annotation
run_annotation() {
    local liver_dataset="$1"
    local ref_dataset="$2"
    local liver_path="$LIVER_DIR/$liver_dataset"
    local ref_path="$REF_DIR/$ref_dataset"
    
    echo "Processing: $liver_dataset with reference $ref_dataset" | tee -a "$LOG_FILE"
    echo "Started at: $(date)" | tee -a "$LOG_FILE"
    
    # Run the cell type annotation with zarr file creation
    uv run python celltype_annotate_cli_v2.py "$liver_path" "$ref_path" \
        --results-dir "$RESULTS_DIR" \
        --min-clusters 5 \
        --max-clusters 15 \
        2>&1 | tee -a "$LOG_FILE"
    
    local exit_code=${PIPESTATUS[0]}
    if [ $exit_code -eq 0 ]; then
        echo "✓ Successfully completed: $liver_dataset with $ref_dataset at $(date)" | tee -a "$LOG_FILE"
    else
        echo "✗ Failed: $liver_dataset with $ref_dataset at $(date) (exit code: $exit_code)" | tee -a "$LOG_FILE"
    fi
    echo "----------------------------------------" | tee -a "$LOG_FILE"
}

# Main execution loop
total_combinations=$((${#LIVER_DATASETS[@]} * ${#REF_DATASETS[@]}))
current_combination=0

echo "Total combinations to process: $total_combinations" | tee -a "$LOG_FILE"
echo "========================================" | tee -a "$LOG_FILE"

for liver_dataset in "${LIVER_DATASETS[@]}"; do
    for ref_dataset in "${REF_DATASETS[@]}"; do
        current_combination=$((current_combination + 1))
        echo "Processing combination $current_combination of $total_combinations" | tee -a "$LOG_FILE"
        
        # Check if liver dataset exists
        if [ ! -d "$LIVER_DIR/$liver_dataset" ]; then
            echo "Warning: Liver dataset not found: $LIVER_DIR/$liver_dataset" | tee -a "$LOG_FILE"
            continue
        fi
        
        # Check if reference dataset exists
        if [ ! -f "$REF_DIR/$ref_dataset" ]; then
            echo "Warning: Reference dataset not found: $REF_DIR/$ref_dataset" | tee -a "$LOG_FILE"
            continue
        fi
        
        # Run the annotation
        run_annotation "$liver_dataset" "$ref_dataset"
        
        # Add a small delay to avoid overwhelming the system
        sleep 2
    done
done

echo "Batch processing completed at $(date)" | tee -a "$LOG_FILE"
echo "Results stored in: $RESULTS_DIR" | tee -a "$LOG_FILE"

# Generate summary
echo "========================================" | tee -a "$LOG_FILE"
echo "SUMMARY:" | tee -a "$LOG_FILE"
successful_runs=$(grep -c "✓ Successfully completed" "$LOG_FILE")
failed_runs=$(grep -c "✗ Failed" "$LOG_FILE")
echo "Successful runs: $successful_runs" | tee -a "$LOG_FILE"
echo "Failed runs: $failed_runs" | tee -a "$LOG_FILE"
echo "Total combinations processed: $((successful_runs + failed_runs))" | tee -a "$LOG_FILE"