#\!/bin/bash

# Base directories
ZARR_DIR="/media/chrism/Share1TB/10x_xenium_mouse_datasets/zarr_xenosplit"
REF_DIR="/media/chrism/Share1TB/10x_xenium_mouse_datasets/xenium_ref"
RESULTS_DIR="/media/chrism/Share1TB/10x_xenium_mouse_datasets/results"

# Maximum number of parallel jobs
MAX_JOBS=8

# Arrays of files
ZARR_FILES=(
    "Xenium_Prime_Mouse_Brain_Coronal_FF.zarr"
    "Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP.zarr"
    "Xenium_V1_FF_Mouse_Brain_Coronal.zarr"
    "Xenium_V1_FF_Mouse_Brain_MultiSection_1.zarr"
    "Xenium_V1_FF_Mouse_Brain_MultiSection_2.zarr"
    "Xenium_V1_FF_Mouse_Brain_MultiSection_3.zarr"
    "Xenium_V1_FFPE_TgCRND8_17_9_months.zarr"
    "Xenium_V1_FFPE_TgCRND8_2_5_months.zarr"
    "Xenium_V1_FFPE_TgCRND8_5_7_months.zarr"
    "Xenium_V1_FFPE_wildtype_13_4_months.zarr"
    "Xenium_V1_FFPE_wildtype_2_5_months.zarr"
    "Xenium_V1_FFPE_wildtype_5_7_months.zarr"
    "Xenium_V1_mFemur_EDTA_3daydecal_section.zarr"
    "Xenium_V1_mFemur_EDTA_PFA_3daydecal_section.zarr"
    "Xenium_V1_mFemur_formic_acid_24hrdecal_section.zarr"
    "Xenium_V1_mouse_Colon_FF.zarr"
)

REF_FILES=(
    "0c44ac00-d0a8-4a5d-93e6-da8df2488ad2.h5ad"
    "2430dcf3-40b8-43dc-8f8c-87856c035db1.h5ad"
    "24b8aa1f-258c-4d47-8f96-e30fa8d4aeb0.h5ad"
    "40320ee1-a655-4cbd-bbca-7608d5fc4eb6.h5ad"
    "69f6ba65-d094-4ee8-96d5-b2754ff1c084.h5ad"
    "756ea4fa-e736-4533-878f-6e6d32339da3.h5ad"
    "fb1eabf9-1fae-4253-97c1-52ecdee0f2df.h5ad"
    "fe173071-fdad-456a-8fb9-46442556b6b9.h5ad"
)

# Function to run a single annotation
run_annotation() {
    local zarr_file="$1"
    local ref_file="$2"
    local job_id="$3"
    
    echo "[Job $job_id] Starting: $zarr_file + $ref_file at $(date)"
    
    uv run python celltype_annotate_cli_v2.py \
        "$ZARR_DIR/$zarr_file" \
        "$REF_DIR/$ref_file" \
        --results-dir "$RESULTS_DIR" \
        --min-clusters 5 \
        --max-clusters 10 \
        2>&1 | sed "s/^/[Job $job_id] /"
    
    local exit_code=$?
    if [ $exit_code -eq 0 ]; then
        echo "[Job $job_id] ✓ Completed: $zarr_file + $ref_file at $(date)"
    else
        echo "[Job $job_id] ✗ Failed: $zarr_file + $ref_file at $(date)"
    fi
    
    return $exit_code
}

# Export function for parallel execution
export -f run_annotation
export ZARR_DIR REF_DIR RESULTS_DIR

# Create list of all combinations
combinations=()
job_id=1
for zarr_file in "${ZARR_FILES[@]}"; do
    for ref_file in "${REF_FILES[@]}"; do
        combinations+=("$zarr_file|$ref_file|$job_id")
        job_id=$((job_id + 1))
    done
done

total_combinations=${#combinations[@]}
echo "Starting parallel cell type annotation for $total_combinations combinations..."
echo "Using $MAX_JOBS parallel jobs"
echo "Results will be stored in: $RESULTS_DIR"
echo ""

# Check if GNU parallel is available
if command -v parallel &> /dev/null; then
    echo "Using GNU parallel for job management"
    printf '%s\n' "${combinations[@]}" | parallel -j $MAX_JOBS --colsep '|' run_annotation {1} {2} {3}
else
    echo "GNU parallel not found, using bash job control"
    # Fallback: use bash job control
    active_jobs=0
    
    for combination in "${combinations[@]}"; do
        IFS='|' read -r zarr_file ref_file job_id <<< "$combination"
        
        # Wait if we've reached max jobs
        while [ $active_jobs -ge $MAX_JOBS ]; do
            wait -n  # Wait for any job to complete
            active_jobs=$((active_jobs - 1))
        done
        
        # Start new job in background
        run_annotation "$zarr_file" "$ref_file" "$job_id" &
        active_jobs=$((active_jobs + 1))
        
        echo "Started job $job_id/$total_combinations (active: $active_jobs)"
    done
    
    # Wait for all remaining jobs to complete
    wait
fi

echo ""
echo "All cell type annotations completed\!"
echo "Results stored in: $RESULTS_DIR"
EOF < /dev/null