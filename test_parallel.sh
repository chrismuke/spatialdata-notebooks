#!/bin/bash

# Base directories
ZARR_DIR="/media/chrism/Share1TB/10x_xenium_mouse_datasets/zarr_xenosplit"
REF_DIR="/media/chrism/Share1TB/10x_xenium_mouse_datasets/xenium_ref"
RESULTS_DIR="/media/chrism/Share1TB/10x_xenium_mouse_datasets/results"

# Arrays of files (just test with first ones)
ZARR_FILES=(
    "Xenium_V1_mouse_pup.zarr"
)

REF_FILES=(
    "0c44ac00-d0a8-4a5d-93e6-da8df2488ad2.h5ad"
    "2430dcf3-40b8-43dc-8f8c-87856c035db1.h5ad"
)

# Function to run a single annotation
run_annotation() {
    local zarr_file="$1"
    local ref_file="$2"  
    local job_id="$3"

    echo "[Job $job_id] Starting: $zarr_file + $ref_file at $(date)"
    echo "[Job $job_id] Paths: $ZARR_DIR/$zarr_file + $REF_DIR/$ref_file"
    
    # Simulated work
    sleep 2
    echo "[Job $job_id] ✓ Completed: $zarr_file + $ref_file at $(date)"
    
    return 0
}

export -f run_annotation
export ZARR_DIR REF_DIR RESULTS_DIR

# Create combinations with bash job control
job_id=1
active_jobs=0
MAX_JOBS=2

for zarr_file in "${ZARR_FILES[@]}"; do
    for ref_file in "${REF_FILES[@]}"; do
        # Wait if we've reached max jobs
        while [ $active_jobs -ge $MAX_JOBS ]; do
            wait -n  # Wait for any job to complete
            active_jobs=$((active_jobs - 1))
        done

        # Start new job in background
        run_annotation "$zarr_file" "$ref_file" "$job_id" &
        active_jobs=$((active_jobs + 1))
        
        echo "Started job $job_id (active: $active_jobs)"
        job_id=$((job_id + 1))
    done
done

# Wait for all remaining jobs to complete
wait

echo "All jobs completed!"