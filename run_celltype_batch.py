#!/usr/bin/env python3

import subprocess
import sys
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

# Configuration
ZARR_DIR = "/media/chrism/Share1TB/10x_xenium_mouse_datasets/zarr_xenosplit"
REF_DIR = "/media/chrism/Share1TB/10x_xenium_mouse_datasets/xenium_ref"
RESULTS_DIR = "/media/chrism/Share1TB/10x_xenium_mouse_datasets/results"
MAX_WORKERS = 8

# File lists
ZARR_FILES = [
    "Xenium_Prime_Mouse_Brain_Coronal_FF.zarr",
    "Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP.zarr", 
    "Xenium_V1_FF_Mouse_Brain_Coronal.zarr",
    "Xenium_V1_FF_Mouse_Brain_MultiSection_1.zarr",
    "Xenium_V1_FF_Mouse_Brain_MultiSection_2.zarr",
    "Xenium_V1_FF_Mouse_Brain_MultiSection_3.zarr",
    "Xenium_V1_FFPE_TgCRND8_17_9_months.zarr",
    "Xenium_V1_FFPE_TgCRND8_2_5_months.zarr",
    "Xenium_V1_FFPE_TgCRND8_5_7_months.zarr",
    "Xenium_V1_FFPE_wildtype_13_4_months.zarr",
    "Xenium_V1_FFPE_wildtype_2_5_months.zarr",
    "Xenium_V1_FFPE_wildtype_5_7_months.zarr",
    "Xenium_V1_mFemur_EDTA_3daydecal_section.zarr",
    "Xenium_V1_mFemur_EDTA_PFA_3daydecal_section.zarr",
    "Xenium_V1_mFemur_formic_acid_24hrdecal_section.zarr",
    "Xenium_V1_mouse_Colon_FF.zarr"
]

REF_FILES = [
    "0c44ac00-d0a8-4a5d-93e6-da8df2488ad2.h5ad",
    "2430dcf3-40b8-43dc-8f8c-87856c035db1.h5ad",
    "24b8aa1f-258c-4d47-8f96-e30fa8d4aeb0.h5ad",
    "40320ee1-a655-4cbd-bbca-7608d5fc4eb6.h5ad",
    "69f6ba65-d094-4ee8-96d5-b2754ff1c084.h5ad",
    "756ea4fa-e736-4533-878f-6e6d32339da3.h5ad",
    "fb1eabf9-1fae-4253-97c1-52ecdee0f2df.h5ad",
    "fe173071-fdad-456a-8fb9-46442556b6b9.h5ad"
]

def run_annotation(zarr_file, ref_file, job_id):
    """Run cell type annotation for a single combination"""
    zarr_path = os.path.join(ZARR_DIR, zarr_file)
    ref_path = os.path.join(REF_DIR, ref_file)
    
    print(f"[Job {job_id}] Starting: {zarr_file} + {ref_file}")
    
    try:
        cmd = [
            "uv", "run", "python", "celltype_annotate_cli_v2.py",
            zarr_path, ref_path,
            "--results-dir", RESULTS_DIR,
            "--min-clusters", "5",
            "--max-clusters", "10"
        ]
        
        result = subprocess.run(cmd, capture_output=True, text=True, cwd=os.getcwd())
        
        if result.returncode == 0:
            print(f"[Job {job_id}] ✓ Completed: {zarr_file} + {ref_file}")
            return True, job_id, zarr_file, ref_file, result.stdout
        else:
            print(f"[Job {job_id}] ✗ Failed: {zarr_file} + {ref_file}")
            print(f"[Job {job_id}] Error: {result.stderr}")
            return False, job_id, zarr_file, ref_file, result.stderr
            
    except Exception as e:
        print(f"[Job {job_id}] Exception: {str(e)}")
        return False, job_id, zarr_file, ref_file, str(e)

def main():
    # Create all combinations
    combinations = []
    job_id = 1
    for zarr_file in ZARR_FILES:
        for ref_file in REF_FILES:
            combinations.append((zarr_file, ref_file, job_id))
            job_id += 1
    
    total_combinations = len(combinations)
    print(f"Starting parallel cell type annotation for {total_combinations} combinations...")
    print(f"Using {MAX_WORKERS} parallel workers")
    print(f"Results will be stored in: {RESULTS_DIR}")
    print()
    
    # Track results
    completed = 0
    failed = 0
    
    # Run in parallel
    with ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futures = {executor.submit(run_annotation, zarr, ref, jid): (zarr, ref, jid) 
                  for zarr, ref, jid in combinations}
        
        for future in as_completed(futures):
            success, job_id, zarr_file, ref_file, output = future.result()
            if success:
                completed += 1
            else:
                failed += 1
            
            print(f"Progress: {completed + failed}/{total_combinations} completed ({completed} success, {failed} failed)")
    
    print()
    print("All cell type annotations completed!")
    print(f"Results: {completed} successful, {failed} failed")
    print(f"Results stored in: {RESULTS_DIR}")

if __name__ == "__main__":
    main()