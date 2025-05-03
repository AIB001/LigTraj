#!/bin/bash

# Save the current directory (starting directory)
start_dir=$(pwd)

# Find all subdirectories in the current directory
find . -type d | while read -r dir; do
    # Skip the current directory (.) itself
    if [ "$dir" != "." ]; then
        # Check if the directory contains a GMX_PROLIG_MD folder
        if [ -d "$dir/GMX_PROLIG_MD" ]; then
            echo "Entering directory: $dir/GMX_PROLIG_MD"
            cd "$dir/GMX_PROLIG_MD" || { echo "Failed to enter $dir/GMX_PROLIG_MD"; continue; }

            # Check if iqb.md.slurm.sh exists in the GMX_PROLIG_MD folder
            if [ -f "iqb.md.slurm.sh" ]; then
                echo "Submitting job: iqb.md.slurm.sh in $dir/GMX_PROLIG_MD"
                sbatch iqb.md.slurm.sh  # Submit the job
            else
                echo "Error: iqb.md.slurm.sh not found in $dir/GMX_PROLIG_MD"
            fi

            # Return to the starting directory
            cd "$start_dir" || { echo "Failed to return to the starting directory"; exit 1; }
            echo $(pwd)
        fi
    fi
done

echo "All jobs submitted successfully."