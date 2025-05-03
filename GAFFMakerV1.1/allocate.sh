#!/bin/bash

# Function to copy iqb.md.slurm.sh into GMX_PROLIG_MD folders
copy_file_to_gmx_prolig_md() {
    local source_dir="${1%/}"  # Source directory containing iqb.md.slurm.sh
    folder_name=$(basename "${source_dir}")
    local seq_name="${2}"
    local time="${3}"
    cp config_files/batch_submit.sh "${source_dir}"
    cp config_files/iqb.md.slurm.sh "${source_dir}"
    cp config_files/localrun.sh "${source_dir}"
    # sed -i "s/SLURM_TASKNAME/${source_dir}/g" "${source_dir}/iqb.md.slurm.sh"
    sed -i "s|SLURM_TASKNAME|${folder_name}|g" "${source_dir}/iqb.md.slurm.sh"
    sed -i "s|SEQ_NAME|${seq_name}|g" "${source_dir}/iqb.md.slurm.sh"
    sed -i "s|TIME|${time}|g" "${source_dir}/iqb.md.slurm.sh"
    local source_file="$source_dir/iqb.md.slurm.sh"  # Path to the file to copy

    # Check if the source file exists
    if [ ! -f "$source_file" ]; then
        echo "Error: File '$source_file' does not exist."
        exit 1
    fi

    # Find all GMX_PROLIG_MD folders under the source directory
    find "$source_dir" -type d -name "GMX_PROLIG_MD" | while read -r gmx_folder; do
        echo "Copying file to: $gmx_folder"
        cp "$source_file" "$gmx_folder"  # Copy the file to the GMX_PROLIG_MD folder
    done
}

# Main script
if [ $# -ne 3 ]; then
    echo "Usage: $0 <source_directory> <serquence_name> <time>"  # Print usage if no directory is provided
    exit 1
fi

source_directory="$1"  # Get the source directory from the command line argument
seq="${2}"
time="${3}"

# Check if the source directory exists
if [ -d "$source_directory" ]; then
    copy_file_to_gmx_prolig_md "$source_directory" "${seq}" "${time}"
    echo "Operation completed."
else
    echo "Error: Directory '$source_directory' does not exist."  # Print error if directory is invalid
    exit 1
fi