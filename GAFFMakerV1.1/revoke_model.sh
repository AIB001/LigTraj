#!/bin/bash

# Function to delete GMX_PROLIG_MD folders
delete_gmx_prolig_md() {
    local dir="$1"  # Directory to search in

    # Find all directories named GMX_PROLIG_MD under the given directory
    find "$dir" -type d -name "GMX_PROLIG_MD" | while read -r folder; do
        echo "Deleting folder: $folder"
        rm -rf "$folder"  # Recursively delete the folder
    done
}

# Main script
if [ $# -ne 1 ]; then
    echo "Usage: $0 <directory>"  # Print usage if no directory is provided
    exit 1
fi

root_directory="$1"  # Get the directory from the command line argument

# Check if the directory exists
if [ -d "$root_directory" ]; then
    delete_gmx_prolig_md "$root_directory"
    echo "Operation completed."
else
    echo "Error: Directory '$root_directory' does not exist."  # Print error if directory is invalid
    exit 1
fi