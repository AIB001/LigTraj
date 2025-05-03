#!/bin/bash

# Usage:
# bash pdb_clean.sh input_directory/ output_directory/

# Check arguments
if [[ $# -ne 2 ]]; then
    echo "Usage: bash pdb_clean.sh input_directory/ output_directory/"
    exit 1
fi

input_dir="${1%/}"    # Remove trailing slash
output_base="${2%/}"

# Check if input directory exists
if [[ ! -d "$input_dir" ]]; then
    echo "Error: Input directory '$input_dir' does not exist."
    exit 1
fi

# Create output base directory if it does not exist
mkdir -p "$output_base"

# Loop over matching PDB files
for pdb_file in "$input_dir"/*_protein_*.pdb; do
    if [[ -f "$pdb_file" ]]; then
        filename="$(basename "$pdb_file")"
        
        # Use pattern matching to extract yyy (pdbid)
        if [[ "$filename" =~ ^[^_]+_([^_]+)_protein_.*\.pdb$ ]]; then
            pdbid="${BASH_REMATCH[1]}"
            output_dir="${output_base}/${pdbid}"
            mkdir -p "$output_dir"
            output_pdb="${output_dir}/${pdbid}_clean.pdb"
            
            # Remove HETATM lines and write cleaned file
            grep -v '^HETATM' "$pdb_file" > "$output_pdb"
            echo "Processed: $filename --> $output_pdb"
        else
            echo "Warning: File '$filename' does not match expected pattern. Skipping."
        fi
    fi
done

echo "All matching files processed."
