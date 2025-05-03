# #!/bin/bash

# #########################################
# # Check if a directory argument is provided
# #########################################
# if [ $# -ne 1 ]; then
#     echo "Usage: $0 <target_directory>"
#     exit 1
# fi

# target_dir="$1"

# # Verify the directory exists
# if [ ! -d "$target_dir" ]; then
#     echo "Error: Directory '$target_dir' does not exist"
#     exit 1
# fi

# #########################################
# # Change to the target directory
# #########################################
# cd "$target_dir" || { echo "Error: Failed to enter directory '$target_dir'"; exit 1; }

# #########################################
# # Detect .amb2gmx directory
# #########################################
# amb2gmx_dirs=( *.amb2gmx )
# if [ ${#amb2gmx_dirs[@]} -eq 0 ]; then
#     echo "Error: No .amb2gmx directory found in '$target_dir'"
#     exit 1
# elif [ ${#amb2gmx_dirs[@]} -gt 1 ]; then
#     echo "Error: Multiple .amb2gmx directories found in '$target_dir'"
#     exit 1
# fi

# original_dir="${amb2gmx_dirs[0]}"
# ligname="${original_dir%.amb2gmx}"

# echo "Detected ligand name: $ligname"

# #########################################
# # If ligand name is already LIG, skip all operations
# #########################################
# if [ "$ligname" == "LIG" ]; then
#     echo "The ligand is already standardized as 'LIG'. Skipping all operations."
#     exit 0
# fi

# #########################################
# # Step 1: Rename the main directory
# #########################################
# echo "Renaming directory: $original_dir -> LIG.amb2gmx"
# mv -v "$original_dir" "LIG.amb2gmx" || { echo "Directory rename failed"; exit 1; }

# cd "LIG.amb2gmx" || { echo "Failed to enter LIG.amb2gmx directory"; exit 1; }

# #########################################
# # Step 2: Rename files
# #########################################
# declare -A file_map=(
#     ["${ligname}_GMX.gro"]="LIG.gro"
#     ["${ligname}_GMX.top"]="LIG.top"
#     ["posre_${ligname}.itp"]="posre_LIG.itp"
# )

# echo "Renaming files:"
# for old_name in "${!file_map[@]}"; do
#     new_name="${file_map[$old_name]}"
#     if [ ! -f "$old_name" ]; then
#         echo "Error: Required file '$old_name' not found"
#         exit 1
#     fi
#     mv -v "$old_name" "$new_name" || { echo "File rename failed"; exit 1; }
# done

# #########################################
# # Step 3: Replace ligand name with LIG in all files
# #########################################

# echo -e "\nUpdating file contents: replacing '$ligname' with 'LIG'"

# # Detect if GNU or BSD sed is available
# if sed --version >/dev/null 2>&1; then
#     sed_command="GNU"
# else
#     sed_command="BSD"
# fi

# # Loop over all regular files and perform in-place replacement
# for file in *; do
#     if [ -f "$file" ]; then
#         echo "Processing $file"
#         if [ "$sed_command" == "GNU" ]; then
#             sed -i.bak "s/$ligname/LIG/g" "$file" && rm -f "${file}.bak"
#         else
#             sed -i "" "s/$ligname/LIG/g" "$file"
#         fi
#     fi
# done

# #########################################
# echo -e "\nOperation completed successfully!"
# echo "All occurrences of '$ligname' have been replaced with 'LIG'"
# #########################################


###########################################
# Version 2.0
###########################################

# #!/bin/bash

# #########################################
# # Usage check
# #########################################
# if [ $# -ne 1 ]; then
#     echo "Usage: $0 <parent_directory>"
#     exit 1
# fi

# parent_dir="$1"

# if [ ! -d "$parent_dir" ]; then
#     echo "Error: Directory '$parent_dir' does not exist"
#     exit 1
# fi

# #########################################
# # Prepare missing amb2gmx list
# #########################################
# missing_amb2gmx=()

# #########################################
# # Loop through all subdirectories
# #########################################
# for folder in "$parent_dir"/*/; do
#     [ -d "$folder" ] || continue  # skip non-directories

#     echo "Processing folder: $folder"

#     #########################################
#     # Step 1: Check and move forcefield content
#     #########################################
#     forcefield_dir="${folder}forcefield"

#     if [ -d "$forcefield_dir" ]; then
#         echo "  Found forcefield directory. Moving contents up..."

#         shopt -s dotglob  # include hidden files
#         mv -v "$forcefield_dir"/* "$folder" 2>/dev/null
#         shopt -u dotglob

#         echo "  forcefield content moved, but directory kept."
#     else
#         echo "  No forcefield directory found."
#     fi

#     #########################################
#     # Step 2: Check for .amb2gmx directory
#     #########################################
#     amb2gmx_dirs=( "$folder"/*.amb2gmx )
#     if [ ${#amb2gmx_dirs[@]} -eq 0 ] || [ ! -d "${amb2gmx_dirs[0]}" ]; then
#         echo "  No .amb2gmx directory found."
#         missing_amb2gmx+=("$(basename "$folder")")
#         continue
#     elif [ ${#amb2gmx_dirs[@]} -gt 1 ]; then
#         echo "  Warning: Multiple .amb2gmx directories found in '$folder'. Skipping."
#         missing_amb2gmx+=("$(basename "$folder") (multiple amb2gmx)")
#         continue
#     fi

#     original_dir="$(basename "${amb2gmx_dirs[0]}")"
#     ligname="${original_dir%.amb2gmx}"

#     echo "  Detected ligand name: $ligname"

#     #########################################
#     # Step 3: If ligand already LIG, skip
#     #########################################
#     if [ "$ligname" == "LIG" ]; then
#         echo "  Ligand already standardized as 'LIG'. Skipping."
#         continue
#     fi

#     #########################################
#     # Step 4: Rename directory to LIG.amb2gmx
#     #########################################
#     echo "  Renaming directory: $original_dir -> LIG.amb2gmx"
#     mv -v "$folder/$original_dir" "$folder/LIG.amb2gmx" || { echo "  Directory rename failed"; continue; }

#     work_dir="$folder/LIG.amb2gmx"

#     #########################################
#     # Step 5: Rename files inside LIG.amb2gmx
#     #########################################
#     declare -A file_map=(
#         ["${ligname}_GMX.gro"]="LIG.gro"
#         ["${ligname}_GMX.top"]="LIG.top"
#         ["posre_${ligname}.itp"]="posre_LIG.itp"
#     )

#     echo "  Renaming files:"
#     for old_name in "${!file_map[@]}"; do
#         new_name="${file_map[$old_name]}"
#         if [ ! -f "$work_dir/$old_name" ]; then
#             echo "    Error: Required file '$old_name' not found in $work_dir"
#             continue 2  # Skip to next folder
#         fi
#         mv -v "$work_dir/$old_name" "$work_dir/$new_name" || { echo "    File rename failed"; continue 2; }
#     done

#     #########################################
#     # Step 6: Replace ligand name with LIG in all files
#     #########################################
#     echo "  Updating file contents: replacing '$ligname' with 'LIG'"

#     # Detect sed type
#     if sed --version >/dev/null 2>&1; then
#         sed_command="GNU"
#     else
#         sed_command="BSD"
#     fi

#     for file in "$work_dir"/*; do
#         [ -f "$file" ] || continue
#         echo "    Processing $file"
#         if [ "$sed_command" == "GNU" ]; then
#             sed -i.bak "s/$ligname/LIG/g" "$file" && rm -f "${file}.bak"
#         else
#             sed -i "" "s/$ligname/LIG/g" "$file"
#         fi
#     done

#     echo "  Completed processing $folder"

# done

# #########################################
# # Summary of missing amb2gmx directories
# #########################################
# echo -e "\n=============================="
# echo "Folders missing .amb2gmx directories:"
# if [ ${#missing_amb2gmx[@]} -eq 0 ]; then
#     echo "None! All folders processed successfully."
# else
#     for miss in "${missing_amb2gmx[@]}"; do
#         echo " - $miss"
#     done
# fi
# echo "=============================="

##########################################
# Version 2.1
##########################################

#!/bin/bash

#########################################
# Usage check
#########################################
if [ $# -ne 1 ]; then
    echo "Usage: $0 <parent_directory>"
    exit 1
fi

parent_dir="$1"

if [ ! -d "$parent_dir" ]; then
    echo "Error: Directory '$parent_dir' does not exist"
    exit 1
fi

#########################################
# Prepare missing amb2gmx list
#########################################
missing_amb2gmx=()

#########################################
# Loop through all subdirectories
#########################################
for folder in "$parent_dir"/*/; do
    [ -d "$folder" ] || continue  # skip non-directories

    echo "Processing folder: $folder"

    #########################################
    # Step 1: Check and move forcefield content
    #########################################
    forcefield_dir="${folder}forcefield"

    if [ -d "$forcefield_dir" ]; then
        echo "  Found forcefield directory. Moving contents up..."

        shopt -s dotglob  # include hidden files
        mv -v "$forcefield_dir"/* "$folder" 2>/dev/null
        shopt -u dotglob

        echo "  forcefield content moved, but directory kept."
    else
        echo "  No forcefield directory found."
    fi

    #########################################
    # Step 2: Check for .amb2gmx directory
    #########################################
    amb2gmx_dirs=( "$folder"/*.amb2gmx )
    if [ ${#amb2gmx_dirs[@]} -eq 0 ] || [ ! -d "${amb2gmx_dirs[0]}" ]; then
        echo "  No .amb2gmx directory found."
        missing_amb2gmx+=("$(basename "$folder")")
        continue
    elif [ ${#amb2gmx_dirs[@]} -gt 1 ]; then
        echo "  Warning: Multiple .amb2gmx directories found in '$folder'. Skipping."
        missing_amb2gmx+=("$(basename "$folder") (multiple amb2gmx)")
        continue
    fi

    original_dir="$(basename "${amb2gmx_dirs[0]}")"
    ligname="${original_dir%.amb2gmx}"

    echo "  Detected ligand name: $ligname"

    #########################################
    # Step 3: If ligand already LIG, skip
    #########################################
    if [ "$ligname" == "LIG" ]; then
        echo "  Ligand already standardized as 'LIG'. Skipping."
        continue
    fi

    #########################################
    # Step 4: Rename directory to LIG.amb2gmx
    #########################################
    echo "  Renaming directory: $original_dir -> LIG.amb2gmx"
    mv -v "$folder/$original_dir" "$folder/LIG.amb2gmx" || { echo "  Directory rename failed"; continue; }

    work_dir="$folder/LIG.amb2gmx"

    #########################################
    # Step 5: Rename files inside LIG.amb2gmx
    #########################################
    declare -A file_map=(
        ["${ligname}_GMX.gro"]="LIG.gro"
        ["${ligname}_GMX.top"]="LIG.top"
        ["posre_${ligname}.itp"]="posre_LIG.itp"
    )

    echo "  Renaming files:"
    for old_name in "${!file_map[@]}"; do
        new_name="${file_map[$old_name]}"
        if [ ! -f "$work_dir/$old_name" ]; then
            echo "    Error: Required file '$old_name' not found in $work_dir"
            missing_amb2gmx+=("$(basename "$folder") (missing $old_name)")
            continue 2  # Skip to next folder
        fi
        mv -v "$work_dir/$old_name" "$work_dir/$new_name" || { echo "    File rename failed"; missing_amb2gmx+=("$(basename "$folder") (rename failed)"); continue 2; }
    done

    #########################################
    # Step 6: Replace ligand name with LIG in all files
    #########################################
    echo "  Updating file contents: replacing '$ligname' with 'LIG'"

    # Detect sed type
    if sed --version >/dev/null 2>&1; then
        sed_command="GNU"
    else
        sed_command="BSD"
    fi

    for file in "$work_dir"/*; do
        [ -f "$file" ] || continue
        echo "    Processing $file"
        if [ "$sed_command" == "GNU" ]; then
            sed -i.bak "s/$ligname/LIG/g" "$file" && rm -f "${file}.bak"
        else
            sed -i "" "s/$ligname/LIG/g" "$file"
        fi
    done

    echo "  Completed processing $folder"

done

#########################################
# Summary of missing amb2gmx directories
#########################################
echo -e "\n=============================="
echo "Folders missing .amb2gmx or with errors:"
if [ ${#missing_amb2gmx[@]} -eq 0 ]; then
    echo "None! All folders processed successfully."
else
    for miss in "${missing_amb2gmx[@]}"; do
        echo " - $miss"
    done
fi
echo "=============================="

#########################################
# Move abnormal folders to ../abnormal
#########################################

if [ ${#missing_amb2gmx[@]} -gt 0 ]; then
    abnormal_dir="$(dirname "$parent_dir")/abnormal"
    mkdir -p "$abnormal_dir"

    echo -e "\nMoving abnormal folders to $abnormal_dir"

    for miss in "${missing_amb2gmx[@]}"; do
        # 取掉可能的 (xxx) 注释
        clean_name="${miss%% *}"

        src="$parent_dir/$clean_name"
        if [ -d "$src" ]; then
            mv -v "$src" "$abnormal_dir/"
        else
            echo "  Warning: $src does not exist or already moved."
        fi
    done
else
    echo "No abnormal folders to move."
fi

echo -e "\nAll operations completed."

