# #!/bin/bash

# # Usage:
# # bash amber2gmx_v2.sh ligand.mol2 output_directory/
# SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
# echo "${SCRIPY_DIR}"
# #########################################
# # Check arguments
# #########################################
# if [[ $# -lt 2 ]]; then
#     echo "Usage: bash amber2gmx_v2.sh ligand.mol2 output_directory/"
#     exit 1
# fi

# input_mol2="$1"
# output_dir="${2%/}"   # Remove trailing slash

# if [[ ! -f "$input_mol2" ]]; then
#     echo "Error: mol2 file '$input_mol2' not found."
#     exit 1
# fi

# mkdir -p "$output_dir"

# ligand_name="$(basename "$input_mol2" .mol2)"

# amber_mol2="${output_dir}/${ligand_name}_amber.mol2"
# prep_file="${output_dir}/${ligand_name}.prep"
# frcmod_file="${output_dir}/${ligand_name}.frcmod"

# #########################################
# # Step 1: Generate Amber-compatible mol2 and prep
# #########################################
# echo "[Step 1] Generating mol2 and prep files..."

# antechamber -i "$input_mol2" -fi mol2 -o "$amber_mol2" -fo mol2 -c bcc -s 2 > "${output_dir}/antechamber_mol2.log" 2>&1
# if [[ ! -f "$amber_mol2" ]]; then
#     echo "Error: antechamber failed to generate Amber mol2."
#     exit 1
# fi

# antechamber -i "$amber_mol2" -fi mol2 -o "$prep_file" -fo prepi -c bcc -s 2 > "${output_dir}/antechamber_prep.log" 2>&1
# if [[ ! -f "$prep_file" ]]; then
#     echo "Error: antechamber failed to generate prep."
#     exit 1
# fi

# #########################################
# # Step 2: Generate frcmod
# #########################################
# echo "[Step 2] Generating frcmod..."

# parmchk2 -i "$prep_file" -f prepi -o "$frcmod_file" > "${output_dir}/parmchk2.log" 2>&1
# if [[ ! -f "$frcmod_file" ]]; then
#     echo "Error: parmchk2 failed to generate frcmod."
#     exit 1
# fi

# #########################################
# # Step 3: Prepare tleap input
# #########################################
# echo "[Step 3] Preparing tleap input..."

# tleap_input="${output_dir}/tleap_${ligand_name}.in"

# cat > "$tleap_input" << EOF
# source leaprc.protein.ff14SB
# source leaprc.gaff

# LIG = loadmol2 $amber_mol2
# loadamberparams $frcmod_file

# saveamberparm LIG ${output_dir}/${ligand_name}.prmtop ${output_dir}/${ligand_name}.rst7

# quit
# EOF

# #########################################
# # Step 4: Run tleap
# #########################################
# echo "[Step 4] Running tleap..."

# tleap -f "$tleap_input" > "${output_dir}/tleap.log" 2>&1

# if [[ ! -f "${output_dir}/${ligand_name}.prmtop" || ! -f "${output_dir}/${ligand_name}.rst7" ]]; then
#     echo "Error: tleap failed to generate prmtop/rst7."
#     exit 1
# fi

# #########################################
# # Step 5: Convert to GROMACS format with acpype
# #########################################
# echo "[Step 5] Running acpype..."

# cd "$output_dir"
# acpype -p "${ligand_name}.prmtop" -x "${ligand_name}.rst7" -d > acpype.log 2>&1

# # Find the actual ITP and GRO files
# itp_file=$(find . -maxdepth 1 -name "*_GMX.itp" | head -n 1)
# gro_file=$(find . -maxdepth 1 -name "*_GMX.gro" | head -n 1)
# top_file=$(find . -maxdepth 1 -name "*_GMX.top" | head -n 1)
# posre_file=$(find . -maxdepth 1 -name "posre_*.itp" | head -n 1)

# if [[ ! -f "$itp_file" || ! -f "$posre_file" ]]; then
#     echo "Hint: acpype may did not produce ITP or posre files, please check manually."
#     # exit 1
# fi

# cd "${SCRIPT_DIR}/../"

# # #########################################
# # # Step 6: Standardize output names
# # #########################################
# # echo "[Step 6] Renaming and organizing output..."

# # cp "$itp_file" LIG.itp
# # cp "$posre_file" posre_LIG.itp
# # cp "$gro_file" LIG_GMX.gro
# # cp "$top_file" LIG_GMX.top

# # #########################################
# # # Step 7: Ensure posre included in ITP
# # #########################################
# # echo "[Step 7] Checking posre include in LIG.itp..."

# # if ! grep -q "posre_LIG.itp" "LIG.itp"; then
# #     echo "" >> LIG.itp
# #     echo "#ifdef POSRES" >> LIG.itp
# #     echo "#include \"posre_LIG.itp\"" >> LIG.itp
# #     echo "#endif" >> LIG.itp
# # fi

# #########################################
# # Cleaning up temporary files...
# #########################################
# echo "[ Cleaning up... ]"
# pwd
# # rm -f sqm.* ANTECHAMBER* ATOMTYPE* leap.log \
# #     "${ligand_name}.acpype/" \
# #     *.ACPYPE* \
# #     *amb2gmx* \
# #     *md.mdp
# rm ANTECHAMBER_*
# rm sqm*
# rm leap*
# rm ATOMTYPE*
# rm NEWPDB*
# rm PREP*

# cd "${output_dir}"
# rm antechamber_*
# rm *frcmod
# rm *prep
# rm *prmtop
# rm *rst7
# rm *log
# rm *in
# cd "${SCRIPT_DIR}/../"

# echo "All done. Files are ready in $output_dir."
# echo "aib remind you: next work can be carried out now."

##########################################################
# Version 2.0
##########################################################

#!/bin/bash

# Usage:
# bash amber2gmx_batch.sh ligand_directory/ protein_output_directory/

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#########################################
# Check arguments
#########################################
if [[ $# -ne 2 ]]; then
    echo "Usage: bash amber2gmx_batch.sh ligand_directory/ protein_output_directory/"
    exit 1
fi

ligand_dir="${1%/}"
protein_base="${2%/}"

if [[ ! -d "$ligand_dir" ]]; then
    echo "Error: Ligand directory '$ligand_dir' does not exist."
    exit 1
fi

if [[ ! -d "$protein_base" ]]; then
    echo "Error: Protein directory '$protein_base' does not exist."
    exit 1
fi

#########################################
# Loop over mol2 files matching pattern
#########################################
for mol2_file in "$ligand_dir"/*_ligand_*.mol2; do
    if [[ -f "$mol2_file" ]]; then
        filename="$(basename "$mol2_file")"

        # Extract yyy (pdbid)
        if [[ "$filename" =~ ^[^_]+_([^_]+)_ligand_.*\.mol2$ ]]; then
            pdbid="${BASH_REMATCH[1]}"
            echo "[INFO] Processing ligand for PDB ID: $pdbid"
        else
            echo "[WARNING] File '$filename' does not match expected pattern. Skipping."
            continue
        fi

        #########################################
        # Prepare output directory
        #########################################
        target_dir="${protein_base}/${pdbid}"

        if [[ ! -d "$target_dir" ]]; then
            echo "[ERROR] Target protein directory '$target_dir' does not exist. Skipping $filename."
            continue
        fi

        forcefield_dir="${target_dir}/forcefield"
        mkdir -p "$forcefield_dir"

        ligand_name="${pdbid}_ligand"
        amber_mol2="${forcefield_dir}/${ligand_name}.mol2"
        prep_file="${forcefield_dir}/${ligand_name}.prep"
        frcmod_file="${forcefield_dir}/${ligand_name}.frcmod"

        #########################################
        # Step 1: Generate Amber-compatible mol2 and prep
        #########################################
        echo "[Step 1] Generating mol2 and prep files for $filename..."

        antechamber -i "$mol2_file" -fi mol2 -o "$amber_mol2" -fo mol2 -c bcc -s 2 > "${forcefield_dir}/antechamber_mol2.log" 2>&1
        if [[ ! -f "$amber_mol2" ]]; then
            echo "Error: antechamber failed to generate Amber mol2 for $filename."
            continue
        fi

        antechamber -i "$amber_mol2" -fi mol2 -o "$prep_file" -fo prepi -c bcc -s 2 > "${forcefield_dir}/antechamber_prep.log" 2>&1
        if [[ ! -f "$prep_file" ]]; then
            echo "Error: antechamber failed to generate prep for $filename."
            continue
        fi

        #########################################
        # Step 2: Generate frcmod
        #########################################
        echo "[Step 2] Generating frcmod for $filename..."

        parmchk2 -i "$prep_file" -f prepi -o "$frcmod_file" > "${forcefield_dir}/parmchk2.log" 2>&1
        if [[ ! -f "$frcmod_file" ]]; then
            echo "Error: parmchk2 failed to generate frcmod for $filename."
            continue
        fi

        #########################################
        # Step 3: Prepare tleap input
        #########################################
        tleap_input="${forcefield_dir}/tleap_${ligand_name}.in"

        cat > "$tleap_input" << EOF
source leaprc.protein.ff14SB
source leaprc.gaff

LIG = loadmol2 $amber_mol2
loadamberparams $frcmod_file

saveamberparm LIG ${forcefield_dir}/${ligand_name}.prmtop ${forcefield_dir}/${ligand_name}.rst7

quit
EOF

        #########################################
        # Step 4: Run tleap
        #########################################
        echo "[Step 4] Running tleap for $filename..."

        tleap -f "$tleap_input" > "${forcefield_dir}/tleap.log" 2>&1

        if [[ ! -f "${forcefield_dir}/${ligand_name}.prmtop" || ! -f "${forcefield_dir}/${ligand_name}.rst7" ]]; then
            echo "Error: tleap failed to generate prmtop/rst7 for $filename."
            continue
        fi

        #########################################
        # Step 5: Convert to GROMACS format with acpype
        #########################################
        echo "[Step 5] Running acpype for $filename..."

        cd "$forcefield_dir"
        acpype -p "${ligand_name}.prmtop" -x "${ligand_name}.rst7" -d > acpype.log 2>&1

        # Check acpype output
        itp_file=$(find . -maxdepth 1 -name "*_GMX.itp" | head -n 1)
        gro_file=$(find . -maxdepth 1 -name "*_GMX.gro" | head -n 1)

        if [[ ! -f "$itp_file" || ! -f "$gro_file" ]]; then
            echo "Warning: acpype did not produce ITP or GRO files for $filename. Please check manually."
        fi

        #########################################
        # Clean up intermediate files
        #########################################
        echo "[Cleanup] Removing intermediate files for $filename..."
        rm -f ANTECHAMBER_* sqm* leap* ATOMTYPE* NEWPDB* PREP*
        rm -f antechamber_*
        rm -f *frcmod *prep *prmtop *rst7 *log *in

        cd "$SCRIPT_DIR"

        echo "[DONE] Finished processing $filename"

    fi
done

echo "All ligands processed."

