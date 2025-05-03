# #!/bin/bash

# # ============================================
# # Script name: prepare_system.sh
# # Purpose: Automate protein-ligand system setup for GROMACS
# # ============================================

# # Usage check
# if [ $# -ne 1 ]; then
#     echo "Usage: $0 <target_directory>"
#     exit 1
# fi

# TARGET_DIR="$1"

# # Step 1: Enter the target directory
# cd "$TARGET_DIR" || { echo "Cannot enter directory $TARGET_DIR"; exit 1; }

# # Step 2: Create GMX_PROLIG_MD folder
# mkdir -p GMX_PROLIG_MD

# # Step 3: Find the PDB file containing 'protein' and .pdb extension
# PROTEIN_PDB=$(ls *.pdb 2>/dev/null | head -n 1)
# if [ -z "$PROTEIN_PDB" ]; then
#     echo "No protein.pdb file found."
#     exit 1
# fi

# echo "Protein PDB file found: $PROTEIN_PDB"

# # Step 3.1: Replace HN1, HN2, HN3 atom names
# sed -i 's/HN1/H1 /g; s/HN2/H2 /g; s/HN3/H3 /g' "$PROTEIN_PDB"

# # Step 3.2: Run pdbfixer to clean the protein
# pdbfixer "$PROTEIN_PDB" --output=GMX_PROLIG_MD/fixed_clean_protein.pdb --add-atoms=heavy --keep-heterogens=none

# # Step 4: Enter GMX_PROLIG_MD and generate topology using gmx pdb2gmx
# cd GMX_PROLIG_MD || { echo "Cannot enter GMX_PROLIG_MD"; exit 1; }
# echo -e "6\n1" | gmx pdb2gmx -f fixed_clean_protein.pdb -o pro_lig.gro -ignh

# # =============================== #
# # Step 4.1: Merge ligand GRO coordinates
# # =============================== #

# # Step 4.1.1: Check ligand GRO exists
# if [ ! -f "../LIG.amb2gmx/LIG.gro" ]; then
#     echo "Ligand GRO file not found in ../LIG.amb2gmx"
#     exit 1
# fi

# # Step 4.1.2: Extract ligand atom count and coordinates
# LIGAND_ATOM_COUNT=$(head -n 2 ../LIG.amb2gmx/LIG.gro | tail -n 1 | awk '{print $1}')
# LIGAND_COORDS=$(sed '1,2d' ../LIG.amb2gmx/LIG.gro | sed '$d')

# # Step 4.1.3: Read protein gro atom count
# PROTEIN_ATOM_COUNT=$(head -n 2 pro_lig.gro | tail -n 1 | awk '{print $1}')

# # Step 4.1.4: Calculate new atom count
# TOTAL_ATOM_COUNT=$((PROTEIN_ATOM_COUNT + LIGAND_ATOM_COUNT))

# # Step 4.1.5: Extract protein gro components
# GRO_HEADER=$(head -n 1 pro_lig.gro)
# GRO_BOX=$(tail -n 1 pro_lig.gro)
# PROTEIN_COORDS=$(sed '1,2d' pro_lig.gro | sed '$d')

# # Step 4.1.6: Write the merged gro file
# {
#     echo "$GRO_HEADER"
#     echo "$TOTAL_ATOM_COUNT"
#     echo "$PROTEIN_COORDS"
#     echo "$LIGAND_COORDS"
#     echo "$GRO_BOX"
# } > pro_lig_merged.gro

# # Backup original
# mv pro_lig.gro pro_lig_backup.gro
# mv pro_lig_merged.gro pro_lig.gro

# echo "Ligand coordinates successfully merged into pro_lig.gro"

# # =============================== #
# # Step 5: Return to the parent directory
# # =============================== #

# cd ..

# # =============================== #
# # Step 6: Process ligand topology
# # =============================== #

# cd LIG.amb2gmx || { echo "Cannot enter LIG.amb2gmx"; exit 1; }

# # Step 6.1: Extract [ moleculetype ] to before [ system ] into LIG.itp
# awk '
#     BEGIN { capture=0 }
#     /^\[ *moleculetype *\]/ { capture=1; print; next }
#     /^\[ *system *\]/ { capture=0 }
#     capture { print }
# ' LIG.top > LIG.itp

# # Step 6.2: Append position restraints to LIG.itp
# cat >> LIG.itp <<EOF

# #ifdef POSRES
# #include "posre_LIG.itp"
# #endif
# EOF

# # Step 6.3: Extract full [ atomtypes ] block
# awk '
#     BEGIN { capture=0 }
#     /^\[ *atomtypes *\]/ { capture=1; print; next }
#     /^\[.*\]/ && capture { capture=0 }
#     capture { print }
# ' LIG.top > atomtypes_block.itp

# cd ..

# # =============================== #
# # Step 7: Modify GMX_PROLIG_MD/topol.top
# # =============================== #

# TOPOL="GMX_PROLIG_MD/topol.top"

# if [ ! -f "$TOPOL" ]; then
#     echo "topol.top does not exist in GMX_PROLIG_MD"
#     exit 1
# fi

# # Step 7.1: Insert [ atomtypes ] and #include "LIG.itp" after forcefield include
# awk -v atomtypes="$(cat LIG.amb2gmx/atomtypes_block.itp)" '
#     BEGIN { inserted=0 }
#     {
#         print $0
#         if ($0 ~ /#include "amber99sb.ff\/forcefield.itp"/ && inserted==0) {
#             print ""
#             print atomtypes
#             print ""
#             print "#include \"../LIG.amb2gmx/LIG.itp\""
#             inserted=1
#         }
#     }
# ' "$TOPOL" > temp_topol.top && mv temp_topol.top "$TOPOL"

# # Step 7.2: Add "LIG 1" to [ molecules ] section (if not already present)
# awk '
#     BEGIN { in_molecules=0; added=0 }
#     /^\[ *molecules *\]/ { in_molecules=1; print; next }
#     /^\[/ && in_molecules { in_molecules=0 }
#     { 
#         if (in_molecules && $1 == "LIG") added=1;
#         print 
#     }
#     END {
#         if (!added) {
#             print "LIG\t1"
#         }
#     }
# ' "$TOPOL" > temp_topol.top && mv temp_topol.top "$TOPOL"

# echo "topol.top successfully updated with ligand parameters."
# echo "All steps completed successfully."

# # =============================== #
# # System Modeling with gmx
# # =============================== #

# echo "Modeling with GROMACS..."

# gmx editconf -f pro_lig.gro -o pro_lig_newbox.gro -c -box 10 10 10 
# gmx solvate -cp pro_lig_newbox.gro -p topol.top -o solv.gro
# gmx grompp -f ../../mdps/ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 10
# echo -e "15" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

###########################################
# Version 2.0
###########################################
#!/bin/bash

# ============================================
# Script name: batch_prepare_system.sh
# Purpose: Batch automate protein-ligand system setup for GROMACS
# ============================================

if [ $# -ne 1 ]; then
    echo "Usage: $0 <parent_directory>"
    exit 1
fi

PARENT_DIR="$1"

if [ ! -d "$PARENT_DIR" ]; then
    echo "Error: Directory '$PARENT_DIR' does not exist"
    exit 1
fi

# 记录哪些子文件夹缺少 LIG.amb2gmx
missing_lig_dirs=()

#########################################
# 遍历所有子文件夹
#########################################

for system_dir in "$PARENT_DIR"/*/; do
    [ -d "$system_dir" ] || continue

    echo "----------------------------------------------"
    echo "Processing system folder: $system_dir"

    cd "$system_dir" || { echo "  Cannot enter $system_dir. Skipping."; continue; }

    #########################################
    # Step 1: 找到 {xxx}.pdb 文件（蛋白质）
    #########################################

    pdb_files=( *.pdb )
    if [ ${#pdb_files[@]} -eq 0 ]; then
        echo "  No PDB file (*.pdb) found. Skipping."
        cd ..
        continue
    fi

    PROTEIN_PDB="${pdb_files[0]}"
    pdbid="${PROTEIN_PDB%.pdb}"

    echo "  Found protein PDB file: $PROTEIN_PDB (pdbid = $pdbid)"

    #########################################
    # Step 2: 检查 LIG.amb2gmx 是否存在
    #########################################

    if [ ! -d "LIG.amb2gmx" ]; then
        echo "  LIG.amb2gmx directory not found. Will skip modeling."
        missing_lig_dirs+=("$(basename "$system_dir")")
        cd ..
        continue
    fi

    #########################################
    # Step 3: 创建 GMX_PROLIG_MD 文件夹
    #########################################

    mkdir -p GMX_PROLIG_MD

    #########################################
    # Step 4: 修正 PDB 文件 HN1/HN2/HN3
    #########################################

    sed -i 's/HN1/H1 /g; s/HN2/H2 /g; s/HN3/H3 /g' "$PROTEIN_PDB"

    #########################################
    # Step 5: 用 pdbfixer 修复蛋白
    #########################################

    pdbfixer "$PROTEIN_PDB" --output=GMX_PROLIG_MD/fixed_clean_protein.pdb --add-atoms=heavy --keep-heterogens=none

    #########################################
    # Step 6: pdb2gmx 生成拓扑
    #########################################

    cd GMX_PROLIG_MD || { echo "  Cannot enter GMX_PROLIG_MD. Skipping."; cd ..; continue; }
    echo -e "6\n1" | gmx pdb2gmx -f fixed_clean_protein.pdb -o pro_lig.gro -ignh

    #########################################
    # Step 7: 合并 LIG 坐标
    #########################################

    if [ ! -f "../LIG.amb2gmx/LIG.gro" ]; then
        echo "  Ligand GRO file not found in ../LIG.amb2gmx. Skipping."
        cd ../..
        missing_lig_dirs+=("$(basename "$system_dir") (missing LIG.gro)")
        continue
    fi

    LIGAND_ATOM_COUNT=$(head -n 2 ../LIG.amb2gmx/LIG.gro | tail -n 1 | awk '{print $1}')
    LIGAND_COORDS=$(sed '1,2d' ../LIG.amb2gmx/LIG.gro | sed '$d')

    PROTEIN_ATOM_COUNT=$(head -n 2 pro_lig.gro | tail -n 1 | awk '{print $1}')
    TOTAL_ATOM_COUNT=$((PROTEIN_ATOM_COUNT + LIGAND_ATOM_COUNT))

    GRO_HEADER=$(head -n 1 pro_lig.gro)
    GRO_BOX=$(tail -n 1 pro_lig.gro)
    PROTEIN_COORDS=$(sed '1,2d' pro_lig.gro | sed '$d')

    {
        echo "$GRO_HEADER"
        echo "$TOTAL_ATOM_COUNT"
        echo "$PROTEIN_COORDS"
        echo "$LIGAND_COORDS"
        echo "$GRO_BOX"
    } > pro_lig_merged.gro

    mv pro_lig.gro pro_lig_backup.gro
    mv pro_lig_merged.gro pro_lig.gro

    echo "  Ligand coordinates merged into pro_lig.gro"

    #########################################
    # Step 8: 准备 LIG.itp 和 atomtypes_block
    #########################################

    cd ../LIG.amb2gmx || { echo "  Cannot enter LIG.amb2gmx. Skipping."; cd ..; continue; }

    awk '
        BEGIN { capture=0 }
        /^\[ *moleculetype *\]/ { capture=1; print; next }
        /^\[ *system *\]/ { capture=0 }
        capture { print }
    ' LIG.top > LIG.itp

    cat >> LIG.itp <<EOF

#ifdef POSRES
#include "posre_LIG.itp"
#endif
EOF

    awk '
        BEGIN { capture=0 }
        /^\[ *atomtypes *\]/ { capture=1; print; next }
        /^\[.*\]/ && capture { capture=0 }
        capture { print }
    ' LIG.top > atomtypes_block.itp

    cd ../GMX_PROLIG_MD

    #########################################
    # Step 9: 修改 topol.top
    #########################################

    TOPOL="topol.top"

    if [ ! -f "$TOPOL" ]; then
        echo "  topol.top does not exist."
        cd ../..
        missing_lig_dirs+=("$(basename "$system_dir") (missing topol.top)")
        continue
    fi

    awk -v atomtypes="$(cat ../LIG.amb2gmx/atomtypes_block.itp)" '
        BEGIN { inserted=0 }
        {
            print $0
            if ($0 ~ /#include "amber99sb.ff\/forcefield.itp"/ && inserted==0) {
                print ""
                print atomtypes
                print ""
                print "#include \"../LIG.amb2gmx/LIG.itp\""
                inserted=1
            }
        }
    ' "$TOPOL" > temp_topol.top && mv temp_topol.top "$TOPOL"

    awk '
        BEGIN { in_molecules=0; added=0 }
        /^\[ *molecules *\]/ { in_molecules=1; print; next }
        /^\[/ && in_molecules { in_molecules=0 }
        { 
            if (in_molecules && $1 == "LIG") added=1;
            print 
        }
        END {
            if (!added) {
                print "LIG\t1"
            }
        }
    ' "$TOPOL" > temp_topol.top && mv temp_topol.top "$TOPOL"

    echo "  topol.top updated."

    #########################################
    # Step 10: Box & solvation & ions
    #########################################

    echo "  Running box, solvate, ions..."

    gmx editconf -f pro_lig.gro -o pro_lig_newbox.gro -c -box 10 10 10 
    gmx solvate -cp pro_lig_newbox.gro -p topol.top -o solv.gro
    gmx grompp -f ../../mdps/ions.mdp -c solv.gro -p topol.top -o ions.tpr -maxwarn 10
    echo -e "15" | gmx genion -s ions.tpr -o solv_ions.gro -p topol.top -pname NA -nname CL -neutral

    echo "  System preparation completed for $system_dir"

    #########################################
    # Step 11: 返回上一级继续下一个
    #########################################

    cd ../..

done

#########################################
# Summary 报告
#########################################

echo ""
echo "=============================================="
echo "Summary: Systems missing LIG.amb2gmx or failed steps:"
if [ ${#missing_lig_dirs[@]} -eq 0 ]; then
    echo "All systems processed successfully!"
else
    for miss in "${missing_lig_dirs[@]}"; do
        echo " - $miss"
    done
fi
echo "=============================================="
