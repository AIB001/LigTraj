#!/bin/bash
# Usage pbc_align.sh <partne_path>
if [[ $# -ne 1 ]]; then
    echo "Parent path needed"
    exit 1
fi

path="${1%/}"  
cd "${path}/GMX_PROLIG_MD/prod"

# echo -e "1\n0" | gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center
gmx trjconv -s md.tpr -f md.xtc -o md_noPBC.xtc -pbc mol -center <<EOF
Protein
System
EOF

# echo -e "1\n0" | gmx trjconv -s md.tpr -f md_noPBC.xtc -o md_aligned.xtc -fit rot+trans
gmx trjconv -s md.tpr -f md_noPBC.xtc -o md_aligned.xtc -fit rot+trans <<EOF
Protein
System
EOF

echo "PBC is canceled and protein is aligned"
