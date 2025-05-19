# LigTraj V1.3 -- By A.I.B.

LigTraj is a Python package for ligand trajectory analysis including RMSD, contact analysis, and covariance analysis etc.

## Installation

```bash
pip install -e .
```

## Usage
### TrajAnalysis

```python
import LigTraj
from LigTraj import TrajAnalysis as ta
import os

base_dir = os.path.dirname(os.path.abspath(__file__))
topol = os.path.join(base_dir, "GMX_PROLIG_MD", "solv_ions.gro")
traj = os.path.join(base_dir, "GMX_PROLIG_MD", "prod", "md_aligned.xtc")
sdf = os.path.join(base_dir, "GMX_PROLIG_MD", "v2020_3tiy_ligand_1746191494002.sdf")

ta.rmsd("topol.gro", "traj.xtc", resname="LIG")
ta.contact("topol.gro", "traj.xtc", "ligand.sdf", resname="LIG", n_frames=50)
ta.covariance("topol.gro", "traj.xtc", resname="LIG")
```

### GAFF_Maker
```python
from LigTraj import GAFF_Maker as GAM
import argparse

if __name__ == "__main__":
    """Main function to parse arguments and run GAFFMaker"""
    parser = argparse.ArgumentParser(description="GAFFMaker: Build protein-ligand systems for GROMACS")
    
    parser.add_argument("protein", help="Path to protein PDB file")
    parser.add_argument("ligand", help="Path to ligand MOL2 file")
    parser.add_argument("--output", "-o", default="output", help="Output directory")
    parser.add_argument("--overwrite", "-f", action="store_true", help="Overwrite existing files")
    
    args = parser.parse_args()
    
    # Create and run GAFFMaker
    gaffmaker = GAM.GAFF_Maker(args.protein, args.ligand, args.output)
    gaffmaker.run()
```


## Requirements
- mdtraj
- numpy
- pandas
- matplotlib
- networkx
- rdkit
- tqdm
- scipy
- scikit-learn
