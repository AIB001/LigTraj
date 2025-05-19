# LigTraj V1.3 -- By A.I.B.

LigTraj is a Python package for ligand trajectory analysis including RMSD, contact analysis, and covariance analysis etc.

## Installation

```bash
pip install -e .
```

## Usage

```python
import LigTraj
from LigTraj import TrajAnalysis as ta

ta.rmsd("topol.gro", "traj.xtc", resname="LIG")
ta.contact("topol.gro", "traj.xtc", "ligand.sdf", resname="LIG", n_frames=50)
ta.covariance("topol.gro", "traj.xtc", resname="LIG")
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
- sklearn
