# LigTraj V1.1 -- By A.I.B.

LigTraj is a Python package for ligand trajectory analysis including RMSD, contact analysis, and covariance analysis etc.

## Installation
Download the package, enter `LigTraj`, then use `pip` for installation
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

## Demo
![ligand_rmsd](https://github.com/user-attachments/assets/d5999e30-b60f-492c-8cf8-27ac240bcecc)
![ligand_residue_contact_network](https://github.com/user-attachments/assets/0d3cec58-48da-474c-9ccc-4673dc5a3d09)
![ligand_contact_frequency](https://github.com/user-attachments/assets/a3117df8-312c-48bb-ac1a-0d8c6227f467)
![covariance_matrix](https://github.com/user-attachments/assets/c7999b11-0992-4ca7-b5f6-9fa410313e31)

