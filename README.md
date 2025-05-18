# LigTraj V1.1 -- By A.I.B.

LigTraj is a Python package for ligand trajectory analysis including RMSD, contact analysis, and covariance analysis etc. Graph module in LigTraj is designed to generate feature of MD traj. Recent updated feature used for embedding includes MaSIF, distance-contact graph, conformationla ensamble etc. 

![726cdd7e513e76386b2c0755fb1ca5c](https://github.com/user-attachments/assets/3a7103c1-4e40-4498-82a1-d6a325c00f18)


And the GAFFMaker tool is used to streamline the process of setting up simulation files for ligand-protein system molecular dynamics(MD) simulation.

## Installation
Download the package, enter `LigTraj`, then use `pip` for installation
```bash
git clone https://github.com/aib001/LigTraj.git #you can use ssh for a quicker fetch
cd LigTrajV1.3
pip install -e .
```

## Usage

```python

from LigTraj import TrajAnalysis as ta
from LigTraj import Graph
import os

base_dir = os.path.dirname(os.path.abspath(__file__))

topol = os.path.join(base_dir, "GMX_PROLIG_MD", "solv_ions.gro")
traj = os.path.join(base_dir, "GMX_PROLIG_MD", "prod", "md_aligned.xtc")
sdf = os.path.join(base_dir, "GMX_PROLIG_MD", "v2020_3tiy_ligand_1746191494002.sdf")

##################################
# Part 1. Traj Analysis Example
##################################

print("Running t-SNE analysis (SE3 invariant, coordinates)...")
ta.tsne(topol, traj, resname="LIG", feature_type="coordinates", se3_invariant=True)

print("Running RMSD analysis...")
ta.rmsd(topol, traj, resname="LIG")

print("Running Contact analysis...")
ta.contact(topol, traj, sdf, resname="LIG", distance_cutoff=0.4, n_frames=50)

print("Running Covariance analysis...")
ta.covariance(topol, traj, resname="LIG")

##################################
# Part 2. Graph Generate Test
##################################

print("Building Graph ensemble...")
Graph.build(topol, traj, sdf, resname="LIG", n_frames=10)

##################################
# Part 3. Graph Feature Embedding Test 
##################################

print("Generating MaSIF embeddings...")
Graph.feature(topol, traj, sdf, resname="LIG", n_frames=10)


# print("Running MaSIF geodesic embedding...")
# Graph.masif_embedding(topol, traj, sdf, resname="LIG", n_frames=20)

##################################
# Part 4. MaSIF Embedding
##################################

print("Running MaSIF-style geodesic embedding (Euclidean distance)...")
Graph.masif_embedding(
    topol, traj, sdf,
    resname="LIG",
    cutoff=0.4,
    n_frames=10,
    distance_mode="euclidean"
)

print("Running MaSIF-style geodesic embedding (Geodesic distance)...")
Graph.masif_embedding(
    topol, traj, sdf,
    resname="LIG",
    cutoff=0.4,
    n_frames=10,
    distance_mode="geodesic"
)


##################################
print("All analyses completed.")

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
- GROMACS recommended version: 2024.2 with CUDA
- AmberTool recommonded version higher than V20

## Demo
### TrajAnalysis
![ligand_rmsd](https://github.com/user-attachments/assets/d5999e30-b60f-492c-8cf8-27ac240bcecc)
![ligand_residue_contact_network](https://github.com/user-attachments/assets/0d3cec58-48da-474c-9ccc-4673dc5a3d09)
![ligand_contact_frequency](https://github.com/user-attachments/assets/a3117df8-312c-48bb-ac1a-0d8c6227f467)
![covariance_matrix](https://github.com/user-attachments/assets/c7999b11-0992-4ca7-b5f6-9fa410313e31)
### Graph
![Curvature_polar_plot_euclidean](https://github.com/user-attachments/assets/4fbb31d2-5909-41be-ba29-73a95af95236)
![Curvature_polar_plot_geodesic](https://github.com/user-attachments/assets/2b67f1e2-f635-44f2-ae48-7c77103aaefa)
![graph_ensemble](https://github.com/user-attachments/assets/c55306dc-2c57-4ea3-9e1f-b8d525763867)



