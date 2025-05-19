# LigTraj V1.3 -- By A.I.B.
---------------------------------------------------------------------
## V1.3 Update Feature -- May 19, 2025
+ Adjust the framework
+ GAFF_Maker is now available through python script:
  ```python
  from LigTraj import GAFF_Maker as GAM
  model=GAM.GAFF_Maker(protein_path, ligand_path, output_path, overwrite=False)
  model.run()
  ```
+ Add 4D-fingerprint module 
---------------------------------------------------------------------

LigTraj is a Python package for ligand trajectory analysis including RMSD, contact analysis, and covariance analysis etc. Graph module in LigTraj is designed to generate feature of MD traj. Recent updated feature used for embedding includes MaSIF, distance-contact graph, conformational ensemble etc. 

![726cdd7e513e76386b2c0755fb1ca5c](https://github.com/user-attachments/assets/3a7103c1-4e40-4498-82a1-d6a325c00f18)

The GAFFMaker tool streamlines the process of setting up simulation files for ligand-protein system molecular dynamics (MD) simulation.

## Environment Setup

### Required Software
- **Python**: 3.8 or higher
- **GROMACS**: Version 2024.2 with CUDA support
- **AmberTools**: Version 20 or higher (recommended V23)
- **RDKit**: For chemical structure handling

### Installation

<details>
<summary><b>Step 1: Install Required Software</b> (click to expand)</summary>

```bash
# Install AmberTools (using conda)
conda install -c conda-forge ambertools=23

# If you don't have GROMACS installed:
# For Ubuntu/Debian:
sudo apt-get install gromacs

# For macOS:
brew install gromacs

# For more options including CUDA support, refer to:
# https://manual.gromacs.org/current/install-guide/index.html
```
</details>

<details>
<summary><b>Step 2: Clone and Install LigTraj</b> (click to expand)</summary>

```bash
# Clone the repository
git clone https://github.com/aib001/LigTraj.git

# Navigate to the directory
cd LigTrajV1.3

# Install the package and dependencies
pip install -e .
```
</details>

### Dependencies
```
mdtraj>=1.9.7
numpy>=1.20.0
pandas>=1.3.0
matplotlib>=3.5.0
networkx>=2.7.0
rdkit>=2022.03.1
tqdm>=4.62.0
scipy>=1.7.0
scikit-learn>=1.0.0
torch>=1.10.0
torch_geometric>=2.0.0
openpyxl>=3.0.0
```

## Usage Examples

### 1. Trajectory Analysis

<div style="background-color:#f8f8f8; padding:10px; border-radius:5px; margin-bottom:20px;">

```python
from LigTraj import TrajAnalysis as ta
import os

# Set paths
topol = "path/to/solv_ions.gro"
traj = "path/to/prod/md_aligned.xtc"
sdf = "path/to/ligand.sdf"

# Run analyses
ta.rmsd(topol, traj, resname="LIG", heavy_only=True, frame_interval_ns=0.5)
ta.contact(topol, traj, sdf, resname="LIG", distance_cutoff=0.4, n_frames=50)
ta.covariance(topol, traj, resname="LIG")
```
</div>

### 2. Protein-Ligand System Setup with GAFF_Maker

<div style="background-color:#f8f8f8; padding:10px; border-radius:5px; margin-bottom:20px;">

```python
from LigTraj import GAFF_Maker as GAM

# Option 1: Simple function call
GAM.GAFF_Maker("protein.pdb", "ligand.mol2", "output_dir", overwrite=False)

# Option 2: More control over the process
model = GAM.GAFFMaker("protein.pdb", "ligand.mol2", "output_dir", overwrite=False)
model.clean_protein()
model.generate_ligand_forcefield()
model.standardize_files()
model.build_model()
model.cleanup()
```
</div>

### 3. Graph and Feature Generation

<div style="background-color:#f8f8f8; padding:10px; border-radius:5px; margin-bottom:20px;">

```python
from LigTraj import Graph
from LigTraj import FourDFingerprint

# Generate contact graphs
Graph.build(topol, traj, sdf, resname="LIG", n_frames=10)

# Generate MaSIF embeddings
Graph.feature(topol, traj, sdf, resname="LIG", n_frames=10)

# Generate geodesic embeddings
Graph.masif_embedding(
    topol, traj, sdf,
    resname="LIG",
    cutoff=0.4,
    n_frames=10,
    distance_mode="geodesic"
)

# Compute 4D fingerprint
FourDFingerprint.compute_4D_fingerprint(
    topol, traj,
    resname="LIG",
    grid_size=0.5,
    cutoff=0.4,
    n_frames=50
)
```
</div>

### 4. Complete Analysis Pipeline

<div style="background-color:#f8f8f8; padding:10px; border-radius:5px; margin-bottom:20px;">

```python
from LigTraj import TrajAnalysis as ta
from LigTraj import Graph
from LigTraj import FourDFingerprint
import os

# Setup paths
base_dir = os.path.dirname(os.path.abspath(__file__))
topol = os.path.join(base_dir, "GMX_PROLIG_MD", "solv_ions.gro")
traj = os.path.join(base_dir, "GMX_PROLIG_MD", "prod", "md_aligned.xtc")
sdf = os.path.join(base_dir, "GMX_PROLIG_MD", "ligand.sdf")

# Part 1: Basic Trajectory Analysis
ta.rmsd(topol, traj, resname="LIG")
ta.contact(topol, traj, sdf, resname="LIG", distance_cutoff=0.4, n_frames=50)
ta.covariance(topol, traj, resname="LIG")

# Part 2: Graph Generation
Graph.build(topol, traj, sdf, resname="LIG", n_frames=10)

# Part 3: Feature Embedding
Graph.feature(topol, traj, sdf, resname="LIG", n_frames=10)

# Part 4: MaSIF Embedding (with both distance modes)
Graph.masif_embedding(topol, traj, sdf, resname="LIG", cutoff=0.4, 
                     n_frames=10, distance_mode="euclidean")
Graph.masif_embedding(topol, traj, sdf, resname="LIG", cutoff=0.4, 
                     n_frames=10, distance_mode="geodesic")

# Part 5: Generate 4D Fingerprint
FourDFingerprint.compute_4D_fingerprint(topol, traj, resname="LIG")

print("All analyses completed successfully!")
```
</div>

## Running Simulations with GAFFMaker

The workflow for setting up and running protein-ligand MD simulations:

<div style="background-color:#e8f4ff; padding:15px; border-radius:5px; margin-bottom:20px;">

1. **Prepare Input Files**:
   - Protein structure in PDB format
   - Ligand structure in MOL2 format

2. **Generate System with GAFFMaker**:
   ```python
   from LigTraj import GAFF_Maker as GAM
   GAM.GAFF_Maker("protein.pdb", "ligand.mol2", "system_dir")
   ```

3. **Run Simulation with GROMACS**:
   ```bash
   # Navigate to the output directory
   cd system_dir/GMX_PROLIG_MD
   
   # Option 1: Run with the provided script
   bash localrun.sh
   
   # Option 2: Run on HPC with SLURM
   sbatch iqb.md.slurm.sh
   
   # Option 3: Submit multiple jobs
   bash batch_submit.sh
   ```

4. **Analyze Trajectory**:
   ```python
   from LigTraj import TrajAnalysis as ta
   
   # Align trajectory first (if not already done)
   # bash pbc_align.sh system_dir  # Script included in GAFFMakerV1.1
   
   ta.rmsd("system_dir/GMX_PROLIG_MD/solv_ions.gro", 
          "system_dir/GMX_PROLIG_MD/prod/md_aligned.xtc",
          resname="LIG")
   ```
</div>

## Demo

<div style="display:flex; flex-wrap:wrap; justify-content:space-between;">

<div style="width:48%; margin-bottom:20px;">
<h3>RMSD Analysis</h3>
<img src="https://github.com/user-attachments/assets/d5999e30-b60f-492c-8cf8-27ac240bcecc" alt="ligand_rmsd" style="width:100%;">
</div>

<div style="width:48%; margin-bottom:20px;">
<h3>Ligand-Protein Contact Network</h3>
<img src="https://github.com/user-attachments/assets/0d3cec58-48da-474c-9ccc-4673dc5a3d09" alt="ligand_residue_contact_network" style="width:100%;">
</div>

<div style="width:48%; margin-bottom:20px;">
<h3>Ligand Contact Frequency</h3>
<img src="https://github.com/user-attachments/assets/a3117df8-312c-48bb-ac1a-0d8c6227f467" alt="ligand_contact_frequency" style="width:100%;">
</div>

<div style="width:48%; margin-bottom:20px;">
<h3>Covariance Matrix</h3>
<img src="https://github.com/user-attachments/assets/c7999b11-0992-4ca7-b5f6-9fa410313e31" alt="covariance_matrix" style="width:100%;">
</div>

<div style="width:48%; margin-bottom:20px;">
<h3>MaSIF Curvature (Euclidean)</h3>
<img src="https://github.com/user-attachments/assets/4fbb31d2-5909-41be-ba29-73a95af95236" alt="Curvature_polar_plot_euclidean" style="width:100%;">
</div>

<div style="width:48%; margin-bottom:20px;">
<h3>MaSIF Curvature (Geodesic)</h3>
<img src="https://github.com/user-attachments/assets/2b67f1e2-f635-44f2-ae48-7c77103aaefa" alt="Curvature_polar_plot_geodesic" style="width:100%;">
</div>

<div style="width:48%; margin-bottom:20px;">
<h3>Ensemble Contact Graph</h3>
<img src="https://github.com/user-attachments/assets/c55306dc-2c57-4ea3-9e1f-b8d525763867" alt="graph_ensemble" style="width:100%;">
</div>

</div>
