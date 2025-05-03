
from LigTraj import TrajAnalysis as ta
import os

base_dir = os.path.dirname(os.path.abspath(__file__))

topol = os.path.join(base_dir, "GMX_PROLIG_MD", "solv_ions.gro")
traj = os.path.join(base_dir, "GMX_PROLIG_MD", "prod", "md_aligned.xtc")
sdf = os.path.join(base_dir, "GMX_PROLIG_MD", "v2020_3tiy_ligand_1746191494002.sdf")


print("Running t-SNE analysis (SE3 invariant, coordinates)...")
ta.tsne(topol, traj, resname="LIG", feature_type="coordinates", se3_invariant=True)

print("Running RMSD analysis...")
ta.rmsd(topol, traj, resname="LIG")

print("Running Contact analysis...")
ta.contact(topol, traj, sdf, resname="LIG", distance_cutoff=0.4, n_frames=50)

print("Running Covariance analysis...")
ta.covariance(topol, traj, resname="LIG")

print("All analyses completed.")
