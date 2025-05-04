# import argparse
# import os
# import mdtraj as md
# import numpy as np
# import pandas as pd
# from rdkit import Chem, RDLogger
# from rdkit.Chem import rdmolops

# # Disable RDKit warnings for cleaner output
# RDLogger.DisableLog('rdApp.*')

# def compute_masif_features(topol, traj, sdf, resname="LIG", cutoff=0.4, n_frames=50):
#     """
#     Compute MaSIF-style surface features for ligand-protein contacts within 4Å.
#     Features per frame:
#       - shape_index_est: approximate shape index (0 ~ flat, 1 ~ highly curved interface)
#       - curvature_est: approximate curvature (relative magnitude of second principal component)
#       - electrostatic_score: count of complementary charged contacts (approximate)
#       - hbond_count: number of potential H-bond ligand-protein atom pairs (N-O within 3.5Å)
#       - hydrophobic_score: fraction of contacting residues that are hydrophobic
#     Saves an Excel file with these features for all frames.
#     Returns a numpy array of shape (n_frames, 5) with the features.
#     """
#     # Load trajectory and limit frames
#     t = md.load(traj, top=topol)
#     total_frames = t.n_frames
#     if n_frames > total_frames:
#         n_frames = total_frames
    
#     # Select heavy atoms for ligand and protein
#     ligand_atoms = t.topology.select(f"resname {resname} and element != H")
#     protein_atoms = t.topology.select("protein and element != H")
#     if len(ligand_atoms) == 0 or len(protein_atoms) == 0:
#         raise ValueError("Ligand or protein heavy atoms not found for feature computation.")
    
#     # Calculate ligand formal charge using RDKit on the ligand SDF
#     mol = Chem.MolFromMolFile(sdf, removeHs=False)
#     if mol is None:
#         raise ValueError("Failed to read ligand from SDF for feature computation.")
#     ligand_charge = rdmolops.GetFormalCharge(mol)
    
#     # Residue type categories for chemical features
#     positive_residues = {"ARG", "LYS"}
#     negative_residues = {"ASP", "GLU"}
#     hydrophobic_residues = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO", "CYS", "TYR", "GLY"}
    
#     # Use mdtraj to find contacting protein atoms for each frame
#     sub_traj = t.slice(range(n_frames))
#     neighbors_list = md.compute_neighbors(sub_traj, cutoff, ligand_atoms, haystack_indices=protein_atoms)
    
#     features = []
#     for frame_idx, prot_indices in enumerate(neighbors_list):
#         prot_indices = np.array(prot_indices, dtype=int)
#         # Geometric features: shape index and curvature
#         shape_score = 0.0
#         curvature_score = 0.0
#         if prot_indices.size >= 3:
#             coords = t.xyz[frame_idx, prot_indices]  # coordinates of contacting protein atoms
#             center = coords.mean(axis=0)
#             centered = coords - center
#             if centered.shape[0] >= 3:
#                 cov = np.cov(centered.T)
#                 try:
#                     eigvals = np.linalg.eigvals(cov)
#                     eigvals = np.sort(np.real(eigvals))
#                     if eigvals.size < 3:
#                         eigvals = np.pad(eigvals, (0, 3 - eigvals.size), 'constant', constant_values=0.0)
#                     e1, e2, e3 = eigvals[0], eigvals[1], eigvals[2]
#                     if abs(e3) > 1e-8:
#                         shape_score = float(e1 / (e3 + 1e-8))
#                         curvature_score = float(e2 / (e3 + 1e-8))
#                         # Clamp values to [0, 1] range for stability
#                         shape_score = max(0.0, min(1.0, shape_score))
#                         curvature_score = max(0.0, min(1.0, curvature_score))
#                 except np.linalg.LinAlgError:
#                     shape_score = curvature_score = 0.0
#         # Chemical features: charged contacts, H-bonds, hydrophobicity
#         contacting_residues = {t.topology.atom(idx).residue for idx in prot_indices}
#         pos_res_count = sum(1 for res in contacting_residues if res.name in positive_residues)
#         neg_res_count = sum(1 for res in contacting_residues if res.name in negative_residues)
#         hydrophobic_res_count = sum(1 for res in contacting_residues if res.name in hydrophobic_residues)
#         # Electrostatic score: count of opposite-charge contacts
#         if ligand_charge > 0:
#             electro_score = neg_res_count
#         elif ligand_charge < 0:
#             electro_score = pos_res_count
#         else:
#             electro_score = pos_res_count + neg_res_count
#         # Hydrogen bond count: count N-O atom pairs < 0.35 nm
#         hbond_count = 0
#         ligand_polar_atoms = [i for i in ligand_atoms if t.topology.atom(i).element.symbol in ('N', 'O')]
#         protein_polar_atoms = [j for j in prot_indices if t.topology.atom(j).element.symbol in ('N', 'O')]
#         polar_pairs = []
#         for i in ligand_polar_atoms:
#             lig_elem = t.topology.atom(i).element.symbol
#             for j in protein_polar_atoms:
#                 prot_elem = t.topology.atom(j).element.symbol
#                 if (lig_elem == 'N' and prot_elem == 'O') or (lig_elem == 'O' and prot_elem == 'N'):
#                     polar_pairs.append((i, j))
#         if polar_pairs:
#             distances = md.compute_distances(t[frame_idx], polar_pairs)[0]
#             hbond_count = int(np.sum(distances < 0.35))
#         # Hydrophobic score: fraction of contacting residues that are hydrophobic
#         hydrophobic_score = 0.0
#         if contacting_residues:
#             hydrophobic_score = hydrophobic_res_count / len(contacting_residues)
#         features.append([shape_score, curvature_score, electro_score, hbond_count, hydrophobic_score])
    
#     features = np.array(features)
#     # Save features to Excel
#     out_dir = os.path.dirname(os.path.abspath(traj))
#     graph_dir = os.path.join(out_dir, "Graph")
#     os.makedirs(graph_dir, exist_ok=True)
#     df = pd.DataFrame(features, columns=["ShapeIndex", "Curvature", "Electrostatics", "HBonds", "Hydrophobicity"])
#     df.index.name = "Frame"
#     feature_path = os.path.join(graph_dir, "ligand_surface_features.xlsx")
#     df.to_excel(feature_path, index=True)
#     print(f"Surface feature embeddings saved to: {feature_path}")
#     return features

# if __name__ == "__main__":
#     parser = argparse.ArgumentParser(description="Compute MaSIF-style surface features for ligand-protein contacts")
#     parser.add_argument('--topol', required=True, help='Topology file (PDB/PRMTOP)')
#     parser.add_argument('--traj', required=True, help='Trajectory file (XTC/DCD)')
#     parser.add_argument('--sdf', required=True, help='Ligand structure file (.sdf)')
#     parser.add_argument('--resname', default="LIG", help='Ligand residue name (default: LIG)')
#     parser.add_argument('--n_frames', type=int, default=10, help='Number of frames to analyze (default: 50)')
#     args = parser.parse_args()
#     compute_masif_features(args.topol, args.traj, args.sdf, resname=args.resname, cutoff=0.4, n_frames=args.n_frames)

import argparse
import os
import mdtraj as md
import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdmolops

# Disable RDKit warnings for cleaner output
RDLogger.DisableLog('rdApp.*')

def compute_masif_features(topol, traj, sdf, resname="LIG", cutoff=0.4, n_frames=50):
    """
    Compute MaSIF-style surface features for ligand-protein contacts within 4Å.
    Features per frame:
      - shape_index_est
      - curvature_est
      - electrostatic_score
      - hbond_count
      - hydrophobic_score
    """
    # Load trajectory
    t = md.load(traj, top=topol)
    total_frames = t.n_frames
    if n_frames > total_frames:
        n_frames = total_frames

    ligand_atoms = t.topology.select(f"resname {resname} and element != H")
    protein_atoms = t.topology.select("protein and element != H")
    if len(ligand_atoms) == 0 or len(protein_atoms) == 0:
        raise ValueError("Ligand or protein heavy atoms not found.")

    mol = Chem.MolFromMolFile(sdf, removeHs=False)
    if mol is None:
        raise ValueError("Failed to read ligand SDF.")
    ligand_charge = rdmolops.GetFormalCharge(mol)

    positive_residues = {"ARG", "LYS"}
    negative_residues = {"ASP", "GLU"}
    hydrophobic_residues = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO", "CYS", "TYR", "GLY"}

    sub_traj = t.slice(range(n_frames))
    neighbors_list = md.compute_neighbors(sub_traj, cutoff, ligand_atoms, haystack_indices=protein_atoms)

    features = []
    for frame_idx, prot_indices in enumerate(neighbors_list):
        prot_indices = np.array(prot_indices, dtype=int)

        shape_score = 0.0
        curvature_score = 0.0
        if prot_indices.size >= 3:
            coords = t.xyz[frame_idx, prot_indices]
            center = coords.mean(axis=0)
            centered = coords - center
            if centered.shape[0] >= 3:
                cov = np.cov(centered.T)
                try:
                    eigvals = np.linalg.eigvals(cov)
                    eigvals = np.sort(np.real(eigvals))
                    if eigvals.size < 3:
                        eigvals = np.pad(eigvals, (0, 3 - eigvals.size), 'constant', constant_values=0.0)
                    e1, e2, e3 = eigvals[0], eigvals[1], eigvals[2]
                    if abs(e3) > 1e-8:
                        shape_score = float(e1 / (e3 + 1e-8))
                        curvature_score = float(e2 / (e3 + 1e-8))
                        shape_score = max(0.0, min(1.0, shape_score))
                        curvature_score = max(0.0, min(1.0, curvature_score))
                except np.linalg.LinAlgError:
                    shape_score = curvature_score = 0.0

        contacting_residues = {t.topology.atom(idx).residue for idx in prot_indices}
        pos_res_count = sum(1 for res in contacting_residues if res.name in positive_residues)
        neg_res_count = sum(1 for res in contacting_residues if res.name in negative_residues)
        hydrophobic_res_count = sum(1 for res in contacting_residues if res.name in hydrophobic_residues)

        if ligand_charge > 0:
            electro_score = neg_res_count
        elif ligand_charge < 0:
            electro_score = pos_res_count
        else:
            electro_score = pos_res_count + neg_res_count

        hbond_count = 0
        ligand_polar_atoms = [i for i in ligand_atoms if t.topology.atom(i).element.symbol in ('N', 'O')]
        protein_polar_atoms = [j for j in prot_indices if t.topology.atom(j).element.symbol in ('N', 'O')]
        polar_pairs = [(i, j) for i in ligand_polar_atoms for j in protein_polar_atoms
                       if (t.topology.atom(i).element.symbol, t.topology.atom(j).element.symbol) in [('N', 'O'), ('O', 'N')]]
        if polar_pairs:
            distances = md.compute_distances(t[frame_idx], polar_pairs)[0]
            hbond_count = int(np.sum(distances < 0.35))

        hydrophobic_score = 0.0
        if contacting_residues:
            hydrophobic_score = hydrophobic_res_count / len(contacting_residues)

        features.append([shape_score, curvature_score, electro_score, hbond_count, hydrophobic_score])

    features = np.array(features)

    out_dir = os.path.dirname(os.path.abspath(traj))
    graph_dir = os.path.join(out_dir, "Graph")
    os.makedirs(graph_dir, exist_ok=True)
    df = pd.DataFrame(features, columns=["ShapeIndex", "Curvature", "Electrostatics", "HBonds", "Hydrophobicity"])
    df.index.name = "Frame"
    feature_path = os.path.join(graph_dir, "ligand_surface_features.xlsx")
    df.to_excel(feature_path, index=True)
    print(f"Surface feature embeddings saved to: {feature_path}")
    return features

# For command line usage
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compute MaSIF-style surface features for ligand-protein contacts")
    parser.add_argument('--topol', required=True, help='Topology file (PDB/PRMTOP)')
    parser.add_argument('--traj', required=True, help='Trajectory file (XTC/DCD)')
    parser.add_argument('--sdf', required=True, help='Ligand structure file (.sdf)')
    parser.add_argument('--resname', default="LIG", help='Ligand residue name (default: LIG)')
    parser.add_argument('--n_frames', type=int, default=10, help='Number of frames to analyze (default: 50)')
    args = parser.parse_args()

    compute_masif_features(args.topol, args.traj, args.sdf,
                           resname=args.resname, cutoff=0.4, n_frames=args.n_frames)
