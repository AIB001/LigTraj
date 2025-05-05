# import os
# import mdtraj as md
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from rdkit import Chem
# from rdkit.Chem import rdmolops
# from tqdm import tqdm

# plt.rcParams.update({'font.size': 21, 'font.family': 'Times New Roman'})

# def compute_masif_geodesic_embedding(topol, traj_file, sdf, resname="LIG", cutoff=0.4, n_frames=10):
#     """
#     Compute MaSIF-like geodesic feature embedding for ligand-protein interface patches.
#     Outputs per-frame geodesic patch embeddings and visualizations for all features.
#     """
#     traj = md.load(traj_file, top=topol)
#     if n_frames > traj.n_frames:
#         n_frames = traj.n_frames

#     # Ligand and protein atoms
#     ligand_atoms = traj.topology.select(f"resname {resname} and element != H")
#     protein_atoms = traj.topology.select("protein and element != H")
#     if len(ligand_atoms) == 0 or len(protein_atoms) == 0:
#         raise ValueError("No ligand or protein heavy atoms found.")

#     mol = Chem.MolFromMolFile(sdf, removeHs=False)
#     if mol is None:
#         raise ValueError("Failed to read ligand SDF.")
#     ligand_charge = rdmolops.GetFormalCharge(mol)

#     pos_residues = {"ARG", "LYS"}
#     neg_residues = {"ASP", "GLU"}
#     hydrophobic_residues = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO", "CYS", "TYR", "GLY"}

#     out_dir = os.path.dirname(os.path.abspath(traj_file))
#     graph_dir = os.path.join(out_dir, "Graph")
#     masif_dir = os.path.join(graph_dir, "masif")
#     os.makedirs(masif_dir, exist_ok=True)

#     all_embeddings = []

#     for frame in tqdm(range(n_frames), desc="Computing MaSIF embeddings"):
#         t_slice = traj.slice(frame)

#         lig_coords = t_slice.xyz[0, ligand_atoms, :]
#         lig_center = lig_coords.mean(axis=0)

#         distances = md.compute_distances(t_slice,
#                         [(i, j) for i in ligand_atoms for j in protein_atoms])[0]
#         contact_pairs = [(ligand_atoms[i // len(protein_atoms)],
#                           protein_atoms[i % len(protein_atoms)])
#                           for i, d in enumerate(distances) if d < cutoff]
#         contacting_prot_atoms = list(set(j for _, j in contact_pairs))

#         if len(contacting_prot_atoms) < 3:
#             continue

#         prot_coords = t_slice.xyz[0, contacting_prot_atoms, :]
#         center = prot_coords.mean(axis=0)

#         r = np.linalg.norm(prot_coords - center, axis=1)

#         direction = lig_center - center
#         direction /= (np.linalg.norm(direction) + 1e-8)

#         vectors = prot_coords - center
#         u_axis = vectors[0] / (np.linalg.norm(vectors[0]) + 1e-8)
#         v_axis = np.cross(direction, u_axis)
#         v_axis /= (np.linalg.norm(v_axis) + 1e-8)

#         theta = []
#         for vec in vectors:
#             proj_u = np.dot(vec, u_axis)
#             proj_v = np.dot(vec, v_axis)
#             angle = np.arctan2(proj_v, proj_u) % (2 * np.pi)
#             theta.append(angle)

#         # ---- Feature per residue ----
#         shape_list = []
#         curvature_list = []
#         electro_list = []
#         hbond_list = []
#         hydro_list = []

#         for idx, prot_idx in enumerate(contacting_prot_atoms):
#             res = t_slice.topology.atom(prot_idx).residue

#             # Local PCA for shape/curvature
#             neighbor_coords = []
#             for jdx, other_idx in enumerate(contacting_prot_atoms):
#                 d = np.linalg.norm(prot_coords[idx] - prot_coords[jdx])
#                 if d < 0.5:
#                     neighbor_coords.append(prot_coords[jdx])
#             neighbor_coords = np.array(neighbor_coords)

#             if len(neighbor_coords) >= 3:
#                 centered = neighbor_coords - neighbor_coords.mean(axis=0)
#                 cov = np.cov(centered.T)
#                 eigvals = np.linalg.eigvalsh(cov)
#                 eigvals = np.sort(eigvals)
#                 if len(eigvals) < 3:
#                     eigvals = np.pad(eigvals, (0, 3 - len(eigvals)), constant_values=0.0)
#                 shape = eigvals[0] / (eigvals[2] + 1e-8)
#                 curvature = eigvals[1] / (eigvals[2] + 1e-8)
#             else:
#                 shape = 0.0
#                 curvature = 0.0

#             shape_list.append(np.clip(shape, 0, 1))
#             curvature_list.append(np.clip(curvature, 0, 1))

#             # Electrostatics
#             if ligand_charge > 0:
#                 electro = int(res.name in neg_residues)
#             elif ligand_charge < 0:
#                 electro = int(res.name in pos_residues)
#             else:
#                 electro = int(res.name in pos_residues or res.name in neg_residues)
#             electro_list.append(electro)

#             # Hydrophobicity
#             hydro = int(res.name in hydrophobic_residues)
#             hydro_list.append(hydro)

#             # HBond
#             lig_polars = [i for i in ligand_atoms if t_slice.topology.atom(i).element.symbol in ('N', 'O')]
#             prot_elem = t_slice.topology.atom(prot_idx).element.symbol
#             hbond = 0
#             if prot_elem in ('N', 'O'):
#                 pairs = [(i, prot_idx) for i in lig_polars]
#                 if pairs:
#                     dists = md.compute_distances(t_slice, pairs)[0]
#                     hbond = int(np.any(dists < 0.35))
#             hbond_list.append(hbond)

#         # ---- DataFrame ----
#         df = pd.DataFrame({
#             'r': r,
#             'theta': theta,
#             'ShapeIndex': shape_list,
#             'Curvature': curvature_list,
#             'Electrostatics': electro_list,
#             'HBonds': hbond_list,
#             'Hydrophobicity': hydro_list
#         })
#         df['Frame'] = frame

#         all_embeddings.append(df)

#     if not all_embeddings:
#         raise ValueError("No valid frames found for MaSIF embedding.")

#     embedding_df = pd.concat(all_embeddings, ignore_index=True)
#     embedding_file = os.path.join(graph_dir, "masif_geodesic_embedding.xlsx")
#     embedding_df.to_excel(embedding_file, index=False)
#     print(f"MaSIF geodesic embedding saved to: {embedding_file}")

#     # Plotting first frame all features
#     print("Saving MaSIF embedding feature maps...")
#     _plot_patch_example(all_embeddings[0], masif_dir)

# def _plot_patch_example(df, save_dir):
#     """Polar plots for all features in one patch."""
#     from matplotlib import colormaps
#     cmap = colormaps['turbo']  # 青 → 橙 → 粉红

#     features = ["ShapeIndex", "Curvature", "Electrostatics", "HBonds", "Hydrophobicity"]
#     for feat in features:
#         fig = plt.figure(figsize=(8, 8))
#         ax = plt.subplot(111, polar=True)
#         c = df[feat]
#         sc = ax.scatter(df['theta'], df['r'], c=c, cmap=cmap, s=120)
#         plt.title(f"MaSIF Patch Embedding - {feat}")
#         cbar = plt.colorbar(sc, label=feat)
#         plt.tight_layout()
#         filename = os.path.join(save_dir, f"{feat}_polar_plot.png")
#         plt.savefig(filename, dpi=300)
#         plt.show()
#         plt.close()

######################################
# Version 1.2
######################################

# import os
# import mdtraj as md
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from rdkit import Chem
# from rdkit.Chem import rdmolops
# from tqdm import tqdm

# plt.rcParams.update({'font.size': 21, 'font.family': 'Times New Roman'})

# def compute_masif_geodesic_embedding(topol, traj_file, sdf, resname="LIG", cutoff=0.4, n_frames=10):
#     """
#     Compute MaSIF-like geodesic feature embedding for ligand-protein interface patches.
#     Output: Per-frame geodesic patch embeddings + visualizations for all MaSIF features.
#     """
#     traj = md.load(traj_file, top=topol)
#     if n_frames > traj.n_frames:
#         n_frames = traj.n_frames

#     ligand_atoms = traj.topology.select(f"resname {resname} and element != H")
#     protein_atoms = traj.topology.select("protein and element != H")
#     if len(ligand_atoms) == 0 or len(protein_atoms) == 0:
#         raise ValueError("No ligand or protein heavy atoms found.")

#     mol = Chem.MolFromMolFile(sdf, removeHs=False)
#     if mol is None:
#         raise ValueError("Failed to read ligand SDF.")
#     ligand_charge = rdmolops.GetFormalCharge(mol)

#     pos_residues = {"ARG", "LYS"}
#     neg_residues = {"ASP", "GLU"}

#     # Kyte-Doolittle hydropathy scale
#     hydropathy_scale = {
#         "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5,
#         "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8,
#         "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2, "GLU": -3.5,
#         "GLN": -3.5, "ASP": -3.5, "ASN": -3.5, "LYS": -3.9, "ARG": -4.5
#     }

#     out_dir = os.path.dirname(os.path.abspath(traj_file))
#     graph_dir = os.path.join(out_dir, "Graph")
#     masif_dir = os.path.join(graph_dir, "masif")
#     os.makedirs(masif_dir, exist_ok=True)

#     all_embeddings = []

#     for frame in tqdm(range(n_frames), desc="Computing MaSIF embeddings"):
#         t_slice = traj.slice(frame)

#         lig_coords = t_slice.xyz[0, ligand_atoms, :]
#         lig_center = lig_coords.mean(axis=0)

#         distances = md.compute_distances(t_slice,
#                         [(i, j) for i in ligand_atoms for j in protein_atoms])[0]
#         contact_pairs = [(ligand_atoms[i // len(protein_atoms)],
#                           protein_atoms[i % len(protein_atoms)])
#                           for i, d in enumerate(distances) if d < cutoff]
#         contacting_prot_atoms = set(j for _, j in contact_pairs)

#         if len(contacting_prot_atoms) < 3:
#             continue

#         prot_coords = t_slice.xyz[0, list(contacting_prot_atoms), :]
#         center = prot_coords.mean(axis=0)

#         r = np.linalg.norm(prot_coords - center, axis=1)

#         direction = lig_center - center
#         direction /= (np.linalg.norm(direction) + 1e-8)

#         vectors = prot_coords - center
#         u_axis = vectors[0] / (np.linalg.norm(vectors[0]) + 1e-8)
#         v_axis = np.cross(direction, u_axis)
#         v_axis /= (np.linalg.norm(v_axis) + 1e-8)

#         theta = []
#         for vec in vectors:
#             proj_u = np.dot(vec, u_axis)
#             proj_v = np.dot(vec, v_axis)
#             angle = np.arctan2(proj_v, proj_u) % (2 * np.pi)
#             theta.append(angle)

#         # ===== Features =====
#         contacting_residues = [t_slice.topology.atom(idx).residue for idx in contacting_prot_atoms]
#         pos_count = sum(res.name in pos_residues for res in contacting_residues)
#         neg_count = sum(res.name in neg_residues for res in contacting_residues)

#         if ligand_charge > 0:
#             electro = neg_count
#         elif ligand_charge < 0:
#             electro = pos_count
#         else:
#             electro = pos_count + neg_count

#         lig_polars = [i for i in ligand_atoms if t_slice.topology.atom(i).element.symbol in ('N', 'O')]
#         prot_polars = [j for j in contacting_prot_atoms if t_slice.topology.atom(j).element.symbol in ('N', 'O')]
#         polar_pairs = [(i, j) for i in lig_polars for j in prot_polars]
#         hbond = 0
#         if polar_pairs:
#             dists = md.compute_distances(t_slice, polar_pairs)[0]
#             hbond = np.sum(dists < 0.35)

#         # ===== Hydrophobicity: continuous =====
#         hydro_values = []
#         for res in contacting_residues:
#             hydro_values.append(hydropathy_scale.get(res.name, 0.0))
#         if hydro_values:
#             hydro_score = np.mean(hydro_values)
#         else:
#             hydro_score = 0.0

#         centered = prot_coords - center
#         cov = np.cov(centered.T)
#         eigvals = np.linalg.eigvalsh(cov)
#         eigvals = np.sort(eigvals)
#         if len(eigvals) < 3:
#             eigvals = np.pad(eigvals, (0, 3 - len(eigvals)), constant_values=0.0)
#         shape_idx = eigvals[0] / (eigvals[2] + 1e-8)
#         curvature = eigvals[1] / (eigvals[2] + 1e-8)

#         shape_idx = np.clip(shape_idx, 0, 1)
#         curvature = np.clip(curvature, 0, 1)

#         df = pd.DataFrame({
#             'r': r,
#             'theta': theta,
#             'ShapeIndex': shape_idx,
#             'Curvature': curvature,
#             'Electrostatics': electro,
#             'HBonds': hbond,
#             'Hydrophobicity': hydro_score
#         })
#         df['Frame'] = frame

#         all_embeddings.append(df)

#     if not all_embeddings:
#         raise ValueError("No valid frames found for MaSIF embedding.")

#     embedding_df = pd.concat(all_embeddings, ignore_index=True)
#     embedding_file = os.path.join(graph_dir, "masif_geodesic_embedding.xlsx")
#     embedding_df.to_excel(embedding_file, index=False)
#     print(f"MaSIF geodesic embedding saved to: {embedding_file}")

#     print("Saving MaSIF embedding feature maps...")
#     _plot_patch_example(all_embeddings[0], masif_dir)

# def _plot_patch_example(df, save_dir):
#     """Generate polar plots for all features of one patch."""
#     from matplotlib import cm
#     from matplotlib.colors import LinearSegmentedColormap

#     features = ["ShapeIndex", "Curvature", "Electrostatics", "HBonds", "Hydrophobicity"]

#     # Custom colormap: cyan → orange → pink
#     cmap = LinearSegmentedColormap.from_list("custom", ["#00ced1", "#ffa500", "#ff69b4"])

#     for feat in features:
#         fig = plt.figure(figsize=(8, 8))
#         ax = plt.subplot(111, polar=True)
#         c = df[feat]
#         sc = ax.scatter(df['theta'], df['r'], c=c, cmap=cmap, s=100)
#         plt.title(f"MaSIF Patch Embedding - {feat}")
#         cbar = plt.colorbar(sc, label=feat)
#         plt.tight_layout()
#         filename = os.path.join(save_dir, f"{feat}_polar_plot.png")
#         plt.savefig(filename, dpi=300)
#         plt.show()
#         plt.close()

# ##########################################
# # Version 1.3
# ##########################################
# import os
# import mdtraj as md
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from rdkit import Chem
# from rdkit.Chem import rdmolops
# from tqdm import tqdm

# plt.rcParams.update({'font.size': 21, 'font.family': 'Times New Roman'})

# def compute_masif_geodesic_embedding(topol, traj_file, sdf, resname="LIG", cutoff=0.4, n_frames=10):
#     """
#     Compute MaSIF-like geodesic feature embedding for ligand-protein interface patches.
#     Output: Per-frame geodesic patch embeddings + visualizations.
#     """
#     traj = md.load(traj_file, top=topol)
#     if n_frames > traj.n_frames:
#         n_frames = traj.n_frames

#     ligand_atoms = traj.topology.select(f"resname {resname} and element != H")
#     protein_atoms = traj.topology.select("protein and element != H")

#     if len(ligand_atoms) == 0 or len(protein_atoms) == 0:
#         raise ValueError("No ligand or protein heavy atoms found.")

#     mol = Chem.MolFromMolFile(sdf, removeHs=False)
#     if mol is None:
#         raise ValueError("Failed to read ligand SDF.")
#     ligand_charge = rdmolops.GetFormalCharge(mol)

#     pos_residues = {"ARG", "LYS"}
#     neg_residues = {"ASP", "GLU"}

#     # Kyte-Doolittle scale
#     kd_scale = {
#         "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5,
#         "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8,
#         "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2,
#         "GLU": -3.5, "GLN": -3.5, "ASP": -3.5, "ASN": -3.5,
#         "LYS": -3.9, "ARG": -4.5
#     }

#     out_dir = os.path.dirname(os.path.abspath(traj_file))
#     graph_dir = os.path.join(out_dir, "Graph")
#     masif_dir = os.path.join(graph_dir, "masif")
#     os.makedirs(masif_dir, exist_ok=True)

#     all_embeddings = []

#     for frame in tqdm(range(n_frames), desc="Computing MaSIF embeddings"):
#         t_slice = traj.slice(frame)
#         lig_coords = t_slice.xyz[0, ligand_atoms, :]
#         lig_center = lig_coords.mean(axis=0)

#         distances = md.compute_distances(t_slice,
#                         [(i, j) for i in ligand_atoms for j in protein_atoms])[0]
#         contact_pairs = [(ligand_atoms[i // len(protein_atoms)],
#                           protein_atoms[i % len(protein_atoms)])
#                           for i, d in enumerate(distances) if d < cutoff]
#         contacting_prot_atoms = list(set(j for _, j in contact_pairs))

#         if len(contacting_prot_atoms) < 3:
#             continue

#         prot_coords = t_slice.xyz[0, contacting_prot_atoms, :]
#         center = prot_coords.mean(axis=0)

#         r = np.linalg.norm(prot_coords - center, axis=1)

#         direction = lig_center - center
#         direction /= (np.linalg.norm(direction) + 1e-8)

#         vectors = prot_coords - center
#         u_axis = vectors[0] / (np.linalg.norm(vectors[0]) + 1e-8)
#         v_axis = np.cross(direction, u_axis)
#         v_axis /= (np.linalg.norm(v_axis) + 1e-8)

#         theta = []
#         for vec in vectors:
#             proj_u = np.dot(vec, u_axis)
#             proj_v = np.dot(vec, v_axis)
#             angle = np.arctan2(proj_v, proj_u) % (2 * np.pi)
#             theta.append(angle)

#         # Per-atom features
#         shape_list = []
#         curvature_list = []
#         electro_list = []
#         hbond_list = []
#         hydro_list = []

#         for idx, prot_idx in enumerate(contacting_prot_atoms):
#             res = t_slice.topology.atom(prot_idx).residue

#             # --- Local shape & curvature ---
#             neighbor_coords = []
#             for jdx, other_idx in enumerate(contacting_prot_atoms):
#                 d = np.linalg.norm(prot_coords[idx] - prot_coords[jdx])
#                 if d < 0.5:
#                     neighbor_coords.append(prot_coords[jdx])
#             neighbor_coords = np.array(neighbor_coords)

#             if len(neighbor_coords) >= 3:
#                 centered = neighbor_coords - neighbor_coords.mean(axis=0)
#                 cov = np.cov(centered.T)
#                 eigvals = np.linalg.eigvalsh(cov)
#                 eigvals = np.sort(eigvals)
#                 if len(eigvals) < 3:
#                     eigvals = np.pad(eigvals, (0, 3 - len(eigvals)), constant_values=0.0)
#                 shape = eigvals[0] / (eigvals[2] + 1e-8)
#                 curvature = eigvals[1] / (eigvals[2] + 1e-8)
#             else:
#                 shape = 0.0
#                 curvature = 0.0

#             shape_list.append(np.clip(shape, 0, 1))
#             curvature_list.append(np.clip(curvature, 0, 1))

#             # --- Electrostatics ---
#             if ligand_charge > 0:
#                 electro = int(res.name in neg_residues)
#             elif ligand_charge < 0:
#                 electro = int(res.name in pos_residues)
#             else:
#                 electro = int(res.name in pos_residues or res.name in neg_residues)
#             electro_list.append(electro)

#             # --- Hydrophobicity ---
#             hydro = kd_scale.get(res.name, 0.0)
#             hydro_list.append(hydro)

#             # --- HBond ---
#             lig_polars = [i for i in ligand_atoms if t_slice.topology.atom(i).element.symbol in ('N', 'O')]
#             prot_elem = t_slice.topology.atom(prot_idx).element.symbol
#             hbond = 0
#             if prot_elem in ('N', 'O'):
#                 pairs = [(i, prot_idx) for i in lig_polars]
#                 if pairs:
#                     dists = md.compute_distances(t_slice, pairs)[0]
#                     hbond = int(np.any(dists < 0.35))
#             hbond_list.append(hbond)

#         df = pd.DataFrame({
#             'r': r,
#             'theta': theta,
#             'ShapeIndex': shape_list,
#             'Curvature': curvature_list,
#             'Electrostatics': electro_list,
#             'HBonds': hbond_list,
#             'Hydrophobicity': hydro_list
#         })
#         df['Frame'] = frame
#         all_embeddings.append(df)

#     if not all_embeddings:
#         raise ValueError("No valid frames found for MaSIF embedding.")

#     embedding_df = pd.concat(all_embeddings, ignore_index=True)
#     embedding_file = os.path.join(graph_dir, "masif_geodesic_embedding.xlsx")
#     embedding_df.to_excel(embedding_file, index=False)
#     print(f"MaSIF geodesic embedding saved to: {embedding_file}")

#     # Plotting first frame for all features
#     print("Saving MaSIF embedding feature maps...")
#     _plot_patch_example(all_embeddings[0], masif_dir)

# def _plot_patch_example(df, save_dir):
#     """Polar plots for all features in one patch."""
#     from matplotlib import colormaps
#     from matplotlib.colors import LinearSegmentedColormap
#     # cmap = colormaps['turbo']  # 青 → 橙 → 粉红
#     cmap = LinearSegmentedColormap.from_list("custom", ["#00ced1", "#ffa500", "#ff69b4"])

#     features = ["ShapeIndex", "Curvature", "Electrostatics", "HBonds", "Hydrophobicity"]
#     for feat in features:
#         fig = plt.figure(figsize=(8, 8))
#         ax = plt.subplot(111, polar=True)
#         c = df[feat]
#         sc = ax.scatter(df['theta'], df['r'], c=c, cmap=cmap, s=120)
#         plt.title(f"MaSIF Patch Embedding - {feat}")
#         cbar = plt.colorbar(sc, label=feat)
#         plt.tight_layout()
#         filename = os.path.join(save_dir, f"{feat}_polar_plot.png")
#         plt.savefig(filename, dpi=300)
#         plt.show()
#         plt.close()

##################################
# Version 1.4
##################################

# import os
# import mdtraj as md
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from rdkit import Chem
# from rdkit.Chem import rdmolops
# from tqdm import tqdm
# import networkx as nx

# plt.rcParams.update({'font.size': 21, 'font.family': 'Times New Roman'})

# def compute_masif_geodesic_embedding(topol, traj_file, sdf, resname="LIG", cutoff=0.4, n_frames=10, distance_mode="euclidean"):
#     """
#     Compute MaSIF-like geodesic feature embedding for ligand-protein interface patches.
#     distance_mode: 'euclidean' or 'geodesic'
#     """
#     traj = md.load(traj_file, top=topol)
#     if n_frames > traj.n_frames:
#         n_frames = traj.n_frames

#     ligand_atoms = traj.topology.select(f"resname {resname} and element != H")
#     protein_atoms = traj.topology.select("protein and element != H")

#     if len(ligand_atoms) == 0 or len(protein_atoms) == 0:
#         raise ValueError("No ligand or protein heavy atoms found.")

#     mol = Chem.MolFromMolFile(sdf, removeHs=False)
#     if mol is None:
#         raise ValueError("Failed to read ligand SDF.")
#     ligand_charge = rdmolops.GetFormalCharge(mol)

#     pos_residues = {"ARG", "LYS"}
#     neg_residues = {"ASP", "GLU"}

#     kd_scale = {
#         "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5,
#         "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8,
#         "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2,
#         "GLU": -3.5, "GLN": -3.5, "ASP": -3.5, "ASN": -3.5,
#         "LYS": -3.9, "ARG": -4.5
#     }

#     out_dir = os.path.dirname(os.path.abspath(traj_file))
#     graph_dir = os.path.join(out_dir, "Graph")
#     masif_dir = os.path.join(graph_dir, "masif")
#     os.makedirs(masif_dir, exist_ok=True)

#     all_embeddings = []

#     for frame in tqdm(range(n_frames), desc="Computing MaSIF embeddings"):
#         t_slice = traj.slice(frame)
#         lig_coords = t_slice.xyz[0, ligand_atoms, :]
#         lig_center = lig_coords.mean(axis=0)

#         distances = md.compute_distances(t_slice,
#                         [(i, j) for i in ligand_atoms for j in protein_atoms])[0]
#         contact_pairs = [(ligand_atoms[i // len(protein_atoms)],
#                           protein_atoms[i % len(protein_atoms)])
#                           for i, d in enumerate(distances) if d < cutoff]
#         contacting_prot_atoms = list(set(j for _, j in contact_pairs))

#         if len(contacting_prot_atoms) < 3:
#             continue

#         prot_coords = t_slice.xyz[0, contacting_prot_atoms, :]
#         center = prot_coords.mean(axis=0)

#         ### ---- Compute r ---- ###
#         if distance_mode == "euclidean":
#             r = np.linalg.norm(prot_coords - center, axis=1)

#         elif distance_mode == "geodesic":
#             # Build protein-protein graph
#             G = nx.Graph()
#             for i, idx_i in enumerate(contacting_prot_atoms):
#                 G.add_node(idx_i, pos=prot_coords[i])

#             for i, idx_i in enumerate(contacting_prot_atoms):
#                 for j, idx_j in enumerate(contacting_prot_atoms):
#                     if i < j:
#                         d = np.linalg.norm(prot_coords[i] - prot_coords[j])
#                         if d < 0.8:
#                             G.add_edge(idx_i, idx_j, weight=d)

#             # Find the protein atom closest to lig_center
#             distances_to_lig = np.linalg.norm(prot_coords - lig_center, axis=1)
#             lig_like_prot_idx = contacting_prot_atoms[np.argmin(distances_to_lig)]

#             # Compute geodesic distances
#             r = []
#             for idx in contacting_prot_atoms:
#                 if nx.has_path(G, idx, lig_like_prot_idx):
#                     path_length = nx.dijkstra_path_length(G, idx, lig_like_prot_idx, weight='weight')
#                 else:
#                     path_length = 999.0
#                 r.append(path_length)
#             r = np.array(r)
#         else:
#             raise ValueError("distance_mode must be 'euclidean' or 'geodesic'")

#         ### ---- Compute theta ---- ###
#         direction = lig_center - center
#         direction /= (np.linalg.norm(direction) + 1e-8)
#         vectors = prot_coords - center
#         u_axis = vectors[0] / (np.linalg.norm(vectors[0]) + 1e-8)
#         v_axis = np.cross(direction, u_axis)
#         v_axis /= (np.linalg.norm(v_axis) + 1e-8)

#         theta = []
#         for vec in vectors:
#             proj_u = np.dot(vec, u_axis)
#             proj_v = np.dot(vec, v_axis)
#             angle = np.arctan2(proj_v, proj_u) % (2 * np.pi)
#             theta.append(angle)

#         # --- Per-atom features ---
#         shape_list = []
#         curvature_list = []
#         electro_list = []
#         hbond_list = []
#         hydro_list = []

#         for idx, prot_idx in enumerate(contacting_prot_atoms):
#             res = t_slice.topology.atom(prot_idx).residue

#             # --- Local shape & curvature ---
#             neighbor_coords = []
#             for jdx, other_idx in enumerate(contacting_prot_atoms):
#                 d = np.linalg.norm(prot_coords[idx] - prot_coords[jdx])
#                 if d < 0.5:
#                     neighbor_coords.append(prot_coords[jdx])
#             neighbor_coords = np.array(neighbor_coords)

#             if len(neighbor_coords) >= 3:
#                 centered = neighbor_coords - neighbor_coords.mean(axis=0)
#                 cov = np.cov(centered.T)
#                 eigvals = np.linalg.eigvalsh(cov)
#                 eigvals = np.sort(eigvals)
#                 if len(eigvals) < 3:
#                     eigvals = np.pad(eigvals, (0, 3 - len(eigvals)), constant_values=0.0)
#                 shape = eigvals[0] / (eigvals[2] + 1e-8)
#                 curvature = eigvals[1] / (eigvals[2] + 1e-8)
#             else:
#                 shape = 0.0
#                 curvature = 0.0

#             shape_list.append(np.clip(shape, 0, 1))
#             curvature_list.append(np.clip(curvature, 0, 1))

#             # --- Electrostatics ---
#             if ligand_charge > 0:
#                 electro = int(res.name in neg_residues)
#             elif ligand_charge < 0:
#                 electro = int(res.name in pos_residues)
#             else:
#                 electro = int(res.name in pos_residues or res.name in neg_residues)
#             electro_list.append(electro)

#             # --- Hydrophobicity ---
#             hydro = kd_scale.get(res.name, 0.0)
#             hydro_list.append(hydro)

#             # --- HBond ---
#             lig_polars = [i for i in ligand_atoms if t_slice.topology.atom(i).element.symbol in ('N', 'O')]
#             prot_elem = t_slice.topology.atom(prot_idx).element.symbol
#             hbond = 0
#             if prot_elem in ('N', 'O'):
#                 pairs = [(i, prot_idx) for i in lig_polars]
#                 if pairs:
#                     dists = md.compute_distances(t_slice, pairs)[0]
#                     hbond = int(np.any(dists < 0.35))
#             hbond_list.append(hbond)

#         df = pd.DataFrame({
#             'r': r,
#             'theta': theta,
#             'ShapeIndex': shape_list,
#             'Curvature': curvature_list,
#             'Electrostatics': electro_list,
#             'HBonds': hbond_list,
#             'Hydrophobicity': hydro_list
#         })
#         df['Frame'] = frame
#         all_embeddings.append(df)

#     if not all_embeddings:
#         raise ValueError("No valid frames found for MaSIF embedding.")

#     embedding_df = pd.concat(all_embeddings, ignore_index=True)
#     embedding_file = os.path.join(graph_dir, f"masif_geodesic_embedding_{distance_mode}.xlsx")
#     embedding_df.to_excel(embedding_file, index=False)
#     print(f"MaSIF geodesic embedding saved to: {embedding_file}")

#     print("Saving MaSIF embedding feature maps...")
#     _plot_patch_example(all_embeddings[0], masif_dir, distance_mode)

# def _plot_patch_example(df, save_dir, distance_mode):
#     from matplotlib.colors import LinearSegmentedColormap
#     cmap = LinearSegmentedColormap.from_list("custom", ["#00ced1", "#ffa500", "#ff69b4"])

#     features = ["ShapeIndex", "Curvature", "Electrostatics", "HBonds", "Hydrophobicity"]
#     for feat in features:
#         fig = plt.figure(figsize=(8, 8))
#         ax = plt.subplot(111, polar=True)
#         c = df[feat]
#         sc = ax.scatter(df['theta'], df['r'], c=c, cmap=cmap, s=120)
#         plt.title(f"MaSIF - {feat} ({distance_mode})")
#         cbar = plt.colorbar(sc, label=feat)
#         plt.tight_layout()
#         filename = os.path.join(save_dir, f"{feat}_polar_plot_{distance_mode}.png")
#         plt.savefig(filename, dpi=300)
#         plt.show()
#         plt.close()

# ############################
# # Version 1.5
# ############################
# import os
# import mdtraj as md
# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from rdkit import Chem
# from rdkit.Chem import rdmolops
# from tqdm import tqdm
# import networkx as nx

# plt.rcParams.update({'font.size': 21, 'font.family': 'Times New Roman'})

# def compute_masif_geodesic_embedding(topol, traj_file, sdf, resname="LIG", cutoff=0.4, n_frames=10, distance_mode="euclidean"):
#     """
#     Compute MaSIF-like embedding for ligand-protein interface patches.
#     For euclidean distance, r is the distance from the pocket center to each contacting protein atom.
#     For geodesic distance, r is the shortest path from the closest ligand atom to each contacting protein atom
#     across a graph that includes ligand and protein contact atoms.
#     """
#     traj = md.load(traj_file, top=topol)
#     if n_frames > traj.n_frames:
#         n_frames = traj.n_frames

#     # Select ligand and protein heavy atoms
#     ligand_atoms = traj.topology.select(f"resname {resname} and element != H")
#     protein_atoms = traj.topology.select("protein and element != H")

#     if len(ligand_atoms) == 0 or len(protein_atoms) == 0:
#         raise ValueError("No ligand or protein heavy atoms found.")

#     # Read ligand SDF to get formal charge
#     mol = Chem.MolFromMolFile(sdf, removeHs=False)
#     if mol is None:
#         raise ValueError("Failed to read ligand SDF.")
#     ligand_charge = rdmolops.GetFormalCharge(mol)

#     # Define residue properties
#     pos_residues = {"ARG", "LYS"}
#     neg_residues = {"ASP", "GLU"}

#     kd_scale = {
#         "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5,
#         "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8,
#         "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2,
#         "GLU": -3.5, "GLN": -3.5, "ASP": -3.5, "ASN": -3.5,
#         "LYS": -3.9, "ARG": -4.5
#     }

#     out_dir = os.path.dirname(os.path.abspath(traj_file))
#     graph_dir = os.path.join(out_dir, "Graph")
#     masif_dir = os.path.join(graph_dir, "masif")
#     os.makedirs(masif_dir, exist_ok=True)

#     all_embeddings = []

#     for frame in tqdm(range(n_frames), desc="Computing MaSIF embeddings"):
#         t_slice = traj.slice(frame)
#         lig_coords = t_slice.xyz[0, ligand_atoms, :]
#         lig_center = lig_coords.mean(axis=0)

#         # Find contacting protein atoms within cutoff distance
#         distances = md.compute_distances(t_slice,
#                         [(i, j) for i in ligand_atoms for j in protein_atoms])[0]
#         contact_pairs = [(ligand_atoms[i // len(protein_atoms)],
#                           protein_atoms[i % len(protein_atoms)])
#                           for i, d in enumerate(distances) if d < cutoff]
#         contacting_prot_atoms = list(set(j for _, j in contact_pairs))

#         if len(contacting_prot_atoms) < 3:
#             continue

#         prot_coords = t_slice.xyz[0, contacting_prot_atoms, :]
#         center = prot_coords.mean(axis=0)

#         ### ---- Compute r ---- ###
#         if distance_mode == "euclidean":
#             # r = Euclidean distance from pocket center to each contacting protein atom
#             r = np.linalg.norm(prot_coords - center, axis=1)

#         elif distance_mode == "geodesic":

#             ### --- 1. Build graph including ligand atoms and contacting protein atoms --- ###
#             G = nx.Graph()

#             # Add ligand atoms as nodes
#             for i, lig_idx in enumerate(ligand_atoms):
#                 G.add_node(lig_idx, pos=t_slice.xyz[0, lig_idx, :])

#             # Add contacting protein atoms as nodes
#             for i, prot_idx in enumerate(contacting_prot_atoms):
#                 G.add_node(prot_idx, pos=t_slice.xyz[0, prot_idx, :])

#             # Add edges between ligand atoms (internal ligand connectivity)
#             for i, idx_i in enumerate(ligand_atoms):
#                 for j, idx_j in enumerate(ligand_atoms):
#                     if i < j:
#                         d = np.linalg.norm(t_slice.xyz[0, idx_i, :] - t_slice.xyz[0, idx_j, :])
#                         if d < 0.5:
#                             G.add_edge(idx_i, idx_j, weight=d)

#             # Add edges between contacting protein atoms
#             for i, idx_i in enumerate(contacting_prot_atoms):
#                 for j, idx_j in enumerate(contacting_prot_atoms):
#                     if i < j:
#                         d = np.linalg.norm(t_slice.xyz[0, idx_i, :] - t_slice.xyz[0, idx_j, :])
#                         if d < 0.8:
#                             G.add_edge(idx_i, idx_j, weight=d)

#             # Add edges between ligand and protein atoms
#             for lig_idx in ligand_atoms:
#                 for prot_idx in contacting_prot_atoms:
#                     d = np.linalg.norm(t_slice.xyz[0, lig_idx, :] - t_slice.xyz[0, prot_idx, :])
#                     if d < 0.5:
#                         G.add_edge(lig_idx, prot_idx, weight=d)

#             ### --- 2. Find ligand atom closest to the protein pocket center --- ###
#             lig_dists = np.linalg.norm(lig_coords - center, axis=1)
#             lig_start_atom = ligand_atoms[np.argmin(lig_dists)]

#             ### --- 3. Compute geodesic distances from ligand atom to each protein contact atom --- ###
#             r = []
#             for prot_idx in contacting_prot_atoms:
#                 try:
#                     path_length = nx.dijkstra_path_length(G, lig_start_atom, prot_idx, weight='weight')
#                 except nx.NetworkXNoPath:
#                     path_length = 999.0  # If no path exists, assign a large distance
#                 r.append(path_length)
#             r = np.array(r)

#         else:
#             raise ValueError("distance_mode must be 'euclidean' or 'geodesic'")

#         ### ---- Compute theta ---- ###
#         # Direction vector from pocket center to ligand center
#         direction = lig_center - center
#         direction /= (np.linalg.norm(direction) + 1e-8)
#         vectors = prot_coords - center

#         # Define reference axes in the plane
#         u_axis = vectors[0] / (np.linalg.norm(vectors[0]) + 1e-8)
#         v_axis = np.cross(direction, u_axis)
#         v_axis /= (np.linalg.norm(v_axis) + 1e-8)

#         theta = []
#         for vec in vectors:
#             proj_u = np.dot(vec, u_axis)
#             proj_v = np.dot(vec, v_axis)
#             angle = np.arctan2(proj_v, proj_u) % (2 * np.pi)
#             theta.append(angle)

#         # --- Compute per-atom features ---
#         shape_list = []
#         curvature_list = []
#         electro_list = []
#         hbond_list = []
#         hydro_list = []

#         for idx, prot_idx in enumerate(contacting_prot_atoms):
#             res = t_slice.topology.atom(prot_idx).residue

#             # Local shape & curvature (PCA eigenvalues of neighboring atoms)
#             neighbor_coords = []
#             for jdx, other_idx in enumerate(contacting_prot_atoms):
#                 d = np.linalg.norm(prot_coords[idx] - prot_coords[jdx])
#                 if d < 0.5:
#                     neighbor_coords.append(prot_coords[jdx])
#             neighbor_coords = np.array(neighbor_coords)

#             if len(neighbor_coords) >= 3:
#                 centered = neighbor_coords - neighbor_coords.mean(axis=0)
#                 cov = np.cov(centered.T)
#                 eigvals = np.linalg.eigvalsh(cov)
#                 eigvals = np.sort(eigvals)
#                 if len(eigvals) < 3:
#                     eigvals = np.pad(eigvals, (0, 3 - len(eigvals)), constant_values=0.0)
#                 shape = eigvals[0] / (eigvals[2] + 1e-8)
#                 curvature = eigvals[1] / (eigvals[2] + 1e-8)
#             else:
#                 shape = 0.0
#                 curvature = 0.0

#             shape_list.append(np.clip(shape, 0, 1))
#             curvature_list.append(np.clip(curvature, 0, 1))

#             # Electrostatics: basic rule based on residue and ligand charge
#             if ligand_charge > 0:
#                 electro = int(res.name in neg_residues)
#             elif ligand_charge < 0:
#                 electro = int(res.name in pos_residues)
#             else:
#                 electro = int(res.name in pos_residues or res.name in neg_residues)
#             electro_list.append(electro)

#             # Hydrophobicity (Kyte-Doolittle scale)
#             hydro = kd_scale.get(res.name, 0.0)
#             hydro_list.append(hydro)

#             # HBond potential: check if any ligand polar atoms are close to the protein atom
#             lig_polars = [i for i in ligand_atoms if t_slice.topology.atom(i).element.symbol in ('N', 'O')]
#             prot_elem = t_slice.topology.atom(prot_idx).element.symbol
#             hbond = 0
#             if prot_elem in ('N', 'O'):
#                 pairs = [(i, prot_idx) for i in lig_polars]
#                 if pairs:
#                     dists = md.compute_distances(t_slice, pairs)[0]
#                     hbond = int(np.any(dists < 0.35))
#             hbond_list.append(hbond)

#         df = pd.DataFrame({
#             'r': r,
#             'theta': theta,
#             'ShapeIndex': shape_list,
#             'Curvature': curvature_list,
#             'Electrostatics': electro_list,
#             'HBonds': hbond_list,
#             'Hydrophobicity': hydro_list
#         })
#         df['Frame'] = frame
#         all_embeddings.append(df)

#     if not all_embeddings:
#         raise ValueError("No valid frames found for MaSIF embedding.")

#     # Save results
#     embedding_df = pd.concat(all_embeddings, ignore_index=True)
#     embedding_file = os.path.join(graph_dir, f"masif_geodesic_embedding_{distance_mode}.xlsx")
#     embedding_df.to_excel(embedding_file, index=False)
#     print(f"MaSIF geodesic embedding saved to: {embedding_file}")

#     # Plot feature maps for the first frame
#     print("Saving MaSIF embedding feature maps...")
#     _plot_patch_example(all_embeddings[0], masif_dir, distance_mode)

# def _plot_patch_example(df, save_dir, distance_mode):
#     """
#     Plot polar feature maps for the first frame.
#     Each point represents a contacting protein atom.
#     """
#     from matplotlib.colors import LinearSegmentedColormap
#     cmap = LinearSegmentedColormap.from_list("custom", ["#00ced1", "#ffa500", "#ff69b4"])

#     features = ["ShapeIndex", "Curvature", "Electrostatics", "HBonds", "Hydrophobicity"]
#     for feat in features:
#         fig = plt.figure(figsize=(8, 8))
#         ax = plt.subplot(111, polar=True)
#         c = df[feat]
#         sc = ax.scatter(df['theta'], df['r'], c=c, cmap=cmap, s=120)
#         plt.title(f"MaSIF - {feat} ({distance_mode})")
#         cbar = plt.colorbar(sc, label=feat)
#         plt.tight_layout()
#         filename = os.path.join(save_dir, f"{feat}_polar_plot_{distance_mode}.png")
#         plt.savefig(filename, dpi=300)
#         plt.show()
#         plt.close()

####################################
# Version 1.6: LigAtom(nearest to pocket center) --> LigAtom(nearest to contact atom) --> contact atom
####################################
import os
import mdtraj as md
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from rdkit import Chem
from rdkit.Chem import rdmolops
from tqdm import tqdm
import networkx as nx

plt.rcParams.update({'font.size': 21, 'font.family': 'Times New Roman'})

def compute_masif_geodesic_embedding(topol, traj_file, sdf, resname="LIG", cutoff=0.4, n_frames=10, distance_mode="euclidean"):
    """
    Compute MaSIF-like embedding for ligand-protein interface patches.
    For euclidean distance: r is distance from pocket center to each contacting protein atom.
    For geodesic distance: r is total path length from ligand start atom to ligand-contact atom to protein atom.
    """
    traj = md.load(traj_file, top=topol)
    if n_frames > traj.n_frames:
        n_frames = traj.n_frames

    ligand_atoms = traj.topology.select(f"resname {resname} and element != H")
    protein_atoms = traj.topology.select("protein and element != H")

    if len(ligand_atoms) == 0 or len(protein_atoms) == 0:
        raise ValueError("No ligand or protein heavy atoms found.")

    mol = Chem.MolFromMolFile(sdf, removeHs=False)
    if mol is None:
        raise ValueError("Failed to read ligand SDF.")
    ligand_charge = rdmolops.GetFormalCharge(mol)

    pos_residues = {"ARG", "LYS"}
    neg_residues = {"ASP", "GLU"}

    kd_scale = {
        "ILE": 4.5, "VAL": 4.2, "LEU": 3.8, "PHE": 2.8, "CYS": 2.5,
        "MET": 1.9, "ALA": 1.8, "GLY": -0.4, "THR": -0.7, "SER": -0.8,
        "TRP": -0.9, "TYR": -1.3, "PRO": -1.6, "HIS": -3.2,
        "GLU": -3.5, "GLN": -3.5, "ASP": -3.5, "ASN": -3.5,
        "LYS": -3.9, "ARG": -4.5
    }

    out_dir = os.path.dirname(os.path.abspath(traj_file))
    graph_dir = os.path.join(out_dir, "Graph")
    masif_dir = os.path.join(graph_dir, "masif")
    os.makedirs(masif_dir, exist_ok=True)

    all_embeddings = []

    for frame in tqdm(range(n_frames), desc="Computing MaSIF embeddings"):
        t_slice = traj.slice(frame)
        lig_coords = t_slice.xyz[0, ligand_atoms, :]
        lig_center = lig_coords.mean(axis=0)

        # Compute all ligand-protein distances
        distances = md.compute_distances(t_slice,
                        [(i, j) for i in ligand_atoms for j in protein_atoms])[0]

        # Identify contact pairs within cutoff
        contact_pairs = [(ligand_atoms[i // len(protein_atoms)],
                          protein_atoms[i % len(protein_atoms)])
                          for i, d in enumerate(distances) if d < cutoff]
        contacting_prot_atoms = list(set(j for _, j in contact_pairs))

        if len(contacting_prot_atoms) < 3:
            continue

        prot_coords = t_slice.xyz[0, contacting_prot_atoms, :]
        center = prot_coords.mean(axis=0)

        ### ---- Compute r ---- ###
        if distance_mode == "euclidean":
            # Simple Euclidean distance from pocket center to each contact protein atom
            r = np.linalg.norm(prot_coords - center, axis=1)

        elif distance_mode == "geodesic":

            ### --- 1. Build Graph: ligand atoms + contacting protein atoms --- ###
            G = nx.Graph()

            # Add ligand atoms
            for i, lig_idx in enumerate(ligand_atoms):
                G.add_node(lig_idx, pos=t_slice.xyz[0, lig_idx, :])

            # Add protein contact atoms
            for i, prot_idx in enumerate(contacting_prot_atoms):
                G.add_node(prot_idx, pos=t_slice.xyz[0, prot_idx, :])

            # Ligand-ligand edges
            for i, idx_i in enumerate(ligand_atoms):
                for j, idx_j in enumerate(ligand_atoms):
                    if i < j:
                        d = np.linalg.norm(t_slice.xyz[0, idx_i, :] - t_slice.xyz[0, idx_j, :])
                        if d < 0.5:
                            G.add_edge(idx_i, idx_j, weight=d)

            # Protein-protein edges
            for i, idx_i in enumerate(contacting_prot_atoms):
                for j, idx_j in enumerate(contacting_prot_atoms):
                    if i < j:
                        d = np.linalg.norm(t_slice.xyz[0, idx_i, :] - t_slice.xyz[0, idx_j, :])
                        if d < 0.8:
                            G.add_edge(idx_i, idx_j, weight=d)

            # Ligand-protein contact edges
            for lig_idx in ligand_atoms:
                for prot_idx in contacting_prot_atoms:
                    d = np.linalg.norm(t_slice.xyz[0, lig_idx, :] - t_slice.xyz[0, prot_idx, :])
                    if d < 0.5:
                        G.add_edge(lig_idx, prot_idx, weight=d)

            ### --- 2. Global ligand start atom: closest to pocket center --- ###
            lig_dists = np.linalg.norm(lig_coords - center, axis=1)
            lig_start_atom = ligand_atoms[np.argmin(lig_dists)]

            ### --- 3. Compute r for each contact protein atom --- ###
            r = []

            for prot_idx in contacting_prot_atoms:
                # Find ligand atoms directly contacting this protein atom
                contact_ligs = []
                for lig_idx in ligand_atoms:
                    d = np.linalg.norm(t_slice.xyz[0, lig_idx, :] - t_slice.xyz[0, prot_idx, :])
                    if d < cutoff:
                        contact_ligs.append((lig_idx, d))

                if len(contact_ligs) == 0:
                    r.append(999.0)
                    continue

                # Pick the closest contacting ligand atom
                lig_contact_atom = sorted(contact_ligs, key=lambda x: x[1])[0][0]

                # Compute two-path geodesic:
                # 1. lig_start_atom -> lig_contact_atom
                # 2. lig_contact_atom -> prot_idx

                try:
                    d1 = nx.dijkstra_path_length(G, lig_start_atom, lig_contact_atom, weight='weight')
                    d2 = nx.dijkstra_path_length(G, lig_contact_atom, prot_idx, weight='weight')
                    total_distance = d1 + d2
                except nx.NetworkXNoPath:
                    total_distance = 999.0

                r.append(total_distance)

            r = np.array(r)

        else:
            raise ValueError("distance_mode must be 'euclidean' or 'geodesic'")

        ### ---- Compute theta ---- ###
        direction = lig_center - center
        direction /= (np.linalg.norm(direction) + 1e-8)
        vectors = prot_coords - center
        u_axis = vectors[0] / (np.linalg.norm(vectors[0]) + 1e-8)
        v_axis = np.cross(direction, u_axis)
        v_axis /= (np.linalg.norm(v_axis) + 1e-8)

        theta = []
        for vec in vectors:
            proj_u = np.dot(vec, u_axis)
            proj_v = np.dot(vec, v_axis)
            angle = np.arctan2(proj_v, proj_u) % (2 * np.pi)
            theta.append(angle)

        ### ---- Compute per-atom features ---- ###
        shape_list = []
        curvature_list = []
        electro_list = []
        hbond_list = []
        hydro_list = []

        for idx, prot_idx in enumerate(contacting_prot_atoms):
            res = t_slice.topology.atom(prot_idx).residue

            # Local shape & curvature (PCA eigenvalues)
            neighbor_coords = []
            for jdx, other_idx in enumerate(contacting_prot_atoms):
                d = np.linalg.norm(prot_coords[idx] - prot_coords[jdx])
                if d < 0.5:
                    neighbor_coords.append(prot_coords[jdx])
            neighbor_coords = np.array(neighbor_coords)

            if len(neighbor_coords) >= 3:
                centered = neighbor_coords - neighbor_coords.mean(axis=0)
                cov = np.cov(centered.T)
                eigvals = np.linalg.eigvalsh(cov)
                eigvals = np.sort(eigvals)
                if len(eigvals) < 3:
                    eigvals = np.pad(eigvals, (0, 3 - len(eigvals)), constant_values=0.0)
                shape = eigvals[0] / (eigvals[2] + 1e-8)
                curvature = eigvals[1] / (eigvals[2] + 1e-8)
            else:
                shape = 0.0
                curvature = 0.0

            shape_list.append(np.clip(shape, 0, 1))
            curvature_list.append(np.clip(curvature, 0, 1))

            # Electrostatics
            if ligand_charge > 0:
                electro = int(res.name in neg_residues)
            elif ligand_charge < 0:
                electro = int(res.name in pos_residues)
            else:
                electro = int(res.name in pos_residues or res.name in neg_residues)
            electro_list.append(electro)

            # Hydrophobicity
            hydro = kd_scale.get(res.name, 0.0)
            hydro_list.append(hydro)

            # Hydrogen bond potential
            lig_polars = [i for i in ligand_atoms if t_slice.topology.atom(i).element.symbol in ('N', 'O')]
            prot_elem = t_slice.topology.atom(prot_idx).element.symbol
            hbond = 0
            if prot_elem in ('N', 'O'):
                pairs = [(i, prot_idx) for i in lig_polars]
                if pairs:
                    dists = md.compute_distances(t_slice, pairs)[0]
                    hbond = int(np.any(dists < 0.35))
            hbond_list.append(hbond)

        df = pd.DataFrame({
            'r': r,
            'theta': theta,
            'ShapeIndex': shape_list,
            'Curvature': curvature_list,
            'Electrostatics': electro_list,
            'HBonds': hbond_list,
            'Hydrophobicity': hydro_list
        })
        df['Frame'] = frame
        all_embeddings.append(df)

    if not all_embeddings:
        raise ValueError("No valid frames found for MaSIF embedding.")

    embedding_df = pd.concat(all_embeddings, ignore_index=True)
    embedding_file = os.path.join(graph_dir, f"masif_geodesic_embedding_{distance_mode}.xlsx")
    embedding_df.to_excel(embedding_file, index=False)
    print(f"MaSIF geodesic embedding saved to: {embedding_file}")

    print("Saving MaSIF embedding feature maps...")
    _plot_patch_example(all_embeddings[0], masif_dir, distance_mode)

def _plot_patch_example(df, save_dir, distance_mode):
    """
    Plot polar feature maps for the first frame.
    Each point represents a contacting protein atom.
    """
    from matplotlib.colors import LinearSegmentedColormap
    cmap = LinearSegmentedColormap.from_list("custom", ["#00ced1", "#ffa500", "#ff69b4"])

    features = ["ShapeIndex", "Curvature", "Electrostatics", "HBonds", "Hydrophobicity"]
    for feat in features:
        fig = plt.figure(figsize=(8, 8))
        ax = plt.subplot(111, polar=True)
        c = df[feat]
        sc = ax.scatter(df['theta'], df['r'], c=c, cmap=cmap, s=120)
        plt.title(f"MaSIF - {feat} ({distance_mode})")
        cbar = plt.colorbar(sc, label=feat)
        plt.tight_layout()
        filename = os.path.join(save_dir, f"{feat}_polar_plot_{distance_mode}.png")
        plt.savefig(filename, dpi=300)
        plt.show()
        plt.close()
