# import os
# import mdtraj as md
# import numpy as np
# import pandas as pd
# from tqdm import tqdm
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D

# plt.rcParams.update({'font.size': 18, 'font.family': 'Times New Roman'})

# # ----------------------- Define IPE Categories -----------------------
# # 'A'    → Any atom (no selection criteria, all atoms included)
# # 'NP'   → Non-polar atoms (typically carbons and sulfurs)
# # 'P+'   → Positively charged groups (e.g., NZ, NH groups, or residues ARG, LYS)
# # 'P-'   → Negatively charged groups (e.g., OE, OD groups, or residues ASP, GLU)
# # 'HBD'  → Hydrogen bond donors (atoms that can donate a hydrogen: nitrogen or hydroxyl groups)
# # 'HBA'  → Hydrogen bond acceptors (atoms that can accept a hydrogen: oxygen or OE/OD groups)


# IPE_CATEGORIES = {
#     'A': lambda atom: True,  # Any atom
#     'NP': lambda atom: atom.element.symbol in {'C', 'S'},  # Non-polar atoms
#     'P+': lambda atom: atom.name.startswith('NZ') or atom.name.startswith('NH') or atom.residue.name in {'ARG', 'LYS'},  # Positively charged
#     'P-': lambda atom: atom.name.startswith('OE') or atom.name.startswith('OD') or atom.residue.name in {'ASP', 'GLU'},  # Negatively charged
#     'HBD': lambda atom: atom.element.symbol == 'N' or atom.name.startswith('OH'),  # Hydrogen bond donor
#     'HBA': lambda atom: atom.element.symbol == 'O' or atom.name.startswith('OD') or atom.name.startswith('OE')  # Hydrogen bond acceptor
# }

# # ----------------------- Compute 4D Fingerprint -----------------------

# def compute_4D_fingerprint(topol, traj_file, resname="LIG", grid_size=1.0, cutoff=4.0, n_frames=50):
#     """
#     Compute the 4D-Fingerprint (GCOD approach).
#     """

#     traj = md.load(traj_file, top=topol)
#     if n_frames > traj.n_frames:
#         n_frames = traj.n_frames

#     out_dir = os.path.dirname(os.path.abspath(traj_file))
#     graph_dir = os.path.join(out_dir, "Graph")
#     fingerprint_dir = os.path.join(graph_dir, "4D_Fingerprint")
#     os.makedirs(fingerprint_dir, exist_ok=True)

#     # Select ligand and protein atoms (excluding hydrogens)
#     ligand_atoms = traj.topology.select(f"resname {resname} and element != H")
#     protein_atoms = traj.topology.select("protein and element != H")
#     atoms = ligand_atoms.tolist() + protein_atoms.tolist()

#     # Compute bounding box for the entire trajectory
#     all_coords = traj.xyz[:, atoms, :].reshape(-1, 3)
#     min_corner = all_coords.min(axis=0) - 1.0
#     max_corner = all_coords.max(axis=0) + 1.0

#     grid_x = np.arange(min_corner[0], max_corner[0], grid_size)
#     grid_y = np.arange(min_corner[1], max_corner[1], grid_size)
#     grid_z = np.arange(min_corner[2], max_corner[2], grid_size)

#     # Initialize occupancy grid
#     fingerprint = {}

#     for i, x in enumerate(grid_x):
#         for j, y in enumerate(grid_y):
#             for k, z in enumerate(grid_z):
#                 fingerprint[(i, j, k)] = {key: 0 for key in IPE_CATEGORIES}

#     # ----------------- Loop through trajectory and count occupancies -----------------

#     print("Building 4D-Fingerprint (GCOD) occupancy maps...")
#     for frame in tqdm(range(n_frames)):
#         t_slice = traj.slice(frame)

#         for atom_idx in atoms:
#             atom = traj.topology.atom(atom_idx)
#             coord = t_slice.xyz[0, atom_idx, :]

#             # Find nearest grid point indices
#             ix = int((coord[0] - min_corner[0]) // grid_size)
#             iy = int((coord[1] - min_corner[1]) // grid_size)
#             iz = int((coord[2] - min_corner[2]) // grid_size)

#             key = (ix, iy, iz)

#             if key not in fingerprint:
#                 continue

#             for category, func in IPE_CATEGORIES.items():
#                 if func(atom):
#                     fingerprint[key][category] += 1

#     # ----------------- Normalize & Convert to DataFrame -----------------

#     data = []
#     for (ix, iy, iz), counts in fingerprint.items():
#         total = sum(counts.values())
#         if total > 0:
#             row = {
#                 'x': grid_x[ix],
#                 'y': grid_y[iy],
#                 'z': grid_z[iz]
#             }
#             for cat in IPE_CATEGORIES:
#                 row[cat] = counts[cat] / n_frames  # Average occupancy per frame
#             data.append(row)

#     df = pd.DataFrame(data)
#     df.to_excel(os.path.join(fingerprint_dir, "4D_fingerprint_GCOD.xlsx"), index=False)
#     np.save(os.path.join(fingerprint_dir, "4D_fingerprint_GCOD.npy"), df.to_records(index=False))

#     print(f"4D-Fingerprint saved to: {fingerprint_dir}")

#     _plot_4D_grid(df, fingerprint_dir)

# # ----------------------- Visualization -----------------------

# def _plot_4D_grid(df, save_dir):
#     """
#     3D grid visualization.
#     Different colors represent different categories of occupancy.
#     """
#     fig = plt.figure(figsize=(10, 8))
#     ax = fig.add_subplot(111, projection='3d')

#     categories = ['A', 'NP', 'P+', 'P-', 'HBD', 'HBA']
#     colors = {
#         'A': 'grey',
#         'NP': 'green',
#         'P+': 'blue',
#         'P-': 'red',
#         'HBD': 'cyan',
#         'HBA': 'magenta'
#     }

#     for cat in categories:
#         mask = df[cat] > 0.01  # Filter sparse grid points
#         ax.scatter(df['x'][mask], df['y'][mask], df['z'][mask],
#                    c=colors[cat], s=df[cat][mask] * 300, label=cat, alpha=0.6)

#     ax.set_xlabel('X (nm)')
#     ax.set_ylabel('Y (nm)')
#     ax.set_zlabel('Z (nm)')
#     plt.title("4D Fingerprint GCOD Occupancy")
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(os.path.join(save_dir, "4D_Fingerprint_GCOD.png"), dpi=300)
#     plt.show()

########################################
# Version 1.2
########################################

# import os
# import mdtraj as md
# import numpy as np
# import pandas as pd
# from tqdm import tqdm
# import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
# from matplotlib.colors import LinearSegmentedColormap
# from scipy.spatial import distance
# import pickle

# plt.rcParams.update({'font.size': 18, 'font.family': 'Times New Roman'})

# # ----------------------- Enhanced IPE Categories with Physics-Based Properties -----------------------
# # More detailed atom typing with physical properties
# # Each atom type now includes information about:
# # - Element type and chemical environment
# # - Partial charge ranges
# # - vdW radius ranges
# # - Polarizability estimates
# # - Typical binding energy contributions

# IPE_CATEGORIES = {
#     # Core interaction types
#     'A': {
#         'description': 'Any atom',
#         'filter': lambda atom: True,
#         'color': 'grey',
#         'vdw_radius': 1.7,  # Average vdW radius in Å
#         'energy_weight': 1.0
#     },
#     'NP': {
#         'description': 'Non-polar atoms (hydrophobic)',
#         'filter': lambda atom: atom.element.symbol in {'C', 'S'} and not any(n in atom.name for n in ['NH', 'OH', 'SH']),
#         'color': 'green',
#         'vdw_radius': 1.9,  # Carbon/sulfur typical radius
#         'energy_weight': 0.75  # Weaker interactions
#     },
#     'P+': {
#         'description': 'Positively charged groups',
#         'filter': lambda atom: (atom.element.symbol == 'N' and 
#                                (atom.name.startswith('NZ') or atom.name.startswith('NH') or 
#                                 atom.residue.name in {'ARG', 'LYS', 'HIS'})),
#         'color': 'blue',
#         'vdw_radius': 1.5,  # Nitrogen typical radius
#         'energy_weight': 2.0  # Strong electrostatic interaction
#     },
#     'P-': {
#         'description': 'Negatively charged groups',
#         'filter': lambda atom: (atom.element.symbol == 'O' and 
#                                (atom.name.startswith('OE') or atom.name.startswith('OD') or 
#                                 atom.residue.name in {'ASP', 'GLU'})),
#         'color': 'red',
#         'vdw_radius': 1.4,  # Oxygen typical radius
#         'energy_weight': 2.0  # Strong electrostatic interaction
#     },
#     'HBD': {
#         'description': 'Hydrogen bond donors',
#         'filter': lambda atom: ((atom.element.symbol == 'N' and not atom.name.startswith('NZ')) or 
#                                 atom.name.startswith('OH') or 
#                                 (atom.residue.name in {'SER', 'THR', 'TYR'} and atom.name.startswith('O'))),
#         'color': 'cyan', 
#         'vdw_radius': 1.6,  # Average of N/O
#         'energy_weight': 1.5  # H-bonds are directional and strong
#     },
#     'HBA': {
#         'description': 'Hydrogen bond acceptors',
#         'filter': lambda atom: (atom.element.symbol == 'O' or 
#                                (atom.element.symbol == 'N' and not any(n in atom.name for n in ['NH', 'NZ']))),
#         'color': 'magenta',
#         'vdw_radius': 1.4,  # Mostly oxygen atoms
#         'energy_weight': 1.5  # H-bonds are directional and strong
#     },
#     # Additional specialized categories
#     'ARO': {
#         'description': 'Aromatic rings',
#         'filter': lambda atom: atom.element.symbol == 'C' and atom.residue.name in {'PHE', 'TYR', 'TRP', 'HIS'},
#         'color': 'orange',
#         'vdw_radius': 1.85,  # Aromatic carbon
#         'energy_weight': 1.2  # Pi-pi interactions and CH-pi
#     },
#     'HAL': {
#         'description': 'Halogens (potential halogen bonds)',
#         'filter': lambda atom: atom.element.symbol in {'F', 'Cl', 'Br', 'I'},
#         'color': 'purple',
#         'vdw_radius': 1.9,  # Average halogen radius
#         'energy_weight': 1.1  # Halogen bonds are directional
#     }
# }

# def get_atom_partial_charge(atom, ligand_charges=None):
#     """
#     Estimate partial charges based on atom type.
#     For ligands, uses provided charge map if available.
#     For proteins, uses heuristics based on amino acid and atom type.
#     """
#     if ligand_charges is not None and atom.index in ligand_charges:
#         return ligand_charges[atom.index]
    
#     # Protein charge estimates
#     if atom.residue.name in {'ARG', 'LYS'} and atom.element.symbol == 'N':
#         return 0.4  # Positively charged groups
#     elif atom.residue.name in {'ASP', 'GLU'} and atom.element.symbol == 'O':
#         return -0.4  # Negatively charged groups
#     elif atom.element.symbol == 'O':
#         return -0.2  # Oxygen typically has partial negative charge
#     elif atom.element.symbol == 'N':
#         return -0.1  # Nitrogen typically has partial negative charge
#     elif atom.element.symbol == 'S':
#         return 0.0  # Sulfur can vary
#     else:
#         return 0.0  # Default for carbon and others
        
# def compute_interaction_energy(atom1_coord, atom2_coord, atom1_type, atom2_type, 
#                                atom1_charge, atom2_charge, cutoff=4.0):
#     """
#     Compute simplified interaction energy between two atoms.
#     Includes:
#     1. Lennard-Jones potential (van der Waals)
#     2. Coulomb electrostatics
#     3. Hydrogen bonding potential (simple distance-based)
    
#     Returns energy in arbitrary units scaled to be in range [-10, 10].
#     """
#     dist = np.linalg.norm(atom1_coord - atom2_coord)
    
#     if dist > cutoff:
#         return 0.0
    
#     # Distance-dependent dielectric constant (approximate solvent effect)
#     dielectric = 1.0 + 78.5 * (1.0 - np.exp(-0.5 * dist))
    
#     # 1. Electrostatic term (Coulomb's law with distance-dependent dielectric)
#     electrostatic = 138.94 * atom1_charge * atom2_charge / (dielectric * dist)
    
#     # 2. van der Waals term (simplified Lennard-Jones)
#     sigma = (IPE_CATEGORIES[atom1_type]['vdw_radius'] + IPE_CATEGORIES[atom2_type]['vdw_radius']) / 2.0
#     vdw = 0.0
#     if dist < 2.5:  # Only calculate for relevant distances
#         vdw = 4.0 * ((sigma/dist)**12 - (sigma/dist)**6)
    
#     # 3. Hydrogen bond term (if applicable)
#     hbond = 0.0
#     if (atom1_type == 'HBD' and atom2_type == 'HBA') or (atom1_type == 'HBA' and atom2_type == 'HBD'):
#         if 2.5 < dist < 3.5:  # Typical H-bond distance range
#             # Gaussian-like function peaking at optimal H-bond distance (~2.9Å)
#             hbond = -5.0 * np.exp(-4.0 * (dist - 2.9)**2)
    
#     # 4. Apply interaction weights
#     weight = IPE_CATEGORIES[atom1_type]['energy_weight'] * IPE_CATEGORIES[atom2_type]['energy_weight']
    
#     # Total energy (scaled to range [-10, 10])
#     energy = weight * (0.2 * electrostatic + 0.1 * vdw + hbond)
    
#     return np.clip(energy, -10.0, 10.0)

# def compute_4D_fingerprint(topol, traj_file, resname="LIG", grid_size=1.0, cutoff=4.0, n_frames=50):
#     """
#     Compute enhanced 4D-Fingerprint with physical energy terms and better grid representation.
    
#     Parameters:
#     -----------
#     topol : str
#         Path to topology file
#     traj_file : str
#         Path to trajectory file
#     resname : str
#         Residue name of the ligand
#     grid_size : float
#         Size of grid cells in Angstroms
#     cutoff : float
#         Maximum distance for considering interactions in Angstroms
#     n_frames : int
#         Number of frames to process from trajectory
    
#     Returns:
#     --------
#     None, saves results to disk
#     """
#     print("Loading trajectory...")
#     traj = md.load(traj_file, top=topol)
#     if n_frames > traj.n_frames:
#         n_frames = traj.n_frames
    
#     # Prepare output directories
#     out_dir = os.path.dirname(os.path.abspath(traj_file))
#     graph_dir = os.path.join(out_dir, "Graph")
#     fingerprint_dir = os.path.join(graph_dir, "4D_Fingerprint")
#     os.makedirs(fingerprint_dir, exist_ok=True)

#     # Select ligand and protein atoms (excluding hydrogens)
#     ligand_atoms = traj.topology.select(f"resname {resname} and element != H")
#     protein_atoms = traj.topology.select("protein and element != H")
    
#     if len(ligand_atoms) == 0:
#         raise ValueError(f"No ligand atoms found with resname {resname}")
    
#     # Estimate partial charges (can be replaced with actual charges if available)
#     ligand_charges = {idx: 0.0 for idx in ligand_atoms}  # Placeholder
    
#     # Create atom category mappings for quick lookup
#     atom_categories = {}
#     for idx in np.concatenate([ligand_atoms, protein_atoms]):
#         atom = traj.topology.atom(idx)
#         atom_categories[idx] = []
#         for cat, props in IPE_CATEGORIES.items():
#             if props['filter'](atom):
#                 atom_categories[idx].append(cat)

#     # Compute bounding box from trajectory's protein-ligand complex
#     complex_atoms = np.concatenate([ligand_atoms, protein_atoms])
#     all_coords = traj.xyz[:, complex_atoms, :].reshape(-1, 3)
#     min_corner = all_coords.min(axis=0) - 2.0  # Extra padding
#     max_corner = all_coords.max(axis=0) + 2.0
    
#     # Create evenly spaced grid
#     grid_x = np.arange(min_corner[0], max_corner[0], grid_size)
#     grid_y = np.arange(min_corner[1], max_corner[1], grid_size)
#     grid_z = np.arange(min_corner[2], max_corner[2], grid_size)
    
#     n_x, n_y, n_z = len(grid_x), len(grid_y), len(grid_z)
    
#     # Initialize data structures
#     # First structure: Grid-based fingerprint (better for CNN inputs)
#     # For each category, we'll have a 3D tensor of grid_shape
#     grid_fingerprint = {
#         cat: np.zeros((n_x, n_y, n_z)) for cat in IPE_CATEGORIES
#     }
    
#     # Add energy grids for physical interpretability
#     grid_fingerprint.update({
#         'total_energy': np.zeros((n_x, n_y, n_z)),
#         'electrostatic': np.zeros((n_x, n_y, n_z)),
#         'vdw': np.zeros((n_x, n_y, n_z)),
#         'hbond': np.zeros((n_x, n_y, n_z)),
#         'persistence': np.zeros((n_x, n_y, n_z)),  # How stable the interaction is over time
#         'solvation': np.zeros((n_x, n_y, n_z))     # Estimate of solvent exposure
#     })
    
#     # Second structure: DataFrame for easier analysis and visualization
#     df_data = []
    
#     print("Building enhanced 4D-Fingerprint with physical energy terms...")
#     for frame in tqdm(range(n_frames)):
#         t_slice = traj.slice(frame)
#         frame_grid = {k: np.zeros((n_x, n_y, n_z)) for k in grid_fingerprint}
        
#         # For each ligand atom, find interacting protein atoms
#         for lig_idx in ligand_atoms:
#             lig_atom = traj.topology.atom(lig_idx)
#             lig_coord = t_slice.xyz[0, lig_idx, :]
#             lig_charge = get_atom_partial_charge(lig_atom, ligand_charges)
            
#             # Get ligand atom categories
#             lig_cats = atom_categories[lig_idx]
            
#             # Find protein atoms within cutoff distance
#             for prot_idx in protein_atoms:
#                 prot_atom = traj.topology.atom(prot_idx)
#                 prot_coord = t_slice.xyz[0, prot_idx, :]
#                 prot_charge = get_atom_partial_charge(prot_atom)
                
#                 # Compute distance
#                 dist = np.linalg.norm(lig_coord - prot_coord)
#                 if dist > cutoff:
#                     continue
                
#                 # Get protein atom categories
#                 prot_cats = atom_categories[prot_idx]
                
#                 # Interaction midpoint - where we'll record the interaction
#                 midpoint = (lig_coord + prot_coord) / 2.0
                
#                 # Find grid indices for midpoint
#                 ix = min(max(0, int((midpoint[0] - min_corner[0]) // grid_size)), n_x-1)
#                 iy = min(max(0, int((midpoint[1] - min_corner[1]) // grid_size)), n_y-1)
#                 iz = min(max(0, int((midpoint[2] - min_corner[2]) // grid_size)), n_z-1)
                
#                 # Calculate physical properties for each interaction type
#                 for lig_cat in lig_cats:
#                     for prot_cat in prot_cats:
#                         # Calculate interaction energy components
#                         energy = compute_interaction_energy(
#                             lig_coord, prot_coord, lig_cat, prot_cat, 
#                             lig_charge, prot_charge, cutoff
#                         )
                        
#                         # Compute distance factor (smooth decay)
#                         dist_factor = 1.0 - (dist / cutoff)**2
                        
#                         # Update grid values for both categories with distance weighting
#                         frame_grid[lig_cat][ix, iy, iz] += dist_factor
#                         frame_grid[prot_cat][ix, iy, iz] += dist_factor
                        
#                         # Update energy grids
#                         frame_grid['total_energy'][ix, iy, iz] += energy * dist_factor
                        
#                         # Simple solvation estimate based on burial depth
#                         # (more negative for deeply buried, approximately 0 for exposed)
#                         depth_estimate = -5.0 * (1.0 - dist/cutoff)
#                         frame_grid['solvation'][ix, iy, iz] += depth_estimate
        
#         # Accumulate frame data into overall fingerprint
#         for k in grid_fingerprint:
#             grid_fingerprint[k] += frame_grid[k]
            
#             # Track persistence - nonzero values that stay consistent
#             if frame > 0:
#                 persistent_mask = (frame_grid[k] > 0) & (grid_fingerprint[k] > 0)
#                 grid_fingerprint['persistence'] += persistent_mask
    
#     # Normalize by number of frames
#     for k in grid_fingerprint:
#         grid_fingerprint[k] /= max(1, n_frames)
    
#     # Convert grid data to dataframe for analysis and visualization
#     for ix in range(n_x):
#         for iy in range(n_y):
#             for iz in range(n_z):
#                 # Only include grid points with significant interactions
#                 total_occupancy = sum(grid_fingerprint[cat][ix, iy, iz] for cat in IPE_CATEGORIES)
#                 if total_occupancy > 0.05:  # Threshold to reduce noise
#                     row = {
#                         'x': grid_x[ix],
#                         'y': grid_y[iy],
#                         'z': grid_z[iz],
#                         'total_energy': grid_fingerprint['total_energy'][ix, iy, iz],
#                         'persistence': grid_fingerprint['persistence'][ix, iy, iz] / n_frames,
#                         'solvation': grid_fingerprint['solvation'][ix, iy, iz]
#                     }
                    
#                     # Add category occupancies
#                     for cat in IPE_CATEGORIES:
#                         row[cat] = grid_fingerprint[cat][ix, iy, iz]
                    
#                     df_data.append(row)
    
#     # Create DataFrame and save results
#     df = pd.DataFrame(df_data)
    
#     # Save all outputs
#     df.to_excel(os.path.join(fingerprint_dir, "4D_fingerprint_physics_v2.xlsx"), index=False)
#     np.save(os.path.join(fingerprint_dir, "4D_fingerprint_grid_v2.npy"), grid_fingerprint)
    
#     # Save compressed tensor for deep learning (better format for ML)
#     with open(os.path.join(fingerprint_dir, "4D_fingerprint_tensor_v2.pkl"), 'wb') as f:
#         pickle.dump({
#             'grid_shape': (n_x, n_y, n_z),
#             'grid_origin': min_corner,
#             'grid_spacing': grid_size,
#             'fingerprint': grid_fingerprint
#         }, f)
    
#     print(f"Enhanced 4D-Fingerprint saved to: {fingerprint_dir}")
    
#     # Create visualizations
#     _plot_4D_grid_enhanced(df, fingerprint_dir)

# def _plot_4D_grid_enhanced(df, save_dir):
#     """
#     Enhanced 3D grid visualization with multiple physical property views.
#     """
#     # Plot 1: Interaction categories with improved coloring
#     fig = plt.figure(figsize=(12, 10))
#     ax = fig.add_subplot(111, projection='3d')
    
#     # Determine point sizing based on persistence and energy
#     point_sizes = 50 + 250 * df['persistence']
    
#     # Plot each category separately
#     for cat, props in IPE_CATEGORIES.items():
#         mask = df[cat] > 0.1  # Filter sparse grid points
#         if mask.sum() > 0:  # Only plot if we have points
#             ax.scatter(df['x'][mask], df['y'][mask], df['z'][mask],
#                       c=props['color'], s=point_sizes[mask], label=props['description'],
#                       alpha=0.7, edgecolors='w', linewidth=0.5)
    
#     ax.set_xlabel('X (Å)')
#     ax.set_ylabel('Y (Å)')
#     ax.set_zlabel('Z (Å)')
#     plt.title("4D Fingerprint - Interaction Types")
#     plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2)
#     plt.tight_layout()
#     plt.savefig(os.path.join(save_dir, "4D_Fingerprint_Categories_v2.png"), dpi=300)
#     plt.show()
#     plt.close()
    
#     # Plot 2: Energy landscape
#     fig = plt.figure(figsize=(12, 10))
#     ax = fig.add_subplot(111, projection='3d')
    
#     # Create custom colormap for energy (blue=favorable, white=neutral, red=unfavorable)
#     energy_cmap = LinearSegmentedColormap.from_list("energy", 
#                                                    [(0, 'blue'), (0.5, 'white'), (1.0, 'red')])
    
#     # Normalize energy values to [-1, 1] range for colormap
#     energy_vals = df['total_energy']
#     vmax = max(abs(energy_vals.min()), abs(energy_vals.max()))
#     norm_energy = energy_vals / (vmax + 1e-10)
    
#     # Plot energy points
#     significant = abs(energy_vals) > 0.1
#     sc = ax.scatter(df['x'][significant], df['y'][significant], df['z'][significant],
#                    c=norm_energy[significant], cmap=energy_cmap, 
#                    s=100 * (abs(norm_energy[significant]) + 0.1),
#                    alpha=0.8, edgecolors='k', linewidth=0.2)
    
#     plt.colorbar(sc, label='Interaction Energy (normalized)', pad=0.1)
#     ax.set_xlabel('X (Å)')
#     ax.set_ylabel('Y (Å)')
#     ax.set_zlabel('Z (Å)')
#     plt.title("4D Fingerprint - Energy Landscape")
#     plt.tight_layout()
#     plt.savefig(os.path.join(save_dir, "4D_Fingerprint_Energy_v2.png"), dpi=300)
#     plt.show()
#     plt.close()
    
#     # Plot 3: Persistence/stability view
#     fig = plt.figure(figsize=(12, 10))
#     ax = fig.add_subplot(111, projection='3d')
    
#     persistence_cmap = plt.cm.viridis
#     persistent = df['persistence'] > 0.1
#     sc = ax.scatter(df['x'][persistent], df['y'][persistent], df['z'][persistent],
#                    c=df['persistence'][persistent], cmap=persistence_cmap, 
#                    s=200 * df['persistence'][persistent],
#                    alpha=0.7)
    
#     plt.colorbar(sc, label='Interaction Persistence', pad=0.1)
#     ax.set_xlabel('X (Å)')
#     ax.set_ylabel('Y (Å)')
#     ax.set_zlabel('Z (Å)')
#     plt.title("4D Fingerprint - Interaction Stability")
#     plt.tight_layout()
#     plt.savefig(os.path.join(save_dir, "4D_Fingerprint_Persistence_v2.png"), dpi=300)
#     plt.show()
#     plt.close()

####################################
# V1.3
####################################
import os
import mdtraj as md
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import LinearSegmentedColormap
from scipy.spatial import distance
import pickle

plt.rcParams.update({'font.size': 18, 'font.family': 'Times New Roman'})

# ----------------------- Enhanced IPE Categories with Physics-Based Properties -----------------------
# More detailed atom typing with physical properties
# Each atom type now includes information about:
# - Element type and chemical environment
# - vdW radius ranges
# - Typical binding energy contributions

IPE_CATEGORIES = {
    # Core interaction types
    'A': {
        'description': 'Any atom',
        'filter': lambda atom: True,
        'color': 'grey',
        'vdw_radius': 1.7,  # Average vdW radius in Å
        'energy_weight': 1.0
    },
    'NP': {
        'description': 'Non-polar atoms (hydrophobic)',
        'filter': lambda atom: atom.element.symbol in {'C', 'S'} and not any(n in atom.name for n in ['NH', 'OH', 'SH']),
        'color': 'green',
        'vdw_radius': 1.9,  # Carbon/sulfur typical radius
        'energy_weight': 0.75  # Weaker interactions
    },
    'P+': {
        'description': 'Positively charged groups',
        'filter': lambda atom: (atom.element.symbol == 'N' and 
                               (atom.name.startswith('NZ') or atom.name.startswith('NH') or 
                                atom.residue.name in {'ARG', 'LYS', 'HIS'})),
        'color': 'blue',
        'vdw_radius': 1.5,  # Nitrogen typical radius
        'energy_weight': 2.0  # Strong electrostatic interaction
    },
    'P-': {
        'description': 'Negatively charged groups',
        'filter': lambda atom: (atom.element.symbol == 'O' and 
                               (atom.name.startswith('OE') or atom.name.startswith('OD') or 
                                atom.residue.name in {'ASP', 'GLU'})),
        'color': 'red',
        'vdw_radius': 1.4,  # Oxygen typical radius
        'energy_weight': 2.0  # Strong electrostatic interaction
    },
    'HBD': {
        'description': 'Hydrogen bond donors',
        'filter': lambda atom: ((atom.element.symbol == 'N' and not atom.name.startswith('NZ')) or 
                                atom.name.startswith('OH') or 
                                (atom.residue.name in {'SER', 'THR', 'TYR'} and atom.name.startswith('O'))),
        'color': 'cyan', 
        'vdw_radius': 1.6,  # Average of N/O
        'energy_weight': 1.5  # H-bonds are directional and strong
    },
    'HBA': {
        'description': 'Hydrogen bond acceptors',
        'filter': lambda atom: (atom.element.symbol == 'O' or 
                               (atom.element.symbol == 'N' and not any(n in atom.name for n in ['NH', 'NZ']))),
        'color': 'magenta',
        'vdw_radius': 1.4,  # Mostly oxygen atoms
        'energy_weight': 1.5  # H-bonds are directional and strong
    },
    # Additional specialized categories
    'ARO': {
        'description': 'Aromatic rings',
        'filter': lambda atom: atom.element.symbol == 'C' and atom.residue.name in {'PHE', 'TYR', 'TRP', 'HIS'},
        'color': 'orange',
        'vdw_radius': 1.85,  # Aromatic carbon
        'energy_weight': 1.2  # Pi-pi interactions and CH-pi
    },
    'HAL': {
        'description': 'Halogens (potential halogen bonds)',
        'filter': lambda atom: atom.element.symbol in {'F', 'Cl', 'Br', 'I'},
        'color': 'purple',
        'vdw_radius': 1.9,  # Average halogen radius
        'energy_weight': 1.1  # Halogen bonds are directional
    }
}

def get_atom_partial_charge(atom, ligand_charges=None):
    """
    Estimate partial charges based on atom type.
    For ligands, uses provided charge map if available.
    For proteins, uses heuristics based on amino acid and atom type.
    """
    if ligand_charges is not None and atom.index in ligand_charges:
        return ligand_charges[atom.index]
    
    # Protein charge estimates
    if atom.residue.name in {'ARG', 'LYS'} and atom.element.symbol == 'N':
        return 0.4  # Positively charged groups
    elif atom.residue.name in {'ASP', 'GLU'} and atom.element.symbol == 'O':
        return -0.4  # Negatively charged groups
    elif atom.element.symbol == 'O':
        return -0.2  # Oxygen typically has partial negative charge
    elif atom.element.symbol == 'N':
        return -0.1  # Nitrogen typically has partial negative charge
    elif atom.element.symbol == 'S':
        return 0.0  # Sulfur can vary
    else:
        return 0.0  # Default for carbon and others
        
def compute_interaction_energy(atom1_coord, atom2_coord, atom1_type, atom2_type, 
                               atom1_charge, atom2_charge, cutoff=4.0):
    """
    Compute simplified interaction energy between two atoms.
    Includes:
    1. Lennard-Jones potential (van der Waals)
    2. Coulomb electrostatics
    3. Hydrogen bonding potential (simple distance-based)
    
    Returns energy in arbitrary units scaled to be in range [-10, 10].
    """
    dist = np.linalg.norm(atom1_coord - atom2_coord)
    
    if dist > cutoff:
        return 0.0
    
    # Distance-dependent dielectric constant (approximate solvent effect)
    dielectric = 1.0 + 78.5 * (1.0 - np.exp(-0.5 * dist))
    
    # 1. Electrostatic term (Coulomb's law with distance-dependent dielectric)
    electrostatic = 138.94 * atom1_charge * atom2_charge / (dielectric * dist)
    
    # 2. van der Waals term (simplified Lennard-Jones)
    sigma = (IPE_CATEGORIES[atom1_type]['vdw_radius'] + IPE_CATEGORIES[atom2_type]['vdw_radius']) / 2.0
    vdw = 0.0
    if dist < 2.5:  # Only calculate for relevant distances
        vdw = 4.0 * ((sigma/dist)**12 - (sigma/dist)**6)
    
    # 3. Hydrogen bond term (if applicable)
    hbond = 0.0
    if (atom1_type == 'HBD' and atom2_type == 'HBA') or (atom1_type == 'HBA' and atom2_type == 'HBD'):
        if 2.5 < dist < 3.5:  # Typical H-bond distance range
            # Gaussian-like function peaking at optimal H-bond distance (~2.9Å)
            hbond = -5.0 * np.exp(-4.0 * (dist - 2.9)**2)
    
    # 4. Apply interaction weights
    weight = IPE_CATEGORIES[atom1_type]['energy_weight'] * IPE_CATEGORIES[atom2_type]['energy_weight']
    
    # Total energy (scaled to range [-10, 10])
    energy = weight * (0.2 * electrostatic + 0.1 * vdw + hbond)
    
    return np.clip(energy, -10.0, 10.0)

def compute_4D_fingerprint(topol, traj_file, resname="LIG", grid_size=1.0, cutoff=4.0, n_frames=50):
    """
    Compute enhanced 4D-Fingerprint with physical energy terms and better grid representation.
    Version 1.2 with improved physical interpretability and machine learning compatibility.
    
    Parameters:
    -----------
    topol : str
        Path to topology file
    traj_file : str
        Path to trajectory file
    resname : str
        Residue name of the ligand
    grid_size : float
        Size of grid cells in Angstroms
    cutoff : float
        Maximum distance for considering interactions in Angstroms
    n_frames : int
        Number of frames to process from trajectory
    """
    print("Loading trajectory...")
    traj = md.load(traj_file, top=topol)
    if n_frames > traj.n_frames:
        n_frames = traj.n_frames
    
    # Prepare output directories
    out_dir = os.path.dirname(os.path.abspath(traj_file))
    graph_dir = os.path.join(out_dir, "Graph")
    fingerprint_dir = os.path.join(graph_dir, "4D_Fingerprint")
    os.makedirs(fingerprint_dir, exist_ok=True)

    # Select ligand and protein atoms (excluding hydrogens)
    ligand_atoms = traj.topology.select(f"resname {resname} and element != H")
    protein_atoms = traj.topology.select("protein and element != H")
    
    if len(ligand_atoms) == 0:
        raise ValueError(f"No ligand atoms found with resname {resname}")
    
    # Estimate partial charges (can be replaced with actual charges if available)
    ligand_charges = {idx: 0.0 for idx in ligand_atoms}  # Placeholder
    
    # Create atom category mappings for quick lookup
    atom_categories = {}
    for idx in np.concatenate([ligand_atoms, protein_atoms]):
        atom = traj.topology.atom(idx)
        atom_categories[idx] = []
        for cat, props in IPE_CATEGORIES.items():
            if props['filter'](atom):
                atom_categories[idx].append(cat)

    # Compute bounding box from trajectory's protein-ligand complex
    complex_atoms = np.concatenate([ligand_atoms, protein_atoms])
    all_coords = traj.xyz[:, complex_atoms, :].reshape(-1, 3)
    min_corner = all_coords.min(axis=0) - 2.0  # Extra padding
    max_corner = all_coords.max(axis=0) + 2.0
    
    # Create evenly spaced grid
    grid_x = np.arange(min_corner[0], max_corner[0], grid_size)
    grid_y = np.arange(min_corner[1], max_corner[1], grid_size)
    grid_z = np.arange(min_corner[2], max_corner[2], grid_size)
    
    n_x, n_y, n_z = len(grid_x), len(grid_y), len(grid_z)
    
    # Initialize data structures
    # First structure: Grid-based fingerprint (better for CNN inputs)
    # For each category, we'll have a 3D tensor of grid_shape
    grid_fingerprint = {
        cat: np.zeros((n_x, n_y, n_z)) for cat in IPE_CATEGORIES
    }
    
    # Add energy grids for physical interpretability
    grid_fingerprint.update({
        'total_energy': np.zeros((n_x, n_y, n_z)),
        'electrostatic': np.zeros((n_x, n_y, n_z)),
        'vdw': np.zeros((n_x, n_y, n_z)),
        'hbond': np.zeros((n_x, n_y, n_z)),
        'persistence': np.zeros((n_x, n_y, n_z)),  # How stable the interaction is over time
        'solvation': np.zeros((n_x, n_y, n_z))     # Estimate of solvent exposure
    })
    
    # Initialize original-style fingerprint format too (for compatibility)
    original_fingerprint = {}
    for i, x in enumerate(grid_x):
        for j, y in enumerate(grid_y):
            for k, z in enumerate(grid_z):
                original_fingerprint[(i, j, k)] = {cat: 0 for cat in ['A', 'NP', 'P+', 'P-', 'HBD', 'HBA']}
    
    # Second structure: DataFrame for easier analysis and visualization
    df_data = []
    
    print("Building enhanced 4D-Fingerprint with physical energy terms...")
    for frame in tqdm(range(n_frames)):
        t_slice = traj.slice(frame)
        frame_grid = {k: np.zeros((n_x, n_y, n_z)) for k in grid_fingerprint}
        
        # For each ligand atom, find interacting protein atoms
        for lig_idx in ligand_atoms:
            lig_atom = traj.topology.atom(lig_idx)
            lig_coord = t_slice.xyz[0, lig_idx, :]
            lig_charge = get_atom_partial_charge(lig_atom, ligand_charges)
            
            # Get ligand atom categories
            lig_cats = atom_categories[lig_idx]
            
            # Find protein atoms within cutoff distance
            for prot_idx in protein_atoms:
                prot_atom = traj.topology.atom(prot_idx)
                prot_coord = t_slice.xyz[0, prot_idx, :]
                prot_charge = get_atom_partial_charge(prot_atom)
                
                # Compute distance
                dist = np.linalg.norm(lig_coord - prot_coord)
                if dist > cutoff:
                    continue
                
                # Get protein atom categories
                prot_cats = atom_categories[prot_idx]
                
                # Interaction midpoint - where we'll record the interaction
                midpoint = (lig_coord + prot_coord) / 2.0
                
                # Find grid indices for midpoint
                ix = min(max(0, int((midpoint[0] - min_corner[0]) // grid_size)), n_x-1)
                iy = min(max(0, int((midpoint[1] - min_corner[1]) // grid_size)), n_y-1)
                iz = min(max(0, int((midpoint[2] - min_corner[2]) // grid_size)), n_z-1)
                
                # Calculate physical properties for each interaction type
                for lig_cat in lig_cats:
                    for prot_cat in prot_cats:
                        # Calculate interaction energy components
                        energy = compute_interaction_energy(
                            lig_coord, prot_coord, lig_cat, prot_cat, 
                            lig_charge, prot_charge, cutoff
                        )
                        
                        # Compute distance factor (smooth decay)
                        dist_factor = 1.0 - (dist / cutoff)**2
                        
                        # Update grid values for both categories with distance weighting
                        frame_grid[lig_cat][ix, iy, iz] += dist_factor
                        frame_grid[prot_cat][ix, iy, iz] += dist_factor
                        
                        # Update original-style fingerprint too (for compatibility)
                        if lig_cat in ['A', 'NP', 'P+', 'P-', 'HBD', 'HBA']:
                            original_fingerprint[(ix, iy, iz)][lig_cat] += 1
                        if prot_cat in ['A', 'NP', 'P+', 'P-', 'HBD', 'HBA']:
                            original_fingerprint[(ix, iy, iz)][prot_cat] += 1
                        
                        # Update energy grids
                        frame_grid['total_energy'][ix, iy, iz] += energy * dist_factor
                        
                        # Simple solvation estimate based on burial depth
                        # (more negative for deeply buried, approximately 0 for exposed)
                        depth_estimate = -5.0 * (1.0 - dist/cutoff)
                        frame_grid['solvation'][ix, iy, iz] += depth_estimate
        
        # Accumulate frame data into overall fingerprint
        for k in grid_fingerprint:
            grid_fingerprint[k] += frame_grid[k]
            
            # Track persistence - nonzero values that stay consistent
            if frame > 0:
                persistent_mask = (frame_grid[k] > 0) & (grid_fingerprint[k] > 0)
                grid_fingerprint['persistence'] += persistent_mask
    
    # Normalize by number of frames
    for k in grid_fingerprint:
        grid_fingerprint[k] /= max(1, n_frames)
    
    # Also normalize original fingerprint
    for key in original_fingerprint:
        for cat in original_fingerprint[key]:
            original_fingerprint[key][cat] /= max(1, n_frames)
    
    # Convert grid data to dataframe for analysis and visualization
    for ix in range(n_x):
        for iy in range(n_y):
            for iz in range(n_z):
                # Only include grid points with significant interactions
                total_occupancy = sum(grid_fingerprint[cat][ix, iy, iz] for cat in IPE_CATEGORIES)
                if total_occupancy > 0.05:  # Threshold to reduce noise
                    row = {
                        'x': grid_x[ix],
                        'y': grid_y[iy],
                        'z': grid_z[iz],
                        'total_energy': grid_fingerprint['total_energy'][ix, iy, iz],
                        'persistence': grid_fingerprint['persistence'][ix, iy, iz] / n_frames,
                        'solvation': grid_fingerprint['solvation'][ix, iy, iz]
                    }
                    
                    # Add category occupancies
                    for cat in IPE_CATEGORIES:
                        row[cat] = grid_fingerprint[cat][ix, iy, iz]
                    
                    df_data.append(row)
    
    # Create original style dataframe too
    original_data = []
    for (ix, iy, iz), counts in original_fingerprint.items():
        total = sum(counts.values())
        if total > 0:
            row = {
                'x': grid_x[ix],
                'y': grid_y[iy],
                'z': grid_z[iz]
            }
            for cat in ['A', 'NP', 'P+', 'P-', 'HBD', 'HBA']:
                row[cat] = counts[cat]
            original_data.append(row)
    
    # Create DataFrames and save results
    df = pd.DataFrame(df_data)
    original_df = pd.DataFrame(original_data)
    
    # Save all outputs
    df.to_excel(os.path.join(fingerprint_dir, "4D_fingerprint_physics_v1_2.xlsx"), index=False)
    original_df.to_excel(os.path.join(fingerprint_dir, "4D_fingerprint_GCOD_v1_2.xlsx"), index=False)
    
    # Save numpy files for both formats
    np.save(os.path.join(fingerprint_dir, "4D_fingerprint_grid_v1_2.npy"), grid_fingerprint)
    
    # Save compressed tensor for deep learning (better format for ML)
    with open(os.path.join(fingerprint_dir, "4D_fingerprint_tensor_v1_2.pkl"), 'wb') as f:
        pickle.dump({
            'grid_shape': (n_x, n_y, n_z),
            'grid_origin': min_corner,
            'grid_spacing': grid_size,
            'fingerprint': grid_fingerprint
        }, f)
    
    print(f"Enhanced 4D-Fingerprint saved to: {fingerprint_dir}")
    
    # Create visualizations
    _plot_4D_grid_enhanced(df, fingerprint_dir)
    _plot_4D_grid_original(original_df, fingerprint_dir)

def _plot_4D_grid_original(df, save_dir):
    """
    3D grid visualization using the original style.
    Different colors represent different categories of occupancy.
    """
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')

    categories = ['A', 'NP', 'P+', 'P-', 'HBD', 'HBA']
    colors = {
        'A': 'grey',
        'NP': 'green',
        'P+': 'blue',
        'P-': 'red',
        'HBD': 'cyan',
        'HBA': 'magenta'
    }

    for cat in categories:
        mask = df[cat] > 0.01  # Filter sparse grid points
        ax.scatter(df['x'][mask], df['y'][mask], df['z'][mask],
                  c=colors[cat], s=df[cat][mask] * 300, label=cat, alpha=0.6)

    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    plt.title("4D Fingerprint GCOD Occupancy")
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, "4D_Fingerprint_GCOD_v1_2.png"), dpi=300)
    plt.show()
    plt.close()

def _plot_4D_grid_enhanced(df, save_dir):
    """
    Enhanced 3D grid visualization with multiple physical property views.
    """
    # Plot 1: Interaction categories with improved coloring
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Determine point sizing based on persistence and energy
    point_sizes = 50 + 250 * df['persistence']
    
    # Plot each category separately
    for cat, props in IPE_CATEGORIES.items():
        mask = df[cat] > 0.1  # Filter sparse grid points
        if mask.sum() > 0:  # Only plot if we have points
            ax.scatter(df['x'][mask], df['y'][mask], df['z'][mask],
                      c=props['color'], s=point_sizes[mask], label=props['description'],
                      alpha=0.7, edgecolors='w', linewidth=0.5)
    
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    plt.title("4D Fingerprint - Interaction Types")
    plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2)
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, "4D_Fingerprint_Categories_v1_2.png"), dpi=300)
    plt.show()
    plt.close()
    
    # Plot 2: Energy landscape
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    # Create custom colormap for energy (blue=favorable, white=neutral, red=unfavorable)
    energy_cmap = LinearSegmentedColormap.from_list("energy", 
                                                   [(0, 'blue'), (0.5, 'white'), (1.0, 'red')])
    
    # Normalize energy values to [-1, 1] range for colormap
    energy_vals = df['total_energy']
    vmax = max(abs(energy_vals.min()), abs(energy_vals.max()))
    norm_energy = energy_vals / (vmax + 1e-10)
    
    # Plot energy points
    significant = abs(energy_vals) > 0.1
    sc = ax.scatter(df['x'][significant], df['y'][significant], df['z'][significant],
                   c=norm_energy[significant], cmap=energy_cmap, 
                   s=100 * (abs(norm_energy[significant]) + 0.1),
                   alpha=0.8, edgecolors='k', linewidth=0.2)
    
    plt.colorbar(sc, label='Interaction Energy (normalized)', pad=0.1)
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    plt.title("4D Fingerprint - Energy Landscape")
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, "4D_Fingerprint_Energy_v1_2.png"), dpi=300)
    plt.show()
    plt.close()
    
    # Plot 3: Persistence/stability view
    fig = plt.figure(figsize=(12, 10))
    ax = fig.add_subplot(111, projection='3d')
    
    persistence_cmap = plt.cm.viridis
    persistent = df['persistence'] > 0.1
    sc = ax.scatter(df['x'][persistent], df['y'][persistent], df['z'][persistent],
                   c=df['persistence'][persistent], cmap=persistence_cmap, 
                   s=200 * df['persistence'][persistent],
                   alpha=0.7)
    
    plt.colorbar(sc, label='Interaction Persistence', pad=0.1)
    ax.set_xlabel('X (Å)')
    ax.set_ylabel('Y (Å)')
    ax.set_zlabel('Z (Å)')
    plt.title("4D Fingerprint - Interaction Stability")
    plt.tight_layout()
    plt.savefig(os.path.join(save_dir, "4D_Fingerprint_Persistence_v1_2.png"), dpi=300)
    plt.show()
    plt.close()