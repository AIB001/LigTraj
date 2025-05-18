import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from tqdm import tqdm
import matplotlib.colors as mcolors

plt.rcParams.update({'font.size': 21, 'font.family': 'Times New Roman'})

def compute_covariance_matrix(positions):
    mean_pos = np.mean(positions, axis=0)
    centered = positions - mean_pos
    cov_matrix = np.cov(centered.T)
    return cov_matrix

def plot_matrix(matrix, title, filename, atom_indices, topology, show=True):
    labels = []
    for atom_idx in atom_indices:
        atom = topology.atom(atom_idx)
        labels.append(f"{atom.element.symbol}_{atom.index}")

    plt.figure(figsize=(10, 8))
    cmap = mcolors.LinearSegmentedColormap.from_list("custom_cmap",
        ["#6A5ACD", "#90EE90", "#ff7979"])
    im = plt.imshow(matrix, cmap=cmap)
    plt.colorbar(im)
    plt.title(title)
    num_atoms = len(atom_indices)
    num = num_atoms * 3
    y_positions = [i * 3 + 1 for i in range(num_atoms)]
    plt.yticks(y_positions, labels, fontsize=12)
    x_positions = [i * 3 + 1 for i in range(num_atoms)]
    plt.xticks(x_positions, labels, rotation=90, fontsize=12)
    for pos in range(3, num, 3):
        plt.axhline(pos - 0.5, color='white', linewidth=0.5)
        plt.axvline(pos - 0.5, color='white', linewidth=0.5)
    plt.tight_layout()
    if show:
        plt.show()
    else:
        plt.savefig(filename, dpi=300)
        plt.close()

def compute_and_plot_covariance(topol, traj, resname="LIG", output_folder=None, save_plots=True, verbose=True):
    steps = [
        "Loading trajectory",
        "Selecting ligand heavy atoms",
        "Extracting atomic coordinates",
        "Computing covariance matrix",
        "Plotting covariance matrix heatmap",
        "Eigen decomposition",
        "Plotting covariance matrices reconstructed from top 5 principal components"
    ]
    with tqdm(total=len(steps), desc="Covariance Analysis", ncols=80) as pbar:
        t = md.load(traj, top=topol)
        pbar.set_description("Step 1/7: Trajectory loaded")
        pbar.update(1)
        selection_string = f"resname {resname} and element != H"
        small_mol_atoms = t.topology.select(selection_string)
        if len(small_mol_atoms) == 0:
            raise ValueError(f"No heavy atoms found for resname '{resname}'.")
        if verbose:
            print(f"Selected {len(small_mol_atoms)} heavy atoms from residue name '{resname}'.")
        pbar.set_description("Step 2/7: Ligand heavy atoms selected")
        pbar.update(1)
        positions = t.xyz[:, small_mol_atoms, :].reshape(t.n_frames, -1)
        pbar.set_description("Step 3/7: Atomic coordinates extracted")
        pbar.update(1)
        cov = compute_covariance_matrix(positions)
        pbar.set_description("Step 4/7: Covariance matrix computed")
        pbar.update(1)
        if output_folder is None:
            output_folder = os.path.dirname(os.path.abspath(traj))
        cov_matrix_path = os.path.join(output_folder, 'covariance_matrix.png')
        plot_matrix(cov, 'Covariance Matrix (Cartesian Coordinates)',
                    cov_matrix_path, atom_indices=small_mol_atoms, topology=t.topology,
                    show=not save_plots)
        pbar.set_description("Step 5/7: Covariance matrix heatmap plotted")
        pbar.update(1)
        eigvals, eigvecs = np.linalg.eigh(cov)
        idx = np.argsort(eigvals)[::-1]
        eigvals = eigvals[idx]
        eigvecs = eigvecs[:, idx]
        eigvecs_top5 = eigvecs[:, :5]
        df = pd.DataFrame(eigvecs_top5, columns=[f'PC{i+1}' for i in range(5)])
        excel_path = os.path.join(output_folder, 'top5_eigenvectors.xlsx')
        df.to_excel(excel_path, index=False)
        if verbose:
            print(f"Top 5 eigenvectors saved to: {excel_path}")
        pbar.set_description("Step 6/7: Eigen decomposition completed")
        pbar.update(1)
        for i in range(5):
            eigval = eigvals[i]
            eigvec = eigvecs[:, i].reshape(-1, 1)
            reconstructed_cov = eigval * (eigvec @ eigvec.T)
            filename = os.path.join(output_folder, f'reconstructed_cov_PC{i+1}.png')
            plot_matrix(reconstructed_cov,
                        f'Reconstructed Cov Matrix PC{i+1}',
                        filename,
                        atom_indices=small_mol_atoms,
                        topology=t.topology,
                        show=not save_plots)
        pbar.set_description("Step 7/7: Reconstructed covariance matrices plotted")
        pbar.update(1)
    if verbose:
        print("Analysis completed ✅")
        print("Top 10 eigenvalues (descending order):")
        for i in range(min(10, len(eigvals))):
            print(f"Eigenvalue {i+1}: {eigvals[i]:.4f}")
