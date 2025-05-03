import argparse
import mdtraj as md
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from tqdm import tqdm
from scipy.signal import savgol_filter

plt.rcParams.update({'font.size': 21, 'font.family': 'Times New Roman'})

def compute_and_plot_rmsd(topol, traj, resname="LIG", heavy_only=True, frame_interval_ns=0.5, verbose=True):
    """
    Compute RMSD of ligand atoms (heavy atoms or all atoms) and plot the result.
    """

    # Step 1: Load trajectory
    t = md.load(traj, top=topol)
    if verbose:
        print("Trajectory loaded.")

    # Step 2: Select ligand atoms
    if heavy_only:
        selection_string = f"resname {resname} and element != H"
    else:
        selection_string = f"resname {resname}"

    small_mol_atoms = t.topology.select(selection_string)
    if len(small_mol_atoms) == 0:
        raise ValueError(f"No atoms found for resname '{resname}' with selection: {selection_string}")

    if verbose:
        atom_type = "heavy atoms" if heavy_only else "all atoms"
        print(f"Selected {len(small_mol_atoms)} {atom_type} from residue name '{resname}'.")

    # Step 3: Compute RMSD (to first frame)
    rmsd = md.rmsd(t, t, 0, atom_indices=small_mol_atoms)

    # Step 4: Smooth RMSD
    window = min(11, len(rmsd) // 2 * 2 + 1)  # odd window length
    smooth_rmsd = savgol_filter(rmsd, window_length=window, polyorder=3)

    # Step 5: Plot RMSD
    output_folder = os.path.dirname(os.path.abspath(traj))
    time_ns = np.arange(len(rmsd)) * frame_interval_ns

    plt.figure(figsize=(8, 6))
    plt.plot(time_ns, rmsd, label='Raw RMSD', color='#8fbcd4', linewidth=1.5)
    plt.plot(time_ns, smooth_rmsd, label='Smoothed RMSD', color='#1f77b4', linewidth=2.5, marker='o', markersize=4)
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD (nm)')
    plt.title('Ligand RMSD')
    plt.legend()
    plt.tight_layout()

    rmsd_plot_path = os.path.join(output_folder, 'ligand_rmsd.png')
    plt.savefig(rmsd_plot_path, dpi=300)
    plt.close()
    if verbose:
        print(f"RMSD plot saved to: {rmsd_plot_path}")

    # Step 6: Save RMSD data
    df_rmsd = pd.DataFrame({'Time (ns)': time_ns, 'RMSD': rmsd, 'Smoothed RMSD': smooth_rmsd})
    rmsd_data_path = os.path.join(output_folder, 'ligand_rmsd_data.xlsx')
    df_rmsd.to_excel(rmsd_data_path, index=False)
    if verbose:
        print(f"RMSD data saved to: {rmsd_data_path}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Ligand RMSD analysis')
    parser.add_argument('--topol', required=True, help='Topology file (e.g., .pdb, .prmtop)')
    parser.add_argument('--traj', required=True, help='Trajectory file (e.g., .dcd, .xtc)')
    parser.add_argument('--resname', default="LIG", help='Ligand residue name (default: LIG)')
    parser.add_argument('--heavy_only', action='store_true', help='Use only heavy atoms (default: True). If not set, use all atoms.')
    args = parser.parse_args()

    compute_and_plot_rmsd(
        args.topol,
        args.traj,
        resname=args.resname,
        heavy_only=args.heavy_only
    )
