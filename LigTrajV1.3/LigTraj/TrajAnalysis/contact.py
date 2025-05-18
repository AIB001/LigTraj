import mdtraj as md
import numpy as np
import pandas as pd
import os
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D
import matplotlib.pyplot as plt
import networkx as nx
from rdkit import RDLogger
from matplotlib import colormaps

RDLogger.DisableLog('rdApp.*')

plt.rcParams.update({'font.size': 21, 'font.family': 'Times New Roman'})

def compute_contact_frequency(topol, traj, resname="LIG", distance_cutoff=0.4, n_frames=50, verbose=True):
    t = md.load(traj, top=topol)

    if n_frames > t.n_frames:
        n_frames = t.n_frames

    if verbose:
        print(f"Trajectory loaded. Total frames: {t.n_frames}, analyzing first {n_frames} frames.")

    ligand_atoms = t.topology.select(f"resname {resname} and element != H")
    protein_atoms = t.topology.select("protein and element != H")

    if len(ligand_atoms) == 0 or len(protein_atoms) == 0:
        raise ValueError("No ligand heavy atoms or protein heavy atoms found.")

    contact_counts = np.zeros(len(ligand_atoms), dtype=int)
    contacts = {}

    for frame in tqdm(range(n_frames), desc="Analyzing contacts"):
        distances = md.compute_distances(
            t.slice(frame),
            atom_pairs=[(i, j) for i in ligand_atoms for j in protein_atoms]
        )[0]

        frame_contacts = set()
        for idx, (lig_atom, prot_atom) in enumerate([(i, j) for i in ligand_atoms for j in protein_atoms]):
            if distances[idx] < distance_cutoff:
                lig_idx = np.where(ligand_atoms == lig_atom)[0][0]
                res = t.topology.atom(prot_atom).residue
                aa_name = f"{res.name[0]}-{res.index}"
                frame_contacts.add((lig_idx, aa_name))

        for lig_idx, aa_name in frame_contacts:
            contacts.setdefault((lig_idx, aa_name), 0)
            contacts[(lig_idx, aa_name)] += 1

        for lig_idx, _ in frame_contacts:
            contact_counts[lig_idx] += 1

    frequency = contact_counts / n_frames
    return ligand_atoms, frequency, contacts, t, n_frames

def draw_molecule_contact(mol, frequencies, filename):
    AllChem.Compute2DCoords(mol)

    max_freq = max(frequencies) if max(frequencies) > 0 else 1
    norm_freq = [f / max_freq for f in frequencies]

    cmap = colormaps['coolwarm']

    drawer = rdMolDraw2D.MolDraw2DCairo(600, 600)
    opts = drawer.drawOptions()
    opts.useBWAtomPalette()
    opts.atomHighlightsAreCircles = True
    opts.circleRadius = 1.8  # 更大的圆圈

    highlight_colors = {}
    for idx, freq in enumerate(norm_freq):
        rgb = cmap(freq)[:3]
        highlight_colors[idx] = (rgb[0], rgb[1], rgb[2], 0.6)

    drawer.ClearDrawing()
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer, mol,
        highlightAtoms=list(range(len(frequencies))),
        highlightAtomColors=highlight_colors
    )

    drawer.FinishDrawing()

    with open(filename, 'wb') as f:
        f.write(drawer.GetDrawingText())

def draw_contact_network(ligand_atoms, contacts, t, resname, filename, n_frames):
    G = nx.Graph()
    lig_min = min(ligand_atoms)
    lig_nodes = []
    prot_nodes = []

    edge_list = []

    lig_atom_names = {}
    for lig_idx in np.unique([i for i, _ in contacts.keys()]):
        atom = t.topology.atom(ligand_atoms[lig_idx])
        new_idx = ligand_atoms[lig_idx] - lig_min
        node_name = f"{atom.element.symbol}-{new_idx}"
        G.add_node(node_name, type='ligand')
        lig_nodes.append(node_name)
        lig_atom_names[lig_idx] = node_name

    for _, aa_name in contacts.keys():
        G.add_node(aa_name, type='protein')
        prot_nodes.append(aa_name)

    for (lig_idx, aa_name), count in contacts.items():
        node_name = lig_atom_names[lig_idx]
        G.add_edge(node_name, aa_name, weight=count)
        edge_list.append(((node_name, aa_name), count))

    top_edges = sorted(edge_list, key=lambda x: x[1], reverse=True)[:10]
    top_edge_set = set([e[0] for e in top_edges])

    node_colors = []
    for node in G.nodes():
        if G.nodes[node]['type'] == 'ligand':
            node_colors.append("#2ca02c")
        else:
            node_colors.append("#1f77b4")

    shells = [lig_nodes, list(set(prot_nodes))]
    pos = nx.shell_layout(G, nlist=shells)

    plt.figure(figsize=(10, 8))
    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=800)

    def interpolate_color(c1, c2, t):
        return tuple(c1[i] + (c2[i] - c1[i]) * t for i in range(3))

    pink_rgb = (0.78, 0.08, 0.52)
    orange_rgb = (1.0, 0.65, 0.0)
    blue_rgb = (0.3, 0.6, 1.0)

    for (u, v, d) in G.edges(data=True):
        freq = d['weight'] / n_frames
        if (u, v) in top_edge_set or (v, u) in top_edge_set:
            color = interpolate_color(pink_rgb, orange_rgb, freq)
            alpha = 1.0
        else:
            color = interpolate_color((0.3, 0.7, 1.0), blue_rgb, freq)
            alpha = 0.1 + 0.4 * freq

        nx.draw_networkx_edges(
            G, pos, edgelist=[(u, v)],
            width=2,
            edge_color=[color],
            alpha=alpha
        )

    nx.draw_networkx_labels(
        G, pos,
        font_size=12,
        font_family="Times New Roman",
        font_color="black"
    )

    edge_labels = {}
    for (u, v, d) in G.edges(data=True):
        freq = d['weight'] / n_frames
        if (u, v) in top_edge_set or (v, u) in top_edge_set:
            edge_labels[(u, v)] = f"{freq:.2f}"

    nx.draw_networkx_edge_labels(
        G, pos,
        edge_labels=edge_labels,
        font_size=12,
        font_family="Times New Roman"
    )

    plt.title("Ligand-Residue Contact Network", fontfamily="Times New Roman")
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

def compute_and_plot_contact(topol, traj, sdf, resname="LIG", distance_cutoff=0.4, n_frames=50):
    ligand_atoms, freq, contacts, t, total_frames = compute_contact_frequency(
        topol, traj, resname=resname, distance_cutoff=distance_cutoff, n_frames=n_frames
    )

    out_dir = os.path.dirname(os.path.abspath(traj))
    df = pd.DataFrame({
        "AtomIndex": ligand_atoms,
        "Frequency": freq
    })
    freq_file = os.path.join(out_dir, "ligand_contact_frequency.xlsx")
    df.to_excel(freq_file, index=False)
    print(f"Contact frequency data saved to: {freq_file}")

    mol = Chem.MolFromMolFile(sdf)
    if mol is None:
        raise ValueError("Failed to read molecule from SDF.")

    mol_fig_path = os.path.join(out_dir, "ligand_contact_frequency.png")
    draw_molecule_contact(mol, freq, mol_fig_path)
    print(f"Molecule contact figure saved to: {mol_fig_path}")

    network_path = os.path.join(out_dir, "ligand_residue_contact_network.png")
    draw_contact_network(ligand_atoms, contacts, t, resname, network_path, total_frames)
    print(f"Contact network figure saved to: {network_path}")
