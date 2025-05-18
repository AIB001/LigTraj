import os
import mdtraj as md
import numpy as np
import networkx as nx
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
import torch
from torch_geometric.data import Data
import pickle
import matplotlib.pyplot as plt
from matplotlib import colormaps, colors

plt.rcParams.update({'font.size': 21, 'font.family': 'Times New Roman'})

class GraphGenerator:
    def __init__(self, cutoff=0.4):
        self.cutoff = cutoff  # 0.4 nm == 4 Å

    def generate_graphs(self, topol, traj_file, sdf, resname="LIG", n_frames=50):
        traj = md.load(traj_file, top=topol)
        if n_frames > traj.n_frames:
            n_frames = traj.n_frames

        ligand_atoms = traj.topology.select(f"resname {resname} and element != H")
        protein_atoms = traj.topology.select("protein and element != H")

        if len(ligand_atoms) == 0 or len(protein_atoms) == 0:
            raise ValueError("No ligand or protein heavy atoms found.")

        mol = Chem.MolFromMolFile(sdf)
        if mol is None:
            raise ValueError("Failed to read ligand SDF.")
        AllChem.Compute2DCoords(mol)

        out_dir = os.path.dirname(os.path.abspath(traj_file))
        graph_dir = os.path.join(out_dir, "Graph")
        os.makedirs(graph_dir, exist_ok=True)

        graphs = []
        ensemble_edges = {}
        distance_sums = {}
        all_protein_nodes = set()

        print(f"Building graphs per frame:")

        lig_min = min(ligand_atoms)
        lig_node_names = {}
        for atom_idx in ligand_atoms:
            atom = traj.topology.atom(atom_idx)
            lig_node_names[atom_idx] = f"{atom.element.symbol}-{atom_idx - lig_min}"

        for frame in tqdm(range(n_frames)):
            t_slice = traj.slice(frame)
            contacts, distances = self._get_contacts(t_slice, ligand_atoms, protein_atoms)

            G = self._build_network(t_slice, ligand_atoms, protein_atoms, contacts, distances, lig_node_names)

            with open(os.path.join(graph_dir, f"graph_frame{frame}.gpickle"), 'wb') as f:
                pickle.dump(G, f)

            graphs.append(G)

            for (lig_idx, prot_idx), dist in distances.items():
                lig_node = lig_node_names[lig_idx]
                res = traj.topology.atom(prot_idx).residue
                prot_node = f"{res.name[0]}-{res.index}"
                edge = (lig_node, prot_node)
                ensemble_edges[edge] = ensemble_edges.get(edge, 0) + 1
                distance_sums[edge] = distance_sums.get(edge, 0.0) + dist
                all_protein_nodes.add(prot_node)

        G_ensemble = nx.Graph()

        for node in lig_node_names.values():
            G_ensemble.add_node(node, type='ligand')

        for prot_node in all_protein_nodes:
            G_ensemble.add_node(prot_node, type='protein')

        for edge, count in ensemble_edges.items():
            avg_distance = distance_sums[edge] / count
            G_ensemble.add_edge(edge[0], edge[1], weight=count / n_frames, distance=avg_distance)

        with open(os.path.join(graph_dir, "graph_ensemble.gpickle"), 'wb') as f:
            pickle.dump(G_ensemble, f)

        print(f"Graphs saved to {graph_dir}")
        self._draw_network(G_ensemble, os.path.join(graph_dir, "graph_ensemble.png"), traj, ligand_atoms, lig_node_names)

        pyg_data_list = [self._to_pyg_data(G) for G in graphs]
        torch.save(pyg_data_list, os.path.join(graph_dir, "pyg_graphs.pt"))
        print("Graph dataset ready for deep learning.")

    def _get_contacts(self, traj, ligand_atoms, protein_atoms):
        pairs = [(i, j) for i in ligand_atoms for j in protein_atoms]
        distances = md.compute_distances(traj, pairs)[0]
        contacts = set()
        dist_dict = {}
        for idx, d in enumerate(distances):
            lig_idx, prot_idx = pairs[idx]
            if d < self.cutoff:
                contacts.add((lig_idx, prot_idx))
                dist_dict[(lig_idx, prot_idx)] = d
        return contacts, dist_dict

    def _build_network(self, traj, ligand_atoms, protein_atoms, contacts, distances, lig_node_names):
        G = nx.Graph()

        for atom_idx, node_name in lig_node_names.items():
            G.add_node(node_name, type='ligand')

        if contacts:
            for lig_idx, prot_idx in contacts:
                lig_node = lig_node_names[lig_idx]
                res = traj.topology.atom(prot_idx).residue
                prot_node = f"{res.name[0]}-{res.index}"
                G.add_node(prot_node, type='protein')
                G.add_edge(lig_node, prot_node, distance=distances[(lig_idx, prot_idx)])
        return G

    def _draw_network(self, G, filename, traj, ligand_atoms, lig_node_names):
        prot_nodes = [n for n, d in G.nodes(data=True) if d.get('type') == 'protein']

        # Ligand真实排布
        lig_positions = {}
        lig_atom_positions = traj.xyz[0][ligand_atoms, :]
        lig_atom_indices = list(lig_node_names.keys())
        for idx, pos in enumerate(lig_atom_positions):
            node_name = lig_node_names[lig_atom_indices[idx]]
            lig_positions[node_name] = pos[:2]

        # Normalize ligand coords
        lig_coords = np.array(list(lig_positions.values()))
        center = lig_coords.mean(axis=0)
        lig_coords -= center
        scale = 1.0 / np.max(np.abs(lig_coords))
        lig_coords *= scale

        for i, node in enumerate(lig_positions):
            lig_positions[node] = lig_coords[i]

        # 蛋白质排布：圆形布局
        angle_step = 2 * np.pi / len(prot_nodes) if prot_nodes else 0
        prot_positions = {}
        for i, node in enumerate(prot_nodes):
            angle = i * angle_step
            prot_positions[node] = (1.5 * np.cos(angle), 1.5 * np.sin(angle))

        pos = {}
        pos.update(lig_positions)
        pos.update(prot_positions)

        plt.figure(figsize=(10, 8))
        nx.draw_networkx_nodes(G, pos,
                               nodelist=list(lig_positions.keys()), node_color="#2ca02c", node_size=800)
        nx.draw_networkx_nodes(G, pos,
                               nodelist=prot_nodes, node_color="#1f77b4", node_size=800)

        # 边颜色根据平均距离渐变
        distances = [d['distance'] for _, _, d in G.edges(data=True)]
        norm = colors.Normalize(vmin=min(distances), vmax=max(distances))
        cmap = colormaps['plasma']
        edge_colors = [cmap(norm(d)) for d in distances]

        nx.draw_networkx_edges(G, pos, edge_color=edge_colors, width=2)

        nx.draw_networkx_labels(G, pos, font_size=12, font_family="Times New Roman", font_color="black")

        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])

        ax = plt.gca()
        cbar = plt.colorbar(sm, ax=ax)
        cbar.set_label("Average Distance (nm)")

        plt.title("Ensemble Contact & Distance Network", fontfamily="Times New Roman")
        plt.tight_layout()
        plt.savefig(filename, dpi=300)
        plt.close()

    def _to_pyg_data(self, G):
        node_mapping = {node: i for i, node in enumerate(G.nodes)}
        edge_index = []
        edge_attr = []
        for u, v, d in G.edges(data=True):
            edge_index.append([node_mapping[u], node_mapping[v]])
            edge_attr.append([d.get('distance', 0.0)])
        if len(edge_index) > 0:
            edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()
            edge_attr = torch.tensor(edge_attr, dtype=torch.float)
        else:
            edge_index = torch.empty((2, 0), dtype=torch.long)
            edge_attr = torch.empty((0, 1), dtype=torch.float)
        x = torch.ones((len(G.nodes), 1))
        return Data(x=x, edge_index=edge_index, edge_attr=edge_attr)
