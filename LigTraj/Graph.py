# import os
# from . import graph_generator
# from . import masif_embedding

# def build(topol, traj, sdf, resname="LIG", n_frames=50):
#     """
#     Top-level function to build ligand-protein graphs and compute interface features.
#     Generates per-frame graph data and MaSIF-style feature embeddings under Graph/ directory.
#     """
#     out_dir = os.path.dirname(os.path.abspath(traj))
#     graph_dir = os.path.join(out_dir, "Graph")
#     os.makedirs(graph_dir, exist_ok=True)
#     # Step 1: Generate per-frame graphs and ensemble graph
#     graph_generator.generate_graphs(topol, traj, sdf, resname=resname, n_frames=n_frames)
#     # Step 2: Compute MaSIF-style surface feature embeddings
#     masif_embedding.compute_masif_features(topol, traj, sdf, resname=resname, n_frames=n_frames)
#     print("Graph construction and feature extraction completed.")

# if __name__ == "__main__":
#     import argparse
#     parser = argparse.ArgumentParser(description="Top-level Graph module to build graphs and compute features")
#     parser.add_argument('--topol', required=True, help='Topology file')
#     parser.add_argument('--traj', required=True, help='Trajectory file')
#     parser.add_argument('--sdf', required=True, help='Ligand structure file (.sdf)')
#     parser.add_argument('--resname', default="LIG", help='Ligand residue name (default: LIG)')
#     parser.add_argument('--n_frames', type=int, default=50, help='Number of frames to process (default: 50)')
#     args = parser.parse_args()
#     build(args.topol, args.traj, args.sdf, resname=args.resname, n_frames=args.n_frames)

from . import graph_generator
from . import graph_feature
# from .graph_masif_emb import compute_masif_geodesic_embedding
from . import graph_masif_emb


"""
    Build per-frame contact graphs and an ensemble graph.

    Parameters:
    - topol: Path to topology file.
    - traj: Path to trajectory file.
    - sdf: Path to ligand SDF file.
    - resname: Ligand residue name (default: "LIG").
    - n_frames: Number of frames to process (default: 50).
"""

def build(topol, traj, sdf, resname="LIG", n_frames=50):
    generator = graph_generator.GraphGenerator()
    generator.generate_graphs(topol, traj, sdf, resname=resname, n_frames=n_frames)

def feature(topol, traj, sdf, resname="LIG", n_frames=10):
    return graph_feature.compute_masif_features(topol, traj, sdf, resname=resname, n_frames=n_frames)

def masif_embedding(topol, traj, sdf, resname="LIG", cutoff=0.4, n_frames=10, distance_mode="euclidean"):
    """
    Compute MaSIF geodesic feature embeddings for the ligand-protein interface.
    """
    graph_masif_emb.compute_masif_geodesic_embedding(
        topol, traj, sdf,
        resname=resname,
        cutoff=cutoff,
        n_frames=n_frames,
        distance_mode=distance_mode
    )