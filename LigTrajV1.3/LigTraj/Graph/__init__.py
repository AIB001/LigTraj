from . import graph_generator
from . import graph_feature
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