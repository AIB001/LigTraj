import numpy as np
import mdtraj as md
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from scipy.stats import gaussian_kde
import matplotlib.pyplot as plt

plt.rcParams.update({
    'font.family' : 'Times New Roman',
    'font.size': 21
})

def compute_and_plot_tsne(topol, traj, resname="LIG", feature_type='distances',
                          se3_invariant=True, perplexity=30, n_iter=2000):
    """
    Compute and visualize t-SNE embedding of ligand conformations.

    Parameters:
    - topol: Topology file path.
    - traj: Trajectory file path.
    - resname: Ligand residue name.
    - feature_type: 'coordinates' or 'distances'.
    - se3_invariant: Whether to apply SE(3)-invariant features.
    - perplexity, n_iter: t-SNE parameters.
    """

    print("Loading trajectory...")
    t = md.load(traj, top=topol)

    sm_atoms = t.topology.select(f"resname {resname} and element != H")
    if len(sm_atoms) == 0:
        raise ValueError(f"No heavy atoms found for resname {resname}.")

    print(f"Selected {len(sm_atoms)} heavy atoms from residue name '{resname}'.")

    if feature_type == 'coordinates':
        coords = t.xyz[:, sm_atoms, :]

        if se3_invariant:
            print("Applying translation and rotation alignment for SE(3) invariance...")
            ref_coords = coords[0]
            ref_centroid = ref_coords.mean(axis=0)
            coords -= ref_centroid

            ref_frame_centered = ref_coords - ref_centroid
            ref_cov = np.cov(ref_frame_centered.T)
            ref_eigvals, ref_evecs = np.linalg.eigh(ref_cov)
            idx = np.argsort(ref_eigvals)[::-1]
            ref_evecs = ref_evecs[:, idx]

            aligned_coords = []
            for frame in coords:
                cov = np.cov(frame.T)
                eigvals, evecs = np.linalg.eigh(cov)
                idx = np.argsort(eigvals)[::-1]
                evecs = evecs[:, idx]

                for i in range(evecs.shape[1]):
                    if np.dot(ref_evecs[:, i], evecs[:, i]) < 0:
                        evecs[:, i] *= -1

                rotation = evecs.T
                aligned_coords.append(frame @ rotation)

            coords = np.array(aligned_coords)

        features = coords.reshape(len(t), -1)

    elif feature_type == 'distances':
        pairs = []
        n_atoms = len(sm_atoms)
        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                pairs.append([sm_atoms[i], sm_atoms[j]])
        features = md.compute_distances(t, pairs)

    else:
        raise ValueError("feature_type must be 'coordinates' or 'distances'")

    scaler = StandardScaler()
    features = scaler.fit_transform(features)

    print("Running t-SNE embedding...")
    tsne = TSNE(n_components=2, perplexity=perplexity, n_iter=n_iter, random_state=42)
    embedding = tsne.fit_transform(features)

    scaler_tsne = MinMaxScaler()
    embedding = scaler_tsne.fit_transform(embedding)

    kde = gaussian_kde(embedding.T)
    density = kde(embedding.T)
    entropy = -np.mean(np.log(density + 1e-10))
    print(f"Conformational Entropy: {entropy:.2f}")

    plt.figure(figsize=(8, 6))
    sc = plt.scatter(embedding[:, 0], embedding[:, 1],
                     c=np.arange(len(t)), cmap='viridis', s=8, alpha=0.8)
    cbar = plt.colorbar(sc)
    cbar.set_label('Frame Number')
    plt.xlabel('t-SNE 1')
    plt.ylabel('t-SNE 2')
    plt.title('Conformational Landscape')
    plt.tight_layout()
    plt.show()
