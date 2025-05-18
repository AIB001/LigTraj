from .FourDFingerprint import compute_4D_fingerprint as _compute_4D_fingerprint
from .tsne import compute_and_plot_tsne

def compute_4D_fingerprint(topol, traj, resname="LIG", grid_size=0.5, cutoff=0.4, n_frames=50):
    """
    Computing 4D Fingerprint
    """
    return _compute_4D_fingerprint(topol, traj, resname=resname, grid_size=grid_size, cutoff=cutoff, n_frames=n_frames)

def tsne(topol, traj, resname="LIG", feature_type="coordinates", se3_invariant=True):
    """
    Performing t-SNE Calculation with SE3 invarient 
    
    """
    return compute_and_plot_tsne(topol, traj, resname="LIG", feature_type="coordinates", se3_invariant=True)