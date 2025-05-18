from .rmsd import compute_and_plot_rmsd
from .contact import compute_and_plot_contact
from .covariance import compute_and_plot_covariance
# from .tsne import compute_and_plot_tsne  

def rmsd(topol, traj, resname="LIG", heavy_only=True, frame_interval_ns=0.5):
    return compute_and_plot_rmsd(topol, traj, resname, heavy_only, frame_interval_ns)

def contact(topol, traj, sdf, resname="LIG", distance_cutoff=0.4, n_frames=10):
    return compute_and_plot_contact(topol, traj, sdf, resname, distance_cutoff, n_frames)

def covariance(topol, traj, resname="LIG"):
    return compute_and_plot_covariance(topol, traj, resname)

# def tsne(topol, traj, resname="LIG", feature_type="coordinates", se3_invariant=True,
#          perplexity=30, n_iter=2000):
#     return compute_and_plot_tsne(topol, traj, resname, feature_type, se3_invariant,
#                                  perplexity, n_iter)
