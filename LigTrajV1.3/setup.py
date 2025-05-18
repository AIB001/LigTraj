from setuptools import setup, find_packages

setup(
    name="LigTraj",
    version="1.3",
    packages=find_packages(),
    install_requires=[
        "mdtraj",
        "numpy",
        "pandas",
        "matplotlib",
        "networkx",
        "rdkit",
        "tqdm",
        "scipy",
        "scikit-learn",
        "torch",
        "torch_geometric",
        "openpyx"
    ],
    author="A.I.B.",
    description="Ligand trajectory analysis toolkit",
    license="Institute of Quantitative Biology, Zhejiang University"
)