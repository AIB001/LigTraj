from setuptools import setup, find_packages

setup(
    name="LigTraj",
    version="1.1.0",
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
        "sklearn"
    ],
    author="A.I.B.",
    description="Ligand trajectory analysis toolkit",
    license="Institute of Quantitative Biology, Zhejiang University"
)
