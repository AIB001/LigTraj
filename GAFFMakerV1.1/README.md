# GAFFMakerV1.1 -- By A.I.B.
This tool is used to streamline the process of modeling protein-ligand molecular dynamics(MD) simulation with GROMACS. This scripts can use for high throughput modeling.

## Usage
````bash
bash GAFFMaker.sh <path_to_protein&ligands> <output_path>
````
Note: All ligands and protein should be out in the same folder; the name should contain `_<pdbid>_ligand_` and `_<pdbid>_protein_`

## Requirement 
GROMACS version:     2024.2

AmberTools Version Higher than V20

````bash
conda install -c conda-forge ambertools=23
````


## Data Access
We strongly recommend you to obtain the data from [PDB Binder](https://www.pdbbind-plus.org.cn/)
