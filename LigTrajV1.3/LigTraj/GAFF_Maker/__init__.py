from .GAFF_Maker import GAFFMaker

def GAFF_Maker(protein, ligand, output_path, overwrite=False):
    """
    Build up protein-ligand MD system with GAFF and AMBER force field
    """
    gaffmaker = GAFFMaker(protein, ligand, output_path, overwrite)
    gaffmaker.run()