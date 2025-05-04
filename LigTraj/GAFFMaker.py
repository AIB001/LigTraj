from . import gaff_maker

def gaff_modeler(protein_ligand_path, output_path):
    return gaff_maker.model(protein_ligand_path, output_path)
