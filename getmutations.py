import os
import Bio.PDB
import Bio.PDB.Polypeptide as pp
import biographs as bp

def MutationsDict(path, prot, pdbs=True):
    """Get dictionary with lists of mutations per position in protein.

    Parameters:
        path (string): path where prot pdb file is located, and mutated
                       pdbs if pdbs is True
        prot (string): name of original pdb file
        pdbs (bool): default True. If mutated pdb files exist, only keep
                    mutations with corresponding pdb file. If false,
                    keeps all mutations

    Returns:
        dict with keys :aa:chain:position, each containing lists with
        :aa:chain:position:mutated_aa for all mutations

    """

    # Sorted list of one letter amino acids
    AA = list(Bio.PDB.Polypeptide.aa1)
    # Generate model of original pdb file
    model = bg.Pmolecule(os.path.join(path, f"{prot}.pdb")).model
    # Dict to store mutations
    mutations = dict()

    for chain in model.get_chains():
        for residue in chain:
            if pp.is_aa(residue):
                code = pp.three_to_one(residue.get_resname())
                chain = residue.parent.id
                position = str(residue.id[1])
                prefix = code+chain+position
                if pdbs is True:
                    first_mutation = os.path.join(path, f"{prot}_{prefix+AA[0]}.pdb")
                    # Assume if first mutation exists, all mutations do
                    # Excludes mutating residue for itself
                    if os.path.exists(first_mutation):
                        mutations[prefix] = [prefix+aa for aa in AA if aa!=prefix[0]]
                else:
                    mutations[prefix] = [prefix+aa for aa in AA if aa!=prefix[0]]

    return mutations
