import os
import subprocess
import shutil
import Bio.PDB
import Bio.PDB.Polypeptide as pp
import biographs as bg
import multiprocessing as mp

def MutationsDict(file, positions=None):
    """Get dictionary with lists of mutations per position in protein, ignore
    positions without residue in pdb file.

    Parameters:
        file (string): pdb file to get mutations from
        positions: list of tuples of the form (chain, first, last) for positions
                   to mutate for all other aminoacids. If None, mutates all
                   positions in all chains

    Returns:
        dict with keys :aa:chain:position, each containing lists with
        :aa:chain:position:mutated_aa for all mutations

    """

    # Sorted list of one letter amino acids
    AA = list(Bio.PDB.Polypeptide.aa1)
    # Generate model of original pdb file
    model = bg.Pmolecule(file).model
    # Dict to store mutations
    mutations = dict()

    for chain_id, first, last in positions:
        # Get chain corresponding to chain_id given
        chain = next(chain for chain in model.get_chains() if chain.id == chain_id)
        for residue in chain:
            if pp.is_aa(residue):
                code = pp.three_to_one(residue.get_resname())
                position = residue.id[1]
                prefix = code+chain_id+str(position)
                # Only save positions between first and last
                if position in range(first, last +1):
                        mutations[prefix] = [prefix+aa for aa in AA if aa!=code]
    return mutations

def GetMutations(path, protein, mutations, foldx_path, out_path=None, use_mp=None):
    """Get mutated pdb files using FoldX.

    Parameters:
        path (string): path where original pdb file is stored, if out_path is None,
                       all resulting files will be stored here
        protein (string): name of original pdb file, without .pdb suffix
        mutations (dict): contains mutations to be made, as returned by MutationsDict
        foldx_path (string): path where FoldX software is stored
        out_path (string): path where resulting files will be stored. Default None
        use_mp (int): number of pools to run if multiprocessing, if None, multiprocessing
                  is not used

    Returns:
        None, only output files

    """
    positions = mutations.keys()
    assert os.path.exists(path), f"{path} does not exist."
    if not out_path:
        out_path = path
    else:
        if not os.path.exists(out_path):
            os.makedirs(out_path)
        # Copy original pdb file to out_path
        original_protein = os.path.join(path, f"{protein}.pdb")
        shutil.copy(original_protein, out_path)

    # Mutations will be made according to :aa:chain:position items in positions,
    # generating an individual list and config file for each, which are necessary
    # to run FoldX BuildModel
    for position in positions:
        IndList(out_path, position, mutations[position])
        ConfigFile(out_path, protein, position)

    if use_mp:
        with mp.Manager() as manager:
            mutations = manager.dict(mutations)
            inputs = [(out_path, protein, position, mutations[position], foldx_path)
                            for position in positions]
            with mp.Pool(use_mp) as pool:
                pool.starmap(Mutate, inputs)
    else:
        for position in positions:
            Mutate(out_path, protein, position, mutations[position], foldx_path)
    return


def Mutate(path, protein, prefix, mutations_list, foldx_path):
    """Call FoldX BuildModel through config_file, rename resulting .pdb files"""
    config_file = os.path.join(path, f"config_{protein}_{prefix}.cfg")
    mutate = f"{foldx_path} -f {config_file}"
    subprocess.check_call(mutate, shell = True)
    for i, x in enumerate(mutations_list):
        source = os.path.join(path, f"{protein}_{i+1}.pdb")
        dest = os.path.join(path, f"{protein}_{x}.pdb")
        os.rename(source, dest)
        source_WT = os.path.join(path, f"WT_{protein}_{i+1}.pdb")
        dest_WT = os.path.join(path, f"WT_{protein}_{x}.pdb")
        os.rename(source_WT, dest_WT)
    return

def IndList(path, prefix, mutation_list):
    """Write individual_list txt file."""
    file = os.path.join(path, f"individual_list_{prefix}.txt")
    indlist = open(file, "w+")
    for mutation in mutation_list:
        indlist.write(f"{mutation};\n")
    indlist.close()
    return

def ConfigFile(path, protein, prefix):
    """Write config file."""
    file = os.path.join(path, f"config_{protein}_{prefix}.cfg")
    config = open(file, "w+")
    config.write("command=BuildModel\n")
    config.write(f"pdb={protein}.pdb\n")
    config.write(f"pdb-dir = {path}\n")
    ind_list = os.path.join(path, f"individual_list_{prefix}.txt")
    config.write(f"mutant-file={ind_list}\n")
    config.write(f"output-dir={path}\n")
    config.write(f"output-file={protein}_{prefix}\n")
    config.close()
    return
