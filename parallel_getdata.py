import os
import multiprocessing as mp
import iterate as it
import numpy as np
import networkx as nx
import Bio.PDB
from ctypes import c_char_p
from multiprocessing import sharedctypes



def GetData(tuple):
    """Get nodes, edges, weight and distance from mutation.

    Parameters:
        tuple containing:
            protein (str): name of original protein pdb file
            mutation (str): :aa:chain:position:mutated_aa
            path (str): path were pdbs are found, assumes mutated file is named
                        'protein_mutation.pdb'
            threshold (float): threshold to generate network

    Uses:
        original_matrix = matrix from original protein
        nodes, edges, weighst, distance = np arrays to write data in
        positions = list of positions, to get corresponding index to write in


    Returns: None
    """
    protein = tuple[0]
    mutation = tuple[1]
    path = tuple[2]
    threshold = tuple[3]

    # Generate network for current mutation
    current_path = os.path.join(path, f"{protein}_{mutation}.pdb")
    current_prot = it.Pmolecule(current_path)
    current = current_prot.network(cutoff=threshold)

    # Obtain the absolute difference in terms of adjacency
    # matrices: the perturbation network.
    current_matrix = nx.adjacency_matrix(current).toarray()
    original_matrix_np = np.ctypeslib.as_array(original_matrix)
    difference = np.abs(original_matrix_np - current_matrix)
    perturbation_network = nx.from_numpy_array(difference)

    # Remove isolates for accurate perturbation network node count
    perturbation_network.remove_nodes_from(
        list(nx.isolates(perturbation_network)))

    # Corresponding row in array according to mutation
    assert mutation[-1] in AA, \
        f"{mutation[-1]} not one of {Bio.PDB.Polypeptide.aa1}"
    aa_index = AA.index(mutation[-1])

    # Corresponding column in array according to position
    assert mutation[:-1] in positions, \
        f"{mutation[:-1]} not one of amino acids in the protein"
    index = positions.index(mutation[:-1])

    # Information obtained from perturbation network
    nodes[aa_index][index] = it.GetNodes(perturbation_network)
    edges[aa_index][index] = it.GetEdges(perturbation_network)
    weight[aa_index][index] = it.GetWeight(perturbation_network)
    distance[aa_index][index] = it.GetDistance(perturbation_network)

    return

if __name__ == '__main__':
    with mp.Manager() as manager:

        # original_data matrix missing here, add when iterating over thresholds

        # Path where all necessary pdbs are stored
        path = "/Users/macbook/Documents/perturbation-networks/proteins/parallel"
        # Path were resultings data csvs will be stored
        csv_path = "/Users/macbook/Documents/perturbation-networks/proteins/data"
        protein = "1nd4"
        threshold = 5

        AA = manager.list(list(Bio.PDB.Polypeptide.aa1)) # List of amino acids
        mutations_positions = it.MutationsList(path, protein)
        mutations = manager.list(mutations_positions[0])
        positions = manager.list(mutations_positions[1])
        N = len(AA) # Number of amino acids
        M = len(positions) # Number of mutated positions

        # Create arrays to store data
        # https://jonasteuwen.github.io/numpy/python/multiprocessing/2017/01/07/multiprocessing-numpy-array.html
        nodes1 = np.ctypeslib.as_ctypes(np.zeros((N, M)))
        edges1 = np.ctypeslib.as_ctypes(np.zeros((N, M)))
        weight1 = np.ctypeslib.as_ctypes(np.zeros((N, M)))
        distance1 = np.ctypeslib.as_ctypes(np.zeros((N, M)))
        #shared arrays
        nodes =  sharedctypes.RawArray(nodes1._type_, nodes1)
        edges =  sharedctypes.RawArray(edges1._type_, edges1)
        weight =  sharedctypes.RawArray(weight1._type_, weight1)
        distance =  sharedctypes.RawArray(distance1._type_, distance1)

        # save original data here

        # Generate molecule of original pdb file
        original_prot = it.Pmolecule(os.path.join(path, f"{protein}.pdb"))
        original = original_prot.network(cutoff=threshold)
        matrix = np.ctypeslib.as_ctypes(nx.adjacency_matrix(original).toarray())
        original_matrix = sharedctypes.RawArray(matrix._type_, matrix)

        # Get list of tuples to map in Pool
        tuples = [(protein, mutation, path, threshold) for mutation in mutations]

        with mp.Pool(4) as pool:
            pool.map(GetData, tuples)

        # Arrays should be full by now, we save as csvs in csv_path
        # positions will be used as header
        it.WriteCSV(csv_path, nodes, positions, f"{protein}_{threshold}_nodes.csv")
        it.WriteCSV(csv_path, edges, positions, f"{protein}_{threshold}_edges.csv")
        it.WriteCSV(csv_path, weight, positions, f"{protein}_{threshold}_weight.csv")
        it.WriteCSV(csv_path, distance, positions, f"{protein}_{threshold}_distance.csv")
