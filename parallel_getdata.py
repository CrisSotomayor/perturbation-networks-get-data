import os
import multiprocessing as mp
import iterate as it
import numpy as np
import networkx as nx
import Bio.PDB


def MultiprocessData(tuple):
    it.GetData(tuple[0], tuple[1], tuple[2], tuple[3])
    return

if __name__ == '__main__':
    # original_data matrix missing here, add when iterating over thresholds

    # Path where all necessary pdbs are stored
    path = "/Users/macbook/Documents/perturbation-networks/proteins/parallel"
    # Path were resultings data csvs will be stored
    csv_path = "/Users/macbook/Documents/perturbation-networks/proteins/data"

    protein = "1nd4"
    threshold = 5

    AA = list(Bio.PDB.Polypeptide.aa1) # List of amino acids
    mutations, positions = it.MutationsList(path, protein)
    N = len(AA) # Number of amino acids
    M = len(positions) # Number of mutated positions

    # Create arrays to store data
    nodes = np.zeros((N, M))
    edges = np.zeros((N, M))
    weight = np.zeros((N, M))
    distance = np.zeros((N, M))

    # save original data here

    # Generate molecule of original pdb file
    original_prot = it.Pmolecule(os.path.join(path, f"{protein}.pdb"))
    original = original_prot.network(cutoff=threshold)
    original_matrix = nx.adjacency_matrix(original).toarray()

    # Get list of tuples to map in Pool
    tuples = [(protein, mutation, path, threshold) for mutation in mutations]

    with mp.Pool(30) as pool:
        pool.map(MultiprocessData, tuples)

    # Arrays should be full by now, we save as csvs in csv_path
    # positions will be used as header
    WriteCSV(csv_path, nodes, positions, f"{protein}_{threshold}_nodes.csv")
    WriteCSV(csv_path, edges, positions, f"{protein}_{threshold}_edges.csv")
    WriteCSV(csv_path, weight, positions, f"{protein}_{threshold}_weight.csv")
    WriteCSV(csv_path, distance, positions, f"{protein}_{threshold}_distance.csv")
