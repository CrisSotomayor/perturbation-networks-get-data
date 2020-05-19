import os
import Bio.PDB
import biographs as bg
import numpy as np
import pandas as pd
import networkx as nx


def GetData(path, prot, mutations, csv_path = None):
    """Get data from amino acid mutation perturbation networks as CSV files.

    Parameters:
        path (string): path where original and mutated pdb files are located
        prot (string): name of original pdb file
        mutations (dict): keys are strings representing positions to mutate
                        (amino acid, chain and index), each contains a list of
                        mutations performed (original aa, chain, index and
                        mutated aa). Mutated pdb files should be found in path,
                        as prot_mutation.pdb according to mutations in said list
        csv_path (string): default None, path where CSV files will be saved,
                        if None, a dir named "perturbation_network_data" will
                        be created in 'path'

    Returns:
        None
    """
    # Sorted list of one letter amino acids
    AA = list(Bio.PDB.Polypeptide.aa1)
    N = len(AA)  # Number of amino acids
    # Generate network of original pdb file
    original_prot = bg.Pmolecule(os.path.join(path, f"{prot}.pdb"))
    original = original_prot.network()
    original_matrix = nx.adjacency_matrix(original).toarray()
    # The range of thresholds will define the networks
    thresholds = [round(i, 1) for i in np.linspace(3, 10, 71)]
    # Create dir to save resulting csv files if not specified
    # error si no existen?
    if csv_path == None:
        csv_path = os.path.join(path, "perturbation_network_data")
        os.makedirs(csv_path)

    # For each threshold we iterate over all mutations
    for threshold in thresholds:
        M = len(mutations.keys())  # Number of mutated positions
        cols = list(mutations.keys())
        nodes = np.zeros((N, M))
        edges = np.zeros((N, M))
        weights = np.zeros((N, M))
        distance = np.zeros((N, M))

        for index, position in enumerate(mutations.keys()):
            for mutation in mutations[position]:
                # Generate network for current mutation
                current_path = os.path.join(path, f"{prot}_{mutation}.pdb")
                current_prot = bg.Pmolecule(current_path)
                current = current_prot.network(cutoff=threshold)

                # Obtain the absolute difference in terms of adjacency
                # matrices: the perturbation network.
                current_matrix = nx.adjacency_matrix(current).toarray()
                difference = np.abs(original_matrix - current_matrix)
                perturbation_network = nx.from_numpy_array(difference)

                # Remove isolates for accurate perturbation network node count
                perturbation_network.remove_nodes_from(
                    list(nx.isolates(perturbation_network)))

                # Corresponding row in array according to mutation
                aa_index = AA.index(mutation[-1])

                # Information obtained from perturbation network
                nodes[aa_index][index]
                          = len(perturbation_network.nodes())
                edges[aa_index][index]
                          = len(perturbation_network.edges())
                weights[aa_index][index]
                          = perturbation_network.size(weight='weight')
                distance[aa_index][index]
                          = nx.diameter(perturbation_network)

        # Save data arrays as csv files in csv_path
        # Usamos DataFrames? creo que es mas rapido con np.savetxt pero esta
        # dificil ponerle los aminoacidos al inicio de las filas
        pd.DataFrame(nodes, columns=cols, index = AA).to_csv(
                    os.path.join(csv_path, f"{prot}_{threshold}_nodes.csv"))
        pd.DataFrame(edges, columns=cols, index = AA).to_csv(
                    os.path.join(csv_path, f"{prot}_{threshold}_edges.csv"))
        pd.DataFrame(distance, columns=cols, index = AA).to_csv(
                    os.path.join(csv_path, f"{prot}_{threshold}_distance.csv"))
        pd.DataFrame(weights, columns=cols, index = AA).to_csv(
                    os.path.join(csv_path, f"{prot}_{threshold}_weights.csv"))
    return
