import os
import Bio.PDB
import biographs as bg
import numpy as np
import pandas as pd
import networkx as nx
import csv


def GetData(path, prot, mutations, csv_path=None):
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
    # Generate molecule of original pdb file
    original_prot = bg.Pmolecule(os.path.join(path, f"{prot}.pdb"))
    # The range of thresholds will define the networks
    thresholds = [round(i, 1) for i in np.linspace(3, 10, 71)]
    # Create dir to save resulting csv files if not specified
    if csv_path is None:
        csv_path = os.path.join(path, "perturbation_network_data")
        os.makedirs(csv_path)
    # Check if path and csv_path exist
    assert os.path.exists(path), "Directory doesn't exist."
    assert os.path.exists(csv_path), "Directory doesn't exist."
    # Array to save data from original protein network
    original_data = np.array((4, len(thresholds)))

    # For each threshold we iterate over all mutations
    for i, threshold in enumerate(thresholds):
        M = len(mutations.keys())  # Number of mutated positions
        cols = list(mutations.keys())
        nodes = np.zeros((N, M))
        edges = np.zeros((N, M))
        weights = np.zeros((N, M))
        distance = np.zeros((N, M))

        # Generate network for original graph with threshold
        original = original_prot.network(cutoff=threshold)
        original_matrix = nx.adjacency_matrix(original).toarray()
        # Saving data from original network
        original_data[0][i] = len(original.nodes())
        original_data[1][i] = len(original.edges())
        original_data[2][i] = original.size(weight='weight')
        original_data[3][i] = nx.diameter(original)

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
                nodes[aa_index][index] = len(perturbation_network.nodes())
                edges[aa_index][index] = len(perturbation_network.edges())
                weights[aa_index][index] = perturbation_network.size(weight='weight')
                distance[aa_index][index] = nx.diameter(perturbation_network)

        # Save data arrays as csv files in csv_path
        with open(os.path.join(csv_path, f"{prot}_{threshold}_nodes.csv"),
         'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(cols)
            writer.writerows(nodes)
        with open(os.path.join(csv_path, f"{prot}_{threshold}_edges.csv"),
         'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(cols)
            writer.writerows(edges)
        with open(os.path.join(csv_path, f"{prot}_{threshold}_weights.csv"),
         'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(cols)
            writer.writerows(weights)
        with open(os.path.join(csv_path, f"{prot}_{threshold}_distance.csv"),
         'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(cols)
            writer.writerows(distance)

        # Save array from original data
        original_data = np.vstack([thresholds, original_data]) # add thresholds
        original_data = np.transpose(original_data) # to add names as header

        with open(os.path.join(csv_path, f"{prot}_original.csv"),
         'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(['threshold', 'nodes', 'edges', 'weights', 'distance'])
            writer.writerows(original_data)
    return
