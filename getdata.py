from Bio import SeqIO
import os
import subprocess
import Bio.PDB
import biographs as bg
import numpy as np
import pandas as pd
import networkx as nx


def GetData(path, prot, mutations):
    """
    Missing docstring
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
    data = {threshold: [] for threshold in thresholds}

    # For each threshold we iterate over all mutations
    for threshold in thresholds:
        M = len(mutations.keys())  # Number of mutated positions
        cols = list(mutations.keys())
        nodes = pd.DataFrame(np.zeros((N, M)),
                             columns=cols, index=AA)
        edges = pd.DataFrame(np.zeros((N, M)),
                             columns=cols, index=AA)
        weights = pd.DataFrame(np.zeros((N, M)),
                               columns=cols, index=AA)
        distance = pd.DataFrame(np.zeros((N, M)),
                                columns=cols, index=AA)

        for position in mutations.keys():
            for mutation in mutations[position]:
                current_path = os.path.join(path, f"{prot}_{mutation}.pdb")
                current_prot = bg.Pmolecule(current_path)
                current = current_prot.network(cutoff=threshold)

                # Obtain the absolute difference in terms of adjacency
                # matrices: the perturbation network.
                current_matrix = nx.adjacency_matrix(current).toarray()
                difference = np.abs(original_matrix - current_matrix)
                perturbation_network = nx.from_numpy_array(difference)

                # quitamos vertices sin vecinos (filas y columnas de solo ceros)
                # Do we really wan to do this?
                perturbation_network.remove_nodes_from(
                    list(nx.isolates(perturbation_network)))

                # informacion que queremos de la grafica
                nodes.at[mutation[-1], mutation[:-1]
                         ] = len(perturbation_network.nodes())
                edges.at[mutation[-1], mutation[:-1]
                         ] = len(perturbation_network.edges())
                weights.at[mutation[-1], mutation[:-1]
                           ] = perturbation_network.size(weight='weight')
                distance.at[mutation[-1], mutation[:-1]
                            ] = nx.diameter(perturbation_network)

        data[threshold] = [nodes, edges, weights, distance]
    return data
