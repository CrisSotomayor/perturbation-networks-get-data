from Bio import SeqIO
import os
import subprocess
import biographs as bg
import numpy as np
import pandas as pd
import networkx as nx



def GetData(path, prot, mutations):
    AA = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    
    #generar grafica de proteina original
    original_prot = bg.Pmolecule(path+prot+".pdb")
    original = original_prot.network()
    original_matrix = nx.adjacency_matrix(original).toarray() #para no calcularla cada vez

    #valores que consideraremos como umbrales, y un diccionario para guardar los datos, indexado por esos umbrales
    thresholds = [round(i, 1) for i in np.linspace(3, 10, 71)]
    data = {i:[] for i in thresholds}

    for i in thresholds: #por cada umbral iteramos sobre todas las mutaciones
        #dataframes donde guardamos los datos
        nodes = pd.DataFrame(np.zeros((20,len(mutations.keys()))), columns=list(mutations.keys()), index = AA)
        edges = pd.DataFrame(np.zeros((20,len(mutations.keys()))), columns=list(mutations.keys()), index = AA)
        weights = pd.DataFrame(np.zeros((20,len(mutations.keys()))), columns=list(mutations.keys()), index = AA)
        distance = pd.DataFrame(np.zeros((20,len(mutations.keys()))), columns=list(mutations.keys()), index = AA)

        for position in mutations.keys():
            for mutation in mutations[position]:
                #generar grafica de proteina correspondiente a esta mutacion con el umbral correspondiente
                current_path = path+prot+"_"+mutation+".pdb"
                current_prot = bg.Pmolecule(current_path)
                current = current_prot.network(cutoff = i)

                #obtenemos la diferencia con las matrices de adyacencia, porque la funcion no considera pesos
                current_matrix = nx.adjacency_matrix(current).toarray()
                difference = np.abs(original_matrix - current_matrix)

                #grafica resultante
                perturbation_network = nx.from_numpy_array(difference)

                #quitamos vertices sin vecinos (filas y columnas de solo ceros)
                perturbation_network.remove_nodes_from(list(nx.isolates(perturbation_network)))

                #informacion que queremos de la grafica
                nodes.at[mutation[-1], mutation[:-1]] = len(perturbation_network.nodes())
                edges.at[mutation[-1], mutation[:-1]] = len(perturbation_network.edges())
                weights.at[mutation[-1], mutation[:-1]] = perturbation_network.size(weight='weight')
                distance.at[mutation[-1], mutation[:-1]] = nx.diameter(perturbation_network)


        data[i] = [nodes, edges, weights, distance]
    return data
