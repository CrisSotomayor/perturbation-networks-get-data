import os
import multiprocessing as mp
from getdata import GetData
from getmutations import MutationsDict


def Multiprocess(path):
    prot = path.split('/')[-1]
    mutations = MutationsDict(path, prot)
    GetData(path, prot, mutations)

    return

if __name__ == '__main__':
    proteins = ["3s4y", "1nd4", "6r5k", "3dqw"]
    main_path = "/Users/macbook/Documents/perturbation-networks/proteins/"
    paths = [os.path.join(main_path, protein) for protein in proteins]

    with mp.Pool(40) as pool:
        pool.map(Multiprocess, paths)
