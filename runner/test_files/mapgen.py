import pandas as pd
import numpy as np
import argparse
from rdkit.Chem.Fingerprints import FingerprintMols
from rdkit.Chem import AllChem
from rdkit import Chem
import networkx as nx
import os

class MapGen():
    def __init__(self, in_sdf, metric):
        self.suppl = Chem.SDMolSupplier(in_sdf)
        self.metric = metric.lower()

    def make_fp(self, mol):
        if self.metric == 'tanimoto':
            fp = FingerprintMols.FingerprintMol(mol)
        if self.metric == 'mfp':
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
        if self.metric == 'smiles':
            fp = Chem.MolToSmiles(mol, isomericSmiles=True)
        return fp

    def get_ligdict(self): #1
        lig_dict = {}
        for mol in self.suppl:
            charge = Chem.rdmolops.GetFormalCharge(mol)
            if charge not in lig_dict.keys():
                lig_dict[charge] = {'Name':[], 'Mol':[], 'FP':[]}
                lig_dict[charge]['Name'].append(mol.GetProp('_Name'))
                lig_dict[charge]['Mol'].append(mol)
                if self.metric != 'mcs':
                    lig_dict[charge]['FP'].append(self.make_fp(mol))
            else:
                lig_dict[charge]['Name'].append(mol.GetProp('_Name'))
                lig_dict[charge]['Mol'].append(mol)
                if self.metric != 'mcs':
                    lig_dict[charge]['FP'].append(self.make_fp(mol))
            
        self.lig_dict = lig_dict
        
    def sim_mx(self): #2
        if self.metric == 'tanimoto' or self.metric == 'mfp': 
            from rdkit import DataStructs
        elif self.metric == 'mcs':
            from rdkit.Chem import rdFMCS
        elif self.metric == 'smiles':
            from Bio import pairwise2
        sim_dfs = {}
        for charge in self.lig_dict.keys():
            df = pd.DataFrame()
            for index, i in enumerate(self.lig_dict[charge]['Name']):
                for jndex, j in enumerate(self.lig_dict[charge]['Name']):
                    if i == j:
                        df.loc[i, j] = 1.0
                        continue
                    elif i != j:
                        if self.metric == 'tanimoto' or self.metric == 'mfp': 
                            fp1,fp2 = self.lig_dict[charge]['FP'][index], self.lig_dict[charge]['FP'][jndex]
                            similarity = DataStructs.FingerprintSimilarity(fp1,fp2)
                            df.loc[i, j] = similarity
                        if self.metric == 'mcs':
                            mcs = rdFMCS.FindMCS([self.lig_dict[charge]['Mol'][index], self.lig_dict[charge]['Mol'][jndex]], atomCompare=rdFMCS.AtomCompare.CompareAny)
                            df.loc[i, j] = mcs.numAtoms + mcs.numBonds
                        if self.metric == 'smiles':
                            fp1,fp2 = self.lig_dict[charge]['FP'][index], self.lig_dict[charge]['FP'][jndex]
                            alignments = pairwise2.align.globalms(fp1, fp2, 1,-1, -0.5, -0.05)
                            df.loc[i, j] = alignments[0][2]

            sim_dfs[charge] = df

        self.sim_dfs = sim_dfs

    def clean_mxs(self):
        for charge in self.lig_dict.keys():
            if self.metric == 'tanimoto' or self.metric == 'mfp' or self.metric == 'smiles':
                m = self.sim_dfs[charge].idxmax()
                for i, j in zip(self.sim_dfs[charge].index, m):
                    vlist = self.sim_dfs[charge].loc[i, :].tolist()
                    index = vlist.index(1.0)
                    vlist[index:] =  [0.0] * len(vlist[index:])
                    self.sim_dfs[charge].loc[i, :] = vlist

                if self.metric == 'tanimoto' or self.metric == 'mfp':
                    self.sim_dfs[charge] = 1 - self.sim_dfs[charge] #get dissimilarity matrix
                self.sim_dfs[charge] = self.sim_dfs[charge].replace(1.0, 0.0) #set diagonal to 0
                self.sim_dfs[charge] = self.sim_dfs[charge].replace(0.0, 1.0) #set zeroes into 1 (in order to search shortest path)

            if self.metric == 'mcs':
                m = self.sim_dfs[charge].idxmax()
                for i, j in zip(self.sim_dfs[charge].index, m):
                    vlist = self.sim_dfs[charge].loc[i, :].tolist()
                    index = vlist.index(1.0)
                    vlist[index:] =  [0.0] * len(vlist[index:])
                    self.sim_dfs[charge].loc[i, :] = vlist
                    
                self.sim_dfs[charge] = 100 - self.sim_dfs[charge] #get dissimilarity matrix
                self.sim_dfs[charge] = self.sim_dfs[charge].replace(1.0, 0.0) #set diagonal to 0
                self.sim_dfs[charge] = self.sim_dfs[charge].replace(0.0, 1.0) #set zeroes into 1 (in order to search shortest path)

    def get_ligpairs(self):
        for charge in self.lig_dict.keys():
            pairs_dict = {}
            for i in self.sim_dfs[charge].index:
                for j in self.sim_dfs[charge].columns:
                    if self.sim_dfs[charge].loc[i, j] == 1.0:
                        continue
                    else:
                        pairs_dict['{} {}'.format(i, j)] = round(self.sim_dfs[charge].loc[i, j], 3)
            
            if self.metric == 'tanimoto' or self.metric == 'mfp' or self.metric == 'mcs':
                pairs_dict = {k: v for k, v in sorted(pairs_dict.items(), key=lambda item: item[1], reverse=False)}
            elif self.metric == 'smiles': 
                pairs_dict = {k: v for k, v in sorted(pairs_dict.items(), key=lambda item: item[1], reverse=True)}

            self.lig_dict[charge]['pairs_dict'] = pairs_dict

    def intersection(self, edge_list, candidate_edge):
        k = False
        r1, r2 = candidate_edge.split()[0], candidate_edge.split()[1]
        for edge in edge_list:
            if r1 == edge[0] or r1 == edge[1] or r2 == edge[0] or r2 == edge[1]:
                k = True
        return k

    def not_ingraph(self, node_list, candidate_edge):
        k = False
        r1, r2 = candidate_edge.split()[0], candidate_edge.split()[1]
        if r1 not in node_list or r2 not in node_list:
            k = True
        return k

    def outer_nodes(self, G):
        node_list = []
        for node in G.nodes:
            if len(G.edges(node)) == 1:
                node_list.append(node)
        return node_list     

    def make_map(self):
        for charge in self.lig_dict.keys():
            H = nx.Graph()
            if len(self.lig_dict[charge]['Name']) == 1: #In case one ligand is found alone in a charge group
                ligcol = self.sim_dfs[self.lig_dict[charge]['Name']].sort_values(by=[self.lig_dict[charge]['Name'][0]]) #complete similarity matrix
                H.add_edge(self.lig_dict[charge]['Name'][0], ligcol.index[1])
                H.add_edge(self.lig_dict[charge]['Name'][0], ligcol.index[2])
                lig_dict[charge]['Graph'] = H
                break

            # 1. Make SPT
            incomplete = True
            while incomplete:
                for pert, score in self.lig_dict[charge]['pairs_dict'].items():
                    if len(H.nodes) == len(self.lig_dict[charge]['Name']):
                        incomplete = False
                        break
                    l1, l2 = pert.split()[0], pert.split()[1]
                    if H.has_edge(l1, l2) or H.has_edge(l2, l1):
                        continue
                    if len(H.nodes) == 0 or self.intersection(H.edges, pert) and self.not_ingraph(H.nodes, pert):
                            H.add_edge(l1, l2, weight=score)
                            break

            # 2. Close Cycles            
            while len(self.outer_nodes(H)) != 0:
                for pert, score in self.lig_dict[charge]['pairs_dict'].items():
                    l1, l2 = pert.split()[0], pert.split()[1]
                    if l1 in self.outer_nodes(H) or l2 in self.outer_nodes(H):
                        if (l1, l2) not in H.edges or (l2, l1) not in H.edges:
                            H.add_edge(l1, l2, weight=score)
                            break
                        else:
                            continue
                    else:
                        continue

            # 3. Add influence edges
            eig_cent = nx.eigenvector_centrality(H, max_iter=1000)
            eig_cent = {k: v for k, v in sorted(eig_cent.items(), key=lambda item: item[1], reverse=True)}
            cent_nodes = [k for k,v in nx.eigenvector_centrality(H, max_iter=1000).items() if v > 0.15]
            per_nodes = [k for k,v in nx.eigenvector_centrality(H, max_iter=1000).items() if v < 0.01]
            per_len = len(per_nodes)
            while per_len > 1:
                for pert, score in self.lig_dict[charge]['pairs_dict'].items():
                    l1, l2 = pert.split()[0], pert.split()[1]
                    if l1 in per_nodes and l2 not in per_nodes or l1 not in per_nodes and l2 in per_nodes:
                        if (l1, l2) not in H.edges or (l2, l1) not in H.edges and intersection(H.edges, pert):
                            H.add_edge(l1, l2, weight=score)
                            nlen = len([v for k,v in nx.eigenvector_centrality(H, max_iter=1000).items() if v < 0.01])
                            per_len = nlen
                            break
                        else:
                            continue
                    else:
                        continue

            self.lig_dict[charge]['Graph'] = H
            for edge in H.edges:
                print('{} {}'.format(edge[0], edge[1]))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog='MapGen',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description='FEP map generator based on selected distance metrics.')
    parser.add_argument('-isdf', '--insdf',
                        dest="isdf",
                        required=True,
                        help=".sdf file name")
    parser.add_argument('-m', '--metric',
                        dest="metric",
                        default='MFP',
                        choices=['MFP', 'Tanimoto', 'MCS', 'SMILES'],
                        required=False,
                        help="Distance metric for ligand pairwairse comparison")

args = parser.parse_args()
mg = MapGen(args.isdf, args.metric)
mg.get_ligdict()
mg.sim_mx()
mg.clean_mxs()
mg.get_ligpairs()
mg.make_map()