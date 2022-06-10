import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sbn
from ete3 import NCBITaxa
import pandas as pd

ncbi = NCBITaxa()

graphpath = '/home/tholstei/repos/PepGM_all/PepGM/results/Adeno_2_auto_taxa/PXD004095_Adenovirus2/refseqViral_PepGM_graph.graphml'
outmax ='/home/tholstei/repos/PepGM_all/PepGM/results/Adeno_2_auto_taxa/PXD004095_Adenovirus2/peptidome_sim_MAX.png'
#outmin = '/home/tholstei/repos/PepGM_all/PepGM/results/Cowpox_nohostfilter/PXD003013_Cowpox_BR/peptidome_sim_MIN.png'
taxonresults = '/home/tholstei/repos/PepGM_all/PepGM/results/Adeno_2_auto_taxa/PXD004095_Adenovirus2/PepGm_Results.csv'


#scripts that compute the peptidome similarity between taxa using only the peptides that were included in the PepGm graph
def GetNhighestTaxa(resultscsv,N):
    IDs = pd.read_csv(resultscsv, names = ['ID','score','type'])
    return IDs.ID.to_list()[0:N]

def GetPeptidesperTaxon(Graphin,Taxa):
    graph = nx.read_graphml(Graphin)
    PeptidomeDict = {}
    for node in graph.nodes(data=True):
        if node[1]['category']=='taxon' and node[0] in Taxa:
            neighbors = graph.neighbors(node[0])
            PeptidomeDict.update({node[0]:[n for n in neighbors]})

    return PeptidomeDict

def ComputeDetectedPeptidomeSimilarity(PeptidomeDict):
    SimMatrixMax = []
    SimMatrixMin = []
    Taxa1 = []
    Taxa2 = []
    for taxon1 in PeptidomeDict.keys():
        Taxa1.append(taxon1)
        SimMatrixMaxRow = []
        SimMatrixMinRow = []
        for taxon2 in PeptidomeDict.keys():
            Taxa2.append(taxon2)
            peptides1 = set(PeptidomeDict[taxon1])
            peptides2 = set(PeptidomeDict[taxon2])
            shared = len(peptides1.intersection(peptides2))
            try:
                SimMax = shared/(max(len(peptides1),len(peptides2)))
            except:
                SimMax = 0
            try:
                SimMin = SimMin = shared/(min(len(peptides1),len(peptides2)))
            except:
                SimMin = 0
            
            SimMatrixMaxRow.append(SimMax)
            SimMatrixMinRow.append(SimMin)
        SimMatrixMax.append(SimMatrixMaxRow)
        SimMatrixMin.append(SimMatrixMinRow)

    return SimMatrixMax,SimMatrixMin,Taxa1,Taxa2


Taxa = GetNhighestTaxa(taxonresults,15)
PeptiDict = GetPeptidesperTaxon(graphpath,Taxa)
SimMax,SimMin,Taxa1,Taxa2 = ComputeDetectedPeptidomeSimilarity(PeptiDict)
Taxid1 = ncbi.get_taxid_translator(Taxa1)
TaxonList1 = [Taxid1[int(tax)] for tax in Taxa1]


sbn.set(rc={'figure.figsize':(13.7,10.27)})

ax1 = sbn.heatmap(SimMax, xticklabels = TaxonList1, yticklabels = TaxonList1,cmap="YlGnBu",annot = [[round(n,2) for n in inner_list] for inner_list in SimMax])
plt.tight_layout()
plt.savefig(outmax)
plt.close()

#ax1 = sbn.heatmap(SimMin, xticklabels = TaxonList, yticklabels = TaxonList,cmap="YlGnBu")
#plt.tight_layout()
#plt.savefig(outmin)
#plt.close()




        
