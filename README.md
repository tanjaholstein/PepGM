# PepGM
for Felix

A snakemake workflow 

 Input: mass spectrum mgf or mzML files, searchGUI parameter file, config file (DB used, target taxa, host)
Output: Organism/ taxonomic classification with a confidence score, visualization of input/output data as phylogenetic tree and barplots<\p>

 The graphical model used was inspired by https://www.openms.de/comp/epifany/ 
 
 
Workflow description:
1. searchDB cleanup : cRaP DB ist added, host is added (if wanted), duplicate entries are removed using seqkit. generation of target-decoy DB using searchCLI
2. peptide search using searchCLI + PeptideShaker. Generation of a a peptide list
3. all descendant strain of the target taxa are queried in the NCBI protein DB (possibility to filter swissprot only/all/refseq onlyetc) through the NCBI API. scripts: CreatePepGMGraph.py and FactorGraphGeneration.py
4. Donwloaded protein recordes are digested using cp-dt and queried again the protein ID list to generate a bipartite taxon-peptide graph scripts: CreatePepGMGraph.py and FactorGraphGeneration.py
5. The bipartite graph is transformed into a factor graph using convolution trees and conditional probability table factors (CPD). scripts: CreatePepGMGraph.py and FactorGraphGeneration.py
6. for different sets of CPD parameters, the belief propagation algorithm is run until convergence to obtain the posterior probabilites of the taxa. scripts: belief_propagation.py and PepGM.py 
7. Through an  empirically deduced metric, the ideal parameter set is inferred. script GridSearchAnalysis.py 
8. For this ideal parameter set, we output a results barchart and phylogenetic tree view showcasing the 15 best scoring tax. scripts: BarPlotResults, PhyloTreeView.py
