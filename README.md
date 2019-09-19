# AOPEDF
## paper "Network-based Prediction of Drug-Target Interactions using an Arbitrary-Order Proximity Embedded Deep Forest"

### 'dataset' directory
- `drug_dict.txt`: list of drug unique identifier and drug names
- `protein_dict.txt`: list of protein unique identifier and protein names
- `disease_dict.txt`: list of disease unique identifier and disease names
- `se_dict.txt`: list of side effect unique identifier and side effect names
- `drugdrug.txt`: Drug-Drug interaction matrix
- `drugDisease.txt`: Drug-Disease association matrix
- `drugsideEffect.txt`: Drug-SideEffect association matrix
- `drugsim1network.txt`: Drug chemical similarity matrix
- `drugsim2network.txt`: Drug therapeutic similarity matrix
- `drugsim3network.txt`: Drug sequence similarity matrix
- `drugsim4network.txt`: Drug biological processes similarity matrix
- `drugsim5network.txt`: Drug cellular component similarity matrix
- `drugsim6network.txt`: Drug molecular function similarity matrix
- `proteinprotein.txt`: Protein-Protein interaction matrix
- `proteinDisease.txt`: Protein-Disease association matrix
- `proteinsim1network.txt`: Protein sequence similarity matrix
- `proteinsim2network.txt`: Protein biological processes similarity matrix
- `proteinsim3network.txt`: Protein cellular component similarity matrix
- `proteinsim4network.txt`: Protein molecular function similarity matrix
- `Sim_drugDisease`: Drug-Disease Jaccard similarity matrix
- `Sim_drugsideEffect.txt`: Drug-SideEffect Jaccard similarity matrix
- `Sim_proteinDisease.txt`: Protein-Disease Jaccard similarity matrix
>Noted: Since drug-disease network, drug-side-effect network and protein-disease network are heterogeneous network, we calculate the corresponding similarity networks based on the Jaccard similarity coefficient for them

### 'AROPE' directory
This directory contains code necessary to use AROPE to extract arbitrary-Order proximity from different network

### 'gcforest' directory
This directory contains library of deep forest classifier

### Basic Usage
```
$ python AOPEDF.py
```


