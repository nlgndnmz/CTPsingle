
## CTPsingle - clonality inference from low coverage single-sample tumors
###### author: Nilgun Donmez

CTPsingle is a tool that aims to infer the subclonal decomposition using low-coverage sequencing data from a single tumor sample. 

---- INSTALLATION

The core functionality of CTPsingle is implemented in R. If your system does not have R, you can install it for free-of-charge from https://www.r-project.org/. In addition, CTPsingle requires the following R packages:

- ggplot2
- RColorBrewer
- DPpackage
- lpSolve

The best way to install these packages is to call the install.packages() function in R. Additionally, CTPsingle comes with several utility scripts written in Python. To run these, you will need Python2.7 and the numpy library.

---- RUNNING CTPsingle

To run CTPsingle, you need a single file containing mutant read counts formatted as follows:

```

Chromosome Position Mutant Reference Mcount Rcount Multiplier Gender
8 116260205 C A 13 8 2 Female
5 106280608 T C 12 29 2 Female
9 62266747 G C 19 26 2 Female
13 64241231 G A 15 34 2 Female
9 34778288 C G 13 16 2 Female

```

Above, Chromosome and Position give the location of the mutation. Mutant and Reference columns give the variant and reference alleles at this location. Mcount and Rcount give the number of reads supporting the mutation and reference alleles respectively. Multiplier column contains the total copy count of the region containing the position in the tumor sample. Note that since CTPsingle does not handle amplified regions, this column should at most be 2. The Gender column should be one of Female, Male or unknown. If it is set to unknown, mutations on non-autosomal chromosomes are ignored.

In addition to a file formatted as above, you will need a folder containing tree structure files. This folder is supplied with the CTPsingle package and is called "GammaAdjMatrices". Typical usage of CTPsingle is as follows:

` Rscript ./CTPsingle.R -f <mutation_file> -o <output_prefix> -m <adjacency_matrices> `

where 'mutation_file' is a file formatted as described and 'output_prefix' is a prefix for CTPsingle's output. 'adjacency_matrices' should give the path of the GammaAdjMatrices folder. For your convenience, we also included a simulation dataset under './data'. To test CTPsingle, simply execute the following command from the installation directory:

` Rscript ./CTPsingle.R -f ./data/simulation15.frq -o ./data/simulation15 -m ./GammaAdjMatrices `

If execution is successful, the following files will be written under ./data:

```
simulation15.png
simulation15_cls.png
simulation15_frq.png
simulation15_cluster_assignments.txt
simulation15_num_3_tree_1.txt
simulation15_num_3_tree_2.txt
```

The simulation15_cluster_assignments.txt contains the assignment of the mutations to subclones. The 'mostLikely' column contains the ID of the subclone, while the 'averageFrequency' column gives the cancer cell fraction of that subcone. The figures visualize various aspects of the data. The '_num_X_Y.txt' files give possible tree topologies on X number of nodes (i.e. the number of subclones found). The last 4 columns of these files denote:

```
... [ parent-node ] [ child-node ] [ cancer cell fraction of child node ] [ objective score ]
```

which can be used the draw the tree structure. Node that the parent-node of the root node is always given as 0. A smaller objective score means that the data fits that tree structure better. In addition to these files, CTPsingle reports its progress and other statistics to the standard output. The ground truth information for this dataset is also given at ./data/groundTruth15.txt. 



