# Installation of phyloscanner_analyse_trees.R
Phyloscanner's analysis of phylogenies with within- and between-host diversity needs R, ideally up-to-date.

Once R is installed, change directory to the "phyloscannerR"
```sh
cd PhyloscannerExtra/phyloscannerR/ 
```
then start an an interactive R session by running R
```sh
R
```
Then inside the interactive R session, run:

install.packages("devtools")
install.packages("BiocManager")
library(devtools)
BiocManager::install("ggtree")
BiocManager::install("RBGL")
install("../phyloscannerR", dependencies = T)

NB! Unfortunately, we never have an up-to-date version of R on our HPC. So, I installed "ggtree", using ggtree_3.8.0.tar.gz which is compatible with our R version (v 4.1.3 ). 
