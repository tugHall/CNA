tugHall version 3.0
====================

[![License](https://img.shields.io/badge/License-GPLv3-orange.svg)](https://github.com/tugHall/Clone-based/blob/master/Documentation/LICENSE)

**tugHall** _(**tu**mor **g**ene-**Hall**mark)_ is a cancer-cell evolution model simulator, wherein gene mutations are linked to the hallmarks of cancer, 
which influence tumor cell behaviors. 
=======
**tugHall** _(**tu**mor **g**ene-**Hall**mark)_ is a cancer-cell evolution model simulator, wherein gene mutations are linked to the hallmarks of cancer,
which influence tumor cell behaviors.


This is an _**R**_-based script to simulate the cancer cell evolution in the framework of the model proposed by _**Prof. Mamoru Kato**_,
_Head of Bioinformatics Department, Research Institute, National Cancer Center, Tokyo, JAPAN_.

Authors and contributor list:
---
_**Iurii Nagornov**_

_**Jo Nishino**_

_**Mamoru Kato**_

_Division of Bioinformatics, Research Institute, National Cancer Center Japan, Tokyo, Japan_

All questions and requests can be sent to inagonov@ncc.go.jp

Project source can be downloaded from websites
---
https://github.com/tugHall/  -  the developing resource

Short description
---
The wide availability of recent cancer genomic data requires a coherent model that can sort out the relevant findings to systematically explain the clonal evolution and resultant intra-tumor heterogeneity (ITH). Here, we present a new mathematical model designed to computationally simulate the evolution of cancer cells. The model connects well-known cancer hallmarks with the specific mutational states of tumor-related genes. The cell behavior phenotypes are stochastically determined and the hallmarks interfere probabilistically with the phenotypic probabilities. In turn, the hallmark variables depend on the mutational states of tumor-related genes. Thus, it is expected our software can be used to deepen our understanding of cancer-cell evolution and generation of ITH.

How to use **tugHall** and how to analyze data, kindly see user-guides in **Documentation** folder in prefered format *Rmd, pdf, html*.


General changes
---

This version is based on the clone consideration instead cell-based version 1.1.
Each clone has one or more cells, that allows to accelerate the calculations when number of clones is much less than number of cells.
Definition of clone: the clone is set of cell with same set of genes, which have same mutated / not mutated sites in genes.

Changes in comparison with tugHall v.2.1 
---

This version can calculate the copy number alterations (CNA) caused by deletions and duplications in comparison with version 2.1. 
CNAs may malfunction genes and change variant allel frequences if the point mutations are located on CNAs.

=======
Changes in comparison with tugHall v.2.1
---

This version can calculate the copy number alterations (CNA) caused by deletions and duplications in comparison with version 2.1.
CNAs may malfunction genes and change variant allel frequences if the point mutations are located on CNAs.

Content of package
---

* **tugHall_3.0.R** is a R-script of simulation, which uses scripts in **/Code/** folder.
* **/Code/** is the directory with scripts to simulate and to analyze data. 
* **/Input/** and **/Output/** are the directories for input and output data during a simulation. 
* **/Documentation/** is a directory with documentation to use tugHall software and to analyze data of simulation.   
* **/Figures/** is the folder with figures of plots of last simulation, which are used in user-guides also. 
* **/TESTs/** is the directory with tests and it's results. Description of tests is in the folder. 


