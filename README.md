### Package QtAC

QtAC provides tools to analyze the maturation process of complex systems in the sense of Holling's adaptive cycle metaphor ([zu Castell and Schrenk, 2020](https://www.nature.com/articles/s41598-020-74888-y)). Using time series of components' abundance data, the dynamics of information transfer among the components are estimated. On the basis of the resulting information networks, potential, connectedness, and resilience can be computed. The development of these three variables deﬁnes the system's course through the adaptive cycle. The package offers several options to visualize the results. The main function requires a transfer entropy estimator being included in the JIDT toolkit ([Lizier, 2014](https://www.frontiersin.org/articles/10.3389/frobt.2014.00011/full)). Besides, functions of the packages pracma, rJava, igraph, and rgl are used.

ESSENTIAL SOFTWARE:

* R >= 3.6.1 (https://www.r-project.org/)
* install_packages_QtAC.R

INSTALLING PROCESS:

__R Studio:__
```R
source("install_packages_QtAC.R")
```
Tools -> Install Packages...
-> Install from: Package Archive File -> Browse: Choose "QtAC_1.0.tar.gz" -> Install
```R
library(QtAC)
```
__Command:__
```R
source("install_packages_QtAC.R")
install.packages("QtAC_1.0.tar.gz", repos = NULL, type = "source")
library(QtAC)
```
EXAMPLE DESCRIPTION:

We provide an example application of QtAC to a SIHUMI community (simplified human intestinal microbiota, Becker et al. 2010). It is based on time series of abundance data of seven bacterial species, which has been simulated using BacArena (Bauer et al. 2017). Each row of the tab-separated file "QtAC_SIHUMI.txt" contains the number of individuals of one of the species at 40 subsequent simulation steps. The name of the respective species is specified in the first column. Following "QtAC_SIHUMI.R", the user exemplarily learns how to estimate a series of information networks from abundance data and how to compute the corresponding systemic variables. Several ways of visualizing the results are illustrated. 

RUN EXAMPLE:
* Access folder "Example"
* Download "QtAC_SIHUMI.txt"
* Follow instructions in "QtAC_SIHUMI.R"

SUPPLEMENTARY MATERIAL:

The QtAC folder contains a manual "QtAC_0.1.0.pdf" describing the main functions. Additionally required packages can be installed via "install_packages_QtAC.R". An example application to a bacterial community can be found in the subfolder "Example". The subfolder "dist" contains the file "MTinfodynamics.jar" which is necessary to run the main function QtAC. 
