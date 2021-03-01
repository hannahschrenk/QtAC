

---Package QtAC---

QtAC provides tools to analyze the maturation process of complex systems in the sense of Holling's adaptive cycle metaphor (see e.g. Holling 2001). Using time series of components' abundance data, the dynamics of information transfer among the components are estimated. On the basis of the resulting information networks, potential, connectedness, and resilience can be computed. The development of these three variables deﬁnes the system's course through the adaptive cycle. The package offers several options to visualize the results. The main function requires a transfer entropy estimator being included in the JIDT toolkit (Lizier 2014). Besides, functions of the packages pracma, rJava, igraph, and rgl are used.

ESSENTIAL SOFTWARE:

    R (https://www.r-project.org/)
    Perl (https://www.perl.org/get.html)
    install_packages_QtAC.R

INSTALLING PROCESS:

===(R Studio)===

    source("install_packages_QtAC.R")
    Tools -> Install Packages...

    -> Install from: Package Archive File -> Browse: Choose "QtAC_1.0.tar.gz" -> Install

    Type library(QtAC)

===Command===

	source("install_packages_QtAC.R")
        install.packages("QtAC_1.0.tar.gz", repos = NULL, type = "source")
	
        library(QtAC)

SUPPLEMENTARY MATERIAL:

The QtAC folder contains a manual "QtAC_0.1.0.pdf" describing the main functions. Required additional packages can be installed via "install_packages_QtAC.R". An example application to a bacterial community can be found in the subfolder "Example". The subfolder "dist" contains the file "MTinfodynamics.jar" which is necessary to run the main function QtAC. 
