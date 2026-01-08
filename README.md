# MORE
MORE (Multi-Omics REgulation) is an R package for the application of Partial Least Squares (PLS) or Multiple Linear Regression (MLR) models with Elastic Net or Iterative Sparse Group Lasso (ISGL) regularizations to multi-omics data. The MORE method applies MLRs or PLS to model a target omic expression as a function of experimental variables, such as diseases or treatments, and the potential regulators of that given target feature. The aim is to obtain specific candidate regulators for the biological system under study.

## Installation

### Prerequisites

Before installation, please ensure you have **Rtools** installed on your system (for Windows users), as it is required for building and running certain dependencies. You can download it from [CRAN](https://cran.r-project.org/bin/windows/Rtools/).

### Installing MORE

Currently, the package can be installed directly from GitHub using the `devtools` R package:

    install.packages("devtools")
    devtools::install_github("BiostatOmics/MORE")


Although all dependencies should install automatically, you can manually install them if the process fails:

* glmnet
* igraph
* MASS
* psych
* car
* ltm
* graph
* ropls
* clusterProfiler

Note that the last three dependencies are Bioconductor packages and must be installed using the BiocManager helper: 

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("graph", "ropls", "clusterProfiler"))
```

## Usage

You can find the User Guide for the package in the vignettes folder or access it directly [here](https://github.com/BiostatOmics/MORE/blob/master/vignettes/UsersGuide.pdf). Furthermore, in the same folder we include a full step-by-step tutorial to run MORE (accessible [here](https://github.com/BiostatOmics/MORE/blob/master/vignettes/tutorial.html))

