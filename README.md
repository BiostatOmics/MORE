# MORE
MORE (Multi-Omics REgulation) is an R package for the application of Partial Least Squares (PLS) or Multiple Linear Regression (MLR) models with Elastic Net or Iterative Sparse Group Lasso (ISGL) regularizations to multi-omics data. The MORE method applies MLRs or PLS to model a target omic expression as a function of experimental variables, such as diseases or treatments, and the potential regulators of that given target feature. The aim is to obtain specific candidate regulators for the biological system under study.


## Installing

Currently, the package can be installed directly from GitHub using the `devtools` R package:

    install.packages("devtools")
    devtools::install_github("BiostatOmics/MORE")

Before installation, it might be necessary to install the required dependencies:

* glmnet
* igraph
* MASS
* psych
* car
* ltm
* ropls

## Usage

You can find the User Guide for the package in the vignettes folder or access it directly [here](https://github.com/BiostatOmics/MORE/blob/master/vignettes/UsersGuide.pdf).

