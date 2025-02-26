#' Data to showcase MORE usage
#'
#' @format List containing data of Gene expression of differentially expressed genes, data from different omic regulators, the potential associations between genes and regulators and the experimental design under study.
#' 
#' \describe{
#' 
#'    \item{GeneExpressionDE}{Matrix with genes in rows and samples in columns. It contains information about the gene expression of the differentially expressed genes.}
#'    \item{data.omics}{List where each element contains a matrix or a data.frame containing data for each regulatory omic. Each element must have a structure similar to gene expression data. Regulators in rows and samples in columns.}
#'    \item{edesign}{Matrix or data.frame containing the experimental conditions. Rows must be same as columns in GeneExpressionDE and in the same order.}
#'    \item{associations}{If prior biological knowledge want to be considered, list containing as many elements as omics considered in data.omics and in the same order. Each element of the list must be a data.frame with genes in the first column and regulators that potentially regulate them in the second column. If there is no biological prior knowledge it should be set equal to NULL.}
#' }
#'
#' @source {Data modified from STATegra datasets to serve as an example}
#' @usage data(TestData)

"TestData"

#' Data to showcase MORE usage
#'
#' @format List containing data of Gene expression of differentially expressed genes, data from different omic regulators, the potential associations between genes and regulators and the experimental design under study.
#' 
#' \describe{
#' 
#'    \item{GeneExpressionDE}{Matrix with genes in rows and samples in columns. It contains information about the gene expression of the differentially expressed genes.}
#'    \item{data.omics}{List where each element contains a matrix or a data.frame containing data for each regulatory omic. Each element must have a structure similar to gene expression data. Regulators in rows and samples in columns.}
#'    \item{edesign}{Matrix or data.frame containing the experimental conditions. Rows must be same as columns in GeneExpressionDE and in the same order.}
#'    \item{associations}{If prior biological knowledge want to be considered, list containing as many elements as omics considered in data.omics and in the same order. Each element of the list must be a data.frame with genes in the first column and regulators that potentially regulate them in the second column. If there is no biological prior knowledge it should be set equal to NULL.}
#' }
#'
#' @source {Data modified from STATegra datasets to serve as an example, extended version of TestData example}
#' @usage data(STATegraReduced)

"STATegraReduced"

