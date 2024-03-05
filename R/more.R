#########################################################################################
######           MORE:      Multi-Omics REgulation                        ######
#########################################################################################


## By Sonia, Monica & Maider
## 15-October-2023
## Last modified: August 2023


options(stringsAsFactors = FALSE)

library(igraph)
library(MASS)
library(glmnet)
library(psych)
library(car)
library(stringr)
library(RColorConesa)

setClass("MORE")

isBin <-function(x){
  if(length(unique(x[!is.na(x[, 1]), 1]))==2 && length(unique(x[1,!is.na(x[1,]),drop=TRUE]))==2){
    return(1)
  }
  else{
    return(0)
  }
}

isBinclinic <-function(x){
  
  if(class(x)=='character'){
    return(1)
  }
  else if(length(unique(x))==2){
    return(1)
  }
  else{
    return(0)
  }
}

#' more: Multi-Omics Regulation
#'
#' \code{more} fits a GLM regression model (when the selected method is GLM) or a PLS model (when the selected method is PLS) for all genes in the dataset to identify
#' the potential regulators that show a significant impact on gene expression under specific experimental conditions.
#'
#' @param GeneExpression Data frame containing gene expression data with genes in rows and
#' experimental samples in columns. Row names must be the gene IDs.
#' @param data.omics List where each element corresponds to a different omic data type to be considered (miRNAs,
#' transcription factors, methylation, etc.). The names of the list will represent the omics, and each element in 
#' the list should be a data matrix with omic regulators in rows and samples in columns.
#' @param associations List where each element corresponds to a different omic data type (miRNAs,
#' transcription factors, methylation, etc.). The names of the list will represent the omics. Each element in 
#' the list should be a data frame with 2 columns (optionally 3), describing the potential interactions between genes
#' and regulators for that omic. First column must contain the genes (or features in
#' GeneExpression object), second column must contain the regulators, and an optional third column can
#' be added to describe the type of interaction (e.g., for methylation, if a CpG site is located in
#' the promoter region of the gene, in the first exon, etc.). If the user lacks prior knowledge of the potential regulators, they can set the parameter to NULL. 
#' In this case, all regulators in \code{\link{data.omics}} will be treated as potential regulators for all genes. In this case, for computational efficiency, it is recommended to use pls2 \code{\link{method}}.
#' Additionally, if the users have prior knowledge for certain omics and want to set other omics to NULL, they can do so.
#' @param edesign Data frame describing the experimental design. Rows must be the samples (columns
#' in \code{\link{GeneExpression}}) and columns must be the experimental variables to be included in the model (e.g. treatment, etc.).
#' @param clinic Data.frame with all clinical variables to consider,with samples in rows and variables in columns.
#' @param clinic.type Vector which indicates the type of data of variables introduced in \code{\link{clinic}}. The user should code as 0 numeric variables and as 1 categorical or binary variables. 
#' By default is set to NULL. In this case, the data type will be predicted automatically. However, the user must verify the prediction and manually input the vector if incorrect.
#' @param center By default TRUE. It determines whether centering is applied to \code{\link{data.omics}}.
#' @param scale By default TRUE. It determines whether scaling is applied to \code{\link{data.omics}}.
#' @param scaletype Type of scaling to be applied. Three options:
#' \itemize{
#' \item auto : Applies the autoscaling. 
#' \item pareto : Applies the pareto scaling. \deqn{\frac{X_k}{s_k \sqrt[4]{m_b}} }
#' \item block : Applies the block scaling. \deqn{ \frac{X_k}{s_k \sqrt{m_b}} }
#' }
#' considering m_b the number of variables of the block. By default, auto.
#' @param epsilon Convergence threshold for coordinate descent algorithm in elasticnet. Default value, 1e-5.
#' @param min.variation  For numerical regulators, it specifies the minimum change required across conditions to retain the regulator in 
#' the regression models. In the case of binary regulators, if the proportion of the most common value is equal to or inferior this value, 
#' the regulator is considered to have low variation and will be excluded from the regression models. The user has the option to set a single 
#' value to apply the same filter to all omics, provide a vector of the same length as omics if they want to specify different levels for each omics, 
#' or use 'NA' when they want to apply a minimum variation filter but are uncertain about the threshold. By default, 0.
#' @param interactions.reg If TRUE, the model includes interactions between regulators and experimental variables. By default, TRUE.
#' @param family.glm Error distribution and link function to be used in the model when \code{\link{method}} glm. By default, gaussian().
#' @param elasticnet.glm ElasticNet mixing parameter. There are three options:
#' \itemize{
#' \item NULL : The parameter is selected from a grid of values ranging from 0 to 1 with 0.1 increments. The chosen value optimizes the mean cross-validated error when optimizing the lambda values.
#' \item A number between 0 and 1 : ElasticNet is applied with this number being the combination between Ridge and Lasso penalization (elasticnet=0 is the ridge penalty, elasticnet=1 is the lasso penalty). 
#' \item A vector with the mixing parameters to try. The one that optimizes the mean cross-validated error when optimizing the lambda values will be used.
#' }
#' By default, NULL.
#' @param col.filter.glm Type of correlation coefficients to use when applying the multicollinearity filter when glm \code{\link{method}} is used. 
#' \itemize{
#' \item cor: Computes the correlation between omics. Pearson correlation between numeric variables, phi coefficient between numeric and binary and biserial correlation between binary variables. 
#' \item pcor : Computes the partial correlation.
#' }
#' @param correlation.glm  Value to determine the presence of collinearity between two regulators when using the glm \code{\link{method}}. By default, 0.7.
#' @param gr.method.isgl Grouping approach to create groups of variables in ISGL penalization. There are two options: 'cor' to cluster variables using correlations and 'pca' to use Principal Component Analysis approach. By default, 'cor'.
#' @param thres.isgl Threshold for the correlation when gr.method.isgl is 'cor' or threshold for the percentage of variability to explain when 'pca'. By default, 0.7.
#' @param alfa.pls Significance level for variable selection in pls1and pls2 \code{\link{method}}. By default, 0.05.
#' @param p.method.pls Type of resampling method to apply for the p-value calculation when pls1 or pls2 \code{\link{method}}. Two options:
#' \itemize{
#' \item jack : Applies Jack-Knife resampling technique.
#' \item perm : Applies a resampling technique in which the response variable is permuted 100 times to obtain the distribution of the coefficients and compute then their associated p-value.
#' }
#' By default, jack.
#' @param vip.pls Value of VIP above which a variable can be considered significant in addition to the computed p-value in \code{\link{p.method}}. By default, 0.8.
#' @param method Model to be fitted. Four options:
#' \itemize{
#' \item glm : Applies a Generalized Linear Model (GLM) with ElasticNet regularization.
#' \item pls1 : Applies a Partial Least Squares (PLS) model, one for each of the genes at \code{\link{GeneExpression}}.
#' \item pls2 : Applies a PLS model to all genes at the same time, only possible when \code{\link{associations}}= NULL.
#' \item isgl : Applies a Generalized Linear Model (GLM) with Iterative Sparse Group Lasso (ISGL) regularization.
#' }
#' By default, glm.
#' @return List containing the following elements:
#' \itemize{
#' \item ResultsPerGene : List with as many elements as genes in \code{\link{GeneExpression}}. For each gene, it includes information about gene values, considered variables, estimated coefficients,
#'                    detailed information about all regulators, and regulators identified as relevant (in glm scenario) or significant (in pls scenarios).
#' \item GlobalSummary : List with information about the fitted models, including model metrics, information about regulators, genes without models, regulators, master regulators and hub genes.
#' \item Arguments : List containing all the arguments used to generate the models.                
#' }
#'
#' @examples
#' 
#' data(TestData)
#' 
#' #Omic type
#' omic.type = c(1,0,0)
#' names(omic.type) = names(TestData$data.omics)
#' SimGLM = more(GeneExpression = TestData$GeneExpressionDE,
#'               associations = TestData$associations, 
#'               data.omics = TestData$data.omics,
#'               omic.type = omic.type,
#'               edesign = TestData$edesign,
#'               center = TRUE, scale = TRUE, 
#'               scaltype = 'auto',
#'               epsilon = 0.00001, family.glm = gaussian(), elasticnet = NULL,
#'               interactions.reg = TRUE,min.variation = 0,  col.filter.glm = 'cor',
#'               correlation.glm = 0.7, method  ='glm')
#' 
#' @export

more <-function(GeneExpression,
                data.omics,
                associations = NULL,
                omic.type = NULL,
                edesign = NULL,
                clinic = NULL,
                clinic.type = NULL,
                center = TRUE, scale = TRUE,
                scaletype = 'auto',
                epsilon = 0.00001,
                min.variation = 0,
                interactions.reg = TRUE,
                family.glm = gaussian(),
                elasticnet.glm = NULL,
                col.filter.glm = 'cor',
                correlation.glm = 0.7,
                thres.isgl = 0.7,
                gr.method.isgl = 'cor',
                alfa.pls = 0.05,
                p.method.pls = 'jack',
                vip.pls = 0.8,
                method  ='glm'){
  
  if(is.null(omic.type)){
    #Create internally omic.type vector
    omic.type = sapply(data.omics, function(x) isBin(x)) #We assume that all regulators are of the same type in the same omic
    
    cat('Considering that we codify by 1 binary omics and by 0 numeric omics, we consider the following nature of the omics: \n')
    print(omic.type)
    
    cat('Please if it is incorrect stop the generation of the models and introduce omic.type manually \n')
  }
  
  if(!is.null(clinic)){
    if(is.null(clinic.type)){
      clinic.type = sapply(clinic, function(x) isBinclinic(x))
      
      cat('Considering that we codify by 1 binary/categorical features and by 0 numeric features, we consider the following nature of the clinical features:\n')
      print(clinic.type)
      
      cat('Please if it is incorrect stop the generation of the models and introduce clinic.type manually \n')
    }
  }
  
  if(method=='glm'){
    
    return(GetGLM(GeneExpression=GeneExpression,
                  data.omics=data.omics,
                  associations=associations,
                  omic.type = omic.type,
                  edesign = edesign,
                  clinic = clinic,
                  clinic.type = clinic.type,
                  center = center, scale = scale,
                  epsilon = epsilon,
                  family = family.glm,
                  elasticnet = elasticnet.glm,
                  interactions.reg = interactions.reg,
                  min.variation = min.variation,
                  col.filter = col.filter.glm,
                  correlation = correlation.glm,
                  scaletype = scaletype))
    
  }
  
  if(method=='isgl'){
    
    return(GetISGL(GeneExpression=GeneExpression,
                  data.omics=data.omics,
                  associations=associations,
                  omic.type = omic.type,
                  edesign = edesign,
                  clinic = clinic,
                  clinic.type = clinic.type,
                  center = center, scale = scale,
                  interactions.reg = interactions.reg,
                  min.variation = min.variation,
                  gr.method = gr.method.isgl,thres = thres.isgl))
    
  }
  
  else {
    
    return(GetPLS(GeneExpression=GeneExpression,
                  data.omics=data.omics,
                  associations=associations,
                  omic.type = omic.type,
                  edesign = edesign,
                  clinic = clinic,
                  clinic.type = clinic.type,
                  center = center, scale = scale,
                  epsilon = epsilon,
                  alfa = alfa.pls, 
                  interactions.reg = interactions.reg,
                  min.variation = min.variation,
                  scaletype = scaletype,
                  p.method =p.method.pls, 
                  vip = vip.pls,
                  method = method))


  }
  
}
