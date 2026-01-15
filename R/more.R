#########################################################################################
######           MORE:      Multi-Omics REgulation                        ######
#########################################################################################


## By Sonia, Monica & Maider
## 15-October-2023
## Last modified: August 2023


options(stringsAsFactors = FALSE)

# library(igraph)
# library(MASS)
# library(glmnet)
# library(psych)
# library(car)
# library(furrr)
#' @import igraph
#' @import MASS
#' @import glmnet
#' @importFrom psych phi
#' @import car
#' @import furrr
#' @import ggplot2
#' @import ropls
#' @import sglfast
#' @import ltm

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
  
  if(inherits(x,'character')){
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
#' \code{more} fits a MLR regression model (when the selected method is MLR) or a PLS model (when the selected method is PLS1 or PLS2) for all target features in the dataset to identify
#' the potential regulators that show a significant impact on the expression of the target omic features under specific experimental conditions.
#'
#' @param targetData  Matrix or data.frame containing data from a target omic with its features in rows and samples in columns. The row names must be the target omic features IDs; e.g. when gene expression is considered as the targetData, gene IDs should be the row names.
#' @param regulatoryData List where each element corresponds to a different omic data type (miRNAs, transcription factors, methylation, etc.). The names of this list will be the omics, and each element of the list is a  matrix or data frame with omic regulators in rows and samples in columns. Attention! We clarify that the user can not use as regulatory omic the data considered as \code{targetData}.
#' @param associations List where each element corresponds to a different omic data type (miRNAs, transcription factors, methylation, etc.). The names of the elements of the list will be the omics (in the same order as in \code{regulatoryData}). Each element is a data.frame with two columns (optionally three) describing the potential interactions between target omic features and potential regulators for that omic. The first column must contain the regulators, the second the target features IDs, and an additional column can be added to describe the type of interaction (for example, in methylation data, if a CpG site is located in the promoter region of the gene, in the first exon, or any other information). Optionally, the user can set the \code{associations} data.frame of an omic equal to NULL if they want to consider all the regulators of that omic as potential regulators for all the target features. They can even set \code{associations} to NULL if they want to consider all regulators of all omics in \code{regulatoryData} as potential regulators to all target features. Even if it can be done, we do not recommend the user to do it as it can be very time-consuming.
#' @param omicType Vector with as many elements as the number of omics, indicating whether the omic values are numeric (0) or binary (1). When NULL is indicated, MORE will estimate which type of omics are provided and display them on the screen. If a single value is provided, the type for all the omics is set to that value. By default, NULL. If the estimated type of omics are incorrect, the user must halt the process and manually specify the correct omic type.
#' @param condition Data.frame or matrix describing the condition to which samples belong. Rows must be the samples, that is, the columns in the \code{targetData}, and columns must be the conditions to be included in the model, such as disease, treatment, etc.
#' @param clinic Data frame or matrix containing clinical variables values where rows must represent samples and columns variables.
#' @param clinicType Vector with as many elements as the number of clinical variables, indicating whether the variables values are numeric (0) or categorical/binary (1). When NULL is indicated, MORE will estimate which type of variables are provided and display them on the screen. If a single value is provided, the type for all the variables is set to that value. By default, NULL. If the estimated type of variables are incorrect, the user must halt the process and manually specify the clinic type.
#' @param minVariation  Vector with as many elements as the number of omics (names of this vector will be the omics), indicating the minimum change in the standard deviation that a regulator must show across conditions in order not to be considered as having low variation and be removed from the regression models, for numerical regulators. For binary regulators, the minimum change in the proportion a regulator must show across conditions. When a single value is given, the minimum change will be considered the same for all omics and equal to this value. The user has the option to set this argument to NA if they do not want to provide a value but want to filter more than constant regulators across conditions. In this case, the value will be calculated as the 10\% of the maximum observed variability across conditions for continuous regulators and as the 10\% of the maximum observed proportion difference across conditions for binary regulators. Additionally, the user can combine both functionalities; indeed, the user has the option to provide a vector containing the minimum change in the standard deviation for some omics and NA others. By default, this argument is set to 0. If the selected threshold is excessively restrictive, resulting in the removal of all regulators, the model will fail to run and a warning message will be shown.
#' @param percNA Maximum percentage of missing values allowed in regulatoryData to construct the models. Only used in PLS models since MLR models do not allow missing values. By default, 0.2.
#' @param scaleType Type of scaling to be applied when adjusting a model, if scaling is requested. It can be:
#' \itemize{
#' \item auto : Applies the auto scaling; so that scales each variable independently. 
#' \item softBlock : Applies the pareto scaling. Type of block-scaling. \deqn{\frac{X_k}{s_k \sqrt[4]{m_b}} }
#' \item hardBlock : Applies the block scaling. Type of block-scaling. \deqn{ \frac{X_k}{s_k \sqrt{m_b}} }
#' \item none : It does not apply any type of scaling. Not recommended if the user does not apply their own scaling.
#' }
#' considering m_b the number of variables of the block. By default, auto.
#' @param epsilon A threshold for the convergence tolerance in the MLR model. By default, 0.00001.
#' @param interactions If TRUE (default), MORE allows for interactions between each regulator and the conditions, or phenotypes under study. 
#' @param varSel Variable selection method to apply. The options are different for MLR or PLS. Four options:
#' \itemize{
#' \item EN : Applies a Multiple Linear Regression (MLR) model with ElasticNet (EN) regularization.
#' \item ISGL : Applies a Multiple Linear Regression (MLR) model with Iterative Sparse Group Lasso (ISGL) regularization.
#' \item Jack : Applies the Jack-Knife resampling technique to compute the significance of the coefficients in Partial Least Squares (PLS) models.
#' \item Perm : Applies a resampling technique to compute the significance of the coefficients in Partial Least Squares (PLS) models in which the response variable is permuted 100 times to obtain the distribution of the coefficients and compute then their associated p-value.
#' }
#' By default, EN.
#' @param alfaEN ElasticNet mixing parameter (\eqn{\alpha}). By default, NULL. These are the values that can be passed to this argument:
#' \itemize{
#' \item NULL : \eqn{\alpha} parameter will be automatically optimized by cross-validation. For computational efficiency, only values ranging from 0 to 1 in increments of 0.1 will be tested.
#' \item A number between 0 and 1 : ElasticNet is applied with this \eqn{\alpha} being the combination between ridge and lasso penalization.
#' \item A vector with the mixing parameters to try: ElasticNet will be applied for each of the \eqn{\alpha} values provided in the vector, and the one that yields the best cross-validation performance will be selected.
#' \item 0 : A ridge penalty will be applied.
#' \item 1: A lasso penalty will be applied.
#' }
#' The shrinkage parameter (\eqn{\lambda}) will be optimized by cross-validation in all cases. By default, NULL.
#' @param correlation  Correlation threshold (in absolute value) for the Multicollinearity Filter (MF) to decide which regulators are correlated, in which case, a/some representative of the group of correlated regulators is chosen to enter the model. By default, 0.7.
#' @param groupingISGL Option to create the groups for the Iterative Sparse Group Lasso regularization. By default, MF_0.7. There are two options:
#' \itemize{
#' \item MF_X: Takes the same approach as in the multicollinearity filter. Grouping the variables which correlate more than X. E.g: MF_0.7 groups variables with correlation higher than 0.7. By default, MF_0.7. 
#' \item PCA_X: Extracts as many PCA components to explain X percent of the variability of the data and assigns the variables to the component in which they had the highest loadings. E.g: PCA_0.8 extracts the minimum number of components required to explain at least 80\% of the variability and groups variables according to the component they obtained the highest loading.
#' }
#' @param alfa Significance level to consider a regulator as significant in PLS models. By default, 0.05.
#' @param vip The Variable Importance in Projection, VIP, score threshold to apply together with the alfa threshold to take a variable as significant; both requirements should be met to take a variable as significant. By default, 0.8.
#' @param method Regression model to be applied:
#' \itemize{
#' \item MLR : Applies a Multiple Linear Regression (MLR).
#' \item PLS1 : Applies a Partial Least Squares (PLS) model, one for each of the features in the target omic of \code{targetData}.
#' \item PLS2 : Applies a PLS model to all features of the target omic simultaneously, only possible when \code{associations}= NULL.
#' }
#' By default, MLR.
#' @param parallel If FALSE, MORE will be run sequentially. If TRUE, MORE will be run using parallelization with as many cores as the available ones minus one so the system is not overcharged. If the user wants to specify how many cores they want to use, they can also provide the number of cores to use in this parameter. Parallelization is only implemented for MLR with ElasticNet regularization and PLS1 methods. By default, FALSE.
#' @param seed Sets the seed to guaranty reproducibility of the results for the methods that require random sampling of the observations (i.e., all MLR approaches and PLS with Perm variable selection). By default, 123.
#' @return List containing the following elements:
#' \itemize{
#' \item ResultsPerTargetF : List with as many elements as features of the target omic in \code{targetData}. For each feature, it includes information about the feature values, considered variables, estimated coefficients,
#'                    detailed information about all regulators, and regulators identified as relevant (in MLR scenario) or significant (in PLS scenarios).
#' \item GlobalSummary : List with information about the fitted models, including model metrics, information about regulators, features of the target omic without models, regulators, master regulators and hub target features.
#' \item Arguments : List containing all the arguments used to generate the models.                
#' }
#'
#' @examples
#' 
#' data(TestData)
#' 
#' #Omic type
#' omicType = c(1,0,0)
#' names(omicType) = names(TestData$data.omics)
#' SimMLR = more(targetData = TestData$GeneExpressionDE,
#'               associations = TestData$associations, 
#'               regulatoryData = TestData$data.omics,
#'               omicType = omicType,
#'               condition = TestData$edesign,
#'               scaleType = 'auto', varSel = 'EN', 
#'               epsilon = 0.00001, alfaEN = NULL,
#'               interactions = TRUE, minVariation = 0,  
#'               correlation = 0.7, method  ='MLR')
#' 
#' @export

more <-function(targetData,
                regulatoryData,
                associations = NULL,
                omicType = NULL,
                condition = NULL,
                clinic = NULL,
                clinicType = NULL,
                minVariation = 0,
                percNA = 0.2,
                scaleType = 'auto',
                epsilon = 0.00001,
                interactions = TRUE,
                varSel = 'EN',
                alfaEN = NULL,
                correlation = 0.7,
                groupingISGL = 'MF_0.7',
                alfa = 0.05,
                vip = 0.8,
                method  ='MLR',
                parallel = FALSE,
                seed = 123){
  
  #Set the seed for the reproducibility
  set.seed(seed)
  
  if(is.null(omicType)){
    #Create internally omicType vector
    omicType = sapply(regulatoryData, function(x) isBin(x)) #We assume that all regulators are of the same type in the same omic
    
    cat('Considering that we codify by 1 binary omics and by 0 numeric omics, we consider the following nature of the omics: \n')
    print(omicType)
    
    cat('Please if it is incorrect stop the generation of the models and introduce omicType manually \n')
  }
  
  if(!is.null(clinic)){
    if(is.null(clinicType)){
      clinicType = sapply(clinic, function(x) isBinclinic(x))
      
      cat('Considering that we codify by 1 binary/categorical features and by 0 numeric features, we consider the following nature of the clinical features:\n')
      print(clinicType)
      
      cat('Please if it is incorrect stop the generation of the models and introduce clinicType manually \n')
    }
  }
  
  
  
  if(varSel=='EN'){
    
    if(method!='MLR') { stop("ERROR: Selected variable selection technique not available for PLS models")}
    
    return(GetMLR(targetData=targetData,
                  regulatoryData=regulatoryData,
                  associations=associations,
                  omicType = omicType,
                  condition = condition,
                  clinic = clinic,
                  clinicType = clinicType,
                  epsilon = epsilon,
                  family = gaussian(),
                  elasticnet = alfaEN,
                  interactions = interactions,
                  minVariation = minVariation,
                  col.filter = 'cor',
                  correlation = correlation,
                  scaleType = scaleType,
                  parallel = parallel))
    
  }
  
  if(varSel=='ISGL'){
    
    if(method!='MLR') { stop("ERROR: Selected variable selection technique not available for PLS models")}
    
    gr.method = strsplit(groupingISGL, "_")[[1]][1]
    thres = as.numeric(strsplit(groupingISGL, "_")[[1]][2])
    
    return(GetISGL(targetData=targetData,
                  regulatoryData=regulatoryData,
                  associations=associations,
                  omicType = omicType,
                  condition = condition,
                  clinic = clinic,
                  clinicType = clinicType,
                  interactions = interactions,
                  minVariation = minVariation,
                  gr.method = gr.method, thres = thres))
    
  }
  
  else {
    
    if(method == 'MLR') { stop("ERROR: Selected variable selection technique not available for MLR models")}
    
    return(GetPLS(targetData=targetData,
                  regulatoryData=regulatoryData,
                  associations=associations,
                  omicType = omicType,
                  condition = condition,
                  clinic = clinic,
                  clinicType = clinicType,
                  alfa = alfa, 
                  interactions = interactions,
                  minVariation = minVariation,
                  percNA = percNA,
                  scaleType = scaleType,
                  varSel =varSel, 
                  vip = vip,
                  method = method,
                  parallel = parallel))


  }
  
}
