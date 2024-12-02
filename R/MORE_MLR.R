#########################################################################################
######           Functions to integrate omics data using MLRs                      ######
#########################################################################################


## By Sonia, Maider & Monica
## 07-July-2016
## Last modified: October 2023


options(stringsAsFactors = FALSE)

library(ltm)

#' Generalized Linear Models
#'
#'\code{GetMLR} fits a regression model for all the target omic features in the data set to identify
#' the conditions and potential regulators that show a significant effect on
#' the expression of each target omic feature.
#'
#' @param targetData Data frame containing data from a target omic with its features in rows and
#' samples (in different conditions) in columns. Row names must be the target omic features IDs; e.g. when Gene Expression is considered gene IDs.
#' @param regulatoryData List where each element corresponds to a different omic data type to be considered (miRNAs,
#' transcription factors, methylation, etc.). The names of the list will represent the omics, and each element in 
#' the list should be a data matrix with omic regulators in rows and samples in columns.
#' @param associations List where each element corresponds to a different omic data type (miRNAs,
#' transcription factors, methylation, etc.). The names of the list will represent the omics. Each element in 
#' the list should be a data frame with 2 columns (optionally 3), describing the potential interactions between target omic features
#' and regulators for that omic. First column must contain the target omic features (the ones of 
#' \code{\link{targetData}} object), second column must contain the regulators, and an optional third column can
#' be added to describe the type of interaction (e.g., for methylation, if a CpG site is located in
#' the promoter region of the gene, in the first exon, etc.). If the user lacks prior knowledge of the potential regulators, they can set the parameter to NULL. 
#' In this case, all regulators in \code{\link{regulatoryData}} will be treated as potential regulators for all target omic features. In this case, for computational efficiency, it is recommended to use PLS2 \code{\link{method}}.
#' Additionally, if the users have prior knowledge for certain omics and want to set other omics to NULL, they can do so.
#' @param omicType Vector which indicates the type of data of omics introduced in \code{\link{regulatoryData}}. The user should code as 0 numeric omics and as 1 categorical or binary omics. 
#' By default is set to NULL. In this case, the data type will be predicted automatically. However, the user must verify the prediction and manually input the vector if incorrect.
#' @param condition Data frame describing the condition or phenotype to which considered samples belong. Rows must be the samples (columns
#' in \code{\link{targetData}}) and columns must be the conditions or phenotypes to be included in the model (e.g., treatment, etc.).
#' @param clinic Data frame with all clinical variables to consider, with samples in rows and variables in columns.
#' @param clinicType Vector which indicates the type of data of variables introduced in \code{\link{clinic}}. The user should code as 0 numeric variables and as 1 categorical or binary variables. 
#' By default is set to NULL. In this case, the data type will be predicted automatically. However, the user must verify the prediction and manually input the vector if incorrect.
#' @param minVariation  For numerical regulators, it specifies the minimum change required across conditions to retain the regulator in 
#' the regression models. In the case of binary regulators, if the proportion of the most common value is equal to or inferior this value, 
#' the regulator is considered to have low variation and will be excluded from the regression models. The user has the option to set a single 
#' value to apply the same filter to all omics, provide a vector of the same length as omics if they want to specify different levels for each omics, 
#' or use 'NA' when they want to apply a minimum variation filter but are uncertain about the threshold. By default, 0.
#' @param scaleType Type of scaling to be applied. Four options:
#' \itemize{
#' \item auto : Applies the auto scaling. 
#' \item softBlock : Applies the pareto scaling. \deqn{\frac{X_k}{s_k \sqrt[4]{m_b}} }
#' \item hardBlock : Applies the block scaling. \deqn{ \frac{X_k}{s_k \sqrt{m_b}} }
#' \item none : It does not apply any type of scaling. Not recommended if the user does not apply their own scaling.
#' }
#' considering m_b the number of variables of the block. By default, auto.
#' @param epsilon Convergence threshold for coordinate descent algorithm in Elasticnet. Default value, 1e-5.
#' @param interactions If TRUE, the model includes interactions between regulators and condition variables. By default, TRUE.
#' @param varSel Type of variable selection method to apply, different options depending on MLR or PLS usage. Four options:
#' \itemize{
#' \item EN : Applies a Multiple Linear Regression (MLR) with ElasticNet regularization.
#' \item ISGL : Applies a Multiple Linear Regression (MLR) with Iterative Sparse Group Lasso (ISGL) regularization.
#' \item Jack : Applies Jack-Knife resampling technique for the calculation of the significance of the coefficients in Partial Least Squares (PLS) models.
#' \item Perm : Applies a resampling technique for the calculation of the significance of the coefficients in Partial Least Squares (PLS) models in which the response variable is permuted 100 times to obtain the distribution of the coefficients and compute then their associated p-value.
#' }
#' By default, Jack.
#' @param alfaEN ElasticNet mixing parameter. There are three options:
#' \itemize{
#' \item NULL : The parameter is selected from a grid of values ranging from 0 to 1 with 0.1 increments. The chosen value optimizes the mean cross-validated error when optimizing the lambda values.
#' \item A number between 0 and 1 : ElasticNet is applied with this number being the combination between Ridge and Lasso penalization (elasticnet=0 is the ridge penalty, elasticnet=1 is the lasso penalty). 
#' \item A vector with the mixing parameters to try. The one that optimizes the mean cross-validated error when optimizing the lambda values will be used.
#' }
#' By default, NULL.
#' @param correlation  Value to determine the presence of collinearity between two regulators when using the MLR \code{\link{method}}. By default, 0.7.
#' @param parallel If FALSE, MORE will be run sequentially. If TRUE, MORE will be run using parallelization with as many cores as the available ones minus one so the system is not overcharged. If the user wants to specify how many cores they want to use, they can also provide the number of cores to use in this parameter. Parallelization is only implemented for MLR with EN variable selection and PLS methods.
#' @return List containing the following elements:
#' \itemize{
#' \item ResultsPerTargetF : List with as many elements as features of the target omic in \code{\link{targetData}}. For each feature, it includes information about the feature values, considered variables, estimated coefficients,
#'                    detailed information about all regulators, and regulators identified as relevant (in MLR scenario) or significant (in PLS scenarios).
#' \item GlobalSummary : List with information about the fitted models, including model metrics, information about regulators, features of the target omic without models, regulators, master regulators and hub target features.
#' \item Arguments : List containing all the arguments used to generate the models.                
#' }
#' @export

GetMLR = function(targetData,
                  regulatoryData,
                  associations = NULL,
                  omicType = 0,
                  condition = NULL,
                  clinic = NULL,
                  clinicType =NULL,
                  epsilon = 0.00001,
                  family = gaussian(),
                  elasticnet = NULL,
                  interactions = TRUE,
                  minVariation = 0,
                  col.filter = 'cor',
                  correlation = 0.7,
                  scaleType = 'auto',
                  parallel = FALSE){


  # Converting matrix to data.frame
  targetData = as.data.frame(targetData)
  regulatoryData = lapply(regulatoryData, as.data.frame)
  
  ## Omic types
  if (length(omicType) == 1) omicType = rep(omicType, length(regulatoryData))
  names(omicType) = names(regulatoryData)
  
  # Creating vector for minVariation
  if (length(minVariation) == 1) minVariation=rep(minVariation,length(regulatoryData))
  names(minVariation)=names(regulatoryData)
  
  if(!is.null(clinic)){
    
    ## Clinic types
    if (length(clinicType) == 1) {clinicType = rep(clinicType, ncol(clinic)); names(clinicType) = colnames(clinic)}
    
    ##Before introducing variables in regulatoryData convert them to numeric type
    ## TO DO: Careful creates k-1 dummies. Is what we want?
    catvar <- which(clinicType == 1)
    dummy_vars <- model.matrix(~ . , data = as.data.frame(clinic[,catvar ,drop=FALSE]))[,-1,drop=FALSE]
    clinic <-clinic[, -catvar,drop=FALSE]
    clinic <- cbind(clinic, dummy_vars)
    
    regulatoryData = c(list(clinic =  as.data.frame(t(clinic))),regulatoryData)
    
    #Add in associations clinic to consider all the clinical variables in all target features
    if(!is.null(associations)){associations = c(list(clinic = NULL),associations)}
    
    #Add information to omicType and minVariation even if it is not relevant
    omicType = c(0,omicType)
    names(omicType)[1] = 'clinic'
    
    minVariation = c(0,minVariation)
    names(minVariation)[1] = 'clinic'
    om= 2
    
  }else{clinicType=NULL; om =1}
  
  # If associations is NULL create a list of associations NULL
  if (is.null(associations)){
    associations=vector('list',length(regulatoryData))
    names(associations)=names(regulatoryData)
  }

  # Checking that the number of samples per omic is equal to number of samples for targetData and the number of samples for condition
  for (i in 1:length(names(regulatoryData))){
    if(!ncol(regulatoryData[[i]]) == ncol(targetData) ) {
      stop("ERROR: Samples in regulatoryData must be the same as in targetData and in condition")
    }
  }
  if(!is.null(condition)){
    if(!ncol(targetData) == nrow(condition) ) {
      stop("ERROR: Samples in regulatoryData must be the same as in targetData and in condition")
    }
  }

  ## Checking that samples are in the same order in targetData, regulatoryData and condition
  orderproblem<-FALSE
  if(is.null(condition)){
    nameproblem<-!all(sapply(regulatoryData, function(x) length(intersect(colnames(x),colnames(targetData))==ncol(targetData))))
    if(nameproblem){
      cat('Warning. targetData and regulatoryData samples have not same names. We assume that they are ordered.\n')
    }else{
      orderproblem<-!all(sapply(regulatoryData, function(x) identical(colnames(x),colnames(targetData))))
      if(orderproblem){
        regulatoryData<-lapply(regulatoryData, function(x) x[,colnames(targetData)])
      }
    }
    
  } else{
    nameproblem<-!all(c(sapply(regulatoryData, function(x) length(intersect(colnames(x),colnames(targetData)))==ncol(targetData)), length(intersect(rownames(condition),colnames(targetData)))==ncol(targetData)))
    if(nameproblem){
      cat('Warning. targetData, condition and regulatoryData samples have not same names. We assume that they are ordered.\n')
    } else{
      orderproblem<-!all(c(sapply(regulatoryData, function(x) identical(colnames(x),colnames(targetData))), identical(colnames(targetData),rownames(condition))))
      if(orderproblem){
        regulatoryData<-lapply(regulatoryData, function(x) x[,colnames(targetData)])
        condition<-condition[colnames(targetData), , drop=FALSE]
      }
    }
  }

  ## Checking if there are regulators with "_R", "_P" or "_N" or with ":" and checking that there are not replicate identifiers compared to targetData
  message = FALSE
  for (i in 1:length(names(regulatoryData))){

    problemas = c(rownames(regulatoryData[[i]])[grep("_R$", rownames(regulatoryData[[i]]))],
                  rownames(regulatoryData[[i]])[grep("_P$", rownames(regulatoryData[[i]]))],
                  rownames(regulatoryData[[i]])[grep("_N$", rownames(regulatoryData[[i]]))])

    problema = c(grep(":", rownames(regulatoryData[[i]]), value = TRUE))
    rownames(regulatoryData[[i]]) = gsub(':', '-', rownames(regulatoryData[[i]]))
    rownames(regulatoryData[[i]]) = gsub('_R$', '-R', rownames(regulatoryData[[i]]))
    rownames(regulatoryData[[i]]) = gsub('_P$', '-P', rownames(regulatoryData[[i]]))
    rownames(regulatoryData[[i]]) = gsub('_N$', '-N', rownames(regulatoryData[[i]]))

    #Change the name in the association matrix only if associations is not NULL
    if(!is.null(associations[[i]])){
      associations[[i]][[2]]=gsub(':', '-', associations[[i]][[2]])
      associations[[i]][[2]] = gsub('_R$', '-R', associations[[i]][[2]])
      associations[[i]][[2]] = gsub('_P$', '-P', associations[[i]][[2]])
      associations[[i]][[2]] = gsub('_N$', '-N', associations[[i]][[2]])
    }

    if(length(problemas) > 0) {
      cat("In",names(regulatoryData)[i], ',', problemas ,"regulators have names that may cause conflict with the algorithm by ending in _R, _P or _N", "\n")
      cat("Endings changed with -R, -P or -N, respectively", "\n")
    }

    if(length(problema) > 0) {
      cat("Some regulators in the omic", names(regulatoryData)[i],  "have names with \":\" that could cause conflict, replaced with \"-\" ", "\n")
      cat("Changed identifiers: ", problema, "\n")
    }
    
    #Checking that there are not replicate identifiers compared to targetData
    
    repeated = intersect(rownames(targetData), rownames(regulatoryData[[i]]))
    
    if(length(repeated) > 0) {
      cat(names(regulatoryData)[i], "omic and target omic have shared identifiers in regulators:", repeated, "\n")
      #Change the name in the association matrix only if is not NULL
      if(!is.null(associations[[i]])){
        associations[[i]][[2]][associations[[i]][[2]]%in%repeated] = paste(names(regulatoryData)[i],'-', associations[[i]][[2]][associations[[i]][[2]]%in%repeated], sep='')
      }
      #Change the name in regulatoryData
      rownames(regulatoryData[[i]])[rownames(regulatoryData[[i]])%in%repeated] = paste(names(regulatoryData)[i],'-',  rownames(regulatoryData[[i]])[rownames(regulatoryData[[i]])%in%repeated], sep='')
    }
  }

  ##Checking that there are no replicates in the identifiers and changing identifiers in case of need
  if(length(names(regulatoryData))>1){
    for (i in 1:(length(names(regulatoryData))-1)){
      for(j in (i+1):(length(names(regulatoryData)))){
        repeated = intersect(rownames(regulatoryData[[i]]), rownames(regulatoryData[[j]]))
        
        if(length(repeated) > 0) {
          cat(names(regulatoryData)[i], "and", names(regulatoryData)[j], "omics have shared identifiers in regulators:", repeated, "\n")
          #Change the name in the association matrix only if is not NULL
          if(!is.null(associations[[i]])){
            associations[[i]][[2]][associations[[i]][[2]]%in%repeated] =  paste(names(regulatoryData)[i],'-', associations[[i]][[2]][associations[[i]][[2]]%in%repeated], sep='')
          }
          if(!is.null(associations[[j]])){
            associations[[j]][[2]][associations[[j]][[2]]%in%repeated] =  paste(names(regulatoryData)[j],'-', associations[[j]][[2]][associations[[j]][[2]]%in%repeated], sep='')
          }
          #Change the name in regulatoryData
          rownames(regulatoryData[[i]])[rownames(regulatoryData[[i]])%in%repeated] =  paste(names(regulatoryData)[i],'-', rownames(regulatoryData[[i]])[rownames(regulatoryData[[i]])%in%repeated],sep='')
          rownames(regulatoryData[[j]])[rownames(regulatoryData[[j]])%in%repeated] =  paste(names(regulatoryData)[j],'-', rownames(regulatoryData[[j]])[rownames(regulatoryData[[j]])%in%repeated],sep = '')
          
        }
      }
    }
  }

  

  # Preparing family for ElasticNet variable selection
  family2 = family$family
  family2 = strsplit(family2, "(", fixed = TRUE)[[1]][1]

  if (family2 %in% c("poisson", "quasipoisson", "Negative Binomial")) {
    family2 = "poisson"
  } else if (family2 %in% c("gaussian", "binomial")) {
    family2 = family2
  } else {
    family2 = NULL
    message(sprintf("Warning message:"))
    message(sprintf("Elasticnet variable selection cannot be applied for family %s", family2))
  }

  #Checking there are not -Inf/Inf values and eliminate target features/regulator that contain them
  infproblemtargetF<-is.infinite(rowSums(targetData))
  infproblemreg<-lapply(regulatoryData[om:length(regulatoryData)], function(x) is.infinite(rowSums(x)))
  if(any(infproblemtargetF)){
    targetFsInf<-rownames(targetData)[infproblemtargetF]
    targetData<-targetData[!infproblemtargetF,]
  }else{targetFsInf <-NULL}
  for (i in 1:(length(names(regulatoryData))-(om-1))){
    if(any(infproblemreg[[i]])){
      cat(rownames(regulatoryData[[i + (om-1)]])[infproblemreg[[i]]], 'regulators of the omic', names(regulatoryData)[i +(om-1)] ,'have been deleted due to -Inf/Inf values. \n')
      regulatoryData[[i + (om-1)]]<-regulatoryData[[i + (om-1)]][!infproblemreg[[i]],]
    }
  }

  ## Removing target features with NAs and keeping track
  min.obs = ncol(targetData)
  targetFsNotNA = apply(targetData, 1, function (x) sum(!is.na(x)))
  targetFsNotNA = names(which(targetFsNotNA >= min.obs))
  targetFsNA = setdiff(rownames(targetData), targetFsNotNA)
  targetData = targetData[targetFsNotNA,]

  ## Removing target features with no regulators only if associations does not have an associations = NULL in any omic
  targetFsNOreg = NULL
  targetFsNOreg = lapply(associations, function(x) if(!is.null(x)) {setdiff( rownames(targetData),x[,1])})
  targetFsNOreg = Reduce(intersect, targetFsNOreg)
  targetData = targetData[!(rownames(targetData) %in% targetFsNOreg),]
  if (length(targetFsNOreg) > 0){
    cat(length(targetFsNOreg), "target features had no initial regulators. Models will be computed for", nrow(targetData), 'target features.\n')
  }

  ## Removing constant target features
  constantTargetFs = apply(targetData, 1, sd, na.rm = TRUE)
  notConstant = names(constantTargetFs)[constantTargetFs > 0]
  constantTargetFs = names(constantTargetFs)[constantTargetFs == 0]
  targetData = targetData[notConstant,]

  AlltargetFs=rownames(targetData)
  ntargetFs = length(AlltargetFs)

  # Experimental groups
  if (is.null(condition)) {
    cat("No experimental covariates were provided.\n")
    Group = 1:ncol(targetData)
    names(Group) = colnames(targetData)
    des.mat = NULL
  } else {
    Group = apply(condition, 1, paste, collapse = "_")
    des.mat = model.matrix(~Group)[, -1, drop = FALSE]
    rownames(des.mat) = colnames(targetData)
    #Change the name to avoid conflicts with RegulationPerCondition
    colnames(des.mat) = sub('Group','Group_',colnames(des.mat))
  }

  ## Remove regulators with NA
  cat("Removing regulators with missing values...\n")

  myregNA = lapply(regulatoryData, rownames)
  regulatoryData = lapply(regulatoryData, na.omit)
  myregNA = lapply(1:length(regulatoryData), function (i) setdiff(myregNA[[i]], rownames(regulatoryData[[i]])))
  names(myregNA)=names(regulatoryData)

  cat("Number of regulators with missing values:\n")
  print(sapply(myregNA, length))
  cat("\n")


  ## Remove regulators with Low Variability
  cat("Removing regulators with low variation...\n")

  tmp = LowVariationRegu(minVariation, regulatoryData, Group, associations, AlltargetFs, omicType, clinicType)
  regulatoryData = tmp[["data.omics"]]
  associations = tmp[["associations"]]
  myregLV = tmp[["myregLV"]]
  rm("tmp"); gc()
  
  if(all(sapply(regulatoryData, function(x)nrow(x)==0))) stop("ERROR: No regulators left after LowVariation filter. Consider being less restrictive.")

  ### Results objects

  ## Global summary for all target features
  GlobalSummary = vector("list", length = 6)
  names(GlobalSummary) = c("GoodnessOfFit", "ReguPerTargetF", "TargetFNOmodel", "TargetFNOregu", "GlobalRegulators", "HubTargetF")

  GlobalSummary$TargetFNOmodel = NULL
  if (length(targetFsNA) > 0) {
    GlobalSummary$TargetFNOmodel = data.frame("targetF" = targetFsNA,
                                            "problem" = rep("Too many missing values", length(targetFsNA)))
  }
  if (length(constantTargetFs) > 0) {
    GlobalSummary$TargetFNOmodel = rbind(GlobalSummary$TargetFNOmodel,
                                       data.frame("targetF" = constantTargetFs,
                                                  "problem" = rep("Response values are constant", length(constantTargetFs))))
  }
  if (length(targetFsInf) > 0){
    GlobalSummary$TargetFNOmodel = rbind(GlobalSummary$TargetFNOmodel,
                                       data.frame("targetF" = targetFsInf,
                                                  "problem" = rep("-Inf/Inf values", length(targetFsInf))))

  }

  GlobalSummary$TargetFNOregu = NULL
  if (length(targetFsNOreg) > 0){
    GlobalSummary$TargetFNOregu = data.frame("targetF" = targetFsNOreg, "problem" = rep("Target feature had no initial regulators", length(targetFsNOreg)))
  }

  GlobalSummary$GoodnessOfFit = matrix(NA, ncol = 4, nrow = ntargetFs)
  rownames(GlobalSummary$GoodnessOfFit) = AlltargetFs
  colnames(GlobalSummary$GoodnessOfFit) = c("Rsquared", "RMSE","NRMSE", "relReg")

  GlobalSummary$ReguPerTargetF = matrix(0, ncol = 3*length(regulatoryData), nrow = ntargetFs)
  rownames(GlobalSummary$ReguPerTargetF) = AlltargetFs
  colnames(GlobalSummary$ReguPerTargetF) = c(paste(names(regulatoryData), "Ini", sep = "-"),
                                          paste(names(regulatoryData), "Mod", sep = "-"),
                                          paste(names(regulatoryData), "Rel", sep = "-"))

  ## Specific results for each target feature
  ResultsPerTargetF=vector("list", length=length(AlltargetFs))
  names(ResultsPerTargetF) = AlltargetFs
  
  ## Specify the scaling type
  scale <- ifelse(scaleType == 'none', FALSE, TRUE)
  center <- ifelse(scaleType == 'none', FALSE, TRUE)

  ### Computing model for each target feature
  cat("Checking multicollinearity, selecting predictors and fitting model for ...\n")
  
  options(future.globals.maxSize = 4000*1024^2)
  if(!isFALSE(parallel)){
    if(is.numeric(parallel)){
      if(.Platform$OS.type == "unix") {
        future::plan("multicore", workers = parallel)
      }else{
        future::plan("multisession", workers = parallel)
      }
    } else{
      nc = future::availableCores() - 1
      if(.Platform$OS.type == "unix") {
        future::plan("multicore", workers = nc)
      }else{
        future::plan("multisession", workers = nc)
      }
      parallel = nc
    }
    ResultsPerTargetF <- furrr::future_map(1:ntargetFs,
                                        ~ResultsPerTargetF.i.mlr(AlltargetFs[.],GlobalSummary,regulatoryData,associations,targetData,omicType,
                                                              condition, des.mat,myregLV,myregNA,col.filter,correlation,epsilon,
                                                              scale,center,scaleType,interactions,elasticnet,family2),
                                        .progress = TRUE,.options = furrr_options(seed = TRUE))
  } else{
    ResultsPerTargetF <- purrr::map(1:ntargetFs,
                                        ~ResultsPerTargetF.i.mlr(AlltargetFs[.],GlobalSummary,regulatoryData,associations,targetData,omicType,
                                                              condition, des.mat,myregLV,myregNA,col.filter,correlation,epsilon,
                                                              scale,center,scaleType,interactions,elasticnet,family2),
                                        .progress = TRUE,.options = furrr_options(seed = TRUE))
  }
  
  
  future::plan("sequential")
  names(ResultsPerTargetF)<-AlltargetFs
  
  
  globalValues<-c('GoodnessOfFit','ReguPerTargetF','TargetFNOmodel', 'TargetFNOregu')
  
  for(i in 1:length(ResultsPerTargetF)){
    
    for (gValue in globalValues){
      if (exists(gValue, ResultsPerTargetF[[i]]) && ! is.null(ResultsPerTargetF[[i]][[gValue]])){
        if(! is.data.frame(ResultsPerTargetF[[i]][[gValue]])) {
          GlobalSummary[[gValue]][AlltargetFs[i],]<-ResultsPerTargetF[[i]][[gValue]]
        } else{
          GlobalSummary[[gValue]]<-rbind(GlobalSummary[[gValue]], ResultsPerTargetF[[i]][[gValue]])
        }
      }
    }
    ResultsPerTargetF[[i]] <- ResultsPerTargetF[[i]][-grep(paste0(globalValues, collapse = "|"), names(ResultsPerTargetF[[i]]))]
  }
    
  
  # Remove from GoodnessOfFit target features with no significant regulators
  
  targetFsNosig = names(which(GlobalSummary$GoodnessOfFit[,4]==0))
  targetFssig = setdiff(rownames(GlobalSummary$GoodnessOfFit), targetFsNosig)
  GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[targetFssig,,drop=FALSE]
  
  targetFsNoreg = rownames(GlobalSummary$GoodnessOfFit)[is.na(rowSums(GlobalSummary$GoodnessOfFit))]
  targetFsreg = setdiff(rownames(GlobalSummary$GoodnessOfFit), targetFsNoreg)
  GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[targetFsreg,]
  
  #Calculate GlobalRegulators
  m_rel_reg<-lapply(ResultsPerTargetF, function(x) x$relevantRegulators)
  m_rel_reg <- unlist(m_rel_reg)
  mrel_vector <- table(m_rel_reg)
  #Calculate third quantile
  q3<-quantile(mrel_vector,0.75)
  if(length(mrel_vector[mrel_vector>q3])<10){
    GlobalSummary$GlobalRegulators = intersect(names(mrel_vector[rev(tail(order(mrel_vector),10))]), names(mrel_vector[mrel_vector>10]) )
  } else{
    GlobalSummary$GlobalRegulators = intersect(names(mrel_vector[mrel_vector>q3]), names(mrel_vector[mrel_vector>10]) ) 
  }
  
  #Calculate HubTargetF
  relevant_regulators<-GlobalSummary$ReguPerTargetF[,c(grep('-Rel$',colnames(GlobalSummary$ReguPerTargetF))),drop=FALSE]
  s_rel_reg<-apply(relevant_regulators, 1, sum)
  #Calculate third quantile
  q3<-quantile(s_rel_reg,0.75)
  if(length(s_rel_reg[s_rel_reg>q3])<10){
    GlobalSummary$HubTargetF = intersect(names(s_rel_reg[rev(tail(order(s_rel_reg),10))]), names(s_rel_reg[s_rel_reg>10]) )
  } else{
    GlobalSummary$HubTargetF = intersect(names(s_rel_reg[s_rel_reg>q3]), names(s_rel_reg[s_rel_reg>10]))
  }
  
  if(all(sapply(associations,is.null))) {associations = NULL}
  
  myarguments = list(condition = condition, finaldesign = des.mat, groups = Group, 
                     elasticnet = elasticnet, minVariation = minVariation, 
                     correlation = correlation, epsilon = epsilon, 
                     associations = associations, targetData = targetData, 
                     regulatoryData = regulatoryData, omicType = omicType,
                     clinic = clinic, clinicType=clinicType, method ='MLR',
                     parallel = parallel)
  
  result <- list("ResultsPerTargetF" = ResultsPerTargetF, "GlobalSummary" = GlobalSummary, "arguments" = myarguments) 
  class(result) <- "MORE"
  return(result)


}

# MORE main function --------------------------

ResultsPerTargetF.i.mlr<-function(targetF,GlobalSummary,regulatoryData,associations,targetData,omicType,
                           condition, des.mat,myregLV,myregNA,col.filter,correlation,epsilon,
                           scale,center,scaleType,interactions,elasticnet,family2){
  
  
  ResultsPerTargetF.i = vector("list", length = 9)
  names(ResultsPerTargetF.i) = c("Y", "X", "coefficients", "allRegulators", "relevantRegulators", "GoodnessOfFit", "ReguPerTargetF","TargetFNOmodel","TargetFNOregu")
  
  #Initialize global summary values
  
  ResultsPerTargetF.i$ReguPerTargetF <- GlobalSummary$ReguPerTargetF[targetF, ,drop=FALSE]
  ResultsPerTargetF.i$TargetFNOmodel <- GlobalSummary$TargetFNOmodel
  
  RetRegul = GetAllReg(targetF=targetF, associations=associations, data.omics = regulatoryData)
  RetRegul.targetF = RetRegul$Results  ## RetRegul$TableTargetF: nr reg per omic
  ## Some of these reg will be removed, because they are not in regulatoryData
  
  # RetRegul.targetF--> target feature/regulator/omic/area
  RetRegul.targetF=RetRegul.targetF[RetRegul.targetF[,"regulator"]!= "No-regulator", ,drop=FALSE] ## Remove rows with no-regulators
  
  ### NO INITIAL REGULATORS
  if(length(RetRegul.targetF)==0){ 
    
    if (is.null(condition)) {
      ResultsPerTargetF.i$X = NULL
      ResultsPerTargetF.i$relevantRegulators = NULL
      ResultsPerTargetF.i$allRegulators = NULL
      isModel = NULL
      
      ResultsPerTargetF.i$TargetFNOregu = rbind(ResultsPerTargetF.i$TargetFNOregu,
                                                data.frame("targetF" = targetF,
                                                           "problem" = 'Target feature had no initial regulators'))
      
    } else {
      des.mat2 = cbind(t(targetData[targetF,]), des.mat)
      colnames(des.mat2)[1] = "response"
      des.mat2 = na.omit(des.mat2)
      
      # Removing predictors with constant values
      sdNo0 = apply(des.mat2, 2, sd)
      sdNo0 = names(sdNo0)[sdNo0 > 0]
      des.mat2 = des.mat2[,sdNo0]
      
      isModel =NULL
      ResultsPerTargetF.i$X = des.mat2[,-1, drop = FALSE]
      ResultsPerTargetF.i$relevantRegulators = NULL
      ResultsPerTargetF.i$allRegulators = NULL
      
      ResultsPerTargetF.i$TargetFNOregu = rbind(ResultsPerTargetF.i$TargetFNOregu,
                                                data.frame("targetF" = targetF,
                                                           "problem" = 'Target feature had no initial regulators'))
    }
    
    
    # GlobalSummary$ReguPerTargetF  # this is initially set to 0 so no need to modify it
    
    
    ### WITH INITIAL REGULATORS
  }
  else { ## There are regulators for this targetF at the beginning
    
    ResultsPerTargetF.i$allRegulators = data.frame(RetRegul.targetF, rep("Model",nrow(RetRegul.targetF)), stringsAsFactors = FALSE)
    colnames(ResultsPerTargetF.i$allRegulators) = c("targetF","regulator","omic","area","filter")
    
    ResultsPerTargetF.i$ReguPerTargetF[1, grep("-Ini", colnames(ResultsPerTargetF.i$ReguPerTargetF))] = as.numeric(RetRegul$TableTargetF[-1])
    # the rest of columns remain 0
    
    ## Identify which regulators where removed because of missing values or low variation
    res = RemovedRegulators(RetRegul.targetF = ResultsPerTargetF.i$allRegulators,
                            myregLV=myregLV, myregNA=myregNA, data.omics=regulatoryData)
    
    if(length(res$RegulatorMatrix)==0){ ## No regulators left after the filtering to compute the model
      
      if (is.null(condition)) {
        ResultsPerTargetF.i$X = NULL
        ResultsPerTargetF.i$relevantRegulators = NULL
        ResultsPerTargetF.i$allRegulators = res$SummaryPerTargetF
        isModel = NULL
        
      } else {
        des.mat2 = cbind(t(targetData[targetF,]), des.mat)
        colnames(des.mat2)[1] = "response"
        des.mat2 = na.omit(des.mat2)
        
        # Removing predictors with constant values
        sdNo0 = apply(des.mat2, 2, sd)
        sdNo0 = names(sdNo0)[sdNo0 > 0]
        des.mat2 = des.mat2[,sdNo0]
        
        isModel = NULL
        
        ResultsPerTargetF.i$TargetFNOmodel = rbind(ResultsPerTargetF.i$TargetFNOmodel,
                                              data.frame("targetF" = targetF,
                                                         "problem" = 'No regulators left after NA/LowVar filtering'))
        
        ResultsPerTargetF.i$X = des.mat2[,-1, drop = FALSE]
        ResultsPerTargetF.i$relevantRegulators = NULL
        ResultsPerTargetF.i$allRegulators = res$SummaryPerTargetF
        
      }
      
    }
    else {  ## Regulators for the model!!
      
      ## Apply multicollinearity filter only if there is more than one regulator for a target feature
      
      if (ncol(res$RegulatorMatrix)>1){
        if(col.filter=='cor'){
          res = CollinearityFilter1(data = res$RegulatorMatrix, reg.table = res$SummaryPerTargetF,
                                    correlation = correlation, omic.type = omicType, scale = scale, center = center)
        }
        if(col.filter=='pcor'){
          res = CollinearityFilter2(data = res$RegulatorMatrix, reg.table = res$SummaryPerTargetF,
                                    correlation = correlation, omic.type = omicType, epsilon = epsilon , scale = scale, center = center)
        }
        
      }
      
      if(is.null(res)){
        des.mat2 = cbind(t(targetData[targetF,]), des.mat)
        colnames(des.mat2)[1] = "response"
        
        colnames(des.mat2) = gsub("\`", "", colnames(des.mat2))
        
        ResultsPerTargetF.i$TargetFNOmodel = rbind(ResultsPerTargetF.i$TargetFNOmodel,
                                              data.frame("targetF" = targetF,
                                                         "problem" = 'Problem with Partial Correlation calculation'))
        isModel =NULL
      } else{
        ResultsPerTargetF.i$allRegulators = res$SummaryPerTargetF
        
        ## Scaling predictors for ElasticNet only in case they were not already scaled
        des.mat2EN = RegulatorsInteractions(interactions, reguValues = res$RegulatorMatrix, reguInfo = res$SummaryPerTargetF,
                                            des.mat, method = 'mlr')
        
        # Removing observations with missing values
        des.mat2EN = lapply(des.mat2EN, na.omit)
        
        #Scale the variables, indispensable for elasticnet application
        des.mat2EN = lapply(des.mat2EN, function(x)scale(x,scale=scale,center=center))
        
        ##Scale if needed for soft or hard block scaling
        res$RegulatorMatrix = Scaling.type(des.mat2EN,scaleType)
        if(is.null(des.mat)){
          des.mat2EN = data.frame(t(targetData[targetF,]),res$RegulatorMatrix, check.names = FALSE)
        } else{
          des.mat2EN = data.frame(t(targetData[targetF,]),scale(des.mat,scale=scale,center=center), res$RegulatorMatrix, check.names = FALSE)
        }
        colnames(des.mat2EN)[1] = "response"
        rm(res); gc()
        
        ###  Variable selection --> Elasticnet
        tmp = ElasticNet(family2, des.mat2EN, epsilon, elasticnet)
        regulatorcoef = tmp[['coefficients']]
        isModel = tmp[['isModel']]
        m = tmp[['m']]
        des.mat2 = as.data.frame(des.mat2EN[,colnames(tmp[["des.mat2"]]),drop = FALSE])
        ResultsPerTargetF.i$X = des.mat2EN[,-1, drop = FALSE]
        rm(des.mat2EN); gc()
        
      }
      
      
      if (ncol(des.mat2) == 1 || is.null(isModel)) {
        
        ## Extracting significant regulators
        ResultsPerTargetF.i$relevantRegulators = NULL
        ResultsPerTargetF.i$allRegulators = data.frame(ResultsPerTargetF.i$allRegulators, "Rel" = 0, stringsAsFactors = FALSE)
        
        ## Counting original regulators in the model per omic
        contando = ResultsPerTargetF.i$allRegulators[which(ResultsPerTargetF.i$allRegulators[,"filter"] == "Model"),]
        contando = table(contando[,"omic"])
        contando = as.numeric(contando[names(regulatoryData)])
        contando[is.na(contando)] = 0
        ResultsPerTargetF.i$ReguPerTargetF[targetF, grep("-Mod", colnames(ResultsPerTargetF.i$ReguPerTargetF))] = contando
        ResultsPerTargetF.i$TargetFNOmodel = rbind(ResultsPerTargetF.i$TargetFNOmodel,
                                                   data.frame("targetF" = targetF,
                                                              "problem" = 'No relevant regulators after variable selection'))
      } else{
        
        isModel = TRUE
        mycoef = colnames(des.mat2[,-1,drop = FALSE])
        myvariables = unlist(strsplit(mycoef, ":", fixed = TRUE))
        mycondi = intersect(myvariables, colnames(des.mat))
        myvariables = intersect(myvariables, rownames(ResultsPerTargetF.i$allRegulators))
        
        ResultsPerTargetF.i$allRegulators = data.frame(ResultsPerTargetF.i$allRegulators, "Rel" = 0, stringsAsFactors = FALSE)
        ResultsPerTargetF.i$allRegulators[myvariables, "Rel"] = 1
        ResultsPerTargetF.i$coefficients = regulatorcoef
        #ResultsPerTargetF[[i]]$coefficients = regulatorcoef[myvariables,, drop =FALSE]
        colnames(ResultsPerTargetF.i$coefficients) = c('coefficient')
        
        ## A las variables significativas le quito "_R", solo quedara omica_mc"num". Luego creo un objeto que contenga a mi tabla de "allRegulators"
        ## para poder modificar los nombres de la misma forma.
        myvariables = sub("_R", "", myvariables)
        ResultsPerTargetF.i$relevantRegulators = myvariables # significant regulators including "new" correlated regulators without _R
        
        contando = ResultsPerTargetF.i$allRegulators[which(ResultsPerTargetF.i$allRegulators[,"filter"] == "Model"),]
        contando = table(contando[,"omic"])
        contando = as.numeric(contando[names(regulatoryData)])
        contando[is.na(contando)] = 0
        ResultsPerTargetF.i$ReguPerTargetF[targetF, grep("-Mod", colnames(ResultsPerTargetF.i$ReguPerTargetF))] = contando
        
        mytable = ResultsPerTargetF.i$allRegulators
        mytable[,"filter"] = sub("_P", "", mytable[,"filter"])
        mytable[,"filter"] = sub("_N", "", mytable[,"filter"])
        mytable[,"filter"] = sub("_R", "", mytable[,"filter"])
        
        collin.regulators = intersect(myvariables, mytable[,"filter"])
        
        if (length(collin.regulators) > 0) {  # there were correlated regulators
          original.regulators = mytable[mytable[,"filter"] %in% collin.regulators, "regulator"]
          
          ResultsPerTargetF.i$allRegulators[original.regulators, "Rel"] = 1
          ResultsPerTargetF.i$relevantRegulators = c(ResultsPerTargetF.i$relevantRegulators, as.character(original.regulators))
          ResultsPerTargetF.i$relevantRegulators = setdiff(ResultsPerTargetF.i$relevantRegulators, collin.regulators)
          
          ## Counting original regulators in the model per omic
          contando = ResultsPerTargetF.i$allRegulators
          quitar = which(contando[,"filter"] == "MissingValue")
          if (length(quitar) > 0) contando = contando[-quitar,]
          quitar = which(contando[,"filter"] == "LowVariation")
          if (length(quitar) > 0) contando = contando[-quitar,]
          contando = contando[-grep("_R", rownames(contando)),]
        } else {
          contando = ResultsPerTargetF.i$allRegulators[which(ResultsPerTargetF.i$allRegulators[,"filter"] == "Model"),]
        }
        contando = table(contando[,"omic"])
        contando = as.numeric(contando[names(regulatoryData)])
        contando[is.na(contando)] = 0
        ResultsPerTargetF.i$ReguPerTargetF[targetF, grep("-Mod", colnames(ResultsPerTargetF.i$ReguPerTargetF))] = contando
        
        ## Counting significant regulators per omic
        if (length(ResultsPerTargetF.i$relevantRegulators) > 0) {
          contando = ResultsPerTargetF.i$allRegulators[ResultsPerTargetF.i$relevantRegulators,, drop=FALSE]
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(regulatoryData)])
          contando[is.na(contando)] = 0
          ResultsPerTargetF.i$ReguPerTargetF[targetF, grep("-Rel", colnames(ResultsPerTargetF.i$ReguPerTargetF))] = contando
        } else{
          ResultsPerTargetF.i$TargetFNOmodel = rbind(ResultsPerTargetF.i$TargetFNOmodel,
                                                     data.frame("targetF" = targetF,
                                                                "problem" = 'No relevant regulators after variable selection'))
        }
        
      }
    }
    
  } ## Close "else" --> None regulators from begining
  
  if (is.null(isModel)) {
    
    ResultsPerTargetF.i$Y = targetData[targetF,]
    ResultsPerTargetF.i$coefficients = NULL
    
  } else {
    ResultsPerTargetF.i$Y = data.frame("y" = des.mat2[,1], "fitted.y" = tmp[['fitted.values']],
                                    "residuals" = des.mat2[,1] - tmp[['fitted.values']], check.names = FALSE)
    colnames(ResultsPerTargetF.i$Y) <- c("y", "fitted.y", "residuals")
    ResultsPerTargetF.i$GoodnessOfFit = c(m$R.squared, m$RMSE, m$NRMSE,length(ResultsPerTargetF.i$relevantRegulators))
    
    
  }  
  
  return(ResultsPerTargetF.i)
}



# Multi-collinearity filter ------------------------------------------------------


## Multicollinearity filter taking into account correlation of different omics. Method 'COR'

correlations<- function(v, data, reg.table, omic.type){
  omic1 = omic.type[[reg.table[v[1], 'omic']]]
  omic2 = omic.type[[reg.table[v[2], 'omic']]]
  
  if(omic1 == 0 & omic2 == 0){
    correlation = cor(data[, v[1]], data[, v[2]])
  } else if(omic1 == 0 & omic2 == 1){
    correlation = ltm::biserial.cor(data[, v[1]], data[, v[2]])
  } else if(omic1 == 1 & omic2 == 0){
    correlation = ltm::biserial.cor(data[, v[2]], data[, v[1]])
  } else{
    contingency.table = table(data[,v[1]], data[,v[2]])
    correlation = psych::phi(contingency.table)
  }
  
  return(correlation)
}

CollinearityFilter1 = function(data, reg.table, correlation = 0.8, omic.type,scale,center) {
  
  ## data = Regulator data matrix for all omics where missing values and regulators with low variation have been filtered out
  #         (regulators must be in columns)
  ## reg.table = Table with "targetF", "regulator", "omic", "area", filter" where omics with no regulators have been removed
  row.names(reg.table) = reg.table[,"regulator"]
  #Scale the data only for correlation calculation
  data2 = scale(data,scale,center)
  myreg = as.character(reg.table[which(reg.table[,"filter"] == "Model"),"regulator"])
  mycorrelations = data.frame(t(combn(myreg,2)),combn(myreg, 2, function(x) correlations(x, data2, reg.table, omic.type)))
  
  ## Compute the correlation between all regulators (even if they are of different omics)
  mycor = mycorrelations[abs(mycorrelations[,3]) >= correlation,]
  
  if (nrow(mycor) == 1) {  ### only 2 regulators are correlated in this omic
    
    correlacionados = unlist(mycor[,1:2])
    regulators = colnames(data)
    keep = sample(correlacionados, 1) # Regulador al azar de la pareja
    
    ## Lo siguiente elimina el no representante de la matriz de reguladores. Al regulador escogido como representante,
    ## le cambia el nombre por "mc_1_R" para que despues pase la seleccion de variables y asi, en reg.table se conserva
    ## la info de que fue escogido como representante.
    
    remove = setdiff(correlacionados, keep)
    regulators = setdiff(regulators, remove)
    data = as.matrix(data[ ,regulators])
    colnames(data) = regulators
    index.reg = which(colnames(data) == as.character(keep))
    colnames(data)[index.reg] = paste(reg.table[keep[[1]], 'omic'], paste("mc", 1, sep = ""), "R", sep = "_")
    
    # Cambio en reg.table. Asignacion de los nombres segun sea representante,
    # correlacion positiva o negativa. Creacion de una nueva fila con el representante
    # para la seleccion de variables y asi, no perder la info del representante.
    
    reg.table = rbind(reg.table, reg.table[keep,])
    reg.table[nrow(reg.table), "regulator"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", 1, sep = ""), "R", sep = "_")
    reg.table[keep, "filter"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", 1, sep = ""), "R", sep = "_")
    rownames(reg.table) = reg.table[ ,"regulator"]
    
    if(mycor[,3] > 0){
      reg.table[remove, "filter"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", 1, sep = ""), "P", sep = "_")
    } else{
      reg.table[remove, "filter"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", 1, sep = ""), "N", sep = "_")
    }
  }
  
  if (nrow(mycor) >= 2) {   ### more than 2 regulators might be correlated in this omic
    
    mygraph = igraph::graph_from_data_frame(mycor, directed=F)
    mycomponents = igraph::components(mygraph)
    mygraph$community<-mycomponents$membership ##save membership information
    
    for (i in 1:mycomponents$no) {
      
      #create the subgraphs of the clusters
      mysubgraph = igraph::subgraph(mygraph,as.numeric(igraph::V(mygraph)[which(mygraph$community==i)]))
      
      nedges = igraph::ecount(mysubgraph)
      
      ## see if it is a fully connected graph
      if (nedges == ((mycomponents$csize[i]*(mycomponents$csize[i]-1))/2)){
        correlacionados = names(mycomponents$membership[mycomponents$membership == i])
        regulators = colnames(data)
        
        ## Escoge un regulador al azar como representante de cada componente conexa. Para cada componente conexa elimina aquellos reguladores
        ## que no han sido escogidos como representante.
        keep = sample(correlacionados, 1)  # mantiene uno al azar
        reg.remove = setdiff(correlacionados, keep) # correlated regulator to remove
        regulators = setdiff(regulators, reg.remove)  # all regulators to keep
        data = as.matrix(data[ ,regulators])
        colnames(data) = regulators
        index.reg = which(colnames(data) == as.character(keep))
        colnames(data)[index.reg] = paste(reg.table[keep[[1]], 'omic'], paste("mc", i, sep = ""), "R", sep = "_")
        
        # Asignacion de nombre al representante y nueva fila para el filtro de
        # seleccion de variables (asi no se pierde la info del representante).
        reg.table = rbind(reg.table, reg.table[keep,])
        reg.table[nrow(reg.table), "regulator"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", i, sep = ""), "R", sep = "_")
        reg.table[keep, "filter"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", i, sep = ""), "R", sep = "_")
        rownames(reg.table) = reg.table[ ,"regulator"]
        
        # Matriz que recoge los reguladores correlacionados con el representante: actual.correlation. Asi se puede ver si la correlacion es
        # positiva o negativa y asignar el nombre. Intente hacer merge(), expand.grid(), pero no daba las mismas combinaciones que combn(), por lo que
        # vi necesario hacer un bucle para quedarme con aquellas parejas que interesan (representante - resto de reguladores).
        actual.couple = data.frame(t(combn(correlacionados,2)), stringsAsFactors = FALSE)
        colnames(actual.couple) = colnames(mycorrelations[,c(1,2)])
        
        actual.correlation = NULL
        for(k in 1:nrow(actual.couple)){
          if (any(actual.couple[k,c(1,2)] == keep)){
            actual.correlation = rbind(actual.correlation, actual.couple[k,])
          }
        }
        
        actual.correlation = merge(actual.correlation[,c(1,2)],mycorrelations)
        
        # Uso la matriz anterior para recorrer las correlaciones y segun sea positiva
        # o negativa, asigno un nombre.
        for(k in 1:nrow(actual.correlation)){
          if(actual.correlation[k,3] > 0){
            index = as.character(actual.correlation[k, which(with(actual.correlation, actual.correlation[k,c(1,2)] != keep))])
            reg.table[index, "filter"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", i, sep = ""), "P", sep = "_")
          } else{
            index = as.character(actual.correlation[k, which(with(actual.correlation, actual.correlation[k,c(1,2)] != keep))])
            reg.table[index, "filter"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", i, sep = ""), "N", sep = "_")
          }
        }
      } else{
        
        j=1
        
        mycomponents2= mycomponents
        
        #Opción 1: Tomar como representante el que más edges tenga y crear subgrafos separando a los que no crean conexión con el
        
        ##Repite the proccess till there are no connected edges in the subgraph
        
        while(sum(mycomponents2$csize)!=mycomponents2$no){
          
          ## Take the regulator(s) with more edges
          mynumedges=table(igraph::as_edgelist(mysubgraph))
          maxcorrelationed = names(which(mynumedges==max(mynumedges)))
          
          if(length(maxcorrelationed)>1){
            
            #Compute the sums (in absolute value) of the correlations and take as a representator the biggest
            sums = sapply(maxcorrelationed, function(x) sum(abs(mycor[which(apply( mycor[,c(1,2)]==c(x), 1, any)),3])))
            if(length(which(sums==max(sums)))>1){
              repre = sample(names(which(sums==max(sums))), 1)
            } else{
              repre = names(which(sums == max(sums)))
            }
            
          } else{
            repre = maxcorrelationed
          }
          correlacionados = names(which(igraph::as_adjacency_matrix(mysubgraph)[,repre]>0))
          regulators = colnames(data)
          
          regulators = setdiff(regulators, correlacionados)  # all regulators to keep
          data = as.matrix(data[ ,regulators])
          colnames(data) = regulators
          index.reg = which(colnames(data) == as.character(repre))
          colnames(data)[index.reg] = paste(reg.table[repre, 'omic'], paste("mc", i, sep = ""), j, "R", sep = "_")
          
          # Asignacion de nombre al representante y nueva fila para el filtro de
          # seleccion de variables (asi no se pierde la info del representante).
          reg.table = rbind(reg.table, reg.table[repre,])
          reg.table[nrow(reg.table), "regulator"] = paste(reg.table[repre, 'omic'], paste("mc", i, sep = ""),j, "R", sep = "_")
          reg.table[repre, "filter"] = paste(reg.table[repre, 'omic'], paste("mc", i, sep = ""),j, "R", sep = "_")
          rownames(reg.table) = reg.table[ ,"regulator"]
          
          # Matriz que recoge los reguladores correlacionados con el representante: actual.correlation. Asi se puede ver si la correlacion es
          # positiva o negativa y asignar el nombre. Intente hacer merge(), expand.grid(), pero no daba las mismas combinaciones que combn(), por lo que
          # vi necesario hacer un bucle para quedarme con aquellas parejas que interesan (representante - resto de reguladores).
          
          actual.correlation = NULL
          for(k in 1:nrow(mycor)){
            if (any(mycor[k,c(1,2)] == repre)){
              actual.correlation = rbind(actual.correlation, mycor[k,])
            }
          }
          
          # Uso la matriz anterior para recorrer las correlaciones y segun sea positiva
          # o negativa, asigno un nombre.
          for(k in 1:nrow(actual.correlation)){
            if(actual.correlation[k,3] > 0){
              index = as.character(actual.correlation[k, which(with(actual.correlation, actual.correlation[k,c(1,2)] != repre))])
              reg.table[index, "filter"] = paste(reg.table[repre, 'omic'], paste("mc", i, sep = ""),j, "P", sep = "_")
            } else{
              index = as.character(actual.correlation[k, which(with(actual.correlation, actual.correlation[k,c(1,2)] != repre))])
              reg.table[index, "filter"] = paste(reg.table[repre, 'omic'], paste("mc", i, sep = ""),j, "N", sep = "_")
            }
          }
          mysubgraph<-igraph::delete_vertices(mysubgraph,correlacionados)
          mycomponents2 = igraph::components(mysubgraph)
          j=j+1
          
        }
        
        
      }
      
    }
  }
  
  resultado = list(RegulatorMatrix = data, SummaryPerTargetF = reg.table)
  rownames(resultado$SummaryPerTargetF) = resultado$SummaryPerTargetF[,"regulator"]
  return(resultado)
}


## Multicollinearity filter: Partial Correlation full order


cor2pcor<-function(m,tol){
  m1 = try(MASS::ginv(m, tol = tol),silent=TRUE)
  if (class(m1)[1]=='try-error'){
    return(matrix(NA, ncol = ncol(m),nrow=nrow(m)))
  } else{
    diag(m1) = diag(m1)
    return(-cov2cor(m1))
  }
}

partialcorrelation <- function(data,reg.table,myreg, omic.type,epsilon){
  
  data<-as.matrix(data)
  pairs<-combn(ncol(data),2)
  
  cvx <- matrix(0,ncol = ncol(data),nrow=ncol(data))
  
  for (i in 1:ncol(pairs)) {
    cvx[pairs[,i][1],pairs[,i][2]] = correlations(pairs[,i], data, reg.table, omic.type)
    cvx[pairs[,i][2],pairs[,i][1]] = cvx[pairs[,i][1],pairs[,i][2]]
  }
  
  diag(cvx)<-1
  
  # partial correlation
  correlation <- cor2pcor(cvx, epsilon)
  
  colnames(correlation)<-colnames(data)
  rownames(correlation)<-colnames(data)
  
  correlation<-data.frame(t(combn(myreg,2)),combn(myreg, 2, function(x) correlation[x[1], x[2]]))
  
  return(correlation)
}

CollinearityFilter2 = function(data, reg.table, correlation = 0.8, omic.type,epsilon,scale,center) {
  
  ## data = Regulator data matrix for all omics where missing values and regulators with low variation have been filtered out
  #         (regulators must be in columns)
  ## reg.table = Table with "targetF", "regulator", "omic", "area", filter" where omics with no regulators have been removed
  row.names(reg.table) = reg.table[,"regulator"]
  #resultado = list(RegulatorMatrix = data, SummaryPerTargetF = reg.table)
  
  myreg = as.character(reg.table[which(reg.table[,"filter"] == "Model"),"regulator"])
  data<-data[,myreg]
  #Scale only the data for correlation calculation
  data2 = scale(data,scale,center)
  mycorrelations = suppressWarnings(partialcorrelation(data2, reg.table,myreg, omic.type, epsilon))
  
  if(any(is.na(mycorrelations[,3]))){
    return(NULL)
  }
  
  ## Compute the correlation between all regulators (even if they are of different omics)
  mycor = mycorrelations[abs(mycorrelations[,3]) >= correlation,]
  
  if (nrow(mycor) == 1) {  ### only 2 regulators are correlated in this omic
    
    correlacionados = unlist(mycor[,1:2])
    regulators = colnames(data)
    keep = sample(correlacionados, 1) # Regulador al azar de la pareja
    
    ## Lo siguiente elimina el no representante de la matriz de reguladores. Al regulador escogido como representante,
    ## le cambia el nombre por "mc_1_R" para que despues pase la seleccion de variables y asi, en reg.table se conserva
    ## la info de que fue escogido como representante.
    
    remove = setdiff(correlacionados, keep)
    regulators = setdiff(regulators, remove)
    data = as.matrix(data[ ,regulators])
    colnames(data) = regulators
    index.reg = which(colnames(data) == as.character(keep))
    colnames(data)[index.reg] = paste(reg.table[keep[[1]], 'omic'], paste("mc", 1, sep = ""), "R", sep = "_")
    
    # Cambio en reg.table. Asignacion de los nombres segun sea representante,
    # correlacion positiva o negativa. Creacion de una nueva fila con el representante
    # para la seleccion de variables y asi, no perder la info del representante.
    
    reg.table = rbind(reg.table, reg.table[keep,])
    reg.table[nrow(reg.table), "regulator"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", 1, sep = ""), "R", sep = "_")
    reg.table[keep, "filter"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", 1, sep = ""), "R", sep = "_")
    rownames(reg.table) = reg.table[ ,"regulator"]
    
    if(mycor[,3] > 0){
      reg.table[remove, "filter"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", 1, sep = ""), "P", sep = "_")
    } else{
      reg.table[remove, "filter"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", 1, sep = ""), "N", sep = "_")
    }
  }
  
  if (nrow(mycor) >= 2) {   ### more than 2 regulators might be correlated in this omic
    
    mygraph = igraph::graph_from_data_frame(mycor, directed=F)
    mycomponents = igraph::clusters(mygraph)
    mygraph$community<-mycomponents$membership ##save membership information
    
    for (i in 1:mycomponents$no) {
      
      #create the subgraphs of the clusters
      mysubgraph = igraph::subgraph(mygraph,as.numeric(igraph::V(mygraph)[which(mygraph$community==i)]))
      
      nedges = igraph::ecount(mysubgraph)
      
      ## see if it is a fully connected graph
      if (nedges == ((mycomponents$csize[i]*(mycomponents$csize[i]-1))/2)){
        correlacionados = names(mycomponents$membership[mycomponents$membership == i])
        regulators = colnames(data)
        
        ## Escoge un regulador al azar como representante de cada componente conexa. Para cada componente conexa elimina aquellos reguladores
        ## que no han sido escogidos como representante.
        keep = sample(correlacionados, 1)  # mantiene uno al azar
        reg.remove = setdiff(correlacionados, keep) # correlated regulator to remove
        regulators = setdiff(regulators, reg.remove)  # all regulators to keep
        data = as.matrix(data[ ,regulators])
        colnames(data) = regulators
        index.reg = which(colnames(data) == as.character(keep))
        colnames(data)[index.reg] = paste(reg.table[keep[[1]], 'omic'], paste("mc", i, sep = ""), "R", sep = "_")
        
        # Asignacion de nombre al representante y nueva fila para el filtro de
        # seleccion de variables (asi no se pierde la info del representante).
        reg.table = rbind(reg.table, reg.table[keep,])
        reg.table[nrow(reg.table), "regulator"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", i, sep = ""), "R", sep = "_")
        reg.table[keep, "filter"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", i, sep = ""), "R", sep = "_")
        rownames(reg.table) = reg.table[ ,"regulator"]
        
        # Matriz que recoge los reguladores correlacionados con el representante: actual.correlation. Asi se puede ver si la correlacion es
        # positiva o negativa y asignar el nombre. Intente hacer merge(), expand.grid(), pero no daba las mismas combinaciones que combn(), por lo que
        # vi necesario hacer un bucle para quedarme con aquellas parejas que interesan (representante - resto de reguladores).
        actual.couple = data.frame(t(combn(correlacionados,2)), stringsAsFactors = FALSE)
        colnames(actual.couple) = colnames(mycorrelations[,c(1,2)])
        
        actual.correlation = NULL
        for(k in 1:nrow(actual.couple)){
          if (any(actual.couple[k,c(1,2)] == keep)){
            actual.correlation = rbind(actual.correlation, actual.couple[k,])
          }
        }
        
        actual.correlation = merge(actual.correlation[,c(1,2)],mycorrelations)
        
        # Uso la matriz anterior para recorrer las correlaciones y segun sea positiva
        # o negativa, asigno un nombre.
        for(k in 1:nrow(actual.correlation)){
          if(actual.correlation[k,3] > 0){
            index = as.character(actual.correlation[k, which(with(actual.correlation, actual.correlation[k,c(1,2)] != keep))])
            reg.table[index, "filter"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", i, sep = ""), "P", sep = "_")
          } else{
            index = as.character(actual.correlation[k, which(with(actual.correlation, actual.correlation[k,c(1,2)] != keep))])
            reg.table[index, "filter"] = paste(reg.table[keep[[1]], 'omic'], paste("mc", i, sep = ""), "N", sep = "_")
          }
        }
      } else{
        
        j=1
        
        mycomponents2= mycomponents
        
        #Opción 1: Tomar como representante el que más edges tenga y crear subgrafos separando a los que no crean conexión con el
        
        ##Repite the proccess till there are no connected edges in the subgraph
        
        while(sum(mycomponents2$csize)!=mycomponents2$no){
          
          ## Take the regulator(s) with more edges
          mynumedges=table(igraph::as_edgelist(mysubgraph))
          maxcorrelationed = names(which(mynumedges==max(mynumedges)))
          
          if(length(maxcorrelationed)>1){
            
            #Compute the sums (in absolute value) of the correlations and take as a representator the biggest
            sums = sapply(maxcorrelationed, function(x) sum(abs(mycor[which(apply( mycor[,c(1,2)]==c(x), 1, any)),3])))
            if(length(which(sums==max(sums)))>1){
              repre = sample(names(which(sums==max(sums))), 1)
            } else{
              repre = names(which(sums == max(sums)))
            }
            
          } else{
            repre = maxcorrelationed
          }
          
          correlacionados = names(which(igraph::as_adjacency_matrix(mysubgraph)[,repre]>0))
          regulators = colnames(data)
          
          regulators = setdiff(regulators, correlacionados)  # all regulators to keep
          data = as.matrix(data[ ,regulators])
          colnames(data) = regulators
          index.reg = which(colnames(data) == as.character(repre))
          colnames(data)[index.reg] = paste(reg.table[repre, 'omic'], paste("mc", i, sep = ""), j, "R", sep = "_")
          
          # Asignacion de nombre al representante y nueva fila para el filtro de
          # seleccion de variables (asi no se pierde la info del representante).
          reg.table = rbind(reg.table, reg.table[repre,])
          reg.table[nrow(reg.table), "regulator"] = paste(reg.table[repre, 'omic'], paste("mc", i, sep = ""),j, "R", sep = "_")
          reg.table[repre, "filter"] = paste(reg.table[repre, 'omic'], paste("mc", i, sep = ""),j, "R", sep = "_")
          rownames(reg.table) = reg.table[ ,"regulator"]
          
          # Matriz que recoge los reguladores correlacionados con el representante: actual.correlation. Asi se puede ver si la correlacion es
          # positiva o negativa y asignar el nombre. Intente hacer merge(), expand.grid(), pero no daba las mismas combinaciones que combn(), por lo que
          # vi necesario hacer un bucle para quedarme con aquellas parejas que interesan (representante - resto de reguladores).
          
          actual.correlation = NULL
          for(k in 1:nrow(mycor)){
            if (any(mycor[k,c(1,2)] == repre)){
              actual.correlation = rbind(actual.correlation, mycor[k,])
            }
          }
          
          # Uso la matriz anterior para recorrer las correlaciones y segun sea positiva
          # o negativa, asigno un nombre.
          for(k in 1:nrow(actual.correlation)){
            if(actual.correlation[k,3] > 0){
              index = as.character(actual.correlation[k, which(with(actual.correlation, actual.correlation[k,c(1,2)] != repre))])
              reg.table[index, "filter"] = paste(reg.table[repre, 'omic'], paste("mc", i, sep = ""),j, "P", sep = "_")
            } else{
              index = as.character(actual.correlation[k, which(with(actual.correlation, actual.correlation[k,c(1,2)] != repre))])
              reg.table[index, "filter"] = paste(reg.table[repre, 'omic'], paste("mc", i, sep = ""),j, "N", sep = "_")
            }
          }
          mysubgraph<-igraph::delete_vertices(mysubgraph,correlacionados)
          mycomponents2 = igraph::clusters(mysubgraph)
          j=j+1
          
        }
        
        
      }
      
    }
  }
  
  resultado = list(RegulatorMatrix = data, SummaryPerTargetF = reg.table)
  rownames(resultado$SummaryPerTargetF) = resultado$SummaryPerTargetF[,"regulator"]
  return(resultado)
}

