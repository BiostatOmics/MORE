#########################################################################################
######           Functions to integrate omics data using PLSs                      ######
#########################################################################################


## By Maider
## 20-July-2023
## Last modified: October 2023


options(stringsAsFactors = FALSE)

library(ropls)

#'
#'\code{GetPLS} fits a PLS model for all the features in the target omic dataset to identify
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
#' @param percNA Maximum percentage of missing values present in regulatoryData regulators and observations to use them to construct the models.
#' @param scaleType Type of scaling to be applied. Four options:
#' \itemize{
#' \item auto : Applies the auto scaling. 
#' \item softBlock : Applies the pareto scaling. \deqn{\frac{X_k}{s_k \sqrt[4]{m_b}} }
#' \item hardBlock : Applies the block scaling. \deqn{ \frac{X_k}{s_k \sqrt{m_b}} }
#' \item none : It does not apply any type of scaling. Not recommended if the user does not apply their own scaling.
#' }
#' considering m_b the number of variables of the block. By default, auto.
#' @param interactions If TRUE, the model includes interactions between regulators and condition variables. By default, TRUE.
#' @param varSel Type of variable selection method to apply, different options depending on MLR or PLS usage. Four options:
#' \itemize{
#' \item EN : Applies a Multiple Linear Regression (MLR) with ElasticNet regularization.
#' \item ISGL : Applies a Multiple Linear Regression (MLR) with Iterative Sparse Group Lasso (ISGL) regularization.
#' \item Jack : Applies Jack-Knife resampling technique for the calculation of the significance of the coefficients in Partial Least Squares (PLS) models.
#' \item Perm : Applies a resampling technique for the calculation of the significance of the coefficients in Partial Least Squares (PLS) models in which the response variable is permuted 100 times to obtain the distribution of the coefficients and compute then their associated p-value.
#' }
#' By default, Jack.
#' @param alfa Significance level for variable selection in PLS1 and PLS2 \code{\link{method}}. By default, 0.05.
#' @param vip Value of VIP above which a variable can be considered significant in addition to the computed p-value in \code{\link{varSel}}. By default, 0.8.
#' @param method Model to be fitted. Two options:
#' \itemize{
#' \item PLS1 : Applies a Partial Least Squares (PLS) model, one for each of the features in the target omic of \code{\link{targetData}}.
#' \item PLS2 : Applies a PLS model to all features of the target omic simultaneously, only possible when \code{\link{associations}}= NULL.
#' }
#' By default, PLS1.
#' @param parallel If FALSE, MORE will be run sequentially. If TRUE, MORE will be run using parallelization with as many cores as the available ones minus one so the system is not overcharged. If the user wants to specify how many cores they want to use, they can also provide the number of cores to use in this parameter. Parallelization is only implemented for MLR with EN variable selection and PLS methods.
#' @return List containing the following elements:
#' \itemize{
#' \item ResultsPerTargetF : List with as many elements as features of the target omic in \code{\link{targetData}}. For each feature, it includes information about the feature values, considered variables, estimated coefficients,
#'                    detailed information about all regulators, and regulators identified as relevant (in MLR scenario) or significant (in PLS scenarios).
#' \item GlobalSummary : List with information about the fitted models, including model metrics, information about regulators, features of the target omic without models, regulators, master regulators and hub target features.
#' \item Arguments : List containing all the arguments used to generate the models.                
#' }
#'
#' @export
GetPLS = function(targetData,
                  regulatoryData,
                  associations =NULL,
                  omicType = 0,
                  condition = NULL,
                  clinic = NULL,
                  clinicType = NULL,
                  alfa = 0.05, 
                  interactions = TRUE,
                  minVariation = 0,
                  percNA = 0.2,
                  scaleType = 'auto',
                  varSel ='Jack',
                  vip = 0.8,
                  method = 'PLS1',
                  parallel = FALSE){
  
  # Converting matrix to data.frame
  targetData = as.data.frame(targetData)
  regulatoryData = lapply(regulatoryData, as.data.frame)
  
  ##Omic types
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
  
  
  # Not possible to apply a PLS2 when associations matrix is provided
  if(!is.null(associations)){
    if(method=='PLS2'){
      warning('WARNING: In PLS2 associations will only be used to reduce regulatoryData to those regulators present in associations, it won’t be used to set specific target feature regulator associations')
      #Reduce regulatoryData to the ones present in associations
      regulatoryData = lapply(names(associations), function(x) regulatoryData[[x]][rownames(regulatoryData[[x]]) %in% associations[[x]][, 2], ,drop=FALSE])
      names(regulatoryData) = names(associations)
      #Set to NULL associations ones its not needed anymore
      associations = NULL
      }
  }
  
  # If associations is NULL create a list of associations NULL
  if (is.null(associations)){
    associations=vector('list',length(regulatoryData))
    names(associations)=names(regulatoryData)
  }
  
  # Checking that the number of samples per omic is equal to number of samples for targetData and the number of samples for condition
  for (i in 1:length(names(regulatoryData))){
    if(!ncol(regulatoryData[[i]]) == ncol(targetData)) {
      stop("ERROR: Samples in regulatoryData must be the same as in targetData and in condition")
    }
  }
  if(!is.null(condition)){
    if(!ncol(targetData) == nrow(condition)) {
      stop("ERROR: Samples in regulatoryData must be the same as in targetData and in condition")
    }
  }
  
  ##Verify that the scaleType selected is a possible selection
  
  if(!scaleType %in% c('auto','softBlock','hardBlock','none')){stop('ERROR: The selected type of scaling is not one of auto, softBlock, hardBlock or none')}
  if(!varSel %in% c('Perm','Jack')){stop('ERROR: The selected method for p value calculation is not one of Perm or Jack')}
  if(!method %in% c('PLS1','PLS2')){stop('ERROR: Not valid pls type. Select one of PLS1 or PLS2')}
  
  ## Checking that samples are in the same order in targetData, regulatoryData and condition
  orderproblem<-FALSE
  if(is.null(condition)){
    nameproblem<-!all(sapply(regulatoryData, function(x) length(intersect(colnames(x),colnames(targetData))==ncol(targetData))))
    if(nameproblem){
      cat('Warning. targetData and regulatoryData samples have not same names. We assume that they are ordered.\n')
      #To avoid problems set in all the cases the same names according to the ones in targetData
      Tnames = colnames(targetData)
      regulatoryData = lapply(regulatoryData, function(x){
        colnames(x) = Tnames
        x
      })
      
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
      #To avoid problems set in all the cases the same names according to the ones in targetData
      Tnames = colnames(targetData)
      regulatoryData = lapply(regulatoryData, function(x){
        colnames(x) = Tnames
        x
      })
      rownames(condition) = Tnames 
    } else{
      orderproblem<-!all(c(sapply(regulatoryData, function(x) identical(colnames(x),colnames(targetData))), identical(colnames(targetData),rownames(condition))))
      if(orderproblem){
        regulatoryData<-lapply(regulatoryData, function(x) x[,colnames(targetData)])
        condition<-condition[colnames(targetData), , drop=FALSE]
      }
    }
  }
  
  ## Checking if there are regulators with "_R", "_P" or "_N" or ":" and checking that there are not replicate identifiers compared to targetData
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
      associations[[i]][[2]] = gsub(':', '-', associations[[i]][[2]])
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
          rownames(regulatoryData[[i]])[rownames(regulatoryData[[i]])%in%repeated] =  paste(names(regulatoryData)[i],'-', rownames(regulatoryData[[i]])[rownames(regulatoryData[[i]])%in%repeated],sep = '')
          rownames(regulatoryData[[j]])[rownames(regulatoryData[[j]])%in%repeated] =  paste(names(regulatoryData)[j],'-', rownames(regulatoryData[[j]])[rownames(regulatoryData[[j]])%in%repeated],sep = '')
          
        }
      }
    }
    rm(repeated)
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
  
  ## ropls does not allow NAs in response variable, so remove target features with missing values
  min.obs = ncol(targetData)
  # target features
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
    des.mat = model.matrix(~0+., data = as.data.frame(Group))
    #Change the name to avoid conflicts with RegulationPerCondition
    colnames(des.mat) = sub('Group','Group_',colnames(des.mat))
  }
  
  ## Removing regulators with more than percNA of NAs and keeping track
  cat("Removing regulators with missing values...\n")
  
  missing_rows <- function(regOmic,percNA) {
    highNA = apply(regOmic, 1, function(x) sum(is.na(x))/length(x)) > percNA
    return(rownames(regOmic)[highNA])
  }
  
  myregNA = lapply(regulatoryData, function(x) missing_rows(x,percNA))
  regulatoryData = lapply(1:length(regulatoryData), function (i) regulatoryData[[i]][setdiff(rownames(regulatoryData[[i]]),myregNA[[i]]),])
  names(regulatoryData)=names(omicType)
  cat("Number of regulators with missing values:\n")
  print(sapply(myregNA, length))
  cat("\n")
  
  #Filter observations with two many missing values
  missing_cols <- function(regOmic,percNA) {
    highNA = apply(regOmic, 2, function(x) sum(is.na(x))/length(x)) > percNA
    return(colnames(regOmic)[highNA])
  }
  
  myobsNA = lapply(regulatoryData, function(x) missing_cols(x,percNA))
  myobsNA = unique(unlist(myobsNA))
  #Remove from the regulatoryData, targetData, Group if available the observations with missing values
  obsNotNA = setdiff(colnames(targetData), myobsNA)
  regulatoryData = lapply(regulatoryData, function(x) x[,obsNotNA])
  cat("Number of observations with missing values:", length(myobsNA))
  cat("\n")
  
  ## Remove observation with more than percNA of NAs
  targetData = targetData[,obsNotNA]
  #Remove from condition and regulators
  des.mat = des.mat[obsNotNA,]
  Group = Group[obsNotNA]
  
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
  
  GlobalSummary$GoodnessOfFit = matrix(NA, ncol = 6, nrow = ntargetFs)
  rownames(GlobalSummary$GoodnessOfFit) = AlltargetFs
  colnames(GlobalSummary$GoodnessOfFit) = c( "RsquaredY", "Qsquared","RMSE","NRMSE","ncomp","sigReg")
  
  GlobalSummary$ReguPerTargetF = matrix(0, ncol = 3*length(regulatoryData), nrow = ntargetFs)
  rownames(GlobalSummary$ReguPerTargetF) = AlltargetFs
  colnames(GlobalSummary$ReguPerTargetF) = c(paste(names(regulatoryData), "Ini", sep = "-"),
                                          paste(names(regulatoryData), "Mod", sep = "-"),
                                          paste(names(regulatoryData), "Sig", sep = "-"))
  rm(infproblemtargetF);rm(infproblemreg);rm(problema);rm(problemas);rm(orderproblem);rm(nameproblem);rm(targetFsNA);rm(targetFsInf);rm(targetFsNOreg);rm(targetFsNotNA);rm(constantTargetFs);rm(notConstant);gc()
  ## Specific results for each target feature
  ResultsPerTargetF=vector("list", length=length(AlltargetFs))
  names(ResultsPerTargetF) = AlltargetFs
  
  ## Specify the scaling type
  scale <- ifelse(scaleType == 'none', FALSE, TRUE)
  center <- ifelse(scaleType == 'none', FALSE, TRUE)
 
  if(method=='PLS1'){
    ### Computing model for each target feature
    cat("Selecting predictors and fitting model for ...\n")
    
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
        parallel =nc
      }
      ResultsPerTargetF <- furrr::future_map(1:ntargetFs,
                                          ~ResultsPerTargetF.i(AlltargetFs[.],GlobalSummary,regulatoryData,associations,targetData,omicType,
                                                            condition, des.mat,myregLV,myregNA,scale,center,scaleType,
                                                            interactions,varSel,vip,alfa),
                                          .progress = TRUE,.options = furrr::furrr_options(seed = TRUE) )
    } else{
      ResultsPerTargetF <- purrr::map(1:ntargetFs,
                                          ~ResultsPerTargetF.i(AlltargetFs[.],GlobalSummary,regulatoryData,associations,targetData,omicType,
                                                            condition, des.mat,myregLV,myregNA,scale,center,scaleType,
                                                            interactions,varSel,vip,alfa),
                                          .progress = TRUE,.options = furrr::furrr_options(seed = TRUE) )
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
    
    
  } else{
    
    cat("Fitting model ...\n")
    
    #When associations = NULL all target features have the same potential regulators
    targetF = AlltargetFs[1]
    RetRegul = GetAllReg(targetF=targetF, associations=associations, data.omics = regulatoryData)
    RetRegul.targetF = RetRegul$Results  ## RetRegul$TableTargetF: nr reg per omic
    ## Some of these reg will be removed, because they are not in regulatoryData
    
    #Add information about initial potential regulators
    GlobalSummary$ReguPerTargetF[1:nrow(GlobalSummary$ReguPerTargetF), grep("-Ini", colnames(GlobalSummary$ReguPerTargetF))]=rep(as.numeric(RetRegul$TableTargetF[-1]),each=length(AlltargetFs))
    
    # RetRegul.targetF--> targetF/regulator/omic/area
    RetRegul.targetF=RetRegul.targetF[RetRegul.targetF[,"regulator"]!= "No-regulator", ,drop=FALSE] ## Remove rows with no-regulators
    
    allRegulators = data.frame(RetRegul.targetF, rep("Model",nrow(RetRegul.targetF)), stringsAsFactors = FALSE)
    colnames(allRegulators) = c("targetF","regulator","omic","area","filter")
    
    ## Identify which regulators where removed because of missing values or low variation
    res = RemovedRegulators(RetRegul.targetF = allRegulators,
                            myregLV=myregLV, myregNA=myregNA, data.omics=regulatoryData)
    rm(myregLV);rm(myregNA);rm(RetRegul.targetF);rm(RetRegul);gc()
    ## Create the interactions between regulators and condition 
    
    des.mat2 = RegulatorsInteractions(interactions, reguValues = res$RegulatorMatrix, reguInfo = res$SummaryPerTargetF,
                                          des.mat, method ='pls')
    
    #Scale the variables, necessary for the models interpretation
    des.mat2 = lapply(des.mat2, function(x)scale(x,scale=scale,center=center))
    
    ##Scale if needed for soft or hard block scaling
    res$RegulatorMatrix = Scaling.type(des.mat2,scaleType)
    if(is.null(des.mat)){
      des.mat2 = res$RegulatorMatrix
    } else{
      des.mat2 = data.frame(scale(des.mat,scale=scale,center=center), res$RegulatorMatrix, check.names = FALSE)
    }
    rm(res); gc()
    
    # Removing predictors with constant values
    sdNo0 = apply(des.mat2, 2, sd)
    sdNo0 = names(sdNo0)[sdNo0 > 0]
    quito = setdiff(colnames(des.mat2), sdNo0)
    des.mat2 = des.mat2[,sdNo0]
    sdNo0 = unique(gsub(".*[:?](.*)$", "\\1", sdNo0))
    quito = unique(gsub(".*[:?](.*)$", "\\1", quito))
    quito = setdiff(quito, sdNo0)
    allRegulators[quito,"filter"] = "Constant"
    
    Y = scale(t(targetData), center = center, scale = scale)
    
    ## Computing PLS model
    if (nrow(des.mat2)<7){cross = nrow(des.mat)-2}else{cross =7}
    myPLS = try(suppressWarnings( ropls::opls(des.mat2, Y, info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0)),silent = TRUE)
    
    if(class(myPLS)=='try-error' || length(myPLS@modelDF)==0){
      myPLS = try(suppressWarnings( ropls::opls(des.mat2, Y, info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0, predI=1)), silent = TRUE)
    }
    
    if (class(myPLS)=='try-error' || length(myPLS@modelDF)==0){
      
      myPLS = NULL
      
      GlobalSummary$TargetFNOmodel = rbind(GlobalSummary$TargetFNOmodel,
                                         data.frame("targetF" = AlltargetFs, "problem" = "No significant components on PLS2"))
      
      ## Extracting significant regulators
      for (i in 1:length(AlltargetFs)) {
        targetF = AlltargetFs[i]
        ResultsPerTargetF[[i]]$significantRegulators = NULL
        ResultsPerTargetF[[i]]$allRegulators = data.frame(ResultsPerTargetF[[i]]$allRegulators, "Sig" = NA, stringsAsFactors = FALSE)
        
        ResultsPerTargetF[[i]]$allRegulators = allRegulators
        ResultsPerTargetF[[i]]$allRegulators[,'targetF']=rep(targetF,nrow(ResultsPerTargetF[[i]]$allRegulators))
        
        ## Counting original regulators in the model per omic
        contando = ResultsPerTargetF[[i]]$allRegulators[which(ResultsPerTargetF[[i]]$allRegulators[,"filter"] == "Model"),]
        contando = table(contando[,"omic"])
        contando = as.numeric(contando[names(regulatoryData)])
        contando[is.na(contando)] = 0
        GlobalSummary$ReguPerTargetF[targetF, grep("-Mod", colnames(GlobalSummary$ReguPerTargetF))] = contando
        
        ## Counting significant regulators per omic
        GlobalSummary$ReguPerTargetF[targetF, grep("-Sig", colnames(GlobalSummary$ReguPerTargetF))] = NA
        
        ResultsPerTargetF[[i]]$Y = targetData[i,]
        ResultsPerTargetF[[i]]$coefficients = NULL
        ResultsPerTargetF[[i]]$X = des.mat2
        
        GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[rownames(GlobalSummary$GoodnessOfFit) != targetF,]
        
      }
      
      
    }else{
      
      myPLS = suppressInnecPLSdata(myPLS)
      if (varSel == 'Jack'){
        pval = p.valuejack.pls2(myPLS, des.mat2, Y, alfa)
      } else {
        pval = p.coef.pls2(myPLS, 100, des.mat2, Y)
      }
      
      #Aunque lo hayamos tratado todo junto ahora separamos la respuesta por target features
      
      for (i in 1:ncol(Y)) {
        
        #Save X matrix
        ResultsPerTargetF[[i]]$X = des.mat2
        
        #Tratar como significativas tan solo las que cumplan ambas condiciones
        sigvariables = intersect(names(myPLS@vipVn[which(myPLS@vipVn>vip)]), rownames(pval[,i,drop=FALSE])[which(pval[,i,drop=FALSE]<alfa)])
        ResultsPerTargetF[[i]]$coefficients = data.frame('coefficient' = myPLS@coefficientMN[sigvariables,i], 'pvalue' = pval[sigvariables,i,drop=FALSE])

        ## Extracting significant regulators and recovering correlated regulators
        myvariables = unlist(strsplit(sigvariables, ":", fixed = TRUE))
        myvariables = intersect(myvariables, allRegulators[,'regulator'])
        
        ResultsPerTargetF[[i]]$allRegulators = allRegulators
        
        targetF = AlltargetFs[i]
        ResultsPerTargetF[[i]]$allRegulators[,'targetF']=rep(targetF,nrow(ResultsPerTargetF[[i]]$allRegulators))
        
        ResultsPerTargetF[[i]]$allRegulators = data.frame(ResultsPerTargetF[[i]]$allRegulators, "Sig" = 0, stringsAsFactors = FALSE)
        ResultsPerTargetF[[i]]$allRegulators[myvariables, "Sig"] = 1
        
        ResultsPerTargetF[[i]]$significantRegulators = myvariables
        
        ## Counting original regulators in the model per omic              
        contando = ResultsPerTargetF[[i]]$allRegulators
        quitar = which(contando[,"filter"] == "MissingValue")
        if (length(quitar) > 0) contando = contando[-quitar,]
        quitar = which(contando[,"filter"] == "LowVariation")
        if (length(quitar) > 0) contando = contando[-quitar,]
        
        contando = table(contando[,"omic"])
        contando = as.numeric(contando[names(regulatoryData)])
        contando[is.na(contando)] = 0
        GlobalSummary$ReguPerTargetF[targetF, grep("-Mod", colnames(GlobalSummary$ReguPerTargetF))] = contando
        
        ## Counting significant regulators per omic
        if (length(ResultsPerTargetF[[i]]$significantRegulators) > 0) {
          contando = ResultsPerTargetF[[i]]$allRegulators[ResultsPerTargetF[[i]]$significantRegulators,]
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(regulatoryData)])
          contando[is.na(contando)] = 0
          GlobalSummary$ReguPerTargetF[targetF, grep("-Sig", colnames(GlobalSummary$ReguPerTargetF))] = contando
        } else{
          GlobalSummary$TargetFNOmodel = rbind(GlobalSummary$TargetFNOmodel,
                                               data.frame("targetF" = targetF, "problem" = "No significant regulators after variable selection"))
        }
        
        ResultsPerTargetF[[i]]$Y = data.frame("y" = myPLS@suppLs$y[,i,drop=FALSE], "fitted.y" = myPLS@suppLs$yPreMN[,i,drop=FALSE],
                                           "residuals" = myPLS@suppLs$y[,i,drop=FALSE]-myPLS@suppLs$yPreMN[,i,drop=FALSE])
        colnames( ResultsPerTargetF[[i]]$Y) = c('y','fitted.y','residuals')
        
        
        GlobalSummary$GoodnessOfFit[targetF,] = c(myPLS@modelDF[,'R2Y(cum)'][myPLS@summaryDF[,'pre']],
                                               myPLS@modelDF[,'Q2(cum)'][myPLS@summaryDF[,'pre']],
                                               myPLS@summaryDF[,'RMSEE'],
                                               round(myPLS@summaryDF[,'RMSEE']/(max(myPLS@suppLs$y)-min(myPLS@suppLs$y)),6),
                                               round(myPLS@summaryDF[,'pre'],digits = 0),
                                               round(length(ResultsPerTargetF[[i]]$significantRegulators),digits = 0))
      }
      
      
    }
    
  } 

  targetFsNosig = names(which(GlobalSummary$GoodnessOfFit[,6]==0))
  targetFssig = setdiff(rownames(GlobalSummary$GoodnessOfFit), targetFsNosig)
  GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[targetFssig,]
  
  targetFsNoreg = rownames(GlobalSummary$GoodnessOfFit)[is.na(rowSums(GlobalSummary$GoodnessOfFit))]
  targetFsreg = setdiff(rownames(GlobalSummary$GoodnessOfFit), targetFsNoreg)
  GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[targetFsreg,]
  
  #Calculate GlobalRegulators
  m_sig_reg<-lapply(ResultsPerTargetF, function(x) x$significantRegulators)
  m_sig_reg <- unlist(m_sig_reg)
  msig_vector <- table(m_sig_reg)
  #Calculate third quantile
  q3<-quantile(msig_vector,0.75)
  if(length(msig_vector[msig_vector>q3])<10){
    GlobalSummary$GlobalRegulators = intersect(names(msig_vector[rev(tail(order(msig_vector),10))]), names(msig_vector[msig_vector>10]) )
  } else{
    GlobalSummary$GlobalRegulators = intersect(names(msig_vector[msig_vector>q3]), names(msig_vector[msig_vector>10]) ) 
  }
  
  #Calculate HubTargetF
  significant_regulators<-GlobalSummary$ReguPerTargetF[,c(grep('-Sig$',colnames(GlobalSummary$ReguPerTargetF))),drop=FALSE]
  s_sig_reg<-apply(significant_regulators, 1, sum)
  #Calculate third quantile
  q3<-quantile(s_sig_reg,0.75)
  if(length(s_sig_reg[s_sig_reg>q3])<10){
    GlobalSummary$HubTargetF = intersect(names(s_sig_reg[rev(tail(order(s_sig_reg),10))]), names(s_sig_reg[s_sig_reg>10]) )
  } else{
    GlobalSummary$HubTargetF = intersect(names(s_sig_reg[s_sig_reg>q3]), names(s_sig_reg[s_sig_reg>10]))
  }
  
  if(all(sapply(associations,is.null))) {associations = NULL}
  
  myarguments = list(condition = condition, finaldesign = des.mat, groups = Group, alfa = alfa, 
                     clinicType = clinicType, minVariation = minVariation, associations = associations, vip = vip,
                     targetData = targetData, regulatoryData = regulatoryData, omicType = omicType,
                     clinic = clinic, scaleType = scaleType, varSel=varSel, method =method, parallel = parallel)

  # Create the results for the scale filter check
  
  result <- list("ResultsPerTargetF" = ResultsPerTargetF, "GlobalSummary" = GlobalSummary, "arguments" = myarguments) 
  class(result) <- "MORE"
  return(result)
}



ResultsPerTargetF.i<-function(targetF,GlobalSummary,regulatoryData,associations,targetData,omicType,
                           condition, des.mat,myregLV,myregNA,scale,center,scaleType,
                           interactions,varSel,vip,alfa){
  
  
    ResultsPerTargetF.i = vector("list", length = 9)
    names(ResultsPerTargetF.i) = c("Y", "X", "coefficients", "allRegulators", "significantRegulators", "GoodnessOfFit", "ReguPerTargetF","TargetFNOmodel","TargetFNOregu")
    
    #Initialize gobal summary values
    ResultsPerTargetF.i$ReguPerTargetF <- GlobalSummary$ReguPerTargetF[targetF, , drop=FALSE]
    
    RetRegul = GetAllReg(targetF=targetF, associations=associations, data.omics = regulatoryData)
    RetRegul.targetF = RetRegul$Results  ## RetRegul$TableTargetF: nr reg per omic
    ## Some of these reg will be removed, because they are not in regulatoryData
    
    
    # RetRegul.targetF--> targetF/regulator/omic/area
    RetRegul.targetF=RetRegul.targetF[RetRegul.targetF[,"regulator"]!= "No-regulator", ,drop=FALSE] ## Remove rows with no-regulators
    
    
    ### NO INITIAL REGULATORS
    if(length(RetRegul.targetF)==0){ ## En el caso de que no hayan INICIALMENTE reguladores -> Calcular modelo con variables experimentales
      
      if (is.null(condition)) {
        ResultsPerTargetF.i$X = NULL
        ResultsPerTargetF.i$significantRegulators = NULL
        ResultsPerTargetF.i$allRegulators = NULL
        myPLS = NULL
        
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
        
        myPLS = NULL
        
        ResultsPerTargetF.i$X = des.mat2[,-1, drop = FALSE]
        ResultsPerTargetF.i$significantRegulators = NULL
        ResultsPerTargetF.i$allRegulators = NULL
        
        ResultsPerTargetF.i$TargetFNOregu = rbind(ResultsPerTargetF.i$TargetFNOregu,
                                                   data.frame("targetF" = targetF,
                                                              "problem" = 'Target feature had no initial regulators'))
      }
      
      # GlobalSummary$ReguPerTargetF  # this is initially set to 0 so no need to modify it
      
      ### WITH INITIAL REGULATORS
    } else { ## There are regulators for this targetF at the beginning
      
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
          ResultsPerTargetF.i$significantRegulators = NULL
          ResultsPerTargetF.i$allRegulators = res$SummaryPerTargetF
          myPLS = NULL
          
          ResultsPerTargetF.i$TargetFNOmodel = rbind(ResultsPerTargetF.i$TargetFNOmodel,
                                                     data.frame("targetF" = targetF,
                                                                "problem" = 'No regulators left after NA/LowVar filtering'))
        } else {
          des.mat2 = cbind(t(targetData[targetF,]), des.mat)
          colnames(des.mat2)[1] = "response"
          des.mat2 = na.omit(des.mat2)
          
          # Removing predictors with constant values
          sdNo0 = apply(des.mat2, 2, sd)
          sdNo0 = names(sdNo0)[sdNo0 > 0]
          des.mat2 = des.mat2[,sdNo0]
          myPLS = NULL
          
          ResultsPerTargetF.i$TargetFNOmodel = rbind(ResultsPerTargetF.i$TargetFNOmodel,
                                                     data.frame("targetF" = targetF,
                                                        "problem" = 'No regulators left after NA/LowVar filtering'))
          ResultsPerTargetF.i$X = des.mat2[,-1, drop = FALSE]
          ResultsPerTargetF.i$significantRegulators = NULL
          ResultsPerTargetF.i$allRegulators = res$SummaryPerTargetF
        }
        
      } else {  ## Regulators for the model!!
        ResultsPerTargetF.i$allRegulators = res$SummaryPerTargetF
        
        ## Create the interactions between regulators and condition 
        des.mat2 = RegulatorsInteractions(interactions, reguValues = res$RegulatorMatrix, reguInfo = res$SummaryPerTargetF,
                                          des.mat, method = 'pls')
        
        #Scale the variables, necessary for the models interpretation
        des.mat2 = lapply(des.mat2, function(x)scale(x,scale=scale,center=center))
        
        ##Scale if needed for soft or hard block scaling
        res$RegulatorMatrix = Scaling.type(des.mat2,scaleType)
        if(is.null(des.mat)){
          des.mat2 = data.frame(t(targetData[targetF,]),res$RegulatorMatrix, check.names = FALSE)
        } else{
          des.mat2 = data.frame(t(targetData[targetF,]),scale(des.mat,scale=scale,center=center), res$RegulatorMatrix, check.names = FALSE)
        }
        colnames(des.mat2)[1] = "response"
        rm(res); gc()
        
        # Removing predictors with constant values
        sdNo0 = apply(des.mat2, 2, function(x) sd(x,na.rm = FALSE))
        # Do not remove a variable for containing NA in its values
        sdNo0 = names(sdNo0)[is.na(sdNo0) | sdNo0 > 0]
        quito = setdiff(colnames(des.mat2), sdNo0)
        des.mat2 = des.mat2[,sdNo0]
        sdNo0 = unique(gsub(".*[:?](.*)$", "\\1", sdNo0))
        quito = unique(gsub(".*[:?](.*)$", "\\1", quito))
        quito = setdiff(quito, sdNo0)
        ResultsPerTargetF.i$allRegulators[quito,"filter"] = "Constant"
        
        #Save X matrix
        
        ResultsPerTargetF.i$X = des.mat2[,-1, drop = FALSE]
        
        ## Computing PLS model
        if (nrow(des.mat2)<7){cross = nrow(des.mat)-2}else{cross =7}
        myPLS = try(suppressWarnings(ropls::opls(des.mat2[,-1,drop=FALSE],scale(des.mat2[,1],scale=scale,center=center) , info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0)),silent = TRUE)
        
        if(class(myPLS)=='try-error' || length(myPLS@modelDF)==0){
          myPLS = try(suppressWarnings(ropls::opls(des.mat2[,-1,drop=FALSE], scale(des.mat2[,1],scale=scale,center=center), info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0, predI=1)), silent = TRUE)
        }
        
        if (class(myPLS)=='try-error' || length(myPLS@modelDF)==0){
          
          myPLS = NULL
          
          ResultsPerTargetF.i$TargetFNOmodel = rbind(ResultsPerTargetF.i$TargetFNOmodel,
                                             data.frame("targetF" = targetF, "problem" = "No significant components on PLS"))
          
          ## Extracting significant regulators
          ResultsPerTargetF.i$significantRegulators = NULL
          ResultsPerTargetF.i$allRegulators = data.frame(ResultsPerTargetF.i$allRegulators, "Sig" = NA, stringsAsFactors = FALSE)
          
          ## Counting original regulators in the model per omic
          contando = ResultsPerTargetF.i$allRegulators[which(ResultsPerTargetF.i$allRegulators[,"filter"] == "Model"),]
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(regulatoryData)])
          contando[is.na(contando)] = 0
          ResultsPerTargetF.i$ReguPerTargetF[targetF, grep("-Mod", colnames(ResultsPerTargetF.i$ReguPerTargetF))] = contando
          
          ## Counting significant regulators per omic
          ResultsPerTargetF.i$ReguPerTargetF[targetF, grep("-Sig", colnames(ResultsPerTargetF.i$ReguPerTargetF))] = NA
          
        }else{
          
          if (varSel == 'Jack'){
            pval = p.valuejack(myPLS, des.mat2, alfa)
          } else {
            pval = p.coef(myPLS, 100, des.mat2)
          }
          
          #Tratar como significativas tan solo las que cumplan ambas condiciones
          sigvariables = intersect(names(myPLS@vipVn[which(myPLS@vipVn>vip)]), rownames(pval)[which(pval<alfa)])
          ResultsPerTargetF.i$coefficients = data.frame('coefficient' = myPLS@coefficientMN[sigvariables,], 'pvalue' = pval[sigvariables,,drop=FALSE])
          #Recalcular los parametros con las variables seleccionadas
          if(length(sigvariables)>0) {
            fPLS = try(suppressWarnings(ropls::opls(des.mat2[,sigvariables,drop=FALSE], scale(des.mat2[,1],scale=scale,center=center), info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI = 0)),silent=TRUE)
            if(class(fPLS)=='try-error' || length(fPLS@modelDF)==0){
              fPLS = try(suppressWarnings(ropls::opls(des.mat2[,sigvariables,drop=FALSE], scale(des.mat2[,1],scale=scale,center=center), info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0, predI=1)), silent = TRUE)
            } 
            if(class(fPLS)!='try-error' && length(fPLS@modelDF)!=0){
              myPLS = fPLS}
          } 
          
          ## Extracting significant regulators and recovering correlated regulators
          myvariables = unlist(strsplit(sigvariables, ":", fixed = TRUE))
          myvariables = intersect(myvariables, rownames(ResultsPerTargetF.i$allRegulators))
          
          ## Remove if only group is significant
          if(length(myvariables)==0){
            ResultsPerTargetF.i$TargetFNOmodel = rbind(ResultsPerTargetF.i$TargetFNOmodel,
                                                       data.frame("targetF" = targetF, "problem" = "No significant regulators after variable selection"))
          }
          
          ResultsPerTargetF.i$allRegulators = data.frame(ResultsPerTargetF.i$allRegulators, "Sig" = 0, stringsAsFactors = FALSE)
          ResultsPerTargetF.i$allRegulators[myvariables, "Sig"] = 1
          
          ResultsPerTargetF.i$significantRegulators = myvariables
          
          ## Counting original regulators in the model per omic              
          contando = ResultsPerTargetF.i$allRegulators
          quitar = which(contando[,"filter"] == "MissingValue")
          if (length(quitar) > 0) contando = contando[-quitar,]
          quitar = which(contando[,"filter"] == "LowVariation")
          if (length(quitar) > 0) contando = contando[-quitar,]
          
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(regulatoryData)])
          contando[is.na(contando)] = 0
          ResultsPerTargetF.i$ReguPerTargetF[targetF, grep("-Mod", colnames(ResultsPerTargetF.i$ReguPerTargetF))] = contando
          
          ## Counting significant regulators per omic
          if (length(ResultsPerTargetF.i$significantRegulators) > 0) {
            contando = ResultsPerTargetF.i$allRegulators[ResultsPerTargetF.i$significantRegulators,]
            contando = table(contando[,"omic"])
            contando = as.numeric(contando[names(regulatoryData)])
            contando[is.na(contando)] = 0
            ResultsPerTargetF.i$ReguPerTargetF[targetF, grep("-Sig", colnames(ResultsPerTargetF.i$ReguPerTargetF))] = contando
          }
          
          
        }
        
        
        
      } 
      
    } ## Close "else" --> None regulators from beginning
    
    if (is.null(myPLS)) {
      
      ResultsPerTargetF.i$Y = targetData[targetF,]
      ResultsPerTargetF.i$coefficients = NULL

    } else {
      ResultsPerTargetF.i$Y = data.frame("y" = myPLS@suppLs$y, "fitted.y" = myPLS@suppLs$yPreMN,
                                      "residuals" = ropls::residuals(myPLS))
      colnames( ResultsPerTargetF.i$Y) = c('y','fitted.y','residuals')
      
      ResultsPerTargetF.i$GoodnessOfFit = c(myPLS@modelDF[,'R2Y(cum)'][myPLS@summaryDF[,'pre']],
                                         myPLS@modelDF[,'Q2(cum)'][myPLS@summaryDF[,'pre']],
                                         myPLS@summaryDF[,'RMSEE'],
                                         round(myPLS@summaryDF[,'RMSEE']/(max(myPLS@suppLs$y)-min(myPLS@suppLs$y)),6),
                                         myPLS@summaryDF[,'pre'],
                                         as.integer(length(ResultsPerTargetF.i$significantRegulators)))
      
      
    }
    
    return(ResultsPerTargetF.i)
    
}

