#########################################################################################
######           Functions to integrate omics data using MLRs                      ######
#########################################################################################


## By Maider
## 07-July-2023
## Last modified: January 2024


options(stringsAsFactors = FALSE)

library(sglfast)
library(ropls)
#'
#'\code{GetISGL} fits a MLR model with Iterative Sparse Group Lasso (ISGL) penalization
#' for all the features in the target omic to identify the conditions and potential
#' regulators that show a relevant effect on the expression of each target omic feature.
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
#' @param interactions If TRUE, the model includes interactions between regulators and condition variables. By default, TRUE.
#' @param thres  Value to determine the presence of collinearity between two regulators and consider them as group. By default, 0.7.
#' @param parallel If FALSE, MORE will be run sequentially. If TRUE, MORE will be run using parallelization with as many cores as the available ones minus one so the system is not overcharged. If the user wants to specify how many cores they want to use, they can also provide the number of cores to use in this parameter. Parallelization is only implemented for MLR with EN variable selection and PLS methods.
#' @return List containing the following elements:
#' \itemize{
#' \item ResultsPerTargetF : List with as many elements as features of the target omic in \code{\link{targetData}}. For each feature, it includes information about the feature values, considered variables, estimated coefficients,
#'                    detailed information about all regulators, and regulators identified as relevant (in MLR scenario) or significant (in PLS scenarios).
#' \item GlobalSummary : List with information about the fitted models, including model metrics, information about regulators, features of the target omic without models, regulators, master regulators and hub target features.
#' \item Arguments : List containing all the arguments used to generate the models.                
#' }
#' 
#' 
#'
#' @export
GetISGL = function(targetData,
                   regulatoryData,
                   associations =NULL,
                   omicType = 0,
                   condition = NULL,
                   clinic = NULL,
                   clinicType = NULL,
                   interactions = TRUE,
                   minVariation = 0,
                   scaleType = 'auto',
                   gr.method = 'MF',
                   thres = 0.7){
  
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
  
  #Verify that the selected grouping type it is a valid option
  if(!gr.method %in% c('MF','PCA')){stop('ERROR: The selected method for grouping is not one of MF or PCA')}
  
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
  
  ## Removing target features with too many NAs and keeping track
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
    colnames(des.mat)= sub('Group','Group_', colnames(des.mat))
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
  
  pap = c(1, 1:round(ntargetFs/100) * 100, ntargetFs)
  
  for (i in 1:ntargetFs) {
    
    targetF=AlltargetFs[i]
    ResultsPerTargetF[[i]] = vector("list", length = 5)
    names(ResultsPerTargetF[[i]]) = c("Y", "X", "coefficients", "allRegulators", "relevantRegulators")
    
    if (is.element(i, pap)) cat(paste("Fitting model for target feature", i, "out of", ntargetFs, "\n"))
    
    RetRegul = GetAllReg(targetF=targetF, associations=associations, data.omics = regulatoryData)
    RetRegul.targetF = RetRegul$Results  ## RetRegul$TableTargetF: nr reg per omic
    ## Some of these reg will be removed, because they are not in regulatoryData
    
    # RetRegul.targetF--> target feature/regulator/omic/area
    RetRegul.targetF=RetRegul.targetF[RetRegul.targetF[,"regulator"]!= "No-regulator", ,drop=FALSE] ## Remove rows with no-regulators
    
    
    ### NO INITIAL REGULATORS
    if(length(RetRegul.targetF)==0){ ## En el caso de que no hayan INICIALMENTE reguladores -> Calcular modelo con variables experiment
      
      if (is.null(condition)) {
        ResultsPerTargetF[[i]]$X = NULL
        ResultsPerTargetF[[i]]$relevantRegulators = NULL
        ResultsPerTargetF[[i]]$allRegulators = NULL
        isModel = NULL
        
        GlobalSummary$TargetFNOregu = rbind(GlobalSummary$TargetFNOregu,
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
        
        isModel = NULL
        
        ResultsPerTargetF[[i]]$X = des.mat2[,-1, drop = FALSE]
        ResultsPerTargetF[[i]]$relevantRegulators = NULL
        ResultsPerTargetF[[i]]$allRegulators = NULL
        
        GlobalSummary$TargetFNOregu = rbind(GlobalSummary$TargetFNOregu,
                                            data.frame("targetF" = targetF,
                                                       "problem" = 'Target feature had no initial regulators'))
      }
      
      
      # GlobalSummary$ReguPerTargetF  # this is initially set to 0 so no need to modify it
      
      ### WITH INITIAL REGULATORS
    } else { ## There are regulators for this target feature at the beginning
      
      ResultsPerTargetF[[i]]$allRegulators = data.frame(RetRegul.targetF, rep("Model",nrow(RetRegul.targetF)), stringsAsFactors = FALSE)
      colnames(ResultsPerTargetF[[i]]$allRegulators) = c("targetF","regulator","omic","area","filter")
      
      GlobalSummary$ReguPerTargetF[targetF, grep("-Ini", colnames(GlobalSummary$ReguPerTargetF))] = as.numeric(RetRegul$TableTargetF[-1])
      # the rest of columns remain 0
      
      ## Identify which regulators where removed because of missing values or low variation
      res = RemovedRegulators(RetRegul.targetF = ResultsPerTargetF[[i]]$allRegulators,
                              myregLV=myregLV, myregNA=myregNA, data.omics=regulatoryData)
      
      if(length(res$RegulatorMatrix)==0){ ## No regulators left after the filtering to compute the model
        
        if (is.null(condition)) {
          ResultsPerTargetF[[i]]$X = NULL
          ResultsPerTargetF[[i]]$relevantRegulators = NULL
          ResultsPerTargetF[[i]]$allRegulators = res$SummaryPerTargetF
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
          
          GlobalSummary$TargetFNOmodel = rbind(GlobalSummary$TargetFNOmodel,
                                             data.frame("targetF" = targetF,
                                                        "problem" = 'No regulators left after NA/LowVar filtering'))
          
          ResultsPerTargetF[[i]]$X = des.mat2[,-1, drop = FALSE]
          ResultsPerTargetF[[i]]$relevantRegulators = NULL
          ResultsPerTargetF[[i]]$allRegulators = res$SummaryPerTargetF
          
        }
        
      } else {  ## Regulators for the model!!
        
        ## Compute only if there is more than one regulator
        ResultsPerTargetF[[i]]$allRegulators = res$SummaryPerTargetF
        if(ncol(res$RegulatorMatrix)>1){
          #Save data needed for running ISGL
          group_index = Creategroups(data = res$RegulatorMatrix, reg.table = res$SummaryPerTargetF, correlation = thres, method = gr.method, omic.type = omicType)
          names(group_index) = colnames(res$RegulatorMatrix)
        } else{
          group_index = c(1)
          names(group_index) = colnames(res$RegulatorMatrix)
        }
        
        ## Create interactions matrix without taking into account which regulators are correlated
        
        des.mat2 = RegulatorsInteractions(interactions, reguValues = res$RegulatorMatrix, reguInfo = res$SummaryPerTargetF,
                                          des.mat, method = 'isgl')
        
        # Input data for the iterative
        train.idx = sample(nrow(des.mat2[[1]]), floor(nrow(des.mat2[[1]])*0.7))
        
        des.mat2 = lapply(des.mat2, function(x){
          data.train = x[train.idx,,drop=FALSE]
          data.validate = x[-train.idx,,drop=FALSE]
          list(data.train = data.train, data.validate = data.validate)
        })
        
        #Add condition data to scale it in the same way
        if(!is.null(des.mat)){
          cond = list(data.train = des.mat[train.idx,,drop=FALSE], data.validate = des.mat[-train.idx,,drop=FALSE])
          des.mat2[[length(des.mat2)+1]] = cond
          names(des.mat2)[length(des.mat2)] = 'condition'
        }
        
        y = t(targetData[targetF,,drop=FALSE])
        y = list(data.train = y[train.idx,,drop=FALSE], data.validate = y[-train.idx,,drop=FALSE])
        des.mat2[[length(des.mat2)+1]] = y
        names(des.mat2)[length(des.mat2)] = 'response'
        
        #Scale the variables, necessary for the models interpretation
        data = scaling.ISGL(des.mat2,des.mat, scaleType = scaleType)
        data.train = data$data.train
        data.validate = data$data.validate
        intercept = data$intercept
        rm(data);gc()
        
        #Find the groups for the interaction variables
        group_index = sapply(colnames(data.train$x), function(x) find_group(x, group_index, colnames(des.mat)))
        
        reguexp = colnames(data.train$x)
        #Call the script to obtain the coeficients
        coefs =  sglfast::isgl(data.train, data.validate, group_index, type = "linear")$beta
        
        regulatorcoef = data.frame(regulators = reguexp, coefficients = coefs)
        mycoef = regulatorcoef[which(regulatorcoef[,2] != 0),1] # selected coefficients
        
        if (length(mycoef) == 0) {
          isModel = NULL
          GlobalSummary$TargetFNOmodel = rbind(GlobalSummary$TargetFNOmodel,
                                             data.frame("targetF" = targetF, "problem" = "No relevant regulators after variable selection"))
          ## Extracting relevant regulators
          ResultsPerTargetF[[i]]$relevantRegulators = NULL
          ResultsPerTargetF[[i]]$allRegulators = data.frame(ResultsPerTargetF[[i]]$allRegulators, "Rel" = 0, stringsAsFactors = FALSE)
          
          ## Counting original regulators in the model per omic
          contando = ResultsPerTargetF[[i]]$allRegulators[which(ResultsPerTargetF[[i]]$allRegulators[,"filter"] == "Model"),]
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(regulatoryData)])
          contando[is.na(contando)] = 0
          GlobalSummary$ReguPerTargetF[targetF, grep("-Mod", colnames(GlobalSummary$ReguPerTargetF))] = contando
          
          ## Counting significant regulators per omic
          GlobalSummary$ReguPerTargetF[targetF, grep("-Sig", colnames(GlobalSummary$ReguPerTargetF))] = NA
          
        } else {
          isModel = TRUE
          ResultsPerTargetF[[i]]$X = data.train$x[,colnames(data.train$x) %in% mycoef, drop = FALSE]
          myvariables = unlist(strsplit(mycoef, ":", fixed = TRUE))
          myvariables = intersect(myvariables, ResultsPerTargetF[[i]]$allRegulators[,'regulator'])
          
          ## Remove if only group is significant
          if(length(myvariables)==0){
            GlobalSummary$TargetFNOmodel = rbind(GlobalSummary$TargetFNOmodel,
                                                 data.frame("targetF" = targetF, "problem" = "No relevant regulators after variable selection"))
          }
          
          ResultsPerTargetF[[i]]$allRegulators = data.frame(ResultsPerTargetF[[i]]$allRegulators, "Rel" = 0, stringsAsFactors = FALSE)
          ResultsPerTargetF[[i]]$allRegulators[myvariables, "Rel"] = 1
          
          ResultsPerTargetF[[i]]$relevantRegulators = ResultsPerTargetF[[i]]$allRegulators[which(ResultsPerTargetF[[i]]$allRegulators[,'Rel']==1),'regulator']
          
          contando = ResultsPerTargetF[[i]]$allRegulators[which(ResultsPerTargetF[[i]]$allRegulators[,"filter"] == "Model"),]
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(regulatoryData)])
          contando[is.na(contando)] = 0
          GlobalSummary$ReguPerTargetF[targetF, grep("-Mod", colnames(GlobalSummary$ReguPerTargetF))] = contando
          
          ## Counting relevant regulators per omic
          if (length(ResultsPerTargetF[[i]]$relevantRegulators) > 0) {
            contando = ResultsPerTargetF[[i]]$allRegulators[ResultsPerTargetF[[i]]$relevantRegulators,]
            contando = table(contando[,"omic"])
            contando = as.numeric(contando[names(regulatoryData)])
            contando[is.na(contando)] = 0
            GlobalSummary$ReguPerTargetF[targetF, grep("-Rel", colnames(GlobalSummary$ReguPerTargetF))] = contando
          }
        }
        
      }
      
      
    } ## Close "else" --> None regulators from begining
    
    if (is.null(isModel)) {
      
      ResultsPerTargetF[[i]]$Y = targetData[i,]
      ResultsPerTargetF[[i]]$coefficients = NULL
      
      GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[rownames(GlobalSummary$GoodnessOfFit) != targetF,, drop = FALSE]
      
      
    } else {
      y.fitted = data.validate$x %*% coefs
      residuals = data.validate$y - y.fitted 
      ResultsPerTargetF[[i]]$Y = data.frame("y" = data.validate$y + intercept, "fitted.y" = y.fitted + intercept, "residuals" = residuals, check.names = FALSE)
      colnames(ResultsPerTargetF[[i]]$Y) <- c("y", "fitted.y", "residuals")
      rownames(regulatorcoef) = regulatorcoef[,1]
      ResultsPerTargetF[[i]]$coefficients = regulatorcoef[mycoef,2, drop = FALSE]
      
      R.squared = round(1-(sum((data.validate$y - y.fitted)^2)/sum((data.validate$y - mean(data.validate$y))^2)),6)
      RMSE = round(sqrt(sum((data.validate$y -y.fitted)^2)/nrow(y.fitted)),6)
      NRMSE = round(RMSE/(max(data.validate$y ) - min(data.validate$y)),6)
      
      GlobalSummary$GoodnessOfFit[targetF,] = c(R.squared, RMSE, NRMSE,length(ResultsPerTargetF[[targetF]]$relevantRegulators))
      
    }
    
  }  ## At this point the loop for all target features is finished
  
  # Remove from GoodnessOfFit target features with no relevant regulators
  
  targetFsNosig = names(which(GlobalSummary$GoodnessOfFit[,4]==0))
  targetFssig = setdiff(rownames(GlobalSummary$GoodnessOfFit), targetFsNosig)
  GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[targetFssig,, drop=FALSE]
  
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
                     scaleType= scaleType, clinicType = clinicType,
                     minVariation = minVariation, associations = associations,
                     thres = thres, gr.method = gr.method,
                     targetData = targetData, regulatoryData = regulatoryData, omicType = omicType,
                     clinic = clinic, method = 'ISGL')
  
  # Create the results for the scale filter check
  
  result <- list("ResultsPerTargetF" = ResultsPerTargetF, "GlobalSummary" = GlobalSummary, "arguments" = myarguments)
  class(result) <- "MORE"
  return(result)
}


# Creating groups for isgl---------------------------------------------

find_group <-function(variable, group_index, des.mat){

  if (variable %in% des.mat){
    return(max(group_index)+1)
  }
  else if(variable %in% names(group_index)){
    return(group_index[[variable]])
  } else {
    return(group_index[[tail(strsplit(variable,':')[[1]],1)]])
  }

}

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

Creategroups = function(data, reg.table, method = 'MF' ,correlation =0.8, omic.type){


  # Take into account only regulators that could enter the model
  myreg = colnames(data)
  # Apply required method (COR: groups based on correlation, PCA: groups based on PCA)

  if (method == 'PCA'){

    # Initially extract all components
    r = min(nrow(data), ncol(data))
    if (nrow(data)<7){cross = nrow(data)-2}else{cross =7}
    respca <- try(suppressWarnings(ropls::opls(scale(data), predI = r, info.txtC='none', fig.pdfC='none',algoC='nipals',permI=0, crossvalI = cross)), silent = TRUE)

    while(class(respca)=='try-error' || length(respca@modelDF)==0 && r>0){
      respca = try(suppressWarnings( ropls::opls(data, info.txtC = 'none', fig.pdfC='none', scaleC = 'none', algoC='nipals',crossvalI = cross, permI=0, predI=r-1)), silent = TRUE)
      r = r-1
    }

    if(r == 0){
      groups = c(1:ncol(data))
    }else{
      #Compute the number of components to extract at least the 80% of the variance
      d = min(which(respca@modelDF[,'R2X(cum)']>correlation))

      if(d<r){
        #Restrict to the first d components that extract the 80% of the variance
        respca@loadingMN <- respca@loadingMN[,1:d, drop=FALSE]
      }

      #Create the groups for the regulators according to their maximal absolute value on the loadings
      groups = max.col(abs(respca@loadingMN))
    }

  }

  if (method == 'MF'){

    #calculate the correlations between the regulator pairs
    mycorrelations = data.frame(t(combn(myreg,2)),combn(myreg, 2, function(x) correlations(x,data,reg.table,omic.type)))
    mycor = mycorrelations[abs(mycorrelations[,3]) >= correlation,]

    #create the graphs for the connected regulators
    mygraph = igraph::graph_from_data_frame(mycor, directed=F)
    mycomponents = igraph::components(mygraph)
    membership<-mycomponents$membership ##save membership information
    groups = NULL

    if ( length(membership)==0){
      groups = c(1:length(myreg))

    } else {
      maxi = max(membership)
      j=1
      for (i in 1:length(myreg)){
        if(myreg[i]%in%names(membership)){
          groups[i] = membership[[myreg[i]]]
        } else{
          groups[i] = maxi+j
          j = j + 1
        }
      }

    }


  }

  return(groups)
}

# Scaling variables specific for ISGL based the scaling methodology of jlaria


scaling.ISGL = function(OmValues,des.mat, scaleType = 'auto'){
  
  
  ## Specify the scaling type
  scale <- ifelse(scaleType == 'none', FALSE, TRUE)
  center <- ifelse(scaleType == 'none', FALSE, TRUE)
  
  for (i in 1:length(OmValues)) {
    
    #Separar data.train y data.validate
    data.train = OmValues[[i]]$data.train
    data.validate = OmValues[[i]]$data.validate
    
    if (names(OmValues)[i] != 'response'){
      
      dataT.scaled = scale(data.train, scale = scale, center = center)
      
      # Centering and scaling factors
      center.at = attr(dataT.scaled, "scaled:center")
      scale.at = attr(dataT.scaled, "scaled:scale")
      
      # Scale validation data acording to scaling parameters of training data
      centered_data = sweep(data.validate, 2, center.at, FUN = "-")
      dataV.scaled = sweep(centered_data, 2, scale.at, FUN = "/")
      
      #Remove regulators lost by low variability and save as new values 
      lowVar =  names(which(apply(dataT.scaled,2,function(x) all(is.nan(x)))))
      OmValues[[i]]$data.train = dataT.scaled[,!colnames(dataT.scaled) %in% lowVar]
      OmValues[[i]]$data.validate = dataV.scaled[,!colnames(dataV.scaled) %in% lowVar]
      
    } else{
      
      if(scaleType != 'none'){
        
        intercept = mean(data.train)
        OmValues[[i]]$data.train = data.train - intercept
        OmValues[[i]]$data.validate = data.validate - intercept
        
      } else{ intercept = 0 }
    }
  }
  
  # Once variables are scaled apply the group scaling if necessary
  data.train = lapply(OmValues, function(x) x$data.train)
  data.validate = lapply(OmValues, function(x) x$data.validate)
  
  dataT.scaled = Scaling.type(data.train[!grepl('condition|response', names(data.train))],scaletype = scaleType)
  dataV.scaled = Scaling.type(data.validate[!grepl('condition|response', names(data.validate))],scaletype = scaleType)
  
  data.train = list( x = as.matrix(if (is.null(des.mat)) dataT.scaled else data.frame(data.train$condition, dataT.scaled, check.names = FALSE)), y = data.train$response )
  data.validate = list( x = as.matrix(if (is.null(des.mat)) dataV.scaled else data.frame(data.validate$condition, dataV.scaled, check.names = FALSE)), y = data.validate$response ) 
  if(!is.null(des.mat)){
    colnames(data.train$x)[1:ncol(des.mat)]=colnames(des.mat)
    colnames(data.validate$x)[1:ncol(des.mat)]=colnames(des.mat)
  }
  
  return(list('data.train'=data.train, 'data.validate'=data.validate, 'intercept'= intercept))
}



