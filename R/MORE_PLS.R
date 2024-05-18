#########################################################################################
######           Functions to integrate omics data using PLSs                      ######
#########################################################################################


## By Maider
## 20-July-2023
## Last modified: October 2023


options(stringsAsFactors = FALSE)

library(ropls)

#'
#'\code{GetPLS} fits a PLS model for all the genes in the data set to identify
#' the experimental variables and potential regulators that show a significant effect on
#' the expression of each gene.
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
#' @param epsilon Convergence threshold for coordinate descent algorithm in elasticnet. Default value, 1e-5.
#' @param alfa Significance level for variable selection in pls1and pls2 \code{\link{method}}. By default, 0.05.
#' @param interactions.reg If TRUE, the model includes interactions between regulators and experimental variables. By default, TRUE.
#' @param min.variation  For numerical regulators, it specifies the minimum change required across conditions to retain the regulator in 
#' the regression models. In the case of binary regulators, if the proportion of the most common value is equal to or inferior this value, 
#' the regulator is considered to have low variation and will be excluded from the regression models. The user has the option to set a single 
#' value to apply the same filter to all omics, provide a vector of the same length as omics if they want to specify different levels for each omics, 
#' or use 'NA' when they want to apply a minimum variation filter but are uncertain about the threshold. By default, 0.
#' @param scaletype Type of scaling to be applied. Three options:
#' \itemize{
#' \item auto : Applies the autoscaling. 
#' \item pareto : Applies the pareto scaling. \deqn{\frac{X_k}{s_k \sqrt[4]{m_b}} }
#' \item block : Applies the block scaling. \deqn{ \frac{X_k}{s_k \sqrt{m_b}} }
#' }
#' considering m_b the number of variables of the block. By default, auto.
#' @param p.method Type of resampling method to apply for the p-value calculation when pls1 or pls2 \code{\link{method}}. Two options:
#' \itemize{
#' \item jack : Applies Jack-Knife resampling technique.
#' \item perm : Applies a resampling technique in which the response variable is permuted 100 times to obtain the distribution of the coefficients and compute then their associated p-value.
#' }
#' #' By default, jack.
#' @param vip Value of VIP above which a variable can be considered significant in addition to the computed p-value in \code{\link{p.method}}. By default, 0.8.
#' @param method Model to be fitted. Two options:
#' \itemize{
#' \item pls1 : Applies a Partial Least Squares (PLS) model, one for each of the genes at \code{\link{GeneExpression}}.
#' \item pls2 : Applies a PLS model to all genes at the same time, only possible when \code{\link{associations}}= NULL.
#' }
#' #' By default, pls1.
#' @return List containing the following elements:
#' \itemize{
#' \item ResultsPerGene : List with as many elements as genes in \code{\link{GeneExpression}}. For each gene, it includes information about gene values, considered variables, estimated coefficients,
#'                    detailed information about all regulators, and regulators identified as relevant (in glm scenario) or significant (in pls scenarios).
#' \item GlobalSummary : List with information about the fitted models, including model metrics, information about regulators, genes without models, regulators, master regulators and hub genes.
#' \item Arguments : List containing all the arguments used to generate the models.                
#' }
#'
#' @export
GetPLS = function(GeneExpression,
                  data.omics,
                  associations =NULL,
                  omic.type = 0,
                  edesign = NULL,
                  clinic = NULL,
                  clinic.type = NULL,
                  center = TRUE, scale = TRUE,
                  epsilon = 0.00001,
                  alfa = 0.05, 
                  interactions.reg = TRUE,
                  min.variation = 0,
                  scaletype = 'auto',
                  p.method ='jack',
                  vip = 0.8,
                  method = 'pls1',
                  parallel = FALSE){
  
  # Converting matrix to data.frame
  GeneExpression = as.data.frame(GeneExpression)
  data.omics = lapply(data.omics, as.data.frame)
  
  ##Omic types
  if (length(omic.type) == 1) omic.type = rep(omic.type, length(data.omics))
  names(omic.type) = names(data.omics)
  
  # Creating vector for min.variation
  if (length(min.variation) == 1) min.variation=rep(min.variation,length(data.omics))
  names(min.variation)=names(data.omics)
  
  if(!is.null(clinic)){
    
    ## Clinic types
    if (length(clinic.type) == 1) {clinic.type = rep(clinic.type, ncol(clinic)); names(clinic.type) = colnames(clinic)}
    
    ##Before introducing variables in data.omics convert them to numeric type
    ## TO DO: Careful creates k-1 dummies. Is what we want?
    catvar <- which(clinic.type == 1)
    dummy_vars <- model.matrix(~ . , data = as.data.frame(clinic[,catvar ,drop=FALSE]))[,-1,drop=FALSE]
    clinic <-clinic[, -catvar,drop=FALSE]
    clinic <- cbind(clinic, dummy_vars)
    
    data.omics = c(list(clinic =  as.data.frame(t(clinic))),data.omics)
    
    #Add in associations clinic to consider all the clinical variables in all genes
    if(!is.null(associations)){associations = c(list(clinic = NULL),associations)}
    
    #Add information to omic.type and min.variation even if it is not relevant
    omic.type = c(0,omic.type)
    names(omic.type)[1] = 'clinic'
    
    min.variation = c(0,min.variation)
    names(min.variation)[1] = 'clinic'
    om= 2
    
  }else{clinic.type=NULL; om =1}
  
  
  # Not possible to apply a PLS2 when associations matrix is provided
  if(!is.null(associations)){
    if(method=='pls2'){stop('ERROR: Not possible to fit a PLS2 when associations is provided')}
  }
  
  # If associations is NULL create a list of associations NULL
  if (is.null(associations)){
    associations=vector('list',length(data.omics))
    names(associations)=names(data.omics)
  }
  
  # Checking that the number of samples per omic is equal to number of samples for gene expression and the number of samples for edesign
  for (i in 1:length(names(data.omics))){
    if(!ncol(data.omics[[i]]) == ncol(GeneExpression)) {
      stop("ERROR: Samples in data.omics must be the same as in GeneExpression and in edesign")
    }
  }
  if(!is.null(edesign)){
    if(!ncol(GeneExpression) == nrow(edesign)) {
      stop("ERROR: Samples in data.omics must be the same as in GeneExpression and in edesign")
    }
  }
  
  ##Verify that the scaletype selected is a possible selection
  
  if(!scaletype %in% c('auto','pareto','block')){stop('ERROR: The selected type of scaling is not one of auto, pareto or block')}
  if(!p.method %in% c('perm','jack')){stop('ERROR: The selected method for p value calculation is not one of perm or jack')}
  if(!method %in% c('pls1','pls2')){stop('ERROR: Not valid pls type. Select one of pls1 or pls2')}
  
  ## Checking that samples are in the same order in GeneExpressionDE, data.omics and edesign
  orderproblem<-FALSE
  if(is.null(edesign)){
    nameproblem<-!all(sapply(data.omics, function(x) length(intersect(colnames(x),colnames(GeneExpression))==ncol(GeneExpression))))
    if(nameproblem){
      cat('Warning. GeneExpression and data.omics samples have not same names. We assume that they are ordered.\n')
    }else{
      orderproblem<-!all(sapply(data.omics, function(x) identical(colnames(x),colnames(GeneExpression))))
      if(orderproblem){
        data.omics<-lapply(data.omics, function(x) x[,colnames(GeneExpression)])
      }
    }
    
  } else{
    nameproblem<-!all(c(sapply(data.omics, function(x) length(intersect(colnames(x),colnames(GeneExpression)))==ncol(GeneExpression)), length(intersect(rownames(edesign),colnames(GeneExpression)))==ncol(GeneExpression)))
    if(nameproblem){
      cat('Warning. GeneExpression, edesign and data.omics samples have not same names. We assume that they are ordered.\n')
    } else{
      orderproblem<-!all(c(sapply(data.omics, function(x) identical(colnames(x),colnames(GeneExpression))), identical(colnames(GeneExpression),rownames(edesign))))
      if(orderproblem){
        data.omics<-lapply(data.omics, function(x) x[,colnames(GeneExpression)])
        edesign<-edesign[colnames(GeneExpression), , drop=FALSE]
      }
    }
  }
  
  ## Checking if there are regulators with "_R", "_P" or "_N" or ":" and checking that there are not replicate identifiers compared to GeneExpression
  message = FALSE
  for (i in 1:length(names(data.omics))){
    
    problemas = c(rownames(data.omics[[i]])[grep("_R$", rownames(data.omics[[i]]))],
                  rownames(data.omics[[i]])[grep("_P$", rownames(data.omics[[i]]))],
                  rownames(data.omics[[i]])[grep("_N$", rownames(data.omics[[i]]))])
    
    problema = c(grep(":", rownames(data.omics[[i]]), value = TRUE))
    rownames(data.omics[[i]]) = gsub(':', '-', rownames(data.omics[[i]]))
    rownames(data.omics[[i]]) = gsub('_R$', '-R', rownames(data.omics[[i]]))
    rownames(data.omics[[i]]) = gsub('_P$', '-P', rownames(data.omics[[i]]))
    rownames(data.omics[[i]]) = gsub('_N$', '-N', rownames(data.omics[[i]]))
    
    #Change the name in the association matrix only if associations is not NULL
    if(!is.null(associations[[i]])){
      associations[[i]][[2]] = gsub(':', '-', associations[[i]][[2]])
      associations[[i]][[2]] = gsub('_R$', '-R', associations[[i]][[2]])
      associations[[i]][[2]] = gsub('_P$', '-P', associations[[i]][[2]])
      associations[[i]][[2]] = gsub('_N$', '-N', associations[[i]][[2]])
    }
    
    if(length(problemas) > 0) {
      cat("In",names(data.omics)[i], ',', problemas ,"regulators have names that may cause conflict with the algorithm by ending in _R, _P or _N", "\n")
      cat("Endings changed with -R, -P or -N, respectively", "\n")
    }
    
    if(length(problema) > 0) {
      cat("Some regulators in the omic", names(data.omics)[i],  "have names with \":\" that could cause conflict, replaced with \"-\" ", "\n")
      cat("Changed identifiers: ", problema, "\n")
    }
    
    #Checking that there are not replicate identifiers compared to GeneExpression
    
    repeated = intersect(rownames(GeneExpression), rownames(data.omics[[i]]))
    
    if(length(repeated) > 0) {
      cat(names(data.omics)[i], "and Gene Expression omics have shared identifiers in regulators:", repeated, "\n")
      #Change the name in the association matrix only if is not NULL
      if(!is.null(associations[[i]])){
        associations[[i]][[2]][associations[[i]][[2]]%in%repeated] = paste(names(data.omics)[i],'-', associations[[i]][[2]][associations[[i]][[2]]%in%repeated], sep='')
      }
      #Change the name in data.omics
      rownames(data.omics[[i]])[rownames(data.omics[[i]])%in%repeated] = paste(names(data.omics)[i],'-',  rownames(data.omics[[i]])[rownames(data.omics[[i]])%in%repeated], sep='')
    }
  }
  
  
  ##Checking that there are no replicates in the identifiers and changing identifiers in case of need
  
  if(length(names(data.omics))>1){
    for (i in 1:(length(names(data.omics))-1)){
      for(j in (i+1):(length(names(data.omics)))){
        repeated = intersect(rownames(data.omics[[i]]), rownames(data.omics[[j]]))
        
        if(length(repeated) > 0) {
          cat(names(data.omics)[i], "and", names(data.omics)[j], "omics have shared identifiers in regulators:", repeated, "\n")
          #Change the name in the association matrix only if is not NULL
          if(!is.null(associations[[i]])){
            associations[[i]][[2]][associations[[i]][[2]]%in%repeated] =  paste(names(data.omics)[i],'-', associations[[i]][[2]][associations[[i]][[2]]%in%repeated], sep='')
          }
          if(!is.null(associations[[j]])){
            associations[[j]][[2]][associations[[j]][[2]]%in%repeated] =  paste(names(data.omics)[j],'-', associations[[j]][[2]][associations[[j]][[2]]%in%repeated], sep='')
          }
          #Change the name in data.omics
          rownames(data.omics[[i]])[rownames(data.omics[[i]])%in%repeated] =  paste(names(data.omics)[i],'-', rownames(data.omics[[i]])[rownames(data.omics[[i]])%in%repeated],sep = '')
          rownames(data.omics[[j]])[rownames(data.omics[[j]])%in%repeated] =  paste(names(data.omics)[j],'-', rownames(data.omics[[j]])[rownames(data.omics[[j]])%in%repeated],sep = '')
          
        }
      }
    }
    rm(repeated)
  }
  
  #Checking there are not -Inf/Inf values and eliminate genes/regulator that contain them
  infproblemgene<-is.infinite(rowSums(GeneExpression))
  infproblemreg<-lapply(data.omics[om:length(data.omics)], function(x) is.infinite(rowSums(x)))
  if(any(infproblemgene)){
    genesInf<-rownames(GeneExpression)[infproblemgene]
    GeneExpression<-GeneExpression[!infproblemgene,]
  }else{genesInf <-NULL}
  for (i in 1:(length(names(data.omics))-(om-1))){
    if(any(infproblemreg[[i]])){
      cat(rownames(data.omics[[i + (om-1)]])[infproblemreg[[i]]], 'regulators of the omic', names(data.omics)[i +(om-1)] ,'have been deleted due to -Inf/Inf values. \n')
      data.omics[[i + (om-1)]]<-data.omics[[i + (om-1)]][!infproblemreg[[i]],]
    }
  }
  
  ## Removing genes with more than 20% of NAs and keeping track
  min.obs = floor( ncol(GeneExpression)*0.8)
  genesNotNA = apply(GeneExpression, 1, function (x) sum(!is.na(x)))
  genesNotNA = names(which(genesNotNA >= min.obs))
  genesNA = setdiff(rownames(GeneExpression), genesNotNA)
  GeneExpression = GeneExpression[genesNotNA,]
  
  ## Removing genes with no regulators only if associations does not have an associations = NULL in any omic
  genesNOreg = NULL
  genesNOreg = lapply(associations, function(x) if(!is.null(x)) {setdiff( rownames(GeneExpression),x[,1])})
  genesNOreg = Reduce(intersect, genesNOreg)
  GeneExpression = GeneExpression[!(rownames(GeneExpression) %in% genesNOreg),]
  if (length(genesNOreg) > 0){
    cat(length(genesNOreg), "genes had no initial regulators. Models will be computed for", nrow(GeneExpression), 'genes.\n')
  }
  
  ## Removing constant genes
  constantGenes = apply(GeneExpression, 1, sd, na.rm = TRUE)
  notConstant = names(constantGenes)[constantGenes > 0]
  constantGenes = names(constantGenes)[constantGenes == 0]
  GeneExpression = GeneExpression[notConstant,]
  
  Allgenes=rownames(GeneExpression)
  nGenes = length(Allgenes)
  
  # Experimental groups
  if (is.null(edesign)) {
    cat("No experimental covariates were provided.\n")
    Group = 1:ncol(GeneExpression)
    names(Group) = colnames(GeneExpression)
    des.mat = NULL
  } else {
    Group = apply(edesign, 1, paste, collapse = "_")
    des.mat = model.matrix(~0+., data = as.data.frame(Group))
    #Change the name to avoid conflicts with RegulationPerCondition
    colnames(des.mat) = sub('Group','Group_',colnames(des.mat))
  }
  
  ## Remove regulators with NA
  cat("Removing regulators with missing values...\n")
  
  myregNA = lapply(data.omics, rownames)
  data.omics = lapply(data.omics, na.omit)
  myregNA = lapply(1:length(data.omics), function (i) setdiff(myregNA[[i]], rownames(data.omics[[i]])))
  names(myregNA)=names(data.omics)
  
  cat("Number of regulators with missing values:\n")
  print(sapply(myregNA, length))
  cat("\n")
  
  
  ## Remove regulators with Low Variability
  cat("Removing regulators with low variation...\n")
  
  tmp = LowVariationRegu(min.variation, data.omics, Group, associations, Allgenes, omic.type, clinic.type)
  data.omics = tmp[["data.omics"]]
  associations = tmp[["associations"]]
  myregLV = tmp[["myregLV"]]
  rm("tmp"); gc()
  
  if(all(sapply(data.omics, function(x)nrow(x)==0))) stop("ERROR: No regulators left after LowVariation filter. Consider being less restrictive.")

  ### Results objects
  
  ## Global summary for all genes
  GlobalSummary = vector("list", length = 6)
  names(GlobalSummary) = c("GoodnessOfFit", "ReguPerGene", "GenesNOmodel", "GenesNOregulators", "GlobalRegulators", "HubGenes")
  
  GlobalSummary$GenesNOmodel = NULL
  if (length(genesNA) > 0) {
    GlobalSummary$GenesNOmodel = data.frame("gene" = genesNA,
                                            "problem" = rep("Too many missing values", length(genesNA)))
  }
  if (length(constantGenes) > 0) {
    GlobalSummary$GenesNOmodel = rbind(GlobalSummary$GenesNOmodel,
                                       data.frame("gene" = constantGenes,
                                                  "problem" = rep("Response values are constant", length(constantGenes))))
  }
  if (length(genesInf) > 0){
    GlobalSummary$GenesNOmodel = rbind(GlobalSummary$GenesNOmodel,
                                       data.frame("gene" = genesInf,
                                                  "problem" = rep("-Inf/Inf values", length(genesInf))))
    
  }
  
  GlobalSummary$GenesNOregulators = NULL
  if (length(genesNOreg) > 0){
    GlobalSummary$GenesNoregulators = data.frame("gene" = genesNOreg, "problem" = rep("Gene had no initial regulators", length(genesNOreg)))
  }
  
  GlobalSummary$GoodnessOfFit = matrix(NA, ncol = 6, nrow = nGenes)
  rownames(GlobalSummary$GoodnessOfFit) = Allgenes
  colnames(GlobalSummary$GoodnessOfFit) = c( "RsquaredY", "Qsquared","RMSE","CV(RMSE)","ncomp","sigReg")
  
  GlobalSummary$ReguPerGene = matrix(0, ncol = 3*length(data.omics), nrow = nGenes)
  rownames(GlobalSummary$ReguPerGene) = Allgenes
  colnames(GlobalSummary$ReguPerGene) = c(paste(names(data.omics), "Ini", sep = "-"),
                                          paste(names(data.omics), "Mod", sep = "-"),
                                          paste(names(data.omics), "Sig", sep = "-"))
  rm(infproblemgene);rm(infproblemreg);rm(problema);rm(problemas);rm(orderproblem);rm(nameproblem);rm(genesNA);rm(genesInf);rm(genesNOreg);rm(genesNotNA);rm(constantGenes);rm(notConstant);gc()
  ## Specific results for each gene
  ResultsPerGene=vector("list", length=length(Allgenes))
  names(ResultsPerGene) = Allgenes
 
  if(method=='pls1'){
    ### Computing model for each gene
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
      }
    }

    ResultsPerGene <- furrr::future_map(1:nGenes,
                                        ~ResultsPerGene.i(Allgenes[.],GlobalSummary,data.omics,associations,GeneExpression,omic.type,
                                                          edesign, des.mat,myregLV,myregNA,epsilon,scale,center,scaletype,
                                                          interactions.reg,p.method,vip,alfa),
                                        .progress = TRUE )
    future::plan("sequential")
    names(ResultsPerGene)<-Allgenes
    
    globalValues<-c('GoodnessOfFit','ReguPerGene','GenesNOmodel')
    
    for(i in 1:length(ResultsPerGene)){
      
      for (gValue in globalValues){
        if (exists(gValue, ResultsPerGene[[i]]) && ! is.null(ResultsPerGene[[i]][[gValue]])){
          if(! is.data.frame(ResultsPerGene[[i]][[gValue]])) {
            GlobalSummary[[gValue]][Allgenes[i],]<-ResultsPerGene[[i]][[gValue]]
          } else{
            GlobalSummary[[gValue]]<-rbind(GlobalSummary[[gValue]], ResultsPerGene[[i]][[gValue]])
          }
        }
      }
      ResultsPerGene[[i]] <- ResultsPerGene[[i]][-grep(paste0(globalValues, collapse = "|"), names(ResultsPerGene[[i]]))]
    }
    
    
  } else{
    
    cat("Fitting model ...\n")
    
    #When associations = NULL all genes have the same potential regulators
    gene = Allgenes[1]
    RetRegul = GetAllReg(gene=gene, associations=associations, data.omics = data.omics)
    RetRegul.gene = RetRegul$Results  ## RetRegul$TableGene: nr reg per omic
    ## Some of these reg will be removed, because they are not in data.omics
    
    #Add information about initial potential regulators
    GlobalSummary$ReguPerGene[1:nrow(GlobalSummary$ReguPerGene), grep("-Ini", colnames(GlobalSummary$ReguPerGene))]=rep(as.numeric(RetRegul$TableGene[-1]),each=length(Allgenes))
    
    # RetRegul.gene--> gene/regulator/omic/area
    RetRegul.gene=RetRegul.gene[RetRegul.gene[,"regulator"]!= "No-regulator", ,drop=FALSE] ## Remove rows with no-regulators
    
    allRegulators = data.frame(RetRegul.gene, rep("Model",nrow(RetRegul.gene)), stringsAsFactors = FALSE)
    colnames(allRegulators) = c("gene","regulator","omic","area","filter")
    
    ## Identify which regulators where removed because of missing values or low variation
    res = RemovedRegulators(RetRegul.gene = allRegulators,
                            myregLV=myregLV, myregNA=myregNA, data.omics=data.omics)
    rm(myregLV);rm(myregNA);rm(RetRegul.gene);rm(RetRegul);gc()
    ## Create the interactions between regulators and edesign 
    
    des.mat2 = RegulatorsInteractionsPLS2(interactions.reg, reguValues = res$RegulatorMatrix,
                                          des.mat)
    
    #Scale the variables, indispensable for elasticnet application
    des.mat2 = scale(des.mat2,scale=scale,center=center)
    
    ##Scale if needed to block scaling or pareto scaling
    if (scaletype!='auto'){
      ## Make the groups of omics
      regupero = lapply(unique(res$SummaryPerGene[,'omic']), function(x) rownames(res$SummaryPerGene)[res$SummaryPerGene[,'omic'] == x & res$SummaryPerGene[,'filter'] == "Model"])
      names(regupero) = unique(res$SummaryPerGene[,'omic'])
      #It does not work in case of really huge amount of data
      regu = try(suppressWarnings( lapply(regupero, function(x) colnames(des.mat2[,grep(paste(x, collapse = "|"), colnames(des.mat2)),drop=FALSE]))),silent = TRUE)
      if(class(regu)=='try-error'){
        #Add the ones related to the interactions
        regupero = filter_columns_by_regexp(regupero, des.mat2,res)
      }else{regupero = regu}
      res$RegulatorMatrix = Scaling.type(des.mat2, regupero, scaletype)
      
      #Use them jointly
      des.mat2 = data.frame(des.mat2[,setdiff(colnames(des.mat2),colnames(res$RegulatorMatrix)),drop=FALSE], res$RegulatorMatrix,check.names = FALSE)
      rm(regupero);rm(res); gc()
    } 
    
    # Removing predictors with constant values
    sdNo0 = apply(des.mat2, 2, sd)
    sdNo0 = names(sdNo0)[sdNo0 > 0]
    quito = setdiff(colnames(des.mat2), sdNo0)
    des.mat2 = des.mat2[,sdNo0]
    sdNo0 = unique(gsub(".*[:?](.*)$", "\\1", sdNo0))
    quito = unique(gsub(".*[:?](.*)$", "\\1", quito))
    quito = setdiff(quito, sdNo0)
    allRegulators[quito,"filter"] = "Constant"
    
    Y = scale(t(GeneExpression), center = center, scale = scale)
    
    ## Computing GLM model
    if (nrow(des.mat2)<7){cross = nrow(des.mat)-2}else{cross =7}
    myPLS = try(suppressWarnings( ropls::opls(des.mat2, Y, info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0)),silent = TRUE)
    
    if(class(myPLS)=='try-error' || length(myPLS@modelDF)==0){
      myPLS = try(suppressWarnings( ropls::opls(des.mat2, Y, info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0, predI=1)), silent = TRUE)
    }
    
    if (class(myPLS)=='try-error' || length(myPLS@modelDF)==0){
      
      myPLS = NULL
      
      GlobalSummary$GenesNOmodel = rbind(GlobalSummary$GenesNOmodel,
                                         data.frame("gene" = Allgenes, "problem" = "No significant components on PLS2"))
      
      ## Extracting significant regulators
      for (i in 1:length(Allgenes)) {
        gene = Allgenes[i]
        ResultsPerGene[[i]]$significantRegulators = NULL
        ResultsPerGene[[i]]$allRegulators = data.frame(ResultsPerGene[[i]]$allRegulators, "Sig" = NA, stringsAsFactors = FALSE)
        
        ResultsPerGene[[i]]$allRegulators = allRegulators
        ResultsPerGene[[i]]$allRegulators[,'gene']=rep(gene,nrow(ResultsPerGene[[i]]$allRegulators))
        
        ## Counting original regulators in the model per omic
        contando = ResultsPerGene[[i]]$allRegulators[which(ResultsPerGene[[i]]$allRegulators[,"filter"] == "Model"),]
        contando = table(contando[,"omic"])
        contando = as.numeric(contando[names(data.omics)])
        contando[is.na(contando)] = 0
        GlobalSummary$ReguPerGene[gene, grep("-Mod", colnames(GlobalSummary$ReguPerGene))] = contando
        
        ## Counting significant regulators per omic
        GlobalSummary$ReguPerGene[gene, grep("-Sig", colnames(GlobalSummary$ReguPerGene))] = NA
        
        ResultsPerGene[[i]]$Y = GeneExpression[i,]
        ResultsPerGene[[i]]$coefficients = NULL
        
        GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[rownames(GlobalSummary$GoodnessOfFit) != gene,]
        
      }
      
      
    }else{
      
      myPLS = suppressInnecPLSdata(myPLS)
      if (p.method == 'jack'){
        pval = p.valuejack.pls2(myPLS, des.mat2, Y, alfa)
      } else {
        pval = p.coef.pls2(myPLS, 100, des.mat2, Y)
      }
      
      #Aunque lo hayamos tratado todo junto ahora separamos la respuesta por genes
      
      for (i in 1:ncol(Y)) {
        
        #Tratar como significativas tan solo las que cumplan ambas condiciones
        sigvariables = intersect(names(myPLS@vipVn[which(myPLS@vipVn>vip)]), rownames(pval[,i,drop=FALSE])[which(pval[,i,drop=FALSE]<alfa)])
        ResultsPerGene[[i]]$coefficients = data.frame('coefficient' = myPLS@coefficientMN[sigvariables,i], 'pvalue' = pval[sigvariables,i,drop=FALSE])

        ## Extracting significant regulators and recovering correlated regulators
        myvariables = unlist(strsplit(sigvariables, ":", fixed = TRUE))
        myvariables = intersect(myvariables, allRegulators[,'regulator'])
        
        ResultsPerGene[[i]]$allRegulators = allRegulators
        
        gene = Allgenes[i]
        ResultsPerGene[[i]]$allRegulators[,'gene']=rep(gene,nrow(ResultsPerGene[[i]]$allRegulators))
        
        ResultsPerGene[[i]]$allRegulators = data.frame(ResultsPerGene[[i]]$allRegulators, "Sig" = 0, stringsAsFactors = FALSE)
        ResultsPerGene[[i]]$allRegulators[myvariables, "Sig"] = 1
        
        ResultsPerGene[[i]]$significantRegulators = myvariables
        
        ## Counting original regulators in the model per omic              
        contando = ResultsPerGene[[i]]$allRegulators
        quitar = which(contando[,"filter"] == "MissingValue")
        if (length(quitar) > 0) contando = contando[-quitar,]
        quitar = which(contando[,"filter"] == "LowVariation")
        if (length(quitar) > 0) contando = contando[-quitar,]
        
        contando = table(contando[,"omic"])
        contando = as.numeric(contando[names(data.omics)])
        contando[is.na(contando)] = 0
        GlobalSummary$ReguPerGene[gene, grep("-Mod", colnames(GlobalSummary$ReguPerGene))] = contando
        
        ## Counting significant regulators per omic
        if (length(ResultsPerGene[[i]]$significantRegulators) > 0) {
          contando = ResultsPerGene[[i]]$allRegulators[ResultsPerGene[[i]]$significantRegulators,]
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(data.omics)])
          contando[is.na(contando)] = 0
          GlobalSummary$ReguPerGene[gene, grep("-Sig", colnames(GlobalSummary$ReguPerGene))] = contando
        }
        
        ResultsPerGene[[i]]$Y = data.frame("y" = myPLS@suppLs$y[,i,drop=FALSE], "fitted.y" = myPLS@suppLs$yPreMN[,i,drop=FALSE],
                                           "residuals" = myPLS@suppLs$y[,i,drop=FALSE]-myPLS@suppLs$yPreMN[,i,drop=FALSE])
        colnames( ResultsPerGene[[i]]$Y) = c('y','fitted.y','residuals')
        
        
        GlobalSummary$GoodnessOfFit[gene,] = c(myPLS@modelDF[,'R2Y(cum)'][myPLS@summaryDF[,'pre']],
                                               myPLS@modelDF[,'Q2(cum)'][myPLS@summaryDF[,'pre']],
                                               myPLS@summaryDF[,'RMSEE'],
                                               round(abs(myPLS@summaryDF[,'RMSEE']/mean(myPLS@suppLs$y)),6),
                                               myPLS@summaryDF[,'pre'],
                                               as.integer(length(ResultsPerGene[[i]]$significantRegulators)))
      }
      
      
    }
    
  } 

  genesNosig = names(which(GlobalSummary$GoodnessOfFit[,6]==0))
  genessig = setdiff(rownames(GlobalSummary$GoodnessOfFit), genesNosig)
  GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[genessig,]
  
  #Calculate GlobalRegulators
  m_sig_reg<-lapply(ResultsPerGene, function(x) x$significantRegulators)
  m_sig_reg <- unlist(m_sig_reg)
  msig_vector <- table(m_sig_reg)
  #Calculate third quantile
  q3<-quantile(msig_vector,0.75)
  if(length(msig_vector[msig_vector>q3])<10){
    GlobalSummary$GlobalRegulators = intersect(names(msig_vector[rev(tail(order(msig_vector),10))]), names(msig_vector[msig_vector>10]) )
  } else{
    GlobalSummary$GlobalRegulators = intersect(names(msig_vector[msig_vector>q3]), names(msig_vector[msig_vector>10]) ) 
  }
  
  #Calculate HubGenes
  significant_regulators<-GlobalSummary$ReguPerGene[,c(grep('-Sig$',colnames(GlobalSummary$ReguPerGene))),drop=FALSE]
  s_sig_reg<-apply(significant_regulators, 1, sum)
  #Calculate third quantile
  q3<-quantile(s_sig_reg,0.75)
  if(length(s_sig_reg[s_sig_reg>q3])<10){
    GlobalSummary$HubGenes = intersect(names(s_sig_reg[rev(tail(order(s_sig_reg),10))]), names(s_sig_reg[s_sig_reg>10]) )
  } else{
    GlobalSummary$HubGenes = intersect(names(s_sig_reg[s_sig_reg>q3]), names(s_sig_reg[s_sig_reg>10]))
  }
  
  if(all(sapply(associations,is.null))) {associations = NULL}
  
  myarguments = list(edesign = edesign, finaldesign = des.mat, groups = Group, alfa = alfa, 
                     center = center, scale = scale, clinic.type = clinic.type,
                     min.variation = min.variation, associations = associations,
                     epsilon = epsilon, vip = vip,
                     GeneExpression = GeneExpression, dataOmics = data.omics, omic.type = omic.type,
                     clinic = clinic, clinic.type = clinic.type, scaletype =scaletype, p.method=p.method, method =method)

  # Create the results for the scale filter check
  
  result <- list("ResultsPerGene" = ResultsPerGene, "GlobalSummary" = GlobalSummary, "arguments" = myarguments) 
  class(result) <- "MORE"
  return(result)
}



ResultsPerGene.i<-function(gene,GlobalSummary,data.omics,associations,GeneExpression,omic.type,
                           edesign, des.mat,myregLV,myregNA,epsilon,scale,center,scaletype,
                           interactions.reg,p.method,vip,alfa){
  
  
    ResultsPerGene.i = vector("list", length = 8)
    names(ResultsPerGene.i) = c("Y", "X", "coefficients", "allRegulators", "significantRegulators", "GoodnessOfFit", "ReguPerGene","GenesNOmodel")
    
    #Initialize gobal summary values
    ResultsPerGene.i$ReguPerGene <- GlobalSummary$ReguPerGene[gene, , drop=FALSE]
    
    RetRegul = GetAllReg(gene=gene, associations=associations, data.omics = data.omics)
    RetRegul.gene = RetRegul$Results  ## RetRegul$TableGene: nr reg per omic
    ## Some of these reg will be removed, because they are not in data.omics
    
    
    # RetRegul.gene--> gene/regulator/omic/area
    RetRegul.gene=RetRegul.gene[RetRegul.gene[,"regulator"]!= "No-regulator", ,drop=FALSE] ## Remove rows with no-regulators
    
    
    ### NO INITIAL REGULATORS
    if(length(RetRegul.gene)==0){ ## En el caso de que no hayan INICIALMENTE reguladores -> Calcular modelo con variables experimentales
      
      if (is.null(edesign)) {
        ResultsPerGene.i$X = NULL
        ResultsPerGene.i$significantRegulators = NULL
        ResultsPerGene.i$allRegulators = NULL
        myPLS = NULL
        
      } else {
        des.mat2 = cbind(t(GeneExpression[gene,]), des.mat)
        colnames(des.mat2)[1] = "response"
        des.mat2 = na.omit(des.mat2)
        
        # Removing predictors with constant values
        sdNo0 = apply(des.mat2, 2, sd)
        sdNo0 = names(sdNo0)[sdNo0 > 0]
        des.mat2 = des.mat2[,sdNo0]
        
        myPLS = NULL
        
        ResultsPerGene.i$X = des.mat2[,-1, drop = FALSE]
        ResultsPerGene.i$significantRegulators = NULL
        ResultsPerGene.i$allRegulators = NULL
      }
      
      # GlobalSummary$ReguPerGene  # this is initially set to 0 so no need to modify it
      
      ### WITH INITIAL REGULATORS
    } else { ## There are regulators for this gene at the beginning
      
      ResultsPerGene.i$allRegulators = data.frame(RetRegul.gene, rep("Model",nrow(RetRegul.gene)), stringsAsFactors = FALSE)
      colnames(ResultsPerGene.i$allRegulators) = c("gene","regulator","omic","area","filter")
      
      ResultsPerGene.i$ReguPerGene[1, grep("-Ini", colnames(ResultsPerGene.i$ReguPerGene))] = as.numeric(RetRegul$TableGene[-1])
      # the rest of columns remain 0
      
      ## Identify which regulators where removed because of missing values or low variation
      res = RemovedRegulators(RetRegul.gene = ResultsPerGene.i$allRegulators,
                              myregLV=myregLV, myregNA=myregNA, data.omics=data.omics)
      
      if(length(res$RegulatorMatrix)==0){ ## No regulators left after the filtering to compute the model
        
        if (is.null(edesign)) {
          ResultsPerGene.i$X = NULL
          ResultsPerGene.i$significantRegulators = NULL
          ResultsPerGene.i$allRegulators = res$SummaryPerGene
          myPLS = NULL
          
        } else {
          des.mat2 = cbind(t(GeneExpression[gene,]), des.mat)
          colnames(des.mat2)[1] = "response"
          des.mat2 = na.omit(des.mat2)
          
          # Removing predictors with constant values
          sdNo0 = apply(des.mat2, 2, sd)
          sdNo0 = names(sdNo0)[sdNo0 > 0]
          des.mat2 = des.mat2[,sdNo0]
          
          myPLS = NULL
          
          ResultsPerGene.i$GenesNOmodel = rbind(ResultsPerGene.i$GenesNOmodel,
                                             data.frame("gene" = gene,
                                                        "problem" = 'No regulators left after NA/LowVar filtering'))
          ResultsPerGene.i$X = des.mat2[,-1, drop = FALSE]
          ResultsPerGene.i$significantRegulators = NULL
          ResultsPerGene.i$allRegulators = res$SummaryPerGene
        }
        
      } else {  ## Regulators for the model!!
        ResultsPerGene.i$allRegulators = res$SummaryPerGene
        
        ## Create the interactions between regulators and edesign 
        des.mat2 = RegulatorsInteractions(interactions.reg, reguValues = res$RegulatorMatrix,
                                          des.mat, GeneExpression, gene)
        
        #Scale the variables, indispensable for elasticnet application
        des.mat2 = data.frame(des.mat2[,1,drop=FALSE], scale(des.mat2[,-1,drop=FALSE],scale=scale,center=center),check.names = FALSE)
        
        ##Scale if needed to block scaling or pareto scaling
        if (scaletype!='auto'){
          ## Make the groups of omics
          regupero = lapply(unique(res$SummaryPerGene[,'omic']), function(x) rownames(res$SummaryPerGene)[res$SummaryPerGene[,'omic'] == x & res$SummaryPerGene[,'filter'] == "Model"])
          names(regupero) = unique(res$SummaryPerGene[,'omic'])
          #Remove empty omics
          regupero = regupero[sapply(regupero, function(x) length(x) > 0)]
          #It does not work in case of really huge amount of data
          regupero1 = try(suppressWarnings( lapply(regupero, function(x) colnames(des.mat2[,grep(paste(x, collapse = "|"), colnames(des.mat2)),drop=FALSE]))),silent = TRUE)
          if(class(regupero1)=='try-error'){
            #Add the ones related to the interactions
            regupero = filter_columns_by_regexp(regupero, des.mat2,res)
          }else{regupero = regupero1}
          res$RegulatorMatrix = Scaling.type(des.mat2[,-1,drop=FALSE], regupero, scaletype)
          
          #Use them jointly
          des.mat2 = data.frame(des.mat2[,1,drop=FALSE], scale(des.mat,scale=scale,center=center), res$RegulatorMatrix,check.names = FALSE)
          rm(regupero);rm(res)
        } 
        
        # Removing predictors with constant values
        sdNo0 = apply(des.mat2, 2, sd)
        sdNo0 = names(sdNo0)[sdNo0 > 0]
        quito = setdiff(colnames(des.mat2), sdNo0)
        des.mat2 = des.mat2[,sdNo0]
        sdNo0 = unique(gsub(".*[:?](.*)$", "\\1", sdNo0))
        quito = unique(gsub(".*[:?](.*)$", "\\1", quito))
        quito = setdiff(quito, sdNo0)
        ResultsPerGene.i$allRegulators[quito,"filter"] = "Constant"
        
        #Save X matrix
        
        ResultsPerGene.i$X = des.mat2[,-1, drop = FALSE]
        
        ## Computing PLS model
        if (nrow(des.mat2)<7){cross = nrow(des.mat)-2}else{cross =7}
        myPLS = try(suppressWarnings(ropls::opls(des.mat2[,-1], scale(des.mat2[,1], center = center,scale=scale), info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0)),silent = TRUE)
        
        if(class(myPLS)=='try-error' || length(myPLS@modelDF)==0){
          myPLS = try(suppressWarnings(ropls::opls(des.mat2[,-1], scale(des.mat2[,1], center =center, scale = scale), info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0, predI=1)), silent = TRUE)
        }
        
        if (class(myPLS)=='try-error' || length(myPLS@modelDF)==0){
          
          myPLS = NULL
          
          ResultsPerGene.i$GenesNOmodel = rbind(ResultsPerGene.i$GenesNOmodel,
                                             data.frame("gene" = gene, "problem" = "No significant components on PLS"))
          
          ## Extracting significant regulators
          ResultsPerGene.i$significantRegulators = NULL
          ResultsPerGene.i$allRegulators = data.frame(ResultsPerGene.i$allRegulators, "Sig" = NA, stringsAsFactors = FALSE)
          
          ## Counting original regulators in the model per omic
          contando = ResultsPerGene.i$allRegulators[which(ResultsPerGene.i$allRegulators[,"filter"] == "Model"),]
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(data.omics)])
          contando[is.na(contando)] = 0
          ResultsPerGene.i$ReguPerGene[gene, grep("-Mod", colnames(ResultsPerGene.i$ReguPerGene))] = contando
          
          ## Counting significant regulators per omic
          ResultsPerGene.i$ReguPerGene[gene, grep("-Sig", colnames(ResultsPerGene.i$ReguPerGene))] = NA
          
        }else{
          
          if (p.method == 'jack'){
            pval = p.valuejack(myPLS, des.mat2, alfa)
          } else {
            pval = p.coef(myPLS, 100, des.mat2)
          }
          
          #Tratar como significativas tan solo las que cumplan ambas condiciones
          sigvariables = intersect(names(myPLS@vipVn[which(myPLS@vipVn>vip)]), rownames(pval)[which(pval<alfa)])
          ResultsPerGene.i$coefficients = data.frame('coefficient' = myPLS@coefficientMN[sigvariables,], 'pvalue' = pval[sigvariables,,drop=FALSE])
          
          ## Extracting significant regulators and recovering correlated regulators
          myvariables = unlist(strsplit(sigvariables, ":", fixed = TRUE))
          myvariables = intersect(myvariables, rownames(ResultsPerGene.i$allRegulators))
          
          ResultsPerGene.i$allRegulators = data.frame(ResultsPerGene.i$allRegulators, "Sig" = 0, stringsAsFactors = FALSE)
          ResultsPerGene.i$allRegulators[myvariables, "Sig"] = 1
          
          ResultsPerGene.i$significantRegulators = myvariables
          
          ## Counting original regulators in the model per omic              
          contando = ResultsPerGene.i$allRegulators
          quitar = which(contando[,"filter"] == "MissingValue")
          if (length(quitar) > 0) contando = contando[-quitar,]
          quitar = which(contando[,"filter"] == "LowVariation")
          if (length(quitar) > 0) contando = contando[-quitar,]
          
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(data.omics)])
          contando[is.na(contando)] = 0
          ResultsPerGene.i$ReguPerGene[gene, grep("-Mod", colnames(ResultsPerGene.i$ReguPerGene))] = contando
          
          ## Counting significant regulators per omic
          if (length(ResultsPerGene.i$significantRegulators) > 0) {
            contando = ResultsPerGene.i$allRegulators[ResultsPerGene.i$significantRegulators,]
            contando = table(contando[,"omic"])
            contando = as.numeric(contando[names(data.omics)])
            contando[is.na(contando)] = 0
            ResultsPerGene.i$ReguPerGene[gene, grep("-Sig", colnames(ResultsPerGene.i$ReguPerGene))] = contando
          }
          
          
        }
        
        
        
      } ## Close "else" --> None regulators from begining
      
      if (is.null(myPLS)) {
        
        ResultsPerGene.i$Y = GeneExpression[gene,]
        ResultsPerGene.i$coefficients = NULL
        
        ResultsPerGene.i$GoodnessOfFit = ResultsPerGene.i$GoodnessOfFit[rownames(ResultsPerGene.i$GoodnessOfFit) != gene,,drop=FALSE]
        
        
      } else {
        ResultsPerGene.i$Y = data.frame("y" = myPLS@suppLs$y, "fitted.y" = myPLS@suppLs$yPreMN,
                                           "residuals" = ropls::residuals(myPLS))
        colnames( ResultsPerGene.i$Y) = c('y','fitted.y','residuals')
        
        ResultsPerGene.i$GoodnessOfFit = c(myPLS@modelDF[,'R2Y(cum)'][myPLS@summaryDF[,'pre']],
                                               myPLS@modelDF[,'Q2(cum)'][myPLS@summaryDF[,'pre']],
                                               myPLS@summaryDF[,'RMSEE'],
                                               round(abs(myPLS@summaryDF[,'RMSEE']/mean(myPLS@suppLs$y)),6),
                                               myPLS@summaryDF[,'pre'],
                                               as.integer(length(ResultsPerGene.i$significantRegulators)))
        
        
      }
      
    }
    
    return(ResultsPerGene.i)
    
}

