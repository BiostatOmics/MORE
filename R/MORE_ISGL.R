#########################################################################################
######           Functions to integrate omics data using GLMs                      ######
#########################################################################################


## By Maider
## 07-July-2023
## Last modified: January 2024


options(stringsAsFactors = FALSE)

library(sglfast)
library(ropls)
#'
#'\code{GetISGL} fits a GLM model with Iterative Sparse Group Lasso (ISGL) penalization
#' for all the genes in the data set to identify the experimental variables and potential
#' regulators that show a relevant effect on the expression of each gene.
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
#' @param interactions.reg If TRUE, the model includes interactions between regulators and experimental variables. By default, TRUE.
#' @param min.variation  For numerical regulators, it specifies the minimum change required across conditions to retain the regulator in
#' the regression models. In the case of binary regulators, if the proportion of the most common value is equal to or inferior this value,
#' the regulator is considered to have low variation and will be excluded from the regression models. The user has the option to set a single
#' value to apply the same filter to all omics, provide a vector of the same length as omics if they want to specify different levels for each omics,
#' or use 'NA' when they want to apply a minimum variation filter but are uncertain about the threshold. By default, 0.
#' @param gr.method Methodology to apply to create the gorups. By default, cor.
#' @param thres Threshold for the correlation when using gr.method ='cor' or threshold for the percentage of variability to explain when gr.method ='pca'. By default, 0.7.
#' @return List containing the following elements:
#' \itemize{
#' \item ResultsPerGene : List with as many elements as genes in \code{\link{GeneExpression}}. For each gene, it includes information about gene values, considered variables, estimated coefficients,
#'                    detailed information about all regulators, and regulators identified as relevant (in glm scenario) or significant (in pls scenarios).
#' \item GlobalSummary : List with information about the fitted models, including model metrics, information about regulators, genes without models, regulators, master regulators and hub genes.
#' \item Arguments : List containing all the arguments used to generate the models.
#' }
#'
#' @export
GetISGL = function(GeneExpression,
                  data.omics,
                  associations =NULL,
                  omic.type = 0,
                  edesign = NULL,
                  clinic = NULL,
                  clinic.type = NULL,
                  center = TRUE, scale = TRUE,
                  interactions.reg = TRUE,
                  min.variation = 0,
                  gr.method = 'cor',
                  thres = 0.7){

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

  # If associations is NULL create a list of associations NULL
  if (is.null(associations)){
    associations=vector('list',length(data.omics))
    names(associations)=names(data.omics)
  }

  # Checking that the number of samples per omic is equal to number of samples for gene expression and the number of samples for edesign
  for (i in 1:length(names(data.omics))){
    if(!ncol(data.omics[[i]]) == ncol(GeneExpression) ) {
      stop("ERROR: Samples in data.omics must be the same as in GeneExpression and in edesign")
    }
  }
  if(!is.null(edesign)){
    if(!ncol(GeneExpression) == nrow(edesign) ) {
      stop("ERROR: Samples in data.omics must be the same as in GeneExpression and in edesign")
    }
  }

  #Verify that the selected grouping type it is a valid option
  if(!gr.method %in% c('cor','pca')){stop('ERROR: The selected method for grouping is not one of cor or pca')}

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

  ## Checking if there are regulators with "_R", "_P" or "_N" or with ":" and checking that there are not replicate identifiers compared to GeneExpression
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
      associations[[i]][[2]]=gsub(':', '-', associations[[i]][[2]])
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
          rownames(data.omics[[i]])[rownames(data.omics[[i]])%in%repeated] =  paste(names(data.omics)[i],'-', rownames(data.omics[[i]])[rownames(data.omics[[i]])%in%repeated],sep='')
          rownames(data.omics[[j]])[rownames(data.omics[[j]])%in%repeated] =  paste(names(data.omics)[j],'-', rownames(data.omics[[j]])[rownames(data.omics[[j]])%in%repeated],sep = '')

        }
      }
    }
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

  ## Removing genes with too many NAs and keeping track
  min.obs = ncol(GeneExpression)
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
    des.mat = model.matrix(~Group)[, -1, drop = FALSE]
    rownames(des.mat) = colnames(GeneExpression)
    #Change the name to avoid conflicts with RegulationPerCondition
    colnames(des.mat)= sub('Group','Group_', colnames(des.mat))
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

  GlobalSummary$GoodnessOfFit = matrix(NA, ncol = 4, nrow = nGenes)
  rownames(GlobalSummary$GoodnessOfFit) = Allgenes
  colnames(GlobalSummary$GoodnessOfFit) = c("Rsquared", "RMSE","CV(RMSE)", "relReg")

  GlobalSummary$ReguPerGene = matrix(0, ncol = 3*length(data.omics), nrow = nGenes)
  rownames(GlobalSummary$ReguPerGene) = Allgenes
  colnames(GlobalSummary$ReguPerGene) = c(paste(names(data.omics), "Ini", sep = "-"),
                                          paste(names(data.omics), "Mod", sep = "-"),
                                          paste(names(data.omics), "Rel", sep = "-"))

  ## Specific results for each gene
  ResultsPerGene=vector("list", length=length(Allgenes))
  names(ResultsPerGene) = Allgenes

  pap = c(1, 1:round(nGenes/100) * 100, nGenes)

  for (i in 1:nGenes) {

    gene=Allgenes[i]
    ResultsPerGene[[i]] = vector("list", length = 5)
    names(ResultsPerGene[[i]]) = c("Y", "X", "coefficients", "allRegulators", "relevantRegulators")

    if (is.element(i, pap)) cat(paste("Fitting model for gene", i, "out of", nGenes, "\n"))

    RetRegul = GetAllReg(gene=gene, associations=associations, data.omics = data.omics)
    RetRegul.gene = RetRegul$Results  ## RetRegul$TableGene: nr reg per omic
    ## Some of these reg will be removed, because they are not in data.omics


    # RetRegul.gene--> gene/regulator/omic/area
    RetRegul.gene=RetRegul.gene[RetRegul.gene[,"regulator"]!= "No-regulator", ,drop=FALSE] ## Remove rows with no-regulators


    ### NO INITIAL REGULATORS
    if(length(RetRegul.gene)==0){ ## En el caso de que no hayan INICIALMENTE reguladores -> Calcular modelo con variables experiment

      if (is.null(edesign)) {
        ResultsPerGene[[i]]$X = NULL
        ResultsPerGene[[i]]$relevantRegulators = NULL
        ResultsPerGene[[i]]$allRegulators = NULL
        isModel = NULL

      } else {
        des.mat2 = cbind(t(GeneExpression[gene,]), des.mat)
        colnames(des.mat2)[1] = "response"
        des.mat2 = na.omit(des.mat2)

        # Removing predictors with constant values
        sdNo0 = apply(des.mat2, 2, sd)
        sdNo0 = names(sdNo0)[sdNo0 > 0]
        des.mat2 = des.mat2[,sdNo0]

        isModel = NULL

        ResultsPerGene[[i]]$X = des.mat2[,-1, drop = FALSE]
        ResultsPerGene[[i]]$relevantRegulators = NULL
        ResultsPerGene[[i]]$allRegulators = NULL
      }


      # GlobalSummary$ReguPerGene  # this is initially set to 0 so no need to modify it

      ### WITH INITIAL REGULATORS
    } else { ## There are regulators for this gene at the beginning

      ResultsPerGene[[i]]$allRegulators = data.frame(RetRegul.gene, rep("Model",nrow(RetRegul.gene)), stringsAsFactors = FALSE)
      colnames(ResultsPerGene[[i]]$allRegulators) = c("gene","regulator","omic","area","filter")

      GlobalSummary$ReguPerGene[gene, grep("-Ini", colnames(GlobalSummary$ReguPerGene))] = as.numeric(RetRegul$TableGene[-1])
      # the rest of columns remain 0

      ## Identify which regulators where removed because of missing values or low variation
      res = RemovedRegulators(RetRegul.gene = ResultsPerGene[[i]]$allRegulators,
                              myregLV=myregLV, myregNA=myregNA, data.omics=data.omics)

      if(length(res$RegulatorMatrix)==0){ ## No regulators left after the filtering to compute the model

        if (is.null(edesign)) {
          ResultsPerGene[[i]]$X = NULL
          ResultsPerGene[[i]]$relevantRegulators = NULL
          ResultsPerGene[[i]]$allRegulators = res$SummaryPerGene
          isModel = NULL

        } else {
          des.mat2 = cbind(t(GeneExpression[gene,]), des.mat)
          colnames(des.mat2)[1] = "response"
          des.mat2 = na.omit(des.mat2)

          # Removing predictors with constant values
          sdNo0 = apply(des.mat2, 2, sd)
          sdNo0 = names(sdNo0)[sdNo0 > 0]
          des.mat2 = des.mat2[,sdNo0]

          isModel = NULL

          GlobalSummary$GenesNOmodel = rbind(GlobalSummary$GenesNOmodel,
                                             data.frame("gene" = gene,
                                                        "problem" = 'No regulators left after NA/LowVar filtering'))

          ResultsPerGene[[i]]$X = des.mat2[,-1, drop = FALSE]
          ResultsPerGene[[i]]$relevantRegulators = NULL
          ResultsPerGene[[i]]$allRegulators = res$SummaryPerGene

        }

      } else {  ## Regulators for the model!!

        ## Compute only if there is more than one regulator
        ResultsPerGene[[i]]$allRegulators = res$SummaryPerGene
        if(ncol(res$RegulatorMatrix)>1){
          #Save data needed for running ASGL
          group_index = Creategroups(data = res$RegulatorMatrix, reg.table = res$SummaryPerGene, correlation = thres, method = gr.method, omic.type = omic.type)
          names(group_index) = colnames(res$RegulatorMatrix)
        } else{
          group_index = c(1)
          names(group_index) = colnames(res$RegulatorMatrix)
        }

        ## Create interactions matrix without taking into account which regulators are correlated

        des.mat2 = RegulatorsInteractions(interactions.reg, reguValues = res$RegulatorMatrix,
                                          des.mat, GeneExpression, gene)

        group_index = sapply(colnames(des.mat2)[-1], function(x) find_group(x, group_index, colnames(des.mat)))
        y = des.mat2[,1,drop=FALSE]
        des.mat2 = scale(des.mat2[,-1, drop = FALSE],scale=scale, center=center)

        train.idx = sample(nrow(des.mat2), floor(nrow(des.mat2)*0.7))

        # Input data for the iterative
        data.train = list(x=as.matrix(des.mat2[train.idx,]), y=y[train.idx, ])
        data.validate = list(x=as.matrix(des.mat2[-train.idx,]), y=y[-train.idx,])

        reguexp = colnames(des.mat2)
        #Call the script to obtain the coeficients
        coefs =  sglfast::isgl(data.train, data.validate, group_index, type = "linear")$beta

        regulatorcoef = data.frame(regulators = reguexp, coefficients = coefs)
        mycoef = regulatorcoef[which(regulatorcoef[,2] != 0),1] # selected coefficients

        if (length(mycoef) == 0) {
          isModel = NULL
          GlobalSummary$GenesNOmodel = rbind(GlobalSummary$GenesNOmodel,
                                                                     data.frame("gene" = gene, "problem" = "No predictors after isgl"))
          ## Extracting relevant regulators
          ResultsPerGene[[i]]$relevantRegulators = NULL
          ResultsPerGene[[i]]$allRegulators = data.frame(ResultsPerGene[[i]]$allRegulators, "Rel" = 0, stringsAsFactors = FALSE)

          ## Counting original regulators in the model per omic
          contando = ResultsPerGene[[i]]$allRegulators[which(ResultsPerGene[[i]]$allRegulators[,"filter"] == "Model"),]
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(data.omics)])
          contando[is.na(contando)] = 0
          GlobalSummary$ReguPerGene[gene, grep("-Mod", colnames(GlobalSummary$ReguPerGene))] = contando

          ## Counting significant regulators per omic
          GlobalSummary$ReguPerGene[gene, grep("-Sig", colnames(GlobalSummary$ReguPerGene))] = NA

        } else {
          isModel = TRUE

          myvariables = unlist(strsplit(mycoef, ":", fixed = TRUE))
          myvariables = intersect(myvariables, ResultsPerGene[[i]]$allRegulators[,'regulator'])


          ResultsPerGene[[i]]$allRegulators = data.frame(ResultsPerGene[[i]]$allRegulators, "Rel" = 0, stringsAsFactors = FALSE)
          ResultsPerGene[[i]]$allRegulators[myvariables, "Rel"] = 1

          ResultsPerGene[[i]]$relevantRegulators = ResultsPerGene[[i]]$allRegulators[which(ResultsPerGene[[i]]$allRegulators[,'Rel']==1),'regulator']

          contando = ResultsPerGene[[i]]$allRegulators[which(ResultsPerGene[[i]]$allRegulators[,"filter"] == "Model"),]
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(data.omics)])
          contando[is.na(contando)] = 0
          GlobalSummary$ReguPerGene[gene, grep("-Mod", colnames(GlobalSummary$ReguPerGene))] = contando

          ## Counting relevant regulators per omic
          if (length(ResultsPerGene[[i]]$relevantRegulators) > 0) {
            contando = ResultsPerGene[[i]]$allRegulators[ResultsPerGene[[i]]$relevantRegulators,]
            contando = table(contando[,"omic"])
            contando = as.numeric(contando[names(data.omics)])
            contando[is.na(contando)] = 0
            GlobalSummary$ReguPerGene[gene, grep("-Rel", colnames(GlobalSummary$ReguPerGene))] = contando
          }
        }

         }


      } ## Close "else" --> None regulators from begining

      if (is.null(isModel)) {

        ResultsPerGene[[i]]$Y = GeneExpression[i,]
        ResultsPerGene[[i]]$coefficients = NULL

        GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[rownames(GlobalSummary$GoodnessOfFit) != gene,, drop = FALSE]


      } else {
        y.fitted = des.mat2%*%coefs
        ResultsPerGene[[i]]$Y = data.frame("y" = y, "fitted.y" = y.fitted, "residuals" = y - des.mat2%*%coefs, check.names = FALSE)
        colnames(ResultsPerGene[[i]]$Y) <- c("y", "fitted.y", "residuals")
        rownames(regulatorcoef) = regulatorcoef[,1]
        ResultsPerGene[[i]]$coefficients = regulatorcoef[mycoef,2, drop = FALSE]

        R.squared = round(1-(sum( y - y.fitted)/sum((y-mean(y[,1]))^2)),6)
        RMSE = round(sqrt(sum((y-y.fitted)^2)/nrow(y)),6)
        cvRMSE = abs(round(sqrt(sum((y-y.fitted)^2)/nrow(y))/mean(y[,1]),6))

        GlobalSummary$GoodnessOfFit[gene,] = c(R.squared, RMSE, cvRMSE,length(ResultsPerGene[[gene]]$relevantRegulators))

      }

  }  ## At this point the loop for all genes is finished

  # Remove from GoodnessOfFit genes with no relevant regulators

  genesNosig = names(which(GlobalSummary$GoodnessOfFit[,1]==0))
  genessig = setdiff(rownames(GlobalSummary$GoodnessOfFit), genesNosig)
  GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[genessig,, drop=FALSE]

  #Calculate GlobalRegulators
  m_rel_reg<-lapply(ResultsPerGene, function(x) x$relevantRegulators)
  m_rel_reg <- unlist(m_rel_reg)
  mrel_vector <- table(m_rel_reg)
  #Calculate third quantile
  q3<-quantile(mrel_vector,0.75)
  if(length(mrel_vector[mrel_vector>q3])<10){
    GlobalSummary$GlobalRegulators = intersect(names(mrel_vector[rev(tail(order(mrel_vector),10))]), names(mrel_vector[mrel_vector>10]) )
  } else{
    GlobalSummary$GlobalRegulators = intersect(names(mrel_vector[mrel_vector>q3]), names(mrel_vector[mrel_vector>10]) )
  }

  #Calculate HubGenes
  relevant_regulators<-GlobalSummary$ReguPerGene[,c(grep('-Rel$',colnames(GlobalSummary$ReguPerGene))),drop=FALSE]
  s_rel_reg<-apply(relevant_regulators, 1, sum)
  #Calculate third quantile
  q3<-quantile(s_rel_reg,0.75)
  if(length(s_rel_reg[s_rel_reg>q3])<10){
    GlobalSummary$HubGenes = intersect(names(s_rel_reg[rev(tail(order(s_rel_reg),10))]), names(s_rel_reg[s_rel_reg>10]) )
  } else{
    GlobalSummary$HubGenes = intersect(names(s_rel_reg[s_rel_reg>q3]), names(s_rel_reg[s_rel_reg>10]))
  }
  if(all(sapply(associations,is.null))) {associations = NULL}

  myarguments = list(edesign = edesign, finaldesign = des.mat, groups = Group,
                     center = center, scale = scale, clinic.type = clinic.type,
                     min.variation = min.variation, associations = associations,
                     thres = thres, gr.method = gr.method,
                     GeneExpression = GeneExpression, dataOmics = data.omics, omic.type = omic.type,
                     clinic = clinic, clinic.type = clinic.type,  method = 'isgl')

  # Create the results for the scale filter check

  result <- list("ResultsPerGene" = ResultsPerGene, "GlobalSummary" = GlobalSummary, "arguments" = myarguments)
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

Creategroups = function(data, reg.table, method = 'cor' ,correlation =0.8, omic.type){


  # Take into account only regulators that could enter the model
  myreg = colnames(data)
  # Apply required method (COR: groups based on correlation, PCA: groups based on PCA)

  if (method == 'pca'){

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
      d = min(which(cumsum(respca@modelDF[,'R2X'])>correlation))

      if(d<r){
        #Restrict to the first d components that extract the 80% of the variance
        respca@loadingMN <- respca@loadingMN[,1:d, drop=FALSE]
      }


      #Create the groups for the regulators according to their maximal absolute value on the loadings
      groups = max.col(abs(respca@loadingMN))
    }

  }

  if (method == 'cor'){

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
