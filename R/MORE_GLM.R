#########################################################################################
######           Functions to integrate omics data using GLMs                      ######
#########################################################################################


## By Sonia, Maider & Monica
## 07-July-2016
## Last modified: October 2023


options(stringsAsFactors = FALSE)

library(ltm)

#' Generalized Linear Models
#'
#'\code{GetGLM} fits a regression model for all the genes in the data set to identify
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
#' @param family Error distribution and link function to be used in the model when \code{\link{method}} glm. By default, gaussian().
#' @param elasticnet ElasticNet mixing parameter. There are three options:
#' \itemize{
#' \item NULL : The parameter is selected from a grid of values ranging from 0 to 1 with 0.1 increments. The chosen value optimizes the mean cross-validated error when optimizing the lambda values.
#' \item A number between 0 and 1 : ElasticNet is applied with this number being the combination between Ridge and Lasso penalization (elasticnet=0 is the ridge penalty, elasticnet=1 is the lasso penalty).
#' \item A vector with the mixing parameters to try. The one that optimizes the mean cross-validated error when optimizing the lambda values will be used.
#' }
#' #' By default, NULL.
#' @param interactions.reg If TRUE, the model includes interactions between regulators and experimental variables. By default, TRUE.
#' @param min.variation  For numerical regulators, it specifies the minimum change required across conditions to retain the regulator in
#' the regression models. In the case of binary regulators, if the proportion of the most common value is equal to or inferior this value,
#' the regulator is considered to have low variation and will be excluded from the regression models. The user has the option to set a single
#' value to apply the same filter to all omics, provide a vector of the same length as omics if they want to specify different levels for each omics,
#' or use 'NA' when they want to apply a minimum variation filter but are uncertain about the threshold. By default, 0.
#' @param col.filter Type of correlation coefficients to use when applying the multicollinearity filter when glm \code{\link{method}} is used.
#' \itemize{
#' \item cor: Computes the correlation between omics. Pearson correlation between numeric variables, phi coefficient between numeric and binary and biserial correlation between binary variables.
#' \item pcor : Computes the partial correlation.
#' }
#' @param correlation  Value to determine the presence of collinearity between two regulators when using the glm \code{\link{method}}. By default, 0.7.
#' @param scaletype Type of scaling to be applied. Three options:
#' \itemize{
#' \item auto : Applies the autoscaling.
#' \item pareto : Applies the pareto scaling. \deqn{\frac{X_k}{s_k \sqrt[4]{m_b}} }
#' \item block : Applies the block scaling. \deqn{ \frac{X_k}{s_k \sqrt{m_b}} }
#' }
#' considering m_b the number of variables of the block. By default, auto.
#' @return List containing the following elements:
#' \itemize{
#' \item ResultsPerGene : List with as many elements as genes in \code{\link{GeneExpression}}. For each gene, it includes information about gene values, considered variables, estimated coefficients,
#'                    detailed information about all regulators, and regulators identified as relevant (in glm scenario) or significant (in pls scenarios).
#' \item GlobalSummary : List with information about the fitted models, including model metrics, information about regulators, genes without models, regulators, master regulators and hub genes.
#' \item Arguments : List containing all the arguments used to generate the models.
#' }
#' @export

GetGLM = function(GeneExpression,
                  data.omics,
                  associations = NULL,
                  omic.type = 0,
                  edesign = NULL,
                  clinic = NULL,
                  clinic.type =NULL,
                  center = TRUE, scale = TRUE,
                  epsilon = 0.00001,
                  family = gaussian(),
                  elasticnet = NULL,
                  interactions.reg = TRUE,
                  min.variation = 0,
                  col.filter = 'cor',
                  correlation = 0.7,
                  scaletype = 'auto'){


  # Converting matrix to data.frame
  GeneExpression = as.data.frame(GeneExpression)
  data.omics = lapply(data.omics, as.data.frame)

  ## Omic types
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
    if(!length(colnames(data.omics[[i]])) == length(colnames(GeneExpression)) ) {
      stop("ERROR: Samples in data.omics must be the same as in GeneExpression and in edesign")
    }
  }
  if(!is.null(edesign)){
    if(!length(colnames(GeneExpression)) == length(rownames(edesign)) ) {
      stop("ERROR: Samples in data.omics must be the same as in GeneExpression and in edesign")
    }
  }

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
      cat(x, "and Gene Expression omics have shared identifiers in regulators:", repeated, "\n")
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

  ## Removing genes with NAs and keeping track
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
    cat(length(genesNOreg), "genes had no initial regulators. Models will be computed for", length(rownames(GeneExpression)), 'genes.\n')
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


  ### Computing model for each gene
  cat("Checking multicollinearity, selecting predictors and fitting model for ...\n")

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
    if(length(RetRegul.gene)==0){ ## En el caso de que no hayan INICIALMENTE reguladores -> Calcular modelo con variables experimentales o clinicas.

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

        isModel =NULL
        ResultsPerGene[[i]]$X = des.mat2[,-1, drop = FALSE]
        ResultsPerGene[[i]]$relevantRegulators = NULL
        ResultsPerGene[[i]]$allRegulators = NULL
      }


      # GlobalSummary$ReguPerGene  # this is initially set to 0 so no need to modify it


      ### WITH INITIAL REGULATORS
    }
    else { ## There are regulators for this gene at the beginning

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

      }
      else {  ## Regulators for the model!!

        ## Apply multicollinearity filter only if there is more than one regulator for a gene

        if (ncol(res$RegulatorMatrix)>1){
          if(col.filter=='cor'){
            res = CollinearityFilter1(data = res$RegulatorMatrix, reg.table = res$SummaryPerGene,
                                      correlation = correlation, omic.type = omic.type, scale = scale, center = center)
          }
          if(col.filter=='pcor'){
            res = CollinearityFilter2(data = res$RegulatorMatrix, reg.table = res$SummaryPerGene,
                                      correlation = correlation, omic.type = omic.type, epsilon = epsilon , scale = scale, center = center)
          }

        }

        if(is.null(res)){
          des.mat2 = cbind(t(GeneExpression[gene,]), des.mat)
          colnames(des.mat2)[1] = "response"

          colnames(des.mat2) = gsub("\`", "", colnames(des.mat2))

          GlobalSummary$GenesNOmodel = rbind(GlobalSummary$GenesNOmodel,
                                             data.frame("gene" = gene,
                                                        "problem" = 'Problem with Partial Correlation calculation'))
          isModel =NULL
        } else{
          ResultsPerGene[[i]]$allRegulators = res$SummaryPerGene

          ## Scaling predictors for ElasticNet only in case they were not already scaled
          des.mat2EN = RegulatorsInteractions(interactions.reg, reguValues = res$RegulatorMatrix,
                                                des.mat, GeneExpression, gene)

          # Removing observations with missing values
          des.mat2EN = na.omit(des.mat2EN)

          #Scale the variables, indispensable for elasticnet application
          des.mat2EN = data.frame(des.mat2EN[,1,drop=FALSE], scale(des.mat2EN[,-1,drop=FALSE],scale=scale,center=center),check.names = FALSE)

          ##Scale if needed to block scaling or pareto scaling
          if (scaletype!='auto'){
            ## Make the groups of omics
            regupero = lapply(unique(res$SummaryPerGene[,'omic']), function(x) rownames(res$SummaryPerGene)[res$SummaryPerGene[,'omic'] == x & res$SummaryPerGene[,'filter'] == "Model"])
            names(regupero) = unique(res$SummaryPerGene[,'omic'])
            #Remove empty omics
            regupero = regupero[sapply(regupero, function(x) length(x) > 0)]
            #It does not work in case of really huge amount of data
            regupero1 = try(suppressWarnings( lapply(regupero, function(x) colnames(des.mat2EN[,grep(paste(x, collapse = "|"), colnames(des.mat2EN)),drop=FALSE]))),silent = TRUE)
            if(class(regupero1)=='try-error'){
              #Add the ones related to the interactions
              regupero = filter_columns_by_regexp(regupero, des.mat2EN,res)
            }else{regupero = regupero1}
            res$RegulatorMatrix = Scaling.type(des.mat2EN[,-1,drop=FALSE], regupero, scaletype)

            #Use them jointly
            des.mat2EN = data.frame(des.mat2EN[,1,drop=FALSE], scale(des.mat,scale=scale,center=center), res$RegulatorMatrix,check.names = FALSE)
            rm(regupero);gc()
          }

          ###  Variable selection --> Elasticnet
          tmp = ElasticNet(family2, des.mat2EN, epsilon, elasticnet)
          regulatorcoef = tmp[['coefficients']]
          isModel = tmp[['isModel']]
          m = tmp[['m']]
          des.mat2 = as.data.frame(des.mat2EN[,colnames(tmp[["des.mat2"]]),drop = FALSE])
          ResultsPerGene[[i]]$X = des.mat2EN[,-1, drop = FALSE]
          rm(des.mat2EN); gc()

        }


        if (ncol(des.mat2) == 1 || is.null(isModel)) {

          ## Extracting significant regulators
          ResultsPerGene[[i]]$relevantRegulators = NULL
          ResultsPerGene[[i]]$allRegulators = data.frame(ResultsPerGene[[i]]$allRegulators, "Rel" = 0, stringsAsFactors = FALSE)

          ## Counting original regulators in the model per omic
          contando = ResultsPerGene[[i]]$allRegulators[which(ResultsPerGene[[i]]$allRegulators[,"filter"] == "Model"),]
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(data.omics)])
          contando[is.na(contando)] = 0
          GlobalSummary$ReguPerGene[gene, grep("-Mod", colnames(GlobalSummary$ReguPerGene))] = contando
        } else{

          isModel = TRUE
          mycoef = colnames(des.mat2[,-1,drop = FALSE])
          myvariables = unlist(strsplit(mycoef, ":", fixed = TRUE))
          mycondi = intersect(myvariables, colnames(des.mat))
          myvariables = intersect(myvariables, rownames(ResultsPerGene[[i]]$allRegulators))

          ResultsPerGene[[i]]$allRegulators = data.frame(ResultsPerGene[[i]]$allRegulators, "Rel" = 0, stringsAsFactors = FALSE)
          ResultsPerGene[[i]]$allRegulators[myvariables, "Rel"] = 1
          ResultsPerGene[[i]]$coefficients = regulatorcoef
          #ResultsPerGene[[i]]$coefficients = regulatorcoef[myvariables,, drop =FALSE]
          colnames(ResultsPerGene[[i]]$coefficients) = c('coefficient')

          ## A las variables significativas le quito "_R", solo quedara omica_mc"num". Luego creo un objeto que contenga a mi tabla de "allRegulators"
          ## para poder modificar los nombres de la misma forma.
          myvariables = sub("_R", "", myvariables)
          ResultsPerGene[[i]]$relevantRegulators = myvariables # significant regulators including "new" correlated regulators without _R

          contando = ResultsPerGene[[i]]$allRegulators[which(ResultsPerGene[[i]]$allRegulators[,"filter"] == "Model"),]
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(data.omics)])
          contando[is.na(contando)] = 0
          GlobalSummary$ReguPerGene[gene, grep("-Mod", colnames(GlobalSummary$ReguPerGene))] = contando

          mytable = ResultsPerGene[[i]]$allRegulators
          mytable[,"filter"] = sub("_P", "", mytable[,"filter"])
          mytable[,"filter"] = sub("_N", "", mytable[,"filter"])
          mytable[,"filter"] = sub("_R", "", mytable[,"filter"])

          collin.regulators = intersect(myvariables, mytable[,"filter"])

          if (length(collin.regulators) > 0) {  # there were correlated regulators
            original.regulators = mytable[mytable[,"filter"] %in% collin.regulators, "regulator"]

            ResultsPerGene[[i]]$allRegulators[original.regulators, "Rel"] = 1
            ResultsPerGene[[i]]$relevantRegulators = c(ResultsPerGene[[i]]$relevantRegulators, as.character(original.regulators))
            ResultsPerGene[[i]]$relevantRegulators = setdiff(ResultsPerGene[[i]]$relevantRegulators, collin.regulators)

            ## Counting original regulators in the model per omic
            contando = ResultsPerGene[[i]]$allRegulators
            quitar = which(contando[,"filter"] == "MissingValue")
            if (length(quitar) > 0) contando = contando[-quitar,]
            quitar = which(contando[,"filter"] == "LowVariation")
            if (length(quitar) > 0) contando = contando[-quitar,]
            contando = contando[-grep("_R", rownames(contando)),]
          } else {
            contando = ResultsPerGene[[i]]$allRegulators[which(ResultsPerGene[[i]]$allRegulators[,"filter"] == "Model"),]
          }
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(data.omics)])
          contando[is.na(contando)] = 0
          GlobalSummary$ReguPerGene[gene, grep("-Mod", colnames(GlobalSummary$ReguPerGene))] = contando

          ## TO DO: Corregir GblobalSummary$ReguPerGene
          ## Counting significant regulators per omic
          if (length(ResultsPerGene[[i]]$relevantRegulators) > 0) {
            contando = ResultsPerGene[[i]]$allRegulators[ResultsPerGene[[i]]$relevantRegulators,, drop=FALSE]
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

      GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[rownames(GlobalSummary$GoodnessOfFit) != gene,,drop=FALSE]


    } else {
      ResultsPerGene[[i]]$Y = data.frame("y" = des.mat2[,1], "fitted.y" = tmp[['fitted.values']],
                                         "residuals" = des.mat2[,1] - tmp[['fitted.values']], check.names = FALSE)
      colnames(ResultsPerGene[[i]]$Y) <- c("y", "fitted.y", "residuals")
      GlobalSummary$GoodnessOfFit[gene,] = c(m$R.squared, m$RMSE, m$cvRMSE,length(ResultsPerGene[[gene]]$relevantRegulators))


    }

  }  ## At this point the loop for all genes is finished

  # Remove from GoodnessOfFit genes with no significant regulators

  genesNosig = names(which(GlobalSummary$GoodnessOfFit[,4]==0))
  genessig = setdiff(rownames(GlobalSummary$GoodnessOfFit), genesNosig)
  GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[genessig,,drop=FALSE]

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
  relevant_regulators<-GlobalSummary$ReguPerGene[,c(grep('-Rel$',colnames(GlobalSummary$ReguPerGene)))]
  s_rel_reg<-apply(relevant_regulators, 1, sum)
  #Calculate third quantile
  q3<-quantile(s_rel_reg,0.75)
  if(length(s_rel_reg[s_rel_reg>q3])<10){
    GlobalSummary$HubGenes = intersect(names(s_rel_reg[rev(tail(order(s_rel_reg),10))]), names(s_rel_reg[s_rel_reg>10]) )
  } else{
    GlobalSummary$HubGenes = intersect(names(s_rel_reg[s_rel_reg>q3]), names(s_rel_reg[s_rel_reg>10]))
  }

  myarguments = list(edesign = edesign, finaldesign = des.mat, groups = Group, family = family,
                     center = center, scale = scale, elasticnet = tmp[['elasticnet']],
                     min.variation = min.variation, correlation = correlation,
                     epsilon = epsilon, associations = associations,
                     GeneExpression = GeneExpression, dataOmics = data.omics, omic.type = omic.type,
                     clinic = clinic, clinic.type=clinic.type,method ='glm')

  result <- list("ResultsPerGene" = ResultsPerGene, "GlobalSummary" = GlobalSummary, "arguments" = myarguments)
  class(result) <- "MORE"
  return(result)


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
  ## reg.table = Table with "gene", "regulator", "omic", "area", filter" where omics with no regulators have been removed
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

        ## Take as representator the one with highest correlations, it case of tie, select it randomly
        sums = sapply(correlacionados, function(x) sum(abs(mycor[which(apply(mycor[,c(1,2)]==c(x),1,any)),3])))
        if(length(which(sums==max(sums)))>1){
          keep = sample(names(which(sums==max(sums))),1)
        } else {
          keep = names(which(sums==max(sums)))
        }

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
          correlacionados = names(which(igraph::as_adj(mysubgraph)[,repre]>0))
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

  resultado = list(RegulatorMatrix = data, SummaryPerGene = reg.table)
  rownames(resultado$SummaryPerGene) = resultado$SummaryPerGene[,"regulator"]
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
  ## reg.table = Table with "gene", "regulator", "omic", "area", filter" where omics with no regulators have been removed
  row.names(reg.table) = reg.table[,"regulator"]
  #resultado = list(RegulatorMatrix = data, SummaryPerGene = reg.table)

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

          correlacionados = names(which(igraph::as_adj(mysubgraph)[,repre]>0))
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

  resultado = list(RegulatorMatrix = data, SummaryPerGene = reg.table)
  rownames(resultado$SummaryPerGene) = resultado$SummaryPerGene[,"regulator"]
  return(resultado)
}

