#########################################
#### Auxiliar functions
#########################################

## By Sonia & Maider
## 15-Oct-2023

# Removing regulators with low variation ----------------------------------

LowVariationRegu = function(min.variation, data.omics, ExpGroups, associations, Allgenes, omic.type, clinic.type) {
  
  if (all(is.na(min.variation))){
    
    for (ov in names(associations)){
      if(!is.null(associations[[ov]])){
        myreg=associations[[ov]][associations[[ov]][,1] %in% Allgenes,2] # removing regulators not associated to our genes
        data.omics[[ov]]=data.omics[[ov]][intersect(myreg, rownames(data.omics[[ov]])),] ## Reduced data.omics
        rm("myreg")
      }
    }
    #### Low variation cutoff is computed automatically
    data.omicsMean = vector("list", length=length(data.omics))
    if(!is.null(clinic.type)){j = 2}else{j=1}
    for(i in j:length(data.omics)){
      data.omicsMean[[i]]=t(apply(data.omics[[i]], 1, tapply, ExpGroups, mean))
    }
    names(data.omicsMean)=names(data.omics)
    
    percVar = c(10, 0.1)
    names(percVar) = 0:1
    percVar = percVar[as.character(omic.type)]
    names(percVar)=names(data.omicsMean)
    
    # Applying Low Variation filter
    LowVar=LowVariatFilter(data=data.omicsMean, method="sd", percVar=percVar, omic.type = omic.type, clinic.type = clinic.type)
    
    # data.omicsMean reduced: without NA and LV
    data.omicsMean=LowVar$data
    
    ## data.omics reduced: only mygenes and without NA and LV
    for (ov in names(data.omics)[j:length(data.omics)]){
      data.omics[[ov]] = data.omics[[ov]][rownames(data.omicsMean[[ov]]),]
      # Remove regulators from associations that have been removed due to LowVariation
      associations[[ov]] = associations[[ov]][associations[[ov]][,2] %in% rownames(data.omics[[ov]]),,drop = FALSE]
      
    }
    
  }
  else {
    
    #### Low variation cutoff is set by the user
    
    # removing regulators not associated to our genes only when there is associations matrix
    for (ov in names(associations)){
      if(!is.null(associations[[ov]])){
        myreg=associations[[ov]][associations[[ov]][,1] %in% Allgenes,2]
        data.omics[[ov]]=data.omics[[ov]][intersect(myreg, rownames(data.omics[[ov]])),]
        rm("myreg")
      }
    }
    
    # Creating vector for min.variation
    if (length(min.variation) == 1) {  ## Including min.variation = 0. I need a vector with omics names
      min.variation=rep(min.variation,length(data.omics))
      names(min.variation)=names(data.omics)
    } 
    
    # computing mean per condition in data.omics
    data.omicsMean=vector("list", length = length(data.omics))
    if(!is.null(clinic.type)){j = 2}else{j=1}
    for(i in j:length(data.omics)){
      data.omicsMean[[i]] = t(apply(data.omics[[i]], 1, tapply, ExpGroups, mean))
    }
    names(data.omicsMean) = names(data.omics)
    
    # Applying Low Variation filter
    LowVar=LowVariatFilter(data = data.omicsMean, method = "user", percVar = min.variation, omic.type = omic.type, clinic.type = clinic.type)
    
    data.omicsMean=LowVar$data  ## data.omicsMean reduced
    
    ## data.omics reduced: only mygenes and without NA and LV
    for (ov in names(data.omics)[j:length(data.omics)]){
      data.omics[[ov]] = data.omics[[ov]][rownames(data.omicsMean[[ov]]),]
      # Remove regulators from associations that have been removed due to LowVariation
      associations[[ov]] = associations[[ov]][associations[[ov]][,2] %in% rownames(data.omics[[ov]]),,drop = FALSE]
      
    }
    
  }
  
  rm("data.omicsMean"); gc()
  
  # Regulators removed due to low variation filter
  myregLV=LowVar$LV.reg
  rm("LowVar"); gc()
  
  cat("Number of regulators with low variation:\n")
  print(sapply(myregLV, length))
  cat("\n")
  
  return(list("myregLV" = myregLV, "data.omics" = data.omics, "associations" = associations))
}





# Low variation filtering -------------------------------------------------
## data: For computing variability. It could be mean values between replicates
## method: One of "sd","range", "IQrange" or "user"
## percVar: percentage of variation defined by the user

LowVariatFilter=function(data, method, percVar, omic.type, clinic.type){
  
  SummaryRes = LV.reg = vector("list", length=length(data))
  names(SummaryRes) = names(LV.reg) = names(data)
  
  if(!is.null(clinic.type)){i = 2}else{i=1}
  
  for (ov in names(data)[i:length(data)]) {
    if(is.na(percVar[ov])){
      method.low="sd"
      if(omic.type[ov]==0){
        percVar[ov]=10
      }else{
        percVar[ov]==0.1
      }
    }else{
      method.low =method
    }
    if (omic.type[ov] == 0) {  # numerical regulators
      
      if (method.low=="sd") {
        met=apply(data[[ov]], 1, sd, na.rm=TRUE)  ## Compute standard deviation between conditions
        maxMet=max(met)*(percVar[ov]/100)  ## Compute minimum variation allowed
        myreg=met[met>maxMet]  # Regulators to be kept
        LV.reg[[ov]]=names(met[met<=maxMet]) ## Keep names of removed regulators
        data[[ov]]=data[[ov]][names(myreg), ,drop=FALSE]
      }
      
      if(method.low=="user") {
        if (min(dim(data[[ov]])) > 0) {
          met = apply(data[[ov]], 1, function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE) )
          maxMet = percVar[ov] ## We don't consider the max, just the percentage defined by the user
          myreg=met[met>maxMet]
          LV.reg[[ov]]=names(met[met<=maxMet]) ## Keep names of removed regulators
          data[[ov]]=data[[ov]][names(myreg), , drop=FALSE]
        }
      }
    }
    
    if (omic.type[ov] == 1) {  # binary categorical regulators
      
      if(method.low=='sd'){
        met = apply(data[[ov]], 1, function (x) { max(x, na.rm = TRUE)-min(x, na.rm = TRUE) })  ## Compute maximun variation between groups
        maxMet=max(met)/10 ## Compute the minimum variation allowed
        myreg = met[met > maxMet]  # Regulators to be kept
        LV.reg[[ov]] = names(met[met <= maxMet]) ## Keep names of removed regulators
        data[[ov]] = data[[ov]][names(myreg), ,drop=FALSE]
      }
      
      if(method.low == 'user'){ 
        met = apply(data[[ov]], 1, function (x) { max(x, na.rm = TRUE)-min(x, na.rm = TRUE) })  ## Compute maximun variation between groups
        myreg = met[met > percVar[ov]]  # Regulators to be kept
        LV.reg[[ov]] = names(met[met <= percVar[ov]]) ## Keep names of removed regulators
        data[[ov]] = data[[ov]][names(myreg), ,drop=FALSE]
      }
    }
  }
  
  results=vector("list", length=2)
  results[[1]]=data
  results[[2]]=LV.reg
  names(results)=c("data", "LV.reg")
  
  return(results)
}

# Get all regulators ------------------------------------------------------


## Auxiliar function to paste different areas

uniquePaste = function (x) {
  x = unique(x)
  
  if (length(x) > 1) {
    x = paste(x, collapse = ";")
  }
  
  return (x)
}





## INPUT
## gene: we are seeking information
## associations: List containing as many elements as original association files between genes and each omic (association$miRNA, association$DNase, etc.)
## data.omics: For removing regulators with no-information

GetAllReg=function(gene, associations, data.omics){
  
  Reg.matrix=NULL
  NrReg=NULL
  myomic=names(data.omics)
  
  for(ov in myomic){
    if (is.null(associations[[ov]])){
      myregulators=rownames(data.omics[[ov]])
      ov.nr=length(myregulators)
      
      if(ov.nr==0){
        myregulators=c("No-regulator")
        myregulators=t(as.matrix(myregulators))
      }
      
      Reg.matrix.temp=cbind(rep(gene,ov.nr),myregulators,ov,"")  ## Lo tengo que dejar igual que las otras omicas
      NrReg.temp=ov.nr
      colnames(Reg.matrix.temp)=c("gene","regulator","omic","area")
      Reg.matrix.temp = t(apply(Reg.matrix.temp,1, as.character))  ## adding this here to avoid factors
      Reg.matrix=rbind(Reg.matrix,Reg.matrix.temp)
      colnames(Reg.matrix)=c("gene","regulator","omic","area")
      NrReg=c(NrReg, NrReg.temp)
      
    } else{
      
      colnames(associations[[ov]])[1]="gene"
      
      ## Regulator with area
      
      if(ncol(associations[[ov]])>2){
        
        myregulators=associations[[ov]][associations[[ov]]$gene==gene, ,drop=FALSE] ## "regulator" with Area--> Matrix
        
        if(nrow(myregulators)==0){
          ov.nr=nrow(myregulators) ## regulators nr per omic --> It could be 0
          myregulators=c("No-regulator","")
          myregulators=t(as.matrix(myregulators))
          
        } else {
          myregulators=aggregate(myregulators,by=list(myregulators[,2]),uniquePaste) ## Agrupado por "region", me devuelve Area separado por comas
          myregulators=myregulators[,-c(1:2)]
          myregulators[,2]=sapply(myregulators[,2], paste, collapse = ";")
          ov.nr=nrow(myregulators)
        }
        
        Reg.matrix.temp = cbind(gene,myregulators,ov)
        Reg.matrix.temp = Reg.matrix.temp[,c(1,2,4,3),drop=FALSE] ## Dejo el mismo orden que tenia
        colnames(Reg.matrix.temp)=c("gene","regulator","omic","area")
        NrReg.temp=ov.nr
        Reg.matrix=rbind(Reg.matrix,Reg.matrix.temp)
        colnames(Reg.matrix)=c("gene","regulator","omic","area")
        NrReg=c(NrReg, NrReg.temp)
        
        ## Regulators without area
        
      } else {
        
        myregulators=associations[[ov]][associations[[ov]]$gene==gene,2]
        myregulators=unique(myregulators) ## Could it be repeated??
        ov.nr=length(myregulators) ## regulators nr per omic --> It could be 0
        
        if(ov.nr==0){
          myregulators=c("No-regulator")
          myregulators=t(as.matrix(myregulators))
        }
        
        Reg.matrix.temp=cbind(gene,myregulators,ov,"")  ## Lo tengo que dejar igual que las otras omicas
        NrReg.temp=ov.nr
        Reg.matrix.temp = t(apply(Reg.matrix.temp,1, as.character))  ## adding this here to avoid factors
        colnames(Reg.matrix.temp)=c("gene","regulator","omic","area")
        Reg.matrix=rbind(Reg.matrix,Reg.matrix.temp)
        colnames(Reg.matrix)=c("gene","regulator","omic","area")
        NrReg=c(NrReg, NrReg.temp)
      }
    }
  }
  Results2=c(gene,NrReg)
  myomic=paste(myomic,"Ini",sep="-")
  names(Results2)=c("gene",myomic)
  
  Results=vector("list", length=2)
  Results[[1]]=Reg.matrix
  Results[[2]]=Results2
  names(Results)=c("Results","TableGene")
  
  return(Results)
  
}






# Obtain which regulators have been removed and why -----------------------
## Input:
## RetRegul.gene: Initial regulator matrix
## myregLV: regulators removed for low variability. List by omics
## myregNA: regulators removed for NA. List by omics

RemovedRegulators = function(RetRegul.gene, myregLV, myregNA, data.omics){
  
  RegulatorsValue=NULL
  RetRegul.geneNEW = NULL
  
  rownames(RetRegul.gene) = RetRegul.gene[,"regulator"]
  myregini = RetRegul.gene[,"regulator"]  ## In our case, one regulator cannot belong to 2 omics
  mygene = RetRegul.gene[1,"gene"]
  
  for(ov in unique(RetRegul.gene[,"omic"])){
    
    ## remove regulators not in data.omics
    regmodel = intersect(RetRegul.gene[,"regulator"], rownames(data.omics[[ov]]))
    RegulatorsValue = cbind(RegulatorsValue, t(data.omics[[ov]][regmodel,,drop=FALSE]))
    RetRegul.geneNEW = rbind(RetRegul.geneNEW, RetRegul.gene[regmodel, , drop=FALSE])
    
    ## NA
    
    if(length(intersect(myregini, myregNA[[ov]]))>0){
      RetRegul.geneNEW = rbind(RetRegul.geneNEW, data.frame(gene = mygene, regulator = intersect(myregini, myregNA[[ov]]),
                                                            omic = ov, area=RetRegul.gene[intersect(myregini, myregNA[[ov]]),"area"],
                                                            filter = "MissingValue", stringsAsFactors = FALSE))
    }
    
    ## LV
    
    if (length(intersect(myregini, myregLV[[ov]]))>0){
      RetRegul.geneNEW = rbind(RetRegul.geneNEW, data.frame(gene = mygene, regulator = intersect(myregini, myregLV[[ov]]),
                                                            omic = ov, area=RetRegul.gene[intersect(myregini, myregLV[[ov]]),"area"],
                                                            filter = "LowVariation", stringsAsFactors = FALSE))
    }
    
    
    ## RetRegul.geneNEW=as.data.frame(RetRegul.geneNEW) ## Si es una matriz con una fila,al hacer el apply me lo convierte en vector
    RetRegul.geneNEW=apply(RetRegul.geneNEW,c(1,2),as.character)
    ##RetRegul.geneNEW[,"filter"] = as.character(RetRegul.geneNEW[,"filter"])
    
    
  }
  
  return(list("SummaryPerGene" = RetRegul.geneNEW, "RegulatorMatrix" = RegulatorsValue))
  
}


# Adding interations with regulators -------------------------------------

filter_columns_by_regexp <- function(regupero, des.mat2, res) {
  filtered_columns <- lapply(regupero, function(x) if (length(x) != 0) {
    pat <- paste(x, collapse = "|")
    logical_indices <- stringr::str_detect(colnames(des.mat2), stringr::regex(pat, ignore_case = TRUE))
    colnames(des.mat2)[logical_indices]
  })
  return(filtered_columns)
}

RegulatorsInteractions = function (interactions.reg, reguValues, des.mat, GeneExpression, gene) {
  
  # Adding regulators
  if (is.null(des.mat)) {
    des.mat2 = data.frame(reguValues, check.names = FALSE)
    des.mat2 = cbind(t(GeneExpression[gene,]), des.mat2)
    colnames(des.mat2)[1] = "response"
  } else {
    des.mat2 = data.frame(des.mat, reguValues, check.names = FALSE)
    ### WITH INTERACTIONS with regulators
    if (interactions.reg > 0) {
      expcond = colnames(des.mat)
      
      if (interactions.reg == 1) {  # Max order of interaction = 2 & Interactions with cont.var not allowed
        if (length(grep(":", expcond, fixed = TRUE)))  expcond = expcond[-grep(":", expcond, fixed = TRUE)]
      }
      
      if (interactions.reg == 2) { # Max order of interaction = 2
        if (length(grep(":", expcond, fixed = TRUE)))  expcond = expcond[-grep(":", expcond, fixed = TRUE)]
      }
      
      if(ncol(reguValues)<=5000){
        
        fff = paste0("~ ",
                     paste(sapply(colnames(reguValues),
                                  function (x) paste(expcond, sprintf("`%s`", x), sep = ":")),
                           collapse = "+"))
        
        fff = as.formula(fff)
        inter.var = model.matrix(fff, des.mat2)[,-1, drop = FALSE]
        colnames(inter.var) = strsplit(as.character(fff),"\\s*\\+\\s*")[-1][[1]]
        des.mat2 = cbind(des.mat2, inter.var)
        
        des.mat2 = cbind(t(GeneExpression[gene,]), des.mat2)
        colnames(des.mat2)[1] = "response"
        
        colnames(des.mat2) = gsub("\`", "", colnames(des.mat2))
        
      } else{
        j = ceiling(ncol(reguValues)/250)
        res.mat = des.mat2
        
        for (k in 1:j){
          
          cols_start = (250 * (k - 1)) + 1
          cols_end = min(250 * k, ncol(reguValues))
          
          fff = paste0("~ ",
                       paste(sapply(colnames(reguValues[, cols_start:cols_end, drop = FALSE]),
                                    function(x) paste(expcond, sprintf("`%s`", x), sep = ":")),
                             collapse = "+"))
          
          fff = as.formula(fff)
          inter.var = model.matrix(fff, des.mat2)[,-1, drop = FALSE]
          colnames(inter.var) = strsplit(as.character(fff), "\\s*\\+\\s*")[-1][[1]]
          res.mat = cbind(res.mat, inter.var)
        }
        
        des.mat2 = res.mat
        
        des.mat2 = cbind(t(GeneExpression[gene,]), des.mat2)
        colnames(des.mat2)[1] = "response"
        
        colnames(des.mat2) = gsub("\`", "", colnames(des.mat2))
        
      }
      
      ## PUEDE OCURRIR QUE AL METER LAS INTERACCIONES TODA LA COLUMNA SEA 0. HAGO UN FILTRO PREVIO
      sd.regulators = apply(des.mat2[,-1,drop=FALSE], 2, sd, na.rm=TRUE)
      regulators0 = names(sd.regulators[sd.regulators==0])
      if (length(regulators0)>0) des.mat2 = des.mat2[, setdiff(colnames(des.mat2), regulators0), drop=FALSE]
      
    } else{  ### WITHOUT INTERACTIONS
      
      des.mat2 = cbind(t(GeneExpression[gene,]), des.mat2)
      colnames(des.mat2)[1] = "response"
      
    }
    
  }
  return(des.mat2)
  
}

# ElasticNet variable selection -------------------------------------------

ElasticNet = function (family2, des.mat2, epsilon, elasticnet) {
  
  if (!is.null(elasticnet) && (ncol(des.mat2) > 2)) {  # Variable selection with ElasticNet
    
    if (length(elasticnet) ==1) { # Variable selection with fixed alpha
      
      if (nrow(des.mat2) > 200) { mynfolds =  ceiling(nrow(des.mat2)/20) } else { mynfolds = nrow(des.mat2) }  # for CV
      
      cvEN = glmnet::cv.glmnet(x = as.matrix(des.mat2[,-1]),
                       y = des.mat2[,1], 
                       nfolds = mynfolds, alpha = elasticnet, standardize = FALSE, thres = epsilon,
                       family = family2, grouped = FALSE)
      
      myS = cvEN$lambda.min  # optimum penalization parameter
      y.fitted = predict(cvEN, s = myS, newx = as.matrix(des.mat2[,-1]))
      
      if (length(myS) > 1) {  # more than 1 values for lambda
        myCVerror = cvEN$cvm[sapply(myS, function (i) which(cvEN$lambda == i))]
        myS = myS[which.min(myCVerror)]
      }
      
      if (length(myS) == 1) {
        x = des.mat2[,-1]
        sel = colnames(x)[which(glmnet::coef.glmnet(cvEN, s = myS)[-1,1] != 0)] # selected coefficients without intercept
        selcoef = as.data.frame(as.matrix(glmnet::coef.glmnet(cvEN, s = myS)[which(glmnet::coef.glmnet(cvEN, s = myS)[,1] != 0),,drop=FALSE]))
        mycoef = c(colnames(des.mat2)[1], sel)
        
        removedCoefs = setdiff(colnames(des.mat2), mycoef)
        des.mat2 = des.mat2[, mycoef, drop = FALSE]
      }
      
      m = modelcharac(cvEN, myS, des.mat2[,1], y.fitted)
      isModel = TRUE
      
    } else { 
      #alpha parameter optimization
      
      alphas = elasticnet
      if (nrow(des.mat2) > 200) { mynfolds =  ceiling(nrow(des.mat2)/20) } else { mynfolds = nrow(des.mat2) }  # for CV
      
      #cv for lambda optimization for each alpha of the grid
      
      cvs <- lapply(alphas, function(x) {
        glmnet::cv.glmnet(x = as.matrix(des.mat2[,-1]), y = des.mat2[,1], nfolds = mynfolds,
                  alpha = x, standardize = FALSE, thres = epsilon,
                  family = family2, grouped = FALSE)})
      
      lambdamin = 10000000
      cvupmin = 1000000000
      cvEN = NULL
      elasticnet = NULL
      for ( i in 1:length(cvs)){
        
        if(cvs[[i]]$cvup[which(cvs[[i]]$lambda == cvs[[i]]$lambda.min)] < cvupmin){
          lambdamin = cvs[[i]]$lambda.min
          cvupmin = cvs[[i]]$cvup[which(cvs[[i]]$lambda == cvs[[i]]$lambda.min)]
          cvEN = cvs[[i]]
          elasticnet = alphas[i]
        }
      }
      
      # optimum penalization parameters
      y.fitted = predict(cvEN, s = lambdamin, newx = as.matrix(des.mat2[,-1]))
      
      if (length(lambdamin) > 1) {  # more than 1 values for lambda
        myCVerror = cvEN$cvm[sapply(lambdamin, function (i) which(cvEN$lambda == i))]
        lambdamin = lambdamin[which.min(myCVerror)]
      }
      
      if (length(lambdamin) == 1) {
        x = des.mat2[,-1]
        sel = colnames(x)[which(glmnet::coef.glmnet(cvEN, s = lambdamin)[-1,1] != 0)] # selected coefficients without intercept
        selcoef = as.data.frame(as.matrix(glmnet::coef.glmnet(cvEN, s = lambdamin)[which(glmnet::coef.glmnet(cvEN, s = lambdamin)[,1] != 0),,drop=FALSE]))
        mycoef = c(colnames(des.mat2)[1], sel)
        
        removedCoefs = setdiff(colnames(des.mat2), mycoef)
        des.mat2 = des.mat2[, mycoef, drop = FALSE]
      }
      
      m = modelcharac(cvEN, lambdamin, des.mat2[,1], y.fitted)
      isModel = TRUE
    }
    
  } else { 
    
    if(is.null(elasticnet)&& ncol(des.mat2) > 2){
      #alpha parameter optimization
      
      alphas = seq(0,1,0.1)
      if (nrow(des.mat2) > 200) { mynfolds =  ceiling(nrow(des.mat2)/20) } else { mynfolds = nrow(des.mat2) }  # for CV
      
      #cv for lambda optimization for each alpha of the grid
      
      cvs <- lapply(alphas, function(x) {
        glmnet::cv.glmnet(x = as.matrix(des.mat2[,-1]), y = des.mat2[,1], nfolds = mynfolds,
                  alpha = x, standardize = FALSE, thres = epsilon,
                  family = family2, grouped = FALSE)})
      
      lambdamin = 10000000
      cvupmin = 1000000000
      cvEN = NULL
      elasticnet = NULL
      for ( i in 1:length(cvs)){
        
        if(cvs[[i]]$cvup[which(cvs[[i]]$lambda == cvs[[i]]$lambda.min)] < cvupmin){
          lambdamin = cvs[[i]]$lambda.min
          cvupmin = cvs[[i]]$cvup[which(cvs[[i]]$lambda == cvs[[i]]$lambda.min)]
          cvEN = cvs[[i]]
          elasticnet = alphas[i]
        }
      }
      
      # optimum penalization parameters
      y.fitted = predict(cvEN, s = lambdamin, newx = as.matrix(des.mat2[,-1]))
      
      if (length(lambdamin) > 1) {  # more than 1 values for lambda
        myCVerror = cvEN$cvm[sapply(lambdamin, function (i) which(cvEN$lambda == i))]
        lambdamin = lambdamin[which.min(myCVerror)]
      }
      
      if (length(lambdamin) == 1) {
        x = des.mat2[,-1]
        sel = colnames(x)[which(glmnet::coef.glmnet(cvEN, s = lambdamin)[-1,1] != 0)] # selected coefficients without intercept
        selcoef = as.data.frame(as.matrix(glmnet::coef.glmnet(cvEN, s = lambdamin)[which(glmnet::coef.glmnet(cvEN, s = lambdamin)[,1] != 0),,drop=FALSE]))
        mycoef = c(colnames(des.mat2)[1], sel)
        
        removedCoefs = setdiff(colnames(des.mat2), mycoef)
        des.mat2 = des.mat2[, mycoef, drop = FALSE]
      }
      
      m = modelcharac(cvEN, lambdamin, des.mat2[,1], y.fitted)
      isModel = TRUE
      
    }
    else{ removedCoefs = NULL; selcoef=NULL;y.fitted = NULL; m = list('AICc'=NULL,'R.squared'=NULL,'cvRMSE'=NULL ); isModel = NULL}}
  
  
  return(list("des.mat2" = des.mat2, "removedCoefs" = removedCoefs, "coefficients" = selcoef, 'fitted.values' = y.fitted, 'm' = m, 'isModel'=isModel, 'elasticnet'=elasticnet))
}

modelcharac = function(fitted.glm,s, y, y.fitted){
  n = fitted.glm$glmnet.fit$nobs
  
  R2 = round(fitted.glm$glmnet.fit$dev.ratio[which(fitted.glm$glmnet.fit$lambda == s)],6)
  RMSE = round(sqrt(sum((y-y.fitted)^2)/n),6)
  cvRMSE= abs(round(sqrt(sum((y-y.fitted)^2)/n)/mean(y),6))
  
  return(list('R.squared'=R2, 'RMSE' = RMSE ,'cvRMSE'=cvRMSE))
}

# Functions needed for PLS


# Adding interations with regulators -------------------------------------

RegulatorsInteractionsPLS2 = function (interactions.reg, reguValues, des.mat) {
  
  # Adding regulators
  if (is.null(des.mat)) {
    des.mat2 = data.frame(reguValues, check.names = FALSE)
  } else {
    des.mat2 = data.frame(des.mat, reguValues, check.names = FALSE)
    ### WITH INTERACTIONS with regulators
    if (interactions.reg > 0) {
      expcond = colnames(des.mat)
      
      if (interactions.reg == 1) {  # Max order of interaction = 2 & Interactions with cont.var not allowed
        if (length(grep(":", expcond, fixed = TRUE)))  expcond = expcond[-grep(":", expcond, fixed = TRUE)]
      }
      
      if (interactions.reg == 2) { # Max order of interaction = 2
        if (length(grep(":", expcond, fixed = TRUE)))  expcond = expcond[-grep(":", expcond, fixed = TRUE)]
      }
      
      if(ncol(reguValues)<=5000){
        
        fff = paste0("~ ",
                     paste(sapply(colnames(reguValues),
                                  function (x) paste(expcond, sprintf("`%s`", x), sep = ":")),
                           collapse = "+"))
        
        fff = as.formula(fff)
        inter.var = model.matrix(fff, des.mat2)[,-1, drop = FALSE]
        colnames(inter.var) = strsplit(as.character(fff),"\\s*\\+\\s*")[-1][[1]]
        des.mat2 = cbind(des.mat2, inter.var)
        
        colnames(des.mat2) = gsub("\`", "", colnames(des.mat2))
        
      } else{
        j = ceiling(ncol(reguValues)/250)
        res.mat = des.mat2
        
        for (k in 1:j){
          
          cols_start = (250 * (k - 1)) + 1
          cols_end = min(250 * k, ncol(reguValues))
          
          fff = paste0("~ ",
                       paste(sapply(colnames(reguValues[, cols_start:cols_end, drop = FALSE]),
                                    function(x) paste(expcond, sprintf("`%s`", x), sep = ":")),
                             collapse = "+"))
          
          fff = as.formula(fff)
          inter.var = model.matrix(fff, des.mat2)[,-1, drop = FALSE]
          colnames(inter.var) = strsplit(as.character(fff), "\\s*\\+\\s*")[-1][[1]]
          res.mat = cbind(res.mat, inter.var)
        }
        
        des.mat2 = res.mat
        
      }
      
      ## PUEDE OCURRIR QUE AL METER LAS INTERACCIONES TODA LA COLUMNA SEA 0. HAGO UN FILTRO PREVIO
      sd.regulators = apply(des.mat2, 2, sd, na.rm=TRUE)
      regulators0 = names(sd.regulators[sd.regulators==0])
      if (length(regulators0)>0) des.mat2 = des.mat2[, setdiff(colnames(des.mat2), regulators0), drop=FALSE]
      
    }
    
  }
  return(des.mat2)
  
}


# Scaling by blocks -------------------------------

Scaling.type= function(reguVal, regu, scaletype){
  res.mat <- NULL
  regu = Filter(function(x) !is.null(x) && length(x) > 0, regu)
  for (ov in names(regu)){
    des.mat2 = as.data.frame(reguVal[, regu[[ov]],drop =FALSE])
    
    if(scaletype=='pareto'){
      #Escalado tipo pareto
      des.mat2 = des.mat2 / (ncol(des.mat2)^(1/4))
    }
    if(scaletype=='block'){
      #Cada bloque de variables tiene el mismo peso total en el modelo
      des.mat2 = des.mat2 / (ncol(des.mat2)^(1/2))
    }
    
    if (is.null(res.mat)) {
      res.mat <- des.mat2
    } else {
      res.mat <- cbind(res.mat, des.mat2)
    }
    
  }
  
  rownames(res.mat) = rownames(reguVal)
  
  return(res.mat)
}

# p-values for PLS models ------------------------------------------------

p.coef<-function(pls,R, datospls){
  #pls: modelo PLS generado con la libreria ropls
  #R: número de veces a repetir la prueba
  #datospls: matriz de datos utilizada para crear el modelo
  k=pls@summaryDF$pre
  coefmod<-pls@coefficientMN
  Y=datospls[,1]
  a<-NULL
  for (i in 1:R){
    Yperm=sample(Y, replace=FALSE)
    plsda.opls<-ropls::opls(datospls[,-1], scale(Yperm), scaleC='none', predI=k,
                     info.txtC='none', fig.pdfC='none', crossvalI=1,permI = 0)
    a<-cbind(a,plsda.opls@coefficientMN)
  }
  p.coefs<-matrix(2, nrow=nrow(a), ncol=1)
  for(i in 1:nrow(a)){
    coefs<-a[i,]
    pvalor = ifelse(coefmod[i] > 0, (sum(coefs > coefmod[i]) + sum(coefs < -coefmod[i]) )/ R, (sum(coefs < coefmod[i]) + sum(coefs > -coefmod[i])) / R)
    p.coefs[i,1]<-pvalor
  }
  rownames(p.coefs)<-rownames(coefmod)
  return(p.coefs)
}

p.valuejack<-function(pls, datospls,alfa){
  
  #pls: modelo PLS generado con la libreria ropls
  #datospls: matriz de datos utilizada para crear el modelo
  
  #Caso de Leave-one-out
  k=pls@summaryDF$pre
  coefmod=pls@coefficientMN
  a=NULL
  
  for (i in 1: nrow(datospls)) {
    
    pls.opls=suppressWarnings(ropls::opls(datospls[-i,-1], scale(datospls[-i,1]), scaleC='none', predI=k,
                                   info.txtC='none', fig.pdfC='none', crossvalI=1, permI = 0))
    
    if(nrow(pls.opls@coefficientMN)<nrow(coefmod)){
      exclvar=setdiff(rownames(coefmod),rownames(pls.opls@coefficientMN))
      b<-matrix(0,ncol=1,nrow = length(exclvar))
      rownames(b)<-exclvar
      pls.opls@coefficientMN=rbind(pls.opls@coefficientMN,b)
      
    }
    order <- match(rownames(coefmod), rownames(pls.opls@coefficientMN))
    plscoefficientMN <- pls.opls@coefficientMN[order, ]
    a=cbind(a,plscoefficientMN)
  }
  
  est = (a-matrix(rep(coefmod, ncol(a)), ncol = ncol(a)))^2
  estjack = sqrt(((nrow(datospls)-1)/nrow(datospls)) * rowSums(est))
  
  pvalor = 2*pt(abs(coefmod/estjack), df = nrow(datospls)-1, lower.tail = FALSE)
  #calculo de los intervalos de confianza para ver que todo va bien, aún y todo podemos ahorrar esto no es necesario en nuestro caso
  #res =cbind(pvalor, paste('(',coefmod- qt(1 - (alfa)/ 2, df = nrow(datospls)-1)* estjack,'-', coefmod+ qt(1 - (alfa)/ 2, df = nrow(datospls)-1)* estjack,')' ))
  #colnames(res) = c('pvalue', 'CI')
  colnames(pvalor) = c('pvalue')
  return(pvalor)
}

# p-values for PLS2 models ------------------------------------------------

p.coef.pls2<-function(pls,R, datospls, Y){
  #pls: modelo PLS generado con la libreria ropls
  #R: número de veces a repetir la prueba
  #datospls: matriz de datos utilizada para crear el modelo
  k=pls@summaryDF$pre
  coefmod<-pls@coefficientMN
  a<-NULL
  for (i in 1:R){
    rows_ord = sample(nrow(Y), replace = FALSE)
    Yperm=Y[rows_ord,]
    pls.opls<-ropls::opls(datospls, Yperm, scaleC='none', predI=k,
                   info.txtC='none', fig.pdfC='none', crossvalI=1,permI = 0)
    a<-cbind(a,pls.opls@coefficientMN)
  }
  p.coefs<-matrix(2, nrow=nrow(a), ncol=ncol(Y))
  for (i in 1:ncol(Y)) {
    b = a[, seq(i, ncol(Y)*R, ncol(Y))]
    
    for(j in 1:nrow(b)){
      coefs<-b[j,]
      pvalor = ifelse(coefmod[j,i] > 0, (sum(coefs > coefmod[j,i]) + sum(coefs < -coefmod[j,i]) )/ R, (sum(coefs < coefmod[j,i]) + sum(coefs > -coefmod[j,i])) / R)
      p.coefs[j,i]<-pvalor
    }
  }
  rownames(p.coefs)<-rownames(coefmod)
  colnames(p.coefs)<-colnames(coefmod)
  
  return(p.coefs)
}

p.valuejack.pls2<-function(pls, datospls, Y,alfa){
  
  #pls: modelo PLS2 generado con la libreria ropls
  #datospls: matriz de datos utilizada para crear el modelo
  
  #Caso de Leave-one-out
  k=pls@summaryDF$pre
  coefmod=pls@coefficientMN
  a=NULL
  pvalores = data.frame()
  pb <- txtProgressBar(min = 0, max = nrow(datospls), style = 3)
  for (i in 1: nrow(datospls)) {
    setTxtProgressBar(pb, value = i)
    pls.opls=suppressWarnings(ropls::opls(datospls[-i, , drop =FALSE], Y[-i,], scaleC='none', predI=k,
                                   info.txtC='none', fig.pdfC='none', crossvalI=1, permI = 0))
    
    if(nrow(pls.opls@coefficientMN)<nrow(coefmod)){
      exclvar=setdiff(rownames(coefmod),rownames(pls.opls@coefficientMN))
      b<-matrix(0,ncol=ncol(pls.opls@coefficientMN),nrow = length(exclvar))
      rownames(b)<-exclvar
      pls.opls@coefficientMN=rbind(pls.opls@coefficientMN,b)
      
    }
    order <- match(rownames(coefmod), rownames(pls.opls@coefficientMN))
    plscoefficientMN <- pls.opls@coefficientMN[order, ]
    rm(pls.opls);gc()
    a=cbind(a,plscoefficientMN)
  }
  close(pb)
  for (i in 1:ncol(pls@coefficientMN)) {
    
    b = a[, seq(i, ncol(pls@coefficientMN)*nrow(datospls), ncol(pls@coefficientMN))]
    
    est = (b-matrix(rep(coefmod[,i], ncol(b)), ncol = ncol(b)))^2
    estjack = sqrt(((nrow(datospls)-1)/nrow(datospls)) * rowSums(est))
    
    pvalor = as.data.frame(2*pt(abs(coefmod[,i]/estjack), df = nrow(datospls), lower.tail = FALSE))
    colnames(pvalor) = colnames(coefmod)[i]
    
    pvalores = c(pvalores,pvalor)
  }
  pvalores =as.data.frame(pvalores)
  rownames(pvalores) = rownames(a)
  
  return(pvalores)
}

