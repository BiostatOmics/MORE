
##  Verify omic type ####

isBin <-function(x){
  if(length(unique(x[!is.na(x[, 1]), 1]))==2 && length(unique(x[1,!is.na(x[1,]),drop=TRUE]))==2){
    return(1)
  }
  else{
    return(0)
  }
}

# Apply for variable type prediction
if(is.null(omicType)){
  #Create internally omicType vector
  omicType = sapply(regulatoryData, function(x) isBin(x)) #We assume that all regulators are of the same type in the same omic
  
  cat('Considering that we codify by 1 binary omics and by 0 numeric omics, we consider the following nature of the omics: \n')
  print(omicType)
  
  cat('Please if it is incorrect stop the generation of the models and introduce omicType manually \n')
}


##  Verify clinic data type ####

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

# Apply for variable type prediction

if(!is.null(clinic)){
  if(is.null(clinicType)){
    clinicType = sapply(clinic, function(x) isBinclinic(x))
    
    cat('Considering that we codify by 1 binary/categorical features and by 0 numeric features, we consider the following nature of the clinical features:\n')
    print(clinicType)
    
    cat('Please if it is incorrect stop the generation of the models and introduce clinicType manually \n')
  }
}

## targetData, regulatoryData, condition dimensions check ####

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

## target, regulatoryData, condition names ####

# Compulsory to have the same names

if(is.null(condition)){
  nameproblem<-!all(sapply(regulatoryData, function(x) length(intersect(colnames(x),colnames(targetData))==ncol(targetData))))
  if(nameproblem){
    stop('ERROR: targetData and regulatoryData samples have not same names.\n')
  }
} else{
  nameproblem<-!all(c(sapply(regulatoryData, function(x) length(intersect(colnames(x),colnames(targetData)))==ncol(targetData)), length(intersect(rownames(condition),colnames(targetData)))==ncol(targetData)))
  if(nameproblem){
    stop('ERROR: targetData, condition and regulatoryData samples have not same names. We assume that they are ordered.\n')
  }
}

## Label restrictions ####

# Only necessary if the Shiny app does not show the console output.

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

## Missing values ####

## Removing regulators with more than percNA of NAs and keeping track
missing_rows <- function(regOmic,percNA) {
  highNA = apply(regOmic, 1, function(x) sum(is.na(x))/length(x)) > percNA
  return(rownames(regOmic)[highNA])
}

myregNA = lapply(regulatoryData, function(x) missing_rows(x,0))
cat("Number of regulators with missing values:\n")
print(sapply(myregNA, length))
cat("\n")

#Filter observations with two many missing values
missing_cols <- function(regOmic,percNA) {
  highNA = apply(regOmic, 2, function(x) sum(is.na(x))/length(x)) > percNA
  return(colnames(regOmic)[highNA])
}

myobsNA = lapply(regulatoryData, function(x) missing_cols(x,0))
myobsNA = unique(unlist(myobsNA))
cat("Number of observations with missing values:", length(myobsNA))
cat("\n")

