
#Apply lowVariation filter

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

tmp = MORE::LowVariationRegu(minVariation, regulatoryData, Group, associations, AlltargetFs, omicType, clinicType)
regulatoryData = tmp[["data.omics"]]
associations = tmp[["associations"]]
myregLV = tmp[["myregLV"]]
rm("tmp"); gc()

if(all(sapply(regulatoryData, function(x)nrow(x)==0))) stop("ERROR: No regulators left after LowVariation filter. Consider being less restrictive.")
