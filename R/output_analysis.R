
#########################################
#### Functions for output analysis
#########################################

## By Sonia, Monica, Maider
## 12-March-2024

# Function to obtain all significant pairs gene-regulator per omic --------

# For only 1 gene
GetPairs1GeneAllReg = function (gene, output) {
  
  if(output$arguments$method=='glm'|| output$arguments$method=='isgl'){
    
    reguSignif = output$ResultsPerGene[[gene]]$relevantRegulators
    
    if (is.null(reguSignif)) {  # NO significant regulators
      return (NULL)
      
    } else {  # Significant regulators
      
      reguSignif = output$ResultsPerGene[[gene]]$allRegulators[reguSignif,]
      reguSignif = reguSignif[,c("gene", "regulator", "omic", "area", "filter")]
      return (reguSignif)
    }
    
  }
  if(output$arguments$method=='pls1' || output$arguments$method=='pls2'){
    
    reguSignif = output$ResultsPerGene[[gene]]$significantRegulators
    
    if (is.null(reguSignif)) {  # NO significant regulators
      return (NULL)
      
    } else {  # Significant regulators
      
      reguSignif = output$ResultsPerGene[[gene]]$allRegulators[reguSignif,]
      reguSignif = reguSignif[,c("gene", "regulator", "omic", "area", "filter")]
      return (reguSignif)
    }
    
  }
  
}


# For all genes
GetPairsGeneRegulator = function (genes = NULL, output) {
  
  if (is.null(genes)) genes = rownames(output$GlobalSummary$ReguPerGene)
  
  myresults = do.call("rbind", lapply(genes, GetPairs1GeneAllReg, output))
  
  #   colnames(myresults) = c("gene", "regulator", "omic", "area")
  return(myresults)
}

#' RegulationPerCondition
#'
#' \code{RegulationPerCondition} Function to be applied to \link{more} main function output.
#' 
#' @param output Output object of MORE main function.
#' 
#' @return Summary table containing all the relevant/significant regulators. Moreover, it provides the regression coefficient that relates the gene and the regulator for each experimental condition after testing if this coefficient is relevant/significant or not.
#'
#' @export


RegulationPerCondition = function(output){
  # output: results of the getGLM/getPLS function.
  method = output$arguments$method
  #Add a progressbar
  pb <- txtProgressBar(min = 0, max = length(rownames(output$GlobalSummary$ReguPerGene)), style = 3)
  if(method =='glm'|| method=='isgl'){
    design = output$arguments$finaldesign
    Group = output$arguments$groups
    
    # Creo a partir de la funcion que ya estaba hecha (linea 1657) la tabla y le anyado los huecos en blanco y cambio el nombre a "representative".
    genes = rownames(output$GlobalSummary$ReguPerGene)
    myresults = do.call("rbind", lapply(genes, GetPairs1GeneAllReg, output))
    colnames(myresults) = c(colnames(myresults)[1:4], "representative")
    myresults[myresults[, "representative"] == "Model", "representative"] = ""
    if (is.null(design)){
      
      # Anyado la columna de coeficientes.
      coeffs = matrix(1, nrow(myresults), 1)
      colnames(coeffs) = "coefficients"
      rownames(coeffs) = rownames(myresults)
      myresults = cbind(myresults, coeffs)
      myresults[grep("_N", myresults[, "representative"]), "coefficients"] = -1  # Para cambiar el signo si pertenece al grupo de correlacionados negativamente
      
      for(k in unique(myresults[,"gene"])){
        setTxtProgressBar(pb, value = which(names(output$ResultsPerGene)==k))
        # Posicion y reguladores que son representantes.
        counts = grep("_R", myresults[myresults[,"gene"] == k, "representative"]) # positions of representatives of mc
        representatives = myresults[myresults[,"gene"] == k, "regulator"][counts]      # Devuelve el nombre real de los reguladores representantes
        omic.representative = myresults[myresults[,"gene"] == k, c("regulator", "representative")][counts,]   # Columna Regulator y Representative
        
        # Necesito el if, si no da error. En caso de entrar, elimino las coletillas para que sea mas sencillo buscar y asignar el representante
        if(length(representatives) != 0){
          norow.nulls = which(myresults[myresults[,"gene"] == k, "representative"] != "")
          myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_P", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])
          myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_N", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])
          myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_R", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])
          
          for(i in 1:length(representatives)){
            # Aquellos que se llamen igual "omica_mc(numero)", se les asignara el representante
            reg.rep = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == representatives[i], "representative"]
            myresults[myresults[,"gene"] == k & myresults[,"representative"] == reg.rep, "representative"] = representatives[i]
          }
          
          # Reguladores significativos del GLM. Pongo gsub() porque haciendo pruebas he visto que hay reguladores que se nombran `nombre regulador`.
          # Las comitas haran que no pueda encontrar el regulador
          # en la tabla. Sin embargo, creo sign.glm para manterner las comitas y poder acceder a la tabla de coeficientes
          significatives = gsub("`", "", names(output$ResultsPerGene[[k]]$coefficients[2:nrow(output$ResultsPerGene[[k]]$coefficients), 1]))
          sign.glm = names(output$ResultsPerGene[[k]]$coefficients[2:nrow(output$ResultsPerGene[[k]]$coefficients), 1])
          
          for(i in 1:length(significatives)){
            if(any(significatives[i] == omic.representative[,2])){
              # index.regul: para saber que regulador es el representante y asi todos los que tengan su nombre en la columna "representative" tendran su coeficiente del modelo GLM.
              index.regul = rownames(omic.representative)[which(omic.representative[,2] == significatives[i])]
              PN = myresults[myresults[,"gene"] == k & myresults[,"representative"] == index.regul, "coefficients"]                        # Sera 1 o -1, segun tenga "_P" o "_N"
              myresults[myresults[,"gene"] == k & myresults[,"representative"] == index.regul, "coefficients"] = PN*output$ResultsPerGene[[k]]$coefficients[sign.glm[i], 1]                                                    # Tendra signo de la tabla si es "_P" y signo opuesto si es "_N".
            } else {
              # En caso de no pertenecer a un grupo de reguladores correlacionados, cogera su coeficiente de la tabla y lo asignara a la posicion correspondiente
              myresults[myresults[,"gene"] == k & myresults[,"regulator"] == significatives[i], "coefficients"] = output$ResultsPerGene[[k]]$coefficients[sign.glm[i], 1]
            }
          }
          
        } else {
          # Si no presenta grupo de reguladores correlacionados, simplemente sacara los coeficientes de la tabla "coefficients"
          myresults[myresults[,"gene"] == k, "coefficients"] = output$ResultsPerGene[[k]]$coefficients[2:nrow(output$ResultsPerGene[[k]]$coefficients), 1]
        }
      }
      
      
    } else {
      
      # Anyado las columnas de las condiciones experimentales. Pongo "Group" porque al hacer model.matrix() siempre coloca "Group" y lo que se almacena en el objeto Group
      index = unique(Group)
      names.groups = paste("Group", index, sep = "_")
      conditions = matrix(0, nrow(myresults), length(names.groups))
      colnames(conditions) = names.groups
      rownames(conditions) = rownames(myresults)
      myresults = cbind(myresults, conditions)
      
      for(k in unique(myresults[,"gene"])){
        setTxtProgressBar(pb, value = which(names(output$ResultsPerGene)==k))
        significant.regulators = output$ResultsPerGene[[k]]$relevantRegulators                    # Reguladores significativos.
        if(method =='glm'){
          model.variables = gsub("`", "", rownames(output$ResultsPerGene[[k]]$coefficients))[-1]       # Reguladores e interacciones en el modelo.
          kc = 2
        } else{
          model.variables = gsub("`", "", rownames(output$ResultsPerGene[[k]]$coefficients))       # Reguladores e interacciones en el modelo.
          kc = 1
        }
        
        # Cojo las interacciones y creo objetos que contengan los reguladores que aparecen con interaccion, solas o ambas.
        interactions.model = gsub("`", "", rownames(output$ResultsPerGene[[k]]$coefficients)[grep(":", rownames(output$ResultsPerGene[[k]]$coefficients))])
        
        inter.variables = unlist(strsplit(interactions.model, ":", fixed = TRUE))
        if(is.null(inter.variables)){
          inter.variables = NULL                                                                            # No hay interacciones.
        } else {
          inter.variables = inter.variables[seq(2, length(inter.variables), by = 2)]                        # Reguladores que presentan interseccion con algun grupo.
        }
        
        variables.only = setdiff(setdiff(model.variables, interactions.model), inter.variables)             # Reguladores solos en el modelo, sin interacciones.
        if(length(grep("Group", variables.only)) != 0){                                                     # No puedo hacer la interseccion con las variables significativas porque me cargo tambien omic_mc: hay que eliminar Group si esta.
          variables.only = variables.only[-grep("Group", variables.only)]
        }
        
        variables.inter.only = intersect(inter.variables, model.variables)                                  # Reguladores con interaccion y solas.
        variables.inter = setdiff(inter.variables, model.variables)                                         # Reguladores con solo interaccion (no aparecen solas en el modelo).
        
        for(j in kc:nrow(output$ResultsPerGene[[k]]$coefficients)){
          regul = unlist(strsplit(gsub("`", "", rownames(output$ResultsPerGene[[k]]$coefficients)[j]), ":"))
          
          # Evaluo en que conjunto se encuentra el regulador correspondiente y segun eso asigno el coeficiente o sumo el nuevo coeficiente a lo que ya habia en esa posicion.
          if(any(regul %in% variables.only)){
            if(any(regul %in% significant.regulators)){
              myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] = output$ResultsPerGene[[k]]$coefficients[j,]
            } else {
              myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul, c(names.groups)] = output$ResultsPerGene[[k]]$coefficients[j,]
            }
          }
          
          if(any(regul %in% variables.inter)){
            if(any(regul %in% significant.regulators)){
              myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul[2], regul[1]] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul[2], regul[1]] + output$ResultsPerGene[[k]]$coefficients[j,]
            } else {
              myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul[2], regul[1]] = myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul[2], regul[1]] + output$ResultsPerGene[[k]]$coefficients[j,]
            }
          }
          
          if(any(regul %in% variables.inter.only)){
            if(any(regul %in% significant.regulators)){
              if(length(regul) == 1){
                myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] + output$ResultsPerGene[[k]]$coefficients[j,]
              } else {
                myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul[2], regul[1]] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul[2], regul[1]] + output$ResultsPerGene[[k]]$coefficients[j,]
              }
            } else {
              if(length(regul) == 1){
                myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul, c(names.groups)] = myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul, c(names.groups)] + output$ResultsPerGene[[k]]$coefficients[j,]
              } else {
                myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul[2], regul[1]] = myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul[2], regul[1]] + output$ResultsPerGene[[k]]$coefficients[j,]
              }
            }
          }
        }
        
        # Veo si hay representantes, en caso de haberlos asignara la misma fila del representante a los reguladores que acaben en "_P" y el opuesto a los que acaban en "_N".
        countsR = grep("_R", myresults[myresults[,"gene"] == k, "representative"])
        
        if(length(countsR) != 0){
          countsR = myresults[myresults[,"gene"] == k, 5:ncol(myresults)][countsR,]
          
          # Para los correlacionados positivamente: mete la misma fila de coeficientes del representante.
          countsP = countsR
          countsP[,"representative"] = sub("_R", "", countsP[,"representative"])
          countsP[,"representative"] = paste(countsP[,"representative"], "_P", sep = "")
          
          for(l in 1:nrow(countsP)){
            myresults[myresults[,"gene"] == k & myresults[,"representative"] == countsP[l,"representative"], 6:ncol(myresults)] = countsP[l,2:ncol(countsP)]
          }
          
          # Para los correlacionados negativamente: mete la fila opuesta de coeficientes del representante.
          countsN = countsR
          countsN[,"representative"] = sub("_R", "", countsN[,"representative"])
          countsN[,"representative"] = paste(countsN[,"representative"], "_N", sep = "")
          
          for(l in 1:nrow(countsN)){
            myresults[myresults[,"gene"] == k & myresults[,"representative"] == countsN[l,"representative"], 6:ncol(myresults)] = -countsN[l,2:ncol(countsN)]
          }
          
          counts = grep("_R", myresults[myresults[,"gene"] == k, "representative"])
          representatives = myresults[myresults[,"gene"] == k, "regulator"][counts]                                    # Devuelve el nombre real de los reguladores representantes
          omic.representative = myresults[myresults[,"gene"] == k, c("regulator", "representative")][counts,]          # Columna Regulator y Representative
          
          # Necesito el if, sino da error. En caso de entrar, elimino las coletillas para que sea mas sencillo buscar y asignar el representante.
          if(length(representatives) != 0){
            norow.nulls = which(myresults[myresults[,"gene"] == k, "representative"] != "")
            myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_P", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])
            myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_N", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])
            myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_R", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])
            
            for(i in 1:length(representatives)){
              # Aquellos que se llamen igual "omica_mc(numero)", se les asignara el representante.
              reg.rep = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == representatives[i], "representative"]
              myresults[myresults[,"gene"] == k & myresults[,"representative"] == reg.rep, "representative"] = representatives[i]
            }
          }
        }
      }
    }
    myresults[,6:ncol(myresults)] = signif(myresults[,6:ncol(myresults)], digits = 4) # Para que no salgan los numeros en diferentes notaciones
    
  }
  if(method=='pls1' || method=='pls2'){
    design = output$arguments$finaldesign
    Group = output$arguments$groups
    
    # Creo a partir de la funcion que ya estaba hecha (linea 1657) la tabla y le anyado los huecos en blanco y cambio el nombre a "representative".
    genes = rownames(output$GlobalSummary$ReguPerGene)
    myresults = do.call("rbind", lapply(genes, GetPairs1GeneAllReg, output))
    colnames(myresults) = c(colnames(myresults)[1:4], "representative")
    myresults[myresults[, "representative"] == "Model", "representative"] = ""
    
    if (is.null(design)){
      
      # Anyado la columna de coeficientes.
      coeffs = matrix(1, nrow(myresults), 1)
      colnames(coeffs) = "coefficients"
      rownames(coeffs) = rownames(myresults)
      myresults = cbind(myresults, coeffs)
      myresults[grep("_N", myresults[, "representative"]), "coefficients"] = -1  # Para cambiar el signo si pertenece al grupo de correlacionados negativamente
      
      for(k in unique(myresults[,"gene"])){
        setTxtProgressBar(pb, value = which(names(output$ResultsPerGene)==k))
        # Posicion y reguladores que son representantes.
        counts = grep("_R", myresults[myresults[,"gene"] == k, "representative"]) # positions of representatives of mc
        representatives = myresults[myresults[,"gene"] == k, "regulator"][counts]      # Devuelve el nombre real de los reguladores representantes
        omic.representative = myresults[myresults[,"gene"] == k, c("regulator", "representative")][counts,]   # Columna Regulator y Representative
        
        # Necesito el if, si no da error. En caso de entrar, elimino las coletillas para que sea mas sencillo buscar y asignar el representante
        if(length(representatives) != 0){
          norow.nulls = which(myresults[myresults[,"gene"] == k, "representative"] != "")
          myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_P", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])
          myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_N", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])
          myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_R", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])
          
          for(i in 1:length(representatives)){
            # Aquellos que se llamen igual "omica_mc(numero)", se les asignara el representante
            reg.rep = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == representatives[i], "representative"]
            myresults[myresults[,"gene"] == k & myresults[,"representative"] == reg.rep, "representative"] = representatives[i]
          }
          
          # Reguladores significativos del GLM. Pongo gsub() porque haciendo pruebas he visto que hay reguladores que se nombran `nombre regulador`.
          # Las comitas haran que no pueda encontrar el regulador
          # en la tabla. Sin embargo, creo sign.glm para manterner las comitas y poder acceder a la tabla de coeficientes
          significatives = gsub("`", "", names(output$ResultsPerGene[[k]]$coefficients[1:nrow(output$ResultsPerGene[[k]]$coefficients), 1]))
          sign.glm = names(output$ResultsPerGene[[k]]$coefficients[1:nrow(output$ResultsPerGene[[k]]$coefficients), 1])
          
          for(i in 1:length(significatives)){
            if(any(significatives[i] == omic.representative[,2])){
              # index.regul: para saber que regulador es el representante y asi todos los que tengan su nombre en la columna "representative" tendran su coeficiente del modelo GLM.
              index.regul = rownames(omic.representative)[which(omic.representative[,2] == significatives[i])]
              PN = myresults[myresults[,"gene"] == k & myresults[,"representative"] == index.regul, "coefficients"]                        # Sera 1 o -1, segun tenga "_P" o "_N"
              myresults[myresults[,"gene"] == k & myresults[,"representative"] == index.regul, "coefficients"] = PN*output$ResultsPerGene[[k]]$coefficients[sign.glm[i], 1]                                                    # Tendra signo de la tabla si es "_P" y signo opuesto si es "_N".
            } else {
              # En caso de no pertenecer a un grupo de reguladores correlacionados, cogera su coeficiente de la tabla y lo asignara a la posicion correspondiente
              myresults[myresults[,"gene"] == k & myresults[,"regulator"] == significatives[i], "coefficients"] = output$ResultsPerGene[[k]]$coefficients[sign.glm[i], 1]
            }
          }
          
        } else {
          # Si no presenta grupo de reguladores correlacionados, simplemente sacara los coeficientes de la tabla "coefficients"
          myresults[myresults[,"gene"] == k, "coefficients"] = output$ResultsPerGene[[k]]$coefficients[1:nrow(output$ResultsPerGene[[k]]$coefficients), 1]
        }
      }
      
      
    } else {
      
      # Añado las columnas de las condiciones experimentales. Pongo "Group" porque al hacer model.matrix() siempre coloca "Group" y lo que se almacena en el objeto Group
      index = unique(Group)
      names.groups = paste("Group", index, sep = "_")
      conditions = matrix(0, nrow(myresults), length(names.groups))
      colnames(conditions) = names.groups
      rownames(conditions) = rownames(myresults)
      myresults = cbind(myresults, conditions)
      
      for(k in unique(myresults[,"gene"])){
        
        setTxtProgressBar(pb, value = which(names(output$ResultsPerGene)==k))

        significant.regulators = output$ResultsPerGene[[k]]$significantRegulators                    # Reguladores significativos.
        model.variables = gsub("`", "", rownames(output$ResultsPerGene[[k]]$coefficients))           # Reguladores e interacciones en el modelo.
        
        # Cojo las interacciones y creo objetos que contengan los reguladores que aparecen con interaccion, solas o ambas.
        interactions.model = gsub("`", "", rownames(output$ResultsPerGene[[k]]$coefficients)[grep(":", rownames(output$ResultsPerGene[[k]]$coefficients))])
        
        inter.variables = unlist(strsplit(interactions.model, ":", fixed = TRUE))
        if(is.null(inter.variables)){
          inter.variables = NULL                                                                            # No hay interacciones.
        } else {
          inter.variables = inter.variables[seq(2, length(inter.variables), by = 2)]                        # Reguladores que presentan interseccion con algun grupo.
        }
        
        variables.only = setdiff(setdiff(model.variables, interactions.model), inter.variables)             # Reguladores solos en el modelo, sin interacciones.
        if(length(grep("Group_", variables.only)) != 0){                                                     # No puedo hacer la interseccion con las variables significativas porque me cargo tambien omic_mc: hay que eliminar Group si esta.
          variables.only = variables.only[-grep("Group_", variables.only)]
        }
        
        variables.inter.only = intersect(inter.variables, model.variables)                                  # Reguladores con interaccion y solas.
        variables.inter = setdiff(inter.variables, model.variables)                                         # Reguladores con solo interaccion (no aparecen solas en el modelo).
        
        for(j in 1:nrow(output$ResultsPerGene[[k]]$coefficients)){
          regul = unlist(strsplit(gsub("`", "", rownames(output$ResultsPerGene[[k]]$coefficients)[j]), ":"))
          
          groups = regul[grepl(paste(paste0('Group_',names(table(Group))),collapse='|'), regul)]
          regula = setdiff(regul,groups)
          # Evaluo en que conjunto se encuentra el regulador correspondiente y segun eso asigno el coeficiente o sumo el nuevo coeficiente a lo que ya habia en esa posicion.
          if(any(regul %in% variables.only)){
            if(any(regul %in% significant.regulators)){
              myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] = output$ResultsPerGene[[k]]$coefficients[j,1]
            } 
          }
          #TO DO: regul[1] se supone que tendría que ser algo sobre los grupos y no un regulador
          if(any(regul %in% variables.inter)){
            if(any(regul %in% significant.regulators)){
              myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regula, groups] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regula, groups] + output$ResultsPerGene[[k]]$coefficients[j,1]
            } 
          }
          
          if(any(regul %in% variables.inter.only)){
            if(any(regul %in% significant.regulators)){
              if(length(regul) == 1){
                myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] + output$ResultsPerGene[[k]]$coefficients[j,1]
              } else {
                myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regula, groups] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regula, groups] + output$ResultsPerGene[[k]]$coefficients[j,1]
              }
            } 
          }
        }
      }
    }
    myresults[,6:ncol(myresults)] = signif(myresults[,6:ncol(myresults)], digits = 4) # Para que no salgan los numeros en diferentes notaciones
    myresults = myresults[,-5,drop=FALSE]
    
  }
  rownames(myresults)=NULL
  close(pb)
  return(myresults)
}

#' RegulationInCondition
#'
#' \code{RegulationInCondition} Function to be applied to \link{RegulationPerCondition} function output.
#' 
#' @param output_regpcond Output object of \link{RegulationPerCondition} function output.
#' @param cond Biological condition from which the user wants to summarize the information.
#' 
#' @return List with the hub genes, global regulators and the regulators with their coefficients specific to the requested condition.
#'
#' @export

RegulationInCondition <- function (output_regpcond, cond){
  
  #Take group column in RegulationPerCondition
  group_col<-grep(cond,colnames(output_regpcond))
  #Use only the group in which we are interested and remove genes that do not affect that group
  output_regpcond = output_regpcond[c(output_regpcond[,group_col]!=0),c(1,2,3,group_col)]
  
  #Calculate the global regulators
  myreg<-table(output_regpcond[,2])
  #Calculate third quantile
  q3<-quantile(myreg,0.75)
  if(length(myreg[myreg>q3])<10){
    GlobalRegulators = intersect(names(myreg[rev(tail(order(myreg),10))]), names(myreg[myreg>10]) )
  } else{
    GlobalRegulators = intersect(names(myreg[myreg>q3]), names(myreg[myreg>10]) ) 
  }
  
  #Calculate the hub genes
  myhub<-table(output_regpcond[,1])
  #Calculate third quantile
  q3<-quantile(myhub,0.75)
  if(length(myhub[myhub>q3])<10){
    HubGenes = names(myhub[rev(tail(order(myhub),10))])
  } else{
    HubGenes = names(myhub[myhub>q3])
  }
  
  return(list('GlobalRegulators'=GlobalRegulators, 'Hubgenes'=HubGenes, 'RegulationInCondition'=output_regpcond))
  
}


# Plot GLM results --------------------------------------------------------

# Plot 1 gene versus 1 regulator ------------------------------------------

plotGeneRegu = function (x.points, geneValues, reguValues, geneErrorValues, reguErrorValues, col = c(1,2),
                         xlab = "", yylab = c("right", "left"), pch = c(16,17), main = "",
                         numLines = NULL, x.names = NULL, yleftlim, yrightlim) {
  
  # Adjust the axis to include the error value
  if (missing(yrightlim)) {
    if (! missing(geneErrorValues) && ! is.null(geneErrorValues)) {
      yrightlim = range(c(geneValues - geneErrorValues, geneValues + geneErrorValues), na.rm = TRUE)
    } else {
      yrightlim = range(geneValues, na.rm = TRUE)
    }
  }
  
  if (missing(yleftlim)) {
    if (! missing(reguErrorValues) && ! is.null(reguErrorValues)) {
      yleftlim = range(c(reguValues - reguErrorValues, reguValues + reguErrorValues), na.rm = TRUE)
    } else {
      yleftlim = range(reguValues, na.rm = TRUE)
    }
  }
  
  plot.y2(x = x.points, yright = geneValues, yleft = reguValues, yleftlim = yleftlim,
          col = col, xlab = xlab, yylab = yylab, pch = pch, main = main, yrightlim = yrightlim,
          yrightErrorValues = geneErrorValues, yleftErrorValues = reguErrorValues)
  
  if (!is.null(numLines)) {
    for (aa in numLines) {
      abline(v = aa, lty = 2, col = 1)
    }
  }
  
  if (!is.null(x.names)) {
    axis(side=1, at = x.points, labels = x.names, cex.axis = 0.8, las=2)
  }
}


# Plot Y2 -----------------------------------------------------------------

# By Ajay Shah (taken from [R] Plot 2 time series with different y axes (left and right),
# in https://stat.ethz.ch/pipermail/r-help/2004-March/047775.html)

# Modified by: Sonia Tarazona

### PARAMETERS (default):
# x: data to be drawn on X-axis
# yright: data to be drawn on Y right axis
# yleft: data to be drawn on Y left axis
# yrightlim (range(yright, na.rm = TRUE)): ylim for rigth Y-axis
# yleftlim (range(yleft, na.rm = TRUE)): ylim for left Y-axis
# xlab (NULL): Label for X-axis
# yylab (c("","")): Labels for right and left Y-axis
# pch (c(1,2)): Type of symbol for rigth and left data
# col (c(1,2)): Color for rigth and left data
# linky (TRUE): If TRUE, points are connected by lines.
# smooth (0): Friedman's super smoothing
# lwds (1): Line width for smoothed line
# length (10): Number of tick-marks to be drawn on axis
# ...: Other graphical parameters to be added by user (such as main, font, etc.)
###


plot.y2 <- function(x, yright, yleft, yrightlim = range(yright, na.rm = TRUE),
                    yleftlim = range(yleft, na.rm = TRUE),
                    xlim = range(x, na.rm = TRUE),
                    xlab = NULL, yylab = c("",""), lwd = c(2,2),
                    pch = c(1,2), col = c(1,2), type = c("o","o"),
                    linky = TRUE, smooth = 0, bg = c("white","white"),
                    lwds = 1, length = 10, ...,
                    x2 = NULL, yright2 = NULL, yleft2 = NULL, col2 = c(3,4),
                    yrightErrorValues, yleftErrorValues
)
{
  #par(mar = c(5,2,4,2), oma = c(0,3,0,3))
  
  ## Plotting RIGHT axis data
  
  plot(x, yright, axes = FALSE, ylab = "", xlab = xlab, ylim = yrightlim,
       xlim = xlim, pch = pch[1], type = type[1], lwd = lwd[1],
       col = col[1], ...)
  
  axis(4, pretty(yrightlim, length), col = 1, col.axis = 1)
  
  if (is.null(yright2) == FALSE) {
    points(x2, yright2, type = type[1], pch = pch[1], lwd = lwd[1], col = col2[1], ...)
  }
  
  #if (linky) lines(x, yright, col = col[1], ...)
  
  if (smooth != 0) lines(supsmu(x, yright, span = smooth), col = col[1], lwd = lwds, ...)
  
  if(yylab[1]=="") {
    mtext(deparse(substitute(yright)), side = 4, outer = FALSE, line = 2,
          col = col[1], cex = 0.9,...)
  } else {
    mtext(yylab[1], side = 4, outer = FALSE, line = 2, col = col[1], cex = 0.9,...)
  }
  
  # Plot arrows showing standard error
  if (!missing(yrightErrorValues) && ! is.null(yrightErrorValues)) {
    arrows(x, yright - yrightErrorValues, x, yright + yrightErrorValues,
           code = 3, length = 0.02, angle = 90, col = col[1])
  }
  
  
  par(new = T)
  
  ## Plotting LEFT axis data
  plot(x, yleft, axes = FALSE, ylab = "" , xlab = xlab, ylim = yleftlim,
       xlim = xlim, bg = bg[1],
       pch = pch[2], type = type[2], lwd = lwd[2], col = col[2], ...)
  
  box()
  
  axis(2, pretty(yleftlim, length), col = 1, col.axis = 1)
  
  if (is.null(yleft2) == FALSE) {
    points(x2, yleft2, type = type[2], pch = pch[2], bg = bg[2],
           lwd = lwd[2], col = col2[2], ...)
  }
  
  
  #if (linky) lines(x, yleft, col = col[2], ...)
  
  if (smooth != 0) lines(supsmu(x, yleft, span = smooth), col = col[2], lwd=lwds, ...)
  
  if(yylab[2] == "") {
    mtext(deparse(substitute(yleft)), side = 2, outer = FALSE, line = 2, col = col[2], cex = 0.9, ...)
  } else {
    mtext(yylab[2], side = 2, outer = FALSE, line = 2, col = col[2], cex = 0.9, ...)
  }
  
  if (!missing(yleftErrorValues) && ! is.null(yleftErrorValues)) {
    arrows(x, yleft - yleftErrorValues, x, yleft + yleftErrorValues,
           code = 3, length = 0.02, angle = 90, col = col[2])
  }
  
  ## X-axis
  ##  axis(1, at = pretty(xlim, length))  ## Comment last line
  
  
}

#' plotmore
#'
#' \code{plotmore} Graphical representation of the relationship between genes and regulators.
#' 
#' @param output Output object of MORE main function.
#' 
#' @param gene ID of the gene to be plotted.
#' @param regulator ID of the regulator to be plotted. If NULL (default), all regulators of the gene are plotted.
#' @param reguValues Vector containing the values of a regulator. If NULL (default), these values are taken from the output object as long as they are available. 
#' @param plotPerOmic If TRUE, all the relevant/significant regulators of the given gene and the same omic are plotted in the same graph. If FALSE (default), each regulator is plotted in a separate plot.
#' @param gene.col Color to plot the gene. By default, 1 (black). 
#' @param regul.col  Color to plot the regulator. If NULL (default), a color will be assigned by the function, that will be different for each regulatory omic.
#' @param order If TRUE (default), the values in X-axis are ordered.
#' @param xlab Label for the X-axis.
#' @param cont.var  Vector with length equal to the number of observations in data, which optionally may contain the values of the numerical variable (e.g. time) to be plotted on the X-axis. By default, NULL.
#' @param cond2plot Vector or factor indicating the experimental group of each value to represent. If NULL (default), the labels are taken from the experimental design matrix. 
#' 
#' @return Graphical representation of the relationship between genes and regulators.
#'
#' @export


plotmore = function(output, gene, regulator = NULL, simplify = FALSE, reguValues = NULL, plotPerOmic = FALSE,
                    gene.col = 1, regu.col = NULL, order = TRUE,
                    xlab = "", cont.var = NULL, cond2plot = NULL,...) {
  
  if(simplify){
    # from which omic is the regulator?
    SigniReguGene = GetPairsGeneRegulator(genes = gene, output = output)
    omic = SigniReguGene[SigniReguGene[,"regulator"] == regulator,'omic']
    
    if (output$arguments$omic.type[omic]==0){
      df <- data.frame(
        gen = unlist(output$arguments$GeneExpression[gene,,drop=TRUE]),
        regulador = unlist(output$arguments$dataOmics[[omic]][regulator,,drop=TRUE]),
        Group = output$arguments$groups)
      
    # Create a scatterplot
      
      num_unique <- length(unique(df$Group))+1
      color_palette <- colorbiostat(num_unique)
      custom_colors <- setNames(color_palette[-1], unique(df$Group))
      ggplot2::ggplot(df, aes(x = regulador, y = gen, color = Group)) +
        geom_point() + scale_color_manual(values = custom_colors)+
        theme_minimal()+
        #geom_abline(intercept = c(coefs[1,2],coefs[2,2],coefs[3,2]),slope = c(coefs[1,3],coefs[2,3],coefs[3,3]),color=c('#15918A','#74CDF0','#EE446F'),linetype=c('solid','solid',"dashed"))+
        #geom_abline(intercept = 0,slope = coefs[2,2],color='#74CDF0')+
        #geom_abline(intercept = 0,slope = coefs[3,2],color='#EE446F',linetype="dashed")+
        labs( x = paste("Regulator\n",regulator), y = paste("Gene Expression\n",gene))
      #geom_smooth(method = "lm", se = FALSE, aes(group = group)) 
      #+geom_abline(intercept = intercept, slope = slope, color="red",  
      #linetype="dashed", size=1.5)+ 
      
    } else {
      
      df <- data.frame(
        gen = unlist(output$arguments$GeneExpression[gene,,drop=TRUE]),
        regulador = unlist(output$arguments$dataOmics[[omic]][regulator,,drop=TRUE]),
        Group = output$arguments$groups)
      df$regulador<-factor(df$regulador)
      
      num_unique <- length(unique(df$Group))+1
      color_palette <- colorbiostat(num_unique)
      custom_colors <- setNames(color_palette[-1], unique(df$Group))
      # Create a scatterplot
      ggplot2::ggplot(df, aes(x = regulador, y = gen,fill=Group)) + theme_minimal()+
        geom_boxplot() + scale_fill_manual(values = custom_colors)+  scale_color_manual(values = custom_colors)+
        scale_x_discrete(labels = c('0','1')) + stat_summary(aes(color = Group),fun='median',geom = 'point', position = position_dodge(width = 0.75))+
      labs( x = paste("Regulator \n",regulator), y = paste("Gene Expression\n",gene))
      
    }
  } else{
    if(output$arguments$method=='glm'||output$arguments$method=='pls2'){
      
      return(plotGLM(output, gene, regulator = regulator, reguValues = reguValues, plotPerOmic = plotPerOmic,
                     gene.col = gene.col, regu.col = regu.col, order = order,
                     xlab = xlab, cont.var = cont.var, cond2plot = cond2plot,...))
    }
    
    if(output$arguments$method=='pls1'||output$arguments$method=='pls2'){
      
      return(plotPLS(output, gene, regulator = regulator, reguValues = reguValues, plotPerOmic = plotPerOmic,
                     gene.col = gene.col, regu.col = regu.col, order = order,
                     xlab = xlab, cont.var = cont.var, cond2plot = cond2plot,...))
    }
  }
  
 
}

# order: Should the experimental groups be ordered for the plot? If TRUE, omic values are also ordered accordingly.
#        If FALSE, the function assumes they were provided in the right order for a meaningful plot.

plotGLM = function (GLMoutput, gene, regulator = NULL, reguValues = NULL, plotPerOmic = FALSE,
                    gene.col = 1, regu.col = NULL, order = TRUE,
                    xlab = "", cont.var = NULL, cond2plot = NULL, verbose =TRUE, ...) {
  
  # Colors for omics
  omic.col = colors()[c(554,89,111,512,17,586,132,428,601,568,86,390,
                        100,200,300,400,500,10,20,30,40,50,60,70,80,90,150,250,350,450,550)]
  
  if (is.null(regu.col)) {
    any.col = omic.col[1:length(GLMoutput$arguments$dataOmics)]
  } else {
    if (length(regu.col) == length(GLMoutput$arguments$dataOmics)) {
      any.col = regu.col
    } else {
      any.col = rep(regu.col, length(GLMoutput$arguments$dataOmics))
    }
  }
  names(any.col) = names(GLMoutput$arguments$dataOmics)
  
  
  # Changing margin
  par(mar = c(5,3,4,3)+0.1)
  
  # Groups to plot
  if (is.null(cond2plot)) {
    if (!is.null(GLMoutput$arguments$edesign)) {
      cond2plot = apply(GLMoutput$arguments$edesign, 1, paste, collapse = "_")
    }
  }
  
  # Replicates
  if (!is.null(cont.var)) {  # we have continuous variable
    if (!is.null(cond2plot)) { # cont.var + cond2plot
      myreplicates = apply(cbind(cond2plot, cont.var), 1, paste, collapse = "_")
      
    } else {  # only continuous variable is provided
      myreplicates = cont.var
    }
    
  } else {   # no cont.var
    
    if (!is.null(cond2plot)) { # only cond2plot
      myreplicates = colnames(GLMoutput$arguments$GeneExpression)
    } else {  # nothing
      myreplicates = colnames(GLMoutput$arguments$GeneExpression)
    }
  }
  
  # Cast myreplicates to character
  myreplicates = as.character(myreplicates)
  names(myreplicates) = colnames(GLMoutput$arguments$GeneExpression)
  if (order) {
    myorder = order(myreplicates)
    myreplicates = sort(myreplicates)
    cond2plot = cond2plot[myorder]
  }
  myrepliUni = unique(myreplicates)
  
  
  if (max(table(myreplicates)) == 1) {
    replicates = FALSE
  } else {
    replicates = TRUE
  }
  
  
  if (is.null(cond2plot)) {
    numLines = NULL
  } else {
    condi1 = unique(cond2plot)
    num1 = 1:length(condi1); names(num1) = condi1
    num2 = num1[as.character(cond2plot)]
    if (replicates) {
      num2 = aggregate(num2, by = list("rep" = myreplicates), unique)$x
    }
    numLines = which(diff(num2) != 0)+0.5
  }
  
  
  # Error values
  getErrorValues = function(realValues, repsInfo) {
    # Disable it
    if (! replicates)
      return(NULL)
    
    out_values = tapply(realValues, repsInfo, function(reps) sd(reps)/sqrt(length(reps)))
    
    return(out_values)
  }
  
  
  myrepliUni = unique(myreplicates)
  
  ## REGULATOR = NULL
  
  if (is.null(regulator)) {  ### Plot all regulators for the given gene
    
    GLMgene = GLMoutput$ResultsPerGene[[gene]]
    
    if (is.null(GLMgene)) {
      stop(paste("No GLM was obtained for gene", gene))
    }
    
    if (is.null(GLMgene$relevantRegulators)) { ## No significant regulators
      
      cat("No relevant regulators were found for this gene.\n")
      
      
    } else
    {  ## Significant regulators:
      
      # Considering multicollinearity
      SigReg = GLMgene$allRegulators
      SigReg = SigReg[SigReg$Rel == 1, c("regulator", "omic", "area", "filter")]
      
      SigReg = SigReg[GLMgene$relevantRegulators,,drop = FALSE]
      
      cat(paste(nrow(SigReg), "relevant regulators are to be plotted for gene", gene)); cat("\n")
      
      # Gene values
      geneValues = GLMgene$Y$y
      if (order) geneValues = geneValues[myorder]
      errorValues = getErrorValues(geneValues, myreplicates)
      geneValues = tapply(geneValues, myreplicates, mean)
      geneValues = geneValues[myrepliUni]
      names(geneValues) = myrepliUni
      
      # X values
      x.points = 1:length(myrepliUni)
      eje = myrepliUni
      
      if (plotPerOmic) { ## All regulators from the same omic in the same plot
        
        myomics = unique(SigReg$omic)
        
        for (oo in myomics) {
          
          SigRegOmic = SigReg[SigReg$omic == oo,]
          
          omicValues = t(GLMoutput$arguments$dataOmics[[oo]])
          omicValues = omicValues[rownames(GLMgene$X),]
          reguValues = omicValues[, colnames(omicValues) == SigRegOmic$regulator[1]]
          if (order) reguValues = reguValues[myorder]
          errorValuesRegu = getErrorValues(reguValues, myreplicates)
          
          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni
          
          mycol = any.col[oo]
          
          if (nrow(SigRegOmic) == 1) {
            leftlab = SigRegOmic$regulator[1]
          } else { leftlab = oo }
          
          yleftlim = range(omicValues[,SigRegOmic$regulator], na.rm = TRUE)
          
          if (! is.null(errorValuesRegu)) {
            yleftlim = range(
              apply(omicValues[, SigRegOmic$regulator, drop = FALSE], 2, function(x) {
                if (order) x = x[myorder]
                errorInd = getErrorValues(x, myreplicates)
                meanValues = tapply(x, myreplicates, mean)
                meanValues = meanValues[myrepliUni]
                
                return(c(meanValues - errorInd, meanValues + errorInd))
              }))
          }
          
          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, mycol), yleftlim = yleftlim,
                       xlab = xlab, yylab = c(gene, leftlab), pch = c(16,16),
                       main = oo, numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)
          
          if (nrow(SigRegOmic) > 1) {
            for (i in 2:nrow(SigRegOmic)) {
              
              reguValues = omicValues[, colnames(omicValues) == SigRegOmic$regulator[i]]
              if (order) reguValues = reguValues[myorder]
              errorValuesRegu = getErrorValues(reguValues, myreplicates)
              
              reguValues = tapply(reguValues, myreplicates, mean)
              reguValues = reguValues[myrepliUni]
              names(reguValues) = myrepliUni
              
              lines(x.points, reguValues, type = "o", lwd = 2, pch = i, col = mycol, lty = i)
              
              if (! is.null(errorValuesRegu)) {
                arrows(x.points, reguValues - errorValuesRegu, x.points, reguValues + errorValuesRegu,
                       code = 3, length = 0.02, angle = 90, col = mycol)
              }
            }
          }
          
        }
        
      } else {  ## Each regulator in a separate plot
        
        for (rr in SigReg$regulator) {
          
          oo = GLMoutput$ResultsPerGene[[gene]]$allRegulators[rr,"omic"]
          omicValues = t(GLMoutput$arguments$dataOmics[[oo]])
          reguValues = omicValues[, colnames(omicValues) == rr]
          if (order) reguValues = reguValues[myorder]
          errorValuesRegu = getErrorValues(reguValues, myreplicates)
          
          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni
          
          mycol = any.col[oo]
          
          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, mycol), xlab = xlab,
                       yylab = c(gene, rr), pch = c(16,16),
                       main = paste(as.character(SigReg[rr, c("omic", "area")]), collapse = " "),
                       numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)
          
        }
        
      }
      
      return(GLMgene$allRegulators[GLMgene$relevantRegulators, -6])
    }
    
  }
  
  
  ## GENE = NULL
  
  if (is.null(gene)) {  ### Plot all genes regulated by the regulator
    
    SigniReguGene = GetPairsGeneRegulator(genes = NULL, output = GLMoutput)
    SigniReguGene = SigniReguGene[SigniReguGene[,"regulator"] == regulator,]
    myomics = SigniReguGene[,"omic"]
    myomic = unique(myomics)
    
    if (nrow(SigniReguGene) > 0) {  # When there are genes regulated by this regulator
      
      if (is.null(reguValues)) {  # User does not provide reguValues
        reguValues = as.numeric(GLMoutput$arguments$dataOmics[[myomic]][regulator,])
      }
      
      numGenes = length(SigniReguGene$gene)
      if(verbose) {cat(paste(numGenes, "genes are regulated by", regulator)); cat("\n")}
      
      if (length(reguValues) > 0) {  # reguValues are available (recovered or given by user)
        
        if (order) reguValues = reguValues[myorder]
        errorValuesRegu = getErrorValues(reguValues, myreplicates)
        reguValues = tapply(reguValues, myreplicates, mean)
        reguValues = reguValues[myrepliUni]
        names(reguValues) = myrepliUni
        
        lapply(1:numGenes, function (i) {
          
          geneValues = GLMoutput$ResultsPerGene[[SigniReguGene[i,"gene"]]]$Y$y
          if (order) geneValues = geneValues[myorder]
          errorValues = getErrorValues(geneValues, myreplicates)
          geneValues = tapply(geneValues, myreplicates, mean)
          geneValues = geneValues[myrepliUni]
          names(geneValues) = myrepliUni
          
          x.points = 1:length(myrepliUni)
          eje = myrepliUni
          
          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, any.col[myomics[i]]), xlab = xlab,
                       yylab = c(SigniReguGene[i,"gene"], regulator), pch = c(16,16),
                       main = paste(as.character(SigniReguGene[1,c("omic", "area")]), collapse = " "),
                       numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)
        })
      } else { cat("Regulator values could not be recovered from GLMoutput. Please provide them in reguValues argument to generate the plot.\n") }
      
      return(SigniReguGene$gene)
      
    } else { cat(paste("There are no genes relevantly regulated by", regulator)); cat("\n") }
    
  }
  
  
  ## GENE + REGULATOR
  
  if (!is.null(gene) && !is.null(regulator)) {  ### Plot only the given gene and the given regulator
    
    geneResults = GLMoutput$ResultsPerGene[[gene]]
    
    if (is.null(geneResults)) {
      stop(paste("No GLM was obtained for gene", gene))
    } else
    {
      myomic = geneResults$allRegulators[regulator, "omic"]
      
      if (is.null(reguValues)) {  # User does not provide reguValues
        reguValues = as.numeric(GLMoutput$arguments$dataOmics[[myomic]][regulator,]) # regulator values
      }
      
      if (length(reguValues) > 0) {  # reguValues are available (recovered or given by user)
        
        if (order)  reguValues = reguValues[myorder]
        errorValuesRegu = getErrorValues(reguValues, myreplicates)
        reguValues = tapply(reguValues, myreplicates, mean, na.rm = TRUE)
        reguValues = reguValues[myrepliUni]
        names(reguValues) = myrepliUni
        
        geneValues = GLMoutput$ResultsPerGene[[gene]]$Y$y
        if (order) geneValues = geneValues[myorder]
        errorValues = getErrorValues(geneValues, myreplicates)
        geneValues = tapply(geneValues, myreplicates, mean)
        geneValues = geneValues[myrepliUni]
        names(geneValues) = myrepliUni
        
        x.points = 1:length(myrepliUni)
        eje = myrepliUni
        
        plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                     col = c(gene.col, any.col[myomic]), xlab = xlab,
                     yylab = c(gene, regulator), pch = c(16,16),
                     main = paste(as.character(geneResults$allRegulators[regulator, c("omic", "area")]),
                                  collapse = " "),
                     numLines = numLines, x.names = eje,
                     geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)
        
      } else {
        
        cat("Regulator values could not be recovered from GLMoutput.\n")
        
        regulator = geneResults$allRegulators[regulator,"filter"]
        
        if (regulator %in% rownames(geneResults$allRegulators)) {
          
          cat(paste(regulator, "values will be plotted instead.")); cat("\n")
          cat(paste(regulator, "summarizes information from the following correlated regulators:")); cat("\n")
          cat(geneResults$allRegulators[geneResults$allRegulators[,"filter"] == regulator,"regulator"]); cat("\n")
          
          reguValues = geneResults$X
          reguValues = reguValues[, colnames(reguValues) == regulator]
          if (order) reguValues = reguValues[myorder]
          errorValuesRegu = getErrorValues(reguValues, myreplicates)
          
          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni
          
          geneValues = GLMoutput$ResultsPerGene[[gene]]$Y$y
          if (order) geneValues = geneValues[myorder]
          errorValues = getErrorValues(geneValues, myreplicates)
          geneValues = tapply(geneValues, myreplicates, mean)
          geneValues = geneValues[myrepliUni]
          names(geneValues) = myrepliUni
          
          x.points = 1:length(myrepliUni)
          eje = myrepliUni
          
          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, any.col[myomic]), xlab = xlab,
                       yylab = c(gene, regulator), pch = c(16,16),
                       main = paste(as.character(geneResults$allRegulators[regulator, c("omic", "area")]),
                                    collapse = " "),
                       numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)
        } else {
          cat("The selected regulator was not declared as relevant by the ElasticNet\n")
          cat("Please, either select another regulator or provide the regulator values.\n")
        }
        
      }
      
    }
    
  }
  
}

plotPLS = function (PLSoutput, gene, regulator = NULL, reguValues = NULL, plotPerOmic = FALSE,
                    gene.col = 1, regu.col = NULL, order = TRUE,
                    xlab = "", cont.var = NULL, cond2plot = NULL, verbose = TRUE,...) {
  
  # Colors for omics
  omic.col = colors()[c(554,89,111,512,17,586,132,428,601,568,86,390,
                        100,200,300,400,500,10,20,30,40,50,60,70,80,90,150,250,350,450,550)]
  
  if (is.null(regu.col)) {
    any.col = omic.col[1:length(PLSoutput$arguments$dataOmics)]
  } else {
    if (length(regu.col) == length(PLSoutput$arguments$dataOmics)) {
      any.col = regu.col
    } else {
      any.col = rep(regu.col, length(PLSoutput$arguments$dataOmics))
    }
  }
  names(any.col) = names(PLSoutput$arguments$dataOmics)
  
  
  # Changing margin
  par(mar = c(5,3,4,3)+0.1)
  
  # Groups to plot
  if (is.null(cond2plot)) {
    if (!is.null(PLSoutput$arguments$edesign)) {
      cond2plot = apply(PLSoutput$arguments$edesign, 1, paste, collapse = "_")
    }
  }
  
  # Replicates
  if (!is.null(cont.var)) {  # we have continuous variable
    if (!is.null(cond2plot)) { # cont.var + cond2plot
      myreplicates = apply(cbind(cond2plot, cont.var), 1, paste, collapse = "_")
      
    } else {  # only continuous variable is provided
      myreplicates = cont.var
    }
    
  } else {   # no cont.var
    
    if (!is.null(cond2plot)) { # only cond2plot
      myreplicates = colnames(PLSoutput$arguments$GeneExpression)
    } else {  # nothing
      myreplicates = colnames(PLSoutput$arguments$GeneExpression)
    }
  }
  
  # Cast myreplicates to character
  myreplicates = as.character(myreplicates)
  names(myreplicates) = colnames(PLSoutput$arguments$GeneExpression)
  if (order) {
    myorder = order(myreplicates)
    myreplicates = sort(myreplicates)
    cond2plot = cond2plot[myorder]
  }
  myrepliUni = unique(myreplicates)
  
  
  if (max(table(myreplicates)) == 1) {
    replicates = FALSE
  } else {
    replicates = TRUE
  }
  
  
  if (is.null(cond2plot)) {
    numLines = NULL
  } else {
    condi1 = unique(cond2plot)
    num1 = 1:length(condi1); names(num1) = condi1
    num2 = num1[as.character(cond2plot)]
    if (replicates) {
      num2 = aggregate(num2, by = list("rep" = myreplicates), unique)$x
    }
    numLines = which(diff(num2) != 0)+0.5
  }
  
  
  # Error values
  getErrorValues = function(realValues, repsInfo) {
    # Disable it
    if (! replicates)
      return(NULL)
    
    out_values = tapply(realValues, repsInfo, function(reps) sd(reps)/sqrt(length(reps)))
    
    return(out_values)
  }
  
  
  myrepliUni = unique(myreplicates)
  
  ## REGULATOR = NULL
  
  if (is.null(regulator)) {  ### Plot all regulators for the given gene
    
    PLSgene = PLSoutput$ResultsPerGene[[gene]]
    
    if (is.null(PLSgene)) {
      stop(paste("No GLM was obtained for gene", gene))
    }
    
    if (is.null(PLSgene$significantRegulators)) { ## No significant regulators
      
      cat("No significant regulators were found for this gene.\n")
      
      
    } else
    {  ## Significant regulators:
      
      # Considering multicollinearity
      SigReg = PLSgene$allRegulators
      SigReg = SigReg[SigReg$Sig == 1, c("regulator", "omic", "area", "filter")]
      
      SigReg = SigReg[PLSgene$significantRegulators,,drop = FALSE]
      
      cat(paste(nrow(SigReg), "significant regulators are to be plotted for gene", gene)); cat("\n")
      
      # Gene values
      geneValues = PLSgene$Y$y
      if (order) geneValues = geneValues[myorder]
      errorValues = getErrorValues(geneValues, myreplicates)
      geneValues = tapply(geneValues, myreplicates, mean)
      geneValues = geneValues[myrepliUni]
      names(geneValues) = myrepliUni
      
      # X values
      x.points = 1:length(myrepliUni)
      eje = myrepliUni
      
      if (plotPerOmic) { ## All regulators from the same omic in the same plot
        
        myomics = unique(SigReg$omic)
        
        for (oo in myomics) {
          
          SigRegOmic = SigReg[SigReg$omic == oo,]
          
          omicValues = t(PLSoutput$arguments$dataOmics[[oo]])
          omicValues = omicValues[rownames(PLSgene$X),]
          reguValues = omicValues[, colnames(omicValues) == SigRegOmic$regulator[1]]
          if (order) reguValues = reguValues[myorder]
          errorValuesRegu = getErrorValues(reguValues, myreplicates)
          
          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni
          
          mycol = any.col[oo]
          
          if (nrow(SigRegOmic) == 1) {
            leftlab = SigRegOmic$regulator[1]
          } else { leftlab = oo }
          
          yleftlim = range(omicValues[,SigRegOmic$regulator], na.rm = TRUE)
          
          if (! is.null(errorValuesRegu)) {
            yleftlim = range(
              apply(omicValues[, SigRegOmic$regulator, drop = FALSE], 2, function(x) {
                if (order) x = x[myorder]
                errorInd = getErrorValues(x, myreplicates)
                meanValues = tapply(x, myreplicates, mean)
                meanValues = meanValues[myrepliUni]
                
                return(c(meanValues - errorInd, meanValues + errorInd))
              }))
          }
          
          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, mycol), yleftlim = yleftlim,
                       xlab = xlab, yylab = c(gene, leftlab), pch = c(16,16),
                       main = oo, numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)
          
          if (nrow(SigRegOmic) > 1) {
            for (i in 2:nrow(SigRegOmic)) {
              
              reguValues = omicValues[, colnames(omicValues) == SigRegOmic$regulator[i]]
              if (order) reguValues = reguValues[myorder]
              errorValuesRegu = getErrorValues(reguValues, myreplicates)
              
              reguValues = tapply(reguValues, myreplicates, mean)
              reguValues = reguValues[myrepliUni]
              names(reguValues) = myrepliUni
              
              lines(x.points, reguValues, type = "o", lwd = 2, pch = i, col = mycol, lty = i)
              
              if (! is.null(errorValuesRegu)) {
                arrows(x.points, reguValues - errorValuesRegu, x.points, reguValues + errorValuesRegu,
                       code = 3, length = 0.02, angle = 90, col = mycol)
              }
            }
          }
          
        }
        
      } else {  ## Each regulator in a separate plot
        
        for (rr in SigReg$regulator) {
          
          oo = PLSoutput$ResultsPerGene[[gene]]$allRegulators[rr,"omic"]
          omicValues = t(PLSoutput$arguments$dataOmics[[oo]])
          reguValues = omicValues[, colnames(omicValues) == rr]
          if (order) reguValues = reguValues[myorder]
          errorValuesRegu = getErrorValues(reguValues, myreplicates)
          
          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni
          
          mycol = any.col[oo]
          
          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, mycol), xlab = xlab,
                       yylab = c(gene, rr), pch = c(16,16),
                       main = paste(as.character(SigReg[rr, c("omic", "area")]), collapse = " "),
                       numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)
          
        }
        
      }
      
      return(PLSgene$allRegulators[PLSgene$significantRegulators, -6])
    }
    
  }
  
  
  ## GENE = NULL
  
  if (is.null(gene)) {  ### Plot all genes regulated by the regulator
    
    SigniReguGene = GetPairsGeneRegulator(genes = NULL, output = PLSoutput)
    SigniReguGene = SigniReguGene[SigniReguGene[,"regulator"] == regulator,]
    myomics = SigniReguGene[,"omic"]
    myomic = unique(myomics)
    
    if (nrow(SigniReguGene) > 0) {  # When there are genes regulated by this regulator
      
      if (is.null(reguValues)) {  # User does not provide reguValues
        reguValues = as.numeric(PLSoutput$arguments$dataOmics[[myomic]][regulator,])
      }
      
      numGenes = length(SigniReguGene$gene)
      if(verbose) {cat(paste(numGenes, "genes are regulated by", regulator)); cat("\n")}
      
      if (length(reguValues) > 0) {  # reguValues are available (recovered or given by user)
        
        if (order) reguValues = reguValues[myorder]
        errorValuesRegu = getErrorValues(reguValues, myreplicates)
        reguValues = tapply(reguValues, myreplicates, mean)
        reguValues = reguValues[myrepliUni]
        names(reguValues) = myrepliUni
        
        lapply(1:numGenes, function (i) {
          
          geneValues = PLSoutput$ResultsPerGene[[SigniReguGene[i,"gene"]]]$Y$y
          if (order) geneValues = geneValues[myorder]
          errorValues = getErrorValues(geneValues, myreplicates)
          geneValues = tapply(geneValues, myreplicates, mean)
          geneValues = geneValues[myrepliUni]
          names(geneValues) = myrepliUni
          
          x.points = 1:length(myrepliUni)
          eje = myrepliUni
          
          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, any.col[myomics[i]]), xlab = xlab,
                       yylab = c(SigniReguGene[i,"gene"], regulator), pch = c(16,16),
                       main = paste(as.character(SigniReguGene[1,c("omic", "area")]), collapse = " "),
                       numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)
        })
      } else { cat("Regulator values could not be recovered from output. Please provide them in reguValues argument to generate the plot.\n") }
      
      return(SigniReguGene$gene)
      
    } else { cat(paste("There are no genes significantly regulated by", regulator)); cat("\n") }
    
  }
  
  
  ## GENE + REGULATOR
  
  if (!is.null(gene) && !is.null(regulator)) {  ### Plot only the given gene and the given regulator
    
    geneResults = PLSoutput$ResultsPerGene[[gene]]
    
    if (is.null(geneResults)) {
      stop(paste("No GLM was obtained for gene", gene))
    } else
    {
      myomic = geneResults$allRegulators[regulator, "omic"]
      
      if (is.null(reguValues)) {  # User does not provide reguValues
        reguValues = as.numeric(PLSoutput$arguments$dataOmics[[myomic]][regulator,]) # regulator values
      }
      
      if (length(reguValues) > 0) {  # reguValues are available (recovered or given by user)
        
        if (order)  reguValues = reguValues[myorder]
        errorValuesRegu = getErrorValues(reguValues, myreplicates)
        reguValues = tapply(reguValues, myreplicates, mean, na.rm = TRUE)
        reguValues = reguValues[myrepliUni]
        names(reguValues) = myrepliUni
        
        geneValues = PLSoutput$ResultsPerGene[[gene]]$Y$y
        if (order) geneValues = geneValues[myorder]
        errorValues = getErrorValues(geneValues, myreplicates)
        geneValues = tapply(geneValues, myreplicates, mean)
        geneValues = geneValues[myrepliUni]
        names(geneValues) = myrepliUni
        
        x.points = 1:length(myrepliUni)
        eje = myrepliUni
        
        plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                     col = c(gene.col, any.col[myomic]), xlab = xlab,
                     yylab = c(gene, regulator), pch = c(16,16),
                     main = paste(as.character(geneResults$allRegulators[regulator, c("omic", "area")]),
                                  collapse = " "),
                     numLines = numLines, x.names = eje,
                     geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)
        
      } else {
        
        cat("Regulator values could not be recovered from output.\n")
        
        regulator = geneResults$allRegulators[regulator,"filter"]
        
        if (regulator %in% rownames(geneResults$allRegulators)) {
          
          cat(paste(regulator, "values will be plotted instead.")); cat("\n")
          cat(paste(regulator, "summarizes information from the following correlated regulators:")); cat("\n")
          cat(geneResults$allRegulators[geneResults$allRegulators[,"filter"] == regulator,"regulator"]); cat("\n")
          
          reguValues = geneResults$X
          reguValues = reguValues[, colnames(reguValues) == regulator]
          if (order) reguValues = reguValues[myorder]
          errorValuesRegu = getErrorValues(reguValues, myreplicates)
          
          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni
          
          geneValues = PLSoutput$ResultsPerGene[[gene]]$Y$y
          if (order) geneValues = geneValues[myorder]
          errorValues = getErrorValues(geneValues, myreplicates)
          geneValues = tapply(geneValues, myreplicates, mean)
          geneValues = geneValues[myrepliUni]
          names(geneValues) = myrepliUni
          
          x.points = 1:length(myrepliUni)
          eje = myrepliUni
          
          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, any.col[myomic]), xlab = xlab,
                       yylab = c(gene, regulator), pch = c(16,16),
                       main = paste(as.character(geneResults$allRegulators[regulator, c("omic", "area")]),
                                    collapse = " "),
                       numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)
        } else {
          cat("The selected regulator was not declared as significant by the ElasticNet\n")
          cat("Please, either select another regulator or provide the regulator values.\n")
        }
        
      }
      
    }
    
  }
  
}

## Plots PLS ----------

#' plotweight
#'
#' \code{plotweight} Function to be applied to \link{more} main function output.
#' 
#' @param output Output object of MORE main function.
#' @param gene Gene of which the user wants visualize the weightings of the model.
#' @param axe1 Number of the PLS component to plot in the X axis
#' @param axe2 Number of the PLS component to plot in the Y axis
#' 
#' @return Plots the weighting star of the regulators identified as significant for the selected gene.
#' @export

plotweight<-function(output, gene,axe1=1,axe2=2){
  
  if(output$GlobalSummary$GoodnessOfFit[gene,'ncomp']==1) {
    warning('The original component extracted a unique component. The visuallization would be hard')
    if (ncol(output$arguments$GeneExpression)<7){cross = output$arguments$GeneExpression -2}else{cross =7}
    pls = ropls::opls(output$ResultsPerGene[[gene]]$X[,output$ResultsPerGene[[gene]]$significantRegulators,drop=FALSE], output$ResultsPerGene[[gene]]$Y$y, info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0, predI=output$GlobalSummary$GoodnessOfFit[gene,'ncomp'])
    
    plot(pls@weightStarMN[,1], rep(0,length(output$ResultsPerGene[[gene]]$significantRegulators)),
         main = "Weights*",
         xlab = paste('w*c', axe1), ylab = paste('w*c', axe2),
         pch = 18, col = "blue", xlim = c(-1,1))
    
    points(pls@cMN[,1], 0, pch = 18, col = "red")
    # Asignamos las etiquetas
    text(pls@weightStarMN[,1], rep(0,length(output$ResultsPerGene[[gene]]$significantRegulators)),
         labels = row.names(pls@weightStarMN),
         cex = 0.6, srt = 60, pos = 2, col = "black")
    text(pls@cMN[,1], 0, labels = gene, cex =
           0.6, srt = 60, pos = 4, col = 'black')
    abline(h=0)
    legend("topleft", legend = c("X", gene),cex = 0.5,
           pch = 18, col = c("blue", "red"))
  } else{
    if(axe1> output$GlobalSummary$GoodnessOfFit[gene,'ncomp'] | axe2>output$GlobalSummary$GoodnessOfFit[gene,'ncomp']) stop('Error: The model did not extracted originally that many components')
    if (ncol(output$arguments$GeneExpression)<7){cross = output$arguments$GeneExpression -2}else{cross =7}
    #Create the PLS model only with the variables that resulted significant in the model
    pls = ropls::opls(output$ResultsPerGene[[gene]]$X[,output$ResultsPerGene[[gene]]$significantRegulators,drop=FALSE], output$ResultsPerGene[[gene]]$Y$y, info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0, predI=output$GlobalSummary$GoodnessOfFit[gene,'ncomp'])

    #Create the weighting plots
    plot(pls@weightStarMN[,axe1], pls@weightStarMN[,axe2],
         main = "Weights*",
         xlab = paste('w*c', axe1), ylab = paste('w*c', axe2),
         pch = 18, col = "blue", asp=1)
    
    points(pls@cMN[,axe1], pls@cMN[,axe2], pch = 18, col = "red")
    # Asignamos las etiquetas
    text(pls@weightStarMN[,axe1], pls@weightStarMN[,axe2],
         labels = row.names(pls@weightStarMN),
         cex = 0.6, pos = 4, col = "black")
    text(pls@cMN[,axe1], pls@cMN[,axe2], labels = gene, cex =
           0.6, pos = 4, col = 'black')
    abline(h=0, v=0)
    legend("topleft", legend = c("X", gene),cex = 0.5,
           pch = 18, col = c("blue", "red"))
  }
  
}

#' plotscores
#'
#' \code{plotscores} Function to be applied to \link{more} main function output.
#' 
#' @param output Output object of MORE main function.
#' @param gene Gene of which the user wants visualize the weightings of the model.
#' @param axe1 Number of the PLS component to plot in the X axis
#' @param axe2 Number of the PLS component to plot in the Y axis
#' 
#' @return Plots the scores of the samples under the PLS model generated by the significant regulators for the selected gene.
#' @export

plotscores<-function(output, gene,axe1=1,axe2=2){
  
  if(output$GlobalSummary$GoodnessOfFit[gene,'ncomp']==1) {
    warning('The original component extracted a unique component. The visuallization would be hard')
    if (ncol(output$arguments$GeneExpression)<7){cross = output$arguments$GeneExpression -2}else{cross =7}
    pls = ropls::opls(output$ResultsPerGene[[gene]]$X[,output$ResultsPerGene[[gene]]$significantRegulators,drop=FALSE], output$ResultsPerGene[[gene]]$Y$y, info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0, predI=output$GlobalSummary$GoodnessOfFit[gene,'ncomp'])
    
    plot(pls@scoreMN[,1], rep(0,length(output$arguments$groups)),
         main = "Scores",
         xlab = paste('t', axe1), ylab = paste('t', axe2),
         pch = 18, col = output$arguments$groups)
    
    # Asignamos las etiquetas
    text(pls@scoreMN[,1], rep(0,length(output$arguments$groups)),
         labels = row.names(pls@scoreMN),
         cex = 0.6, srt = 60, pos = 2, col = 'black')
    abline(h=0)
  } else{
    if(axe1> output$GlobalSummary$GoodnessOfFit[gene,'ncomp'] | axe2>output$GlobalSummary$GoodnessOfFit[gene,'ncomp']) stop('Error: The model did not extracted originally that many components')
    if (ncol(output$arguments$GeneExpression)<7){cross = output$arguments$GeneExpression -2}else{cross =7}
    #Create the PLS model only with the variables that resulted significant in the model
    pls = ropls::opls(output$ResultsPerGene[[gene]]$X[,output$ResultsPerGene[[gene]]$significantRegulators,drop=FALSE], output$ResultsPerGene[[gene]]$Y$y, info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0, predI=output$GlobalSummary$GoodnessOfFit[gene,'ncomp'])
    #Create the weighting plots
    plot(pls@scoreMN[,axe1], pls@scoreMN[,axe2],
         main = "Scores",
         xlab = paste('t', axe1), ylab = paste('t', axe2),
         pch = 18, col = output$arguments$groups)
    
    # Asignamos las etiquetas
    text(pls@scoreMN[,axe1], pls@scoreMN[,axe2],
         labels = row.names(pls@scoreMN),
         cex = 0.6, pos = 4, col = "black")
    abline(h=0, v=0)
  }
  
}


## Summary ------------


#' summary
#'
#' \code{summary.MORE} Function to be applied to MORE object.
#' 
#' @param object MORE object obtained from applying \link{more} function. 
#' @param plot.more If TRUE top 10 global regulators are plotted against the genes they regulate. By default, FALSE.
#' 
#' @return Summary of more analysis.
#' @export

summary.MORE <-function(object, plot.more=FALSE){
  
  cat('A model was computed for',length(object$ResultsPerGene), 'genes.' ,'\n')
  cat(ifelse(is.null(object$GlobalSummary$GenesNoregulators),0,nrow(object$GlobalSummary$GenesNoregulators)), 'genes had no intial regulators.' ,'\n')
  
  if(object$arguments$method == 'glm'||object$arguments$method=='isgl'){
    cat('For', ifelse(is.null(object$GlobalSummary$GenesNOmodel),0,length(object$GlobalSummary$GenesNOmodel)), 'genes, the final GLM model could not be obtained.','\n')
    cat('Genes presented a mean of ',mean(na.omit(object$GlobalSummary$GoodnessOfFit[,'relReg'])),'relevant regulators.','\n')
    
    #Top hub genes
    relevant_regulators<-object$GlobalSummary$ReguPerGene[,c(grep('-Rel$',colnames(object$GlobalSummary$ReguPerGene))),drop=FALSE]
    #globally
    
    s_rel_reg<-apply(relevant_regulators, 1, sum)
    
    cat('These are the top 10 hub genes and the number of relevant regulators for each:\n')
    print(s_rel_reg[rev(tail(order(s_rel_reg),10))])
    
    #Global regulators
    
    m_rel_reg<-lapply(object$ResultsPerGene, function(x) x$relevantRegulators)
    m_rel_reg <- unlist(m_rel_reg)
    
    ## Count occurrences
    mrel_vector <- table(m_rel_reg)
    #Ask to regulate at least 10 genes
    mrel_vector<-mrel_vector[mrel_vector>10]
    if(length(mrel_vector!=0)){
      cat('These are the top 10 global regulators and the number of genes that they regulate:\n')
      print(mrel_vector[rev(tail(order(mrel_vector),10))])
    } else{
      cat('There were not global regulators (regulators that regulate more than 10 genes).')
    }
    
    
    if(plot.more){
      mreg<-mrel_vector[rev(tail(order(mrel_vector),10))]
      for (i in 1:10) {
        par(mfrow=c(2,4))
        plotGLM(object, gene = NULL, regulator = names(mreg)[i], plotPerOmic = FALSE ,order = FALSE, gene.col = 'skyblue', regu.col = 'tan1', verbose = FALSE)
      }
    }
  }
  else{
    cat('Genes presented a mean of ',mean(na.omit(object$GlobalSummary$GoodnessOfFit[,'sigReg'])),'significant regulators.','\n')
    
    #Top hub genes
    significant_regulators<-object$GlobalSummary$ReguPerGene[,c(grep('-Sig$',colnames(object$GlobalSummary$ReguPerGene))),drop=FALSE]
    #globally
    
    s_sig_reg<-apply(significant_regulators, 1, sum)
    cat('These are the top 10 hub genes and the number of significant regulators for each:\n')
    print(s_sig_reg[tail(order(s_sig_reg),10)])
    
    #Global regulators
    
    m_sig_reg<-lapply(object$ResultsPerGene, function(x) x$significantRegulators)
    m_sig_reg <- unlist(m_sig_reg)
    
    ## Count occurrences
    msig_vector <- table(m_sig_reg)
    #Ask to regulate at least 10 genes
    msig_vector<-msig_vector[msig_vector>10]
    cat('These are the top 10 global regulators and the number of genes that they regulate:\n')
    print(msig_vector[rev(tail(order(msig_vector),10))])
    
    if(plot.more){
      msig<-msig_vector[rev(tail(order(msig_vector),10))]
      for (i in 1:10) {
        par(mfrow=c(2,4))
        plotPLS(object, gene = NULL, regulator = names(msig)[i], plotPerOmic = FALSE ,order = FALSE, gene.col = 'skyblue', regu.col = 'tan1', verbose = FALSE)
        
      }
    }
    
  }
  
}

#' summary_plot
#'
#' \code{summary_plot} Function to be applied to MORE object and the object obtained from RegulationPerCondition.
#' 
#' @param output MORE object obtained from applying \link{more} function. 
#' @param output_regpcond Object generated by the function \link{RegulationPerCondition} when applied to a \link{more} object. 
#' @param by_genes If TRUE, the function plots the percentage of genes with significant regulators globally and per omic. If FALSE, it plots the percentage of significant regulations per omic. By default, TRUE.
#' 
#' @return Summary plot of the MORE analysis.
#' @export 


summary_plot<-function(output, output_regpcond, by_genes =TRUE){
  
  #output should by a MORE object
  #output_regpcond should by the output of calculating RegulationPerCondition to a MORE object
  #by_genes by default TRUE calculates the percentage of genes with significant regulators globally and per omic. FALSE to calculate the percentage of significant regulations globally and per omic.
  
  if(by_genes){
    #Calculate the vector with % of significant regulators by condition and globally
    
    ngroups = length(unique(output$arguments$groups))
    omics = names(output$arguments$dataOmics)
    totalgenes = length(output$ResultsPerGene)+ifelse(is.null(output$GlobalSummary$GenesNoregulators),0,length(output$GlobalSummary$GenesNoregulators))+ifelse(is.null(output$GlobalSummary$GenesNomodel),0,length(output$GlobalSummary$GenesNomodel))
    
    pos = grep('Group',colnames(output_regpcond))[1]
    #Create all the counts needed globally and per groups
    
    cts = matrix(NA, nrow=(ngroups)+1,ncol=length(omics)+1)
    
    for (i in 1:((ngroups)+1)){
      #Create the global values
      for (j in 1:((length(omics))+1)){
        if (i==1 && j==1){
          cts[i,j]= length(unique(output_regpcond$gene))
        }else if(i==1){
          cts[i,j]= length(unique(output_regpcond[output_regpcond$omic==omics[j-1],]$gene))
        }else if(j==1){
          cts[i,j] = length(unique(output_regpcond[output_regpcond[,pos+(i-2)]!=0,]$gene))
        }else{
          cts[i,j] = length(intersect(unique(output_regpcond[output_regpcond[,pos+(i-2)]!=0,]$gene), unique(output_regpcond[output_regpcond$omic==omics[j-1],]$gene)))
        }
        
      }
    }
    
    #Convert to percentage
    
    cts<-cts/totalgenes*100
    
    group_levels <- c('Global', unique(output$arguments$groups))
    #Create a df with the percentage of genes with significant regulators by omic and condition
    df <- data.frame(Group=factor(rep(group_levels, times = length(omics) + 1), levels = group_levels),
                     omic=rep(c('Any',names(output$arguments$dataOmics)),each = ngroups+1),
                     genes=as.vector(cts))
    
    num_unique <- ngroups+1
    color_palette <- colorbiostat(num_unique)
    custom_colors <- setNames(color_palette, unique(df$Group))
    ggplot2::ggplot(data=df, aes(x=omic, y=genes, fill=Group)) +
      geom_bar(stat="identity", position=position_dodge()) +
      theme_minimal()+
      scale_fill_manual(values = custom_colors)+ 
      labs(x="Omic", y = "% genes with significant regulators") +
      theme(legend.text = element_text(size = 12),panel.grid = element_line(color = "black",linewidth = 0.5,linetype = 1)) 
    
  } else{
    #Calculate the vector with % of significant regulations by condition in each omic
    
    ngroups = length(unique(output$arguments$groups))
    omics = names(output$arguments$dataOmics)
    pos = grep('Group',colnames(output_regpcond))[1]
    #Create all the counts needed globally and per groups
    
    cts = matrix(NA, nrow=(ngroups),ncol=length(omics))
    
    total_reg_omic <- if (is.null(output$arguments$associations)) {
      sapply(output$arguments$dataOmics, nrow)
    } else {
      sapply(omics, function(x) length(intersect(rownames(output$arguments$dataOmics[[x]]) ,output$arguments$associations[[x]][,2])))
    }
    
    for (i in 1:ngroups){
      #Create the global values
      for (j in 1:length(omics)){
        
        cts[i,j] = length(intersect(output_regpcond[output_regpcond[,pos+i-1]!=0,]$regulator,output_regpcond[output_regpcond$omic==omics[j],]$regulator))/total_reg_omic[j]*100
        
      }
    }

    #Create a df with the percentage of genes with significant regulators by omic and condition
    df <- data.frame(Group=rep(unique(output$arguments$groups), times=length(omics)),
                     omic=rep(names(output$arguments$dataOmics),each = ngroups),
                     genes=as.vector(cts))
    
    num_unique <- ngroups+1
    color_palette <- colorbiostat(num_unique)
    custom_colors <- setNames(color_palette[-1], unique(df$Group))
    ggplot2::ggplot(data=df, aes(x=omic, y=genes, fill=Group)) +
      geom_bar(stat="identity", position=position_dodge()) +
      theme_minimal()+scale_x_discrete(labels = paste(unique(df$omic),'\n',total_reg_omic,'regulators')) +
      scale_fill_manual(values = custom_colors)+
      scale_y_continuous(limits = c(0, max(df$genes) + 1)) +
      labs(x="Omic", y = "% significant regulations") +
      theme(legend.text = element_text(size = 12),panel.grid = element_line(color = "black",linewidth = 0.5,linetype = 1)) 
    
  }
  
}

getallreg <- function(x, gene) {
  reg <- x[x[, 1] == gene, 2]
  return(as.character(reg))
}

#' globalreg_plot
#'
#' \code{globalreg_plot} Function to be applied to \link{RegulationInCondition} function output.
#' 
#' @param output_regincond Output object of running \link{RegulationInCondition} function.
#' @param by_network If TRUE, information would be plotted on a network instead of a corrplot. By default, FALSE.
#' 
#' @return Graphical visualization between global regulators and genes in a specific condition.
#' @export

globalreg_plot<-function(output_regincond, by_network=FALSE){
  
  #output_regincond: Output object of running RegulationInCondition function
  #by_network: By faulta, FALSE. If TRUE plots the results in a network

  regulators<-output_regincond$GlobalRegulators
  if (length(regulators)==0){
    stop("ERROR:No global regulator was identified for this condition")
  }
  if (by_network){
    
    #Take group column in RegulationPerCondition
    df<-output_regincond$RegulationInCondition[grepl(paste(regulators, collapse = "|"), output_regincond$RegulationInCondition$regulator),]
    
    #Create the graph
    mygraph = igraph::graph_from_data_frame(df, directed=F)
    
    odf<-df[,c(2,3)]
    odf<-rbind(odf, data.frame('regulator'=unique(df$gene),'omic'=rep('gene',length(unique(df$gene)))))
    odf<-unique(odf)
    rownames(odf)<-odf$regulator
    
    mygraph<-igraph::set.vertex.attribute(mygraph,'omic', index = igraph::V(mygraph), value = odf[igraph::V(mygraph)$name,]$omic)
    mygraph<-igraph::set.edge.attribute(mygraph, 'sign', index = igraph::E(mygraph), value = df[,4])
    igraph::E(mygraph)$sign<-df[,4]
    
    igraph::plot.igraph(mygraph,vertex.label.cex= 0.4, vertex.size = 4,vertex.color=as.factor(igraph::V(mygraph)$omic),edge.color = ifelse(igraph::E(mygraph)$sign >0, "blue", "red"), main='Gene - Global regulators network')
    legend("topright", legend = unique(igraph::V(mygraph)$omic), col = categorical_pal(length(unique(as.factor(igraph::V(mygraph)$omic)))), pch = 16, cex = 1.5, bty = "n")
    
    
  }else{
    
    # Identify the genes they regulate
    genes<-unique(output_regincond$RegulationInCondition[grepl(paste(regulators, collapse = "|"), output_regincond$RegulationInCondition$regulator),1])
    
    #Create the matrix 
    gen_reg<-matrix(0, nrow= length(genes),ncol=length(regulators))
    colnames(gen_reg)=regulators
    rownames(gen_reg)=genes
    
    for ( i in 1:nrow(gen_reg)){
      #Get all the potential regulators of the gene
      potential_regulator <- unlist(sapply(output_regincond$arguments$associations, function(x) getallreg(x,genes[i])),use.names = FALSE)
      
      #Use NA for any potential regulator
      gen_reg[i,intersect(potential_regulator,regulators)]<-NA
      for ( j in 1:ncol(gen_reg)){
        if ( length(intersect(which(output_regincond$RegulationInCondition$gene==genes[i]),which(output_regincond$RegulationInCondition$regulator==regulators[j])) )!=0){
          gen_reg[i,j] = output_regincond$RegulationInCondition[intersect(which(output_regincond$RegulationInCondition$gene==genes[i]),which(output_regincond$RegulationInCondition$regulator==regulators[j])),4]
        }
      }
    }
    
    df<-data.frame(genes=rep(rownames(gen_reg), times=ncol(gen_reg)),
                   regulators=rep(colnames(gen_reg),each = nrow(gen_reg)),
                   value=as.vector(gen_reg),
                   color= ifelse(as.vector(gen_reg)>0, '#5577FF',ifelse(as.vector(gen_reg)<0, '#FF7755','#FFFFFF')))
    df$color[is.na(df$color)] <- '#aaaaaa'
    
    ggplot2::ggplot(data = df, aes(x = regulators, y = genes, fill = color)) +
      geom_tile(color = "black",lwd = 0.5,  linetype = 1) +
      scale_fill_manual(values = c("#5577FF", "#aaaaaa", "#FF7755",  "#FFFFFF"), 
                        name = "Legend",
                        labels = c('Activator','Potential','Repressor','Not potential')) +
      labs(title = paste0("Gene - Global Regulators \n correlation plot in ",colnames(output_regincond$RegulationInCondition)[4]," condition  \n"), 
           x = "Global regulators", y = "Genes") +
      theme(plot.title = element_text(hjust = 0.5, colour = "black"), 
            axis.title.x = element_text(face="bold", colour="black", size = 8),
            axis.title.y = element_text(face="bold", colour="black", size = 8),
            axis.text.x = element_text(size = 6, angle = 45, hjust = 1),  # Adjust size for x-axis text
            axis.text.y = element_text(size = 6),
            legend.title = element_text(face="bold", colour="black", size = 10),legend.position = 'right') 
    
  }
  
}

## Network creation -------

#' network_more
#'
#' \code{network_more} Function to be applied to RegulationPerConidtion function output.
#' 
#' @param output_regpcond Output object of RegulationPerCondition applied to MORE main function.
#' @param cytoscape TRUE for plotting the network in Cytoscape. FALSE to plot the network in R. 
#' @param group1 Name of the group to take as reference in the differential network creation.
#' @param group2 Name of the group to compare to the reference in the differential network creation. 
#' 
#' @return Plot of the network induced from more.
#' @export



network_more <- function(output_regpcond, cytoscape = TRUE, group1 = NULL, group2 = NULL) {
  
  create_graph <- function(df) {
    #Remove rows with 0 coef
    df <- df[df[,4] != 0, ]
    mygraph <- igraph::graph_from_data_frame(df, directed = FALSE)
    mygraph <- igraph::simplify(mygraph)
    #Add atributtes to the edge and create df with the omic
    odf <- unique(df[, c(2, 3)])
    odf <- rbind(odf, data.frame('regulator' = unique(df$gene), 'omic' = rep('gene', length(unique(df$gene)))))
    rownames(odf) <- odf$regulator
    
    mygraph <- igraph::set_vertex_attr(mygraph, 'omic', index = igraph::V(mygraph), value = odf[igraph::V(mygraph)$name,]$omic)
    mygraph <- igraph::set_edge_attr(mygraph, 'sign', index = igraph::E(mygraph), value = df[, 4])
    
    return(list('mygraph'=mygraph,'df'=df,'odf'=odf))
  }
  
  create_network <- function(mygraph, df,odf, prefix, group_names,diff) {
    RCy3::setEdgeLineWidthDefault(10)
    cy_network <- RCy3::createNetworkFromIgraph(mygraph, paste0(prefix, group_names))
    
    edge_names <- gsub(" \\(interacts with\\) ", "--", RCy3::getAllEdges(cy_network))
    edges_graph <- apply(df[, c(1, 2)], 1, function(row) paste(row, collapse = "--"))
    #Order modified vector based on the order_index
    order_index <- match(edge_names, edges_graph)
    edge_colors <- ifelse(df[order_index, 4] > 0, '#5577FF', '#FF3333')
    RCy3::setEdgeColorBypass(network = cy_network, edge.names = RCy3::getAllEdges(cy_network), edge_colors)
    
    if(diff){
      edge_lines<-ifelse(df[order_index,5] == 0, 'SOLID', ifelse(df[order_index,5] == 1,'ZIGZAG','DOT'))
      RCy3::setEdgeLineStyleBypass(network= cy_network, edge.names = RCy3::getAllEdges(cy_network),  edge_lines)
    }
    #Set node color and generate a color palette
    omic_c <- factor(odf[RCy3::getAllNodes(cy_network), ]$omic)
    num_unique <- length(unique(omic_c))
    color_palette <- colorbiostat(num_unique)
    node_colors <- color_palette[as.integer(omic_c)]
    
    nshaps <-setdiff(RCy3::getNodeShapes(), c("TRIANGLE", "DIAMOND","RECTANGLE"))[1:num_unique]
    node_shapes <- nshaps[as.integer(omic_c)]
    if('TF'%in% omic_c){
      i=grep('TF', omic_c)
      node_shapes[i]<-'TRIANGLE'
    }
    if('miRNA'%in% omic_c){
      i=grep('miRNA', omic_c)
      node_shapes[i]<-'DIAMOND'
    }
    if('gene'%in% omic_c){
      i=grep('gene', omic_c)
      node_shapes[i]<-'RECTANGLE'
    }
    RCy3::setNodeColorBypass(network = cy_network, node.names = RCy3::getAllNodes(cy_network), node_colors)
    RCy3::setNodeShapeBypass(network = cy_network, node.names = RCy3::getAllNodes(cy_network), node_shapes)
  }
  
  if (cytoscape) {
    
    if (is.null(group1) && is.null(group2)) {
      
      ngroups <- grep('Group', colnames(output_regpcond))
      #Create as many networks as groups
      for (i in 1:length(ngroups)) {
        #Data.frame of that network
        df <- output_regpcond[, c(1, 2, 3, ngroups[i])]
        my_graph <- create_graph(df)
        create_network(my_graph$mygraph, my_graph$df,my_graph$odf, 'mynet', colnames(output_regpcond)[ngroups[i]],diff = FALSE)
      }
      
    } else {
      #Look for the groups to consider
      gr1 <- grep(group1, colnames(output_regpcond))
      gr2 <- grep(group2, colnames(output_regpcond))
      
      if (length(gr1) != 1 || length(gr2) != 1 || gr1 == gr2){stop("ERROR: group1 and group2 should be different names of groups to compare")}
      #Create the differential coefficient and the indicator of sign change
      df <- output_regpcond[, c(1,2,3,gr2,gr1)]
      df[, 6] = df[, 4] - df[, 5]
      df[, 7] = ifelse(df[,4]==0 | df[,5]==0,2, ifelse(sign(df[, 6]) == sign(df[, 5]), 0, 1))
      df <- df[, -c(4, 5)]
      names(df)[5] = 'line'
      
      my_graph <- create_graph(df)
      create_network(my_graph$mygraph,my_graph$df,my_graph$odf, 'mynet', paste0(group2,'-', group1),diff = TRUE)
    }
    
  } else {
    
    if (is.null(group1) && is.null(group2)) {
      
      ngroups <- grep('Group', colnames(output_regpcond))
      #Create as many networks as groups
      for (i in 1:length(ngroups)) {
        #Data.frame of that network
        df <- output_regpcond[, c(1, 2, 3, ngroups[i])]
        my_graph <- create_graph(df)
        
        igraph::E(my_graph$mygraph)$sign<-my_graph$df[,4]
        
        igraph::plot.igraph(my_graph$mygraph, vertex.label.cex = 0.3, vertex.size = 3, 
             vertex.color = as.factor(igraph::V(my_graph$mygraph)$omic), 
             edge.color = ifelse(igraph::E(my_graph$mygraph)$sign > 0, "blue", "red"))
        
        igraph::write_graph(my_graph$mygraph, format = 'gml', file = paste0('mynet', colnames(output_regpcond)[ngroups[i]], '.gml'))
      }
      
    } else {
      #Look for the groups to consider
      gr1 <- grep(group1, colnames(output_regpcond))
      gr2 <- grep(group2, colnames(output_regpcond))
      
      if (length(gr1) != 1 || length(gr2) != 1 || gr1 == gr2){stop("ERROR: group1 and group2 should be different names of groups to compare")}
      #Create the differential coefficient and the indicator of sign change
      df <- output_regpcond[, c(1,2,3,gr2,gr1)]
      df[, 6] = df[, 4] - df[, 5]
      df[, 7] = ifelse(df[,4]==0 | df[,5]==0, 2, ifelse(sign(df[, 6]) == sign(df[, 5]), 0, 1))
      df <- df[, -c(4, 5)]
      names(df)[5] = 'line'
      
      my_graph <- create_graph(df)
      
      igraph::E(my_graph$mygraph)$sign<-my_graph$df[,4]
      my_graph$mygraph<-igraph::set.edge.attribute(my_graph$mygraph, 'line', index = igraph::E(my_graph$mygraph), value = my_graph$df[,5])
      igraph::E(my_graph$mygraph)$line<-my_graph$df[,5]
      
      igraph::plot.igraph(my_graph$mygraph, vertex.label.cex = 0.3, vertex.size = 3, 
           vertex.color = as.factor(igraph::V(my_graph$mygraph)$omic), 
           edge.color = ifelse(igraph::E(my_graph$mygraph)$sign > 0, "blue", "red"), 
           edge.lty = ifelse(igraph::E(my_graph$mygraph)$line == 0, "solid", "dashed"))
      
      igraph::write_graph(my_graph$mygraph, format = 'gml', file = paste0('mynet', group2, '-', group1, '.gml'))
    }
  }
}

## Downstream analysis -------


ReguEnrich1regu1function = function(term, test, notTest, annotation) {
  annotTest = length(intersect(test, annotation[annotation[,2] == term,1]))
  if ((annotTest) > 0) {
    annotNOTtest = length(intersect(notTest, annotation[annotation[,2] == term,1]))
    mytest = matrix(c(annotTest, length(test)-annotTest, annotNOTtest, length(notTest)-annotNOTtest), ncol = 2)
    resultat = c(term, annotTest, length(test), annotNOTtest, length(notTest),
                 fisher.test(mytest, alternative = "greater")$p.value)
    names(resultat) = c("term", "annotTest", "test", "annotNotTest", "notTest", "pval")
  } else {
    resultat = c(term, 0, 0, 0, 0, 100)
    names(resultat) = c("term", "annotTest", "test", "annotNotTest", "notTest", "pval")
  }
  return(resultat)
}

ReguEnrich1regu = function(test, notTest, annotation, p.adjust.method = "fdr") {
  annot2test = unique(annotation[,2])
  resultat = t(sapply(annot2test, ReguEnrich1regu1function, test = test, notTest = notTest, annotation = annotation))
  resultat = resultat[-which(as.numeric(resultat[,"annotTest"]) == 0),]
  return(data.frame(resultat,
                    "adjPval" = p.adjust(as.numeric(resultat[,"pval"]), method = p.adjust.method),
                    stringsAsFactors = F))
}


#' ORA_more
#'
#' \code{ORA_more} Function to be applied to RegulationInCondition function output.
#' 
#' @param output Output object of running more
#' @param output_globregincond Output object of running RegulationInCondition function
#' @param annotation Annotation matrix with genes in the first column, GO terms in the second 
#' @param alpha The adjusted pvalue cutoff to consider
#' @param p.adjust.method One of holm, hochberg, hommel, bonferroni, BH, BY, fdr or none
#' 
#' @return Plot of the network induced from more.
#' @export

ORA_more = function(output, output_globregincond, annotation, alpha = 0.05,
                    p.adjust.method = "fdr") {
  
  # output: Output object of running more function 
  # output_globregincond: Output object of running RegulationInCondition function
  # annotation: Annotation matrix with genes in the first column, GO terms in the second and a description in the third
  
  #Take the GlobalRegulators to which we want to apply the object and against which we test it
  regulators = output_globregincond$GlobalRegulators
  reference = rownames(output$arguments$GeneExpression)
  
  #Take only the ones that were DE
  annotation = annotation[annotation[,1] %in% reference,]
  annotDescr = unique(annotation[,2:3])
  rownames(annotDescr) = annotDescr[,1]
  myresults = vector("list", length = length(regulators))
  names(myresults) = regulators
  for (rrr in regulators) {
    print(rrr)
    test = unique(as.character(output_globregincond$RegulationInCondition[which(output_globregincond$RegulationInCondition[,"regulator"] == rrr),"gene"]))
    notTest = setdiff(reference, test)
    resuRegu = ReguEnrich1regu(test, notTest, annotation, p.adjust.method = p.adjust.method)
    resuRegu = resuRegu[which(resuRegu$adjPval < alpha),]
    if (nrow(resuRegu) > 0) {
      resuRegu = data.frame("regulator" = rrr, "termDescr" = annotDescr[resuRegu[,1],2],
                            resuRegu, stringsAsFactors = FALSE)
    } else {
      
      resuRegu = NULL
    }
    myresults[[rrr]] = resuRegu
    rm(resuRegu); gc()
  }
  return(myresults)
}

library(clusterProfiler)

#' GSEA_more
#'
#' \code{GSEA_more} Function to be applied to RegulationInCondition function output.
#' 
#' @param output_globregincond Output object of running RegulationInCondition function
#' @param output_globregincond2 Output object of running RegulationInCondition function for other group different to the previous. By default, NULL.
#' @param annotation Annotation matrix with genes in the first column, GO terms in the second 
#' @param alpha The adjusted pvalue cutoff to consider
#' @param p.adjust.method One of holm, hochberg, hommel, bonferroni, BH, BY, fdr or none
#' 
#' @return Plot of the network induced from more.
#' @export

GSEA_more<-function(output_globregincond, output_globregincond2 = NULL, annotation, alpha = 0.05,
                    p.adjust.method = "fdr"){
  
  #output_globregincond: Output object of running RegulationInCondition
  #output_globregincond2: Output object of running RegulationInCondition for other group to which compare the reference. By default, NULL.
  #annotation: Annotation matrix with genes in the first column, GO terms in the second 
  #alpha: The adjusted pvalue cutoff to consider
  #p.adjust.method: One of holm, hochberg, hommel, bonferroni, BH, BY, fdr or none
  term2gene_bp<-annotation[,c(,1)]
  if(is.null(output_globregincond2)){
    #Store the genes in decreasing order

    geneList<-table(output_globregincond$RegulationInCondition$gene)
    geneList = sort(geneList, decreasing = TRUE)
    
    selected_genes <- names(geneList)
    counts <- as.numeric(geneList)
    geneList <- setNames(counts, selected_genes)
    
  } else{
    
    geneList<-as.data.frame(table(output_globregincond$RegulationInCondition$gene))
    
    geneList2<-as.data.frame(table(output_globregincond2$RegulationInCondition$gene))
    
    #Create a data frame 
    all_genes <- union(geneList[,1], geneList2[,1])
    merged_df <- data.frame(Gene = all_genes, Count1 = 0, Count2 = 0, stringsAsFactors = FALSE)
    merged_df$Count1[merged_df$Gene %in% geneList[,1]] <- geneList[, 2]
    merged_df$Count2[merged_df$Gene %in% geneList2[,1]] <- geneList2[, 2]
    rownames(merged_df)<-merged_df[,1]
    merged_df<-merged_df[,-1]
    #Create the score for the GSEA
    merged_df[,3]<-merged_df[,2]-merged_df[,1]
    
    geneList <- setNames(merged_df[,3], rownames(merged_df))
    geneList <- sort(geneList, decreasing = TRUE)
    
  }
  
  y <- GSEA(geneList = geneList, TERM2GENE = term2gene_bp, pvalueCutoff = alpha, pAdjustMethod = p.adjust.method)
  
  return(y)
}
