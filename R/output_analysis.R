
#########################################
#### Functions for output analysis
#########################################

## By Sonia, Monica, Maider
## 12-March-2024

# Function to obtain all significant pairs targetF-regulator per omic --------

# For only 1 targetF
GetPairs1targetFAllReg = function (targetF, output) {
  
  if(output$arguments$method=='MLR'|| output$arguments$method=='ISGL'){
    
    reguSignif = output$ResultsPerTargetF[[targetF]]$relevantRegulators
    
    if (is.null(reguSignif)) {  # NO significant regulators
      return (NULL)
      
    } else {  # Significant regulators
      
      reguSignif = output$ResultsPerTargetF[[targetF]]$allRegulators[reguSignif,]
      reguSignif = reguSignif[,c("targetF", "regulator", "omic", "area", "filter")]
      return (reguSignif)
    }
    
  }
  if(output$arguments$method=='PLS1' || output$arguments$method=='PLS2'){
    
    reguSignif = output$ResultsPerTargetF[[targetF]]$significantRegulators
    
    if (is.null(reguSignif)) {  # NO significant regulators
      return (NULL)
      
    } else {  # Significant regulators
      
      reguSignif = output$ResultsPerTargetF[[targetF]]$allRegulators[reguSignif,]
      reguSignif = reguSignif[,c("targetF", "regulator", "omic", "area", "filter")]
      return (reguSignif)
    }
    
  }
  
}


# For all targetFs
GetPairstargetFRegulator = function (targetFs = NULL, output) {
  
  if (is.null(targetFs)) targetFs = rownames(output$GlobalSummary$ReguPerTargetF)
  
  myresults = do.call("rbind", lapply(targetFs, GetPairs1targetFAllReg, output))
  
  #   colnames(myresults) = c("targetF", "regulator", "omic", "area")
  return(myresults)
}

#Make the coefficients comparable

ComparableBetas = function(myresults, output){
  
  for(k in unique(myresults[,'targetF'])){
    sd_targetF = sd(output$arguments$targetData[k,])
    myresults[myresults[,"targetF"] == k, grep('Group',colnames(myresults)) ] = myresults[myresults[,"targetF"] == k, grep('Group',colnames(myresults)) ]/sd_targetF
  }
  return(myresults)
}

#' RegulationPerCondition
#'
#' \code{RegulationPerCondition} Function to be applied to \link{more} main function output.
#' 
#' @param output Output object of MORE main function.
#' @param filterR2 Filters out genes with less R2 than the specified. By default, 0.
#' 
#' @return Summary table containing all the relevant/significant regulators. Moreover, it provides the regression coefficient that relates the targetF and the regulator for each experimental condition after testing if this coefficient is relevant/significant or not.
#'
#' @export


RegulationPerCondition = function(output, filterR2 = 0){
  
  #Filter to only those targetF of interest
  filtered_targetF = rownames(output$GlobalSummary$GoodnessOfFit)[which(output$GlobalSummary$GoodnessOfFit[,1]>filterR2)]
  output$ResultsPerTargetF = output$ResultsPerTargetF[names(output$ResultsPerTargetF) %in% filtered_targetF]
  
  # output: results of the getMLR/getPLS function.
  method = output$arguments$method
  #Add a progressbar
  pb <- txtProgressBar(min = 0, max = length(rownames(output$GlobalSummary$ReguPerTargetF)), style = 3)
  if(method =='MLR'|| method=='ISGL'){
    design = output$arguments$finaldesign
    Group = output$arguments$groups
    
    # Creo a partir de la funcion que ya estaba hecha (linea 1657) la tabla y le anyado los huecos en blanco y cambio el nombre a "representative".
    targetFs = rownames(output$GlobalSummary$ReguPerTargetF)
    myresults = do.call("rbind", lapply(targetFs, GetPairs1targetFAllReg, output))
    colnames(myresults) = c(colnames(myresults)[1:4], "representative")
    myresults[myresults[, "representative"] == "Model", "representative"] = ""
    if (is.null(design)){
      
      # Anyado la columna de coeficientes.
      coeffs = matrix(1, nrow(myresults), 1)
      colnames(coeffs) = "coefficients"
      rownames(coeffs) = rownames(myresults)
      myresults = cbind(myresults, coeffs)
      myresults[grep("_N", myresults[, "representative"]), "coefficients"] = -1  # Para cambiar el signo si pertenece al grupo de correlacionados negativamente
      
      for(k in unique(myresults[,"targetF"])){
        setTxtProgressBar(pb, value = which(names(output$ResultsPerTargetF)==k))
        # Posicion y reguladores que son representantes.
        counts = grep("_R", myresults[myresults[,"targetF"] == k, "representative"]) # positions of representatives of mc
        representatives = myresults[myresults[,"targetF"] == k, "regulator"][counts]      # Devuelve el nombre real de los reguladores representantes
        omic.representative = myresults[myresults[,"targetF"] == k, c("regulator", "representative")][counts,]   # Columna Regulator y Representative
        
        # Necesito el if, si no da error. En caso de entrar, elimino las coletillas para que sea mas sencillo buscar y asignar el representante
        if(length(representatives) != 0){
          norow.nulls = which(myresults[myresults[,"targetF"] == k, "representative"] != "")
          myresults[myresults[,"targetF"] == k, "representative"][norow.nulls] = sub("_P", "", myresults[myresults[,"targetF"] == k, "representative"][norow.nulls])
          myresults[myresults[,"targetF"] == k, "representative"][norow.nulls] = sub("_N", "", myresults[myresults[,"targetF"] == k, "representative"][norow.nulls])
          myresults[myresults[,"targetF"] == k, "representative"][norow.nulls] = sub("_R", "", myresults[myresults[,"targetF"] == k, "representative"][norow.nulls])
          
          for(i in 1:length(representatives)){
            # Aquellos que se llamen igual "omica_mc(numero)", se les asignara el representante
            reg.rep = myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == representatives[i], "representative"]
            myresults[myresults[,"targetF"] == k & myresults[,"representative"] == reg.rep, "representative"] = representatives[i]
          }
          
          # Reguladores significativos del MLR. Pongo gsub() porque haciendo pruebas he visto que hay reguladores que se nombran `nombre regulador`.
          # Las comitas haran que no pueda encontrar el regulador
          # en la tabla. Sin embargo, creo sign.glm para manterner las comitas y poder acceder a la tabla de coeficientes
          significatives = gsub("`", "", names(output$ResultsPerTargetF[[k]]$coefficients[2:nrow(output$ResultsPerTargetF[[k]]$coefficients), 1]))
          sign.glm = names(output$ResultsPerTargetF[[k]]$coefficients[2:nrow(output$ResultsPerTargetF[[k]]$coefficients), 1])
          
          for(i in 1:length(significatives)){
            if(any(significatives[i] == omic.representative[,2])){
              # index.regul: para saber que regulador es el representante y asi todos los que tengan su nombre en la columna "representative" tendran su coeficiente del modelo MLR.
              index.regul = rownames(omic.representative)[which(omic.representative[,2] == significatives[i])]
              PN = myresults[myresults[,"targetF"] == k & myresults[,"representative"] == index.regul, "coefficients"]                        # Sera 1 o -1, segun tenga "_P" o "_N"
              myresults[myresults[,"targetF"] == k & myresults[,"representative"] == index.regul, "coefficients"] = PN*output$ResultsPerTargetF[[k]]$coefficients[sign.glm[i], 1]                                                    # Tendra signo de la tabla si es "_P" y signo opuesto si es "_N".
            } else {
              # En caso de no pertenecer a un grupo de reguladores correlacionados, cogera su coeficiente de la tabla y lo asignara a la posicion correspondiente
              myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == significatives[i], "coefficients"] = output$ResultsPerTargetF[[k]]$coefficients[sign.glm[i], 1]
            }
          }
          
        } else {
          # Si no presenta grupo de reguladores correlacionados, simplemente sacara los coeficientes de la tabla "coefficients"
          myresults[myresults[,"targetF"] == k, "coefficients"] = output$ResultsPerTargetF[[k]]$coefficients[2:nrow(output$ResultsPerTargetF[[k]]$coefficients), 1]
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
      
      for(k in unique(myresults[,"targetF"])){
        setTxtProgressBar(pb, value = which(names(output$ResultsPerTargetF)==k))
        significant.regulators = output$ResultsPerTargetF[[k]]$relevantRegulators                    # Reguladores significativos.
        if(method =='MLR'){
          model.variables = gsub("`", "", rownames(output$ResultsPerTargetF[[k]]$coefficients))[-1]       # Reguladores e interacciones en el modelo.
          kc = 2
        } else{
          model.variables = gsub("`", "", rownames(output$ResultsPerTargetF[[k]]$coefficients))       # Reguladores e interacciones en el modelo.
          kc = 1
        }
        
        # Cojo las interacciones y creo objetos que contengan los reguladores que aparecen con interaccion, solas o ambas.
        interactions.model = gsub("`", "", rownames(output$ResultsPerTargetF[[k]]$coefficients)[grep(":", rownames(output$ResultsPerTargetF[[k]]$coefficients))])
        
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
        
        for(j in kc:nrow(output$ResultsPerTargetF[[k]]$coefficients)){
          regul = unlist(strsplit(gsub("`", "", rownames(output$ResultsPerTargetF[[k]]$coefficients)[j]), ":"))
          
          # Evaluo en que conjunto se encuentra el regulador correspondiente y segun eso asigno el coeficiente o sumo el nuevo coeficiente a lo que ya habia en esa posicion.
          if(any(regul %in% variables.only)){
            if(any(regul %in% significant.regulators)){
              myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == regul, c(names.groups)] = output$ResultsPerTargetF[[k]]$coefficients[j,]
            } else {
              myresults[myresults[,"targetF"] == k & myresults[,"representative"] == regul, c(names.groups)] = output$ResultsPerTargetF[[k]]$coefficients[j,]
            }
          }
          
          if(any(regul %in% variables.inter)){
            if(any(regul %in% significant.regulators)){
              myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == regul[2], regul[1]] = myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == regul[2], regul[1]] + output$ResultsPerTargetF[[k]]$coefficients[j,]
            } else {
              myresults[myresults[,"targetF"] == k & myresults[,"representative"] == regul[2], regul[1]] = myresults[myresults[,"targetF"] == k & myresults[,"representative"] == regul[2], regul[1]] + output$ResultsPerTargetF[[k]]$coefficients[j,]
            }
          }
          
          if(any(regul %in% variables.inter.only)){
            if(any(regul %in% significant.regulators)){
              if(length(regul) == 1){
                myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == regul, c(names.groups)] = myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == regul, c(names.groups)] + output$ResultsPerTargetF[[k]]$coefficients[j,]
              } else {
                myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == regul[2], regul[1]] = myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == regul[2], regul[1]] + output$ResultsPerTargetF[[k]]$coefficients[j,]
              }
            } else {
              if(length(regul) == 1){
                myresults[myresults[,"targetF"] == k & myresults[,"representative"] == regul, c(names.groups)] = myresults[myresults[,"targetF"] == k & myresults[,"representative"] == regul, c(names.groups)] + output$ResultsPerTargetF[[k]]$coefficients[j,]
              } else {
                myresults[myresults[,"targetF"] == k & myresults[,"representative"] == regul[2], regul[1]] = myresults[myresults[,"targetF"] == k & myresults[,"representative"] == regul[2], regul[1]] + output$ResultsPerTargetF[[k]]$coefficients[j,]
              }
            }
          }
        }
        
        # Veo si hay representantes, en caso de haberlos asignara la misma fila del representante a los reguladores que acaben en "_P" y el opuesto a los que acaban en "_N".
        countsR = grep("_R", myresults[myresults[,"targetF"] == k, "representative"])
        
        if(length(countsR) != 0){
          countsR = myresults[myresults[,"targetF"] == k, 5:ncol(myresults)][countsR,]
          
          # Para los correlacionados positivamente: mete la misma fila de coeficientes del representante.
          countsP = countsR
          countsP[,"representative"] = sub("_R", "", countsP[,"representative"])
          countsP[,"representative"] = paste(countsP[,"representative"], "_P", sep = "")
          
          for(l in 1:nrow(countsP)){
            myresults[myresults[,"targetF"] == k & myresults[,"representative"] == countsP[l,"representative"], 6:ncol(myresults)] = countsP[l,2:ncol(countsP)]
          }
          
          # Para los correlacionados negativamente: mete la fila opuesta de coeficientes del representante.
          countsN = countsR
          countsN[,"representative"] = sub("_R", "", countsN[,"representative"])
          countsN[,"representative"] = paste(countsN[,"representative"], "_N", sep = "")
          
          for(l in 1:nrow(countsN)){
            myresults[myresults[,"targetF"] == k & myresults[,"representative"] == countsN[l,"representative"], 6:ncol(myresults)] = -countsN[l,2:ncol(countsN)]
          }
          
          counts = grep("_R", myresults[myresults[,"targetF"] == k, "representative"])
          representatives = myresults[myresults[,"targetF"] == k, "regulator"][counts]                                    # Devuelve el nombre real de los reguladores representantes
          omic.representative = myresults[myresults[,"targetF"] == k, c("regulator", "representative")][counts,]          # Columna Regulator y Representative
          
          # Necesito el if, sino da error. En caso de entrar, elimino las coletillas para que sea mas sencillo buscar y asignar el representante.
          if(length(representatives) != 0){
            norow.nulls = which(myresults[myresults[,"targetF"] == k, "representative"] != "")
            myresults[myresults[,"targetF"] == k, "representative"][norow.nulls] = sub("_P", "", myresults[myresults[,"targetF"] == k, "representative"][norow.nulls])
            myresults[myresults[,"targetF"] == k, "representative"][norow.nulls] = sub("_N", "", myresults[myresults[,"targetF"] == k, "representative"][norow.nulls])
            myresults[myresults[,"targetF"] == k, "representative"][norow.nulls] = sub("_R", "", myresults[myresults[,"targetF"] == k, "representative"][norow.nulls])
            
            for(i in 1:length(representatives)){
              # Aquellos que se llamen igual "omica_mc(numero)", se les asignara el representante.
              reg.rep = myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == representatives[i], "representative"]
              myresults[myresults[,"targetF"] == k & myresults[,"representative"] == reg.rep, "representative"] = representatives[i]
            }
          }
        }
      }
    }
    myresults = ComparableBetas(myresults, output)
    myresults[,6:ncol(myresults)] = signif(myresults[,6:ncol(myresults)], digits = 4) # Para que no salgan los numeros en diferentes notaciones
    
  }
  if(method=='PLS1' || method=='PLS2'){
    design = output$arguments$finaldesign
    Group = output$arguments$groups
    
    # Creo a partir de la funcion que ya estaba hecha (linea 1657) la tabla y le anyado los huecos en blanco y cambio el nombre a "representative".
    targetFs = rownames(output$GlobalSummary$ReguPerTargetF)
    myresults = do.call("rbind", lapply(targetFs, GetPairs1targetFAllReg, output))
    colnames(myresults) = c(colnames(myresults)[1:4], "representative")
    myresults[myresults[, "representative"] == "Model", "representative"] = ""
    
    if (is.null(design)){
      
      # Anyado la columna de coeficientes.
      coeffs = matrix(1, nrow(myresults), 1)
      colnames(coeffs) = "coefficients"
      rownames(coeffs) = rownames(myresults)
      myresults = cbind(myresults, coeffs)
      myresults[grep("_N", myresults[, "representative"]), "coefficients"] = -1  # Para cambiar el signo si pertenece al grupo de correlacionados negativamente
      
      for(k in unique(myresults[,"targetF"])){
        setTxtProgressBar(pb, value = which(names(output$ResultsPerTargetF)==k))
        # Posicion y reguladores que son representantes.
        counts = grep("_R", myresults[myresults[,"targetF"] == k, "representative"]) # positions of representatives of mc
        representatives = myresults[myresults[,"targetF"] == k, "regulator"][counts]      # Devuelve el nombre real de los reguladores representantes
        omic.representative = myresults[myresults[,"targetF"] == k, c("regulator", "representative")][counts,]   # Columna Regulator y Representative
        
        # Necesito el if, si no da error. En caso de entrar, elimino las coletillas para que sea mas sencillo buscar y asignar el representante
        if(length(representatives) != 0){
          norow.nulls = which(myresults[myresults[,"targetF"] == k, "representative"] != "")
          myresults[myresults[,"targetF"] == k, "representative"][norow.nulls] = sub("_P", "", myresults[myresults[,"targetF"] == k, "representative"][norow.nulls])
          myresults[myresults[,"targetF"] == k, "representative"][norow.nulls] = sub("_N", "", myresults[myresults[,"targetF"] == k, "representative"][norow.nulls])
          myresults[myresults[,"targetF"] == k, "representative"][norow.nulls] = sub("_R", "", myresults[myresults[,"targetF"] == k, "representative"][norow.nulls])
          
          for(i in 1:length(representatives)){
            # Aquellos que se llamen igual "omica_mc(numero)", se les asignara el representante
            reg.rep = myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == representatives[i], "representative"]
            myresults[myresults[,"targetF"] == k & myresults[,"representative"] == reg.rep, "representative"] = representatives[i]
          }
          
          # Reguladores significativos del MLR. Pongo gsub() porque haciendo pruebas he visto que hay reguladores que se nombran `nombre regulador`.
          # Las comitas haran que no pueda encontrar el regulador
          # en la tabla. Sin embargo, creo sign.glm para manterner las comitas y poder acceder a la tabla de coeficientes
          significatives = gsub("`", "", names(output$ResultsPerTargetF[[k]]$coefficients[1:nrow(output$ResultsPerTargetF[[k]]$coefficients), 1]))
          sign.glm = names(output$ResultsPerTargetF[[k]]$coefficients[1:nrow(output$ResultsPerTargetF[[k]]$coefficients), 1])
          
          for(i in 1:length(significatives)){
            if(any(significatives[i] == omic.representative[,2])){
              # index.regul: para saber que regulador es el representante y asi todos los que tengan su nombre en la columna "representative" tendran su coeficiente del modelo MLR.
              index.regul = rownames(omic.representative)[which(omic.representative[,2] == significatives[i])]
              PN = myresults[myresults[,"targetF"] == k & myresults[,"representative"] == index.regul, "coefficients"]                        # Sera 1 o -1, segun tenga "_P" o "_N"
              myresults[myresults[,"targetF"] == k & myresults[,"representative"] == index.regul, "coefficients"] = PN*output$ResultsPerTargetF[[k]]$coefficients[sign.glm[i], 1]                                                    # Tendra signo de la tabla si es "_P" y signo opuesto si es "_N".
            } else {
              # En caso de no pertenecer a un grupo de reguladores correlacionados, cogera su coeficiente de la tabla y lo asignara a la posicion correspondiente
              myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == significatives[i], "coefficients"] = output$ResultsPerTargetF[[k]]$coefficients[sign.glm[i], 1]
            }
          }
          
        } else {
          # Si no presenta grupo de reguladores correlacionados, simplemente sacara los coeficientes de la tabla "coefficients"
          myresults[myresults[,"targetF"] == k, "coefficients"] = output$ResultsPerTargetF[[k]]$coefficients[1:nrow(output$ResultsPerTargetF[[k]]$coefficients), 1]
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
      
      for(k in unique(myresults[,"targetF"])){
        
        setTxtProgressBar(pb, value = which(names(output$ResultsPerTargetF)==k))

        significant.regulators = output$ResultsPerTargetF[[k]]$significantRegulators                    # Reguladores significativos.
        model.variables = gsub("`", "", rownames(output$ResultsPerTargetF[[k]]$coefficients))           # Reguladores e interacciones en el modelo.
        
        # Cojo las interacciones y creo objetos que contengan los reguladores que aparecen con interaccion, solas o ambas.
        interactions.model = gsub("`", "", rownames(output$ResultsPerTargetF[[k]]$coefficients)[grep(":", rownames(output$ResultsPerTargetF[[k]]$coefficients))])
        
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
        
        for(j in 1:nrow(output$ResultsPerTargetF[[k]]$coefficients)){
          regul = unlist(strsplit(gsub("`", "", rownames(output$ResultsPerTargetF[[k]]$coefficients)[j]), ":"))
          
          groups = regul[grepl(paste(paste0('Group_',names(table(Group))),collapse='|'), regul)]
          regula = setdiff(regul,groups)
          # Evaluo en que conjunto se encuentra el regulador correspondiente y segun eso asigno el coeficiente o sumo el nuevo coeficiente a lo que ya habia en esa posicion.
          if(any(regul %in% variables.only)){
            if(any(regul %in% significant.regulators)){
              myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == regul, c(names.groups)] = output$ResultsPerTargetF[[k]]$coefficients[j,1]
            } 
          }
          #TO DO: regul[1] se supone que tendría que ser algo sobre los grupos y no un regulador
          if(any(regul %in% variables.inter)){
            if(any(regul %in% significant.regulators)){
              myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == regula, groups] = myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == regula, groups] + output$ResultsPerTargetF[[k]]$coefficients[j,1]
            } 
          }
          
          if(any(regul %in% variables.inter.only)){
            if(any(regul %in% significant.regulators)){
              if(length(regul) == 1){
                myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == regul, c(names.groups)] = myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == regul, c(names.groups)] + output$ResultsPerTargetF[[k]]$coefficients[j,1]
              } else {
                myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == regula, groups] = myresults[myresults[,"targetF"] == k & myresults[,"regulator"] == regula, groups] + output$ResultsPerTargetF[[k]]$coefficients[j,1]
              }
            } 
          }
        }
      }
    }
    myresults = ComparableBetas(myresults, output)
    myresults[,6:ncol(myresults)] = signif(myresults[,6:ncol(myresults)], digits = 4) # Para que no salgan los numeros en diferentes notaciones
    myresults = myresults[,-5,drop=FALSE]
    
  }
  rownames(myresults)=NULL
  close(pb)
  return(myresults)
}

#' FilterRegulationPerCondition
#'
#' \code{FilterRegulationPerCondition} Function to be applied to \link{more} and \link{RegulationPerCondition} function outputs.
#' 
#' @param output Output object of \link{more} function.
#' @param outputRegpcond Output object of \link{RegulationPerCondition} function.
#' @param filterR2 Filters the results for the genes that showed a R2 above the one indicated. By default, 0.
#' 
#' @return RegulationPerCondition matrix filtered to the genes with a R2 above the indicated.
#'
#' @export

FilterRegulationPerCondition <- function(output, outputRegpcond, filterR2 = 0){
  filtered_targetF <- rownames(output$GlobalSummary$GoodnessOfFit)[which(output$GlobalSummary$GoodnessOfFit[, 1] > filterR2)]
  filtered_outputRegpcond <- outputRegpcond[outputRegpcond$targetF %in% filtered_targetF,,drop=FALSE]
  return(filtered_outputRegpcond)
}

#' RegulationInCondition
#'
#' \code{RegulationInCondition} Function to be applied to \link{RegulationPerCondition} function output.
#' 
#' @param outputRegpcond Output object of \link{RegulationPerCondition} function.
#' @param cond Biological condition from which the user wants to summarize the information.
#' 
#' @return List with the hub target features, global regulators and the regulators with their coefficients specific to the requested condition.
#'
#' @export

RegulationInCondition <- function (outputRegpcond, cond){
  
  #Take group column in RegulationPerCondition
  group_col<-grep(cond,colnames(outputRegpcond))
  #Use only the group in which we are interested and remove target features that do not affect that group
  outputRegpcond = outputRegpcond[c(outputRegpcond[,group_col]!=0),c(1,2,3,group_col)]
  
  #Calculate the global regulators
  myreg<-table(outputRegpcond[,2])
  #Calculate third quantile
  q3<-quantile(myreg,0.75)
  if(length(myreg[myreg>q3])<10){
    GlobalRegulators = intersect(names(myreg[rev(tail(order(myreg),10))]), names(myreg[myreg>10]) )
  } else{
    GlobalRegulators = intersect(names(myreg[myreg>q3]), names(myreg[myreg>10]) ) 
  }
  
  #Calculate the hub target features
  myhub<-table(outputRegpcond[,1])
  #Calculate third quantile
  q3<-quantile(myhub,0.75)
  if(length(myhub[myhub>q3])<10){
    HubTargetF = names(myhub[rev(tail(order(myhub),10))])
  } else{
    HubTargetF = names(myhub[myhub>q3])
  }
  
  return(list('GlobalRegulators'=GlobalRegulators, 'HubTargetF'=HubTargetF, 'RegulationInCondition'=outputRegpcond))
  
}


# Plot MLR results --------------------------------------------------------

# Plot 1 targetF versus 1 regulator ------------------------------------------

plotTargetFRegu = function (x.points, targetFValues, reguValues, targetFErrorValues, reguErrorValues, col = c(1,2),
                            xlab = "", yylab = c("right", "left"), type = c(16,16), lwd = c(1,1), main = "",
                            numLines = NULL, x.names = NULL, yleftlim, yrightlim, group_names, smooth, size = 3, breakby = c(0.2,0.3)) {
  
  
  if(smooth){
    #Apply smoothing
    splineTargetF <- smooth.spline(x.points, targetFValues, spar = 0.5)
    splineRegu <- smooth.spline(x.points, reguValues, spar = 0.5)
    # Extract smoothed values
    smoothedTargetF <- splineTargetF$y
    smoothedRegu <- splineRegu$y
  } else{
    smoothedTargetF = NULL
    smoothedRegu = NULL
  }
  
  #Plot
  p = plot.y2(x = x.points, yright = targetFValues, yleft = reguValues, 
              yright2 = smoothedTargetF, yleft2 = smoothedRegu, numLines = numLines,
              xlab = xlab, yylab = yylab, col = col, type = type, lwd = lwd, main = main, 
              yrightErrorValues = targetFErrorValues, yleftErrorValues = reguErrorValues, 
              group_names = group_names,size = size, breakby = breakby)
  
  return(p)
  
}


# Plot Y2 -----------------------------------------------------------------

# Inspired by Ajay Shah Plot 2 time series with different y axes (left and right),
# in https://stat.ethz.ch/pipermail/r-help/2004-March/047775.html)



plot.y2 <- function(x, yright, yleft, 
                    yright2 = NULL, yleft2 = NULL, numLines = NULL,
                    xlab = NULL, yylab = c("", ""), col = c(1, 2),
                    type = c(16, 16), lwd = c(1, 1), main = NULL,
                    #Standard error en series temporales
                    yrightErrorValues = NULL, yleftErrorValues = NULL,
                    group_names = NULL, size = 3, breakby = c(0.2,0.3)) {
  
  # Creating a data frame inside the function
  data <- data.frame(x = x, yright = yright, yleft = yleft)
  
  if (!is.null(yrightErrorValues)) data$yrightErrorValues <- yrightErrorValues
  if (!is.null(yleftErrorValues)) data$yleftErrorValues <- yleftErrorValues
  if (!is.null(yright2)) data$yright2 <- yright2
  if (!is.null(yleft2)) data$yleft2 <- yleft2
  
  p <- ggplot(data, aes(x = x)) +
    
    # Left axis (Primary y-axis)
    geom_point(aes(y = yleft), color = col[2], shape = type[1], size = 2, alpha=0.4) +
    geom_line(aes(y = yleft), color = col[2], linewidth = lwd[1], alpha=0.2) +
    
    # Right axis (Secondary y-axis, rescaled to match yright values)
    geom_point(aes(y = (yright - min(yright)) / (max(yright) - min(yright)) * (max(yleft) - min(yleft)) + min(yleft)), color = col[1], shape = type[2], size = 2, alpha=0.4) +
    geom_line(aes(y = (yright - min(yright)) / (max(yright) - min(yright)) * (max(yleft) - min(yleft)) + min(yleft)), color = col[1], linewidth = lwd[2], alpha=0.2) +
    
    # Labels and Theme
    labs(x = "", y = yylab[2], title = main) +
    scale_y_continuous(
      limits = range(data$yleft),
      breaks = seq(min(data$yleft), max(data$yleft), by = breakby[1]),
      labels =  function(x) sprintf("%.1f", x),
      sec.axis = sec_axis(
        transform = ~ . * (max(data$yright) - min(data$yright)) / (max(data$yleft) - min(data$yleft)) + min(data$yright) - min(data$yleft) * (max(data$yright) - min(data$yright)) / (max(data$yleft) - min(data$yleft)),
        breaks = seq(min(data$yright), max(data$yright), by = breakby[2]),  # Right y-axis ticks (adjust the step size if needed)
        labels = function(x) sprintf("%.1f", x), 
        name = yylab[1]
      )
    ) +
    # Modify x-axis to show custom labels
    scale_x_continuous(
      breaks = data$x, # Set x breaks to be the data points
      labels = rownames(data) # Use combined labels for x-axis
    ) +
    theme_minimal() +
    theme(
      axis.title.y = element_text(color = col[2]), # Set left y-axis label color
      axis.title.y.right = element_text(
        color = col[1], 
        angle = 90    ), # Set right y-axis label color
      panel.grid.major = element_blank(), # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.5),
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = size),
      axis.text.y = element_text(angle = 90, hjust = 1, vjust = 0.5),
      axis.ticks.x = element_line(color = "black", linewidth  = 0.5),
      axis.ticks.y = element_line(color = "black", linewidth  = 0.5),
      plot.title = element_text(hjust = 0.5)
    )
  
  if(!is.null(yright2)){
    # Plot the spline for yleft (Primary y-axis)
    p = p + 
      geom_line(aes(y = yleft2), color = col[2], linewidth = 1) +
      
      # Plot the spline for yright (Secondary y-axis)
      geom_line(aes(y = (yright2 - min(yright2)) / (max(yright2) - min(yright2)) * (max(yleft) - min(yleft)) + min(yleft)), color = col[1], linewidth = lwd[1])
    
  } else{
    p = p + geom_point(aes(y = yleft), color = col[2], shape = type[1], size = 2) +
      geom_line(aes(y = yleft), color = col[2], linewidth = lwd[1]) +
      # Right axis (Secondary y-axis, rescaled to match yright values)
      geom_point(aes(y = (yright - min(yright)) / (max(yright) - min(yright)) * (max(yleft) - min(yleft)) + min(yleft)), color = col[1], shape = type[2], size = 2, alpha=0.4) +
      geom_line(aes(y = (yright - min(yright)) / (max(yright) - min(yright)) * (max(yleft) - min(yleft)) + min(yleft)), color = col[1], linewidth = lwd[2]) 
  }
  
  if (!is.null(group_names)) {
    num_groups = length(group_names)
    group_position = seq(min(data$x), max(data$x), length.out = num_groups + 1)
    group_centers = head(group_position, -1) + diff(group_position) / 2
    
    usr <- range(data$yleft)  # usr[3] es el límite inferior del eje Y, usr[4] es el superior
    y_text <- usr[2] + (usr[2] - usr[1]) * 0.1
    
    group_centers <- group_centers[group_centers >= min(data$x) & group_centers <= max(data$x)]
    
    if (y_text > max(data$yleft)) {
      y_text <- max(data$yleft) - (max(data$yleft) - min(data$yleft)) * 0.005
    }
    
    group_labels_df <- data.frame(
      x = group_centers,
      y = rep(y_text, length(group_centers)),
      label = group_names
    )
    
    p = p + geom_vline(xintercept = numLines, color = "black", linetype = "dashed", linewidth = 0.5) + # Dynamically add the group names above the plot
      geom_text(data = group_labels_df, aes(x = x, y = y, label = label), color = "black", size = 3) 
    
  }
  
  if(!is.null(yrightErrorValues) && !is.null(yleftErrorValues)){
    
    
    y_min = min(targetFValues - yrightErrorValues, na.rm = TRUE)
    y_max = max(targetFValues + yrightErrorValues, na.rm = TRUE)
    
    y_minreg = min(reguValues - yleftErrorValues, na.rm = TRUE)
    y_maxreg = max(reguValues + yleftErrorValues, na.rm = TRUE)
    
    p = p + geom_point(aes(y = yleft), color = col[2], shape = type[1], size = 2) +
      geom_line(aes(y = yleft), color = col[2], linewidth = lwd[1]) +
      
      # Right axis (Secondary y-axis, rescaled to match yright values)
      geom_point(aes(y = (yright - min(yright)) / (max(yright) - min(yright)) * (max(yleft) - min(yleft)) + min(yleft)), color = col[1], shape = type[2], size = 2) +
      geom_line(aes(y = (yright - min(yright)) / (max(yright) - min(yright)) * (max(yleft) - min(yleft)) + min(yleft)), color = col[1], linewidth = lwd[2]) +
      
      # Adjust the axis to include the error value
      geom_errorbar(aes(ymin = yleft - yleftErrorValues, ymax = yleft + yleftErrorValues), 
                    width = 0.2, color = col[2]) +
      geom_errorbar(aes(
        ymin = (yright - yrightErrorValues - min(yright, na.rm = TRUE)) / (max(yright,na.rm = TRUE) - min(yright,na.rm = TRUE)) * (max(yleft,na.rm = TRUE) - min(yleft,na.rm = TRUE)) + min(yleft,na.rm = TRUE),
        ymax = (yright + yrightErrorValues - min(yright, na.rm = TRUE)) / (max(yright,na.rm = TRUE) - min(yright,na.rm = TRUE)) * (max(yleft,na.rm = TRUE) - min(yleft, na.rm = TRUE)) + min(yleft,na.rm = TRUE)
      ), width = 0.2, color = col[1]) +
      
      scale_y_continuous(
        limits = c(y_minreg, y_maxreg),   # Expand the limits
        breaks = seq(y_minreg, y_maxreg, by = breakby[1]),
        labels = function(x) sprintf("%.1f", x),
        sec.axis = sec_axis(
          transform = ~ . * (y_max - y_min) / (y_maxreg - y_minreg) + y_min - y_minreg * (y_max - y_min) / (y_maxreg - y_minreg),
          breaks = seq(y_min, y_max, by = breakby[2]),
          labels = function(x) sprintf("%.1f", x), 
          name = targetF
        )
      )
    
  }
  
  return(p)
}

#' plotMORE
#'
#' \code{plotMORE} Graphical representation of the relationship between target features and regulators.
#' 
#' @param output Output object of MORE main function.
#' 
#' @param targetF ID of the target feature to be plotted.
#' @param regulator ID of the regulator to be plotted. If NULL (default), all regulators of the target feature are plotted.
#' @param simplify If TRUE, a boxplot (if the regulator is binary) or a Scatterplot (otherwise) is plotted to represent the relationship between the target feature and the regulator provided to the function. If FALSE (default), the target feature and the regulator profiles will be plotted.
#' @param reguValues Vector containing the values of a regulator. If NULL (default), these values are taken from the output object as long as they are available. 
#' @param plotPerOmic If TRUE, all the relevant/significant regulators of the given target feature and the same omic are plotted in the same graph. If FALSE (default), each regulator is plotted in a separate plot.
#' @param targetF.col Color to plot the target feature. By default, biostatomic colors. 
#' @param regul.col  Color to plot the regulator. If NULL (default), a color will be assigned by the function, that will be different for each regulatory omic.
#' @param order If TRUE (default), the values in X-axis are ordered.
#' @param xlab Label for the X-axis.
#' @param cont.var  Vector with length equal to the number of observations in data, which optionally may contain the values of the numerical variable (e.g. time) to be plotted on the X-axis. By default, NULL.
#' @param cond2plot Vector or factor indicating the experimental group of each value to represent. If NULL (default), the labels are taken from the experimental design matrix. 
#' @param smooth If TRUE (default), smoothing is applied via splines so that the profiles could be more easily seen.
#' @param size Size of the X-axis labels. By default, 3.
#' @param breakby Range in which values should be displayed in Y-axis for regulators and target features, respectively. By default, c(0.3,0.3).
#' 
#' @return Graphical representation of the relationship between target features and regulators.
#'
#' @export


plotMORE = function(output, targetF, regulator = NULL, simplify = FALSE, reguValues = NULL, plotPerOmic = FALSE,
                    targetF.col = NULL, regu.col = NULL, order = TRUE,
                    xlab = "", cont.var = NULL, cond2plot = NULL,smooth =TRUE, size = 3, breakby = c(0.3,0.3),...) {
  
  if(simplify){
    # from which omic is the regulator?
    SigniReguTargetF = GetPairstargetFRegulator(targetFs = targetF, output = output)
    omic = SigniReguTargetF[SigniReguTargetF[,"regulator"] == regulator,'omic']
    
    if(is.null(output$arguments$condition)){
      
      df = data.frame(
        gen = unlist(output$arguments$targetData[targetF,,drop=TRUE]),
        regulador = unlist(output$arguments$regulatoryData[[omic]][regulator,,drop=TRUE]))
      
      color_palette = colorbiostat(1)
      
      if (output$arguments$omicType[omic]==0){
        
        ggplot2::ggplot(df, aes(x = regulador, y = gen)) +
          theme_minimal()+
          geom_point(color = color_palette, fill = NA, shape=21) +
          labs( x = paste("Regulator\n",regulator), y = paste("Target Feature\n",targetF))
        
      } else {
        
        df$regulador<-factor(df$regulador)
        ggplot2::ggplot(df, aes(x = regulador, y = gen)) + theme_minimal()+
          geom_boxplot(color = color_palette) +  
          scale_x_discrete(labels = c('0','1')) + 
          labs( x = paste("Regulator \n",regulator), y = paste("Target feature\n",targetF))
        
      }
      
    } else{
      
      df = data.frame(
        gen = unlist(output$arguments$targetData[targetF,,drop=TRUE]),
        regulador = unlist(output$arguments$regulatoryData[[omic]][regulator,,drop=TRUE]),
        Group = output$arguments$groups)
      
      num_unique = length(unique(df$Group))+1
      color_palette = colorbiostat(num_unique)
      custom_colors = setNames(color_palette[-1], unique(df$Group))
      
      if (output$arguments$omicType[omic]==0){
        
        ggplot2::ggplot(df, aes(x = regulador, y = gen, color = Group)) +
          geom_point() + scale_color_manual(values = custom_colors)+
          theme_minimal()+
          labs( x = paste("Regulator\n",regulator), y = paste("Target Feature\n",targetF))
        
      } else {
        
        df$regulador<-factor(df$regulador)
        ggplot2::ggplot(df, aes(x = regulador, y = gen,fill=Group)) + theme_minimal()+
          geom_boxplot() + scale_fill_manual(values = custom_colors)+  scale_color_manual(values = custom_colors)+
          scale_x_discrete(labels = c('0','1')) + stat_summary(aes(color = Group),fun='median',geom = 'point', position = position_dodge(width = 0.75))+
          labs( x = paste("Regulator \n",regulator), y = paste("Target feature\n",targetF))
        
      }
      
    }
    
  } else{
    if(output$arguments$method=='MLR'||output$arguments$method=='ISGL'){
      
      return(plotMLR(output, targetF, regulator = regulator, reguValues = reguValues, plotPerOmic = plotPerOmic,
                     targetF.col = targetF.col, regu.col = regu.col, order = order,
                     xlab = xlab, cont.var = cont.var, cond2plot = cond2plot,smooth = smooth,size = size, breakby = breakby, ...))
    }
    
    if(output$arguments$method=='PLS1'||output$arguments$method=='PLS2'){
      
      return(plotPLS(output, targetF, regulator = regulator, reguValues = reguValues, plotPerOmic = plotPerOmic,
                     targetF.col = targetF.col, regu.col = regu.col, order = order,
                     xlab = xlab, cont.var = cont.var, cond2plot = cond2plot,smooth = smooth,size = size, breakby = breakby, ...))
    }
  }
  
  
}


plotMLR = function (MLRoutput, targetF, regulator = NULL, reguValues = NULL, plotPerOmic = FALSE,
                    targetF.col = NULL, regu.col = NULL, order = TRUE,
                    xlab = "", cont.var = NULL, cond2plot = NULL, verbose =TRUE,smooth = TRUE,size = 3, breakby = c(0.2,0.3), ...) {
  
  # Colors for omics
  #omic.col = colors()[c(554,89,111,512,17,586,132,428,601,568,86,390,100,200,300,400,500,10,20,30,40,50,60,70,80,90,150,250,350,450,550)]
  
  if (is.null(regu.col)) {
    omic_names = names(MLRoutput$arguments$regulatoryData)
    num_unique = length(omic_names) + 1
    color_palette = colorbiostat(num_unique)
    if(is.null(targetF.col)){
      targetF.col = color_palette[1]
    }
    any.col = color_palette[-1]
  } else {
    if (length(regu.col) == length(MLRoutput$arguments$regulatoryData)) {
      any.col = regu.col
    } else {
      any.col = rep(regu.col, length(MLRoutput$arguments$regulatoryData))
    }
  }
  names(any.col) = names(MLRoutput$arguments$regulatoryData)
  
  
  # Changing margin
  par(mar = c(5,3,4,3)+0.1)
  
  # Groups to plot
  if (is.null(cond2plot)) {
    if (!is.null(MLRoutput$arguments$condition)) {
      cond2plot = apply(MLRoutput$arguments$condition, 1, paste, collapse = "_")
    }
  }
  
  # Replicates
  if (!is.null(cont.var)) {  # we have continuous variable
    if (!is.null(cond2plot)) { # cont.var + cond2plot
      myreplicates = apply(cbind(cond2plot, cont.var), 1, paste, collapse = "_")
      
    } else {  # only continuous variable is provided
      myreplicates = cont.var
    }
    smooth = FALSE
    order = FALSE
  } else {   # no cont.var
    
    if (!is.null(cond2plot)) { # only cond2plot
      myreplicates = colnames(MLRoutput$arguments$targetData)
    } else {  # nothing
      myreplicates = colnames(MLRoutput$arguments$targetData)
    }
  }
  
  # Cast myreplicates to character
  myreplicates = as.character(myreplicates)
  names(myreplicates) = colnames(MLRoutput$arguments$targetData)
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
  group_names = NULL
  
  # Error values
  getErrorValues = function(realValues, repsInfo) {
    # Disable it
    if (! replicates)
      return(NULL)
    
    out_values = tapply(realValues, repsInfo, function(reps) sd(reps)/sqrt(length(reps)))
    
    return(out_values)
  }
  
  ## REGULATOR = NULL
  
  if (is.null(regulator)) {  ### Plot all regulators for the given target feature
    
    MLRtargetF = MLRoutput$ResultsPerTargetF[[targetF]]
    
    if (is.null(MLRtargetF)) {
      stop(paste("No MLR was obtained for target feature", targetF))
    }
    
    if (is.null(MLRtargetF$relevantRegulators)) { ## No significant regulators
      
      cat("No relevant regulators were found for this target feature.\n")
      
      
    } else
    {  ## Significant regulators:
      
      # Considering multicollinearity
      SigReg = MLRtargetF$allRegulators
      SigReg = SigReg[SigReg$Rel == 1, c("regulator", "omic", "area", "filter")]
      
      SigReg = SigReg[MLRtargetF$relevantRegulators,,drop = FALSE]
      
      cat(paste(nrow(SigReg), "relevant regulators are to be plotted for target feature", targetF)); cat("\n")
      
      # Target feature values
      targetFValues = MLRtargetF$Y$y
      if (order) {
        if(is.null(MLRoutput$arguments$condition)){
          myorder = order(targetFValues)
        } else{
          myorder = order(unlist(MLRoutput$arguments$condition),targetFValues)
          group_names = unique(sort(MLRoutput$arguments$groups))
        }
        cond2plot = cond2plot[myorder]
        targetFValues = targetFValues[myorder]
        myreplicates = myreplicates[myorder]
        myrepliUni = myrepliUni[myorder]
        
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
      }
      
      errorValues = getErrorValues(targetFValues, myreplicates)
      targetFValues = tapply(targetFValues, myreplicates, mean)
      targetFValues = targetFValues[myrepliUni]
      names(targetFValues) = myrepliUni
      
      # X values
      x.points = 1:length(myrepliUni)
      eje = myrepliUni
      
      if (plotPerOmic) { ## All regulators from the same omic in the same plot
        
        myomics = unique(SigReg$omic)
        
        for (oo in myomics) {
          
          SigRegOmic = SigReg[SigReg$omic == oo,]
          
          omicValues = t(MLRoutput$arguments$regulatoryData[[oo]])
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
          
          p = plotTargetFRegu(x.points = x.points, targetFValues = targetFValues, reguValues = reguValues,
                              col = c(targetF.col, mycol), xlab = xlab, yylab = c(targetF, leftlab), type = c(16,16),
                              main = oo, numLines = numLines, x.names = eje,
                              targetFErrorValues = errorValues, reguErrorValues = errorValuesRegu, group_names = group_names, smooth = smooth,
                              size = size, breakby = breakby)
          
          if (nrow(SigRegOmic) > 1) {
            for (i in 2:nrow(SigRegOmic)) {
              
              reguValues = omicValues[, colnames(omicValues) == SigRegOmic$regulator[i]]
              if (order) reguValues = reguValues[myorder]
              errorValuesRegu = getErrorValues(reguValues, myreplicates)
              
              reguValues = tapply(reguValues, myreplicates, mean)
              reguValues = reguValues[myrepliUni]
              names(reguValues) = myrepliUni
              
              regulator_data = data.frame(x=x.points, y = reguValues)
              
              y_min = min(c(p$data$yleft, regulator_data$y), na.rm = TRUE)
              y_max = max(c(p$data$yleft, regulator_data$y), na.rm = TRUE)
              
              if (! is.null(errorValuesRegu)) {
                
                y_min = min(targetFValues - targetFErrorValues, na.rm = TRUE)
                y_max = max(targetFValues + targetFErrorValues, na.rm = TRUE)
                
                y_minreg = min(reguValues - reguErrorValues, na.rm = TRUE)
                y_maxreg = max(reguValues + reguErrorValues, na.rm = TRUE)
                
                p = p + geom_point(data = regulator_data, aes(x = x, y = y), color = mycol, shape = i, size = 2) +
                  geom_line(data = regulator_data, aes(y = y), color = mycol, linewidth = 1, linetype = i) +
                  
                  # Right axis (Secondary y-axis, rescaled to match yright values)
                  geom_point(aes(y = (yright - min(yright)) / (max(yright) - min(yright)) * (max(yleft) - min(yleft)) + min(yleft)), color = col[1], shape = type[2], size = 2) +
                  geom_line(aes(y = (yright - min(yright)) / (max(yright) - min(yright)) * (max(yleft) - min(yleft)) + min(yleft)), color = col[1], linewidth = lwd[2]) +
                  
                  # Adjust the axis to include the error value
                  geom_errorbar(aes(ymin = yleft - yleftErrorValues, ymax = yleft + yleftErrorValues), 
                                width = 0.2, color = col[2]) +
                  geom_errorbar(aes(
                    ymin = (yright - yrightErrorValues - min(yright, na.rm = TRUE)) / (max(yright,na.rm = TRUE) - min(yright,na.rm = TRUE)) * (max(yleft,na.rm = TRUE) - min(yleft,na.rm = TRUE)) + min(yleft,na.rm = TRUE),
                    ymax = (yright + yrightErrorValues - min(yright, na.rm = TRUE)) / (max(yright,na.rm = TRUE) - min(yright,na.rm = TRUE)) * (max(yleft,na.rm = TRUE) - min(yleft, na.rm = TRUE)) + min(yleft,na.rm = TRUE)
                  ), width = 0.2, color = col[1]) +
                  
                  scale_y_continuous(
                    limits = c(y_minreg, y_maxreg),   # Expand the limits
                    breaks = seq(y_minreg, y_maxreg, by = breakby[1]),
                    labels = function(x) sprintf("%.1f", x),
                    sec.axis = sec_axis(
                      transform = ~ . * (y_max - y_min) / (y_maxreg - y_minreg) + y_min - y_minreg * (y_max - y_min) / (y_maxreg - y_minreg),
                      breaks = seq(y_min, y_max, by = breakby[2]),
                      labels = function(x) sprintf("%.1f", x), 
                      name = targetF
                    )
                  )
                
              } else{
                p = p + 
                  scale_y_continuous(
                    limits = c(y_min, y_max),   # Expand the limits
                    breaks = seq(y_min, y_max, by = breakby[1]),
                    labels = function(x) sprintf("%.1f", x),
                    sec.axis = sec_axis(
                      transform = ~ . * (max(p$data$yright) - min(p$data$yright)) / (max(p$data$yleft) - min(p$data$yleft)) + 
                        min(p$data$yright) - min(p$data$yleft) * (max(p$data$yright) - min(p$data$yright)) / (max(p$data$yleft) - min(p$data$yleft)),
                      breaks = seq(min(p$data$yright), max(p$data$yright), by = breakby[2]),
                      labels = function(x) sprintf("%.1f", x), 
                      name = targetF
                    )
                  ) + geom_point(data = regulator_data, aes(x = x, y = y), color = mycol, shape = i, size = 2) +
                  geom_line(data = regulator_data, aes(y = y), color = mycol, linewidth = 1, linetype = i) 
              }
            }
          }
          print(p)
        }
        
      } else {  ## Each regulator in a separate plot
        
        for (rr in SigReg$regulator) {
          
          oo = MLRoutput$ResultsPerTargetF[[targetF]]$allRegulators[rr,"omic"]
          omicValues = t(MLRoutput$arguments$regulatoryData[[oo]])
          reguValues = omicValues[, colnames(omicValues) == rr]
          if (order) reguValues = reguValues[myorder]
          errorValuesRegu = getErrorValues(reguValues, myreplicates)
          
          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni
          
          mycol = any.col[oo]
          
          p = plotTargetFRegu(x.points = x.points, targetFValues = targetFValues, reguValues = reguValues,
                              col = c(targetF.col, mycol), xlab = xlab,
                              yylab = c(targetF, rr), type = c(16,16),
                              main = paste(as.character(SigReg[rr, c("omic", "area")]), collapse = " "),
                              numLines = numLines, x.names = eje,
                              targetFErrorValues = errorValues, reguErrorValues = errorValuesRegu, group_names = group_names, smooth = smooth,
                              size = size, breakby = breakby)
          print(p)
          
        }
        
      }
      
      return(MLRtargetF$allRegulators[MLRtargetF$relevantRegulators, -6])
    }
    
  }
  
  
  ## targetF = NULL
  
  if (is.null(targetF)) {  ### Plot all target features regulated by the regulator
    
    SigniReguTargetF = GetPairstargetFRegulator(targetFs = NULL, output = MLRoutput)
    SigniReguTargetF = SigniReguTargetF[SigniReguTargetF[,"regulator"] == regulator,]
    myomics = SigniReguTargetF[,"omic"]
    myomic = unique(myomics)
    
    if (nrow(SigniReguTargetF) > 0) {  # When there are target features regulated by this regulator
      
      if (is.null(reguValues)) {  # User does not provide reguValues
        reguValues = as.numeric(MLRoutput$arguments$regulatoryData[[myomic]][regulator,])
      }
      
      numtargetFs = length(SigniReguTargetF$targetF)
      if(verbose) {cat(paste(numtargetFs, "target features are regulated by", regulator)); cat("\n")}
      
      if (length(reguValues) > 0) {  # reguValues are available (recovered or given by user)
        
        lapply(1:numtargetFs, function (i) {
          
          targetFValues = MLRoutput$ResultsPerTargetF[[SigniReguTargetF[i,"targetF"]]]$Y$y
          if (order) {
            if(is.null(MLRoutput$arguments$condition)){
              myorder = order(targetFValues)
            } else{
              myorder = order(unlist(MLRoutput$arguments$condition),targetFValues)
              group_names = unique(sort(MLRoutput$arguments$groups))
            }
            targetFValues = targetFValues[myorder]
            reguValues = reguValues[myorder]
            myreplicates = myreplicates[myorder]
            myrepliUni = myrepliUni[myorder]
            if (is.null(cond2plot)) {
              numLines = NULL
            } else {
              cond2plot = cond2plot[myorder]
              condi1 = unique(cond2plot)
              num1 = 1:length(condi1); names(num1) = condi1
              num2 = num1[as.character(cond2plot)]
              if (replicates) {
                num2 = aggregate(num2, by = list("rep" = myreplicates), unique)$x
              }
              numLines = which(diff(num2) != 0)+0.5
            }
          }
          errorValuesRegu = getErrorValues(reguValues, myreplicates)
          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni
          
          errorValues = getErrorValues(targetFValues, myreplicates)
          targetFValues = tapply(targetFValues, myreplicates, mean)
          targetFValues = targetFValues[myrepliUni]
          names(targetFValues) = myrepliUni
          
          
          x.points = 1:length(myrepliUni)
          eje = myrepliUni
          
          p = plotTargetFRegu(x.points = x.points, targetFValues = targetFValues, reguValues = reguValues,
                              col = c(targetF.col, any.col[myomics[i]]), xlab = xlab,
                              yylab = c(SigniReguTargetF[i,"targetF"], regulator), type = c(16,16),
                              main = paste(as.character(SigniReguTargetF[1,c("omic", "area")]), collapse = " "),
                              numLines = numLines, x.names = eje,
                              targetFErrorValues = errorValues, reguErrorValues = errorValuesRegu, group_names = group_names, smooth = smooth,
                              size = size, breakby = breakby)
          print(p)
          
        })
        
      } else { cat("Regulator values could not be recovered from MLRoutput. Please provide them in reguValues argument to generate the plot.\n") }
      
      return(SigniReguTargetF$targetF)
      
    } else { cat(paste("There are no target features relevantly regulated by", regulator)); cat("\n") }
    
  }
  
  
  ## targetF + REGULATOR
  
  if (!is.null(targetF) && !is.null(regulator)) {  ### Plot only the given target feature and the given regulator
    
    targetFResults = MLRoutput$ResultsPerTargetF[[targetF]]
    
    if (is.null(targetFResults)) {
      stop(paste("No MLR was obtained for target feature", targetF))
    } else
    {
      myomic = targetFResults$allRegulators[regulator, "omic"]
      
      if (is.null(reguValues)) {  # User does not provide reguValues
        reguValues = as.numeric(MLRoutput$arguments$regulatoryData[[myomic]][regulator,]) # regulator values
      }
      
      if (length(reguValues) > 0) {  # reguValues are available (recovered or given by user)
        
        targetFValues = MLRoutput$ResultsPerTargetF[[targetF]]$Y$y
        if (order) {
          if(is.null(MLRoutput$arguments$condition)){
            myorder = order(targetFValues)
          } else{
            myorder = order(unlist(MLRoutput$arguments$condition),targetFValues)
            group_names = unique(sort(MLRoutput$arguments$groups))
          }
          targetFValues = targetFValues[myorder]
          reguValues = reguValues[myorder]
          myreplicates = myreplicates[myorder]
          myrepliUni = myrepliUni[myorder]
          if (is.null(cond2plot)) {
            numLines = NULL
          } else {
            cond2plot = cond2plot[myorder]
            condi1 = unique(cond2plot)
            num1 = 1:length(condi1); names(num1) = condi1
            num2 = num1[as.character(cond2plot)]
            if (replicates) {
              num2 = aggregate(num2, by = list("rep" = myreplicates), unique)$x
            }
            numLines = which(diff(num2) != 0)+0.5
          }
        }
        
        errorValuesRegu = getErrorValues(reguValues, myreplicates)
        reguValues = tapply(reguValues, myreplicates, mean, na.rm = TRUE)
        reguValues = reguValues[myrepliUni]
        names(reguValues) = myrepliUni
        
        errorValues = getErrorValues(targetFValues, myreplicates)
        targetFValues = tapply(targetFValues, myreplicates, mean, na.rm = TRUE)
        targetFValues = targetFValues[myrepliUni]
        names(targetFValues) = myrepliUni
        
        x.points = 1:length(myrepliUni)
        eje = myrepliUni
        
        p = plotTargetFRegu(x.points = x.points, targetFValues = targetFValues, reguValues = reguValues,
                            col = c(targetF.col, any.col[myomic]), xlab = xlab,
                            yylab = c(targetF, regulator), type = c(16,16),
                            main = paste(as.character(targetFResults$allRegulators[regulator, c("omic", "area")]),
                                         collapse = " "), numLines = numLines, x.names = eje,
                            targetFErrorValues = errorValues, reguErrorValues = errorValuesRegu,
                            group_names = group_names, smooth = smooth, size = size, breakby = breakby)
        print(p)
        
      } else {
        
        cat("Regulator values could not be recovered from MLRoutput.\n")
        
        regulator = targetFResults$allRegulators[regulator,"filter"]
        
        if (regulator %in% rownames(targetFResults$allRegulators)) {
          
          cat(paste(regulator, "values will be plotted instead.")); cat("\n")
          cat(paste(regulator, "summarizes information from the following correlated regulators:")); cat("\n")
          cat(targetFResults$allRegulators[targetFResults$allRegulators[,"filter"] == regulator,"regulator"]); cat("\n")
          
          reguValues = targetFResults$X
          reguValues = reguValues[, colnames(reguValues) == regulator]
          
          targetFValues = MLRoutput$ResultsPerTargetF[[gene]]$Y$y
          if (order) {
            if(is.null(MLRoutput$arguments$condition)){
              myorder = order(targetFValues)
            } else{
              myorder = order(unlist(MLRoutput$arguments$condition),targetFValues)
              group_names = unique(sort(MLRoutput$arguments$groups))
            }
            targetFValues = targetFValues[myorder]
            reguValues = reguValues[myorder]
            myreplicates = myreplicates[myorder]
            myrepliUni = myrepliUni[myorder]
            if (is.null(cond2plot)) {
              numLines = NULL
            } else {
              cond2plot = cond2plot[myorder]
              condi1 = unique(cond2plot)
              num1 = 1:length(condi1); names(num1) = condi1
              num2 = num1[as.character(cond2plot)]
              if (replicates) {
                num2 = aggregate(num2, by = list("rep" = myreplicates), unique)$x
              }
              numLines = which(diff(num2) != 0)+0.5
            }
          }
          
          errorValuesRegu = getErrorValues(reguValues, myreplicates)
          
          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni
          
          errorValues = getErrorValues(targetFValues, myreplicates)
          targetFValues = tapply(targetFValues, myreplicates, mean)
          targetFValues = targetFValues[myrepliUni]
          names(targetFValues) = myrepliUni
          
          x.points = 1:length(myrepliUni)
          eje = myrepliUni
          
          p = plotTargetFRegu(x.points = x.points, targetFValues = targetFValues, reguValues = reguValues,
                              col = c(targetF.col, any.col[myomic]), xlab = xlab,
                              yylab = c(targetF, regulator), type = c(16,16),
                              main = paste(as.character(targetFResults$allRegulators[regulator, c("omic", "area")]),
                                           collapse = " "), numLines = numLines, x.names = eje,
                              targetFErrorValues = errorValues, reguErrorValues = errorValuesRegu, group_names = group_names, smooth = smooth,
                              size = size, breakby = breakby)
          print(p)
        } else {
          cat("The selected regulator was not declared as relevant by the ElasticNet\n")
          cat("Please, either select another regulator or provide the regulator values.\n")
        }
        
      }
      
    }
    
  }
}

plotPLS = function (PLSoutput, targetF, regulator = NULL, reguValues = NULL, plotPerOmic = FALSE,
                    targetF.col = NULL, regu.col = NULL, order = TRUE,
                    xlab = "", cont.var = NULL, cond2plot = NULL, verbose = TRUE, smooth=TRUE, size = 3, breakby = c(0.2,0.3),...) {
  
  # Colors for omics
  #omic.col = colors()[c(554,89,111,512,17,586,132,428,601,568,86,390,100,200,300,400,500,10,20,30,40,50,60,70,80,90,150,250,350,450,550)]
  
  if (is.null(regu.col)) {
    omic_names = names(PLSoutput$arguments$regulatoryData)
    num_unique = length(omic_names) + 1
    color_palette = colorbiostat(num_unique)
    if(is.null(targetF.col)){
      targetF.col = color_palette[1]
    }
    any.col = color_palette[-1]
  } else {
    if (length(regu.col) == length(PLSoutput$arguments$regulatoryData)) {
      any.col = regu.col
    } else {
      any.col = rep(regu.col, length(PLSoutput$arguments$regulatoryData))
    }
  }
  names(any.col) = names(PLSoutput$arguments$regulatoryData)
  
  
  # Changing margin
  par(mar = c(5,3,4,3)+0.1)
  
  # Groups to plot
  if (is.null(cond2plot)) {
    if (!is.null(PLSoutput$arguments$condition)) {
      cond2plot = apply(PLSoutput$arguments$condition, 1, paste, collapse = "_")
    }
  }
  
  # Replicates
  if (!is.null(cont.var)) {  # we have continuous variable
    if (!is.null(cond2plot)) { # cont.var + cond2plot
      myreplicates = apply(cbind(cond2plot, cont.var), 1, paste, collapse = "_")
      
    } else {  # only continuous variable is provided
      myreplicates = cont.var
    }
    smooth = FALSE
    order = FALSE
  } else {   # no cont.var
    
    if (!is.null(cond2plot)) { # only cond2plot
      myreplicates = colnames(PLSoutput$arguments$targetData)
    } else {  # nothing
      myreplicates = colnames(PLSoutput$arguments$targetData)
    }
  }
  
  # Cast myreplicates to character
  myreplicates = as.character(myreplicates)
  names(myreplicates) = colnames(PLSoutput$arguments$targetData)
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
  group_names = NULL
  
  # Error values
  getErrorValues = function(realValues, repsInfo) {
    # Disable it
    if (! replicates)
      return(NULL)
    
    out_values = tapply(realValues, repsInfo, function(reps) sd(reps)/sqrt(length(reps)))
    
    return(out_values)
  }
  
  
  ## REGULATOR = NULL
  
  if (is.null(regulator)) {  ### Plot all regulators for the given target feature
    
    PLStargetF = PLSoutput$ResultsPerTargetF[[targetF]]
    
    if (is.null(PLStargetF)) {
      stop(paste("No PLS was obtained for target feature", targetF))
    }
    
    if (is.null(PLStargetF$significantRegulators)) { ## No significant regulators
      
      cat("No significant regulators were found for this target feature.\n")
      
      
    } else
    {  ## Significant regulators:
      
      # Considering multicollinearity
      SigReg = PLStargetF$allRegulators
      SigReg = SigReg[SigReg$Sig == 1, c("regulator", "omic", "area", "filter")]
      
      SigReg = SigReg[PLStargetF$significantRegulators,,drop = FALSE]
      
      cat(paste(nrow(SigReg), "significant regulators are to be plotted for target feature", targetF)); cat("\n")
      
      # Target feature values
      targetFValues = PLStargetF$Y$y
      if (order) {
        if(is.null(PLSoutput$arguments$condition)){
          myorder = order(targetFValues)
        } else{
          myorder = order(unlist(PLSoutput$arguments$condition),targetFValues)
          group_names = unique(sort(PLSoutput$arguments$groups))
        }
        cond2plot = cond2plot[myorder]
        targetFValues = targetFValues[myorder]
        myreplicates = myreplicates[myorder]
        myrepliUni = myrepliUni[myorder]
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
      }
      
      errorValues = getErrorValues(targetFValues, myreplicates)
      targetFValues = tapply(targetFValues, myreplicates, mean)
      targetFValues = targetFValues[myrepliUni]
      names(targetFValues) = myrepliUni
      
      # X values
      x.points = 1:length(myrepliUni)
      eje = myrepliUni
      
      if (plotPerOmic) { ## All regulators from the same omic in the same plot
        
        myomics = unique(SigReg$omic)
        
        for (oo in myomics) {
          
          SigRegOmic = SigReg[SigReg$omic == oo,]
          
          omicValues = t(PLSoutput$arguments$regulatoryData[[oo]])
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
          
          p = plotTargetFRegu(x.points = x.points, targetFValues = targetFValues, reguValues = reguValues,
                              col = c(targetF.col, mycol), xlab = xlab, yylab = c(targetF, leftlab), type = c(16,16),
                              main = oo, numLines = numLines, x.names = eje,
                              targetFErrorValues = errorValues, reguErrorValues = errorValuesRegu, group_names = group_names, smooth = smooth,
                              size = size, breakby = breakby)
          
          if (nrow(SigRegOmic) > 1) {
            for (i in 2:nrow(SigRegOmic)) {
              
              reguValues = omicValues[, colnames(omicValues) == SigRegOmic$regulator[i]]
              if (order) reguValues = reguValues[myorder]
              errorValuesRegu = getErrorValues(reguValues, myreplicates)
              
              reguValues = tapply(reguValues, myreplicates, mean)
              reguValues = reguValues[myrepliUni]
              names(reguValues) = myrepliUni
              
              regulator_data = data.frame(x=x.points, y = reguValues)
              
              y_min = min(c(p$data$yleft, regulator_data$y), na.rm = TRUE)
              y_max = max(c(p$data$yleft, regulator_data$y), na.rm = TRUE)
              
              if (! is.null(errorValuesRegu)) {
                
                y_min = min(targetFValues - targetFErrorValues, na.rm = TRUE)
                y_max = max(targetFValues + targetFErrorValues, na.rm = TRUE)
                
                y_minreg = min(reguValues - reguErrorValues, na.rm = TRUE)
                y_maxreg = max(reguValues + reguErrorValues, na.rm = TRUE)
                
                p = p + geom_point(data = regulator_data, aes(x = x, y = y), color = mycol, shape = i, size = 2) +
                  geom_line(data = regulator_data, aes(y = y), color = mycol, linewidth = 1, linetype = i) +
                  
                  # Right axis (Secondary y-axis, rescaled to match yright values)
                  geom_point(aes(y = (yright - min(yright)) / (max(yright) - min(yright)) * (max(yleft) - min(yleft)) + min(yleft)), color = col[1], shape = type[2], size = 2) +
                  geom_line(aes(y = (yright - min(yright)) / (max(yright) - min(yright)) * (max(yleft) - min(yleft)) + min(yleft)), color = col[1], linewidth = lwd[2]) +
                  
                  # Adjust the axis to include the error value
                  geom_errorbar(aes(ymin = yleft - yleftErrorValues, ymax = yleft + yleftErrorValues), 
                                width = 0.2, color = col[2]) +
                  geom_errorbar(aes(
                    ymin = (yright - yrightErrorValues - min(yright, na.rm = TRUE)) / (max(yright,na.rm = TRUE) - min(yright,na.rm = TRUE)) * (max(yleft,na.rm = TRUE) - min(yleft,na.rm = TRUE)) + min(yleft,na.rm = TRUE),
                    ymax = (yright + yrightErrorValues - min(yright, na.rm = TRUE)) / (max(yright,na.rm = TRUE) - min(yright,na.rm = TRUE)) * (max(yleft,na.rm = TRUE) - min(yleft, na.rm = TRUE)) + min(yleft,na.rm = TRUE)
                  ), width = 0.2, color = col[1]) +
                  
                  scale_y_continuous(
                    limits = c(y_minreg, y_maxreg),   # Expand the limits
                    breaks = seq(y_minreg, y_maxreg, by = breakby[1]),
                    labels = function(x) sprintf("%.1f", x),
                    sec.axis = sec_axis(
                      transform = ~ . * (y_max - y_min) / (y_maxreg - y_minreg) + y_min - y_minreg * (y_max - y_min) / (y_maxreg - y_minreg),
                      breaks = seq(y_min, y_max, by = breakby[2]),
                      labels = function(x) sprintf("%.1f", x), 
                      name = targetF
                    )
                  )
                
              } else{
                p = p + 
                  scale_y_continuous(
                    limits = c(y_min, y_max),   # Expand the limits
                    breaks = seq(y_min, y_max, by = breakby[1]),
                    labels = function(x) sprintf("%.1f", x),
                    sec.axis = sec_axis(
                      transform = ~ . * (max(p$data$yright) - min(p$data$yright)) / (max(p$data$yleft) - min(p$data$yleft)) + 
                        min(p$data$yright) - min(p$data$yleft) * (max(p$data$yright) - min(p$data$yright)) / (max(p$data$yleft) - min(p$data$yleft)),
                      breaks = seq(min(p$data$yright), max(p$data$yright), by = breakby[2]),
                      labels = function(x) sprintf("%.1f", x), 
                      name = targetF
                    )
                  ) + geom_point(data = regulator_data, aes(x = x, y = y), color = mycol, shape = i, size = 2) +
                  geom_line(data = regulator_data, aes(y = y), color = mycol, linewidth = 1, linetype = i) 
              }
            }
          }
          print(p)
        }
        
      } else {  ## Each regulator in a separate plot
        
        for (rr in SigReg$regulator) {
          
          oo = PLSoutput$ResultsPerTargetF[[targetF]]$allRegulators[rr,"omic"]
          omicValues = t(PLSoutput$arguments$regulatoryData[[oo]])
          reguValues = omicValues[, colnames(omicValues) == rr]
          if (order) reguValues = reguValues[myorder]
          errorValuesRegu = getErrorValues(reguValues, myreplicates)
          
          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni
          
          mycol = any.col[oo]
          
          p = plotTargetFRegu(x.points = x.points, targetFValues = targetFValues, reguValues = reguValues,
                              col = c(targetF.col, mycol), xlab = xlab,
                              yylab = c(targetF, rr), type = c(16,16),
                              main = paste(as.character(SigReg[rr, c("omic", "area")]), collapse = " "),
                              numLines = numLines, x.names = eje,
                              targetFErrorValues = errorValues, reguErrorValues = errorValuesRegu, group_names = group_names, smooth = smooth,
                              size = size, breakby = breakby)
          print(p)
        }
        
      }
      
      return(PLStargetF$allRegulators[PLStargetF$significantRegulators, -6])
    }
    
  }
  
  ## targetF = NULL
  
  if (is.null(targetF)) {  ### Plot all target features regulated by the regulator
    
    SigniReguTargetF = GetPairstargetFRegulator(targetFs = NULL, output = PLSoutput)
    SigniReguTargetF = SigniReguTargetF[SigniReguTargetF[,"regulator"] == regulator,]
    myomics = SigniReguTargetF[,"omic"]
    myomic = unique(myomics)
    
    if (nrow(SigniReguTargetF) > 0) {  # When there are targetFs regulated by this regulator
      
      if (is.null(reguValues)) {  # User does not provide reguValues
        reguValues = as.numeric(PLSoutput$arguments$regulatoryData[[myomic]][regulator,])
      }
      
      numtargetFs = length(SigniReguTargetF$targetF)
      if(verbose) {cat(paste(numtargetFs, "target features are regulated by", regulator)); cat("\n")}
      
      if (length(reguValues) > 0) {  # reguValues are available (recovered or given by user)
        
        lapply(1:numtargetFs, function (i) {
          
          targetFValues = PLSoutput$ResultsPerTargetF[[SigniReguTargetF[i,"targetF"]]]$Y$y
          if (order) {
            if(is.null(PLSoutput$arguments$condition)){
              myorder = order(targetFValues)
            } else{
              myorder = order(unlist(PLSoutput$arguments$condition),targetFValues)
              group_names = unique(sort(PLSoutput$arguments$groups))
            }
            targetFValues = targetFValues[myorder]
            reguValues = reguValues[myorder]
            myreplicates = myreplicates[myorder]
            myrepliUni = myrepliUni[myorder]
            if (is.null(cond2plot)) {
              numLines = NULL
            } else {
              cond2plot = cond2plot[myorder]
              condi1 = unique(cond2plot)
              num1 = 1:length(condi1); names(num1) = condi1
              num2 = num1[as.character(cond2plot)]
              if (replicates) {
                num2 = aggregate(num2, by = list("rep" = myreplicates), unique)$x
              }
              numLines = which(diff(num2) != 0)+0.5
            }
          }
          errorValuesRegu = getErrorValues(reguValues, myreplicates)
          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni
          
          errorValues = getErrorValues(targetFValues, myreplicates)
          targetFValues = tapply(targetFValues, myreplicates, mean)
          targetFValues = targetFValues[myrepliUni]
          names(targetFValues) = myrepliUni
          
          x.points = 1:length(myrepliUni)
          eje = myrepliUni
          
          p = plotTargetFRegu(x.points = x.points, targetFValues = targetFValues, reguValues = reguValues,
                              col = c(targetF.col, any.col[myomics[i]]), xlab = xlab,
                              yylab = c(SigniReguTargetF[i,"targetF"], regulator), type = c(16,16),
                              main = paste(as.character(SigniReguTargetF[1,c("omic", "area")]), collapse = " "),
                              numLines = numLines, x.names = eje,
                              targetFErrorValues = errorValues, reguErrorValues = errorValuesRegu, group_names = group_names, smooth = smooth,
                              size = size, breakby = breakby)
          print(p)
        })
      } else { cat("Regulator values could not be recovered from output. Please provide them in reguValues argument to generate the plot.\n") }
      
      return(SigniReguTargetF$targetF)
      
    } else { cat(paste("There are no target features significantly regulated by", regulator)); cat("\n") }
    
  }
  
  
  ## targetF + REGULATOR
  
  if (!is.null(targetF) && !is.null(regulator)) {  ### Plot only the given target feature and the given regulator
    
    targetFResults = PLSoutput$ResultsPerTargetF[[targetF]]
    
    if (is.null(targetFResults)) {
      stop(paste("No PLS was obtained for target feature", targetF))
    } else
    {
      myomic = targetFResults$allRegulators[regulator, "omic"]
      
      if (is.null(reguValues)) {  # User does not provide reguValues
        reguValues = as.numeric(PLSoutput$arguments$regulatoryData[[myomic]][regulator,]) # regulator values
      }
      
      if (length(reguValues) > 0) {  # reguValues are available (recovered or given by user)
        
        targetFValues = PLSoutput$ResultsPerTargetF[[targetF]]$Y$y
        if (order) {
          if(is.null(PLSoutput$arguments$condition)){
            myorder = order(targetFValues)
          } else{
            myorder = order(unlist(PLSoutput$arguments$condition),targetFValues)
            group_names = unique(sort(PLSoutput$arguments$groups))
          }
          targetFValues = targetFValues[myorder]
          reguValues = reguValues[myorder]
          myreplicates = myreplicates[myorder]
          myrepliUni = myrepliUni[myorder]
          if (is.null(cond2plot)) {
            numLines = NULL
          } else {
            cond2plot = cond2plot[myorder]
            condi1 = unique(cond2plot)
            num1 = 1:length(condi1); names(num1) = condi1
            num2 = num1[as.character(cond2plot)]
            if (replicates) {
              num2 = aggregate(num2, by = list("rep" = myreplicates), unique)$x
            }
            numLines = which(diff(num2) != 0)+0.5
          }
        }
        
        errorValuesRegu = getErrorValues(reguValues, myreplicates)
        reguValues = tapply(reguValues, myreplicates, mean, na.rm = TRUE)
        reguValues = reguValues[myrepliUni]
        names(reguValues) = myrepliUni
        
        errorValues = getErrorValues(targetFValues, myreplicates)
        targetFValues = tapply(targetFValues, myreplicates, mean)
        targetFValues = targetFValues[myrepliUni]
        names(targetFValues) = myrepliUni
        
        x.points = 1:length(myrepliUni)
        eje = myrepliUni
        
        p = plotTargetFRegu(x.points = x.points, targetFValues = targetFValues, reguValues = reguValues,
                            col = c(targetF.col, any.col[myomic]), xlab = xlab,
                            yylab = c(targetF, regulator), type = c(16,16),
                            main = paste(as.character(targetFResults$allRegulators[regulator, c("omic", "area")]),
                                         collapse = " "),
                            numLines = numLines, x.names = eje,
                            targetFErrorValues = errorValues, reguErrorValues = errorValuesRegu, group_names = group_names, smooth = smooth,
                            size = size, breakby = breakby)
        print(p)
        
      } else {
        
        cat("Regulator values could not be recovered from output.\n")
        
        regulator = targetFResults$allRegulators[regulator,"filter"]
        
        if (regulator %in% rownames(targetFResults$allRegulators)) {
          
          cat(paste(regulator, "values will be plotted instead.")); cat("\n")
          cat(paste(regulator, "summarizes information from the following correlated regulators:")); cat("\n")
          cat(targetFResults$allRegulators[targetFResults$allRegulators[,"filter"] == regulator,"regulator"]); cat("\n")
          
          reguValues = targetFResults$X
          reguValues = reguValues[, colnames(reguValues) == regulator]
          targetFValues = PLSoutput$ResultsPerTargetF[[targetF]]$Y$y
          if (order) {
            if(is.null(PLSoutput$arguments$condition)){
              myorder = order(targetFValues)
            } else{
              myorder = order(unlist(PLSoutput$arguments$condition),targetFValues)
              group_names = unique(sort(PLSoutput$arguments$groups))
            }
            targetFValues = targetFValues[myorder]
            reguValues = reguValues[myorder]
            myreplicates = myreplicates[myorder]
            myrepliUni = myrepliUni[myorder]
            
            if (is.null(cond2plot)) {
              numLines = NULL
            } else {
              cond2plot = cond2plot[myorder]
              condi1 = unique(cond2plot)
              num1 = 1:length(condi1); names(num1) = condi1
              num2 = num1[as.character(cond2plot)]
              if (replicates) {
                num2 = aggregate(num2, by = list("rep" = myreplicates), unique)$x
              }
              numLines = which(diff(num2) != 0)+0.5
            }
          }
          
          errorValuesRegu = getErrorValues(reguValues, myreplicates)
          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni
          
          errorValues = getErrorValues(targetFValues, myreplicates)
          targetFValues = tapply(targetFValues, myreplicates, mean)
          targetFValues = targetFValues[myrepliUni]
          names(targetFValues) = myrepliUni
          
          x.points = 1:length(myrepliUni)
          eje = myrepliUni
          
          p = plotTargetFRegu(x.points = x.points, targetFValues = targetFValues, reguValues = reguValues,
                              col = c(targetF.col, any.col[myomic]), xlab = xlab,
                              yylab = c(targetF, regulator), type = c(16,16),
                              main = paste(as.character(targetFResults$allRegulators[regulator, c("omic", "area")]),
                                           collapse = " "),
                              numLines = numLines, x.names = eje,
                              targetFErrorValues = errorValues, reguErrorValues = errorValuesRegu, group_names = group_names, smooth = smooth,
                              size = size, breakby = breakby)
          print(p)
        } else {
          cat("The selected regulator was not declared as significant by the PLS model\n")
          cat("Please, either select another regulator or provide the regulator values.\n")
        }
        
      }
      
    }
    
  }
  
}


## Plots PLS ----------

#' plotWeight
#'
#' \code{plotWeight} Function to be applied to \link{more} main function output.
#' 
#' @param output Output object of MORE main function.
#' @param targetF Target feature of which the user wants visualize the weightings of the model.
#' @param axe1 Number of the PLS component to plot in the X axis
#' @param axe2 Number of the PLS component to plot in the Y axis
#' 
#' @return Plots the weighting star of the regulators identified as significant for the selected target feature.
#' @export

plotWeight<-function(output, targetF,axe1=1,axe2=2){
  
  if(output$GlobalSummary$GoodnessOfFit[targetF,'ncomp']==1) {
    warning('The original PLS model extracted a unique component. The visuallization would be hard')
    if (ncol(output$arguments$targetData)<7){cross = output$arguments$targetData -2}else{cross =7}
    pls = ropls::opls(output$ResultsPerTargetF[[targetF]]$X[,output$ResultsPerTargetF[[targetF]]$significantRegulators,drop=FALSE], output$ResultsPerTargetF[[targetF]]$Y$y, info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0, predI=1)
    
    plot(pls@weightStarMN[,1], rep(0,length(output$ResultsPerTargetF[[targetF]]$significantRegulators)),
         main = "Weights*",
         xlab = paste0('w*c', axe1), ylab = paste0('w*c', axe2),
         pch = 18, col = "blue", xlim = c(-1,1))
    
    points(pls@cMN[,1], 0, pch = 18, col = "red")
    # Asignamos las etiquetas
    text(pls@weightStarMN[,1], rep(0,length(output$ResultsPerTargetF[[targetF]]$significantRegulators)),
         labels = row.names(pls@weightStarMN),
         cex = 0.6, srt = 60, pos = 2, col = "black")
    text(pls@cMN[,1], 0, labels = targetF, cex =
           0.6, srt = 60, pos = 4, col = 'black')
    abline(h=0)
    legend("topleft", legend = c("X", targetF),cex = 0.5,
           pch = 18, col = c("blue", "red"))
  } else{
    if(axe1> output$GlobalSummary$GoodnessOfFit[targetF,'ncomp'] | axe2>output$GlobalSummary$GoodnessOfFit[targetF,'ncomp']) stop('Error: The model did not extracted originally that many components')
    if (ncol(output$arguments$targetData)<7){cross = output$arguments$targetData -2}else{cross =7}
    #Create the PLS model only with the variables that resulted significant in the model
    pls = ropls::opls(output$ResultsPerTargetF[[targetF]]$X[,output$ResultsPerTargetF[[targetF]]$significantRegulators,drop=FALSE], output$ResultsPerTargetF[[targetF]]$Y$y, info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0, predI=output$GlobalSummary$GoodnessOfFit[targetF,'ncomp'])

    #Create the weighting plots
    plot(pls@weightStarMN[,axe1], pls@weightStarMN[,axe2],
         main = "Weights*",
         xlab = paste0('w*c', axe1), ylab = paste0('w*c', axe2),
         pch = 18, col = "blue", asp=1)
    
    points(pls@cMN[,axe1], pls@cMN[,axe2], pch = 18, col = "red")
    # Asignamos las etiquetas
    text(pls@weightStarMN[,axe1], pls@weightStarMN[,axe2],
         labels = row.names(pls@weightStarMN),
         cex = 0.6, pos = 4, col = "black")
    text(pls@cMN[,axe1], pls@cMN[,axe2], labels = targetF, cex =
           0.6, pos = 4, col = 'black')
    abline(h=0, v=0)
    legend("topleft", legend = c("X", targetF),cex = 0.5,
           pch = 18, col = c("blue", "red"))
  }
  
}

#' plotScores
#'
#' \code{plotScores} Function to be applied to \link{more} main function output.
#' 
#' @param output Output object of MORE main function.
#' @param targetF Target feature of which the user wants visualize the weightings of the model.
#' @param axe1 Number of the PLS component to plot in the X axis
#' @param axe2 Number of the PLS component to plot in the Y axis
#' 
#' @return Plots the scores of the samples under the PLS model generated by the significant regulators for the selected target feature.
#' @export

plotScores<-function(output, targetF,axe1=1,axe2=2){
  
  if(output$GlobalSummary$GoodnessOfFit[targetF,'ncomp']==1) {
    warning('The original PLS model extracted a unique component. The visuallization would be hard')
    if (ncol(output$arguments$targetData)<7){cross = output$arguments$targetData -2}else{cross =7}
    pls = ropls::opls(output$ResultsPerTargetF[[targetF]]$X[,output$ResultsPerTargetF[[targetF]]$significantRegulators,drop=FALSE], output$ResultsPerTargetF[[targetF]]$Y$y, info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0, predI=output$GlobalSummary$GoodnessOfFit[targetF,'ncomp'])
    
    num_unique <- length(unique(output$arguments$groups)) + 1
    color_palette <- colorbiostat(num_unique)
    custom_colors <- setNames(color_palette[-1], unique(output$arguments$groups))
    
    plot(pls@scoreMN[,1], rep(0,length(output$arguments$groups)),
         main = "Scores",
         xlab = paste0('t', axe1), ylab = paste0('t', axe2),
         pch = 18, col = custom_colors)
    
    # Asignamos las etiquetas
    text(pls@scoreMN[,1], rep(0,length(output$arguments$groups)),
         labels = row.names(pls@scoreMN),
         cex = 0.6, srt = 60, pos = 2, col = 'black')
    abline(h=0)
    legend("topleft", legend = names(custom_colors),cex = 0.5,
           pch = 18, col = custom_colors)
  } else{
    if(axe1> output$GlobalSummary$GoodnessOfFit[targetF,'ncomp'] | axe2>output$GlobalSummary$GoodnessOfFit[targetF,'ncomp']) stop('Error: The model did not extracted originally that many components')
    if (ncol(output$arguments$targetData)<7){cross = output$arguments$targetData -2}else{cross =7}
    #Create the PLS model only with the variables that resulted significant in the model
    pls = ropls::opls(output$ResultsPerTargetF[[targetF]]$X[,output$ResultsPerTargetF[[targetF]]$significantRegulators,drop=FALSE], output$ResultsPerTargetF[[targetF]]$Y$y, info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0, predI=output$GlobalSummary$GoodnessOfFit[targetF,'ncomp'])
    
    num_unique <- length(unique(output$arguments$groups)) + 1
    color_palette <- colorbiostat(num_unique)
    custom_colors <- setNames(color_palette[-1], unique(output$arguments$groups))
    
    #Create the weighting plots
    plot(pls@scoreMN[,axe1], pls@scoreMN[,axe2],
         main = "Scores",
         xlab = paste0('t', axe1), ylab = paste0('t', axe2),
         pch = 18, col = custom_colors)
    
    # Asignamos las etiquetas
    text(pls@scoreMN[,axe1], pls@scoreMN[,axe2],
         labels = row.names(pls@scoreMN),
         cex = 0.6, pos = 4, col = "black")
    abline(h=0, v=0)
    legend("topleft", legend = names(custom_colors),cex = 0.5,
           pch = 18, col = custom_colors)
  }
  
}


## Summary ------------


#' summary
#'
#' \code{summary.MORE} Function to be applied to MORE object.
#' 
#' @param object MORE object obtained from applying \link{more} function. 
#' @param plot.more If TRUE top 10 global regulators are plotted against the target features they regulate. By default, FALSE.
#' 
#' @return Summary of more analysis.
#' @export

summary.MORE <-function(object, plot.more=FALSE){
  
  cat('A model was computed for',length(object$ResultsPerTargetF), 'target features.' ,'\n')
  cat(ifelse(is.null(object$GlobalSummary$TargetFNOregu),0,nrow(object$GlobalSummary$TargetFNOregu)), 'target features had no intial regulators.' ,'\n')
  
  if(object$arguments$method == 'MLR'||object$arguments$method=='ISGL'){
    cat('For', ifelse(is.null(object$GlobalSummary$TargetFNOmodel),0,nrow(object$GlobalSummary$TargetFNOmodel)), 'target features, the final MLR model could not be obtained.','\n')
    cat('Target features presented a mean of ',mean(na.omit(object$GlobalSummary$GoodnessOfFit[,'relReg'])),'relevant regulators.','\n')
    
    #Top hub target features
    relevant_regulators<-object$GlobalSummary$ReguPerTargetF[,c(grep('-Rel$',colnames(object$GlobalSummary$ReguPerTargetF))),drop=FALSE]
    #globally
    
    s_rel_reg<-apply(relevant_regulators, 1, sum)
    
    cat('These are the top 10 hub target features and the number of relevant regulators for each:\n')
    print(s_rel_reg[rev(tail(order(s_rel_reg),10))])
    
    #Global regulators
    
    m_rel_reg<-lapply(object$ResultsPerTargetF, function(x) x$relevantRegulators)
    m_rel_reg <- unlist(m_rel_reg)
    
    ## Count occurrences
    mrel_vector <- table(m_rel_reg)
    #Ask to regulate at least 10 target features
    mrel_vector<-mrel_vector[mrel_vector>10]
    if(length(mrel_vector!=0)){
      cat('These are the top 10 global regulators and the number of target features that they regulate:\n')
      print(mrel_vector[rev(tail(order(mrel_vector),10))])
    } else{
      cat('There were not global regulators (regulators that regulate more than 10 target features).')
    }
    
    
    if(plot.more){
      mreg<-mrel_vector[rev(tail(order(mrel_vector),10))]
      for (i in 1:10) {
        par(mfrow=c(2,4))
        plotMLR(object, targetF = NULL, regulator = names(mreg)[i], plotPerOmic = FALSE ,order = FALSE, targetF.col = 'skyblue', regu.col = 'tan1', verbose = FALSE)
      }
    }
  }
  else{
    cat('Target features presented a mean of ',mean(na.omit(object$GlobalSummary$GoodnessOfFit[,'sigReg'])),'significant regulators.','\n')
    
    #Top hub target features
    significant_regulators<-object$GlobalSummary$ReguPerTargetF[,c(grep('-Sig$',colnames(object$GlobalSummary$ReguPerTargetF))),drop=FALSE]
    #globally
    
    s_sig_reg<-apply(significant_regulators, 1, sum)
    cat('These are the top 10 hub target features and the number of significant regulators for each:\n')
    print(s_sig_reg[tail(order(s_sig_reg),10)])
    
    #Global regulators
    
    m_sig_reg<-lapply(object$ResultsPerTargetF, function(x) x$significantRegulators)
    m_sig_reg <- unlist(m_sig_reg)
    
    ## Count occurrences
    msig_vector <- table(m_sig_reg)
    #Ask to regulate at least 10 target features
    msig_vector<-msig_vector[msig_vector>10]
    cat('These are the top 10 global regulators and the number of target features that they regulate:\n')
    print(msig_vector[rev(tail(order(msig_vector),10))])
    
    if(plot.more){
      msig<-msig_vector[rev(tail(order(msig_vector),10))]
      for (i in 1:10) {
        par(mfrow=c(2,4))
        plotPLS(object, targetF = NULL, regulator = names(msig)[i], plotPerOmic = FALSE ,order = FALSE, targetF.col = 'skyblue', regu.col = 'tan1', verbose = FALSE)
        
      }
    }
    
  }
  
}

#' summaryPlot
#'
#' \code{summaryPlot} Function to be applied to MORE object and the object obtained from RegulationPerCondition.
#' 
#' @param output MORE object obtained from applying \link{more} function. 
#' @param outputRegpcond Object generated by the function \link{RegulationPerCondition} when applied to a \link{more} object. 
#' @param filterR2 Highlights the results for the genes that showed a R2 above the one indicated. By default, 0.
#' @param byTargetF If TRUE, the function plots the percentage of target features with significant regulators globally and per omic. If FALSE, it plots the percentage of significant regulations per omic. By default, TRUE.
#' 
#' @return Summary plot of the MORE analysis.
#' @export 


summaryPlot <- function(output, outputRegpcond, filterR2 = 0, byTargetF = TRUE) {
  
  
  if (byTargetF) {
    # Calculate percentages and prepare data for plotting
    ngroups <- length(unique(output$arguments$groups))
    omics <- names(output$arguments$regulatoryData)
    totaltargetFs <- length(output$ResultsPerTargetF) +
      ifelse(is.null(output$GlobalSummary$TargetFNOregu), 0, length(output$GlobalSummary$TargetFNOregu)) +
      ifelse(is.null(output$GlobalSummary$TargetFNOmodel), 0, length(output$GlobalSummary$TargetFNOmodel))
    
    # If R2filter is applied calculate the ones that are specific for the genes that fulfill the threshold
    filtered_targetF <- rownames(output$GlobalSummary$GoodnessOfFit)[which(output$GlobalSummary$GoodnessOfFit[, 1] > filterR2)]
    filtered_outputRegpcond <- outputRegpcond[outputRegpcond$targetF %in% filtered_targetF,,drop=FALSE]
    
    pos <- grep("Group", colnames(outputRegpcond))[1]
    if(is.na(pos)){
      
      cts <- matrix(NA, nrow = 1, ncol = length(omics) + 1)
      cts2 <- matrix(NA, nrow = 1, ncol = length(omics) + 1)
      
      for (j in 1:((length(omics)) + 1)) {
        if ( j == 1) {
          cts[1, j] <- length(unique(outputRegpcond$targetF))
          cts2[1, j] <- length(unique(filtered_outputRegpcond$targetF))
        } else {
          cts[1, j] <- length(unique(outputRegpcond[outputRegpcond$omic == omics[j - 1], ]$targetF))
          cts2[1, j] <- length(unique(filtered_outputRegpcond[filtered_outputRegpcond$omic == omics[j - 1], ]$targetF))
        } 
      }
      
      cts <- cts / totaltargetFs * 100
      cts2 <- cts2 / totaltargetFs * 100
      
      df <- data.frame(
        omic = c("Any", names(output$arguments$regulatoryData)),
        targetFs = as.vector(cts),
        filteredR2 =  as.vector(cts2) )
      
      num_unique <- length(omics)+1
      color_palette <- colorbiostat(num_unique)
      custom_colors <- setNames(color_palette, c('Any',omics))
      
      ggplot2::ggplot() +
        # First layer: filteredR2
        geom_bar(data = df, aes(x = omic, y = filteredR2, fill = omic, alpha = "R2"),
                 stat = "identity", position = position_dodge()) +
        # Second layer: targetFs
        geom_bar(data = df, aes(x = omic, y = targetFs, fill = omic, alpha = "None"),
                 stat = "identity", position = position_dodge()) +
        theme_minimal() +
        scale_fill_manual(values = custom_colors) +
        scale_alpha_manual(values = c("R2" = 1, "None" = 0.3), labels = c("R2" = paste0("R2>", filterR2), "None" = "None"), guide = guide_legend(title = "Model filtering")) +
        theme(legend.text = element_text(size = 12), panel.grid = element_line(color = "black", linewidth = 0.5, linetype = 1),
              axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11)) +
        labs( x = "", y = "% targetFs with significant regulators")
      
    } else{
      cts <- matrix(NA, nrow = (ngroups) + 1, ncol = length(omics) + 1)
      cts2 <- matrix(NA, nrow = (ngroups) + 1, ncol = length(omics) + 1)
      for (i in 1:((ngroups) + 1)) {
        for (j in 1:((length(omics)) + 1)) {
          if (i == 1 && j == 1) {
            cts[i, j] <- length(unique(outputRegpcond$targetF))
            cts2[i, j] <- length(unique(filtered_outputRegpcond$targetF))
          } else if (i == 1) {
            cts[i, j] <- length(unique(outputRegpcond[outputRegpcond$omic == omics[j - 1], ]$targetF))
            cts2[i, j] <- length(unique(filtered_outputRegpcond[filtered_outputRegpcond$omic == omics[j - 1], ]$targetF))
          } else if (j == 1) {
            cts[i, j] <- length(unique(outputRegpcond[outputRegpcond[, pos + (i - 2)] != 0, ]$targetF))
            cts2[i, j] <- length(unique(filtered_outputRegpcond[filtered_outputRegpcond[, pos + (i - 2)] != 0, ]$targetF))
          } else {
            sub_gr <- outputRegpcond[outputRegpcond$omic == omics[j - 1], ]
            cts[i, j] <- length(unique(sub_gr[sub_gr[, pos + (i - 2)] != 0, ]$targetF))
            sub_gr <- filtered_outputRegpcond[filtered_outputRegpcond$omic == omics[j - 1], ]
            cts2[i, j] <- length(unique(sub_gr[sub_gr[, pos + (i - 2)] != 0, ]$targetF))
          }
        }
      }
      
      cts <- cts / totaltargetFs * 100
      cts2 <- cts2 / totaltargetFs * 100
      group_levels <- c("Global", gsub("Group_", "", colnames(outputRegpcond)[pos:ncol(outputRegpcond)]))
      
      df <- data.frame(
        Group = factor(rep(group_levels, times = length(omics) + 1), levels = group_levels),
        omic = rep(c("Any", names(output$arguments$regulatoryData)), each = ngroups + 1),
        targetFs = as.vector(cts),
        filteredR2 =  as.vector(cts2) )
      
      num_unique <- ngroups + 1
      color_palette <- colorbiostat(num_unique)
      custom_colors <- setNames(color_palette, unique(df$Group))
      
      ggplot2::ggplot() +
        # First layer: filteredR2
        geom_bar(data = df, aes(x = omic, y = filteredR2, fill = Group, alpha = "R2"),
                 stat = "identity", position = position_dodge()) +
        # Second layer: targetFs
        geom_bar(data = df, aes(x = omic, y = targetFs, fill = Group, alpha = "None"),
                 stat = "identity", position = position_dodge()) +
        theme_minimal() +
        scale_fill_manual(values = custom_colors) +
        scale_alpha_manual(values = c("R2" = 1, "None" = 0.3), labels = c("R2" = paste0("R2>", filterR2), "None" = "None"), guide = guide_legend(title = "Model filtering")) +
        theme(legend.text = element_text(size = 12), panel.grid = element_line(color = "black", linewidth = 0.5, linetype = 1),
              axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11)) +
        labs( x = "", y = "% targetFs with significant regulators") 
    }
    
    
  } else {
    
    ngroups = length(unique(output$arguments$groups))
    omics = names(output$arguments$regulatoryData)
    
    total_reg_omic <- if (is.null(output$arguments$associations)) {
      sapply(output$arguments$regulatoryData, nrow)
    } else {
      sapply(omics, function(x) {
        if (!is.null(output$arguments$associations[[x]])) {
          temp =output$arguments$associations[[x]][output$arguments$associations[[x]][,1] %in% names(output$ResultsPerTargetF), ]
          nrow(temp[temp[,2] %in% rownames(output$arguments$regulatoryData[[x]]), ])
        } else {
          nrow(output$arguments$regulatoryData[[x]])
        }
      })
    }
    
    # If R2filter is applied calculate the ones that are specific for the genes that fulfill the threshold
    filtered_targetF <- rownames(output$GlobalSummary$GoodnessOfFit)[which(output$GlobalSummary$GoodnessOfFit[, 1] > filterR2)]
    filtered_outputRegpcond <- outputRegpcond[outputRegpcond$targetF %in% filtered_targetF,,drop=FALSE]
    
    
    pos = grep('Group',colnames(outputRegpcond))[1]
    #Create all the counts needed globally and per groups
    
    if(is.na(pos)){
      
      cts = matrix(NA, nrow=1,ncol=length(omics))
      cts2 <- matrix(NA, nrow = 1, ncol = length(omics))
      
      for (j in 1:length(omics)){
        cts[1,j] = nrow(outputRegpcond[outputRegpcond$regulator %in% outputRegpcond[outputRegpcond$omic==omics[j],]$regulator,])/total_reg_omic[j]*100
        cts2[1,j] = nrow(filtered_outputRegpcond[filtered_outputRegpcond$regulator %in% filtered_outputRegpcond[filtered_outputRegpcond$omic==omics[j],]$regulator,])/total_reg_omic[j]*100
      }
      
      #Create a df with the percentage of target features with significant regulators by omic and condition
      df <- data.frame(omic=names(output$arguments$regulatoryData),
                       targetFs=as.vector(cts), filteredR2 =  as.vector(cts2))
      
      num_unique <- length(omics)+1
      color_palette <- colorbiostat(num_unique)
      custom_colors <- setNames(color_palette[-1], omics)
      
      ggplot2::ggplot() +
        # First layer: filteredR2
        geom_bar(data = df, aes(x = omic, y = filteredR2, fill = omic, alpha = "R2"),
                 stat = "identity", position = position_dodge()) +
        # Second layer: targetFs
        geom_bar(data = df, aes(x = omic, y = targetFs, fill = omic, alpha = "None"),
                 stat = "identity", position = position_dodge()) +
        theme_minimal() + scale_x_discrete(labels = paste(unique(df$omic),'\n',total_reg_omic)) +
        scale_fill_manual(values = custom_colors) +
        scale_alpha_manual(values = c("R2" = 1, "None" = 0.3), labels = c("R2" = paste0("R2>", filterR2), "None" = "None"),  guide = guide_legend(title = "Model filtering")) +
        scale_y_continuous(limits = c(0, max(df$targetFs) + 1)) +
        theme(legend.text = element_text(size = 12), panel.grid = element_line(color = "black", linewidth = 0.5, linetype = 1),
              axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11)) +
        labs( x = "Number of associations per omic", y = "% significant regulations") 
      
    } else{
      cts = matrix(NA, nrow=(ngroups),ncol=length(omics))
      cts2 <- matrix(NA, nrow = (ngroups), ncol = length(omics))
      
      for (i in 1:ngroups){
        #Create the global values
        for (j in 1:length(omics)){
          temp =outputRegpcond[outputRegpcond[,pos+i-1]!=0,]
          cts[i,j] = nrow(temp[temp$regulator %in% outputRegpcond[outputRegpcond$omic==omics[j],]$regulator,])/total_reg_omic[j]*100
          temp =filtered_outputRegpcond[filtered_outputRegpcond[,pos+i-1]!=0,]
          cts2[i,j] = nrow(temp[temp$regulator %in% filtered_outputRegpcond[filtered_outputRegpcond$omic==omics[j],]$regulator,])/total_reg_omic[j]*100
          
        }
      }
      group_levels <- gsub('Group_','',colnames(outputRegpcond)[pos:ncol(outputRegpcond)])
      #Create a df with the percentage of target features with significant regulators by omic and condition
      df <- data.frame(Group=factor(rep(group_levels, times=length(omics)), levels = group_levels),
                       omic=rep(names(output$arguments$regulatoryData),each = ngroups),
                       targetFs=as.vector(cts), filteredR2 =  as.vector(cts2))
      
      num_unique <- ngroups+1
      color_palette <- colorbiostat(num_unique)
      custom_colors <- setNames(color_palette[-1], unique(df$Group))
      
      ggplot2::ggplot() +
        # First layer: filteredR2
        geom_bar(data = df, aes(x = omic, y = filteredR2, fill = Group, alpha = "R2"),
                 stat = "identity", position = position_dodge()) +
        # Second layer: targetFs
        geom_bar(data = df, aes(x = omic, y = targetFs, fill = Group, alpha = "None"),
                 stat = "identity", position = position_dodge()) +
        theme_minimal() + scale_x_discrete(labels = paste(unique(df$omic),'\n',total_reg_omic)) +
        scale_fill_manual(values = custom_colors) +
        scale_alpha_manual(values = c("R2" = 1, "None" = 0.3), labels = c("R2" = paste0("R2>", filterR2), "None" = "None"),  guide = guide_legend(title = "Model filtering")) +
        scale_y_continuous(limits = c(0, max(df$targetFs) + 1)) +
        theme(legend.text = element_text(size = 12), panel.grid = element_line(color = "black", linewidth = 0.5, linetype = 1),
              axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11)) +
        labs( x = "Number of associations per omic", y = "% significant regulations") 
      
    }
    
    
  }
}

getallreg <- function(x, targetF) {
  reg <- x[x[, 1] == targetF, 2]
  return(as.character(reg))
}

#' globalregPlot
#'
#' \code{globalregPlot} Function to be applied to \link{RegulationInCondition} function output.
#' 
#' @param outputRegincond Output object of running \link{RegulationInCondition} function.
#' @param byNetwork If TRUE, information would be plotted on a network instead of a corrplot. By default, FALSE.
#' 
#' @return Graphical visualization between global regulators and target features in a specific condition.
#' @export

globalregPlot<-function(outputRegincond, byNetwork=FALSE){
  
  #outputRegincond: Output object of running RegulationInCondition function
  #byNetwork: By faulta, FALSE. If TRUE plots the results in a network

  regulators<-outputRegincond$GlobalRegulators
  if (length(regulators)==0){
    stop("ERROR:No global regulator was identified for this condition")
  }
  if (byNetwork){
    
    #Take group column in RegulationPerCondition
    df<-outputRegincond$RegulationInCondition[grepl(paste(regulators, collapse = "|"), outputRegincond$RegulationInCondition$regulator),]
    
    #Create the graph
    mygraph = igraph::graph_from_data_frame(df, directed=F)
    
    odf<-df[,c(2,3)]
    odf<-rbind(odf, data.frame('regulator'=unique(df$targetF),'omic'=rep('targetF',length(unique(df$targetF)))))
    odf<-unique(odf)
    rownames(odf)<-odf$regulator
    
    mygraph<-igraph::set.vertex.attribute(mygraph,'omic', index = igraph::V(mygraph), value = odf[igraph::V(mygraph)$name,]$omic)
    mygraph<-igraph::set.edge.attribute(mygraph, 'sign', index = igraph::E(mygraph), value = df[,4])
    igraph::E(mygraph)$sign<-df[,4]
    
    igraph::plot.igraph(mygraph,vertex.label.cex= 0.4, vertex.size = 4,vertex.color=as.factor(igraph::V(mygraph)$omic),edge.color = ifelse(igraph::E(mygraph)$sign >0, "blue", "red"), main='targetF - Global regulators network')
    legend("topright", legend = unique(igraph::V(mygraph)$omic), col = categorical_pal(length(unique(as.factor(igraph::V(mygraph)$omic)))), pch = 16, cex = 1.5, bty = "n")
    
    
  }else{
    
    # Identify the target features they regulate
    targetFs<-unique(outputRegincond$RegulationInCondition[grepl(paste(regulators, collapse = "|"), outputRegincond$RegulationInCondition$regulator),1])
    
    #Create the matrix 
    gen_reg<-matrix(0, nrow= length(targetFs),ncol=length(regulators))
    colnames(gen_reg)=regulators
    rownames(gen_reg)=targetFs
    
    for ( i in 1:nrow(gen_reg)){
      #Get all the potential regulators of the target feature
      potential_regulator <- unlist(sapply(outputRegincond$arguments$associations, function(x) getallreg(x,targetFs[i])),use.names = FALSE)
      
      #Use NA for any potential regulator
      gen_reg[i,intersect(potential_regulator,regulators)]<-NA
      for ( j in 1:ncol(gen_reg)){
        if ( length(intersect(which(outputRegincond$RegulationInCondition$targetF==targetFs[i]),which(outputRegincond$RegulationInCondition$regulator==regulators[j])) )!=0){
          gen_reg[i,j] = outputRegincond$RegulationInCondition[intersect(which(outputRegincond$RegulationInCondition$targetF==targetFs[i]),which(outputRegincond$RegulationInCondition$regulator==regulators[j])),4]
        }
      }
    }
    
    df<-data.frame(targetFs=rep(rownames(gen_reg), times=ncol(gen_reg)),
                   regulators=rep(colnames(gen_reg),each = nrow(gen_reg)),
                   value=as.vector(gen_reg),
                   color= ifelse(as.vector(gen_reg)>0, '#5577FF',ifelse(as.vector(gen_reg)<0, '#FF7755','#FFFFFF')))
    df$color[is.na(df$color)] <- '#aaaaaa'
    
    ggplot2::ggplot(data = df, aes(x = regulators, y = targetFs, fill = color)) +
      geom_tile(color = "black",lwd = 0.5,  linetype = 1) +
      scale_fill_manual(values = c("#5577FF", "#aaaaaa", "#FF7755",  "#FFFFFF"), 
                        name = "Legend",
                        labels = c('Activator','Potential','Repressor','Not potential')) +
      labs(title = paste0("targetF - Global Regulators \n correlation plot in ",colnames(outputRegincond$RegulationInCondition)[4]," condition  \n"), 
           x = "Global regulators", y = "targetFs") +
      theme(plot.title = element_text(hjust = 0.5, colour = "black"), 
            axis.title.x = element_text(face="bold", colour="black", size = 8),
            axis.title.y = element_text(face="bold", colour="black", size = 8),
            axis.text.x = element_text(size = 6, angle = 45, hjust = 1),  # Adjust size for x-axis text
            axis.text.y = element_text(size = 6),
            legend.title = element_text(face="bold", colour="black", size = 10),legend.position = 'right') 
    
  }
  
}

## Network creation -------


#' networkMORE
#'
#' \code{networkMORE} Function to be applied to RegulationPerConidtion function output.
#' 
#' @param outputRegpcond Output object of RegulationPerCondition applied to MORE main function.
#' @param cytoscape TRUE for plotting the network in Cytoscape. FALSE to plot the network in R. 
#' @param group1 Name of the group to take as reference in the differential network creation. It also can be used for creating networks of a specific group. If it is not provided the networks of all conditions will be plotted. By default, NULL.
#' @param group2 Name of the group to compare to the reference in the differential network creation. By default, NULL.
#' @param pc Percentile to consider to plot the most affecting regulators into the target omic. It must be a value comprissed between 0 and 1. By default, 0.
#' @param pathway If provided, the function will print the regulatory network involved in the specified pathway instead of the entire regulatory network. By default, NULL.
#' @param annotation Annotation matrix with target features in the first column, GO terms in the second and GO term description in the third. Only necessary when a specific pathway has to be plotted. By default, NULL.
#' @param save If TRUE a gml extension network is saved when cytoscape = FALSE. By default, FALSE.
#' @return Plot of the network induced from more.
#' @export


networkMORE <- function(outputRegpcond, cytoscape = TRUE, group1 = NULL, group2 = NULL, pc = 0, pathway= NULL, annotation = NULL, save=FALSE) {
  
  create_graph <- function(df,pc) {
    
    #Remove rows with 0 coef
    df = df[df[,4] != 0, ]
    
    #Select regulations with effects on response in the selected percentile
    qc = quantile(abs(df[,4]),pc)[[1]]
    df = df[which(abs(df[,4])>qc),,drop=FALSE]
    
    #Data.frame of that network
    nodes = data.frame(id = c(unique(df[,'targetF']),unique(df[,'regulator'])),
                       omic = c(rep('targetF',length(unique(df[,'targetF']))),unique(df[,c('regulator','omic')])[,2]))
    
    #Save only four digits as it cannot be loaded greater in Cytoscape
    interactions = data.frame(from = df[,'targetF'], to = df[,'regulator'],
                              coef = round(df[,4],digits = 5), sign = ifelse(df[,4]>0,'p','n'))
    
    ig = igraph::graph_from_data_frame(interactions, vertices = nodes, directed = FALSE)
    
    return(ig)
  }
  
  create_network <- function(mygraph, group_names,diff) {
    
    cy_network <- RCy3::createNetworkFromIgraph(mygraph, group_names)
    
    #Set node color and generate a color palette
    omic_c <- unique(factor(igraph::vertex_attr(mygraph)$omic))
    num_unique <- length(omic_c)
    color_palette <- colorbiostat(num_unique)
    RCy3::setNodeColorMapping('omic', table.column.values = omic_c ,colors = color_palette, mapping.type='d')
    
    nshaps <-setdiff(RCy3::getNodeShapes(), c("TRIANGLE", "DIAMOND","RECTANGLE"))[1:num_unique]
    if(any(grepl("tf", omic_c, ignore.case = TRUE))){
      i=grep('tf', omic_c, ignore.case = TRUE)
      nshaps[i]<-'TRIANGLE'
    }
    if(any(grepl("mirna", omic_c, ignore.case = TRUE))){
      i=grep('mirna', omic_c, ignore.case = TRUE)
      nshaps[i]<-'DIAMOND'
    }
    if('targetF'%in% omic_c){
      i=grep('targetF', omic_c)
      nshaps[i]<-'RECTANGLE'
    }
    RCy3::setNodeShapeMapping('omic', table.column.values = omic_c, shapes = nshaps )
    RCy3::setEdgeColorMapping('sign', c('n','p'), c('#FF3333','#5577FF'),mapping.type='d')
    if(diff){
      RCy3::setEdgeLineStyleMapping('line',c('s','e','d','v','p'),c('SOLID','EQUAL_DASH','DOT','VERTICAL_SLASH','PARALLEL_LINES'))
    } 
  }
  
  DifLineType <- function(df) {
    df[, 7:8] <- NA
    
    # Assign column 8 based on the sign of Group2
    df[, 8] <- ifelse(sign(df[, 5]) == -1, 'n', 'p')
    
    # Assign line type
    df[, 7] <- ifelse(df[, 4] == 0, 'v', 
                      ifelse(df[, 5] == 0, 'd', 
                             ifelse(sign(df[, 4]) == sign(df[, 5]), 
                                    ifelse(abs(df[, 4]) > abs(df[, 5]), 'e', 's'), 
                                    'p')))
    
    # Assign column 8 based on the sign of Reference if Group2 = 0
    df[df[, 5] == 0, 8] <- ifelse(sign(df[df[, 5] == 0, 4]) == -1, 'n', 'p')
    
    return(df)
  }
  
  if(!is.null(pathway)){
    if(is.null(annotation)){stop('No annotation matrix was provided')}
    
    targetFs=annotation[which(annotation[,3] == pathway),1]
    if(length(targetFs)==0){stop('No target feature was considered in that pathway')}
    outputRegpcond = outputRegpcond[c(outputRegpcond$targetF %in% targetFs),,drop=FALSE]
  }
  
  if (cytoscape) {
    if (is.null(group1) && is.null(group2)) {
      
      ngroups <- grep('Group', colnames(outputRegpcond))
      #Create as many networks as groups
      for (i in 1:length(ngroups)) {
        
        df = outputRegpcond[, c(1, 2, 3, ngroups[i])]
        
        ig = create_graph(df,pc)
        if(i == 1){
          create_network(ig, colnames(outputRegpcond)[ngroups[i]],diff = FALSE)
        }  else{
          RCy3::createNetworkFromIgraph(ig, colnames(outputRegpcond)[ngroups[i]])
        }
        
      }
      
    } else if (is.null(group2)){
      ngroup <- grep(group1, colnames(outputRegpcond))
      df = outputRegpcond[, c(1, 2, 3, ngroup)]
      ig = create_graph(df,pc)
      create_network(ig, colnames(outputRegpcond)[ngroup],diff = FALSE)
    } else {
      #Look for the groups to consider
      gr1 <- grep(group1, colnames(outputRegpcond))
      gr2 <- grep(group2, colnames(outputRegpcond))
      
      if (length(gr1) != 1 || length(gr2) != 1 || gr1 == gr2){stop("ERROR: group1 and group2 should be different names of groups to compare")}
      #Create the differential coefficient and the indicator of sign change
      df <- outputRegpcond[, c(1,2,3,gr1,gr2)]
      df[, 6] = df[, 5] - df[, 4]
      #Remove rows with same effect
      df = df[df[,6] != 0, ]
      #Select regulations with effects on response in the selected percentile
      qc = quantile(abs(df[,6]),pc)[[1]]
      df = df[which(abs(df[,6])>qc),,drop=FALSE]
      #Add the line type and sign
      df = DifLineType(df)
      
      #Data.frame of that network
      nodes = data.frame(id = c(unique(df[,'targetF']),unique(df[,'regulator'])),
                         omic = c(rep('targetF',length(unique(df[,'targetF']))),unique(df[,c('regulator','omic')])[,2]))
      
      interactions = data.frame(from = df[,'targetF'], to = df[,'regulator'],
                                sign = df[,8], line = df[,7])
      
      ig = igraph::graph_from_data_frame(interactions, vertices = nodes, directed = FALSE)
      
      create_network(ig, paste0(colnames(outputRegpcond)[gr2], '-', colnames(outputRegpcond)[gr1]) ,diff = TRUE)
    }
    
  } else {
    
    if (is.null(group1) && is.null(group2)) {
      
      ngroups <- grep('Group', colnames(outputRegpcond))
      #Create as many networks as groups
      for (i in 1:length(ngroups)) {
        #Data.frame of that network
        df = outputRegpcond[, c(1, 2, 3, ngroups[i])]
        ig = create_graph(df,pc)
        
        if(save){
          igraph::write_graph(ig, format = 'gml', file = paste0(colnames(outputRegpcond)[ngroups[i]], '.gml'))
        } else{
          igraph::plot.igraph(ig, vertex.label.cex = 0.3, vertex.size = 3,
                              vertex.color = as.factor(igraph::vertex_attr(ig)$omic),
                              edge.color = ifelse(igraph::edge_attr(ig)$sign == 'p', 'blue','red'))
        }
        
      }
      
    } else if(is.null(group2)){
      ngroup <- grep(group1, colnames(outputRegpcond))
      df = outputRegpcond[, c(1, 2, 3, ngroup)]
      ig = create_graph(df,pc)
      if(save){
        igraph::write_graph(ig, format = 'gml', file = paste0(colnames(outputRegpcond)[ngroup], '.gml'))
      } else{
        igraph::plot.igraph(ig, vertex.label.cex = 0.3, vertex.size = 3,
                            vertex.color = as.factor(igraph::vertex_attr(ig)$omic),
                            edge.color = ifelse(igraph::edge_attr(ig)$sign == 'p', 'blue','red'))
      }
      
    } else {
      gr1 <- grep(group1, colnames(outputRegpcond))
      gr2 <- grep(group2, colnames(outputRegpcond))
      
      if (length(gr1) != 1 || length(gr2) != 1 || gr1 == gr2){stop("ERROR: group1 and group2 should be different names of groups to compare")}
      #Create the differential coefficient and the indicator of sign change
      df <- outputRegpcond[, c(1,2,3,gr2,gr1)]
      df[, 6] = df[, 4] - df[, 5]
      #Remove rows with same effect
      df = df[df[,6] != 0, ]
      #Select regulations with effects on response in the selected percentile
      qc = quantile(abs(df[,6]),pc)[[1]]
      df = df[which(abs(df[,6])>qc),,drop=FALSE]
      #Add the line type and sign
      df = DifLineType(df)
      
      #Data.frame of that network
      nodes = data.frame(id = c(unique(df[,'targetF']),unique(df[,'regulator'])),
                         omic = c(rep('targetF',length(unique(df[,'targetF']))),unique(df[,c('regulator','omic')])[,2]))
      
      interactions = data.frame(from = df[,'targetF'], to = df[,'regulator'],
                                sign = df[,8], line = df[,7])
      
      ig = igraph::graph_from_data_frame(interactions, vertices = nodes, directed = FALSE)
      
      if(save){
        igraph::write_graph(ig, format = 'gml', file = paste0(colnames(outputRegpcond)[gr2],'-',colnames(outputRegpcond)[gr1], '.gml'))
      } else{
        igraph::plot.igraph(ig, vertex.label.cex = 0.3, vertex.size = 3,
                            vertex.color = as.factor(igraph::vertex_attr(ig)$omic),
                            edge.color = ifelse(igraph::edge_attr(ig)$sign == 'p', 'blue','red'),
                            edge.lty = ifelse(igraph::edge_attr(ig)$line=='s','solid','dashed'))
      }
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
  if(any(class(resultat)=='matrix')){
    return(data.frame(resultat,
                      "adjPval" = p.adjust(as.numeric(resultat[,"pval"]), method = p.adjust.method),
                      stringsAsFactors = F))
  } else {
    return(data.frame(t(c(resultat,
                      "adjPval" = as.numeric(resultat["pval"])))))
  }
  
}


#' oraMORE
#'
#' \code{oraMORE} Function to be applied to RegulationInCondition function output.
#' 
#' @param output Output object of running more
#' @param outputRegincond Output object of running RegulationInCondition function
#' @param byHubs Indicates whether to perform the ORA for the Hub target features, TRUE, or for the target features regulated by the global regulators, FALSE. By default, TRUE.
#' @param byOmic If provided (it must follow the same nomenclature that in regulatoryData), it performs the ORA to the regulators of the specified omic. Incompatible with other methodologies specified in byHubs parameter. By default, NULL.
#' @param annotation Annotation matrix with target features in the first column, GO terms in the second and GO term description in the third
#' @param alpha The adjusted pvalue cutoff to consider
#' @param p.adjust.method One of holm, hochberg, hommel, bonferroni, BH, BY, fdr or none
#' @param parallel parallel If FALSE, MORE will be run sequentially. If TRUE, MORE will be run using parallelization with as many cores as the available ones minus one so the system is not overcharged. If the user wants to specify how many cores they want to use, they can also provide the number of cores to use in this parameter.
#' 
#' @return Plot of the network induced from more.
#' @export

oraMORE = function(output, outputRegincond, byHubs = FALSE, byOmic = NULL, annotation, alpha = 0.05,
                   p.adjust.method = "fdr", parallel = FALSE) {
  
  # output: Output object of running more function 
  # outputRegincond: Output object of running RegulationInCondition function
  # annotation: Annotation matrix with target features in the first column, GO terms in the second and a description in the third
  
  #Take the reference to compare the set of genes
  reference = rownames(output$arguments$targetData)
  
  annotation = annotation[annotation[,1] %in% reference,]
  annotDescr = unique(annotation[,2:3])
  rownames(annotDescr) = annotDescr[,1]
  
  if(byHubs){
    
    test = outputRegincond$HubTargetF
    notTest = setdiff(reference, test)
    resuRegu = ReguEnrich1regu(test, notTest, annotation, p.adjust.method = p.adjust.method)
    myresults = resuRegu[which(resuRegu$adjPval < alpha),,drop=FALSE]
    myresults = data.frame( "termDescr" = annotDescr[myresults[,1],2],
                            myresults, stringsAsFactors = FALSE)
  } else{
    
    if(!is.null(byOmic)){
      regulators = outputRegincond$RegulationInCondition$regulator[which(outputRegincond$RegulationInCondition$omic==byOmic)]
    } else{
      regulators = outputRegincond$GlobalRegulators
    }
    
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
      myresults <- furrr::future_map(1:length(regulators),
                                     ~ORA.i(regulators[.],outputRegincond, reference, annotation, p.adjust.method, annotDescr,alpha),
                                     .progress = TRUE )
    } else{
      myresults <- purrr::map(1:length(regulators),
                              ~ORA.i(regulators[.],outputRegincond, reference, annotation, p.adjust.method, annotDescr,alpha),
                              .progress = TRUE )
    }
    future::plan('sequential')
    names(myresults) = regulators
    
  }
  
  return(myresults)
}


ORA.i = function(regulator, outputRegincond, reference, annotation, p.adjust.method, annotDescr,alpha){
  
  test = unique(as.character(outputRegincond$RegulationInCondition[which(outputRegincond$RegulationInCondition[,"regulator"] == regulator),"targetF"]))
  notTest = setdiff(reference, test)
  resuRegu = ReguEnrich1regu(test, notTest, annotation, p.adjust.method = p.adjust.method)
  resuRegu = resuRegu[which(resuRegu$adjPval < alpha),,drop=FALSE]
  if (nrow(resuRegu) > 0) {
    resuRegu = data.frame("regulator" = regulator, "termDescr" = annotDescr[resuRegu[,1],2],
                          resuRegu, stringsAsFactors = FALSE)
  } else {
    
    resuRegu = NULL
  }
  
  return(resuRegu)
  
}

#' gseaMORE
#'
#' \code{gseaMORE} Function to be applied to RegulationInCondition function output.
#' 
#' @param outputRegincond Output object of running RegulationInCondition function
#' @param outputRegincond2 Output object of running RegulationInCondition function for other group different to the previous. By default, NULL.
#' @param annotation Annotation matrix with target features in the first column, GO terms in the second and GO term description in the third
#' @param alpha The adjusted pvalue cutoff to consider
#' @param p.adjust.method One of holm, hochberg, hommel, bonferroni, BH, BY, fdr or none
#' 
#' @return Plot of the network induced from more.
#' @export

gseaMORE<-function(outputRegincond, outputRegincond2 = NULL, annotation, alpha = 0.05,
                    p.adjust.method = "fdr"){
  
  #outputRegincond: Output object of running RegulationInCondition
  #outputRegincond2: Output object of running RegulationInCondition for other group to which compare the reference. By default, NULL.
  #annotation: Annotation matrix with target features in the first column, GO terms in the second and GO term description in the third
  #alpha: The adjusted pvalue cutoff to consider
  #p.adjust.method: One of holm, hochberg, hommel, bonferroni, BH, BY, fdr or none
  term2gene_bp<-annotation[,c(2,1)]
  term2name_bp <- unique(annotation[,c(2,3)])
  if(is.null(outputRegincond2)){
    #Store the target features in decreasing order

    geneList<-table(outputRegincond$RegulationInCondition$targetF)
    geneList = sort(geneList, decreasing = TRUE)
    
    selected_genes <- names(geneList)
    counts <- as.numeric(geneList)
    geneList <- setNames(counts, selected_genes)
    
    y <- clusterProfiler::GSEA(geneList = geneList, TERM2GENE = term2gene_bp, TERM2NAME = term2name_bp , pvalueCutoff = alpha, pAdjustMethod = p.adjust.method)
    clusterProfiler::dotplot(y, split = '.sign') 
    
  } else{
    
    geneList<-as.data.frame(table(outputRegincond$RegulationInCondition$targetF))
    
    geneList2<-as.data.frame(table(outputRegincond2$RegulationInCondition$targetF))
    
    #Create a data frame 
    merged_df = merge(geneList, geneList2, by = 'Var1', all = TRUE)
    merged_df[is.na(merged_df)] = 0
    rownames(merged_df)<-merged_df[,1]
    merged_df<-merged_df[,-1]
    #Create the score for the GSEA
    merged_df[,3]<-merged_df[,2]-merged_df[,1]
    
    geneList <- setNames(merged_df[,3], rownames(merged_df))
    geneList <- sort(geneList, decreasing = TRUE)
    
    y <- clusterProfiler::GSEA(geneList = geneList, TERM2GENE = term2gene_bp, TERM2NAME = term2name_bp , pvalueCutoff = alpha, pAdjustMethod = p.adjust.method)
    clusterProfiler::dotplot(y, split = '.sign') + facet_grid(.~.sign, labeller = as_labeller(c(activated = gsub('Group_','',colnames(outputRegincond2$RegulationInCondition)[4]), suppressed = gsub('Group_','',colnames(outputRegincond$RegulationInCondition)[4]))))
    
  }
  
  
  return(y)
}

#' upsetMORE
#'
#' \code{upsetMORE} Function to be applied to a list of RegulationInCondition function outputs.
#' 
#' @param listRegincond List with the output of RegulationInCondition for the conditions to be compared. It must be a named list with the condition as names. 
#' @param byHubs Indicates whether to create the UpSet for the Hub target features, TRUE, or for the Global Regulators, FALSE. By default, TRUE.
#' 
#' @return UpSet plot of the hub target features/global regulators induced from more.
#' @export


upsetMORE = function(listRegincond, byHubs = TRUE){
  
  #Verify required packages are installed
  
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop("Package 'UpSetR' is required but not installed. Please install it to use this function.")
  }
  if (!requireNamespace("ComplexUpset", quietly = TRUE)) {
    stop("Package 'ComplexUpset' is required but not installed. Please install it to use this function.")
  }
  
  #Create the lists according to the plots to create
  if(byHubs){
    x = lapply(listRegincond, function(x) x$HubTargetF)
  } else {
    x = lapply(listRegincond, function(x) x$GlobalRegulators)
  }
  
  data = UpSetR::fromList(x)
  
  #Ask for biostatomics colors
  conditions = names(data)
  num_unique = length(conditions)+1
  color_palette = colorbiostat(num_unique)
  custom_colors = setNames(color_palette[-1], conditions)
  
  queries = lapply(conditions, function(cond) {
    ComplexUpset::upset_query(set = cond, fill = custom_colors[cond]) })
  
  #Create the Upset
  ComplexUpset::upset(
    data,
    conditions,
    queries=queries,
    base_annotations=list(
      'Intersection size'=(
        ComplexUpset::intersection_size(
          bar_number_threshold=1,  # show all numbers on top of bars
          width=0.5,   # reduce width of the bars
        )
        # add some space on the top of the bars
        + scale_y_continuous(expand=expansion(mult=c(0, 0.05)))
        + theme(
          # hide grid lines
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          # show axis lines
          axis.line=element_line(colour='black')
        )
      )
    ),
    stripes=ComplexUpset::upset_stripes(
      geom=geom_segment(size=12),  # make the stripes larger
      colors=c('grey95', 'white')
    ),
    matrix=ComplexUpset::intersection_matrix(
      geom=geom_point(
        shape='circle filled',
        size=3,
        stroke=0
      )
    ),
    set_sizes=(
      ComplexUpset::upset_set_size(geom=geom_bar(width=0.4))
      + theme(
        axis.line.x=element_line(colour='black'),
        axis.ticks.x=element_line()
      )
    ),
    sort_sets=FALSE,
    sort_intersections='descending'
  )
  
}


