#------------------------------------------------------------------------------#
# funcion de nombres de las irf
nombres.irf <- function(nombres){
  aux <- c()
  for(j in 1: length(nombres)){
    for(k in 1: length(nombres)){
      aux <- c(aux,paste0(nombres[j],"-->",nombres[k]))
    }
  }
  aux.1 <-  paste("Inf.",aux,sep="")          # incluye inf, med y sup          
  aux.2 <-  paste("Med.",aux,sep="")
  aux.3 <-  paste("Sup.",aux,sep="")
  aux.4 <-  cbind(aux.1,aux.2,aux.3)
  return(aux.4)
}
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
# funcion de nombres de las fevd
nombres.fevd <- function(nombres){
  aux <- c()
  for(j in 1: length(nombres)){
    for(k in 1: length(nombres)){
      aux <- c(aux,paste0(nombres[j],"-->",nombres[k]))
    }
  }
  aux.1 <-  paste("Inf.",aux,sep="")          # incluye inf, med y sup          
  aux.2 <-  paste("Med.",aux,sep="")
  aux.3 <-  paste("Sup.",aux,sep="")
  aux.4 <-  paste("sd." ,aux,sep="")
  aux.5 <-  cbind(aux.1,aux.2,aux.3,aux.4)
  return(aux.5)
}
#------------------------------------------------------------------------------#
