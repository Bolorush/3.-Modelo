#------------------------------------------------------------------------------#
## implementar fevd en contexto bayesiano ----
### matriz de coeficientes, matrices companion ----
fevd.bayes <- function(datos, p, nirf,inf,sup){
  ## implementar fevd en contexto bayesiano ----
  ### matriz de coeficientes, matrices companion ----
  # parametros necesarios
  p               = p                                       # rezagos
  k               = dim(datos)[2]                           # variables
  tt              = dim(datos)[1]                           # periodos
  nirf            = nirf                                    # periodos h irf y fevd
  s               = dim(posterior[["posterior"]][["A"]])[3] # samplin posterior 
  mat.A           = posterior[["posterior"]][["A"]]         # matrices A
  mat.B           = posterior[["posterior"]][["B"]]         # matrices B
  matriz.den.irf  = matrix(0,nrow = k*k*nirf,ncol = s)      # matriz densidades irf
  matriz.den.fevd = matrix(0,nrow = k*k*nirf,ncol = s)      # matriz densidades fevd
  
  for(draw in 1:s){
    # construccion de matriz companion
    matrices.A <- mat.A[,,draw][,-(k+1)]               # matrices de coeficientes
    
    # matriz companion Y = AYt-1 + V
    matriz.i  <- diag(ncol(matrices.A) -nrow(matrices.A))
    matriz.0  <- matrix(0,nrow(matriz.i),(ncol(matrices.A)-ncol(matriz.i)))
    aux.1     <- cbind(matriz.i,matriz.0)
    companion <- rbind(matrices.A,aux.1)            # matriz companion
    
    # matriz estructural S para determinar funciones Phi
    S = rbind(diag(nrow(matrices.A)),matriz.0)      # matriz S = I
    
    ### cálculo de las matrices impulso respuesta no ortogonales Phi----
    matriz.irf <- array(0,dim = c(k,k,nirf))        # matrices Phi
    
    # inicializar la matriz companion y las matrices irf
    Bpowtm1               <- diag(nrow(companion))
    matriz.irf0           <- Bpowtm1%*%S 
    matriz.irf[1:k,1:k,1] <- matriz.irf0[1:k,1:k]   # primera Phi
    
    # bucle para cálculo de cada Phi
    for (t in 2:nirf) {
      Yimpt   <- Bpowtm1%*%S           # irf en el tiempo t del choque i
      Bpowtm1 <- Bpowtm1%*%companion   # actualizo companion = companion t * companion t-1
      matriz.irf[1:k,1:k,t] <-  Bpowtm1[1:k,1:k]
    }
    
    #______________________________________________________________________________#
    # matrices impulso respuesta ortogonales Psi, punto de partida es B normalizada
    B <- mat.B[,,draw]
    
    # estandarizar la matriz B para que los coeficientes de la diagonal sean 1----
    for(l in 1:k){
      diagonal <- B[l,l]
      B[l,] <- B[l,]/diagonal
    }
    
    # matrices impulso respuesta ortogonales
    matriz.irf.o <- array(0,dim = c(k,k,nirf))
    
    # inicializar las matrices irf ortogonales, la primera es B
    matriz.irf.o[1:k,1:k,1] <- B
    
    # bucle, Psi = Phi*S
    for (t in 2:nirf) {
      matriz.irf.o[1:k,1:k,t] <- matriz.irf[1:k,1:k,t]%*%B  
    }
    
    #-----------------CALCULO DE DESCOMPOSICIÓN DE LA VARIANZA-----------------#
    # calculo de las MSE inicial, la diagonal de la multiplicacion Psi*Psi'
    mse.1 <- diag(matriz.irf.o[,,1]%*%t(matriz.irf.o[,,1]))
    
    # el numerador esta compuesto por los elementos de Psi al cuadrado
    num.1 <- matriz.irf.o[,,1]^2
    
    # la FEVD es numerador sobre mse
    fevd.1 <- num.1/mse.1
    # fevd en formato tabla
    fevd.tabla <- matrix(0,nirf,k*k)
    # primer elemento 
    fevd.tabla[1,] <- matrix(t(fevd.1), nrow=1)
    # las fevd se obtiene acumulando el valor previo
    for(t in 2:nirf){
      mse.1 <- mse.1 + diag(matriz.irf.o[,,t]%*%t(matriz.irf.o[,,t]))
      num.1 <- num.1 + matriz.irf.o[,,t]^2
      fevd.tabla[t,] <-  matrix(t(num.1/mse.1), nrow=1)
    }
    
    #------------------------------------------------------------------------------# 
    ### densidades de las fevd ----
    matriz.den.fevd[,draw] <- matrix(fevd.tabla, ncol = 1)
  }
  
  ### obtener mediana, intervalos creibles y sd por fila para todas las fevd ----
  
  matriz.est.fevd <- t(apply(matriz.den.fevd, 1, quantile, probs=c(inf,0.5,sup)))
  matriz.est.fevd <- cbind(matriz.est.fevd, apply(matriz.den.fevd, 1, sd) )
  
  ### obtener las fevd de cada variable ----
  nombres.1 <- nombres.fevd(nombres)
  
  matrices.fevd <- list()                                      
  for(m in 1:(k*k)){
    matrices.fevd[[m]] <- matrix(0,nrow = nirf,ncol = 4)      # inicializo las fevd
  }
  
  # extraccion de cada fevd, se pone el nombre de cada fevd
  fin  <- 0
  for (m in 1:(k*k)){
    fin  <- fin + nirf
    inicio <- fin - nirf +1
    matrices.fevd[[m]] <- matriz.est.fevd[inicio:fin,]
    colnames(matrices.fevd[[m]]) <- nombres.1[m,]
  }
  
  ### tablas de cada fevd ----
  
  # crear lista para almacenar tablas
  tablas.fevd = list()
  
  # bucle para iterar por cada fevd estimada
  fin  <- 0
  for (m in 1:k){
    fin  <- fin + k
    inicio <- fin - k +1
    tablas.fevd[[m]] <- matrices.fevd[inicio:fin]
  }
  
  return(tablas.fevd)
  
}
#------------------------------------------------------------------------------#