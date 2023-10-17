#------------------------------------------------------------------------------#
## implementar irf en contexto bayesiano ----
### matriz de coeficientes, matrices companion ----
irf.bayes <- function(datos, p, nirf,mog,inf,sup){
  # parametros necesarios
  p              = p                                       # rezagos
  k              = dim(datos)[2]                           # variables
  tt             = dim(datos)[1]                           # periodos
  nirf           = nirf                                    # periodos h irf y fevd
  s              = dim(posterior[["posterior"]][["A"]])[3] # samplin posterior 
  mat.A          = posterior[["posterior"]][["A"]]         # matrices A
  mat.B          = posterior[["posterior"]][["B"]]         # matrices B
  matriz.den.irf = matrix(0,nrow = k*k*nirf,ncol = s)      # matriz densidades irf
  
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
    
    #------------------------------------------------------------------------------# 
    ### irf en formato de tabla ----
    irf.tabla <- matrix(0,nirf,k*k)
    for(t in 1:nirf){
      irf.tabla[t,] <-  matrix(matriz.irf.o[,,t],nrow = 1)
    }
    
    matriz.den.irf[,draw] <- matrix(irf.tabla, ncol = 1)
  }
  
  ### obtener mediana e intervalos creibles por fila para todas las irf ----
  matriz.est.irf <- t(apply(matriz.den.irf, 1, quantile, probs=c(inf,0.5,sup)))
  
  ### obtener las irf de cada variable ----
  nombres <- colnames(us_fiscal_lsuw)
  nombres <- nombres.irf(nombres)
  
  matrices.irf <- list()                                      
  for(m in 1:(k*k)){
    matrices.irf[[m]] <- matrix(0,nrow = nirf,ncol = 3)      # inicializo las irf
  }
  
  # extraccion de cada irf, se pone el nombre de cada irf
  fin  <- 0
  for (m in 1:(k*k)){
    fin  <- fin + nirf
    inicio <- fin - nirf +1
    matrices.irf[[m]] <- matriz.est.irf[inicio:fin,]
    colnames(matrices.irf[[m]]) <- nombres[m,]
  }
  
  ### gráficos de cada irf ----
  
  # crear lista para almacenar graficos
  graficos.irf = list()
  
  # crear variable de fecha para grafico temporal
  periodo <- as.Date(1:nirf)
  
  # bucle para iterar por cada irf estimada
  for(g in 1: (k*k)){
    # uno con variable de fechas
    grafico.irf <- as.data.frame(cbind(periodo,matrices.irf[[g]])) 
    
    # modifico el dataframe para el grafico
    df <- grafico.irf %>%
      dplyr:: select(colnames(grafico.irf)[1], colnames(grafico.irf)[2], colnames(grafico.irf)[3],colnames(grafico.irf)[4]) %>%
      gather(key = "variable", value = "value", -periodo)
    
    plot.irf <- ggplot(df, aes(x = periodo, y = value)) + 
      geom_line(aes(color = variable, linetype = variable), size = 0.8, alpha=0.8) +
      scale_color_manual(values = c("#00696b", "#E7B800","#00696b")) +
      scale_linetype_manual(values = c("dotted","solid","dotted"))  +
      theme_minimal() + ggtitle(paste("Impulso Respuesta ",sub("Med.", "",colnames(grafico.irf[3])), sep="")) + 
      xlab("") + ylab("") + labs(color = "") + theme(legend.position = "none")
    #+ scale_y_continuous(breaks = seq(-0.2, 1, 0.1))
    
    graficos.irf[[g]] <-  plot.irf
  
    }
  
  if(mog==T){return(matrices.irf)} else {return(graficos.irf)}
  
}
#------------------------------------------------------------------------------#