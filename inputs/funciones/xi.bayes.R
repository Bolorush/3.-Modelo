#------------------------------------------------------------------------------#
## implementar proabilidades de transición de estado en contexto bayesiano ----
xi.bayes <- function(xi){
  
  # dimensiones de cada array
  m.xi   <- dim(xi)[1]                    # número de estados 
  t.xi   <- dim(xi)[2]                    # periodos de tiempo
  s.xi   <- dim(xi)[3]                    # número de draws posterior
  den.xi <- matrix(0,nrow = m.xi*t.xi, ncol = s.xi)   #matriz de densidades xi
  
  # llenar la matriz de densidades
  for(draw in 1: s.xi){
    den.xi[,draw] <- matrix(t(xi[,,draw]),ncol = 1)
  }
  # extraer estadísticos requeridos
  inf=0.16; sup=0.84
  mean.xi <- apply(den.xi, 1, mean )
  sd.xi   <- apply(den.xi, 1, sd )
  matriz.est.xi <- cbind(mean.xi,sd.xi)
  # matriz de estados con estadisticos
  estados.xi <- vector()
  
  fin  <- 0
  for(m in 1:m.xi){
    fin  <- fin + t.xi
    inicio <- fin - t.xi +1
    estados.xi <- cbind(estados.xi,matriz.est.xi[inicio:fin,1])
  }
  
  # asignar nombres por cada estado  
  nombres.xi <- vector()
  for(m in 1:m.xi){
    nombres.xi <- c(nombres.xi,paste("Estado",m)) 
  }
  colnames(estados.xi) <- nombres.xi
  
  # dataframe para grafico
  estados.xi <- data.frame(estados.xi,periodo.q)
  df_c <- estados.xi %>%
    dplyr:: select(periodo.q, colnames(estados.xi)[1:m]) %>%
    gather(key = "variable", value = "value", -periodo.q)
  
  # 
  grafico <- ggplot(data = df_c, aes(x = periodo.q, y = value))+
    geom_line() + theme_minimal() +
    scale_y_continuous(limits=c(0,1),breaks = seq(0, 1, 0.1))+
    facet_wrap(~variable, ncol = 1) + ylab("") +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(colour = "gray", size = 1)) 
  
  return(grafico)
}
#------------------------------------------------------------------------------#