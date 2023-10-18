# DESARROLLO DE UN MODELO VAR CON REGIMEN CAMBIANTE ----

rm(list = ls())

#------------------------------------------------------------------------------#
## Cargar los paquetes y funciones definidas ----
library(bsvars)
library(matlib)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(lubridate)
library(plyr)
library(readxl)

# funciones personalizadas
source("./inputs/funciones/nombres.R")     # nombres para las irf y fevd
source("./inputs/funciones/irf.bayes.R")   # irf en contexto bayesiano
source("./inputs/funciones/fevd.bayes.R")  # fevd en contexto bayesiano
source("./inputs/funciones/xi.bayes.R")    # prob transicion en contexto bayesiano

#------------------------------------------------------------------------------#
## Apertura y procesamiento de la base de datos ----

# apertura, empieza en 2006q1 y termina en 2023q2
bdd <- read_excel("inputs/bdd.xlsx", range = "trimestral!A1:N71",col_names = T)










