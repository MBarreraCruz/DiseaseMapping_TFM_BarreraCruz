# =============================================================================
# Coregionalization.R
# Modelo de corregionalización multivariante — ordenación original
# -----------------------------------------------------------------------------
# TFM en Bioestadística — María Barrera Cruz
# Universidad de Valencia, Facultad de Matemáticas, 2026
# Tutor: Miguel Ángel Martínez-Beneito (Dept. Estadística e Investigación Operativa, UV)
# -----------------------------------------------------------------------------
# Descripción:
#   Ajuste del modelo de corregionalización bayesiano para 78 causas de
#   mortalidad prematura en 542 municipios de la Comunitat Valenciana
#   (2000–2023), con la ordenación original de enfermedades (referencia).
#
#   El mismo código, descomentando el bloque correspondiente en la sección
#   "ORDEN DE LAS ENFERMEDADES", sirve para las ordenaciones directa, inversa
#   y las 50 aleatorias (53 versiones en total).
#
#
# NOTA: este script no es autoejecutable en ausencia de los datos de
#   mortalidad (O_ij, E_ij) y los objetos de estructura espacial, que no se
#   distribuyen por razones de protección de datos (Conselleria de Sanitat /
#   INE). Véase el README del repositorio para más información.
# -----------------------------------------------------------------------------
# Dependencia externa:
#   pNimble — Martínez-Beneito, M.A. (Universitat de València)
#   https://github.com/MigueBeneito/pNimble
# =============================================================================


# ── Paquetes ──────────────────────────────────────────────────────────────────

# En caso necesario, instalar nimble desde el branch requerido:
# remotes::install_github("nimble-dev/nimble", ref = "transform_user_dist",
#                         subdir = "packages/nimble")

library(nimble)
library(nimbleHMC)
library(coda)
library(MCMCvis)
library(spdep)
library(RColorBrewer)


# ── pNimble ───────────────────────────────────────────────────────────────────
# Cargado directamente desde el repositorio original (no se distribuye aquí).
# Descomentar la línea local solo si se trabaja sin conexión a internet.

source("https://raw.githubusercontent.com/MigueBeneito/pNimble/main/RutinasNimble.0.2.R")
# source("./RutinasNimble.0.2.R")   # alternativa local
load.leroux()


# ── Datos [no disponibles públicamente] ───────────────────────────────────────
# Los objetos cargados aquí contienen conteos observados y esperados de
# mortalidad prematura por municipio, enfermedad y sexo. No se incluyen en
# el repositorio por motivos de protección de datos.

load("MortCV_102_2000_2023.RData", verbose = TRUE)
load("CartoMunis542.RData", verbose = TRUE)


# =============================================================================
# ESTRUCTURA ESPACIAL (CAR de Leroux)
# =============================================================================
W <- nb2mat(Veci, style="B")
D.W <- diag(apply(W,1,sum)) - W
valores.propios <- eigen(D.W)$values
aux <- apply(D.W, 1, function(x){which(x==-1)})
from.to <- as.matrix(rep(1:length(aux), sapply(aux,length)), ncol=1)
from.to <- cbind(from.to, unlist(aux))
from.to <- from.to[from.to[,1] < from.to[,2], ]

# =============================================================================
# ESPECIFICACIÓN DEL MODELO NIMBLE
# =============================================================================

ModeloCoregionalizacion_Original <- nimbleCode({
  
  for(j in 1:NDis){
    for(i in 1:N) {
      y[i,j] ~ dpois(E[i,j]*lambda[i,j])
      log(lambda[i,j]) <- alpha[j] + theta[i,j]
    }
    rho[j] ~ dunif(0, 1)
    alpha[j] ~ dflat()
    beta[1:N,j] ~ dcar_leroux(rho = rho[j], sd.theta = 1,  
                              Lambda = Lambda[1:N], 
                              from.to = from.to[1:n_vec, 1:2])
    
    zero.beta[j] ~ dnorm(mean(beta[1:N,j]), sd = 0.01)
    sd.theta[j] ~  dunif(0,10) 
  }
  
  # =============================================
  # ESTRUCTURA DE M: Triangular superior
  # =============================================
  
  # Diagonal: M[j,j] = 1
  for(j in 1:NDis){
    M[j,j] <- sd.theta[j] # NEW-IN
  }
  
  # Triángulo inferior: M[j,k] = 0 para j > k
  for(j in 2:NDis){
    for(k in 1:(j-1)){
      M[j,k] <- 0
    }
  }
  
  # Triángulo superior: M[j,k] con prior
  for(j in 1:(NDis-1)){
    for(k in (j+1):NDis){
      M[j,k] ~ dflat() # NEW-IN
    }
  }
  
  # ============================================
  # Cálculo de theta y Sigma 
  # ============================================
  
  theta[1:N, 1:NDis] <- beta[1:N,1:NDis] %*% M[1:NDis, 1:NDis]
  # Sigma[1:NDis, 1:NDis] <- t(M[1:NDis,1:NDis]) %*% M[1:NDis, 1:NDis]
})

# =============================================================================
# ORDEN DE LAS ENFERMEDADES
# =============================================================================
# Descomentando el bloque correspondiente se obtienen las distintas versiones
# del modelo. Solo uno de los bloques debe estar activo en cada ejecución.

cubo_esp2 <- apply(cubo_esp, c(1,3,4), sum)
cubo_esp_prem2 <- apply(cubo_esp_prem, c(1,3,4), sum)
cubo_obs2 <- apply(cubo_obs, c(1,3,4), sum)
cubo_obs_prem2 <- apply(cubo_obs_prem, c(1,3,4), sum)

whichHigh <- apply(cubo_obs2[,,1:102], c(1,3), sum)
whichHighM <- which(whichHigh[1,] > 542*2)
whichHighW <- which(whichHigh[2,] > 542*2)

whichHigh_prem <- apply(cubo_obs_prem2[,,1:102], c(1,3), sum)
whichHighM_prem <- which(whichHigh_prem[1,] > 542*2)
whichHighW_prem <- which(whichHigh_prem[2,] > 542*2)


## ── Ordenación original (referencia) ────────────────────────────────────────
y_full <- cbind(cubo_obs_prem2[1,,whichHighM_prem], 
                cubo_obs_prem2[2,,whichHighW_prem])
E_full <- cbind(cubo_esp_prem2[1,,whichHighM_prem], 
                cubo_esp_prem2[2,,whichHighW_prem])

names_causes <- colnames(y_full)
nombres_enf <- c(paste0("H_", names_causes[1:47]),
                 paste0("M_", names_causes[48:78]))

y_ordered       <- y_full
E_ordered       <- E_full
nombres_ordered <- nombres_enf

## ── Ordenación directa (menor a mayor casos totales) ────────────────────────
# total_casos    <- colSums(y_full)
# orden_directo  <- order(total_casos, decreasing = FALSE)
# y_ordered       <- y_full[, orden_directo]
# E_ordered       <- E_full[, orden_directo]
# nombres_ordered <- nombres_enf[orden_directo]

## ── Ordenación inversa (mayor a menor casos totales) ────────────────────────
# total_casos    <- colSums(y_full)
# orden_inverso  <- order(total_casos, decreasing = TRUE)
# y_ordered       <- y_full[, orden_inverso]
# E_ordered       <- E_full[, orden_inverso]
# nombres_ordered <- nombres_enf[orden_inverso]

## ── Ordenación aleatoria (repetir con set.seed(1) a set.seed(50)) ───────────
# set.seed(1)   # cambiar de 1 a 50 para las 50 permutaciones
# orden_aleatorio <- sample(1:78)
# y_ordered       <- y_full[, orden_aleatorio]
# E_ordered       <- E_full[, orden_aleatorio]
# nombres_ordered <- nombres_enf[orden_aleatorio]

# =============================================================================
# CONSTANTES, DATOS E INICIALES
# =============================================================================

constantes <- list(
  N = 542, 
  E = E_full,
  Lambda = valores.propios, 
  from.to = from.to, 
  n_vec = nrow(from.to), 
  NDis = 78
)

datos <- list(
  y = y_full,
  zero.beta = rep(0, 78) 
)
 
parametros <- c("M", "theta", "sd.theta", "alpha")

# Función para crear iniciales
iniciales <- function(){
  # M_init <- diag(78)
  M_init <- matrix(NA, 78, 78)
  for(j in 1:77){
    for(k in (j+1):78){
      M_init[j,k] <- rnorm(1, 0, 0.01)
    }
  }
  list(
    alpha = rnorm(78, 0, 0.05), 
    rho = runif(78), 
    beta = matrix(rnorm(542*78, 0, 1), ncol = 78), 
    M = M_init,
    sd.theta = runif(78, 0, 0.1)  
  )
}

# =============================================================================
# AJUSTE DEL MODELO
# =============================================================================
  
Coreg.MW.prem.Final_OrdenOriginal <- pNimble(
    code = ModeloCoregionalizacion_Original,  
    data = datos, 
    constants = constantes, 
    inits = iniciales, 
    seeds = 1:3, 
    niter = 5000,     
    nburnin = 2000,    
    thin = 3, 
    nchains = 3, 
    summary = TRUE, 
    WAIC = TRUE,  
    monitors = parametros, 
    ntfyAccount = "mariabarrera_nimble", 
    HMC = TRUE)


# =============================================================================
# GUARDAR RESULTADOS
# =============================================================================
# Ajustar la ruta de destino según el entorno de ejecución.

save(Coreg.MW.prem.Final_OrdenOriginal, names_causes, nombres_enf, # Original
     # orden_por_casos, nombres_ordered,  # Orden por casos
     # orden_aleatorio, nombres_ordered,   # Orden aleatorio
     file = "AA_Evaluacion_Leroux_Final/Coreg.MW.prem.Final_OrdenOriginal.Rdata")

summary_Coreg_reparam_Original <- 
  round(MCMCvis::MCMCsummary(Coreg.MW.prem.Final_OrdenOriginal$samples),3)

save(summary_Coreg_reparam_Original, 
     file = "summary_Coreg_reparam2a_Original.RData")

# =============================================================================
# DIAGNÓSTICO RÁPIDO (opcional)
# =============================================================================
# Descomentar para inspección visual de cadenas.

# MCMCvis::MCMCtrace(Coreg.MW.prem.OrdenOriginal$samples,
#                   pdf = FALSE, type = "trace")
# MCMCvis::MCMCtrace(Coreg.MW.prem.OrdenOriginal$samples,
#                   params = c("rho", "sd.theta", "alpha"),
#                   pdf = FALSE)
