# =============================================================================
# Univariante_Leroux.R
# Modelos univariantes CAR de Leroux — 78 causas de mortalidad prematura
# -----------------------------------------------------------------------------
# TFM en Bioestadística — María Barrera Cruz
# Universidad de Valencia, Facultad de Matemáticas, 2026
# Tutor: Miguel Ángel Martínez-Beneito (Dept. Estadística e Investigación Operativa, UV)
# -----------------------------------------------------------------------------
# Descripción:
#   Ajuste de 78 modelos bayesianos univariantes independientes, uno por cada
#   causa de mortalidad prematura en los 542 municipios de la Comunitat
#   Valenciana (2000–2023). 
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

library(nimble)
library(nimbleHMC)
library(coda)
library(MCMCvis)
library(sf)
library(spdep)


# ── pNimble ───────────────────────────────────────────────────────────────────
# Cargado directamente desde el repositorio original (no se distribuye aquí).
# Descomentar la línea local solo si se trabaja sin conexión a internet.

source("https://raw.githubusercontent.com/MigueBeneito/pNimble/main/RutinasNimble.0.2.R")
# source("./RutinasNimble.0.2.R")   # alternativa local
load.leroux()


# ── Datos [no disponibles públicamente] ───────────────────────────────────────
# Los objetos cargados contienen conteos observados y esperados de mortalidad
# prematura por municipio, enfermedad y sexo. No se incluyen en el repositorio
# por motivos de protección de datos.

load("MortCV_102_2000_2023.RData", verbose = TRUE)
load("CartoMunis542.RData", verbose = TRUE)


# =============================================================================
# PREPARACIÓN DE DATOS
# =============================================================================

cubo_esp2      <- apply(cubo_esp,      c(1, 3, 4), sum)
cubo_esp_prem2 <- apply(cubo_esp_prem, c(1, 3, 4), sum)
cubo_obs2      <- apply(cubo_obs,      c(1, 3, 4), sum)
cubo_obs_prem2 <- apply(cubo_obs_prem, c(1, 3, 4), sum)

# Selección de enfermedades con suficientes casos (> 2 por municipio en media)
whichHigh      <- apply(cubo_obs2[,, 1:102],      c(1, 3), sum)
whichHighM     <- which(whichHigh[1, ] > 542 * 2)
whichHighW     <- which(whichHigh[2, ] > 542 * 2)
whichHigh_prem <- apply(cubo_obs_prem2[,, 1:102], c(1, 3), sum)
whichHighM_prem <- which(whichHigh_prem[1, ] > 542 * 2)
whichHighW_prem <- which(whichHigh_prem[2, ] > 542 * 2)

# Matrices de observados y esperados: 542 municipios × 78 enfermedades
# (columnas 1:47 = hombres, 48:78 = mujeres)
O <- cbind(cubo_obs_prem2[1, , whichHighM_prem],
           cubo_obs_prem2[2, , whichHighW_prem])
E <- cbind(cubo_esp_prem2[1, , whichHighM_prem],
           cubo_esp_prem2[2, , whichHighW_prem])


# =============================================================================
# ESTRUCTURA ESPACIAL (CAR de Leroux)
# =============================================================================
# Construida a partir de la lista de vecindad `Veci`.

Q       <- diag(sapply(Veci, length)) - nb2mat(Veci, style = "B")
aux     <- apply(Q, 1, function(x) which(x == -1))
from.to <- as.matrix(rep(1:length(aux), sapply(aux, length)), ncol = 1)
from.to <- cbind(from.to, unlist(aux))
from.to <- from.to[from.to[, 1] < from.to[, 2], ]
eig.Q   <- eigen(Q)


# =============================================================================
# ESPECIFICACIÓN DEL MODELO NIMBLE
# =============================================================================
# Modelo univariante CAR de Leroux para una única enfermedad j:

Leroux <- nimbleCode({
  for (i in 1:N) {
    O[i]          ~ dpois(E[i] * R[i])
    log(R[i])     <- alfa + sigma.tot * theta[i]
    theta.tot[i]  <- sigma.tot * theta[i]
  }
  
  theta[1:N] ~ dcar_leroux(
    rho      = rho,
    sd.theta = 1,
    Lambda   = Lambda[1:N],
    from.to  = from.to[1:n_vec, 1:2]
  )
  
  # Restricción suma-cero (estabilidad del intercepto)
  zero.theta ~ dnorm(mean(theta[1:N]), sd = 0.001)
  
  # Distribuciones a priori
  rho       ~ dunif(0, 1)
  sigma.tot ~ dunif(0, 2)
  alfa      ~ dflat()
})


# =============================================================================
# AJUSTE: BUCLE SOBRE LAS 78 ENFERMEDADES
# =============================================================================
# Cada modelo se ajusta, se guarda en disco y se elimina de memoria antes de
# pasar a la siguiente enfermedad. El vector `waics.uni` acumula el WAIC de
# cada modelo para la comparación posterior.

iniciales.uni <- function() {
  list(
    alfa      = rnorm(1, 0, 0.1),
    rho       = runif(1, 0.3, 0.7),
    sigma.tot = runif(1, 0.1, 1),
    theta     = rnorm(542, 0, 1)
  )
}

parametros.uni <- c("alfa", "rho", "sigma.tot", "theta.tot")

waics.uni <- rep(NA_real_, 78)

for (i in 1:78) {
  
  constantes.uni <- list(
    N       = 542,
    E       = E[, i],
    Lambda  = eig.Q$values,
    from.to = from.to,
    n_vec   = nrow(from.to)
  )
  
  datos.uni <- list(
    O          = O[, i],
    zero.theta = 0
  )
  
  Resul.UniModel <- pNimble(
    code      = Leroux,
    data      = datos.uni,
    constants = constantes.uni,
    inits     = iniciales.uni,
    seeds     = 1:3,
    niter     = 5000,
    nburnin   = 2000,
    thin      = 3,
    nchains   = 3,
    summary   = TRUE,
    WAIC      = TRUE,
    monitors  = parametros.uni,
    HMC       = TRUE
  )
  
  waics.uni[i] <- Resul.UniModel$WAIC$WAIC
  
  # Ajustar la ruta de destino según el entorno de ejecución
  save(Resul.UniModel,
       file = paste0("Resul.UniModel_", i, ".RData"))
  
  rm(Resul.UniModel)
}

# WAIC acumulado de los 78 modelos (usado en la comparación global)
save(waics.uni, file = "waics_uni.RData")