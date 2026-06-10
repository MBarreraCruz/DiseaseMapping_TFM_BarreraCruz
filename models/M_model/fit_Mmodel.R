# =============================================================================
# Mmodel.R
# M-modelo multivariante conjunto — 78 causas de mortalidad prematura
# -----------------------------------------------------------------------------
# TFM en Bioestadística — María Barrera Cruz
# Universidad de Valencia, Facultad de Matemáticas, 2026
# Tutor: Miguel Ángel Martínez-Beneito (Dept. Estadística e Investigación Operativa, UV)
# -----------------------------------------------------------------------------
# Descripción:
#   Ajuste del M-modelo bayesiano multivariante (Botella-Rocamora et al., 2015)
#   para 78 causas de mortalidad prematura en los 542 municipios de la
#   Comunitat Valenciana (2000–2023), considerando conjuntamente hombres
#   (47 causas) y mujeres (31 causas).

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
# Los objetos cargados contienen conteos observados y esperados de mortalidad
# prematura por municipio, enfermedad y sexo. No se incluyen en el repositorio
# por motivos de protección de datos.

load("MortCV_102_2000_2023.RData", verbose = TRUE)
load("CartoMunis542.RData", verbose = TRUE)


# =============================================================================
# PREPARACIÓN DE DATOS
# =============================================================================

cubo_esp_prem2 <- apply(cubo_esp_prem, c(1, 3, 4), sum)
cubo_obs_prem2 <- apply(cubo_obs_prem, c(1, 3, 4), sum)

# Selección de enfermedades con suficientes casos (> 2 por municipio en media)
whichHigh_prem  <- apply(cubo_obs_prem2[,, 1:102], c(1, 3), sum)
whichHighM_prem <- which(whichHigh_prem[1, ] > 542 * 2)   # 47 causas en hombres
whichHighW_prem <- which(whichHigh_prem[2, ] > 542 * 2)   # 31 causas en mujeres

# Matrices de observados y esperados: 542 municipios × 78 enfermedades
# (columnas 1:47 = hombres, 48:78 = mujeres)
y_prem <- cbind(cubo_obs_prem2[1, , whichHighM_prem],
                cubo_obs_prem2[2, , whichHighW_prem])
E_prem <- cbind(cubo_esp_prem2[1, , whichHighM_prem],
                cubo_esp_prem2[2, , whichHighW_prem])


# =============================================================================
# ESTRUCTURA ESPACIAL (CAR de Leroux)
# =============================================================================
# Construida a partir de la lista de vecindad `Veci`. 

W               <- nb2mat(Veci, style = "B")
D.W             <- diag(apply(W, 1, sum)) - W
valores.propios <- eigen(D.W)$values
aux             <- apply(D.W, 1, function(x) which(x == -1))
from.to         <- as.matrix(rep(1:length(aux), sapply(aux, length)), ncol = 1)
from.to         <- cbind(from.to, unlist(aux))
from.to         <- from.to[from.to[, 1] < from.to[, 2], ]


# =============================================================================
# ESPECIFICACIÓN DEL MODELO NIMBLE
# =============================================================================

Mmodelo <- nimbleCode({
  for (j in 1:NDis) {
    for (i in 1:N) {
      y[i, j]           ~ dpois(E[i, j] * lambda[i, j])
      log(lambda[i, j]) <- alpha[j] + theta[i, j]
    }
    rho[j]   ~ dunif(0, 1)
    alpha[j] ~ dflat()
    beta[1:N, j] ~ dcar_leroux(
      rho      = rho[j],
      sd.theta = 1,
      Lambda   = Lambda[1:N],
      from.to  = from.to[1:n_vec, 1:2]
    )
    
    zero.beta[j] ~ dnorm(mean(beta[1:N, j]), sd = 0.01)
    
    for (k in 1:NDis) { M[j, k] ~ dnorm(0, sd = sigma.m[k]) }
    sigma.m[j] ~ dhalfflat()
  }
  
  theta[1:N, 1:NDis] <- beta[1:N, 1:NDis] %*% M[1:NDis, 1:NDis]
  
  Sigma[1:NDis, 1:NDis] <- t(M[1:NDis, 1:NDis]) %*% M[1:NDis, 1:NDis]
})


# =============================================================================
# CONSTANTES, DATOS E INICIALES
# =============================================================================

constantes <- list(
  N       = 542,
  E       = E_prem,
  Lambda  = valores.propios,
  from.to = from.to,
  n_vec   = nrow(from.to),
  NDis    = 78
)

datos <- list(
  y         = y_prem,
  zero.beta = rep(0, 78)
)

parametros <- c("alpha", "theta", "Sigma", "sigma.m")

iniciales <- function() {
  list(
    alpha   = rnorm(78, 0, 0.05),
    rho     = runif(78),
    beta    = matrix(rnorm(542 * 78, 0, 1), ncol = 78),
    M       = matrix(rnorm(78 * 78), ncol = 78),
    sigma.m = runif(78)
  )
}


# =============================================================================
# AJUSTE DEL MODELO
# =============================================================================

Mmodel.MW.prem <- pNimble(
  code      = Mmodelo,
  data      = datos,
  constants = constantes,
  inits     = iniciales,
  seeds     = 1:3,
  niter     = 5000,
  nburnin   = 2000,
  thin      = 3,
  nchains   = 3,
  summary   = TRUE,
  WAIC      = TRUE,
  monitors  = parametros,
  HMC       = TRUE
)


# =============================================================================
# GUARDAR RESULTADOS
# =============================================================================
# Ajustar la ruta de destino según el entorno de ejecución.

save(Mmodel.MW.prem, file = "Mmodel.MW.prem.RData")

summary_Mmodel.MW.prem <-
  round(MCMCvis::MCMCsummary(Mmodel.MW.prem$samples), 3)

save(summary_Mmodel.MW.prem, file = "summary_Mmodel.MW.prem.RData")