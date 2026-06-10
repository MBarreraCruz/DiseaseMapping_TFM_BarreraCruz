# Mapas de riesgo sin artefactos: el M-modelo como solución bayesiana multivariante para la toma de decisiones en salud pública a gran escala. 

**Multivariate Bayesian disease mapping of 78 causes of premature mortality in 542 municipalities of the Comunitat Valenciana (2000–2023)**

Trabajo Final de Máster en Bioestadística — Universidad de Valencia, 2026  
Presentado por: **María Barrera Cruz**  
Dirigido por: Miguel Ángel Martínez Beneito (Dept. de Estadística e Investigación Operativa, UV)  
Licencia de código: [MIT]

---

## Descripción

Este repositorio contiene el código en R asociado al Trabajo de Fin de Máster (TFM) en Bioestadística titulado *"Mapas de riesgo sin artefactos: el M-modelo como solución bayesiana multivariante para la toma de decisiones en salud pública a gran escala"*. El objetivo del trabajo es analizar el riesgo de mortalidad prematura por 78 causas de enfermedad en los 542 municipios de la Comunitat Valenciana durante el período 2000–2023, mediante la comparación de tres familias de modelos bayesianos jerárquicos con estructura espacial:

1. **Modelos univariantes independientes** — 78 modelos CAR de Leroux ajustados de forma separada para cada enfermedad.
2. **Modelo de corregionalización** — 53 versiones (50 ordenaciones aleatorias + ordenación original + directa + inversa) que permiten estudiar el artefacto de dependencia de orden introducido por la parametrización triangular.
    - **Ordenación original** (referencia respecto al orden empleado en el M-modelo, según sexo y la propia ordenación de las causas de mortalidad del INE)
    - **Ordenación directa** (ordenación de mayor a menor número de casos observados por enfermedad)
    - **Ordenación inversa** (ordenación de menor a mayor número de casos observados por enfermedad)
3. **M-modelo** — modelo multivariante conjunto (Botella-Rocamora et al., 2015) que modela simultáneamente la estructura espacial y la dependencia entre enfermedades mediante una descomposición matricial invariante al orden.

La principal contribución del trabajo es demostrar, a partir de datos reales, las ventajas del enfoque multivariante. Asimismo, se pone de manifiesto que la parametrización basada en el modelo de corregionalización induce una dependencia respecto al orden de las enfermedades en el vector de respuesta, mientras que el M-modelo es invariante a dicho orden y presenta una mayor eficiencia computacional. Además de la comparación metodológica entre modelos, el trabajo demuestra la viabilidad computacional del análisis conjunto de 78 causas de mortalidad en 542 municipios durante un periodo de 23 años, constituyendo, hasta donde alcanza nuestro conocimiento, uno de los estudios de riesgo espacial multivariante de mayor dimensión realizados con datos reales. Los resultados muestran que, gracias al uso de herramientas modernas de programación probabilística, este tipo de análisis es actualmente factible y puede incorporarse de forma sistemática a tareas de vigilancia epidemiológica y apoyo a la toma de decisiones en salud pública.

## Contenido del repositorio

Este repositorio consta de los siguientes directorios:

1. **Datos**. Los datos de mortalidad proceden del Registro de Mortalidad de la Comunitat Valenciana (Conselleria de Sanitat) y del INE. Por motivos de confidencialidad, **los microdatos no se incluyen en este repositorio**. La carpeta `data/` contiene únicamente ... DEFINIR.
2. **Modelos**. Esta carpeta contiene los

    - **Análisis comparativos incluidos:**
      - **WAIC**: comparación de bondad de ajuste predictivo entre los tres enfoques. El WAIC del M-modelo se contrasta con la suma de los 78 WAICs univariantes y con los WAICs de todas las versiones del modelo de corregionalización.
      - **Matrices de correlación**: entre efectos espaciales estimados por cada modelo, para los dos sexos.
      - **Inflación de varianza**: correlación de Spearman entre la posición de cada enfermedad en el vector de respuesta y el ratio de varianza espacial estimada respecto al modelo de referencia.
      - **Mapas de riesgo**: mapas coropléticos para enfermedades seleccionadas, comparando los tres modelos.
      - **Test de Welch**: comparación de correlaciones cruce-sexo vs. correlaciones intra-sexo.

--

## Estructura del repositorio

```
DiseaseMapping_TFM_BarreraCruz/
│
├── README.md
├── LICENSE
├── .gitignore
│
├── data/                        
│   ├── RDatas/                  # 
│   ├── processed/               # 

│
├── models/                      # Scripts de especificación y ajuste de modelos
│   ├── univariate/
│   │   ├── fit_univariate_all.R          # Bucle principal: ajusta los 78 modelos univariantes
│   │   └── model_leroux_univariate.R     # Especificación NIMBLE del CAR de Leroux univariante
│   │
│   ├── mmodel/
│   │   ├── fit_mmodel.R                  # Ajuste del M-modelo conjunto
│   │   └── model_mmodel_nimble.R         # Especificación NIMBLE del M-modelo
│   │
│   └── coregionalization/
│       ├── fit_coreg_single.R            # Ajuste de una ordenación concreta
│       ├── fit_coreg_permutations.R      # Bucle: 50 ordenaciones aleatorias + directa + inversa
│       └── model_coreg_nimble.R          # Especificación NIMBLE del modelo de corregionalización
│
├── analysis/                    # Scripts de análisis comparativo
│   ├── waic/
│   │   └── compare_waic.R               # Cálculo y comparación de WAIC entre los tres enfoques
│   ├── correlation/
│   │   └── correlation_matrices.R       # Matrices de correlación entre efectos espaciales
│   ├── variance_inflation/
│   │   └── variance_inflation_coreg.R   # Análisis de inflación de varianza por posición (Spearman)
│   └── maps/
│       └── risk_maps.R                  # Mapas coropléticos de riesgo relativo (carto_munis)
│
├── results/                     # Salidas generadas (no versionadas, generadas localmente)
│   ├── figures/                 # Figuras (.pdf, .png)
│   ├── tables/                  # Tablas (.csv, .tex)
│   └── diagnostics/             # Diagnósticos MCMC (Rhat, trazas, etc.)
│
│   [Dependencia externa: pNimble — https://github.com/MigueBeneito/pNimble/blob/main/RutinasNimble.0.2.R
│   [Se carga vía source() desde cada script de modelo; no se incluye en este repositorio]
```

--

## Dependencias:

### pNimble (Martínez-Beneito, M.A.)

Los scripts de ajuste de modelos dependen de **pNimble**, una librería auxiliar que implementa `pNimble()` (wrapper paralelo para NIMBLE/HMC) y `load.leroux()` (distribución CAR de Leroux como función NIMBLE personalizada). El código se encuentra disponible en su repositorio público:

> Repositorio original: <https://github.com/MigueBeneito/pNimble>  
> Autor: Martínez-Beneito, M.A. (Universitat de València / FISABIO)
