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
├── docs/                                # HTML renderizado — visualización directa
|   └── index.html            # Resultados con figuras (GitHub Pages)
|
├── models/                              # Scripts de ajuste de modelos (NIMBLE/HMC)
│   ├── Univariate/
│   │   └── fit_Univariate.R             # Bucle sobre las 78 enfermedades (CAR de Leroux)
│   ├── M_model/
│   │   └── fit_Mmodel.R                 # M-modelo multivariante conjunto
│   └── Coregionalization/
│       └── fit_Coregionalization.R      # Modelo de corregionalización
│                                        # Descomentando el bloque correspondiente
│                                        # se ejecutan las distintas ordenaciones
│                                        # (original / directa / inversa / aleatoria)
│
└── analysis/                            # Análisis comparativo y generación de figuras
    ├── Analysis_Results.Rmd             # Documento reproducible — ejecutar desde aquí
    │
    ├── data/                            # Objetos precalculados (sin datos brutos de mortalidad)
    │   ├── RData/                       # Resúmenes de las distribuciones posteriores
    │   │   ├── Theta_mean_UNI.RData         # Medias posteriores θ — univariante (542 × 78)
    │   │   ├── theta_multi.RData            # Medias posteriores θ — M-modelo (542 × 78)
    │   │   ├── theta_coreg_original.RData   # Medias posteriores θ — coreg. original (542 × 78)
    │   │   ├── prob_pos_matrix.RData        # P(θ_ij > 0) — M-modelo (542 × 78)
    │   │   ├── prob_pos_matrix_coreg.RData  # P(θ_ij > 0) — coreg. original (542 × 78)
    │   │   ├── waics_uni.RData              # Vector de 78 WAICs univariantes
    │   │   ├── waic_mmodel.RData            # WAIC escalar del M-modelo
    │   │   ├── waic_coreg_Original.RData    # WAIC — coreg. ordenación original
    │   │   ├── waic_coreg_directa_inversa.RData  # WAICs ordenaciones directa e inversa
    │   │   ├── waic_coreg_vec.RData         # Vector de 50 WAICs (ordenaciones aleatorias)
    │   │   ├── var_theta_coreg.RData        # Lista de varianzas de θ por permutación
    │   │   ├── orden_aleatorio_list.RData   # Lista de 50 vectores de permutación
    │   │   ├── df_cor.RData                 # Correlaciones entre enfermedades (formato largo)
    │   │   └── variabilidad_por_enfermedad.RData  # Rango de inflación de varianza por enfermedad
    │   └── spatial/
    │       └── CartoMunis542.RData          # Cartografía de los 542 municipios (IGN/INE)
    │
    └── results/
        └── output/                      # Figuras generadas por Analysis_Results.Rmd
            ├── boxplot_waic.png
            ├── correlacion_plot.png
            ├── varianza_plot.png
            ├── varianza_enfermedad.png
            ├── I_Moran.png
            ├── mapa_H012.png  ··· mapa_M090.png   # θ estimado (3 modelos, escala común)
            └── probH012.png   ··· probM090.png     # P(θ_ij > 0) (coreg. vs M-modelo)--
```


## Dependencias:

### pNimble (Martínez-Beneito, M.A.)

Los scripts de ajuste de modelos dependen de **pNimble**, una librería auxiliar que implementa `pNimble()` (wrapper paralelo para NIMBLE/HMC) y `load.leroux()` (distribución CAR de Leroux como función NIMBLE personalizada). El código se encuentra disponible en su repositorio público:

> Repositorio original: <https://github.com/MigueBeneito/pNimble>  
> Autor: Martínez-Beneito, M.A. (Universitat de València)
