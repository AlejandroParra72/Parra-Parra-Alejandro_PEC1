

################### IMPORTACIÓN Y DEPURACIÓN DE LOS DATOS ################### 

# Comenzamos importando los datos metabólicos y los metadatos procedentes del estudio MTBLS2452 a R:

library(readr)
library(dplyr)

# Importación de metadatos: 
metadatos <- read.table("C:/Users/aleja/Downloads/s_MTBLS2452.txt", header = TRUE, sep = "\t") # Importamos los 
# metadoatos del archivo txt

metadatos <- metadatos %>%
  select(Factor.Value.Treatment., everything()) # Reorganizamos las columnas del dataframe "metadatos", dejando
# la columna que se corresponde con el tipo de muestra como la primera

metadatos$Factor.Value.Treatment. <- as.factor(metadatos$Factor.Value.Treatment.) # Transformamos los datos de 
# la columna Factor.Value.Treatment. en factores.

str(metadatos) # Comprobamos que la variable Factor.Value.Treatment. es un factor

colnames(metadatos)[colnames(metadatos) == "Factor.Value.Treatment."] <- "groups" # Cambiamos el nombre de la
# variable Factor.Value.Treatment. a "groups", lo cual facilitará algunos aspectos del análisis posterior. 

metadatos <- metadatos[ , -c(17:18)] # Eliminamos las últimas dos columnas del dataframe "metadatos", ya que 
# solo presentan valores NA. 

rownames(metadatos) <- metadatos$Sample.Name #Asignamos los identificadores de las muestras "FREI-XXX" a cada fila.

# Importación de datos metabólicos: 
datos_metabolicos <- read.table("C:/Users/aleja/Downloads/m_MTBLS2452_metabolite_profiling_mass_spectrometry_v2_maf.tsv", 
                                header = TRUE, sep = "\t", quote = "", comment.char = "", stringsAsFactors = FALSE, 
                                fill = TRUE, check.names = FALSE)

# A continucaión, depuramos este dataframe para favorecer el posterior análisis: 
niveles_metabolitos <- datos_metabolicos[, grep("FREI", colnames(datos_metabolicos))] # Extraemos las columnas que nos interesan 

rownames(niveles_metabolitos) <- datos_metabolicos$metabolite_identification # Asignamos los metabolitos a cada fila.

niveles_metabolitos <- niveles_metabolitos[-1, ] # Eliminamos la primera fila, que se correponde con otro ID de las muestras
# y no es información sobre los metabolitos



################### SummarizedExperiment ################### 

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("SummarizedExperiment")
library(SummarizedExperiment)

# A continuación ya podemos crear el objeto de clase SummarizedExperiment que contenga los datos y los metadatos: 
Objeto_SE <- SummarizedExperiment(assays = list(counts = as.matrix(niveles_metabolitos)), colData = metadatos)  

# Tal y como menciona el tutorial propuesto en el enunciado: "SummarizedExperiment es muy similar  a ExpressionSet, 
# pero la principal diferencia es que SummarizedExperiment es más flexible en su información de filas, permitiendo 
# tanto GRanges como aquellos descritos por DataFrames arbitrarios. Esto lo hace ideal para una variedad de 
# experimentos, particularmente experimentos basados en secuenciación como RNA-Seq y ChIp-Seq."



################### ANÁLISIS EXPLORATORIO INICIAL ################### 

# Mostramos la cabecera del objeto: 
head(assay(Objeto_SE))

# ¿Cuántos metabolitos y muestras hay en el objeto?
dim(assay(Objeto_SE)) 
# Podemos identificar 662 metabolitos (filas) y 12 columnas (muestras).

# Mostramos cuales son las muestras del estudio (por motivos de extensión no mostramos los metabolitos
# puesto que son 662): 
colnames(assay(Objeto_SE))

# Mostramos los metadatos
colData(Objeto_SE)

# Identificamos cuantos valores NA hay en el objeto:
sum(is.na(assay(Objeto_SE)))
# Hay 939 valores faltantes en el objeto. 



################### PREPROCESADO DE LOS DATOS ################### 

# 1) IMPUTACIÓN DE VALORES NA: 

# Instalamos POMA de Bioconductor, para poder usar la función PomaImpute. Esta función nos permite 
# imputar los valores NA, ya que como hemos visto, son 939 valores N. 

BiocManager::install("POMA")
library(POMA)
BiocManager::install("ggtext")
library(ggtext)

rownames(Objeto_SE) <- make.names(rownames(Objeto_SE), unique = TRUE)
Objeto_SE_NA <- Objeto_SE %>%
  PomaImpute(method = "knn", zeros_as_na = TRUE, remove_na = TRUE, cutoff = 20)

dim(assay(Objeto_SE_NA)) # Podemos ver que ahora pasamos a tener 598 metabolitos (se han eliminado
# 64 de ellos al imputar los valores NA). 

# 2) NORMALIZACIÓN

# Llevamos a cabo una normalización de los datos, de esta forma evitamos el sesgo de valores muy grandes
# o muy pequeños de las intensidades de señal producidas por MS:  

Objeto_SE_norm <- Objeto_SE_NA %>% 
  PomaNorm(method = "log_pareto")

# A continuación. mostramos como varía la intensidad de señal por muestra antes y despues de la 
# normalización:
PomaBoxplots(Objeto_SE_NA) # Antes de la normalización
PomaBoxplots(Objeto_SE_norm) # Después de la normalización

# Procedemos de la misma forma con la densidad: 
library(ggplot2)
PomaDensity(Objeto_SE_NA, x = "features") + 
  theme(legend.position = "none")

PomaDensity(Objeto_SE_norm, x = "features") + 
  theme(legend.position = "none")

# 3) Detección de valores atípicos

PomaOutliers(Objeto_SE_norm)$polygon_plot

# Podemos ver que las muestras se separan entre controles y tratadas perfectamente. .  

PomaOutliers(Objeto_SE_norm, labels = TRUE)$polygon_plot # En este gráfico, podemos identificar 
# que el valor atípico es el FREI-01737. 

# Realizamos el análisis de outliers:
Objeto_SE_Out<- PomaOutliers(Objeto_SE_norm,  method = "euclidean",
                             type = "median",
                             coef = 1.5,
                             labels = FALSE)
Objeto_SE_Out$data # Podemos ver que detecta la muestra FREI-01741 como outlier, sin embargo, en la 
# gráfica anterior vemos que se separa de los controles aunque esté situado muy abajo en la gráfica. 
# Por lo tanto, se ha decidido mantenerlo en el estudio y prestar atención a como se comporta en 
# el posterior análisis exploratorio. 


################### ANÁLISIS EXPLORATORIO ###################

# Hasta ete punto, hemos realizado una descripción de los datos y los hemos transformado
# (eliminación de valores NA y normalización). A continuación, se muestra el código necesario
# para el análisis exploratorio. Este análisis consta de dos partes: 
# 1) Análisis de componentes principales (PCA)
# 2) Análisis de patrones y clustering (HeatMap y Clustering)

# 1) PCA

# Utilizamos la función PomaPCA del paque POMA para realiza el análisis de componentes
# principales sobre nuestro Objeto_SE normalizado. 

# Esta función nos permite llevar a cabo una reducción de la dimensiones
PomaPCA(Objeto_SE_norm, 
        outcome = "groups", # Escogemos que las muestras se coloreen en función de 
        # si están tratadas o son controles
        center = TRUE, # Cada variable será centrada restándole su media
        scale = FALSE, # No escalamos los datos ya que las variables tienen unidades similares
        ncomp = 10, # Número de componentes principales
        labels = FALSE, 
        ellipse = TRUE, # Añadimos una elipse que representa el área en la que se espera que caigan
        # el 95% de los puntos si los datos siguen una distribución normal. Como hemos hecho
        # una normalizacion antes, es normal que todos los datos caigan dentro de la elise
        load_length = NA) # Eliminamos del gráfico los metabolitos

# Se ha ajustado inicialmente ncomp = 10 para poder visualizar el "scree plot" y aplicar el método
# del codo. Se ha podido observar que a partir del PCA3 ya se supera más del 80% de la varianza explicada, 
# por lo que con un ncomp = 3 es suficiente. 

# De este modo, el PCA quedaría de la siguiente forma: 
PomaPCA(Objeto_SE_norm, 
        outcome = "groups", # Escogemos que las muestras se coloreen en función de 
        # si están tratadas o son controles
        center = TRUE, # Cada variable será centrada restándole su media
        scale = FALSE, # No escalamos los datos ya que las variables tienen unidades similares
        ncomp = 3, # Número de componentes principales
        labels = FALSE, 
        ellipse = TRUE, # Añadimos una elipse que representa el área en la que se espera que caigan
        # el 95% de los puntos si los datos siguen una distribución normal. Como hemos hecho
        # una normalizacion antes, es normal que todos los datos caigan dentro de la elise
        load_length = NA) # Eliminamos del gráfico los metabolitos

# 2.1) Heatmap

# Heatmap de la intensidad de señal de cada metabolito (Método de búsqueda de patrones)
PomaHeatmap(Objeto_SE_norm,
            covs = "groups",
            sample_names = TRUE,
            feature_names = FALSE,
            show_legend = TRUE)

# Podemos identficar que los metabolitos con menor intensidad en controles, son los de mayor 
# intensidad en las muestras tratadas, y viceversa. Por lo tanto, observamos dos patrones
# diferentes de abundancia de metabolitos. Cada uno de estos patrones, se corresponde con 
# un grupo de pacientes (controles y tratados)

# 2.1) Clustering

PomaClust(
  Objeto_SE_norm,
  method = "euclidean",
  k = NA,
  k_max = floor(min(dim(data))/2),
  show_clusters = TRUE,
  labels = TRUE)

# Podemos identficar que las muestras se clusterizan perfectamente en las tratadas y las controles


################### ANÁLISIS BIOLÓGICO ###################

# Finalmente, representamos la abundancia relativa (normalizada) de algunos metabolitos
# concretos interesantes para la interpretación biológica

piruvato <- assay(Objeto_SE_norm)["pyruvate",]
lactate <- assay(Objeto_SE_norm)["lactate",]
NAD_mas <- assay(Objeto_SE_norm)["nicotinamide_adenine_dinucleotide", ]
acetyl_CoA <- assay(Objeto_SE_norm)["acetyl_CoA", ]

boxplot(piruvato ~ colData(Objeto_SE_norm)$groups, 
        main = "Abundancia relativa del piruvato", 
        ylab = "Abundancia", 
        xlab = "Grupo de Tratamiento",
        col = c("slateblue1", "goldenrod1")) 

boxplot(lactate ~ colData(Objeto_SE_norm)$groups, 
        main = "Abundancia relativa del lactato", 
        ylab = "Abundancia", 
        xlab = "Grupo de Tratamiento",
        col = c("slateblue1", "goldenrod1"))

boxplot(NAD_mas ~ colData(Objeto_SE_norm)$groups, 
        main = "Abundancia relativa del NAD+", 
        ylab = "Abundancia", 
        xlab = "Grupo de Tratamiento",
        col = c("slateblue1", "goldenrod1"))

boxplot(acetyl_CoA ~ colData(Objeto_SE_norm)$groups, 
        main = "Abundancia relativa del acetil CoA", 
        ylab = "Abundancia", 
        xlab = "Grupo de Tratamiento",
        col = c("slateblue1", "goldenrod1")) 


