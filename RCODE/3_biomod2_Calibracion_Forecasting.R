#############################################################
## Modelos de prediccion de habitat en ciencias forestales ##
#############################################################

###############################################
# Ejemplo de un modelo de hábitat con biomod2 #
###############################################

## --------------------------------------------------------- ##
## --------------------------------------------------------- ##

## PUNTOS IMPORTANTES ##
########################

# Normalmente, `biomod2´ se suele configurar en varios pasos en base a diferentes funciones que contiene el paquete:
# 1.- BIOMOD_FormatingData -> Datos de entrada de la(s) especie(s)
# 2.- BIOMOD_ModelingOptions -> Parametrizacion de los algoritmos
# 3.- BIOMOD_Modeling -> Calibracion de modelos individuales
# 4.- Evaluacion de los modelos parciales
# 5.- BIOMOD_EnsembleModeling -> Obtencion de modelos de conjunto (ensemble models)
# 6.- Evaluacion de los modelos de conjunto
# 7.- BIOMOD_Projection -> Proyeccion de modelos individuales
# 8.- BIOMOD_EnsembleForecasting -> Proyeccion de modelos de conjunto

## INPUTS NECESARIOS ##
#######################

# Para ejecutar el framework de modelado de conjuntos `biomod2´ se necesitan las siguientes entradas:
# 
# (i) una tabla *csv con ubicaciones de presencia conocidas (x/y o lat/lon) para la(s) especie(s) objetivo de interes
#
# (ii) Un conjunto de variables predictoras que se sabe que influyen o restringen la idoneidad ambiental de la(s) especie(s) objetivo. Estos datos deben estar basados en raster y representar las condiciones actuales temporalmente coherentes con las ubicaciones de presencia. Se esperan archivos separados para cada variable (banda unica). ¡No utilizar archivos multibanda!
# Los formatos de archivo permitidos incluyen ESRI ASCII, GeoTIFF, ... (cualquiera compatible con GDAL)
#
# 

## --------------------------------------------------------- ##
## --------------------------------------------------------- ##

# 1.- Cargar librerias basicas en el espacio de trabajo
#######################################################
library(raster)
library(magrittr) # Un Forward-Pipe Operator para R
library(dismo)
library(rgeos) # Operaciones geometricas con informacion geografica
library(biomod2) # Modelos de distribucion de especies y ensamblado
library(ggplot2)
library(RColorBrewer) # Crea bonitas paletas de colores especialmente para mapas tematicos
library(dplyr) # Para manipular datos y trabajar con objetos, tanto en memoria como fuera de ella
library(sp) # Clases y métodos para datos espaciales
library(rasterVis) # Metodos de visualizacion para datos raster 
# 2.- Definir el directorio de trabajo
######################################
setwd("C:/GitHub_repo")

# 3.- Rutas para importar/exportar archivos
###########################################
path2sp<-"./DATA/ResponseVariable/filtering_records/ pb_thinned/thinned_data_thin3.csv" # Path para las variables bioclimaticas
path2layers<- "./DATA/PredictorVariables/2_test_multicolinealidad_seleccion_variables/2_3_Variables_definitivas" # Path para el raster de human footprint
path2work<- "./OUTPUTS/" # Path para la salida de resultados

# 4.- Cargar los datos de la(s) especie(s) y las variables predictoras
##########################################################

# NOTA: Estos parametros deben ser configurados antes de computar los modelos

# Cargar los datos de la especie
DataSpecies <- read.csv(path2sp,sep=',')

# Explorar los datos
head(DataSpecies, n = 5)

# Plot
qplot(x=lon,y=lat,data=DataSpecies)

# Convertir los datos en un objeto espacial
spPointData <- SpatialPoints(DataSpecies[,c("lon","lat")], 
                             proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# Nombre de la(s) especie(s) de estudio
myRespName <- 'Pb'

# Datos de presencia/ausencia para nuestra especie 
myResp <- spPointData
myResp

# Cargar las variables ambientales (predictores)
myExpl<-stack(list.files(path=path2layers,pattern='asc',full.names=TRUE))
crs(myExpl) <- CRS("+init=epsg:4326")
plot(myExpl)
myExpl

# 5.- Formato de los datos
##########################

# Configurar `biomod2´ para la calibracion
# Numero de conjuntos de PA (pseudoausencias) = 10. Para este ejemplo se configura en 1 conjunto de pseudo-absences, sin embargo, se recomienda que sean más (e.g., 3,5,10...).
# Numero de PA por conjunto = Numero de presencias (RECOMENDABLE). Se recomienda no menos de 30.
# Estrategia de seleccion de PA = 'random'. Ver ayuda para otros metodos.

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp, # Datos de la(s) especie(s)
                                     expl.var = myExpl, # Predictores
                                     resp.name = myRespName, # Nombre de la(s) especie(s)
                                     PA.nb.rep = 1, # Replicas de pseudoausencias - Se recomienda minimo 10
                                     PA.nb.absences = length(myResp),
                                     myResp, # Nº de pseudoausencias - Se recomienda = presencias
                                     PA.strategy = 'random', # Ver ayuda para otros metodos
                                     PA.dist.min = 0,
                                     PA.dist.max = NULL,
                                     PA.sre.quant = 0.025,
                                     PA.table = NULL,
                                     na.rm = TRUE) # Eliminar posibles NAs

# Print formato de los datos
myBiomodData

# Plot Formato de los datos
plot(myBiomodData)

# 6.- Opciones de calibracion
#############################

# Consultar help
?BIOMOD_ModelingOptions

# GAM: cambiar k=4 para evitar modelos demasiado complejos
# Maxent: para trabajar con Maxent en biomod2, necesitamos especificar el archivo *.jar. Para eso, antes debe estar instalado JAVA en el ordenador, además de instalar rjava en la version de R en la que se este trabajando
myBiomodOption <- BIOMOD_ModelingOptions(GAM = list(k = 4),
                                         MAXENT.Phillips = list(path_to_maxent.jar="./DATA/maxent/maxent.jar"))
print(myBiomodOption)

# 7.- Calibracion de los modelos
################################

# Numero de rounds de evaluacion = minimo 20
# Data train/test porcentaje = depende del numero de ocurrencias de la especie de estudio
# Prevalence = 0.5 - Ya sea NULL (predeterminado) o un numero entre 0-1
# Importancia de las variables - Nº de permutaciones = 3 (RECOMENDABLE)
# Metricas para la evaluacion y para crear subconjuntos de modelos parciales para construir el conjunto = ROC, TSS, kappa
# Full models = TRUE - Por defecto

# Consultar help
?BIOMOD_Modeling

# Cambiar el directorio de trabajo al path de la salida de los modelos
setwd(path2work)

myBiomodModelOut <- BIOMOD_Modeling( 
  myBiomodData, # Los datos incluyen los registros de la especie, el nombre de la especie y las variables
  models = c('GLM','GBM','GAM','CTA','ANN',
             'FDA','MARS','RF','MAXENT.Phillips.2'), # Metodos (algoritmos)
  models.options = myBiomodOption, # Configurado en el paso anterior
  NbRunEval=1,
  DataSplit=80, # % de los datos utilizados para training
  Yweights=NULL,
  Prevalence=0.5, # Suma ponderada de presencias es igual a la suma ponderada de ausencias
  VarImport=1,
  models.eval.meth = c('TSS','ROC'),
  BinRoc=TRUE, # Los resultados tambien se produciran como valores binarios utilizando el threshold de la curva ROC para todos los modelos
  BinTSS=TRUE, # Los resultados tambien se produciran como valores binarios utilizando el threshold de TSS para todos los modelos
  FiltRoc=TRUE, # Los resultados tambien se produciran como valores filtrados utilizando el threshold de la curva ROC para todos los modelos
  FiltTSS=TRUE, # Los resultados tambien se produciran como valores filtrados utilizando el threshold de TSS para todos los modelos
  SaveObj = TRUE, # Mantener todos los resultados y salidas en el disco duro
  do.full.models = TRUE, # Modelos calibrados y evaluados con todo el conjunto de datos
  modeling.id = "allmodels") # ID (= nombre) del procedimiento de modelado

# Resumen de la modelizacion
myBiomodModelOut

# 8.- Evaluacion de los modelos
###############################

options(max.print=999999) # Dar mas memoria RAM

# Obtener los valores de evaluacion de los modelos sobre el objeto creado anteriormente
myBiomodModelEval <- get_evaluations(myBiomodModelOut)
myBiomodModelEval

# Print los dimnames de este objeto
dimnames(myBiomodModelEval)

# Print ROC scores de todos los modelos seleccionados
print(myBiomodModelEval["ROC","Testing.data",,,])

# Print TSS scores de todos los modelos seleccionados
print(myBiomodModelEval["TSS","Testing.data",,,])

# Obtener importancia de las variables (predictores)
var.imp.mod.out <- get_variables_importance(myBiomodModelOut, as.data.frame=T)
head(var.imp.mod.out, n = 5)

# Guardar la importancia de las variables en formato *.csv
capture.output(get_variables_importance(myBiomodModelOut),
               file=file.path(myRespName, 
                              paste(myRespName,"_variables_importance.csv", sep=",")))

# Guardar las metricas de la evaluacion en formato *.csv
evalDF.ROC <- as.data.frame(myBiomodModelEval["ROC","Testing.data",,,])
evalDF.TSS <- as.data.frame(myBiomodModelEval["TSS","Testing.data",,,])

write.csv(evalDF.ROC, file = paste(getwd(),"/",myRespName,"/",myRespName,"_evalDF_ROC.csv",sep=""))
write.csv(evalDF.TSS, file = paste(getwd(),"/",myRespName,"/",myRespName,"_evalDF_TSS.csv",sep=""))

# 9.- Ensemble modeling - Obtencion de los modelos de conjunto
##############################################################

# Ensamblar todos los modelos parciales-individuales
# Usamos ROC > 0.9 como la metrica del threshold para la seleccion de modelos parciales y constuccion de los modelos de conjunto. Pueden usarse el TSS y el kappa igualmente.
# Los metodos de ensamblaje son: media, mediana y media ponderada
myBiomodEM <- BIOMOD_EnsembleModeling( 
  modeling.output = myBiomodModelOut, # Resultados de los modelos
  chosen.models = 'all', # Modelos que se incluiran al ensamblar - Seleccionamos todos
  em.by='all', # Indicador que define la forma en que se combinaran los modelos para construir los modelos de conjunto
  eval.metric = c('ROC'), # Metrica de evaluacion utilizada para construir los modelos de conjunto
  eval.metric.quality.threshold = 0.9,
  prob.mean = TRUE, # Estimar las probabilidades medias a traves de predicciones
  prob.cv = FALSE, # Coeficiente de variacion
  prob.ci = FALSE, # Intervalos de confianza de la prob. media
  prob.ci.alpha = 0.05, # Nivel de significancia para estimar el intervalo de confianza. Por defecto = 0.05
  prob.median = TRUE, # Estimar la mediana
  prob.mean.weight = TRUE, # Estimar la media ponderada de las probabilidades
  prob.mean.weight.decay = 'proportional' ) # Definir la importancia relativa de los pesos

# Print summary                     
myBiomodEM

# Obtener evaluaciones para los modelos de conjunto
emEvalDF <- as.data.frame(get_evaluations(myBiomodEM))
head(emEvalDF, n = 2)

# Guardar las metricas de la evaluacion en formato *.csv
write.csv(emEvalDF, file = paste(getwd(),"/",myRespName,"/",myRespName,"_EnsMod_evalDF_AllMetrics.csv",sep=""))

# 10.- Proyeccion de modelos "actuales" al espacio geografico
#############################################################

# Proyeccion sobre el area de estudio bajo condiciones actuales
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut, # Resultados de los modelos
  new.env = myExpl, # Variables ambientales
  proj.name = 'current', # Nombre de las proyecciones
  selected.models = 'all', # Modelos para ser proyectados - Todos
  compress = 'gzip',
  output.format = '.grd', # Formato de archivos GIS (tambien *.img)
  do.stack=TRUE)

# Resumen del objeto creado
myBiomodProj

# Archivos creados en el disco duro
list.files("C:/GitHub_repo/OUTPUTS/Pb/proj_current")

# Obtener las predicciones
mod_projPres <- get_predictions(myBiomodProj)

# Calcular una metrica para las predicciones
presentResult <- calc(mod_projPres,fun = median) # Elegir cualquier estadistica descriptiva que se desee

# Plot de los resultados
plot(presentResult)

# Guardar en formato raster la proyeccion espacial de las predicciones
writeRaster(presentResult, filename = "Pb_Present", format = "GTiff", overwrite = T)

# 11.- Proyeccion de modelos de conjunto
########################################

myBiomodEF <- BIOMOD_EnsembleForecasting( 
  EM.output = myBiomodEM,
  projection.output = myBiomodProj,
  binary.meth = c('ROC', 'TSS'),
  filtered.meth = c('ROC', 'TSS'),
  output.format = '.grd')

# Print summary
myBiomodEF

# Plot del objeto
plot(myBiomodEF)

# Archivos creados en el disco duro
list.files("C:/GitHub_repo/OUTPUTS/Pb/proj_current")

# Calcular una metrica para las predicciones
mod_projPresEnsemble <- get_predictions(myBiomodEF);
presentEnsembleResult <- mod_projPresEnsemble[[2]] # Este es la mediana del conjunto de modelos

# # Guardar en formato raster la proyeccion espacial del modelo de conjunto
writeRaster(presentEnsembleResult, filename = "Pb_PresentEnsemble", format = "GTiff", overwrite = T)

## --------------------------------------------------------- ##
## --------------------------------------------------------- ##

## OBTENCION DE GRAFICOS INTERESANTES ##
########################################

# EVALUACION DE LOS MODELOS
###########################

# Grafico de evaluacion por modelo
gg1 <- models_scores_graph( myBiomodModelOut,
                            by = 'models',
                            metrics = c('ROC','TSS') )
# Se observa el rendimiento de cada modelo segun las metricas utilizadas

# CONTRIBUCION DE LAS VARIABLES A LOS MODELOS
#############################################

# Cargar modelos (en este caso RF)
RF.mod <- BIOMOD_LoadModels(myBiomodModelOut,models='RF')

# Obtener importancia de las variables en este modelo
RF.vi <- variables_importance(get(RF.mod[1]), getModelsInputData(myBiomodModelOut,'expl.var'), nb_rand=1)

# Comprobar los resultados
RF.vi$mat

# Convertir este objeto en tabla = data.frame
contributions <- as.data.frame(RF.vi$mat)

# Obtener grafica de barras
barplot(height = t(contributions),
        beside = TRUE,
        horiz = TRUE,
        xlab = "Variable Importance",
        legend = c("RF"))


# GRAFICO DE CURVAS DE RESPUESTA
################################

# Cargar los modelos para los cuales queremos extraer las curvas de respuesta predicha - en nuestro caso de los modelos de conjunto
myModel <- BIOMOD_LoadModels(myBiomodEM, models='mergedData')

# Plot 2D
myRespPlot2D.Pb <- response.plot2(models = myModel[1],
                                  Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                                  show.variables= get_formal_data(myBiomodModelOut,'expl.var.names')[c(1,2,3,4,5,6)], # Especificamos las variables que queremos mostrar
                                  do.bivariate = FALSE,
                                  col = c("red"),
                                  legend = TRUE,
                                  data_species = get_formal_data(myBiomodModelOut,'resp.var'))

# plot 3D response plots
# Aquí sólo para un modelo
myRespPlot3D <- response.plot2(models  = myModel[1],
                               Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                               show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
                               do.bivariate = TRUE,
                               fixed.var.metric = 'median',
                               data_species = get_formal_data(myBiomodModelOut,'resp.var'),
                               legend = F,
                               display_title = FALSE)

# Clicar en Zoom para ver el plot ampliado

# OBTENCION DE MAPAS PREDICTIVOS DE HABITAT
###########################################

# Plot de probabilidad de presencia de la especie
# Cargar el raster correspondiente
Pb_proj_current <- raster("C:/GitHub_repo/OUTPUTS/Pb/proj_current/proj_current_Pb.grd")

dev.off() # Elimina cualquier formato de plot que se haya generado anteriormente y asi evitar incompatibilidades entre paquetes

# Plot
plot(Pb_proj_current)

# Configurar paleta de colores
blues <- brewer.pal('Blues', n = 5)
reds <- rev(brewer.pal('Reds', n = 5))

# Combinacion de colores
myTheme <- rasterTheme(region = c(blues, reds))

# Plot
levelplot(Pb_proj_current, par.settings = BuRdTheme, margin = FALSE, main = expression("Probabilidad de presencia | Pb"))

# Plot de habitat predicho
# Cargar el raster correspondiente
Pb_binary <- raster("C:/GitHub_repo/OUTPUTS/Pb/proj_current/proj_current_Pb_ensemble_TSSbin.grd")

# Plot
plot(Pb_binary)

# Orden alfabetico
Pb.bn <- as.factor(Pb_binary)

# Add una columna con el nuevo atributo
rat <- levels(Pb.bn)[[1]]
rat[["Pb.bn"]] <- c("Unsuitable", "Suitable")
levels(Pb.bn) <- rat

# Plot
Pb.hab <- levelplot(Pb.bn, col.regions=rev(terrain.colors(2)), xlab="", ylab="", main = expression("Habitat Suitability | Pb"))

plot(Pb.hab)

# Guardar como geotiff
if (require(rgdal)) {
  rf <- writeRaster(Pb_binary, filename="C:/GitHub_repo/OUTPUTS/bin_Pb.tif", format="GTiff", overwrite=TRUE)
}
