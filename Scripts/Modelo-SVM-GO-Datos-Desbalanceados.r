# Cargamos los datos

#setwd("D:/OneDrive - fbioyf.unr.edu.ar/Especializ_Bioinf/7-Procesamiento_inteligente_de_datos/Trabajo-Final")

IDS <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/id_tp.txt",
                as.is = TRUE, sep = "\t", header= FALSE)
IDS <- IDS$V1

Clases <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/clases.txt",
                   as.is = TRUE, sep = "\t", header= TRUE)

rownames(Clases) <- Clases$ID

Clases <- Clases[IDS,]

FeaturesPept <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/FeaturesPeptides.txt",
                         as.is = TRUE, sep = "\t", header= TRUE)

FeaturesPfam <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/FeaturesPfam.txt",
                         as.is = TRUE, sep = "\t", header= TRUE)


#----------------------------------------------------------------------------------------
##### Hacemos un SVM para predecir las clases GO3 ####

FeatPeptGO <- merge(Clases[, - ncol(Clases)], FeaturesPept, by="ID")
FeatPfamGO <- merge(Clases[, - ncol(Clases)], FeaturesPfam, by="ID")

#vamos a probar predecir el GO3 que es de lipid metabolism
FeatPeptGO <- FeatPeptGO[,-c(2,3,5,6,7,8,9)]
FeatPfamGO <- FeatPfamGO[,-c(2,3,5,6,7,8,9)]

FeatPeptGO[,2] <- as.factor(FeatPeptGO[,2])
rownames(FeatPeptGO) <- FeatPeptGO$ID
FeatPeptGO <-FeatPeptGO[,-1]

FeatPfamGO[,2] <- as.factor(FeatPfamGO[,2])
rownames(FeatPfamGO) <- FeatPfamGO$ID
FeatPfamGO <- FeatPfamGO[,-1]


#Estructura de los datos
str(FeatPeptGO) 
str(FeatPfamGO)

# Distribución variable respuesta
library(ggplot2)

ggplot(data = FeatPeptGO, aes(x = GO3, y = ..count.., fill = GO3)) +
  geom_bar(color= "black") +
  labs(title = "Distribución de GO3: lipid metabolism") +
  scale_fill_manual(values = c("white", "black"), 
                    labels = c("-", "+")) +
  theme_bw() + theme(plot.title = element_text(hjust = 0.5))

# Tabla frecuencias variable respuesta
table(FeatPeptGO$GO3)


# Particionamos la data en 80% TRAIN y 20% TEST
library(caret)

set.seed(123)

indexPeptGO <- createFolds(t(FeatPeptGO[, "GO3"]), k = 5)
PeptGOTest <- FeatPeptGO[indexPeptGO[[3]], ]
PeptGOTrain <- FeatPeptGO[setdiff(seq(1:dim(FeatPeptGO)[1]), indexPeptGO[[3]]), ]
dim(PeptGOTrain) #885 observaciones en el set de train
dim(PeptGOTest) #222 observaciones en el set de test


indexPfamGO <- createFolds(t(FeatPfamGO[, "GO3"]), k = 5)
PfamGOTest <- FeatPfamGO[indexPfamGO[[3]], ]
PfamGOTrain <- FeatPfamGO[setdiff(seq(1:dim(FeatPfamGO)[1]), indexPfamGO[[3]]), ]
dim(PfamGOTrain) #885 observaciones en el set de train
dim(PfamGOTest) #222 observaciones en el set de test

##### SVM para Features Peptides y clase GO3 #####
library(e1071) 

## SVM lineal ###
## Kernel lineal ##

#Costo = 1
svmLineal <- svm(x= PeptGOTrain[,-1], y = PeptGOTrain[,1], kernel="linear", cost=1,
                 type="C-classification", scale =TRUE) 

summary(svmLineal)
# Number of Support Vectors:  339 ( 225 114 )

#Error de train
confusionMatrix(predict(svmLineal, PeptGOTrain[,-1], type = "class"), PeptGOTrain[,1])
# 'Positive' Class : 0 
#           Reference
#Prediction   0   1
#         0 771 114
#         1   0   0

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.8712
# parece alto, pero en realidad predice todas las muestras como clase 0
# no logra separar las clases.

# Sensitivity (Recall) = TP/(TP+FN) = 1 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) = 0 (seria cuantos predije bien del total de negativos)

#Error de test
confusionMatrix(predict(svmLineal, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])
#            Reference
#Prediction  0   1
#         0 193 29
#         1   0  0

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.8694
# Sensitivity (Recall) = TP/(TP+FN) = 1 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) = 0 (seria cuantos predije bien del total de negativos)

#### Optimización de hiperparámetros KERNEL LINEAL mediante validación cruzada 10-fold#####
costs <- c(0.001, 0.01, 0.1, 1, 10)
set.seed(123)
tuningLineal <- tune(svm, GO3 ~ ., data = PeptGOTrain, cross=10,
                     kernel = "linear", 
                     ranges = list(cost = costs), 
                     scale = TRUE)

summary(tuningLineal)
# - best parameters: cost 10
# - best performance: 0.2684874 

ggplot(data = tuningLineal$performances, aes(x = cost, y = error)) +
  geom_line() +
  geom_point() +
  labs(title = "Error de validación ~ hiperparámetro C",
       x= "Costo", y= "Error de validación") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(hjust = 1,size=12),
        axis.text.y = element_text(size=12))

# Como vemos, el error se mantiene constante para los distintos
# valores de Costo... es decir que no mejora el modelo.

### PLOT DEL SVM###
# Creo una lista con el valor media de cada una de las variables
medias <- list()

for (i in 1:ncol(PeptGOTest[,-1])){
  medias[[colnames(PeptGOTest[,-1][i])]] <- mean(PeptGOTest[,-1][,i])
}

plot(best_linear_model, PeptGOTest, Length ~ MW)
plot(best_linear_model, PeptGOTest, Length ~ pI)

plot(best_linear_model, PeptGOTest, Length ~ MW, slice = medias[-c(1,2)])
# Con slice, fijo los otros features en su valor medio para graficarlos
#  Las X son los vectores soportes. Los 0 son los otros puntos que no afectan
# en el calculo del margen

plot(best_linear_model, PeptGOTest, Length ~ Charge, slice = medias[-c(1,3)])


#### CALCULO DEL ERROR MEDIO PARA LAS 5 PARTICIONES###
matrixError <- matrix(0,5,2) #genero una fila de 5 filas y 2 columnas
colnames(matrixError) <- c("Train", "Test")


for (i in 1:5){
  PeptGOTestKFold <- FeatPeptGO[indexPeptGO[[i]], ] 
  #ME QUEDO CON LA LISTA i COMO SUBSET DE DATOS DE TEST
  PeptGOTrainKFold <- FeatPeptGO[setdiff(seq(1:dim(FeatPeptGO)[1]), indexPeptGO[[i]]), ]
  #ME QUEDO CON EL RESTO COMO SUBSET DE DATOS PARA TRAIN
  fit1 <- svm(GO3 ~ ., data = PeptGOTrainKFold, kernel = "linear", cost = 1, scale = TRUE)
  #AJUSTO EL MODELO CON MIS DATOS DE TRAIN
  predictGO <- predict(fit1, PeptGOTrainKFold[,-1], type = "class")
  #PREDIGO LAS CLASES PARA LOS DATOS DE TRAIN
  matrixError[i,1] <- (1- confusionMatrix(predictGO, PeptGOTrainKFold[,1])$overall['Accuracy'])*100
  #confusionMatrix es una funcion del paquete caret
  #le tengo q pasar las clases predichas y las clases reales.
  #CALCULO EL ERROR DE TRAINING Y LO GUARDO EN LA MATRIZ DE ERROR
  predictGO <- predict(fit1, PeptGOTestKFold[,-1], type = "class")
  #PREDIGO LAS CLASES CON LOS DATOS DE TEST
  matrixError[i,2] <- (1- confusionMatrix(predictGO, PeptGOTestKFold[,1])$overall['Accuracy'])*100
  #CALCULO EL ERROR DE TEST Y LO GUARDO EN LA MATRIZ DE ERROR
}

# Error de train promedio
mean(matrixError[,1]) #12.91776

#Error de test promedio
mean(matrixError[,2]) #12.91723


#### Probamos kernel radial #####

# Costo=1 y Gamma Default
svmRadial <- svm(x= PeptGOTrain[,-1], y = PeptGOTrain[,1], kernel="radial", cost=1,
                 type="C-classification", scale =TRUE) 

summary(svmRadial)
# Number of Support Vectors:  338 ( 224 114)

confusionMatrix(predict(svmRadial, PeptGOTrain[,-1], type = "class"), PeptGOTrain[,1])
#          Reference
#Prediction   0   1
#         0 771 106
#         1   0   8

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.8802
# Sensitivity (Recall) = TP/(TP+FN) =  1 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) =  0.07018 (seria cuantos predije bien del total de negativos)


# Probamos con los datos de test
confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])
#          Reference
#Prediction   0   1
#         0 193  29
#         1   0   0

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.8694
# Sensitivity (Recall) = TP/(TP+FN) =  1 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) =  0 (seria cuantos predije bien del total de negativos)


confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$overall['Accuracy']
confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$byClass['Sensitivity']
confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$byClass['Specificity']

##### Optimizamos el costo y el gamma en simultaneo ####

costs <- c(0.001, 0.01, 0.1, 1 ,10)
gammas <- c(0.001, 0.01, 0.1, 1 ,10, 100)

# sampling method: 10-fold cross validation 
set.seed(123)
tuningRadial <- tune(svm, GO3 ~ ., data = PeptGOTrain, cross= 10,
                     kernel = "radial", 
                     ranges = list(cost = costs, 
                                   gamma = gammas), 
                     scale = TRUE)

summary(tuningRadial)
# best parameters: cost= 10 gamma=1
# best performance: 0.1219356 
# Con un kernel radial, los hiperparámetros que reducen el error de 
# validación son coste = 10, gamma = 1

ggplot(data = tuningRadial$performances, aes(x = cost, y = error, color = factor(gamma))) +
  geom_line() +
  geom_point() +
  labs(title = "Error de validación ~ hiperparámetro C y gamma",
       x= "Costo", y= "Error de validación",
       color="Gamma") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(hjust = 1,size=12),
        axis.text.y = element_text(size=12),
        legend.position = "bottom")

best_radial_model <- tuningRadial$best.model

summary(best_radial_model)
# Number of Support Vectors:  878 ( 764 114 )

# Para los datos de train
confusionMatrix(predict(best_radial_model, PeptGOTrain[,-1], type = "class"), PeptGOTrain[,1])
#            Reference
#Prediction   0   1
#         0 771   0
#         1   0 114
# Accuracy : 1
# Sensitivity : 1
# Specificity : 1

# Para los datos de test
confusionMatrix(predict(best_radial_model, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])
#            Reference
#Prediction   0  1
#         0 192 28
#         1   1  1
# Accuracy : 0.8694
# Sensitivity : 0.99482
# Specificity : 0.03448


# CALCULO DEL ERROR MEDIO PARA LAS 5 PARTICIONES#
matrixErrorRadial <- matrix(0,5,2) #genero una fila de 5 filas y 2 columnas
colnames(matrixErrorRadial) <- c("Train", "Test")


for (i in 1:5){
  PeptGOTestKFold <- FeatPeptGO[indexPeptGO[[i]], ] 
  #ME QUEDO CON LA LISTA i COMO SUBSET DE DATOS DE TEST
  PeptGOTrainKFold <- FeatPeptGO[setdiff(seq(1:dim(FeatPeptGO)[1]), indexPeptGO[[i]]), ]
  #ME QUEDO CON EL RESTO COMO SUBSET DE DATOS PARA TRAIN
  fit1 <- svm(GO3 ~ ., data = PeptGOTrainKFold, kernel = "radial",
              cost = 10, gamma=1, scale = TRUE)
  #AJUSTO EL MODELO CON MIS DATOS DE TRAIN
  predictGO <- predict(fit1, PeptGOTrainKFold[,-1], type = "class")
  #PREDIGO LAS CLASES PARA LOS DATOS DE TRAIN
  matrixErrorRadial[i,1] <- (1- confusionMatrix(predictGO, PeptGOTrainKFold[,1])$overall['Accuracy'])*100
  #confusionMatrix es una funcion del paquete caret
  #le tengo q pasar las clases predichas y las clases reales.
  #CALCULO EL ERROR DE TRAINING Y LO GUARDO EN LA MATRIZ DE ERROR
  predictGO <- predict(fit1, PeptGOTestKFold[,-1], type = "class")
  #PREDIGO LAS CLASES CON LOS DATOS DE TEST
  matrixErrorRadial[i,2] <- (1- confusionMatrix(predictGO, PeptGOTestKFold[,1])$overall['Accuracy'])*100
  #CALCULO EL ERROR DE TEST Y LO GUARDO EN LA MATRIZ DE ERROR
}

# Error de train promedio
mean(matrixErrorRadial[,1]) #0

#Error de test promedio
mean(matrixErrorRadial[,2]) #12.4655

#### Probamos kernel sigmoideal #####

# Costo=1 y Gamma Default
svmSigmoid <- svm(x= PeptGOTrain[,-1], y = PeptGOTrain[,1], kernel="sigmoid", cost=1,
                 type="C-classification", scale =TRUE) 

summary(svmSigmoid)
# Number of Support Vectors:  204 ( 102 102)

confusionMatrix(predict(svmSigmoid, PeptGOTrain[,-1], type = "class"), PeptGOTrain[,1])
#          Reference
#Prediction   0   1
#         0 725  87
#         1  46  27

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.8497
# Sensitivity (Recall) = TP/(TP+FN) =  0.9403 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) =  0.2368 (seria cuantos predije bien del total de negativos)


# Probamos con los datos de test
confusionMatrix(predict(svmSigmoid, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])
#          Reference
#Prediction   0   1
#         0 179  24
#         1  14   5

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.8288
# Sensitivity (Recall) = TP/(TP+FN) =  0.9275 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) =  0.1724 (seria cuantos predije bien del total de negativos)


confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$overall['Accuracy']
confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$byClass['Sensitivity']
confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$byClass['Specificity']

##### Optimizamos el costo y el gamma en simultaneo ####

costs <- c(0.001, 0.01, 0.1, 1 ,10)
gammas <- c(0.001, 0.01, 0.1, 1 ,10, 100)

# sampling method: 10-fold cross validation 
set.seed(123)
tuningSigmoid <- tune(svm, GO3 ~ ., data = PeptGOTrain, cross= 10,
                     kernel = "sigmoid", 
                     ranges = list(cost = costs, 
                                   gamma = gammas), 
                     scale = TRUE)

summary(tuningSigmoid)
# best parameters: cost= 0.1 gamma=0.1
# best performance: 0.1275792 
# Con un kernel sigmoidio, los hiperparámetros que reducen el error de 
# validación son coste = 0.1, gamma = 0.1

ggplot(data = tuningSigmoid$performances, aes(x = cost, y = error, color = factor(gamma))) +
  geom_line() +
  geom_point() +
  labs(title = "Error de validación ~ hiperparámetro C y gamma",
       x= "Costo", y= "Error de validación",
       color="Gamma") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size=12),
        axis.text.x = element_text(hjust = 1,size=12),
        axis.text.y = element_text(size=12),
        legend.position = "bottom")

best_sigmoid_model <- tuningSigmoid$best.model

summary(best_sigmoid_model)
# Number of Support Vectors:  229 ( 115 114 )

# Para los datos de train
confusionMatrix(predict(best_sigmoid_model, PeptGOTrain[,-1], type = "class"), PeptGOTrain[,1])
#            Reference
#Prediction   0   1
#         0 768 111
#         1   3   3
# Accuracy : 0.8712
# Sensitivity : 0.99611
# Specificity : 0.02632

# Para los datos de test
confusionMatrix(predict(best_sigmoid_model, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])
#            Reference
#Prediction   0  1
#         0 191 28
#         1   2  1
# Accuracy : 0.8694
# Sensitivity : 0.98964
# Specificity : 0.03448


# CALCULO DEL ERROR MEDIO PARA LAS 5 PARTICIONES#
matrixErrorRadial <- matrix(0,5,2) #genero una fila de 5 filas y 2 columnas
colnames(matrixErrorRadial) <- c("Train", "Test")


for (i in 1:5){
  PeptGOTestKFold <- FeatPeptGO[indexPeptGO[[i]], ] 
  #ME QUEDO CON LA LISTA i COMO SUBSET DE DATOS DE TEST
  PeptGOTrainKFold <- FeatPeptGO[setdiff(seq(1:dim(FeatPeptGO)[1]), indexPeptGO[[i]]), ]
  #ME QUEDO CON EL RESTO COMO SUBSET DE DATOS PARA TRAIN
  fit1 <- svm(GO3 ~ ., data = PeptGOTrainKFold, kernel = "radial",
              cost = 10, gamma=1, scale = TRUE)
  #AJUSTO EL MODELO CON MIS DATOS DE TRAIN
  predictGO <- predict(fit1, PeptGOTrainKFold[,-1], type = "class")
  #PREDIGO LAS CLASES PARA LOS DATOS DE TRAIN
  matrixErrorRadial[i,1] <- (1- confusionMatrix(predictGO, PeptGOTrainKFold[,1])$overall['Accuracy'])*100
  #confusionMatrix es una funcion del paquete caret
  #le tengo q pasar las clases predichas y las clases reales.
  #CALCULO EL ERROR DE TRAINING Y LO GUARDO EN LA MATRIZ DE ERROR
  predictGO <- predict(fit1, PeptGOTestKFold[,-1], type = "class")
  #PREDIGO LAS CLASES CON LOS DATOS DE TEST
  matrixErrorRadial[i,2] <- (1- confusionMatrix(predictGO, PeptGOTestKFold[,1])$overall['Accuracy'])*100
  #CALCULO EL ERROR DE TEST Y LO GUARDO EN LA MATRIZ DE ERROR
}

# Error de train promedio
mean(matrixErrorRadial[,1]) #0

#Error de test promedio
mean(matrixErrorRadial[,2]) #12.4655






#### OPTIMIZACION A MANO #####
# Optimizamos el costo y el gamma en simultaneo

costs <- c(0.001, 0.01, 0.1, 1 ,10)
gammas <- c(0.001, 0.01, 0.1, 1 ,10, 100)


#ajustarModelos
#Recibe (parametro1= costo, parametro2=gamma) y ajusta un modelo de svm con kernel
#radial y esos parámetros a los datos PeptGOTrain.
#Devuelve un vector numérico, con las métricas precisión, sensibilidad y specificidad
#del modelo ajustado, para el set de datos de test PeptGOTest.

ajustarModelos <- function(parametro1, parametro2){
  svmRadial <- svm(x= PeptGOTrain[,-c(1:8)], y = PeptGOTrain[,3], kernel="radial",
                   gamma = parametro2, 
                   cost= parametro1, type="C-classification", scale =TRUE)
  return (c(confusionMatrix(predict(svmRadial, PeptGOTest[,-c(1:8)], type = "class"), PeptGOTest[,3])$overall['Accuracy'],
            confusionMatrix(predict(svmRadial, PeptGOTest[,-c(1:8)], type = "class"), PeptGOTest[,3])$byClass['Sensitivity'],
            confusionMatrix(predict(svmRadial, PeptGOTest[,-c(1:8)], type = "class"), PeptGOTest[,3])$byClass['Specificity']))
}

metricas <- matrix(data=NA, nrow = 36, ncol= 5, 
                   dimnames = list(c(), 
                                   c("Costo", "Gamma", "Acuraccy", "Sensib", "Specif")))

k = 1
for (i in costs){
  for (j in gammas){
    metricas[k,"Costo"] <- i
    metricas[k,"Gamma"] <- j
    metricas[k,"Acuraccy"] <- ajustarModelos(i, j)[1]
    metricas[k,"Sensib"] <- ajustarModelos(i, j)[2]
    metricas[k,"Specif"] <- ajustarModelos(i, j)[3]
    k = k+1
  }
}

metricas <- as.data.frame(metricas)
max(metricas$Acuraccy) # Mayor Acuraccy = 0.8823529
which(metricas$Acuraccy == max(metricas$Acuraccy)) #28 y #34 (indice de la fila con mayor Acc)
metricas[28,]
metricas[34,]
#   Costo Gamma  Acuraccy    Sensib    Specif
#28   10   1    0.8823529     1     0.07142857
#34   100  1    0.8823529      1     0.07142857

# Vemos que el mejor Acuraccy lo obtenemos con cost = 10 o cost=100, gamma = 1

# Ploteamos acuraccy vs Gamma, para costo = parametro1 = 10
acuraccyPept <- c()
for (i in gammas){
  acuraccyPept <- c(acuraccyPept, ajustarModelos(parametro1=10, parametro2=i)[1])
}

plot(log10(gammas), acuraccyPept)
acuraccyPept[4] #0.8823529

# Vemos que el acuracy es mayor para un valor de gamma=1

svmRadialBest <- svm(x= PeptGOTrain[,-c(1:8)], y = PeptGOTrain[,3], kernel="radial", 
                     gamma = 1, cost= 10, type="C-classification", scale =TRUE)
confusionMatrix(predict(svmRadialBest, PeptGOTest[,-c(1:8)], type = "class"), PeptGOTest[,3])
#'            Reference
#' Prediction  0  1
#'          0 193 26
#'          1  0  2

#### Probamos kernel sigmoid #####

# Optimizamos el costo 

ajustarModelos2 <- function(parametro1){
  svmSigmoid <- svm(x= PeptGOTrain[,-c(1:8)], y = PeptGOTrain[,3], kernel="sigmoid", cost= parametro1, type="C-classification", scale =TRUE)
  return (c(confusionMatrix(predict(svmSigmoid, PeptGOTest[,-c(1:8)], type = "class"), PeptGOTest[,3])$overall['Accuracy'],
            confusionMatrix(predict(svmSigmoid, PeptGOTest[,-c(1:8)], type = "class"), PeptGOTest[,3])$byClass['Sensitivity'],
            confusionMatrix(predict(svmSigmoid, PeptGOTest[,-c(1:8)], type = "class"), PeptGOTest[,3])$byClass['Specificity']))
}

metricas2 <- matrix(data=NA, nrow = 6, ncol= 4, dimnames = list(c(), c("Costo", "Acuraccy", "Sensib", "Specif")))

#Completamos la matriz
k = 1
for (i in costs){
    metricas2[k, "Costo"] <- i
    metricas2[k, "Acuraccy"] <- ajustarModelos2(i)[1]
    metricas2[k, "Sensib"] <- ajustarModelos2(i)[2]
    metricas2[k, "Specif"] <- ajustarModelos2(i)[3]
    k = k + 1
  }


metricas2 <- as.data.frame(metricas2)
max(metricas2$Acuraccy) # Mayor Acuraccy = 0.8733032
which(metricas2$Acuraccy == max(metricas2$Acuraccy)) #1,2,3
metricas2[1,]
#   Costo   Acuraccy Sensib Specif
#1 0.001  0.8733032      1      0
metricas2[3,]
#   Costo   Acuraccy Sensib Specif
#3   0.1  0.8733032      1      0

#usamos costo=0.1

#Chequeo con el costo 0.1

svmSigmoid <- svm(x= PeptGOTrain[,-c(1:8)], y = PeptGOTrain[,3], kernel="sigmoid", cost=0.1,
                  type="C-classification", scale =TRUE) 

summary(svmSigmoid)
#Number of Support Vectors:  232 ( 117 115 )

#test
confusionMatrix(predict(svmSigmoid, PeptGOTest[,-c(1:8)], type = "class"), PeptGOTest[,3])
#Accuracy : 0.8733

#          Reference
#Prediction   0   1
#          0 193  28
#          1   0   0

#nuevamente vemos que no puede diferenciar las clase negativa mientras que acierta en su totalidad con la clase positiva

#-------------------------------SVM kernel radial para Pfam y clase GO3---------------------------------

ajustarModelos_Pfam <- function(parametro1, parametro2){
  svmRadial_Pfam <- svm(x= PfamGOTrain[,-c(1:8)], y = PfamGOTrain[,3], kernel="radial",
                   gamma = parametro2, 
                   cost= parametro1, type="C-classification", scale =TRUE)
  return (c(confusionMatrix(predict(svmRadial_Pfam, PfamGOTest[,-c(1:8)], type = "class"), PfamGOTest[,3])$overall['Accuracy'],
            confusionMatrix(predict(svmRadial_Pfam, PfamGOTest[,-c(1:8)], type = "class"), PfamGOTest[,3])$byClass['Sensitivity'],
            confusionMatrix(predict(svmRadial_Pfam, PfamGOTest[,-c(1:8)], type = "class"), PfamGOTest[,3])$byClass['Specificity']))
}

metricas3 <- matrix(data=NA, nrow = 36, ncol= 5, 
                   dimnames = list(c(), 
                                   c("Costo", "Gamma", "Acuraccy", "Sensib", "Specif")))

k = 1
for (i in costs){
  for (j in gammas){
    metricas3[k,"Costo"] <- i
    metricas3[k,"Gamma"] <- j
    metricas3[k,"Acuraccy"] <- ajustarModelos_Pfam(i, j)[1]
    metricas3[k,"Sensib"] <- ajustarModelos_Pfam(i, j)[2]
    metricas3[k,"Specif"] <- ajustarModelos_Pfam(i, j)[3]
    k = k+1
  }
}

metricas3 <- as.data.frame(metricas3)
max(metricas3$Acuraccy) # Mayor Acuraccy = 0.8828829
which(metricas3$Acuraccy == max(metricas3$Acuraccy))#27
metricas3[27,]

    #Costo Gamma  Acuraccy    Sensib    Specif                     
#27    10   0.1 0.8828829 0.9948187 0.3103448  


# Vemos que el mejor Acuraccy lo obtenemos con cost = 10,gamma = 0.1

# Ploteamos acuraccy vs Gamma, para costo = parametro1 = 10
acuraccyPfam <- c()  
for (i in gammas){
  acuraccyPfam <- c(acuraccyPfam, ajustarModelos_Pfam(parametro1=10, parametro2=i)[1])
}

plot(log10(gammas), acuraccyPfam)
acuraccyPfam[3] #0.8828829

# Vemos que el acuracy es mayor para un valor de gamma=0.1

#chequeamos con los parametros optimizados y obtenemos la matriz de confucion

svmRadialpfam <- svm(x= PfamGOTrain[,-c(1:8)], y = PfamGOTrain[,3], kernel="radial", 
                     gamma = 0.1, cost= 10, type="C-classification", scale =TRUE)
summary(svmRadialpfam)


confusionMatrix(predict(svmRadialpfam, PfamGOTrain[,-c(1:8)], type = "class"), PfamGOTrain[,3])
#Accuracy : 0.9955
#Prediction   0   1
#         0 768   1
#         1   3 113

#Test
confusionMatrix(predict(svmRadialpfam, PfamGOTest[,-c(1:8)], type = "class"), PfamGOTest[,3])

#Accuracy : 0.8829

#'            Reference
#' Prediction  0  1
#'          0 192 20
#'          1  6  9

'''
# Gamma default
svmRadial1 <- svm(x= PeptGOTrain[,-c(1:8)], y = PeptGOTrain[,3], kernel="radial", cost=1,
                  type="C-classification", scale =TRUE) 

summary(svmRadial1)
# Number of Support Vectors:  342 ( 227 115 )

confusionMatrix(predict(svmRadial1, PeptGOTrain[,-c(1:8)], type = "class"), PeptGOTrain[,3])
# Accuracy : 0.8736 ; solo logró diferenciar 3 individuos como clase 1

# Probamos con los datos de test
confusionMatrix(predict(svmRadial1, PeptGOTest[,-c(1:8)], type = "class"), PeptGOTest[,3])
# No logro predecir ninguna como clase 1


# Gamma=1
svmRadial2 <- svm(x= PeptGOTrain[,-c(1:8)], y = PeptGOTrain[,3], kernel="radial", gamma = 1, cost=1,
                  type="C-classification", scale =TRUE) 

summary(svmRadial2) 
# Number of Support Vectors:  882 ( 767 115 )

confusionMatrix(predict(svmRadial2, PeptGOTrain[,-c(1:8)], type = "class"), PeptGOTrain[,3])
# Accuracy : 0.9977  

# Aplicamos este modelo con los datos de Test 
confusionMatrix(predict(svmRadial2, PeptGOTest[,-c(1:8)], type = "class"), PeptGOTest[,3])
# Nuevamente no predijo ninguna clase 1

# Gamma=0.3
svmRadial3 <- svm(x= PeptGOTrain[,-c(1:8)], y = PeptGOTrain[,3], kernel="radial", gamma = 0.3, cost=1,
                  type="C-classification", scale =TRUE) 

summary(svmRadial3) 
# Number of Support Vectors:  649 ( 534 115 )

confusionMatrix(predict(svmRadial3, PeptGOTrain[,-c(1:8)], type = "class"), PeptGOTrain[,3])
# Accuracy : 0.9594  

# Aplicamos este modelo con los datos de Test 
confusionMatrix(predict(svmRadial3, PeptGOTest[,-c(1:8)], type = "class"), PeptGOTest[,3])
# Solo logra predecir 2 de la clase 1


gammas <- seq(0.1, 1, 0.05)

ajustarModelos <- function(parametro){
  svmRadial <- svm(x= PeptGOTrain[,-c(1:8)], y = PeptGOTrain[,3], kernel="radial", gamma = parametro, cost=1,
                   type="C-classification", scale =TRUE)
  return (confusionMatrix(predict(svmRadial, PeptGOTest[,-c(1:8)], type = "class"), PeptGOTest[,3])$overall['Accuracy'])
}

ajustarModelos(0.5)

acuraccy <- c()

for (i in gammas){
  acuraccy <- c(acuraccy, ajustarModelos(i))
}


plot(gammas, acuraccy)
# Vemos que el mejor gamma es 0.3
'''

