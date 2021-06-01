# BALANCEAMOS LOS DATOS

# Cargamos los datos

setwd("D:/OneDrive - fbioyf.unr.edu.ar/Especializ_Bioinf/7-Procesamiento_inteligente_de_datos/Trabajo-Final")

IDS <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/id_tp.txt",
                as.is = TRUE, sep = "\t", header= FALSE)

IDS <- IDS$V1 # Creamos un array con los IDS

Clases <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/clases.txt",
                   as.is = TRUE, sep = "\t", header= TRUE)
head(Clases)

rownames(Clases) <- Clases$ID

Clases <- Clases[IDS,] # Nos quedamos solo con los IDS que tenían Pfam

FeaturesPept <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/FeaturesPeptides.txt",
                         as.is = TRUE, sep = "\t", header= TRUE)

rownames(FeaturesPept) <- FeaturesPept$ID

FeaturesPfam <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/FeaturesPfam.txt",
                         as.is = TRUE, sep = "\t", header= TRUE)

rownames(FeaturesPfam) <- FeaturesPfam$ID

##### BALANCEO DE DATOS #####
# Para comenzar, dejamos los 145 datos de la clase GO3, 
# y dejamos solo 290 negativos al azar de las otras clases.

set.seed(123)

IDS_GO3 <- IDS[Clases$GO3 == 1]
IDS_notGO3 <- IDS[Clases$GO3 != 1]
IDS_notGO3_reducidos <- sample(IDS_notGO3, 290, replace=FALSE)
IDS_balanceados <- c(IDS_GO3, IDS_notGO3_reducidos)

Clases_Red <- Clases[IDS_balanceados,] # Nos quedamos solo con los IDS que tenían Pfam
FeaturesPept_Red <- FeaturesPept[IDS_balanceados,]
FeaturesPfam_Red <- FeaturesPfam[IDS_balanceados,]

# Guardamos un archivo con los IDs balanceados

write.table(IDS_balanceados, file = 'id_balanceados.txt', append = FALSE, 
            sep=" ", dec=".", row.names = FALSE, col.names=FALSE)

##### Hacemos un SVM para predecir la clase GO3 = lipid metabolism -----

FeatPeptGO <- merge(Clases_Red[, c("ID", "GO3")], FeaturesPept_Red, by="ID")
FeatPfamGO <- merge(Clases_Red[, c("ID", "GO3")], FeaturesPfam_Red, by="ID")


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
dim(PeptGOTrain) #347 observaciones en el set de train
dim(PeptGOTest) #86 observaciones en el set de test

indexPfamGO <- createFolds(t(FeatPfamGO[, "GO3"]), k = 5)
PfamGOTest <- FeatPfamGO[indexPfamGO[[3]], ]
PfamGOTrain <- FeatPfamGO[setdiff(seq(1:dim(FeatPfamGO)[1]), indexPfamGO[[3]]), ]
dim(PfamGOTrain) #346 observaciones en el set de train
dim(PfamGOTest) #87 observaciones en el set de test


#### COMENZAMOS AJUSTANDO MODELOS DE SVM CON LOS FEATURES EXTRAIDOS DE PEPTIDES ####
library(e1071)

## Kernel lineal ##
svmLineal <- svm(x= PeptGOTrain[,-1], y = PeptGOTrain[,1], kernel="linear", cost=1,
                  type="C-classification", scale =TRUE) 

summary(svmLineal)
# Number of Support Vectors:  210 ( 108 102 )

#Error de train
confusionMatrix(predict(svmLineal, PeptGOTrain[,-1], type = "class"), PeptGOTrain[,1])
# 'Positive' Class : 0 
#           Reference
#Prediction   0   1
#         0 204  49
#         1  28  66

# Accuracy (Precision) : Los predichos bien/ TOTAL (204+66)/(204+28+49+66) = 0.7781
# Sensitivity (Recall) = TP/(TP+FN) = 204/(204+28) = 0.8793 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) = 66/(49+66)= 0.5739 (seria cuantos predije bien del total de negativos)

#Error de test
confusionMatrix(predict(svmLineal, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])
#            Reference
#Prediction  0  1
#         0 48 18
#         1 10 10

# Accuracy (Precision) : Los predichos bien/ TOTAL (58+9)/(58+9+20+0) = 0.6744
# Sensitivity (Recall) = TP/(TP+FN) = 0.8276 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) = 0.3571 (seria cuantos predije bien del total de negativos)


#### OPTIMIZACION DEL PARAMETRO COSTO PARA SVM LINEAL, A MANO...####
costs <- c(0.001, 0.01, 0.1, 1 ,10, 100)

ajustarModelosLineal <- function(costo){
  svmLineal <- svm(x= PeptGOTrain[,-1], y = PeptGOTrain[,1], kernel="linear",
                   cost= costo, type="C-classification", scale =TRUE)
  return (c(confusionMatrix(predict(svmLineal, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$overall['Accuracy'],
            confusionMatrix(predict(svmLineal, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$byClass['Sensitivity'],
            confusionMatrix(predict(svmLineal, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$byClass['Specificity']))
}

metricasCosto <- matrix(data=NA, nrow = 6, ncol= 4, 
                   dimnames = list(c(), 
                                   c("Costo", "Acuraccy", "Sensib", "Specif")))

k = 1
for (i in costs){
  metricasCosto[k,"Costo"] <- i
  metricasCosto[k,"Acuraccy"] <- ajustarModelosLineal(i)[1]
  metricasCosto[k,"Sensib"] <- ajustarModelosLineal(i)[2]
  metricasCosto[k,"Specif"] <- ajustarModelosLineal(i)[3]
  k = k+1
}


metricasCosto <- as.data.frame(metricasCosto)
max(metricasCosto$Acuraccy) # Mayor Acuraccy = 0.6744186
which(metricasCosto$Acuraccy == max(metricasCosto$Acuraccy)) #1 3 4 6 (indice de la fila con mayor Acc)
metricasCosto[c(1,3,4,6),]
# Costo  Acuraccy     Sensib    Specif
#1 1e-03 0.6744186 1.0000000 0.0000000
#3 1e-01 0.6744186 0.8793103 0.2500000
#4 1e+00 0.6744186 0.8275862 0.3571429
#6 1e+02 0.6744186 0.7931034 0.4285714
# Vemos que el mejor Acuraccy lo obtenemos con cost = 0.001, 0.1, 1, 100

# Ploteamos acuraccy vs costo
acuraccyPeptLineal <- c()
for (i in costs){
  acuraccyPeptLineal <- c(acuraccyPeptLineal, ajustarModelosLineal(costo=i)[1])
}

plot(log10(costs), acuraccyPeptLineal)
acuraccyPept[6] #0.7356322 
# Vemos que el mejor cost es 100


#### Optimización de hiperparámetros KERNEL LINEAL mediante validación cruzada 10-fold#####
costs <- c(0.001, 0.01, 0.1, 1 ,10, 100)
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

#El valor de coste que resulta en el menor error de validación (0.2625210) es 1.

tuningLineal$best.model
  
# Almacenamos el modelo optimo obtenido y accedemos a su información
best_linear_model <- tuningLineal$best.model

confusionMatrix(predict(best_linear_model, PeptGOTrain[,-1], type = "class"), PeptGOTrain[,1])
#  Accuracy train: 0.7608
#            Reference
#Prediction   0  1
#         0 197 48
#         1  35 67

confusionMatrix(predict(best_linear_model, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])
#            Reference
#Prediction  0  1
#         0 46 17
#         1 12 11
# Accuracy test : 0.6628  
# Sensitivity : 0.7931          
# Specificity : 0.3929

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
matrixError <- matrix(0,5,6) #genero una fila de 5 filas y 6 columnas
colnames(matrixError) <- c("Error-Train", "Especif-Train", "Sensib-Train",
                           "ErrorTest", "Especif-Test", "Sensib-Test")


for (i in 1:5){
  PeptGOTestKFold <- FeatPeptGO[indexPeptGO[[i]], ] 
  #ME QUEDO CON LA LISTA i COMO SUBSET DE DATOS DE TEST
  PeptGOTrainKFold <- FeatPeptGO[setdiff(seq(1:dim(FeatPeptGO)[1]), indexPeptGO[[i]]), ]
  #ME QUEDO CON EL RESTO COMO SUBSET DE DATOS PARA TRAIN
  fit1 <- svm(GO3 ~ ., data = PeptGOTrainKFold, kernel = "linear", cost = 10, scale = TRUE)
  #AJUSTO EL MODELO CON MIS DATOS DE TRAIN
  predictGO <- predict(fit1, PeptGOTrainKFold[,-1], type = "class")
  #PREDIGO LAS CLASES PARA LOS DATOS DE TRAIN
  matrixError[i,1] <- (1- confusionMatrix(predictGO, PeptGOTrainKFold[,1])$overall['Accuracy'])*100
  #confusionMatrix es una funcion del paquete caret
  #le tengo q pasar las clases predichas y las clases reales.
  #CALCULO EL ERROR DE TRAINING Y LO GUARDO EN LA MATRIZ DE ERROR
  matrixError[i,2] <- confusionMatrix(predictGO, PeptGOTrainKFold[,1])$byClass['Specificity']
  matrixError[i,3] <- confusionMatrix(predictGO, PeptGOTrainKFold[,1])$byClass['Sensitivity']
  predictGO <- predict(fit1, PeptGOTestKFold[,-1], type = "class")
  #PREDIGO LAS CLASES CON LOS DATOS DE TEST
  matrixError[i,4] <- (1- confusionMatrix(predictGO, PeptGOTestKFold[,1])$overall['Accuracy'])*100
  #CALCULO EL ERROR DE TEST Y LO GUARDO EN LA MATRIZ DE ERROR
  matrixError[i,5] <- confusionMatrix(predictGO, PeptGOTestKFold[,1])$byClass['Specificity']
  matrixError[i,6] <- confusionMatrix(predictGO, PeptGOTestKFold[,1])$byClass['Sensitivity']
}

# Error de train promedio
mean(matrixError[,1]) #24.42205

# Especificidad de train promedio
mean(matrixError[,2]) #0.5191

#Error de test promedio
mean(matrixError[,4]) #27.0355; Acurracy promedio : 72.97%

# Especificidad de test promedio
mean(matrixError[,5]) #0.467734



### Hago la matriz de error LINEAL promedio usando el costo 100 como optimo####

matrixError2 <- matrix(0,5,2) #genero una fila de 5 filas y 2 columnas
colnames(matrixError2) <- c("Train", "Test")


for (i in 1:5){
  PeptGOTestKFold <- FeatPeptGO[indexPeptGO[[i]], ] 
  #ME QUEDO CON LA LISTA i COMO SUBSET DE DATOS DE TEST
  PeptGOTrainKFold <- FeatPeptGO[setdiff(seq(1:dim(FeatPeptGO)[1]), indexPeptGO[[i]]), ]
  #ME QUEDO CON EL RESTO COMO SUBSET DE DATOS PARA TRAIN
  fit1 <- svm(GO3 ~ ., data = PeptGOTrainKFold, kernel = "linear", cost = 100, scale = TRUE)
  #AJUSTO EL MODELO CON MIS DATOS DE TRAIN
  predictGO <- predict(fit1, PeptGOTrainKFold[,-1], type = "class")
  #PREDIGO LAS CLASES PARA LOS DATOS DE TRAIN
  matrixError2[i,1] <- (1- confusionMatrix(predictGO, PeptGOTrainKFold[,1])$overall['Accuracy'])*100
  #confusionMatrix es una funcion del paquete caret
  #le tengo q pasar las clases predichas y las clases reales.
  #CALCULO EL ERROR DE TRAINING Y LO GUARDO EN LA MATRIZ DE ERROR
  predictGO <- predict(fit1, PeptGOTestKFold[,-1], type = "class")
  #PREDIGO LAS CLASES CON LOS DATOS DE TEST
  matrixError2[i,2] <- (1- confusionMatrix(predictGO, PeptGOTestKFold[,1])$overall['Accuracy'])*100
  #CALCULO EL ERROR DE TEST Y LO GUARDO EN LA MATRIZ DE ERROR
}

# Error de train promedio
mean(matrixError2[,1]) #23.96079

#Error de test promedio
mean(matrixError2[,2]) #27.95242




####### Kernel Radial#####
# Gamma Default = 1
svmRadial <- svm(x= PeptGOTrain[,-1], y = PeptGOTrain[,1], kernel="radial", cost=1,
                  type="C-classification", scale =TRUE) 

summary(svmRadial)
# Number of Support Vectors:  245 ( 138 107)

confusionMatrix(predict(svmRadial, PeptGOTrain[,-1], type = "class"), PeptGOTrain[,1])
#          Reference
#Prediction   0   1
#         0 219  57
#         1  13  58

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.7983
# Sensitivity (Recall) = TP/(TP+FN) =  0.9440 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) =  0.5043 (seria cuantos predije bien del total de negativos)


# Probamos con los datos de test
confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])

confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$overall['Accuracy']
confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$byClass['Sensitivity']
confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$byClass['Specificity']

#           Reference
#Prediction  0  1
#         0 52 18
#         1  6 10

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.7209
# Sensitivity (Recall) = TP/(TP+FN) = 0.8966 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) = 0.3571 (seria cuantos predije bien del total de negativos)

##### Optimizamos Kerner Radial: costo y el gamma en simultaneo ####

costs <- c(0.001, 0.01, 0.1, 1 ,10, 100)
gammas <- c(0.001, 0.01, 0.1, 1 ,10, 100)

# sampling method: 10-fold cross validation 
set.seed(123)
tuningRadial <- tune(svm, GO3 ~ ., data = PeptGOTrain, cross= 10,
               kernel = "radial", 
               ranges = list(cost = costs, 
                             gamma = gammas), 
               scale = TRUE)

summary(tuningRadial)
# best parameters: cost= 1 gamma=0.1
# best performance: 0.2565546 
# Con un kernel radial, los hiperparámetros que reducen el error de 
# validación son coste = 1, gamma = 0.1.

ggplot(data = tuningRadial$performances, aes(x = cost, y = error, color = factor(gamma))) +
  geom_line() +
  geom_point() +
  labs(title = "Error de validación ~ hiperparámetro C y gamma",
       x= "Costo", y= "Error de validación") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(hjust = 1,size=12),
        axis.text.y = element_text(size=12),
        legend.position = "bottom")

best_radial_model <- svm(GO3 ~ . , data = PeptGOTrain, 
                   kernel = "radial", 
                   cost = 1, 
                   gamma = 0.1, 
                   scale = TRUE)

summary(best_radial_model)
# Number of Support Vectors:  274 ( 164 110 )

# Para los datos de train
confusionMatrix(predict(best_radial_model, PeptGOTrain[,-1], type = "class"), PeptGOTrain[,1])
#            Reference
#Prediction   0  1
#         0 222 41
#         1  10 74
# Accuracy : 0.853
# Sensitivity : 0.9569
# Specificity : 0.6435

# Para los datos de test
confusionMatrix(predict(best_radial_model, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])
#            Reference
#Prediction  0  1
#         0 52 17
#         1  6 11
# Accuracy : 0.7326
# Sensitivity : 0.8966
# Specificity : 0.3929


# CALCULO DEL ERROR MEDIO PARA LAS 5 PARTICIONES#
matrixErrorRadial <- matrix(0,5,6) #genero una fila de 5 filas y 6 columnas
colnames(matrixErrorRadial) <- c("Error-Train", "Especif-Train", "Sensib-Train",
                           "ErrorTest", "Especif-Test", "Sensib-Test")

for (i in 1:5){
  PeptGOTestKFold <- FeatPeptGO[indexPeptGO[[i]], ] 
  #ME QUEDO CON LA LISTA i COMO SUBSET DE DATOS DE TEST
  PeptGOTrainKFold <- FeatPeptGO[setdiff(seq(1:dim(FeatPeptGO)[1]), indexPeptGO[[i]]), ]
  #ME QUEDO CON EL RESTO COMO SUBSET DE DATOS PARA TRAIN
  fit1 <- svm(GO3 ~ ., data = PeptGOTrainKFold, kernel = "radial",
              cost = 1, gamma=0.1, scale = TRUE)
  #AJUSTO EL MODELO CON MIS DATOS DE TRAIN
  predictGO <- predict(fit1, PeptGOTrainKFold[,-1], type = "class")
  #PREDIGO LAS CLASES PARA LOS DATOS DE TRAIN
  matrixErrorRadial[i,1] <- (1- confusionMatrix(predictGO, PeptGOTrainKFold[,1])$overall['Accuracy'])*100
  #confusionMatrix es una funcion del paquete caret
  #le tengo q pasar las clases predichas y las clases reales.
  #CALCULO EL ERROR DE TRAINING Y LO GUARDO EN LA MATRIZ DE ERROR
  matrixErrorRadial[i,2] <- confusionMatrix(predictGO, PeptGOTrainKFold[,1])$byClass['Specificity']
  matrixErrorRadial[i,3] <- confusionMatrix(predictGO, PeptGOTrainKFold[,1])$byClass['Sensitivity']
  predictGO <- predict(fit1, PeptGOTestKFold[,-1], type = "class")
  #PREDIGO LAS CLASES CON LOS DATOS DE TEST
  matrixErrorRadial[i,4] <- (1- confusionMatrix(predictGO, PeptGOTestKFold[,1])$overall['Accuracy'])*100
  #CALCULO EL ERROR DE TEST Y LO GUARDO EN LA MATRIZ DE ERROR
  matrixErrorRadial[i,5] <- confusionMatrix(predictGO, PeptGOTestKFold[,1])$byClass['Specificity']
  matrixErrorRadial[i,6] <- confusionMatrix(predictGO, PeptGOTestKFold[,1])$byClass['Sensitivity']
}

# Error de train promedio
mean(matrixErrorRadial[,1]) #13.79837

# Especificidad de train promedio
mean(matrixErrorRadial[,2]) #0.6591457

#Error de test promedio
mean(matrixErrorRadial[,4]) #26.54905

# Especificidad de test promedio
mean(matrixErrorRadial[,5]) #0.3849754


#### Probamos kernel sigmoideal #####

# Costo=1 y Gamma Default
svmSigmoid <- svm(x= PeptGOTrain[,-1], y = PeptGOTrain[,1], kernel="sigmoid", cost=1,
                  type="C-classification", scale =TRUE) 

summary(svmSigmoid)
# Number of Support Vectors:  157 ( 80 77 )

confusionMatrix(predict(svmSigmoid, PeptGOTrain[,-1], type = "class"), PeptGOTrain[,1])
#          Reference
#Prediction   0   1
#         0 185  65
#         1  47  50

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.6772
# Sensitivity (Recall) = TP/(TP+FN) =  0.7974 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) =  0.4348 (seria cuantos predije bien del total de negativos)


# Probamos con los datos de test
confusionMatrix(predict(svmSigmoid, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])
#          Reference
#Prediction   0   1
#         0  43  21
#         1  15   7

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.5814
# Sensitivity (Recall) = TP/(TP+FN) =  0.5814 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) =  0.2500 (seria cuantos predije bien del total de negativos)


confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$overall['Accuracy']
confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$byClass['Sensitivity']
confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$byClass['Specificity']

##### Optimizamos Kernel sigmoideal: costo y el gamma en simultaneo ####

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
# best performance: 0.2884034 
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
# Number of Support Vectors:  228 ( 114 114 )

# Para los datos de train
confusionMatrix(predict(best_sigmoid_model, PeptGOTrain[,-1], type = "class"), PeptGOTrain[,1])
#            Reference
#Prediction   0   1
#         0 209  83
#         1  23  32
# Accuracy : 0.6945
# Sensitivity : 0.9009
# Specificity : 0.2783

# Para los datos de test
confusionMatrix(predict(best_sigmoid_model, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])
#            Reference
#Prediction   0  1
#         0  53 22
#         1   5  6
# Accuracy : 0.686
# Sensitivity : 0.9138
# Specificity : 0.2143


#### CALCULO DEL ERROR MEDIO PARA LAS 5 PARTICIONES###
matrixErrorSigmoid <- matrix(0,5,6) #genero una fila de 5 filas y 6 columnas
colnames(matrixErrorSigmoid) <- c("Error-Train", "Especif-Train", "Sensib-Train",
                                 "ErrorTest", "Especif-Test", "Sensib-Test")

for (i in 1:5){
  PeptGOTestKFold <- FeatPeptGO[indexPeptGO[[i]], ] 
  #ME QUEDO CON LA LISTA i COMO SUBSET DE DATOS DE TEST
  PeptGOTrainKFold <- FeatPeptGO[setdiff(seq(1:dim(FeatPeptGO)[1]), indexPeptGO[[i]]), ]
  #ME QUEDO CON EL RESTO COMO SUBSET DE DATOS PARA TRAIN
  fit1 <- svm(GO3 ~ ., data = PeptGOTrainKFold, kernel = "sigmoid",
              cost = 0.1, gamma=0.1, scale = TRUE)
  #AJUSTO EL MODELO CON MIS DATOS DE TRAIN
  predictGO <- predict(fit1, PeptGOTrainKFold[,-1], type = "class")
  #PREDIGO LAS CLASES PARA LOS DATOS DE TRAIN
  matrixErrorSigmoid[i,1] <- (1- confusionMatrix(predictGO, PeptGOTrainKFold[,1])$overall['Accuracy'])*100
  #confusionMatrix es una funcion del paquete caret
  #le tengo q pasar las clases predichas y las clases reales.
  #CALCULO EL ERROR DE TRAINING Y LO GUARDO EN LA MATRIZ DE ERROR
  matrixErrorSigmoid[i,2] <- confusionMatrix(predictGO, PeptGOTrainKFold[,1])$byClass['Specificity']
  matrixErrorSigmoid[i,3] <- confusionMatrix(predictGO, PeptGOTrainKFold[,1])$byClass['Sensitivity']
  predictGO <- predict(fit1, PeptGOTestKFold[,-1], type = "class")
  #PREDIGO LAS CLASES CON LOS DATOS DE TEST
  matrixErrorSigmoid[i,4] <- (1- confusionMatrix(predictGO, PeptGOTestKFold[,1])$overall['Accuracy'])*100
  #CALCULO EL ERROR DE TEST Y LO GUARDO EN LA MATRIZ DE ERROR
  matrixErrorSigmoid[i,5] <- confusionMatrix(predictGO, PeptGOTestKFold[,1])$byClass['Specificity']
  matrixErrorSigmoid[i,6] <- confusionMatrix(predictGO, PeptGOTestKFold[,1])$byClass['Sensitivity']
}

# Error de train promedio
mean(matrixErrorSigmoid[,1]) #31.06

# Especificidad de train promedio
mean(matrixErrorSigmoid[,2]) #0.1886

#Error de test promedio
mean(matrixErrorSigmoid[,4]) #30.25394

# Especificidad de test promedio
mean(matrixErrorSigmoid[,5]) #0.2098522

##### SVM POLINOMICO #####
# Probamos Costo = 1, Grado = 2
svmPolinomico <- svm(x= PeptGOTrain[,-1], y = PeptGOTrain[,1], kernel="polynomial", cost=1,
                     degree = 2, type="C-classification", scale =TRUE) 

summary(svmPolinomico)
# Number of Support Vectors:  236 ( 124 112)

#Error de traim
confusionMatrix(predict(svmPolinomico, PeptGOTrain[,-1], type = "class"), PeptGOTrain[,1])
#          Reference
#Prediction   0   1
#         0 231 100
#         1   1  15

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.7089
# Sensitivity (Recall) = TP/(TP+FN) =  0.9957 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) =  0.1304 (seria cuantos predije bien del total de negativos)


# Probamos con los datos de test
confusionMatrix(predict(svmPolinomico, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])

#           Reference
#Prediction  0  1
#         0 56 26
#         1  2  2

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.6744
# Sensitivity (Recall) = TP/(TP+FN) = 0.96552 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) = 0.07143 (seria cuantos predije bien del total de negativos)

##### Optimizamos el costo y el grado en simultaneo ####

costs <- c(0.001, 0.01, 0.1, 1 ,10, 100)
grados <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)

# sampling method: 10-fold cross validation 
set.seed(123)
tuningPolinomico <- tune(svm, GO3 ~ ., data = PeptGOTrain, cross= 10,
                     kernel = "polynomial", 
                     ranges = list(cost = costs, 
                                   degrees = grados), 
                     scale = TRUE)

summary(tuningPolinomico)
# best parameters: cost= 10 degrees= 0.5
# best performance: 0.2828571
# Con un kernel polinomico, los hiperparámetros que reducen el error de 
# validación son coste = 10, degrees = 0.5.

ggplot(data = tuningPolinomico$performances, aes(x = cost, y = error, color = factor(degrees))) +
  geom_line() +
  geom_point() +
  labs(title = "Error de validación ~ hiperparámetro C y grado",
       x= "Costo", y= "Error de validación") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(legend.position = "bottom")

# VEMOS EN EL GRAFICO Q CUALQUIER GRADO ES EXACTAMENTE LO MISMO...
best_polynomial_model <- tuningPolinomico$best.model
summary(best_polynomial_model)
# Costo = 10, degree = 3
# Number of Support Vectors:  221 ( 164 110 )

# Para los datos de train
confusionMatrix(predict(best_polynomial_model, PeptGOTrain[,-1], type = "class"), PeptGOTrain[,1])
#            Reference
#Prediction   0  1
#         0 230 50
#         1   2 65
# Accuracy : 0.8501
# Sensitivity : 0.9914
# Specificity : 0.5652

# Para los datos de test
confusionMatrix(predict(best_polynomial_model, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])
#            Reference
#Prediction  0  1
#         0 52 19
#         1  6  9
# Accuracy : 0.7093
# Sensitivity : 0.8966
# Specificity : 0.3214


# CALCULO DEL ERROR MEDIO PARA LAS 5 PARTICIONES#
matrixErrorPolinomico <- matrix(0,5,2) #genero una fila de 5 filas y 2 columnas
colnames(matrixErrorPolinomico) <- c("Train", "Test")


for (i in 1:5){
  PeptGOTestKFold <- FeatPeptGO[indexPeptGO[[i]], ] 
  #ME QUEDO CON LA LISTA i COMO SUBSET DE DATOS DE TEST
  PeptGOTrainKFold <- FeatPeptGO[setdiff(seq(1:dim(FeatPeptGO)[1]), indexPeptGO[[i]]), ]
  #ME QUEDO CON EL RESTO COMO SUBSET DE DATOS PARA TRAIN
  fit1 <- svm(GO3 ~ ., data = PeptGOTrainKFold, kernel = "polynomial",
              cost = 10, degree=3, scale = TRUE)
  #AJUSTO EL MODELO CON MIS DATOS DE TRAIN
  predictGO <- predict(fit1, PeptGOTrainKFold[,-1], type = "class")
  #PREDIGO LAS CLASES PARA LOS DATOS DE TRAIN
  matrixErrorPolinomico[i,1] <- (1- confusionMatrix(predictGO, PeptGOTrainKFold[,1])$overall['Accuracy'])*100
  #confusionMatrix es una funcion del paquete caret
  #le tengo q pasar las clases predichas y las clases reales.
  #CALCULO EL ERROR DE TRAINING Y LO GUARDO EN LA MATRIZ DE ERROR
  predictGO <- predict(fit1, PeptGOTestKFold[,-1], type = "class")
  #PREDIGO LAS CLASES CON LOS DATOS DE TEST
  matrixErrorPolinomico[i,2] <- (1- confusionMatrix(predictGO, PeptGOTestKFold[,1])$overall['Accuracy'])*100
  #CALCULO EL ERROR DE TEST Y LO GUARDO EN LA MATRIZ DE ERROR
}

# Error de train promedio
mean(matrixErrorPolinomico[,1]) #15.82016

#Error de test promedio
mean(matrixErrorPolinomico[,2]) #29.55627; Acurracy promedio : 73.45%



##### TABLA CON ERRORES PROMEDIO Pept SVM DISTINTOS KERNEL PARAMETROS OPTIMIZADOS ####
ResumenErrores <- matrix(data=NA, nrow = 3, ncol= 6, 
                         dimnames = list(c("Lineal", "Radial", "Sigmoideo"), 
                                         c("ErrorTrain-promedio", "Especificidad-Train-Promedio",
                                           "Sensibilidad-Train-Promedio",
                                           "ErrorTest-promedio", "Especificidad-Test-Promedio",
                                           "Sensibilidad-Test-Promedio")))

ResumenErrores["Lineal","ErrorTrain-promedio"] <- mean(matrixError[,1])
ResumenErrores["Lineal","Especificidad-Train-Promedio"] <- mean(matrixError[,2])
ResumenErrores["Lineal","Sensibilidad-Train-Promedio"] <- mean(matrixError[,3])
ResumenErrores["Lineal","ErrorTest-promedio"] <- mean(matrixError[,4])
ResumenErrores["Lineal","Especificidad-Test-Promedio"] <- mean(matrixError[,5])
ResumenErrores["Lineal","Sensibilidad-Test-Promedio"] <- mean(matrixError[,6])

ResumenErrores["Radial","ErrorTrain-promedio"] <- mean(matrixErrorRadial[,1])
ResumenErrores["Radial","Especificidad-Train-Promedio"] <- mean(matrixErrorRadial[,2])
ResumenErrores["Radial","Sensibilidad-Train-Promedio"] <- mean(matrixErrorRadial[,3])
ResumenErrores["Radial","ErrorTest-promedio"] <- mean(matrixErrorRadial[,4])
ResumenErrores["Radial","Especificidad-Test-Promedio"] <- mean(matrixErrorRadial[,5])
ResumenErrores["Radial","Sensibilidad-Test-Promedio"] <- mean(matrixErrorRadial[,6])

ResumenErrores["Sigmoideo","ErrorTrain-promedio"] <- mean(matrixErrorSigmoid[,1])
ResumenErrores["Sigmoideo","Especificidad-Train-Promedio"] <- mean(matrixErrorSigmoid[,2])
ResumenErrores["Sigmoideo","Sensibilidad-Train-Promedio"] <- mean(matrixErrorSigmoid[,3])
ResumenErrores["Sigmoideo","ErrorTest-promedio"] <- mean(matrixErrorSigmoid[,4])
ResumenErrores["Sigmoideo","Especificidad-Test-Promedio"] <- mean(matrixErrorSigmoid[,5])
ResumenErrores["Sigmoideo","Sensibilidad-Test-Promedio"] <- mean(matrixErrorSigmoid[,6])

ResumenErrores  
  

###### HECHO A MANO #####
#'ajustarModelos
#' Recibe (parametro1= costo, parametro2=gamma) y ajusta un modelo de svm con kernel
#' radial y esos parámetros a los datos PeptGOTrain.
#' Devuelve un vector numérico, con las métricas precisión, sensibilidad y specificidad
#' del modelo ajustado, para el set de datos de test PeptGOTest.

ajustarModelos <- function(parametro1, parametro2){
  svmRadial <- svm(x= PeptGOTrain[,-1], y = PeptGOTrain[,1], kernel="radial",
                   gamma = parametro2, 
                   cost= parametro1, type="C-classification", scale =TRUE)
  return (c(confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$overall['Accuracy'],
            confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$byClass['Sensitivity'],
            confusionMatrix(predict(svmRadial, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])$byClass['Specificity']))
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
max(metricas$Acuraccy) # Mayor Acuraccy = 0.7816092
which(metricas$Acuraccy == max(metricas$Acuraccy)) #21 (indice de la fila con mayor Acc)
metricas[21,]
#   Costo Gamma  Acuraccy    Sensib    Specif
#21     1   0.1 0.7816092 0.9482759 0.4482759
# Vemos que el mejor Acuraccy lo obtenemos con cost = 1, gamma = 0.1

# Ploteamos acuraccy vs Gamma, para costo = 1
acuraccyPept <- c()
for (i in gammas){
  acuraccyPept <- c(acuraccyPept, ajustarModelos(parametro1=1, parametro2=i)[1])
}

plot(log10(gammas), acuraccyPept)
acuraccyPept[3] #0.7816092
# Vemos que el mejor gamma es 0.1

svmRadialBest <- svm(x= PeptGOTrain[,-1], y = PeptGOTrain[,1], kernel="radial", 
                     gamma = 0.1, cost= 1, type="C-classification", scale =TRUE)
confusionMatrix(predict(svmRadialBest, PeptGOTest[,-1], type = "class"), PeptGOTest[,1])
#'            Reference
#' Prediction  0  1
#'          0 55 16
#'          1  3 13

#### AJUSTAMOS MODELOS CON LOS FEATURES EXTRAIDOS DE PFam ####
### Kernel lineal Pfam ####
svmLinealPfam <- svm(x= PfamGOTrain[,-1], y = PfamGOTrain[,1], kernel="linear", cost=1,
                 type="C-classification", scale =TRUE) 

summary(svmLinealPfam)
# Number of Support Vectors:  282 ( 185 97 )

confusionMatrix(predict(svmLinealPfam, PfamGOTrain[,-1], type = "class"), PfamGOTrain[,1])
# 'Positive' Class : 0 
#           Reference
#Prediction   0   1
#         0 232   0
#         1   0 114

# Accuracy (Precision) : Los predichos bien/ TOTAL = 1
# Sensitivity (Recall) = TP/(TP+FN) = 1(seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) =  1(seria cuantos predije bien del total de negativos)

confusionMatrix(predict(svmLinealPfam, PfamGOTest[,-1], type = "class"), PfamGOTest[,1])
#            Reference
#Prediction  0  1
#         0 56  9
#         1  2 20

# Accuracy (Precision) : Los predichos bien/ TOTAL  = 0.8736
# Sensitivity (Recall) = TP/(TP+FN) = 0.9655 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) = 0.6897 (seria cuantos predije bien del total de negativos)

#### Optimización de hiperparámetros KERNEL LINEAL mediante validación cruzada 10-fold#####
costs <- c(0.001, 0.01, 0.1, 1 ,10)
set.seed(123)
tuningLinealPfam <- tune(svm, GO3 ~ ., data = PfamGOTrain, cross=10,
                     kernel = "linear", 
                     ranges = list(cost = costs), 
                     scale = TRUE)

summary(tuningLinealPfam)
# - best parameters: cost 10
# - best performance: 0.2080672

ggplot(data = tuningLinealPfam$performances, aes(x = cost, y = error)) +
  geom_line() +
  geom_point() +
  labs(title = "Error de validación ~ hiperparámetro C",
       x= "Costo", y= "Error de validación") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(hjust = 1,size=12),
        axis.text.y = element_text(size=12))

tuningLinealPfam$best.model

# Almacenamos el modelo optimo obtenido y accedemos a su información
best_linear_model_Pfam <- tuningLinealPfam$best.model

confusionMatrix(predict(best_linear_model_Pfam, PfamGOTrain[,-1], type = "class"), PfamGOTrain[,1])
#  Accuracy train: 1
#            Reference
#Prediction   0   1
#         0 232   0
#         1   0 114

confusionMatrix(predict(best_linear_model_Pfam, PfamGOTest[,-1], type = "class"), PfamGOTest[,1])
#            Reference
#Prediction  0  1
#         0 57  9
#         1  1 20
# Accuracy test : 0.8851  
# Sensitivity : 0.9828          
# Specificity : 0.6897

matrixErrorPfam <- matrix(0,5,6) #genero una fila de 5 filas y 6 columnas
colnames(matrixErrorPfam) <- c("Error-Train", "Especif-Train", "Sensib-Train",
                           "ErrorTest", "Especif-Test", "Sensib-Test")


for (i in 1:5){
  PfamGOTestKFold <- FeatPfamGO[indexPfamGO[[i]], ] 
  #ME QUEDO CON LA LISTA i COMO SUBSET DE DATOS DE TEST
  PfamGOTrainKFold <- FeatPfamGO[setdiff(seq(1:dim(FeatPfamGO)[1]), indexPfamGO[[i]]), ]
  #ME QUEDO CON EL RESTO COMO SUBSET DE DATOS PARA TRAIN
  fit1 <- svm(GO3 ~ ., data = PfamGOTrainKFold, kernel = "linear", cost = 10, scale = TRUE)
  #AJUSTO EL MODELO CON MIS DATOS DE TRAIN
  predictGO <- predict(fit1, PfamGOTrainKFold[,-1], type = "class")
  #PREDIGO LAS CLASES PARA LOS DATOS DE TRAIN
  matrixErrorPfam[i,1] <- (1- confusionMatrix(predictGO, PfamGOTrainKFold[,1])$overall['Accuracy'])*100
  #confusionMatrix es una funcion del paquete caret
  #le tengo q pasar las clases predichas y las clases reales.
  #CALCULO EL ERROR DE TRAINING Y LO GUARDO EN LA MATRIZ DE ERROR
  matrixErrorPfam[i,2] <- confusionMatrix(predictGO, PfamGOTrainKFold[,1])$byClass['Specificity']
  matrixErrorPfam[i,3] <- confusionMatrix(predictGO, PfamGOTrainKFold[,1])$byClass['Sensitivity']
  predictGO <- predict(fit1, PfamGOTestKFold[,-1], type = "class")
  #PREDIGO LAS CLASES CON LOS DATOS DE TEST
  matrixErrorPfam[i,4] <- (1- confusionMatrix(predictGO, PfamGOTestKFold[,1])$overall['Accuracy'])*100
  #CALCULO EL ERROR DE TEST Y LO GUARDO EN LA MATRIZ DE ERROR
  matrixErrorPfam[i,5] <- confusionMatrix(predictGO, PfamGOTestKFold[,1])$byClass['Specificity']
  matrixErrorPfam[i,6] <- confusionMatrix(predictGO, PfamGOTestKFold[,1])$byClass['Sensitivity']
}

# Error de train promedio
mean(matrixErrorPfam[,1]) #0

# Especificidad de train promedio
mean(matrixErrorPfam[,2]) #1

#Error de test promedio
mean(matrixErrorPfam[,4]) #18.01925

# Especificidad de test promedio
mean(matrixErrorPfam[,5]) #0.5310345



### Kernel Radial Pfam####
# Gamma Default = 1
svmRadialPfam <- svm(x= PfamGOTrain[,-1], y = PfamGOTrain[,1], kernel="radial", cost=1,
                 type="C-classification", scale =TRUE) 

summary(svmRadialPfam)
# Number of Support Vectors:  299 ( 185 114)

confusionMatrix(predict(svmRadialPfam, PfamGOTrain[,-1], type = "class"), PfamGOTrain[,1])
#          Reference
#Prediction   0   1
#         0 232 114
#         1   0   0

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.6705
# Sensitivity (Recall) = TP/(TP+FN) = 1 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) = 0/(115)= 0 (seria cuantos predije bien del total de negativos)
# No logro predecir ninguna como clase 1


##### Optimizamos Kerner Radial Pfam: costo y el gamma en simultaneo ####

costs <- c(0.001, 0.01, 0.1, 1 ,10, 10)
gammas <- c(0.001, 0.01, 0.1, 1 ,10, 100)

# sampling method: 10-fold cross validation 
set.seed(123)
tuningRadialPfam <- tune(svm, GO3 ~ ., data = PfamGOTrain, cross= 10,
                     kernel = "radial", 
                     ranges = list(cost = costs, 
                                   gamma = gammas), 
                     scale = TRUE)

summary(tuningRadialPfam)
# best parameters: cost= 10 gamma=0.1
# best performance: 0.222605
# Con un kernel radial, los hiperparámetros que reducen el error de 
# validación son coste = 1, gamma = 0.1.

ggplot(data = tuningRadialPfam$performances, aes(x = cost, y = error, color = factor(gamma))) +
  geom_line() +
  geom_point() +
  labs(title = "Error de validación ~ hiperparámetro C y gamma",
       x= "Costo", y= "Error de validación") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(hjust = 1,size=12),
        axis.text.y = element_text(size=12),
        legend.position = "bottom")


best_radial_model_Pfam <- tuningRadialPfam$best.model

summary(best_radial_model_Pfam)
# Number of Support Vectors:  300 ( 201 99 )

# Para los datos de train
confusionMatrix(predict(best_radial_model_Pfam, PfamGOTrain[,-1], type = "class"), PfamGOTrain[,1])
#            Reference
#Prediction   0   1
#         0 232   0
#         1   0 114
# Accuracy : 1
# Sensitivity : 1
# Specificity : 1

# Para los datos de test
confusionMatrix(predict(best_radial_model_Pfam, PfamGOTest[,-1], type = "class"), PfamGOTest[,1])
#            Reference
#Prediction  0  1
#         0 57 12
#         1  1 17
# Accuracy : 0.8506
# Sensitivity : 0.9828
# Specificity : 0.5862


# CALCULO DEL ERROR MEDIO PARA LAS 5 PARTICIONES#
matrixErrorRadialPfam <- matrix(0,5,6) #genero una fila de 5 filas y 6 columnas
colnames(matrixErrorRadialPfam) <- c("Error-Train", "Especif-Train", "Sensib-Train",
                                 "ErrorTest", "Especif-Test", "Sensib-Test")

for (i in 1:5){
  PfamGOTestKFold <- FeatPfamGO[indexPfamGO[[i]], ] 
  #ME QUEDO CON LA LISTA i COMO SUBSET DE DATOS DE TEST
  PfamGOTrainKFold <- FeatPfamGO[setdiff(seq(1:dim(FeatPfamGO)[1]), indexPfamGO[[i]]), ]
  #ME QUEDO CON EL RESTO COMO SUBSET DE DATOS PARA TRAIN
  fit1 <- svm(GO3 ~ ., data = PfamGOTrainKFold, kernel = "radial",
              cost = 10, gamma=0.1, scale = TRUE)
  #AJUSTO EL MODELO CON MIS DATOS DE TRAIN
  predictGO <- predict(fit1, PfamGOTrainKFold[,-1], type = "class")
  #PREDIGO LAS CLASES PARA LOS DATOS DE TRAIN
  matrixErrorRadialPfam[i,1] <- (1- confusionMatrix(predictGO, PfamGOTrainKFold[,1])$overall['Accuracy'])*100
  #confusionMatrix es una funcion del paquete caret
  #le tengo q pasar las clases predichas y las clases reales.
  #CALCULO EL ERROR DE TRAINING Y LO GUARDO EN LA MATRIZ DE ERROR
  matrixErrorRadialPfam[i,2] <- confusionMatrix(predictGO, PfamGOTrainKFold[,1])$byClass['Specificity']
  matrixErrorRadialPfam[i,3] <- confusionMatrix(predictGO, PfamGOTrainKFold[,1])$byClass['Sensitivity']
  predictGO <- predict(fit1, PfamGOTestKFold[,-1], type = "class")
  #PREDIGO LAS CLASES CON LOS DATOS DE TEST
  matrixErrorRadialPfam[i,4] <- (1- confusionMatrix(predictGO, PfamGOTestKFold[,1])$overall['Accuracy'])*100
  #CALCULO EL ERROR DE TEST Y LO GUARDO EN LA MATRIZ DE ERROR
  matrixErrorRadialPfam[i,5] <- confusionMatrix(predictGO, PfamGOTestKFold[,1])$byClass['Specificity']
  matrixErrorRadialPfam[i,6] <- confusionMatrix(predictGO, PfamGOTestKFold[,1])$byClass['Sensitivity']
}

# Error de train promedio
mean(matrixErrorRadialPfam[,1]) #0

# Especificidad de train promedio
mean(matrixErrorRadialPfam[,2]) #1

#Error de test promedio
mean(matrixErrorRadialPfam[,4]) #19.17402

# Especificidad de test promedio
mean(matrixErrorRadialPfam[,5]) #0.4749


#### Probamos kernel sigmoideal Pfam #####

# Costo=1 y Gamma Default
svmSigmoidPfam <- svm(x= PfamGOTrain[,-1], y = PfamGOTrain[,1], kernel="sigmoid", cost=1,
                  type="C-classification", scale =TRUE) 

summary(svmSigmoidPfam)
# Number of Support Vectors:  288 ( 174 114 )

confusionMatrix(predict(svmSigmoidPfam, PfamGOTrain[,-1], type = "class"), PfamGOTrain[,1])
#          Reference
#Prediction   0    1
#         0 232  114
#         1   0    0

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.6705
# Sensitivity (Recall) = TP/(TP+FN) =  1 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) =  0 (seria cuantos predije bien del total de negativos)


# Probamos con los datos de test
confusionMatrix(predict(svmSigmoidPfam, PfamGOTest[,-1], type = "class"), PfamGOTest[,1])
#          Reference
#Prediction   0   1
#         0  58  29
#         1   0   0

# Accuracy (Precision) : Los predichos bien/ TOTAL = 0.667
# Sensitivity (Recall) = TP/(TP+FN) =  1 (seria cuantos predije bien del total de positivos)
# Specificity = TN/(TN + FP) =  0 (seria cuantos predije bien del total de negativos)

##### Optimizamos Kernel sigmoideal Pfam: costo y el gamma en simultaneo ####

costs <- c(0.001, 0.01, 0.1, 1 ,10)
gammas <- c(0.001, 0.01, 0.1, 1 ,10, 100)

# sampling method: 10-fold cross validation 
set.seed(123)
tuningSigmoidPfam <- tune(svm, GO3 ~ ., data = PfamGOTrain, cross= 10,
                      kernel = "sigmoid", 
                      ranges = list(cost = costs, 
                                    gamma = gammas), 
                      scale = TRUE)

summary(tuningSigmoidPfam)
# best parameters: cost= 0.1 gamma=0.1
# best performance: 0.2884034 
# Con un kernel sigmoidio, los hiperparámetros que reducen el error de 
# validación son coste = 0.1, gamma = 0.1

ggplot(data = tuningSigmoidPfam$performances, aes(x = cost, y = error, color = factor(gamma))) +
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

best_sigmoid_model_Pfam <- tuningSigmoidPfam$best.model

summary(best_sigmoid_model_Pfam)
# Number of Support Vectors:  261 ( 161 100 )

# Para los datos de train
confusionMatrix(predict(best_sigmoid_model_Pfam, PfamGOTrain[,-1], type = "class"), PfamGOTrain[,1])
#            Reference
#Prediction   0   1
#         0 230  15
#         1   2  99
# Accuracy : 0.9914
# Sensitivity : 0.8684
# Specificity : 0.8684

# Para los datos de test
confusionMatrix(predict(best_sigmoid_model_Pfam, PfamGOTest[,-1], type = "class"), PfamGOTest[,1])
#            Reference
#Prediction   0  1
#         0  53 10
#         1   5 19
# Accuracy : 0.8276
# Sensitivity : 0.9138
# Specificity : 0.6552


### CALCULO DEL ERROR MEDIO PARA LAS 5 PARTICIONES###
matrixErrorSigmoidPfam <- matrix(0,5,6) #genero una fila de 5 filas y 6 columnas
colnames(matrixErrorSigmoidPfam) <- c("Error-Train", "Especif-Train", "Sensib-Train",
                                  "ErrorTest", "Especif-Test", "Sensib-Test")

for (i in 1:5){
  PfamGOTestKFold <- FeatPfamGO[indexPfamGO[[i]], ] 
  #ME QUEDO CON LA LISTA i COMO SUBSET DE DATOS DE TEST
  PfamGOTrainKFold <- FeatPfamGO[setdiff(seq(1:dim(FeatPfamGO)[1]), indexPfamGO[[i]]), ]
  #ME QUEDO CON EL RESTO COMO SUBSET DE DATOS PARA TRAIN
  fit1 <- svm(GO3 ~ ., data = PfamGOTrainKFold, kernel = "sigmoid",
              cost = 1, gamma=1, scale = TRUE)
  #AJUSTO EL MODELO CON MIS DATOS DE TRAIN
  predictGO <- predict(fit1, PfamGOTrainKFold[,-1], type = "class")
  #PREDIGO LAS CLASES PARA LOS DATOS DE TRAIN
  matrixErrorSigmoidPfam[i,1] <- (1- confusionMatrix(predictGO, PfamGOTrainKFold[,1])$overall['Accuracy'])*100
  #confusionMatrix es una funcion del paquete caret
  #le tengo q pasar las clases predichas y las clases reales.
  #CALCULO EL ERROR DE TRAINING Y LO GUARDO EN LA MATRIZ DE ERROR
  matrixErrorSigmoidPfam[i,2] <- confusionMatrix(predictGO, PfamGOTrainKFold[,1])$byClass['Specificity']
  matrixErrorSigmoidPfam[i,3] <- confusionMatrix(predictGO, PfamGOTrainKFold[,1])$byClass['Sensitivity']
  predictGO <- predict(fit1, PfamGOTestKFold[,-1], type = "class")
  #PREDIGO LAS CLASES CON LOS DATOS DE TEST
  matrixErrorSigmoidPfam[i,4] <- (1- confusionMatrix(predictGO, PfamGOTestKFold[,1])$overall['Accuracy'])*100
  #CALCULO EL ERROR DE TEST Y LO GUARDO EN LA MATRIZ DE ERROR
  matrixErrorSigmoidPfam[i,5] <- confusionMatrix(predictGO, PfamGOTestKFold[,1])$byClass['Specificity']
  matrixErrorSigmoidPfam[i,6] <- confusionMatrix(predictGO, PfamGOTestKFold[,1])$byClass['Sensitivity']
}

# Error de train promedio
mean(matrixErrorSigmoidPfam[,1]) #5.311089

# Especificidad de train promedio
mean(matrixErrorSigmoidPfam[,2]) #0.8828833

#Error de test promedio
mean(matrixErrorSigmoidPfam[,4]) #19.40658

# Especificidad de test promedio
mean(matrixErrorSigmoidPfam[,5]) #0.5167488


##### TABLA CON ERRORES PROMEDIO Pfam SVM DISTINTOS KERNEL PARAMETROS OPTIMIZADOS ####
ResumenErroresPfam <- matrix(data=NA, nrow = 3, ncol= 6, 
                         dimnames = list(c("Lineal", "Radial", "Sigmoideo"), 
                                         c("ErrorTrain-promedio", "Especificidad-Train-Promedio",
                                           "Sensibilidad-Train-Promedio",
                                           "ErrorTest-promedio", "Especificidad-Test-Promedio",
                                           "Sensibilidad-Test-Promedio")))

ResumenErroresPfam["Lineal","ErrorTrain-promedio"] <- mean(matrixErrorPfam[,1])
ResumenErroresPfam["Lineal","Especificidad-Train-Promedio"] <- mean(matrixErrorPfam[,2])
ResumenErroresPfam["Lineal","Sensibilidad-Train-Promedio"] <- mean(matrixErrorPfam[,3])
ResumenErroresPfam["Lineal","ErrorTest-promedio"] <- mean(matrixErrorPfam[,4])
ResumenErroresPfam["Lineal","Especificidad-Test-Promedio"] <- mean(matrixErrorPfam[,5])
ResumenErroresPfam["Lineal","Sensibilidad-Test-Promedio"] <- mean(matrixErrorPfam[,6])

ResumenErroresPfam["Radial","ErrorTrain-promedio"] <- mean(matrixErrorRadialPfam[,1])
ResumenErroresPfam["Radial","Especificidad-Train-Promedio"] <- mean(matrixErrorRadialPfam[,2])
ResumenErroresPfam["Radial","Sensibilidad-Train-Promedio"] <- mean(matrixErrorRadialPfam[,3])
ResumenErroresPfam["Radial","ErrorTest-promedio"] <- mean(matrixErrorRadialPfam[,4])
ResumenErroresPfam["Radial","Especificidad-Test-Promedio"] <- mean(matrixErrorRadialPfam[,5])
ResumenErroresPfam["Radial","Sensibilidad-Test-Promedio"] <- mean(matrixErrorRadialPfam[,6])

ResumenErroresPfam["Sigmoideo","ErrorTrain-promedio"] <- mean(matrixErrorSigmoidPfam[,1])
ResumenErroresPfam["Sigmoideo","Especificidad-Train-Promedio"] <- mean(matrixErrorSigmoidPfam[,2])
ResumenErroresPfam["Sigmoideo","Sensibilidad-Train-Promedio"] <- mean(matrixErrorSigmoidPfam[,3])
ResumenErroresPfam["Sigmoideo","ErrorTest-promedio"] <- mean(matrixErrorSigmoidPfam[,4])
ResumenErroresPfam["Sigmoideo","Especificidad-Test-Promedio"] <- mean(matrixErrorSigmoidPfam[,5])
ResumenErroresPfam["Sigmoideo","Sensibilidad-Test-Promedio"] <- mean(matrixErrorSigmoidPfam[,6])

ResumenErroresPfam  





##### Ajustamos costos y gammas en simultaneo de forma manual#####
costs <- c(0.001, 0.01, 0.1, 1 ,10, 100)
gammas <- c(0.001, 0.01, 0.1, 1 ,10, 100)

#'ajustarModelosPfam
#' Recibe (parametro1= costo, parametro2=gamma) y ajusta un modelo de svm con kernel
#' radial y esos parámetros a los datos PfamGOTrain.
#' Devuelve un vector numérico, con las métricas precisión, sensibilidad y specificidad
#' del modelo ajustado, para el set de datos de test PfamGOTest.

ajustarModelosPfam <- function(parametro1, parametro2){
  svmRadial <- svm(x= PfamGOTrain[,-1], y = PfamGOTrain[,1], kernel="radial",
                   gamma = parametro2, 
                   cost= parametro1, type="C-classification", scale =TRUE)
  return (c(confusionMatrix(predict(svmRadial, PfamGOTest[,-1], type = "class"), PfamGOTest[,1])$overall['Accuracy'],
            confusionMatrix(predict(svmRadial, PfamGOTest[,-1], type = "class"), PfamGOTest[,1])$byClass['Sensitivity'],
            confusionMatrix(predict(svmRadial, PfamGOTest[,-1], type = "class"), PfamGOTest[,1])$byClass['Specificity']))
}

metricasPfam <- matrix(data=NA, nrow = 36, ncol= 5, 
                   dimnames = list(c(), 
                                   c("Costo", "Gamma", "Acuraccy", "Sensib", "Specif")))

k = 1
for (i in costs){
  for (j in gammas){
    metricasPfam[k,"Costo"] <- i
    metricasPfam[k,"Gamma"] <- j
    metricasPfam[k,"Acuraccy"] <- ajustarModelosPfam(i, j)[1]
    metricasPfam[k,"Sensib"] <- ajustarModelosPfam(i, j)[2]
    metricasPfam[k,"Specif"] <- ajustarModelosPfam(i, j)[3]
    k = k+1
  }
}

metricasPfam <- as.data.frame(metricasPfam)
max(metricasPfam$Acuraccy) # Mayor Acuraccy = 0.7906977
max(metricasPfam$Specif) # Mayor Acuraccy = 0.7816092
which(metricasPfam$Acuraccy == max(metricasPfam$Acuraccy)) #32 (indice de la fila con mayor Acc)
metricasPfam[32,]
#   Costo Gamma  Acuraccy    Sensib    Specif
#32   100  0.01 0.7906977 0.9482759 0.4642857
# Vemos que el mejor Acuraccy lo obtenemos con cost = 1, gamma = 0.1

# Ploteamos Acuraccy vs Gamma, para Costo = 100
acuraccyPfam <- c()
for (i in gammas){
  acuraccyPfam <- c(acuraccyPfam, ajustarModelosPfam(parametro1=100, parametro2=i)[1])
}

plot(log10(gammas), acuraccyPfam)
acuraccyPfam[2] #0.7906977
# Vemos que el mejor gamma es 0.01

svmRadialBest <- svm(x= PfamGOTrain[,-1], y = PfamGOTrain[,1], kernel="radial", 
                     gamma = 0.01, cost= 100, type="C-classification", scale =TRUE)
confusionMatrix(predict(svmRadialBest, PfamGOTest[,-1], type = "class"), PfamGOTest[,1])
#'            Reference
#' Prediction  0  1
#'          0 55 15
#'          1  3 13
