# BALANCEAMOS LOS DATOS

# Cargamos los datos

setwd("D:/OneDrive - fbioyf.unr.edu.ar/Especializ_Bioinf/7-Procesamiento_inteligente_de_datos/Trabajo-Final")

IDS <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/id_tp.txt",
                as.is = TRUE, sep = "\t", header= FALSE)#Lista de IDS

IDS <- IDS$V1 # Creamos un array con los IDS

Clases <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/clases.txt",
                   as.is = TRUE, sep = "\t", header= TRUE)#Matriz de clases

rownames(Clases) <- Clases$ID#convierto la columna ID en el indice

Clases <- Clases[IDS,] #Nos quedamos solo con los IDS que tenían Pfam

FeaturesPept <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/FeaturesPeptides.txt",
                         as.is = TRUE, sep = "\t", header= TRUE)#Tabla Features propiedades FQ

rownames(FeaturesPept) <- FeaturesPept$ID#convierto la columna ID en el indice

FeaturesPfam <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/FeaturesPfam.txt",
                         as.is = TRUE, sep = "\t", header= TRUE)#Tabla Features Pfam

rownames(FeaturesPfam) <- FeaturesPfam$ID#convierto la columna ID en el indice

##### BALANCEO DE DATOS #####
# Para comenzar, dejamos los 145 datos de la clase GO3, 
# y dejamos solo 290 negativos al azar de las otras clases.

IDS_GO3 <- IDS[Clases$GO3 == 1]#me quedo con los ID de la clase GO3 (+)
IDS_notGO3 <- IDS[Clases$GO3 != 1]#me quedo con los ID que no son de la clase GO3 (-)

set.seed(123)#fijo una semilla
IDS_notGO3_reducidos <- sample(IDS_notGO3, 290, replace=FALSE)#sampleo 290 datos de los (-)
IDS_balanceados <- c(IDS_GO3, IDS_notGO3_reducidos)#uno los (-) con los que me quede con los (+)

Clases_Red <- Clases[IDS_balanceados,] # Nos quedamos solo con los IDS que tenían Pfam y las clases de los nuevos datos balanceados
FeaturesPept_Red <- FeaturesPept[IDS_balanceados,]#reduzco la tabla de FeaturesPep con los nuevos datos
FeaturesPfam_Red <- FeaturesPfam[IDS_balanceados,]#reduzco la tabla de FeaturesPfam con los nuevos datos


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

library(caret)#para generar las particiones

set.seed(123)

indexPeptGO <- createFolds(t(FeatPeptGO[, "GO3"]), k = 5)
PeptGOTest <- FeatPeptGO[indexPeptGO[[3]], ]
PeptGOTrain <- FeatPeptGO[setdiff(seq(1:dim(FeatPeptGO)[1]), indexPeptGO[[3]]), ]

dim(PeptGOTrain) #346 observaciones en el set de train
dim(PeptGOTest) #87 observaciones en el set de test

indexPfamGO <- createFolds(t(FeatPfamGO[, "GO3"]), k = 5)
PfamGOTest <- FeatPfamGO[indexPfamGO[[3]], ]
PfamGOTrain <- FeatPfamGO[setdiff(seq(1:dim(FeatPfamGO)[1]), indexPfamGO[[3]]), ]

######----------------------Redes Neuronales Pfam Datos Balanceados-----------------###########
library(neuralnet)

### Arquitectura red: 2 capas ocultas, 1 con 3 neuronas y la otra con 2 ####
Pfam_GO3ANN_red <- neuralnet(GO3 == 1 ~ ., data = PfamGOTrain, hidden = c(3,2))
summary(Pfam_GO3ANN_red)

# Para los datos de train
predictPfam_GO3ANN.Train_red <- predict(Pfam_GO3ANN_red, PfamGOTrain)
View(predictPfam_GO3ANN.Train_red)

table(Realidad = PfamGOTrain$GO3 == 1, Predicho = predictPfam_GO3ANN.Train_red[, 1] > 0.5)
#         Predicho
#Realidad FALSE(0) TRUE(1)
#FALSE(0)   232       0      ; De los negativos, uno solo clasifico mal
#TRUE(1)      0      114        ;  Los positivos los clasifica todos bien 


# Para los datos de Test
predictPfam_GO3ANN_Test_red<- predict(Pfam_GO3ANN_red, PfamGOTest)
table(Realidad = PfamGOTest$GO3 == 1, Predicho = predictPfam_GO3ANN_Test_red[, 1] > 0.5)
#         Predicho
#Realidad FALSE TRUE
#FALSE   40     18
#TRUE     5     24

1-((40+24)/(40+24+5+18)) # 0.2644 Error Test
24/29 # 0.827586 Especificidad 
40/58 # 0.68965 Sensitivity


### Arquitectura red: 3 capas ocultas, (4,3,2) ####
Pfam_GO3ANN_red2 <- neuralnet(GO3 == 1 ~ ., data = PfamGOTrain, hidden = c(4,3,2))
summary(Pfam_GO3ANN_red2)

# Para los datos de train
predictPfam_GO3ANN.Train_red2 <- predict(Pfam_GO3ANN_red2, PfamGOTrain)
#View(predictPfam_GO3ANN.Train_red2)

table(Realidad = PfamGOTrain$GO3 == 1, Predicho = predictPfam_GO3ANN.Train_red2[, 1] > 0.5)
#         Predicho
#Realidad FALSE(0) TRUE(1)
#FALSE(0)   232       0      ; De los negativos, uno solo clasifico mal
#TRUE(1)      0      114        ;  Los positivos los clasifica todos bien 


# Para los datos de Test
predictPfam_GO3ANN_Test_red2<- predict(Pfam_GO3ANN_red2, PfamGOTest)
table(Realidad = PfamGOTest$GO3 == 1, Predicho = predictPfam_GO3ANN_Test_red2[, 1] > 0.5)
#         Predicho
#Realidad FALSE TRUE
#FALSE    34     24
#TRUE      4     25

1-((34+25)/(34+4+24+25)) # 0.3218 Error Test
25/29 # 0.862069 Especificidad 
34/58 # 0.5862 Sensitivity

