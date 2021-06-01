#setwd("D:/OneDrive - fbioyf.unr.edu.ar/Especializ_Bioinf/7-Procesamiento_inteligente_de_datos/Trabajo-Final")

# ---------- CREATE OF CLASES MATRIX -----------------------

GO1 <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/GO0016570_ID_GO_Hs_histone%20modification%20(%2B).tsv", as.is = TRUE, sep = "\t", header= FALSE)
GO2 <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/GO0044249_ID_GO_At_cellular%20biosynthetic%20process(%2B).tsv", as.is = TRUE, sep = "\t", header= FALSE)  
GO3 <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/GO0044255_ID_GO_HS_cellular%20lipid%20metabolic%20process(%2B).tsv", as.is = TRUE, sep = "\t", header= FALSE)  
GO4 <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/GO0044257_ID_GO_At_cellular%20protein%20catabolic%20process(%2B).tsv", as.is = TRUE, sep = "\t", header= FALSE)  
GO5 <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/GO0002682_ID_GO_Hs_regulation%20of%20immune%20system%20process(-).tsv",as.is = TRUE, sep = "\t", header= FALSE)
GO6 <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/GO0009737_ID_GO_At_responce%20to%20abscisic%20acid(-).tsv", as.is = TRUE, sep = "\t", header= FALSE)
GO7 <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/GO0015979_ID_GO_At_photosynthesis%20(-).tsv", as.is = TRUE, sep = "\t", header= FALSE)
GO8 <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/GO0071555_ID_GO_At_cell%20wall%20organization%20(-).tsv", as.is = TRUE, sep = "\t", header= FALSE)

#LE pongo as.is = TRUE, para que no me convierta los caracteres a factor.
#Esto es importante, xq si tengo atributos que son caracteres, me los va a convertir a factor.
#Yo lo unico que quiero es que la columna donde este la clase sea un factor (y esa la convierto a factor
#mas adelante)

listaGO = list(GO1, GO2, GO3, GO4, GO5, GO6, GO7, GO8)

IDs = vector(mode = "list", length = 8)

for(i in 1:length(listaGO)){
  listaGO[[i]] <- listaGO[[i]][-1,]
  # Nos quedamos solo con 145 IDs de cada GO
  IDs[[i]] <- unique(listaGO[[i]][,2])[1:145]
}


data = vector(mode = "list", length = 8)

for(i in 1:length(data)){
  data[[i]] <- data.frame(ID = IDs_prueba[[i]], 
                          GO1 = if(i==1) {1} else {as.integer(IDs[[i]] %in% listaGO[[1]]$V2)},
                          GO2 = if(i==2) {1} else {as.integer(IDs[[i]] %in% listaGO[[2]]$V2)}, 
                          GO3 = if(i==3) {1} else {as.integer(IDs[[i]] %in% listaGO[[3]]$V2)}, 
                          GO4 = if(i==4) {1} else {as.integer(IDs[[i]] %in% listaGO[[4]]$V2)}, 
                          GO5 = if(i==5) {1} else {as.integer(IDs[[i]] %in% listaGO[[5]]$V2)}, 
                          GO6 = if(i==6) {1} else {as.integer(IDs[[i]] %in% listaGO[[6]]$V2)}, 
                          GO7 = if(i==7) {1} else {as.integer(IDs[[i]] %in% listaGO[[7]]$V2)}, 
                          GO8 = if(i==8) {1} else {as.integer(IDs[[i]] %in% listaGO[[8]]$V2)},
                          Tetrahymena = if(i%in%c(1:4)) {1} else{0})
}


tabla_final = rbind(data[[1]], data[[2]], data[[3]], data[[4]], data[[5]], data[[6]], data[[7]],data[[8]])

# Eliminamos ids que tienen - o : en su nombre.
tabla_final = tabla_final[- c(grep("-", tabla_final$ID),grep(":", tabla_final$ID)),]

#Eliminamos IDs duplicados
tabla_final$ID <- as.character(tabla_final$ID)
tabla_final <- tabla_final[!duplicated(tabla_final$ID),]

rownames(tabla_final) <- tabla_final$ID

# Tenemos 574 IDs que son de Tetrahymena
sum(tabla_final$Tetrahymena)


#-------------DOWNLOAD OF FASTA SEQUENCES ------------------

#Descargar secuencias de AA de uniprot

install.packages("UniprotR")
library(UniprotR)

# '''
# #me pide Biostring:
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("Biostrings")
# 
# #me pide GenomicaAlignments:
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("GenomicAlignments")'''

#Retrieve the names/identifiers for the sequences
#UniprotNames3 <- GetNamesTaxa(tabla_final$ID)

#Retrieve the sequences
UniprotSeqs <- GetSequences(tabla_final$ID)

#the row names are the accesion number. i want them as a column
UniprotSeqs$Acc_No <- row.names(UniprotSeqs)

#Keep only the Acc_No and sequences from UniprotSeqs
seqs <- c("Acc_No", "Sequence")
UniprotSeq2 <- UniprotSeqs[seqs]

UniprotSeq2$Sequence<-as.character(UniprotSeq2$Sequence)

write.table(UniprotSeq2, file = 'SeqFASTA_tp.txt', append = FALSE, sep="\t", dec=".", 
            row.names = FALSE, col.names=TRUE)


#------------- RETRIEVE OF FEATURES --------------------------
install.packages("Peptides", dependencies=TRUE)
library(Peptides)

help('Peptides')

aaList() #Da la lista de 20 AA

# Las sequencias 5, 269, 587 y 629 tiran error
# Por ejemplo, la secuencia 5 tiene una X de AA
# Decidimos eliminar esas proteinas
UniprotSeq2 <- UniprotSeq2[-c(5,269,587,629),]

ID2 <- UniprotSeq2[,1]
FASTAS <- UniprotSeq2[,2]

### Proteins length ###
LongProteinas<- lengthpep(FASTAS)

### Molecular Weight in Dalton ###
MWProteinas <- mw(FASTAS)
'?'(mw) # Para ver la ayuda de esa funcion

### Amino acid composition ###
CompAAProteinas <- aaComp(seq = FASTAS)

#Devuelve una matriz con el numero de cada clase de AA en una columna, y en la otra columna,
# el % del total de AAs que representa esa clase de AA. Ejemplo:

128/LongProteinas[1]*100 #% del total que representan los AA Tiny


#Creamos un DataFrame vacio con la clase de cada AA como columna
CompAAFinal <- data.frame(matrix(nrow = length(ID2), ncol=10, 
                                      dimnames = list(c(), 
                                                      c("ID", "TinyAA", 
                                                        "SmallAA", 
                                                        "AliphaticAA", 
                                                        "AromaticAA",
                                                        "NonPolarAA",
                                                        "PolarAA",
                                                        "ChargedAA",
                                                        "BasicAA",
                                                        "AcidicAA"))))

#Asignamos los ID
CompAAFinal$ID = ID2
row.names(CompAAFinal) <- CompAAFinal$ID


for (i in 1:nrow(CompAAFinal)){
  CompAAFinal$TinyAA[i] <- as.data.frame(CompAAProteinas[[i]])$'Mole%'[1]
  CompAAFinal$SmallAA[i] <- as.data.frame(CompAAProteinas[[i]])$'Mole%'[2]
  CompAAFinal$AliphaticAA[i] <- as.data.frame(CompAAProteinas[[i]])$'Mole%'[3]
  CompAAFinal$AromaticAA[i] <- as.data.frame(CompAAProteinas[[i]])$'Mole%'[4]
  CompAAFinal$NonPolarAA[i] <- as.data.frame(CompAAProteinas[[i]])$'Mole%'[5]
  CompAAFinal$PolarAA[i] <- as.data.frame(CompAAProteinas[[i]])$'Mole%'[6]
  CompAAFinal$ChargedAA[i] <- as.data.frame(CompAAProteinas[[i]])$'Mole%'[7]
  CompAAFinal$BasicAA[i] <- as.data.frame(CompAAProteinas[[i]])$'Mole%'[8]
  CompAAFinal$AcidicAA[i] <- as.data.frame(CompAAProteinas[[i]])$'Mole%'[9]
}

### Net Charge ###
CargaNetaProteinas <- charge(FASTAS, pH = 7, pKscale = "EMBOSS")
# The net charge of a protein can be calculated specifying the pH value and 
# one of the nine pKa scales availables 
#(Bjellqvist, Dawson, EMBOSS, Lehninger, Murray, Rodwell, Sillero, Solomon or Stryer).


### Isoelectric point ###
PuntoIsoelecProteinas <- pI(FASTAS, pKscale = "EMBOSS")
# the isoelectric point of a peptide may be performed through the function pI
#specifying one of the nine pKa scales available 
#(Bjellqvist, Dawson, EMBOSS, Lehninger, Murray, Rodwell, Sillero, Solomon or Stryer).

### Aliphatic index ###
AliphaticIndexProteinas<- aIndex(FASTAS)
# El indice alifatico se calcula como una combinacion lineal de la composicion 
# de AA alifaticos, A, L, I y V... Puede ser redundante con la composicion de AA...
# VER SI LO DEJAMOS

### Instability index ###
# This index predicts the stability of a protein based on its amino acid composition.
InstabilityIndexProteinas <- instaIndex(FASTAS)
# Peptids are considered stable with index values less than 40.


### Boman index ###
# Diferencia si una proteina interacciona con otra proteina o con la membrana
# This function predicts the potential peptide interaction with another protein.
# It is calculated by adding each amino acid solubilities divided by the sequence length. 
# Por ejemplo, antimicrobial peptides tend to not interact with other proteins (the
# proposed mechanism of action is based on the interaction with membranes), so the values for the
# Boman index are usually negative or nearby to 0.
BomanProteinas <- boman(seq = FASTAS)

### Hydrophobicity index ###
#The hydrophobicity is an important stabilization force in protein folding;
# The hydrophobicity index is calculated adding the hydrophobicity of individual 
#amino acids and dividing this value by the length of the sequence.
# Highly expected transmembrane peptides generally have higher hydrophobicity values than 0.5 using
# Eisenberg scale.
#Peptides includes thirty-eight scales of hydrophobicity (Aboderin, AbrahamLeo, Argos, BlackMould,
#BullBreese, Casari, Chothia, Cid, Cowan3.4, Cowan7.5, Eisenberg, Engelman, Fasman, Fauchere,
# Goldsack, Guy, HoppWoods, Janin, Jones, Juretic, Kidera, Kuhn, KyteDoolittle, Levitt, Manavalan,
# Miyazawa, Parker, Ponnuswamy, Prabhakaran, Rao, Rose, Roseman, Sweet, Tanford, Welling, Wilson,
#Wolfenden or Zimmerman)
HidrofIndexProteinas <- hydrophobicity(FASTAS, scale = "Eisenberg")


### Hydrophobic moment index ### --> NO ENTENDEMOS BIEN Q CARAJO ES --> VER SI LO DEJAMOS
# HAY QUE PASARLE COMO PARAMETRO EL VALOR DE UN ANGULO, Q NI IDEA
# The hydrophobic moment is computed using the standardized Eisenberg (1984) scale, 
# windows (fragment of sequence) of eleven AA and specifying the 
# rotational angle at which it should be calculated.
HidrofMomentProteinas <- hmoment(FASTAS, angle = 100, window = 11)
# Dejamos angle = 100 y window = 11 xq es lo q figura en el paper de ejemplo
# Suggested: a-helix = 100, b-sheet=160
'?'(hmoment)


### Membrane position ###
#Eisenberg et al. (1982) found a correlation between hydrophobicity and hydrophobic moment that
#defines the protein section as globular, transmembrane or superficial.

'?'(membpos)

# ''' EJEMPLO
# data(pepdata)
# membpos(seq="DAEFRHDSGYEVHHQKLVFFAEDVGSNK", angle=100) #100% globular
# membpos(seq="SLDRSSCFTGSLDSIRAQSGLGCNSFRY", angle= 100) ##100% globular
# 
# membpos(seq=pepdata$sequence[8], angle= 100)
# 
# 3/15*100 #el % de cada tipo 
# 1/15*100
# 11/15*100
# '''

# Calculamos este parametro para todas nuestras proteinas con el angulo de 100
MembranePosProteinas <- membpos(FASTAS, angle = 100)

#Creamos un Data Frame vacio 
MembranePosFinal <- data.frame(matrix(nrow = length(ID2), ncol=4, 
                                      dimnames = list(c(), c("ID", "Globular", "Surface", "Transmembrane"))))

#Asignamos los ID
MembranePosFinal$ID = ID2
row.names(MembranePosFinal) <- MembranePosFinal$ID

#Calculamos el % de "Globular"/"Surface"/"Transmembrane" para cada proteina y lo ponemos
# en una columna para cada ID
for (i in 1:nrow(MembranePosFinal)){
  MembranePosFinal$Globular[i] <- sum(MembranePosProteinas[[i]][4] =="Globular")/nrow(MembranePosProteinas[[i]])
  MembranePosFinal$Surface[i] <- sum(MembranePosProteinas[[i]][4] =="Surface")/nrow(MembranePosProteinas[[i]])
  MembranePosFinal$Transmembrane[i] <- sum(MembranePosProteinas[[i]][4] =="Transmembrane")/nrow(MembranePosProteinas[[i]])
}



## DataFrame FINAL FINAL

Features <- data.frame(ID = ID2, 
                       Length = LongProteinas,
                       MW = MWProteinas,
                       Charge = CargaNetaProteinas,
                       pI = PuntoIsoelecProteinas,
                       aIndex = AliphaticIndexProteinas,
                       InstabilityIndex = InstabilityIndexProteinas,
                       Boman = BomanProteinas,
                       HidrophobicityIndex = HidrofIndexProteinas,
                       HMoment = HidrofMomentProteinas)

Features<- cbind(Features, CompAAFinal[,-1], MembranePosFinal[,-1])

write.table(Features, file = 'Features.txt', append = FALSE, sep="\t", dec=".", 
            row.names = FALSE, col.names=TRUE)

#------------- Creamos archivo MULTIFASTA para hacer Pfam-------

MultiFasta<-"D:/OneDrive - fbioyf.unr.edu.ar/Especializ_Bioinf/7-Procesamiento_inteligente_de_datos/Trabajo-Joaquin/MultiFASTA.fasta"

for (i in 1:nrow(UniprotSeq2)){
  cat(paste(">", UniprotSeq2[i,1],sep = " "),"\n", file= MultiFasta, append=TRUE)
  cat(UniprotSeq2[i,2],"\n",file=MultiFasta, append=TRUE)
}

# With this file, execute the following UNIX commands in terminal

# wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz 
#Descarga base de datos Pfam

# hmmpress Pfam-A.hmm #comprime base de datos para ejecutar comando hmmscan

#hmmscan --tblout PfamTP-DataMining Pfam-A.hmm MultiFASTA.fasta
# Crea el archivo MultiFASTA.fasta, que es sencillo de parsear

# Ahora hacemos unas transformaciones para q el resultado quede en columnas
# sed 's/ \+ /\t/g' PfamTP-DataMining | sed 's/ /\t/g' | cut -f 2,4 > PfamTP-DataMining-formatted.txt

# sed 's/ \+ /\t/g' PfamTP-DataMining
# siempre que haya mas de un espacio, lo reemplazamos por una tabulacion
# sed 's/ /\t/g'
# este segundo sed lo hacemos xq quedan algunas columnas separadas solo x un espacio
# entonces reemplazamos estos espacios x una tabulacion
# cut -f 2,4
#Nos quedamos solo con la fila 2 y 4, donde estan los Pfam y los ID respectivamente


# ------------- Creamos un data frame con los PFAM para cada proteina ------

Pfam <- read.csv("https://raw.githubusercontent.com/ferbracalente/ML/main/Data/PfamTP-DataMining-formatted.txt",
                 as.is = TRUE, sep = "\t", header= FALSE, skip = 3)
Pfam <- Pfam[-c(5336:nrow(Pfam)),]

length(unique(Pfam$V2)) #Vemos que hay algunos IDs que no se encontró Pfam

sum(! ID2 %in% unique(Pfam$V2)) # 19 IDs sin Pfam

IDSinPfam <- which(! ID2 %in% unique(Pfam$V2)) #Para obtener los indices de los ID2 q no tienen Pfam

ID3 <- unique(Pfam$V2) #Armo un nuevo vector con los IDs definitivos

PFAMs <- unique(Pfam$V1) #Armo un vector con los Pfam unicos

# Armamos un data frame con todos 0
# nrow = cantidad de ID; ncol = cantidad de Pfam unicos + 1 (columna para el ID)
TablaPfam <- data.frame(matrix(0, nrow = length(ID3), ncol = (1+length(unique(Pfam$V1))),
                               dimnames = list(ID3, c("ID", PFAMs))))

TablaPfam$ID <- ID3

# PfamsParaUnID : String -> DataFrame
# Recibe un ID y filtra el dataframe Pfam por ese ID
PfamsParaUnID <- function(ID){
  return(subset(Pfam, V2==ID))
}

for (i in 1:length(ID3)){
  # j itera en la cantidad de Pfam que hay para cada compuesto ID3[i]
  for (j in 1:nrow(PfamsParaUnID(ID3[i]))){
    # En la TablaPfam, fila = ID3[i], columna = Pfam[j], le asignamos un 1
    TablaPfam[ID3[i], PfamsParaUnID(ID3[i])[j,1]] = 1
  }
}


#------------ Guardamos los archivos------------------------

# Guardamos un archivo con los IDs que no tuvieron problemas ni con Peptides ni con Pfam
write.table(ID3, file = 'id_tp.txt', append = FALSE, sep=" ", dec=".", 
            row.names = FALSE, col.names=FALSE)

# Tabla con Features extraidos con Peptides

Features <- Features[-IDSinPfam,] #Le quito estos IDs a la tabla de Features

write.table(Features, file = 'FeaturesPeptides.txt', append = FALSE, sep="\t", dec=".", 
            row.names = FALSE, col.names=TRUE)


#Tabla con Features extraidos con Pfam

write.table(TablaPfam, file = 'FeaturesPfam.txt', append = FALSE, sep="\t", dec=".", 
            row.names = FALSE, col.names=TRUE)


# Tabla con las CLASES

#Dejamos solo los IDs que no tuvieron problemas ni con Peptides ni con Pfam
tabla_final <- tabla_final[ID3,] 
write.table(tabla_final, file = 'clases.txt', append = FALSE, sep="\t", dec=".", 
            row.names = FALSE, col.names=TRUE)


# Guardamos un archivo con FeaturesPeptides y clases
FeaturesAndClases <- merge(Features, tabla_final, by = 'ID')

write.table(FeaturesAndClases, file = 'FeaturesPeptidesAndClases.txt', append = FALSE, sep="\t", dec=".", 
            row.names = FALSE, col.names=TRUE)

# Guardamos un archivo con FeaturesPeptides y clases
FeaturesPfamAndClases <- merge(TablaPfam, tabla_final, by = 'ID')

write.table(FeaturesPfamAndClases, file = 'FeaturesPfamAndClases.txt', append = FALSE, sep="\t", dec=".", 
            row.names = FALSE, col.names=TRUE)

