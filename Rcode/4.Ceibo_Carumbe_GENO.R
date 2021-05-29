#####################################################################
# Datos para la construccion del mapa geneticos de Cebada x Carumbe #
# datos enviados por Ariel Castro                                   #
# de:	vontruch@fagro.edu.uy
# para:	gastonquero@gmail.com
# Cc:	lviega@fagro.edu.uy
# fecha:	28 de junio de 2018, 15:26
# asunto:	Marcadores Ceibo carumbe
# enviado por:	fagro.edu.uy
# Gaston Quero
# 19/7/2018 
############################################################


# fijar Directorio
getwd()
setwd("E:/Proyecto_Cebada")


# Cargar Paquetes
library("lmem.gwaser")
library("lmem.qtler")
library("qtl")
qtlversion()
library("ggplot2")
library("FactoMineR")
library("dplyr")
library("car")
library (ggjoy)
#library (hrbrthemes)
library(tidyverse)
#library(forcats)
library("viridis")
library("ASMap")

# Cargar datos como data frame
# Estos son los datos originales 

Ceibo.Carumbe.G1 <- read.cross (format="csv",
           dir="./genetic_map/Data/rawdata",
           file="CeiboCarumbe_1.csv", 
           na.strings= "NA",sep=";",dec=",",
           genotypes=c("AA","BB"),
           #alleles=c("A","B"),
           crosstype ="riself",
           estimate.map= FALSE, convertXdata=TRUE, error.prob=0.0001)

summary (Ceibo.Carumbe.G1)

nind   (Ceibo.Carumbe.G1)
nchr   (Ceibo.Carumbe.G1)
totmar (Ceibo.Carumbe.G1)
nmar   (Ceibo.Carumbe.G1)
nphe   (Ceibo.Carumbe.G1)

geno.image  (Ceibo.Carumbe.G1, main="CeiboxCarumbe_original")
plotMissing (Ceibo.Carumbe.G1, main="CeiboxCarumbe_original")

###### vamos a sacar los misssing data ##########
par(mfrow=c(1,2), las=1)
plot(ntyped(Ceibo.Carumbe.G1), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(Ceibo.Carumbe.G1, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")


# Eliminar marcadores con mas del 50 % de datos faltantes
n.missing <- nmissing (Ceibo.Carumbe.G1, what="mar")[(nmissing(Ceibo.Carumbe.G1, what="mar"))/sum(summary(Ceibo.Carumbe.G1)$n.ind)> 0.50]  
names.marker <- c (names(n.missing))
length(names.marker)

Ceibo.Carumbe.G1a <- drop.markers (Ceibo.Carumbe.G1, names.marker)


plotMissing (Ceibo.Carumbe.G1a, main="Ceibo.Carumbe.G1a")
summary (Ceibo.Carumbe.G1a)
nmar(Ceibo.Carumbe.G1a)
totmar(Ceibo.Carumbe.G1a)

#par(mfrow=c(1,2), las=1)
#plot(ntyped(DM.SO.cross.1), ylab="No. typed markers", main="No. genotypes by individual")
#plot(ntyped(DM.SO.cross.1, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")

# eliminar los individuos con mas de 50 no genotipado
Ceibo.Carumbe.G2 <- subset(Ceibo.Carumbe.G1a, 
                        ind=(ntyped(Ceibo.Carumbe.G1a) > ((totmar(Ceibo.Carumbe.G1a) * 70)/100)))

indiv <- subset(Ceibo.Carumbe.G1a, 
                ind=(ntyped(Ceibo.Carumbe.G1a)  < ((totmar(Ceibo.Carumbe.G1a) * 70)/100)))

indiv$pheno$id

plotMissing (Ceibo.Carumbe.G2, main="Ceibo.Carumbe.G2")
geno.image  (Ceibo.Carumbe.G2, main="CeiboxCarumbe_filtrada")

nind(Ceibo.Carumbe.G2)
par(mfrow=c(1,2), las=1)
plot(ntyped(Ceibo.Carumbe.G2), ylab="No. typed markers", main="No. genotypes by individual")
plot(ntyped(Ceibo.Carumbe.G2, "mar"), ylab="No. typed individuals",main="No. genotypes by marker")

nmar(Ceibo.Carumbe.G2)
totmar(Ceibo.Carumbe.G2)


###### voy a comparar los genotipos 
cg <- comparegeno(Ceibo.Carumbe.G2)
hist(cg[lower.tri(cg)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg[lower.tri(cg)])

wh <- which(cg > 0.9, arr=TRUE)
wh <- wh[wh[,1] < wh[,2],]
wh

########## distorcion de segregacion 
gt <- geno.table (Ceibo.Carumbe.G2)
head(gt)
x <- gt[gt$P.value < 0.05/totmar(Ceibo.Carumbe.G2),]
class(x)
dim(x)

### calculo de las frecuencias 
total <- nind (Ceibo.Carumbe.G2)
gt <- geno.table(Ceibo.Carumbe.G2)
missing <- gt[,2]
num.obs <- c(total - missing)
obs.alelo <- gt[, c(3,4,5)]
Geno.freq <- as.vector(obs.alelo/num.obs)
head(Geno.freq)
head(obs.alelo)

##### No polimorficos para AA
no.polimorficos.AA <- Geno.freq[,"AA"] >= 0.8
no.poli.AA <- rownames(Geno.freq [no.polimorficos.AA,])
length(no.poli.AA)
names.marker.AA <- c (rownames(Geno.freq [no.poli.AA,]))
5405

##### No polimorficos para BB
no.polimorficos.BB <- Geno.freq[,"BB"] >= 0.8
no.poli.BB <- rownames(Geno.freq [no.polimorficos.BB,])
names.marker.BB <- c (rownames(Geno.freq [no.poli.BB,]))

length(no.poli.BB)


Ceibo.Carumbe.G2a <- drop.markers (Ceibo.Carumbe.G2, names.marker.AA)
totmar(Ceibo.Carumbe.G2a)

Ceibo.Carumbe.G2b <- drop.markers (Ceibo.Carumbe.G2a, names.marker.BB)
totmar (Ceibo.Carumbe.G2b)

summary (Ceibo.Carumbe.G2b)



####################################################
### code chunk number 21: genofreqbyind (eval = FALSE)
###################################################
g <- pull.geno(Ceibo.Carumbe.G2b )
 gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:2)))
 gfreq <- t(t(gfreq) / colSums(gfreq))
 par(mfrow=c(1,2), las=1)
for(i in 1:2)
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "BB")[i],
       ylim=c(0,1))

 ###################################################
 ### code chunk number 23: triangleplot
 ###################################################
 source("holmans_triangle.R")
 par(mar=rep(0.1,4), pty="s")
 triplot(labels=c("AA","AB","BB"))
 tripoints(gfreq, cex=0.8)
 tripoints(c(0.25, 0.5, 0.25), col="red", lwd=2, cex=1, pch=4)
 
 
 ###################################################
 ### code chunk number 24: pairwiselinkage
 ###################################################
 Ceibo.Carumbe.G2b <- est.rf(Ceibo.Carumbe.G2b)
 
 
 ###################################################
 ### code chunk number 25: checkAlleles
 ###################################################
 checkAlleles ( Ceibo.Carumbe.G2b, threshold=5)
 
 
 ###################################################
 ### code chunk number 26: lodvrf (eval = FALSE)
 ###################################################
rf <- pull.rf (Ceibo.Carumbe.G2b)
lod <- pull.rf (Ceibo.Carumbe.G2b, what="lod")
plot (as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
 
 
 ###################################################
 ### code chunk number 27: lodvrfplot
 ###################################################
 rf <- pull.rf(Ceibo.Carumbe.G2b)
 lod <- pull.rf(Ceibo.Carumbe.G2b, what="lod")
 par(mar=c(4.1,4.1,0.6,0.6), las=1, cex=0.8)
 plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")
 
 
 
 
 ###################################################
 ### code chunk number 30: plotrf (eval = FALSE)
 ###################################################
  plotRF (Ceibo.Carumbe.G2b, alternate.chrid=TRUE)
 
 
 ###################################################
 ### code chunk number 28: forminitialgroups
 ###################################################
 lg <- formLinkageGroups(Ceibo.Carumbe.G2b, max.rf=0.35, min.lod=6)
 table(lg[,2])
 
 
 ###################################################
 ### code chunk number 29: reorganizemarkers
 ###################################################
 Ceibo.Carumbe.G2b <- formLinkageGroups(Ceibo.Carumbe.G2b, max.rf=0.35, min.lod=6, reorgMarkers=TRUE)
 
 
 ###################################################
 ### code chunk number 30: plotrf (eval = FALSE)
 ###################################################
 plotRF (Ceibo.Carumbe.G2b, alternate.chrid=TRUE)
 
 
 #############3 AsMAP
 Ceibo.Carumbe.G2b.asmap.1 <- mstmap ( Ceibo.Carumbe.G2b, chr=1, id = "id")
 pull.map(Ceibo.Carumbe.G2b.asmap.1)
  

write.cross(Ceibo.Carumbe.G2b, format="csv",
            filestem="./genetic_map/Data/procdata/Ceibo.Carumbe.G2b")


######################3
library(onemap)
map_data <- read_mapmaker (dir="./genetic_map/Data/rawdata",file="datos_prueba_onemap2.txt")

Ceibo.Carumbe.G3 <- read_mapmaker (dir="./genetic_map/Data/procdata", 
                                      file = "Ceibo.Carumbe.G3a.txt")

onemap_example_f2 <- read_onemap (inputfile= system.file("extdata/onemap_example_f2.raw", 
                                                        package = "onemap"))

onemap_Ceibo.Carumbe.G3 <- read_onemap(dir="./genetic_map/Data/procdata", 
                                 inputfile =  "Ceibo.Carumbe.G3b.txt")

summary (datos_prueba_onemap)

plot (Ceibo.Carumbe.G3)
plot (onemap_Ceibo.Carumbe.G3)

Ceibo.Carumbe.G3_test <- test_segregation (Ceibo.Carumbe.G3)
#DM.SO_f2_test.0.8.1 <- test_segregation (onemap_DM.SO_f2.0.8.1)

print (Ceibo.Carumbe.G3_test)


# la función Bonferroni_alpha muestra el valor alfa que se debe considerar 
# para este número de loci si se aplica la 
# corrección de Bonferroni con alfa global de 0,05

Bonferroni_alpha(Ceibo.Carumbe.G3_test)


# ver qué marcadores están distorsionados bajo el criterio de Bonferroni
class(DM.SO.G3_test)
plot (Ceibo.Carumbe.G3_test)


# El gráfico se explica por sí mismo: 
# los valores de p se transformaron utilizando -log10 (valores p)
# para una mejor visualización. 
# Una línea vertical muestra el umbral para las pruebas si se aplica
# la corrección de Bonferroni. 
# Se identifican pruebas significativas y no significativas.

# Por favor, recuerde que la corrección de Bonferroni 
# es conservadora, y también que el descarte de datos de marcador 
# puede no ser un buen enfoque para su análisis. 
# Este gráfico es solo para sugerir un criterio, 
# así que úsalo con precaución


# list of markers with non-distorted segregation
lnds <- select_segreg (Ceibo.Carumbe.G3_test, distorted = FALSE, numbers = FALSE) 
length(lnds)

# list of markers with distorted segregation
lds <- select_segreg (Ceibo.Carumbe.G3_test , distorted = TRUE, numbers = TRUE) 
lds # estos me los tengo que sacar


write.table(lds, file = "./genetic_map/Data/procdata/marcadores_distorsionados.txt", 
            append = FALSE, quote = TRUE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE) 

############# 
Ceibo.Carumbe.G4 <- read_mapmaker (dir="./genetic_map/Data/procdata", 
                                   file = "Ceibo.Carumbe.G4.txt")

DM.SO.G4 <- read_onemap (dir="./Data/procdata", 
                         inputfile = "DMSO_G4.txt")

plot(Ceibo.Carumbe.G4)

Ceibo.Carumbe.G4_test <- test_segregation (Ceibo.Carumbe.G4)
#DM.SO_f2_test.0.8.1 <- test_segregation (onemap_DM.SO_f2.0.8.1)

print(Ceibo.Carumbe.G4_test)


# la función Bonferroni_alpha muestra el valor alfa que se debe considerar 
# para este número de loci si se aplica la 
# corrección de Bonferroni con alfa global de 0,05

Bonferroni_alpha(Ceibo.Carumbe.G4_test)


# ver qué marcadores están distorsionados bajo el criterio de Bonferroni
class(DM.SO.G4_test)
plot (Ceibo.Carumbe.G4_test)

## Estimating two-point recombination fractions
###############
(LOD_sug_Ceibo.Carumbe.G4 <- suggest_lod(Ceibo.Carumbe.G4))

twopts_Ceibo.Carumbe.G4 <- rf_2pts (Ceibo.Carumbe.G4, LOD = LOD_sug_Ceibo.Carumbe.G4  )

# Using only recombinations informations

mark_all_Ceibo.Carumbe.G4 <- make_seq (twopts_Ceibo.Carumbe.G4, "all")

#str(mark_all_DM.SO.G4)

#Forming the groups
#You can assign markers to linkage groups using the
# function group

LGs_Ceibo.Carumbe.G4 <- group (mark_all_Ceibo.Carumbe.G4)

LGs_DM.SO.G4


# LOD = 4.5 
#Group 2 : 6 markers
#M67 M68 M69 M70 M71 M72 

#Group 13 : 2 markers
#M1003 M1004 


#Group 18 : 4 markers
#M1329 M1330 M1331 M1332

#LOD = LOD_sug_DM.SO.G4

#Group 2 : 6 markers
#M67 M68 M69 M70 M71 M72 

#Group 16 : 2 markers
#M1003 M1004 

#Group 22 : 4 markers
#M1329 M1330 M1331 M1332 


# Ordering markers within linkage groups

#set_map_fun(type = "haldane")
set_map_fun(type = "kosambi")


#### CHR. 1

LG1          <- make_seq (LGs_DM.SO.G4, 1)
LG1_ord      <- order_seq (LG1)
LG1_frame    <- make_seq(LG1_ord, "force")
LG1_test_map <- map(LG1_frame)



svg (filename="./Figures/Figures_GENO/DM.SO.LG1.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG1_test_map, main = "LG1", inter = FALSE)

dev.off()


#### LG. 2
#LG2          <- make_seq (LGs_DM.SO.G4, 2)
#LG2_ord      <- order_seq (LG2)
#LG2_frame    <- make_seq(LG2_ord, "force")
#LG2_test_map <- map(LG2_frame)

#svg (filename="./Figures/Figures_GENO/DM.SO.LG2.svg", 
    # width=7, 
     #height=5, 
     #pointsize=12)

#rf_graph_table(LG2_test_map, main = "LG2", inter = FALSE)

#dev.off()

#### LG. 3
LG3          <- make_seq (LGs_DM.SO.G4, 3)
LG3_ord      <- order_seq (LG3)
LG3_frame    <- make_seq(LG3_ord, "force")
LG3_test_map <- map(LG3_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG3.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG3_test_map, main = "LG3", inter = FALSE)

dev.off()


#### LG. 4
LG4          <- make_seq (LGs_DM.SO.G4, 4)
LG4_ord      <- order_seq (LG4)
LG4_frame    <- make_seq(LG4_ord, "force")
LG4_test_map <- map(LG4_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG4.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG4_test_map, main = "LG4", inter = FALSE)

dev.off()

#### LG. 5
LG5          <- make_seq (LGs_DM.SO.G4, 5)
LG5_ord      <- order_seq (LG5)
LG5_frame    <- make_seq(LG5_ord, "force")
LG5_test_map <- map(LG5_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG5.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG5_test_map, main = "LG5", inter = FALSE)

dev.off()

#### LG. 6
LG6          <- make_seq (LGs_DM.SO.G4, 6)
LG6_ord      <- order_seq (LG6)
LG6_frame    <- make_seq(LG6_ord, "force")
LG6_test_map <- map(LG6_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG6.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG6_test_map, main = "LG6", inter = FALSE)

dev.off()

#### LG. 7
LG7          <- make_seq (LGs_DM.SO.G4, 7)
LG7_ord      <- order_seq (LG7)
LG7_frame    <- make_seq(LG7_ord, "force")
LG7_test_map <- map(LG7_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG7.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG7_test_map, main = "LG7", inter = FALSE)

dev.off()


#### LG. 8
LG8          <- make_seq (LGs_DM.SO.G4, 8)
LG8_ord      <- order_seq (LG8)
LG8_frame    <- make_seq(LG8_ord, "force")
LG8_test_map <- map(LG8_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG8.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG8_test_map, main = "LG8", inter = FALSE)

dev.off()

#### LG. 9
LG9          <- make_seq (LGs_DM.SO.G4, 9)
LG9_ord      <- order_seq (LG9)
LG9_frame    <- make_seq(LG9_ord, "force")
LG9_test_map <- map(LG9_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG9.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG9_test_map, main = "LG9", inter = FALSE)

dev.off()

#### LG. 10
LG10          <- make_seq (LGs_DM.SO.G4, 10)
LG10_ord      <- order_seq (LG10)
LG10_frame    <- make_seq(LG10_ord, "force")
LG10_test_map <- map(LG10_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG10.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG10_test_map, main = "LG10", inter = FALSE)

dev.off()

#### LG. 11
LG11          <- make_seq (LGs_DM.SO.G4, 11)
LG11_ord      <- order_seq (LG11)
LG11_frame    <- make_seq(LG11_ord, "force")
LG11_test_map <- map(LG11_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG11.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG11_test_map, main = "LG11", inter = FALSE)

dev.off()

#### LG. 12
LG12          <- make_seq (LGs_DM.SO.G4, 12)
LG12_ord      <- order_seq (LG12)
LG12_frame    <- make_seq(LG12_ord, "force")
LG12_test_map <- map(LG12_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG12.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG12_test_map, main = "LG12", inter = FALSE)

dev.off()


#### LG. 13
#LG13          <- make_seq (LGs_DM.SO.G4, 13)
#LG13_ord      <- order_seq (LG13)
#LG13_frame    <- make_seq(LG13_ord, "force")
#LG13_test_map <- map(LG13_frame)

#svg (filename="./Figures/Figures_GENO/DM.SO.LG13.svg", 
 #    width=7, 
  #   height=5, 
   #  pointsize=12)

#rf_graph_table(LG13_test_map, main = "LG13", inter = FALSE)

#dev.off()

#### LG. 14
LG14          <- make_seq (LGs_DM.SO.G4, 14)
LG14_ord      <- order_seq (LG14)
LG14_frame    <- make_seq(LG14_ord, "force")
LG14_test_map <- map(LG14_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG14.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG14_test_map, main = "LG14", inter = FALSE)

dev.off()

#### LG. 15
LG15          <- make_seq (LGs_DM.SO.G4, 15)
LG15_ord      <- order_seq (LG15)
LG15_frame    <- make_seq(LG15_ord, "force")
LG15_test_map <- map(LG15_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG15.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG15_test_map, main = "LG15", inter = FALSE)

dev.off()

#### LG. 16
LG16          <- make_seq (LGs_DM.SO.G4, 16)
LG16_ord      <- order_seq (LG16)
LG16_frame    <- make_seq(LG16_ord, "force")
LG16_test_map <- map(LG16_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG16.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG16_test_map, main = "LG16", inter = FALSE)

dev.off()

#### LG. 17
LG17          <- make_seq (LGs_DM.SO.G4, 17)
LG17_ord      <- order_seq (LG17)
LG17_frame    <- make_seq(LG17_ord, "force")
LG17_test_map <- map(LG17_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG17.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG17_test_map, main = "LG17", inter = FALSE)

dev.off()


#### LG. 18
#LG18          <- make_seq (LGs_DM.SO.G4, 18)
#LG18_ord      <- order_seq (LG18)
#LG18_frame    <- make_seq(LG18_ord, "force")
#LG18_test_map <- map(LG18_frame)

#svg (filename="./Figures/Figures_GENO/DM.SO.LG18.svg", 
 #    width=7, 
  #   height=5, 
   #  pointsize=12)

#rf_graph_table(LG18_test_map, main = "LG18", inter = FALSE)

#dev.off()

#### LG. 19
LG19          <- make_seq (LGs_DM.SO.G4, 19)
LG19_ord      <- order_seq (LG19)
LG19_frame    <- make_seq(LG19_ord, "force")
LG19_test_map <- map(LG19_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG19.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG19_test_map, main = "LG19", inter = FALSE)

dev.off()

#### LG. 20
LG20          <- make_seq (LGs_DM.SO.G4, 20)
LG20_ord      <- order_seq (LG20)
LG20_frame    <- make_seq(LG20_ord, "force")
LG20_test_map <- map(LG20_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG20.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG20_test_map, main = "LG20", inter = FALSE)

dev.off()

str(LGs_DM.SO.G4)
#### LG. 21
LG21          <- make_seq (LGs_DM.SO.G4, 21)
LG21_ord      <- order_seq (LG21)
LG21_frame    <- make_seq(LG21_ord, "force")
LG21_test_map <- map(LG21_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG21.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG21_test_map, main = "LG21", inter = FALSE)

dev.off()

#### LG. 22
LG22          <- make_seq (LGs_DM.SO.G4, 22)
LG22_ord      <- order_seq (LG22)
LG22_frame    <- make_seq(LG22_ord, "force")
LG22_test_map <- map(LG22_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG22.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG22_test_map, main = "LG22", inter = FALSE)

dev.off()

#### LG. 23
LG23          <- make_seq (LGs_DM.SO.G4, 23)
LG23_ord      <- order_seq (LG23)
LG23_frame    <- make_seq(LG23_ord, "force")
LG23_test_map <- map(LG23_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.LG23.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(LG23_test_map, main = "LG23", inter = FALSE)

dev.off()

#### LG. 24
#LG24          <- make_seq (LGs_DM.SO.G4, 24)
#LG24_ord      <- order_seq (LG24)
#LG24_frame    <- make_seq(LG24_ord, "force")
#LG24_test_map <- map(LG24_frame)

#svg (filename="./Figures/Figures_GENO/DM.SO.LG24.svg", 
     #width=7, 
     #height=5, 
     #pointsize=12)

#rf_graph_table(LG24_test_map, main = "LG24", inter = FALSE)

#dev.off()


#### LG. 25
#LG25          <- make_seq (LGs_DM.SO.G4, 25)
#LG25_ord      <- order_seq (LG25)
#LG25_frame    <- make_seq(LG25_ord, "force")
#LG25_test_map <- map(LG25_frame)

#svg (filename="./Figures/Figures_GENO/DM.SO.LG25.svg", 
     #width=7, 
     #height=5, 
     #pointsize=12)

#rf_graph_table(LG25_test_map, main = "LG25", inter = FALSE)

#dev.off()

#### LG. 26
#LG26          <- make_seq (LGs_DM.SO.G4, 26)
#LG26_ord      <- order_seq (LG26)
#LG26_frame    <- make_seq(LG26_ord, "force")
#LG26_test_map <- map(LG26_frame)

#svg (filename="./Figures/Figures_GENO/DM.SO.LG26.svg", 
     #width=7, 
     #height=5, 
     #pointsize=12)

#rf_graph_table(LG26_test_map, main = "LG26", inter = FALSE)

#dev.off()

#### LG. 27
#LG27          <- make_seq (LGs_DM.SO.G4, 27)
#LG27_ord      <- order_seq (LG27)
#LG27_frame    <- make_seq(LG27_ord, "force")
#LG27_test_map <- map(LG27_frame)

#svg (filename="./Figures/Figures_GENO/DM.SO.LG27.svg", 
     #width=7, 
     #height=5, 
     #pointsize=12)

#rf_graph_table(LG27_test_map, main = "LG27", inter = FALSE)

#dev.off()



################


########### riple
LG1_final <- LG1_test_map
#LG2_final <- LG2_test_map
LG3_final <- LG3_test_map
LG4_final <- LG4_test_map
LG5_final <- LG5_test_map
LG6_final <- LG6_test_map
LG7_final <- LG7_test_map
LG8_final <- LG8_test_map
LG9_final <- LG9_test_map
LG10_final <- LG10_test_map
LG11_final <- LG11_test_map
LG12_final <- LG12_test_map
#LG13_final <- LG13_test_map
LG14_final <- LG14_test_map
LG15_final <- LG15_test_map
LG16_final <- LG16_test_map
LG17_final <- LG17_test_map
#LG18_final <- LG18_test_map
LG19_final <- LG19_test_map
LG20_final <- LG20_test_map
LG21_final <- LG21_test_map
LG22_final <- LG22_test_map
LG23_final <- LG23_test_map
#LG24_final <- LG24_test_map
#LG25_final <- LG25_test_map
#LG26_final <- LG26_test_map
#LG27_final <- LG27_test_map


ripple_seq(LG1_final)
#ripple_seq(LG2_final)
ripple_seq(LG3_final)
ripple_seq(LG4_final)
ripple_seq(LG5_final)
ripple_seq(LG6_final)
ripple_seq(LG7_final)
ripple_seq(LG8_final)
ripple_seq(LG9_final)
ripple_seq(LG10_final)
ripple_seq(LG11_final)
ripple_seq(LG12_final)
#ripple_seq(LG13_final)
ripple_seq(LG14_final)
ripple_seq(LG15_final)
ripple_seq(LG16_final)
ripple_seq(LG17_final)
#ripple_seq(LG18_final)
ripple_seq(LG19_final)
ripple_seq(LG20_final)
ripple_seq(LG21_final)
ripple_seq(LG22_final)
ripple_seq(LG23_final)
#ripple_seq(LG24_final)
#ripple_seq(LG25_final)
#ripple_seq(LG26_final)
#ripple_seq(LG27_final)


DMSO.4_map <- list(LG1_final
                  #,LG2_final
                  ,LG3_final
                  ,LG4_final
                  ,LG5_final
                  ,LG6_final
                  ,LG7_final
                  ,LG8_final
                  ,LG9_final
                  ,LG10_final
                  ,LG11_final
                  ,LG12_final
                  #,LG13_final
                  ,LG14_final
                  ,LG15_final
                  ,LG16_final
                  ,LG17_final
                  #,LG18_final
                  ,LG19_final
                  ,LG20_final
                  ,LG21_final
                  ,LG22_final
                  ,LG23_final
                  #,LG24_final
                  #,LG25_final
                  #,LG26_final
                  #,LG27_final
                  )

draw_map (DMSO.4_map, names = FALSE, grid = FALSE, cex.mrk = 0.7) 

write_map (DMSO.4_map, "./Data/procdata/DMSO.5.map.txt")

###############################################################3
# LG por cromosoma asignado 
#set_map_fun(type = "haldane")
set_map_fun(type = "kosambi")


#### CHR. 1
twopts_DM.SO.G4
CHR1          <- make_seq (twopts_DM.SO.G4, "1")
CHR1_ord      <- order_seq (CHR1)
CHR1_frame    <- make_seq(CHR1_ord, "force")
CHR1_test_map <- map(CHR1_frame)



svg (filename="./Figures/Figures_GENO/DM.SO.CHR1.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR1_test_map, main = "CHR1", inter = FALSE)

dev.off()


#### CHR. 2
CHR2          <- make_seq (twopts_DM.SO.G4, "2")
CHR2_ord      <- order_seq (CHR2)
CHR2_frame    <- make_seq(CHR2_ord, "force")
CHR2_test_map <- map(CHR2_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR2.svg", 
 width=7, 
height=5, 
pointsize=12)

rf_graph_table(CHR2_test_map, main = "CHR2", inter = FALSE)

dev.off()

#### CHR. 3
CHR3          <- make_seq (twopts_DM.SO.G4, "3")
CHR3_ord      <- order_seq (CHR3)
CHR3_frame    <- make_seq(CHR3_ord, "force")
CHR3_test_map <- map(CHR3_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR3.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR3_test_map, main = "CHR3", inter = FALSE)

dev.off()


#### CHR. 4
CHR4          <- make_seq (twopts_DM.SO.G4, "4")
CHR4_ord      <- order_seq (CHR4)
CHR4_frame    <- make_seq(CHR4_ord, "force")
CHR4_test_map <- map(CHR4_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR4.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR4_test_map, main = "CHR4", inter = FALSE)

dev.off()

#### CHR. 5
CHR5          <- make_seq (twopts_DM.SO.G4, "5")
CHR5_ord      <- order_seq (CHR5)
CHR5_frame    <- make_seq(CHR5_ord, "force")
CHR5_test_map <- map(CHR5_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR5.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR5_test_map, main = "CHR5", inter = FALSE)

dev.off()

#### CHR. 6
CHR6          <- make_seq (twopts_DM.SO.G4, "6")
CHR6_ord      <- order_seq (CHR6)
CHR6_frame    <- make_seq(CHR6_ord, "force")
CHR6_test_map <- map(CHR6_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR6.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR6_test_map, main = "CHR6", inter = FALSE)

dev.off()

#### CHR. 7
CHR7          <- make_seq (twopts_DM.SO.G4, "7")
CHR7_ord      <- order_seq (CHR7)
CHR7_frame    <- make_seq(CHR7_ord, "force")
CHR7_test_map <- map(CHR7_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR7.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR7_test_map, main = "CHR7", inter = FALSE)

dev.off()


#### CHR. 8
CHR8          <- make_seq (twopts_DM.SO.G4, "8")
CHR8_ord      <- order_seq (CHR8)
CHR8_frame    <- make_seq(CHR8_ord, "force")
CHR8_test_map <- map(CHR8_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR8.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR8_test_map, main = "CHR8", inter = FALSE)

dev.off()

#### CHR. 9
CHR9          <- make_seq (twopts_DM.SO.G4, "9")
CHR9_ord      <- order_seq (CHR9)
CHR9_frame    <- make_seq(CHR9_ord, "force")
CHR9_test_map <- map(CHR9_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR9.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR9_test_map, main = "CHR9", inter = FALSE)

dev.off()

#### CHR. 10
CHR10          <- make_seq (twopts_DM.SO.G4, "10")
CHR10_ord      <- order_seq (CHR10)
CHR10_frame    <- make_seq(CHR10_ord, "force")
CHR10_test_map <- map(CHR10_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR10.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR10_test_map, main = "CHR10", inter = FALSE)

dev.off()

#### CHR. 11
CHR11          <- make_seq (twopts_DM.SO.G4, "11")
CHR11_ord      <- order_seq (CHR11)
CHR11_frame    <- make_seq(CHR11_ord, "force")
CHR11_test_map <- map(CHR11_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR11.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR11_test_map, main = "CHR11", inter = FALSE)

dev.off()

#### CHR. 12
CHR12          <- make_seq (twopts_DM.SO.G4, "12")
CHR12_ord      <- order_seq (CHR12)
CHR12_frame    <- make_seq(CHR12_ord, "force")
CHR12_test_map <- map(CHR12_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR12.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR12_test_map, main = "CHR12", inter = FALSE)

dev.off()


#### CHR. 13
CHR13          <- make_seq (twopts_DM.SO.G4, "13")
CHR13_ord      <- order_seq (CHR13)
CHR13_frame    <- make_seq(CHR13_ord, "force")
CHR13_test_map <- map(CHR13_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR13.svg", 
    width=7, 
   height=5, 
  pointsize=12)

rf_graph_table(CHR13_test_map, main = "CHR13", inter = FALSE)

dev.off()

#### CHR. 14
CHR14          <- make_seq (twopts_DM.SO.G4, "14")
CHR14_ord      <- order_seq (CHR14)
CHR14_frame    <- make_seq(CHR14_ord, "force")
CHR14_test_map <- map(CHR14_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR14.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR14_test_map, main = "CHR14", inter = FALSE)

dev.off()

#### CHR. 15
CHR15          <- make_seq (twopts_DM.SO.G4, "15")
CHR15_ord      <- order_seq (CHR15)
CHR15_frame    <- make_seq(CHR15_ord, "force")
CHR15_test_map <- map(CHR15_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR15.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR15_test_map, main = "CHR15", inter = FALSE)

dev.off()

#### CHR. 16
CHR16          <- make_seq (twopts_DM.SO.G4, "16")
CHR16_ord      <- order_seq (CHR16)
CHR16_frame    <- make_seq(CHR16_ord, "force")
CHR16_test_map <- map(CHR16_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR16.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR16_test_map, main = "CHR16", inter = FALSE)

dev.off()

#### CHR. 17
CHR17          <- make_seq (twopts_DM.SO.G4, "17")
CHR17_ord      <- order_seq (CHR17)
CHR17_frame    <- make_seq(CHR17_ord, "force")
CHR17_test_map <- map(CHR17_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR17.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR17_test_map, main = "CHR17", inter = FALSE)

dev.off()


#### CHR. 18
CHR18          <- make_seq (twopts_DM.SO.G4, "18")
CHR18_ord      <- order_seq (CHR18)
CHR18_frame    <- make_seq(CHR18_ord, "force")
CHR18_test_map <- map(CHR18_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR18.svg", 
    width=7, 
   height=5, 
  pointsize=12)

rf_graph_table(CHR18_test_map, main = "CHR18", inter = FALSE)

dev.off()

#### CHR. 19
CHR19          <- make_seq (twopts_DM.SO.G4, "19")
CHR19_ord      <- order_seq (CHR19)
CHR19_frame    <- make_seq(CHR19_ord, "force")
CHR19_test_map <- map(CHR19_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR19.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR19_test_map, main = "CHR19", inter = FALSE)

dev.off()

#### CHR. 20
CHR20          <- make_seq (twopts_DM.SO.G4, "20")
CHR20_ord      <- order_seq (CHR20)
CHR20_frame    <- make_seq(CHR20_ord, "force")
CHR20_test_map <- map(CHR20_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR20.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR20_test_map, main = "CHR20", inter = FALSE)

dev.off()

str(twopts_DM.SO.G4)
#### CHR. 21
CHR21          <- make_seq (twopts_DM.SO.G4, "21")
CHR21_ord      <- order_seq (CHR21)
CHR21_frame    <- make_seq(CHR21_ord, "force")
CHR21_test_map <- map(CHR21_frame)

svg (filename="./Figures/Figures_GENO/DM.SO.CHR21.svg", 
     width=7, 
     height=5, 
     pointsize=12)

rf_graph_table(CHR21_test_map, main = "CHR21", inter = FALSE)

dev.off()

########### riple
CHR1_final <- CHR1_test_map
CHR2_final <- CHR2_test_map
CHR3_final <- CHR3_test_map
CHR4_final <- CHR4_test_map
CHR5_final <- CHR5_test_map
CHR6_final <- CHR6_test_map
CHR7_final <- CHR7_test_map
CHR8_final <- CHR8_test_map
CHR9_final <- CHR9_test_map
CHR10_final <- CHR10_test_map
CHR11_final <- CHR11_test_map
CHR12_final <- CHR12_test_map
CHR13_final <- CHR13_test_map
CHR14_final <- CHR14_test_map
CHR15_final <- CHR15_test_map
CHR16_final <- CHR16_test_map
CHR17_final <- CHR17_test_map
CHR18_final <- CHR18_test_map
CHR19_final <- CHR19_test_map
CHR20_final <- CHR20_test_map



ripple_seq(CHR1_final)
ripple_seq(CHR2_final)
ripple_seq(CHR3_final)
ripple_seq(CHR4_final)
ripple_seq(CHR5_final)
ripple_seq(CHR6_final)
ripple_seq(CHR7_final)
ripple_seq(CHR8_final)
ripple_seq(CHR9_final)
ripple_seq(CHR10_final)
ripple_seq(CHR11_final)
ripple_seq(CHR12_final)
ripple_seq(CHR13_final)
ripple_seq(CHR14_final)
ripple_seq(CHR15_final)
ripple_seq(CHR16_final)
ripple_seq(CHR17_final)
ripple_seq(CHR18_final)
ripple_seq(CHR19_final)
ripple_seq(CHR20_final)



DMSO.4.1_map <- list(CHR1_final
                   ,CHR2_final
                   ,CHR3_final
                   ,CHR4_final
                   ,CHR5_final
                   ,CHR6_final
                   ,CHR7_final
                   ,CHR8_final
                   ,CHR9_final
                   ,CHR10_final
                   ,CHR11_final
                   ,CHR12_final
                   ,CHR13_final
                   ,CHR14_final
                   ,CHR15_final
                   ,CHR16_final
                   ,CHR17_final
                   ,CHR18_final
                   ,CHR19_final
                   ,CHR20_final
                   )

draw_map (DMSO.4.1_map, names = FALSE, grid = FALSE, cex.mrk = 0.7) 

write_map (DMSO.4.1_map, "./Data/procdata/DMSO.5.1.map.txt")













#########################################
getwd()

DM.SO.G5 <- read.cross (format="csv",
                        dir="./Data/procdata",
                        file="DM.SO_G5.csv", 
                        na.strings= "NA",sep=";",dec=".",
                        genotypes=c("AA","AB","BB"),
                        #alleles=c("A","B"),
                        estimate.map=FALSE, convertXdata=TRUE, error.prob=0.0001)
str(DM.SO.G5)
summary(DM.SO.G5)
plot(DM.SO.G5)
str(DM.SO.G5)
nind(DM.SO.G5)
nchr(DM.SO.G5)
totmar(DM.SO.G5)
nmar(DM.SO.G5)
nphe(DM.SO.G5)

geno.image(DM.SO.G5)
plotMissing(DM.SO.G5)


#################
DM.SO.G5 <- est.rf(DM.SO.G5)

plotRF (DM.SO.G5, alternate.chrid=TRUE)

pq.diagnostics (crossobj=DM.SO.G5)
mq.diagnostics (crossobj=DM.SO.G5)

# QTL_SMA
QTL.result <- qtl.analysis (crossobj=DM.SO.G5,step=0,
                            method='SIM', trait="k",
                            threshold="Li&Ji",distance=30,
                            cofactors=NULL,window.size=30)
# QTL_SIM
QTL.result <- qtl.analysis ( crossobj=DM.SO.G5, step=5,
                             method='SIM',trait="k", threshold="Li&Ji",
                             distance=30,cofactors=NULL,window.size=30)

# QTL CIM
cofactors <- as.vector (QTL.result$selected$marker)

QTL.result <- qtl.analysis ( crossobj=DM.SO.G5, step=5,
                             method='CIM', trait="k", threshold="Li&Ji", distance=30,
                             cofactors=cofactors, window.size=30)

#############3 B 

# QTL_SMA
QTL.result <- qtl.analysis (crossobj=DM.SO.G5,step=0,
                            method='SIM', trait="beta",
                            threshold="Li&Ji",distance=30,
                            cofactors=NULL,window.size=30)
# QTL_SIM
QTL.result <- qtl.analysis ( crossobj=DM.SO.G5, step=5,
                             method='SIM',trait="beta", threshold="Li&Ji",
                             distance=30,cofactors=NULL,window.size=30)

# QTL CIM
cofactors <- as.vector (QTL.result$selected$marker)

QTL.result <- qtl.analysis ( crossobj=DM.SO.G5, step=5,
                             method='CIM', trait="beta", threshold="Li&Ji", distance=30,
                             cofactors=cofactors, window.size=30)


#############3 tm

# QTL_SMA
QTL.result <- qtl.analysis (crossobj=DM.SO.G5,step=0,
                            method='SIM', trait="tm",
                            threshold="Li&Ji",distance=30,
                            cofactors=NULL,window.size=30)
# QTL_SIM
QTL.result <- qtl.analysis ( crossobj=DM.SO.G5, step=5,
                             method='SIM',trait="tm", threshold="Li&Ji",
                             distance=30,cofactors=NULL,window.size=30)

# QTL CIM
cofactors <- as.vector (QTL.result$selected$marker)

QTL.result <- qtl.analysis ( crossobj=DM.SO.G5, step=5,
                             method='CIM', trait="tm", threshold="Li&Ji", distance=30,
                             cofactors=cofactors, window.size=30)


## Not run: 
data (SxM_geno)
data (SxM_map)
data (SxMxE_pheno)

P.data <- SxMxE_pheno
G.data <- SxM_geno
map.data <- SxM_map

cross.data <- qtl.cross (P.data, G.data, map.data, cross='dh',
                         heterozygotes=FALSE)

summary (cross.data)

## Pheno Quality
pq.diagnostics (crossobj=cross.data, boxplot =FALSE)

## Marker Quality
mq.diagnostics (crossobj=cross.data,I.threshold=0.1,
                p.val=0.01,na.cutoff=0.1)

# QTL_SIM
QTL.result <- qtl.memq (crossobj = cross.data, P.data = P.data,
                        env.label = c('ID91','ID92','MAN92','MTd91',
                                      'MTd92','MTi91','MTi92','SKs92','WA91','WA92'),
                        trait = 'yield', step = 10, method = 'SIM',
                        threshold = 'Li&Ji', distance = 50, cofactors = NULL,
                        window.size = 50)

## QTL_CIM
QTL.result <- qtl.memq (crossobj = cross.data, P.data = P.data,
                        env.label = c('ID91','ID92','MAN92','MTd91','MTd92',
                                      'MTi91','MTi92','SKs92','WA91','WA92'),
                        trait = 'yield', step = 10, method = 'CIM',
                        threshold = 'Li&Ji', distance = 50,
                        cofactors = QTL.result$selected$marker, window.size = 50)

## End(Not run)










QTL.result <- qtl.analysis ( crossobj=DM.SO.cross_onemap_qtl, step=5,
                             method='CIM', trait="em.k", threshold="Li&Ji", distance=30,
                             cofactors=NULL, window.size=30)





all_mark <- make_seq(twopts_DM.SO.G3a,"all")
all_mark1 <- drop_marker(all_mark, lds )


# si max.rf <0.5 podemos afirmar que los marcadores están vinculados.
# La puntuación LOD es la estadística utilizada para evaluar 
# la importancia de la prueba para max.rf = 0,50.

twopts_DM.SO.G3 <- rf_2pts (DM.SO.G3)



print(twopts_DM.SO.G3)


###############
(LOD_sug_DM.SO.G3<- suggest_lod(DM.SO.G3))

(LOD_sug_DM.SO_f2.1 <- suggest_lod(onemap_DM.SO_f2.0.8.1))

twopts_DM.SO.G3a <- rf_2pts (DM.SO.G3, LOD =  6.630826)
mark_all_DM.SO.G3a <- make_seq (twopts_DM.SO.G3a, "all")
class(mark_all_DM.SO.G3a)
LGs_DM.SO.G3a <- group(mark_all_DM.SO.G3a)





LGs_DM.SO.G3a
set_map_fun(type = "haldane")

LGs <- map()
rf_graph_table (mark_all_DM.SO.G3a, inter = FALSE)

LG2_f2_ord <- order_seq(input.seq = LG2_f2, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3)

print(twopts_DM.SO_f2a)

print(twopts_DM.SO_f2.1a)

# Se va a usar dos estrategias 

# A : Usando solo la frecuencia de  recombinacion 
## A.1 Assigning markers to linkage groups


mark_all_DM.SO_f2.1 <- make_seq (twopts_DM.SO_f2.1, "1")
class(mark_all_DM.SO_f2)

mark_all_DM.SO_f2a <- make_seq (twopts_DM.SO_f2a, "all")
class(mark_all_DM.SO_f2a)
## A.2 Forming the groups

LGs_DM.SO_f2 <- group(mark_all_DM.SO_f2)

LGs_DM.SO_f2a <- group(mark_all_DM.SO_f2a)

LGs_DM.SO_f2

LGs_DM.SO_f2a

# B : Usando la frecuencia de  recombinacion y la referencia del genoma
# B.1 
twopts_DM.SO_f2.1 <- rf_2pts (onemap_DM.SO_f2.0.8.1)
twopts_DM.SO_f2.1a <- rf_2pts (onemap_DM.SO_f2.0.8.1, LOD =  6.701726)




##################################3
DM.SO.cross_onemap_qtl <- read.cross (format="csv",
                           dir="./Data/procdata",
                           file="DM.SO.onemap0.8_qtl.csv", 
                           na.strings= "NA",sep=";",dec=".",
                           genotypes=c("AA","BB","AB"), alleles=c("A","B"),
                           estimate.map=FALSE, convertXdata=TRUE, error.prob=0.0001)


#################
DM.SO.cross_onemap_qtl <- est.rf(DM.SO.cross_onemap_qtl)

plotRF (DM.SO.cross_onemap_qtl, alternate.chrid=TRUE)

pq.diagnostics (crossobj=DM.SO.cross_onemap_qtl)
mq.diagnostics (crossobj=DM.SO.cross_onemap_qtl)


#############
phenames (DM.SO.cross_onemap_qtl)



# QTL_SIM




DMSO_sug <- calc.genoprob(DM.SO.cross_onemap_qtl, step=1)

out.em.DMSO_sug <- scanone(DMSO_sug, pheno.col=4)

summary(out.em.DMSO_sug)

summary(out.em, threshold=3)

plot(out.em.DMSO_sug)

### hasta aca todos los mardadores no son monomorficos absolutos
DM.SO.cross.2b <- est.rf(DM.SO.cross.2b)

checkAlleles(DM.SO.cross.2b, threshold=5)
plotRF (DM.SO.cross.2b, alternate.chrid=TRUE)

#############
phenames(DM.SO.cross.2b)
QTL_SMA
QTL.result <- qtl.analysis (crossobj=DM.SO.cross.2b,step=0,
                            method='SIM', trait="random", threshold="Li&Ji",
                            distance=30, cofactors=NULL,
                            window.size=30)
phenames(DM.SO.cross.2b)

write.cross (DM.SO.cross.2b, format="csv",
            filestem="./Data/rawdata/DM.SO.cross.2b")


write.cross()

### voy a usar ASMap ################33







####  Voy a probar si el analisis de QTL corre ##############
summary(DM.SO.cross.2b)


## Not run: 
QTL_SMA
QTL.result <- qtl.analysis (crossobj=DM.SO.cross.2b,step=0,
                            method='SIM', trait="random", threshold="Li&Ji",
                            distance=30, cofactors=NULL,
                            window.size=30)

## End(Not run)

# QTL_SIM
QTL.result <- qtl.analysis ( crossobj=DM.SO.cross.2b, step=5,
                             method='SIM',trait="random", threshold="Li&Ji",
                             distance=30,cofactors=NULL,window.size=30)

# QTL CIM
cofactors <- as.vector (QTL.result$selected$marker)

QTL.result <- qtl.analysis ( crossobj=cross.data, step=5,
                             method='CIM', trait="height", threshold="Li&Ji", distance=30,
                             cofactors=cofactors, window.size=30)






rf <- pull.rf(DM.SO.cross.2b)
lod <- pull.rf(DM.SO.cross.2b, what="lod")
plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score")

lg <- formLinkageGroups(DM.SO.cross.2b, max.rf=0.3, min.lod=9)
table(lg[,2])

DM.SO.cross.2c <- orderMarkers(DM.SO.cross.2b, chr = 1)

pull.map(DM.SO.cross.2b, chr=1)
pull.map(DM.SO.cross.2c, chr=1)

DM.SO.cross.2c  <- formLinkageGroups(DM.SO.cross.2b, max.rf=0.3, min.lod=9, 
                                     reorgMarkers=TRUE)

plotRF (DM.SO.cross.2c, alternate.chrid=TRUE)
plotMap (DM.SO.cross.2b)


mapthis <- orderMarkers(mapthis, chr=5)


###### vuelvo a calcular la frecuencia de segregacion ###########

gt.2 <- geno.table(DM.SO.cross.2b)
gt.2[gt.2$P.value < 0.05/totmar(DM.SO.cross.2b),]

g.2 <- pull.geno(DM.SO.cross.2b)
gfreq.2 <- apply(g.2, 1, function(a) table(factor(a, levels=1:3)))
gfreq.2 <- t(t(gfreq.2) / colSums(gfreq.2))
par(mfrow=c(1,3), las=1)
for(i in 1:3)
  plot(gfreq.2[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i],
       ylim=c(0,1))

todrop.2 <- rownames(gt.2[gt.2$P.value < 1e-10,])
length(todrop.2)



5959 - 1278 - 1187
head(gt)
class(gt)
str(gt)
chr.2 <- gt %>% filter (chr == 2) %>%
          arrange(desc(P.value))
row.names(chr.2)


Gm02_13476678_G_A <-

View(chr.2)
class(gt)


arrange(gt, P.value)


write.table(gt, file = "./Data/rawdata/gt_DM.SO.txt", 
            append = FALSE, quote = TRUE, sep = ",",
            eol = "\n", na = "NA", dec = ".")

g <- pull.geno(DM.SO.cross.2)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 1:3)
   plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i],
          ylim=c(0,1))

############### 

todrop <- rownames(gt[gt$P.value < 1e-10,])
length(todrop)


todropa <- rownames(gt[gt$P.value > 1e+8,])
length(todropa)
5822 + 137

length(todrop)
Gm07_8591371_G_A <- pull.markers(DM.SO.cross.2, "Gm07_8591371_G_A" )
View(todrop)
> mapthis <- drop.markers(mapthis, todrop)






cross.data.3 <- drop.markers (cross.data.2, names.marker)

summary (cross.data.3)

plot.missing (cross.data.3, main="cross.data.3.MAF")

geno.image (cross.data.3 , reorder=FALSE, main="cross.data.3",
            alternate.chrid=FALSE)









DM.SO.cross <- read.cross (format="csv",
                           dir="./Data/rawdata",
                           file="previo_DM.SO_cross.csv", 
                           na.strings= "NA",sep=";",dec=",",
                           genotypes=c("AA","BB","AB"), alleles=c("A","B"),
                           estimate.map=FALSE, convertXdata=TRUE, error.prob=0.0001)



### pvalue
todrop <- rownames(gt[gt$P.value < 1e-30,])
length(todrop)
DM.SO.cross.3 <- drop.markers(DM.SO.cross.2, todrop)

nmar(DM.SO.cross.3)
totmar(DM.SO.cross.3)

g <- pull.geno(DM.SO.cross.3)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,3), las=1)
for(i in 1:3)
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i],
       ylim=c(0,1))





par(mfrow=c(1,3), las=1)
for(i in 1:3)
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i],
       ylim=c(0,1))







total <- nind (DM.SO.cross.5)
gt <- geno.table(DM.SO.cross.5)
missing <- gt[,2]
num.obs <- c(total - missing)
obs.alelo <- gt[, c(3,4)]
Geno.freq <- as.vector(obs.alelo/num.obs)

SNPs.MAF <-  Geno.freq[,2] < 0.05 
rownames(Geno.freq [SNPs.MAF,])
names.marker <- c (rownames(Geno.freq [SNPs.MAF,]))
length(names.marker)

todrop <- rownames(gt[gt$P.value < 1e-10,])
length(todrop)
mapthis <- drop.markers(mapthis, todrop)



g <- pull.geno(DM.SO.cross.2)
table(g[10,], g[101,])

print(dup <- findDupMarkers(DM.SO.cross.2, exact.only=FALSE))

class (dup)
str(dup)


# Estos son con el mapa estimado estimado
DM.SO.cross.3 <- read.cross (format="csv",
                           dir="./Data/procdata",
                           file="DM.SO.cross.3.csv", 
                           na.strings= "-",sep=",",dec=".",
                           genotypes=c("AA","BB","AB"), alleles=c("A","B"),
                           estimate.map=FALSE, convertXdata=TRUE, error.prob=0.0001)

jittermap(DM.SO.cross.3)

summary(DM.SO.cross.3)

nind(DM.SO.cross.3)
nchr(DM.SO.cross.3)
totmar(DM.SO.cross.3)
nmar(DM.SO.cross.3)
nphe(DM.SO.cross.3)

plot(DM.SO.cross.3)

plotMissing(DM.SO.cross.3)
nullmarkers(DM.SO.cross.3)
DM.SO.cross.3 <- drop.nullmarkers(DM.SO.cross.3)

# Eliminar marcadores con mas del 50 % de datos faltantes
n.missing <- nmissing(DM.SO.cross.3, what="mar")[(nmissing(DM.SO.cross.3, what="mar"))/sum(summary(DM.SO.cross.3)$n.ind)> 0.50]  
names.marker <- c (names(n.missing))
length(names.marker)

DM.SO.cross.4 <- drop.markers (DM.SO.cross.3, names.marker)

# Ver la info de la poblacion sin los missing markers
summary (DM.SO.cross.4)
plotMissing(DM.SO.cross.4)

geno.image (DM.SO.cross.4, reorder=T,
            alternate.chrid=FALSE)

############################################################
# Single-QTL analysis
############################################################
DM.SO.cross.4 <- calc.genoprob(DM.SO.cross.4, step=0)

out.em <- scanone(DM.SO.cross.4)
out <- scanone(DM.SO.cross.4)
out.cim <- cim(DM.SO.cross.4)
, n.marcovar=3)
plot(out, out.cim, chr=c(1,4,6,15), col=c("blue", "red"))

add.cim.covar(out.cim, chr=c(1,4,6,15))
summary(out.em)

summary(out.em, threshold=3)

plot(out.em)

out.hk <- scanone(DM.SO.cross.4, method="hk")

plot(out.em, out.hk, col=c("blue", "red"))

plot(out.em, col="blue")
plot(out.hk, col="red", add=TRUE)

plot(out.hk - out.em, ylim=c(-0.3, 0.3), ylab="LOD(HK)-LOD(EM)")

sug <- sim.geno(sug, step=1, n.draws=64)
out.imp <- scanone(sug, method="imp")

plot(out.em, out.hk, out.imp, col=c("blue", "red", "green"))

plot(out.em, out.hk, out.imp, col=c("blue", "red", "green"), chr=c(7,15))

plot(out.imp - out.em, out.hk - out.em, col=c("green", "red"), ylim=c(-1,1))

############################################################
# Permutation tests
############################################################
load(url("http://www.rqtl.org/various.RData"))

operm <- scanone(sug, method="hk", n.perm=1000)

plot(operm)

summary(operm)
summary(operm, alpha=c(0.05, 0.2))

summary(out.hk, perms=operm, alpha=0.2, pvalues=TRUE)

############################################################
# Interval estimates of QTL location
############################################################
lodint(out.hk, chr=7)
bayesint(out.hk, chr=7)

lodint(out.hk, chr=7, expandtomarkers=TRUE)
bayesint(out.hk, chr=7, expandtomarkers=TRUE)

lodint(out.hk, chr=7, drop=2)
bayesint(out.hk, chr=7, prob=0.99)

lodint(out.hk, chr=15)
bayesint(out.hk, chr=15)

############################################################
# QTL effects
############################################################
max(out.hk)
mar <- find.marker(sug, chr=7, pos=47.7)
plotPXG(sug, marker=mar)

effectplot(sug, mname1=mar)

effectplot(sug, mname1="7@47.7")

max(out.hk, chr=15)
mar2 <- find.marker(sug, chr=15, pos=12)
plotPXG(sug, marker=mar2)
effectplot(sug, mname1="15@12")

plotPXG(sug, marker=c(mar, mar2))
plotPXG(sug, marker=c(mar2, mar))

effectplot(sug, mname1="7@47.7", mname2="15@12")
effectplot(sug, mname2="7@47.7", mname1="15@12")

############################################################
# Other phenotypes
############################################################
out.hr <- scanone(sug, pheno.col=2, method="hk")

out.bw <- scanone(sug, pheno.col="bw", method="hk")

out.logbw <- scanone(sug, pheno.col=log(sug$pheno$bw), method="hk")

out.all <- scanone(sug, pheno.col=1:4, method="hk")

summary(out.all, threshold=3)

summary(out.all, threshold=3, lodcolumn=4)

summary(out.all, threshold=3, format="allpeaks")

summary(out.all, threshold=3, format="allpheno")

summary(out.all, threshold=3, format="tabByCol")
summary(out.all, threshold=3, format="tabByChr")

############################################################
# Two-dimensional, two-QTL scans
############################################################
sug <- calc.genoprob(sug, step=2)

out2 <- scantwo(sug, method="hk")

plot(out2)

plot(out2, lower="fv1")

plot(out2, lower="fv1", upper="av1")

operm2 <- scantwo(sug, method="hk", n.perm=5)

summary(out2, perms=operm2, alpha=0.2, pvalues=TRUE)

############################################################
# Multiple-QTL analyses
############################################################
sug <- calc.genoprob(sug, step=1)

qtl <- makeqtl(sug, chr=c(7,15), pos=c(47.7, 12), what="prob")

out.fq <- fitqtl(sug, qtl=qtl, method="hk")
summary(out.fq)

summary(fitqtl(sug, qtl=qtl, method="hk", get.ests=TRUE, dropone=FALSE))

out.fqi <- fitqtl(sug, qtl=qtl, method="hk", formula=y~Q1*Q2)
out.fqi <- fitqtl(sug, qtl=qtl, method="hk", formula=y~Q1+Q2+Q1:Q2)
summary(out.fqi)

addint(sug, qtl=qtl, method="hk")

rqtl <- refineqtl(sug, qtl=qtl, method="hk")
rqtl

summary(out.fqr <- fitqtl(sug, qtl=rqtl, method="hk"))

plotLodProfile(rqtl)

plot(out.hk, chr=c(7,15), col="red", add=TRUE)

out.aq <- addqtl(sug, qtl=rqtl, method="hk")

plot(out.aq)

print(pen <- calc.penalties(operm2))

out.sq <- stepwiseqtl(sug, max.qtl=5, penalties=pen, method="hk", verbose=2)
out.sq

# end of rqtltour2.R





dupmar <- findDupMarkers(DM.SO.cross.4)
dupmar.adjonly <- findDupMarkers(DM.SO.cross.4, adjacent.only=TRUE) # finds 4 pairs
# one might consider dropping the extra markers
totmar(DM.SO.cross.4) # 173 markers
DM.SO.cross.5 <- drop.markers(DM.SO.cross.4, unlist(dupmar.adjonly))
totmar(DM.SO.cross.5) # 169 markers

#  MAF se sacan los marcadores con un frecuencia menor al 5% y mayor al 95

total <- nind (DM.SO.cross.5)
gt <- geno.table(DM.SO.cross.5)
missing <- gt[,2]
num.obs <- c(total - missing)
obs.alelo <- gt[, c(3,4)]
Geno.freq <- as.vector(obs.alelo/num.obs)

SNPs.MAF <-  Geno.freq[,2] < 0.05 
rownames(Geno.freq [SNPs.MAF,])
names.marker <- c (rownames(Geno.freq [SNPs.MAF,]))
length(names.marker)

DM.SO.cross.6 <- drop.markers (DM.SO.cross.5, names.marker)
summary (DM.SO.cross.6)
plotMissing(DM.SO.cross.6)

geno.image (DM.SO.cross.6, reorder=T,
            alternate.chrid=FALSE)




QTL.result <- qtl.analysis (
  crossobj=DM.SO.cross.6
  step=0
  method='SIM'
  trait="em.k" ,
  threshold="Li&Ji" 
  distance=30 
  cofactors=NULL
  window.size=30
  

plotMap(sug)
plotPheno(sug, pheno.col=1)
plotPheno(sug, pheno.col=2)
plotPheno(sug, pheno.col=3)
plotPheno(sug, pheno.col=4)
plotPheno(sug, pheno.col=5)
plotPheno(sug, pheno.col=6)




summary(DM.SO.cross)

(G.soja.ET.data.2$Genotype %in% P.soja.ET.data.1$Genotype)
setdiff(G.soja.ET.data.2$Genotype, P.soja.ET.data.1$Genotype)


map.soja.ET.data <- read.table ("./Data/procdata/DM_SO.map.txt",
                           header = F, sep = "\t",
                           dec = ".", na.strings = "-")

dim(map.soja.ET.data)
head(map.soja.ET.data)

P.data    <- P.soja.ET.data.1
#P.data.1 <-na.omit(P.soja.ET.data.1)

#P.E1.data    <- P.E1
#P.E2.data    <- P.E2
#P.E3.data    <- P.E3

G.data    <- G.soja.ET.data.2
str(G.data)

map.data  <- map.soja.ET.data


cross.data <- qtl.cross (P.data, G.data, map.data,
                         cross='f2', heterozygotes = TRUE)
summary (cross.data)

pq.diagnostics (DM.SO.cross, boxplot = TRUE, qqplot = FALSE,
                scatterplot = TRUE,heatplot = TRUE)

# data(badorder)
DM.SO.cross <- est.rf(DM.SO.cross)
plotRF(DM.SO.cross)

genetic.map <- est.map(DM.SO.cross)

logliks <- sapply(genetic.map, attr, "loglik")
plotMap(genetic.map)
#plotMap(DM.SO.cross)

DM.SO.cross.2 <- replace.map(DM.SO.cross, genetic.map)

write.cross(DM.SO.cross.2, format="csv",
            filestem="DM.SO.cross.2")

pq.diagnostics (DM.SO.cross.2, boxplot = TRUE, qqplot = FALSE,
                scatterplot = TRUE,heatplot = TRUE)

DM.SO.cross.6 <- calc.genoprob(DM.SO.cross.6, step=0)

sug <- calc.genoprob(sug, step=1)

out.em <- scanone(DM.SO.cross.6)



summary(out.em)

summary(out.em, threshold=3)

plot(out.em)

out.hk <- scanone(sug, method="hk")

plot(out.em, out.hk, col=c("blue", "red"))

plot(out.em, col="blue")
plot(out.hk, col="red", add=TRUE)

plot(out.hk - out.em, ylim=c(-0.3, 0.3), ylab="LOD(HK)-LOD(EM)")

sug <- sim.geno(sug, step=1, n.draws=64)
out.imp <- scanone(sug, method="imp")

plot(out.em, out.hk, out.imp, col=c("blue", "red", "green"))

plot(out.em, out.hk, out.imp, col=c("blue", "red", "green"), chr=c(7,15))

plot(out.imp - out.em, out.hk - out.em, col=c("green", "red"), ylim=c(-1,1))



out.em <- scanone(DM.SO.cross.4, method="em")


out.hk <- scanone(hyper, method="hk")
summary(out.em)

plot(qtl


## Not run: 
QTL_SMA


QTL.result <- qtl.analysis (
  crossobj=DM.SO.cross.3
  step=0
  method='SIM'
  trait="em.k" 
  threshold="Li&Ji" 
  distance=30 
  cofactors=NULL
  window.size=30
  




mq.diagnostics (DM.SO.cross,I.threshold=0.1,
                p.val=0.01,na.cutoff=0.1)

test1 <- mstmap(cross.data , bychr = FALSE, 
                dist.fun = "kosambi",
                id = "Genotype",trace = FALSE)


pull.map(test1)

Examples

data(fake.f2)

newmap <- est.map(cross.data)
logliks <- sapply(newmap, attr, "loglik")
plotMap(fake.f2, newmap)
fake.f2 <- replace.map(fake.f2, newmap)



## Not run: 
QTL_SMA
QTL.result <- qtl.analysis (crossobj=cross.data,step=0,
                            method='SIM', trait="em.tm.nls", threshold="Li&Ji", distance=30, cofactors=NULL,
                            window.size=30)


## End(Not run)

# QTL_SIM
QTL.result <- qtl.analysis ( crossobj=cross.data, step=5,
                             method='SIM',trait="height", threshold="Li&Ji",
                             distance=30,cofactors=NULL,window.size=30)

# QTL CIM
cofactors <- as.vector (QTL.result$selected$marker)

QTL.result <- qtl.analysis ( crossobj=cross.data, step=5,
                             method='CIM', trait="height", threshold="Li&Ji", distance=30,
                             cofactors=cofactors, window.size=30)


qtl.cross



##### loop para datos faltantes
na_count <- apply(G.soja.prot.data, 1, function(x) sum(is.na(x)))
porcentaje <- 70 # menos del % de missing data

hp1 <- (porcentaje * (ncol(G.soja.prot.data)-1))/ 100
hp2 <- NULL
hp3 <- NULL # son los que quedan 
length(hp3)

for(i in 1:length(na_count)){
  if (na_count [i] < hp1 ){
    hp2 <- cbind(hp2, na_count[i])
    hp3 <- cbind (hp3, i)
  }
}

G.soja.prot.data.1 <- G.soja.prot.data [hp3, ]
dim(G.soja.prot.data.1)
G.data.1 <- G.soja.prot.data.1

G.data.1 [is.na(G.data.1)] <- "-"
dim(G.data.1)
str(G.data.1)


id.g1 <-G.data.1$Genotype
id.g <-G.data$Genotype

(id.g %in% id.g1)
setdiff(id.g, id.g1)

id.diff <-setdiff(id.g, id.g1)

write.table(id.diff, file = "./Data/procdata/id.diff.txt"
            ,append = FALSE, quote = TRUE, 
            sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = TRUE,
            col.names = TRUE)

# Cargar los datos lmem.gwaser 
# este es la matriz con los ambientes completos
cross.data.1 <- gwas.cross (P.data, G.data.1, map.data,
                          cross='gwas', heterozygotes=FALSE)
summary (cross.data.1)

plotMissing (cross.data.1)

geno.image (cross.data.1, reorder=T,
            alternate.chrid=FALSE)

plotMap (cross.data.1)


# genero las matrices por ambiente
cross.data.E1 <- gwas.cross (P.E1.data, G.data.1, map.data,
                            cross='gwas', heterozygotes=FALSE)

plotPheno(cross.data.E1, pheno.col=3)


cross.data.E2 <- gwas.cross (P.E2.data, G.data.1, map.data,
                             cross='gwas', heterozygotes=FALSE)

plotPheno(cross.data.E2, pheno.col=3)


cross.data.E3 <- gwas.cross (P.E3.data, G.data.1, map.data,
                             cross='gwas', heterozygotes=FALSE)

plotPheno(cross.data.E3, pheno.col=3)

summary (cross.data.E1)
summary (cross.data.E2)
summary (cross.data.E3)

# Marker Quality
mq.g.diagnostics (crossobj=cross.data.E1, I.threshold=0.1, p.val=0.01, na.cutoff=0.1)

#PCA
pca.E1 <- pca.analysis (crossobj=cross.data.E1, p.val=0.05)

#pca1 <- pca.analysis (crossobj=cross.data.1, p.val=0.05)

# LD.plots
linkdis.plots(crossobj = cross.data.E1, heterozygotes = FALSE, chr = c('1'))

### ambiente E1  ########################
# Mixed model: Eigenanalysis (PCA as random component)
(pca.E1.GWAS <- gwas.analysis(crossobj=cross.data.E1, method="eigenstrat",
                            provide.K=FALSE, covariates=pca.E1$scores,
                            trait="proteina", threshold="Li&Ji",
                            p=0.05, 
                            out.file="GWAS PCA as Random model"))$selected

#Mixed model: Kinship model
(k.E1.GWAS <- gwas.analysis (crossobj=cross.data.E1, method="kinship",
                          provide.K=FALSE, covariates=FALSE, trait="proteina",
                          threshold="Li&Ji", p=0.05, out.file ="GWAS K as Random model "))$selected

# Naive
(naive.E1.GWAS <- gwas.analysis(crossobj=cross.data.E1, method="naive",
                             provide.K=FALSE, covariates=FALSE, 
                             trait="proteina", threshold="Li&Ji",
                             p=0.05, out.file="GWAS naive model"))$selected


# Mixed model: Q+K
(qk.E1.GWAS <- gwas.analysis (crossobj=cross.data.E1, method="QK", provide.K=FALSE,
                           covariates=pca.E1$scores, trait="proteina", 
                           threshold="Li&Ji", p=0.05,
                           out.file="GWAS Q + K model"))$selected

##### ambiente E2 ##############
#PCA
pca.E2 <- pca.analysis (crossobj=cross.data.E2, p.val=0.05)


# Mixed model: Eigenanalysis (PCA as random component)
(pca.E2.GWAS <- gwas.analysis(crossobj=cross.data.E2, method="eigenstrat",
                              provide.K=FALSE, covariates=pca.E2$scores,
                              trait="proteina", threshold="Li&Ji",
                              p=0.05, 
                              out.file="GWAS PCA as Random model"))$selected

#Mixed model: Kinship model
(k.E2.GWAS <- gwas.analysis (crossobj=cross.data.E2, method="kinship",
                             provide.K=FALSE, covariates=FALSE, trait="proteina",
                             threshold="Li&Ji", p=0.05, out.file ="GWAS K as Random model "))$selected

# Naive
(naive.E2.GWAS <- gwas.analysis(crossobj=cross.data.E2, method="naive",
                                provide.K=FALSE, covariates=FALSE, 
                                trait="proteina", threshold="Li&Ji",
                                p=0.05, out.file="GWAS naive model"))$selected


# Mixed model: Q+K
(qk.E2.GWAS <- gwas.analysis (crossobj=cross.data.E2, method="QK", provide.K=FALSE,
                              covariates=pca.E2$scores, trait="proteina", 
                              threshold="Li&Ji", p=0.05,
                              out.file="GWAS Q + K model"))$selected

##### ambiente E3 ##############
pca.E3 <- pca.analysis (crossobj=cross.data.E3, p.val=0.05)

# Mixed model: Eigenanalysis (PCA as random component)
(pca.E3.GWAS <- gwas.analysis(crossobj=cross.data.E3, method="eigenstrat",
                              provide.K=FALSE, covariates=pca.E3$scores,
                              trait="proteina", threshold="FDR",
                              p=0.05, 
                              out.file="GWAS PCA as Random model"))$selected

#Mixed model: Kinship model
(k.E3.GWAS <- gwas.analysis (crossobj=cross.data.E3, method="kinship",
                             provide.K=FALSE, covariates=FALSE, trait="proteina",
                             threshold="Li&Ji", p=0.05, out.file ="GWAS K as Random model "))$selected

# Naive
(naive.E3.GWAS <- gwas.analysis(crossobj=cross.data.E3, method="naive",
                                provide.K=FALSE, covariates=FALSE, 
                                trait="proteina", threshold="Li&Ji",
                                p=0.05, out.file="GWAS naive model"))$selected


# Mixed model: Q+K
(qk.E3.GWAS <- gwas.analysis (crossobj=cross.data.E3, method="QK", provide.K=FALSE,
                              covariates=pca.E3$scores, trait="proteina", 
                              threshold="Li&Ji", p=0.05,
                              out.file="GWAS Q + K model"))$selected
