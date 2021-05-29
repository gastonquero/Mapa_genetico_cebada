##############################################################################
# analis de GWAS tesis 
# datos FPTA 
# Gaston Quero - Lorena Cammarota
# 21-6-2019 
##############################################
getwd()
#setwd("D:/Documentos/1-Material para la maestría")
setwd("C:/Users/Usuario/Dropbox/Tesis_Lorena")



library("Matrix")
library("RColorBrewer")
library(lme4)
library(nlme)
library(car)
library("ggplot2")       
library("lattice")
library("latticeExtra")
library(multcompView)
library(dplyr)
library(xtable)
library(tidyverse)
library (emmeans)
library("qtl")
library(stringr)
library(data.table)
library(svMisc)
library(ggpubr)
library("ggsci")
library("FactoMineR")
library("factoextra")
library("corrplot")
library(gameofthrones)
#library(lmerTest)
# Cargar Paquetes
# Paquetes 
library(lme4)
#library(lmerTest)
library(nlme)
library(car)
library("ggplot2")       
library("lattice")
library("latticeExtra")
library(multcompView)
library(dplyr)
library(xtable)
library(tidyverse)
library (emmeans)
#library("qtl")
library(stringr)
library(data.table)
library(svMisc)
#library(mmQTLv0.8)
library(LDheatmap)
library(genetics)
library(qqman)
library("qtl")


EEMAC07.cross <- read.cross (format="csv",
                            dir="./Data/rawdata", file= "cross_EEMAC_07.csv",
                            na.strings="NA", 
                            genotypes=c("0","1"),
                            #alleles=c("0","1"),
                            estimate.map=FALSE, 
                            convertXdata=TRUE, error.prob=0.0001)



EEMAC07.cross.chr3 <- subset(EEMAC07.cross,  chr="3")

summary(EEMAC07.cross.chr3 )


crossobj =EEMAC07.cross.chr3
heterozygotes = FALSE
chr= 3

data <- NULL
for (i in 1:length(chr)) {
  a <- paste("crossobj$geno$'", chr[i], "'$data", sep = "")
  p1 <- eval(parse(text = a))
  data <- cbind(data, p1)
}

if (heterozygotes == "FALSE") {
  data[data == 1] <- "A/A"
  data[data == 2] <- "B/B"
}

###### Check this part!!!!
if (heterozygotes == "TRUE") {
  data[data == 1] <- "A/A"
  data[data == 2] <- "B/B"
  data[data == 3] <- "A/B"
}

genos <- genotype(data[, 1])
for (i in 2:dim(data)[2]) {
  g <- genotype(data[, i])
  genos <- data.frame(genos, g)
}

#### add marker names calculate LD
ld <- LD(genos)

# plot LD heatmap
plot.hm <- LDheatmap (ld$"R^2", title="Pairwise LD_chr3",
                      color = colorRampPalette(c("red4", "red","orangered", "orange","yellow1",  "blue4"))(60))





#plot.hm.r <- LDheatmap(ld$"r", 
#                    color = colorRampPalette(c("red4", "red","orangered", "orange","yellow1",  "blue4"))(60))

#plot.hm.Dprim <- LDheatmap(ld$"D'", 
#                      color = colorRampPalette(c("red4", "red","orangered", "orange","yellow1",  "blue4"))(60))


LD.EEMAC.07.cross.chr3.matrix <- plot.hm$LDmatrix

write.table (LD.EEMAC.07.cross.chr3.matrix , 
             file = "./Data/procdata/LD.EEMAC.07.cross.chr3.matrix.txt",
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = TRUE,
             col.names = TRUE)

### 
EEMAC.07.cross.tibble <- read_delim (file="./Data/rawdata/cross_EEMAC_07.csv",
                                      delim=",", quote = "\"", 
                                      na = "-", col_names = TRUE)





head(EEMAC.07.cross.chr1.tibble)

EEMAC.07.cross.chr1.tibble <- EEMAC.07.cross.chr1.tibble %>%
  dplyr::select(-S)


EEMAC.07.cross.chr1.tibble_snp  <- EEMAC.07.cross.chr1.tibble %>%
  dplyr::select (starts_with("S"))

head (EEMAC.07.cross.chr1.tibble_snp)

EEMAC.07.cross.chr1.tibble_snp.1 <- EEMAC.07.cross.chr1.tibble_snp [-c(1,2),]


EEMAC.07.cross.chr1.tibble_bp <- EEMAC.07.cross.chr1.tibble_snp [2,]

mkr.chr1 <- colnames (EEMAC.07.cross.chr1.tibble_snp)

x.mk <- data.frame (mrks = mkr.chr1)

pos.chr1    <- EEMAC.07.cross.chr1.tibble_snp [2,]
pos.chr1.1  <- data.frame (pos =pos.chr1[1,])

pos.chr1.2  <- as.numeric(t(pos.chr1.1))

EEMAC.07.cross.chr1.tibble_bp.1 <- x.mk %>%
  dplyr::mutate (pos=pos.chr1.2)

### 
dt <- EEMAC.07.cross.chr1.tibble_bp.1 
list.pos <- (unique (dt$pos))
list.mrks <- (unique (dt$mrks))

Z <- matrix (0, nrow=length(list.pos), ncol=length(list.pos)) 
Z <- as.data.frame(Z)
colnames(Z) <-list.mrks
rownames(Z) <-list.mrks

#### aca empieza la distancia en bp
x.bp <- lapply (list.pos, function (filtro.x1) { 
  lapply (list.pos, function (filtro.x2) {
    dt.x1 <- dt %>% 
      dplyr::filter (pos == filtro.x1)
    x1    <- dt.x1 [,2]
    id.x1 <- dt.x1 [,1]
    dt.x2 <- dt %>% 
      dplyr::filter (pos == filtro.x2)
    x2    <- dt.x2 [,2]
    id.x2 <- dt.x2 [,1]
    dt.z <- abs (x1 - x2)
    Z ["id.x1","id.x2"] <- dt.z
  })
}) 

#View (x.bp)
names(x.bp) <- list.mrks
XX.x.bp <- as.data.frame (do.call (cbind, x.bp))
is.list(XX.x.bp)

dt.LD.decay <- lapply (list.mrks, function (filtro) {
  XX.x.bp.1 <- XX.x.bp %>% 
    dplyr::select (filtro) 
  delta.bp <- unlist (XX.x.bp.1 , use.names=FALSE)
  
  colnames (LD.EEMAC.07.cross.chr1.matrix) <- list.mrks
  rownames (LD.EEMAC.07.cross.chr1.matrix) <- list.mrks
  LD_EEMAC.07.cross.chr1 <-  as_tibble(LD.EEMAC.07.cross.chr1.matrix)
  LD_EEMAC.07.cross.chr1 <- LD_EEMAC.07.cross.chr1 %>%
    dplyr::mutate (mrk.id =list.mrks)%>%
    dplyr::select (mrk.id, everything())
  
  id <- colnames (XX.x.bp.1)
  
  LD_EEMAC.07.cross.chr1.1 <- LD_EEMAC.07.cross.chr1%>%
    dplyr::filter (mrk.id==id)
  
  id.gather <- colnames (LD_EEMAC.07.cross.chr1.1) [-1]
  
  LD_EEMAC.07.cross.chr1.2 <- LD_EEMAC.07.cross.chr1.1 %>%
    gather (id.gather, key="mrk.2" , value = "R2") %>%
    dplyr::rename (mrk.1 = mrk.id) 
  
  LD_EEMAC.07.cross.chr1.2$mrk.1 <- as.character(LD_EEMAC.07.cross.chr1.2$mrk.1)
  
  
  colnames (XX.x.bp.1) == unique ( LD_EEMAC.07.cross.chr1.2$mrk.1)
  
  LD_EEMAC.07.cross.chr1.3 <- LD_EEMAC.07.cross.chr1.2 %>%
    dplyr::mutate (delta.bp = delta.bp )
  
  print (LD_EEMAC.07.cross.chr1.3)
  return ( LD_EEMAC.07.cross.chr1.3)
  
})

plot.LD.decay.chr1 <- as.data.frame (do.call (rbind, dt.LD.decay))
unique(plot.LD.decay.chr1$mrk.1)


write_delim (plot.LD.decay.chr1 , path="./Data/procdata/plot.LD.decay.chr1.R2.txt",
             delim = ",", na = "NA")


#########################################














summary (EEMAC07.cross)

write.cross (EEMAC07.cross , format="csv",
            filestem="./Data/rawdata/chr.3.LD", chr=3, digits=NULL)


EEMAC07.chr3.cross <- read.cross (format="csv",
                             dir="./Data/rawdata", file= "chr.3.LD.csv",
                             na.strings="NA", 
                             genotypes=c("0","1"),
                             #alleles=c("0","1"),
                             estimate.map=FALSE, 
                             convertXdata=TRUE, error.prob=0.0001)



chr.3.LD <- pull.map(EEMAC07.cross, chr=1, as.table=FALSE)



M_1111223 
M_1110443
#############33


chr.1.LD.cross.tibble <-  read_delim ("./Data/rawdata/chr.1.LD.csv", 
                                     delim = ",", na = "NA", col_types = "f")


LD.chr1 <- EEMAC07.cross.tibble %>%
                           dplyr::select(c(genotypes,M_1111223:M_1110443))


write_delim (LD.qtl.chr4.kinship.EX.1, path="./Data/procdata/LD.qtl.chr4.kinship.EX.1.csv", 
             delim = ",", na = " ", 
             append = FALSE, quote_escape = "double")

qtl.chr4.kinship.EX.1.LD.cross <- read.cross (format="csv",
                                     dir="./Data/procdata", file= "LD.qtl.chr4.kinship.EX.1.csv",
                                     na.strings="NA", 
                                     genotypes=c("0","1"),
                                     crosstype = 'dh',
                                     #alleles=c("0","1"),
                                     estimate.map=FALSE,
                                     convertXdata=TRUE, error.prob=0.0001)

file = qtl.chr4.kinship.EX.1.LD.cross

#construct matrix in correct format
data <- NULL
for(i in 4:nchr(file)){
  a<-paste("file$geno$'", i, "'$data", sep="")
  p1<-eval(parse(text=a))
  data <-cbind(data, p1)
}
dim(data)
data[data==1]<-"A/A"
data[data==2]<-"B/B"

genos <- genotype(data[,1])
for(i in 2:dim(data)[2]){
  g <-genotype(data[,i])
  genos <-data.frame(genos,g)
}

names(genos) <- markernames (qtl.chr4.kinship.EX.1.LD.cross)

####add marker names

#calculate LD
ld <- LD(genos)

#ld$"D'"
rgb.palette <- colorRampPalette(rev(c("blue", "orange", "red")), space = "rgb")


#plot LD heatmap

map.chr4 <- pull.map(qtl.chr4.kinship.EX.1.LD.cross, as.table=TRUE)






LD.qtl.chr4.kinship.EX.1 <-LDheatmap(ld$"D'",
                                color= rgb.palette(50),
                                #color=colorRampPalette(c("red","yellow"))(50),
                                add.map = TRUE,add.key=TRUE,geneMapLocation = 0.1,
                                pop=TRUE,flip=TRUE,
                                genetic.distances= map.chr4$pos,
                                SNP.name=c("M_1120145","M_1121228", "M_1110490"), 
                                title="qtl.chr4.kinship.EX.1")



LD.qtl.chr4.kinship.EX.1.matrix  <- LD.qtl.chr4.kinship.EX.1$LDmatrix

write.table (LD.qtl.chr4.kinship.EX.1.matrix, 
             file = "./Data/procdata/LD.qtl.chr4.kinship.EX.1.matrix.txt",
             append = FALSE, quote = TRUE, sep = ",",
             eol = "\n", na = "NA", dec = ".", row.names = TRUE,
             col.names = TRUE)




LD.qtl.chr4.kinship.EX.1.dt <- as.data.frame (LD.qtl.chr4.kinship.EX.1.matrix)


LD.qtl.chr4.kinship.EX.1.dt <-  LD.qtl.chr4.kinship.EX.1.dt %>%
                               dplyr::mutate (mrks = row.names (LD.qtl.chr4.kinship.EX.1.dt)) %>%
                               dplyr::select (mrks, everything())


r.markers.max <- LD.qtl.chr4.kinship.EX.1.dt %>%
                dplyr::filter(mrks == "M_1121228")


#### cargar los datos segun R/qtl
#en este caso solo se usan los datos Geno, NO  feno
EEMAC07.cross <- read.cross(format="csv",
                           dir="./Data/rawdata", file= "cross_EEMAC_07.csv",
                           na.strings="NA", 
                           genotypes=c("0","1"),
                           #alleles=c("0","1"),
                           estimate.map=FALSE, convertXdata=TRUE, error.prob=0.0001)

# Pide la función jittermap porque hay marcadores con la misma localización (cM)
EEMAC07.cross <- jittermap(EEMAC07.cross)
summary(EEMAC07.cross)
summary.map (EEMAC07.cross)
plotMap(EEMAC07.cross)
#plot(EEMAC07.cross)
plotMissing(EEMAC07.cross)

crossobj = EEMAC07.cross
I.threshold = 0.1
I.quant = FALSE
p.val = 0.01
na.cutoff = 0.1

#mq_comparegenotypes_plot <- function(crossobj) {
  jittermap(crossobj)
  par(mfrow = c(1, 1))
  output <- comparegeno (crossobj)
  n.ind <- nind (crossobj)
  
  pal <- got(10, alpha = 1, option = "Targaryen2", direction = -1)
  pal1 <- got(100, option = "Stark")
  pal12 <- got(100, option = "Tully")
  pal2 <- got(100, option = "Daenerys", direction = 1)
  pal3 <- got(250, option = "Lannister", direction = -1)
  pal4 <- got(100, option = "Martell", direction = 1)

#Permite acceder a la matriz de propocrción
#comparegeno(EEMAC07.cross, what=c("proportion","number","both"))
              
  image(1:n.ind, 1:n.ind, output,
        col = pal3,
        main = "Pairwise comparation of genotypes")
  box()

  image(1:n.ind, 1:n.ind, output,
        main = "Pairwise comparation of genotypes")
  
  image(1:n.ind, 1:n.ind, output,
        col = gray((0:99) / 99), breaks = seq(0, 1, len = 101),
        main = "Pairwise comparation of genotypes")

  #devtools::install_github("aljrico/gameofthrones")
  #install.packages("gameofthrones")
  

## Esto hay que hacerlo con una codigo antes del r/qtl
# MAF se sacan los marcadores con un frecuencia menor al 5% y mayor al 95%

total <- nind (EEMAC07.cross)
gt <- geno.table (EEMAC07.cross)
missing <- gt[,2]
num.obs <- c(total - missing)
obs.alelo <- gt[, c(3,4)]
Geno.freq <- as.vector(obs.alelo/num.obs)

#saco los marcadores que sacó con frecuencia menores a 5%, no saco los mayores a 95%
SNPs.MAF <-  Geno.freq[,2] < 0.05 
rownames(Geno.freq [SNPs.MAF,])
#no hay marcador con frecuancia menor
names.marker <- c (rownames(Geno.freq [SNPs.MAF,]))
length(names.marker)
#control de cuanto saqué, en este caso 0

#quito los marcadores anteriores, nada en este caso
EEMAC07.cross.1 <- drop.markers (EEMAC07.cross, names.marker)

summary (EEMAC07.cross.1)
# si la longitud da 20, el drop markers me tiene que sacar los marcadores mencionados.
plotMissing (EEMAC07.cross.1 )

geno.image (EEMAC07.cross.1 , reorder=FALSE, main="Genotype data_FPTA",
            alternate.chrid=FALSE)

mq.g.diagnostics (EEMAC07.cross.1, I.threshold = 0.1,
                  I.quant = FALSE, p.val = 0.01, na.cutoff = 0.1) 

#PCA Geno, usa todos los datos de los marcadores, SOLO
pca.geno <- pca.analysis (crossobj=EEMAC07.cross.1, p.val=0.05)
pca.x<-pca.geno$scores
class(pca.x)
pca.x<-as.data.frame(pca.x)
write_delim(pca.x, path = "./Data/procdata/pca_matrix.txt", delim = ",", 
            na = "NA", append = FALSE,
             quote_escape = "double")

class(EEMAC07.cross.1)

##### ACA COMIENZO A USAR LOS DATOS FENOTÍPICOS QUE ESTÁN EN LA MATRIZ EEMAC07.cross.1###
# lo que hace comienza a usar los datos Pheno, columna EX, una variable por vez, localidad EEMAC07
plotPheno (EEMAC07.cross.1, pheno.col="EX")

#Mixed model: Eigenanalysis (PCA as random component)
(pcaR.GWAS <- gwas.analysis (crossobj= EEMAC07.cross.1,
                             method="eigenstrat",
                             provide.K=FALSE, 
                             covariates=pca.geno$scores, 
                             trait="EX",
                             threshold="Li&Ji", 
                             p=0.05, 
                             out.file="GWAS PCA as Random model"))$selected

#'#Mixed model: Kinship model
(k.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                         method="kinship",
                         provide.K=FALSE, covariates=FALSE, 
                         trait="EX",
                         threshold="Li&Ji", p=0.05,
                         out.file =" GWAS K as Random model "))$selected

#'#Mixed model: Q+K
(qk.GWAS <- gwas.analysis (crossobj=EEMAC07.cross.1, 
                           method="QK", 
                           provide.K=FALSE,
                           covariates=pca.geno$scores, 
                           trait="EX", 
                           threshold="Li&Ji", p=0.05,
                           out.file="GWAS Q + K model"))$selected

# Naive
(naive.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                             method="naive",
                             provide.K=FALSE,
                             covariates=FALSE, 
                             trait="EX", 
                             threshold="Li&Ji",
                             p=0.05, 
                             out.file="GWAS naive model"))$selected

######PM######
plotPheno (EEMAC07.cross.1, pheno.col="PM")



#Mixed model: Eigenanalysis (PCA as random component)
(pcaR.GWAS <- gwas.analysis (crossobj= EEMAC07.cross.1,
                             method="eigenstrat",
                             provide.K=FALSE, 
                             covariates=pca.geno$scores, 
                             trait="PM",
                             threshold="Li&Ji", 
                             p=0.05, 
                             out.file="GWAS PCA as Random model"))$selected

#'#Mixed model: Kinship model
(k.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                         method="kinship",
                         provide.K=FALSE, covariates=FALSE, 
                         trait="PM",
                         threshold="Li&Ji", p=0.05,
                         out.file =" GWAS K as Random model "))$selected

#'#Mixed model: Q+K
(qk.GWAS <- gwas.analysis (crossobj=EEMAC07.cross.1, 
                           method="QK", 
                           provide.K=FALSE,
                           covariates=pca.geno$scores, 
                           trait="PM", 
                           threshold="Li&Ji", p=0.05,
                           out.file="GWAS Q + K model"))$selected

# Naive
(naive.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                             method="naive",
                             provide.K=FALSE,
                             covariates=FALSE, 
                             trait="PM", 
                             threshold="Li&Ji",
                             p=0.05, 
                             out.file="GWAS naive model"))$selected

######NS######
plotPheno (EEMAC07.cross.1, pheno.col="NS")

#Mixed model: Eigenanalysis (PCA as random component)
(pcaR.GWAS <- gwas.analysis (crossobj= EEMAC07.cross.1,
                             method="eigenstrat",
                             provide.K=FALSE, 
                             covariates=pca.geno$scores, 
                             trait="NS",
                             threshold="Li&Ji", 
                             p=0.05, 
                             out.file="GWAS PCA as Random model"))$selected

#'#Mixed model: Kinship model
(k.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                         method="kinship",
                         provide.K=FALSE, covariates=FALSE, 
                         trait="NS",
                         threshold="Li&Ji", p=0.05,
                         out.file =" GWAS K as Random model "))$selected

#'#Mixed model: Q+K
(qk.GWAS <- gwas.analysis (crossobj=EEMAC07.cross.1, 
                           method="QK", 
                           provide.K=FALSE,
                           covariates=pca.geno$scores, 
                           trait="NS", 
                           threshold="Li&Ji", p=0.05,
                           out.file="GWAS Q + K model"))$selected

# Naive
(naive.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                             method="naive",
                             provide.K=FALSE,
                             covariates=FALSE, 
                             trait="NS", 
                             threshold="Li&Ji",
                             p=0.05, 
                             out.file="GWAS naive model"))$selected

########IK########
plotPheno (EEMAC07.cross.1, pheno.col="IK")

#Mixed model: Eigenanalysis (PCA as random component)
(pcaR.GWAS <- gwas.analysis (crossobj= EEMAC07.cross.1,
                             method="eigenstrat",
                             provide.K=FALSE, 
                             covariates=pca.geno$scores, 
                             trait="IK",
                             threshold="Li&Ji", 
                             p=0.05, 
                             out.file="GWAS PCA as Random model"))$selected

#'#Mixed model: Kinship model
(k.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                         method="kinship",
                         provide.K=FALSE, covariates=FALSE, 
                         trait="IK",
                         threshold="Li&Ji", p=0.05,
                         out.file =" GWAS K as Random model "))$selected

#'#Mixed model: Q+K
(qk.GWAS <- gwas.analysis (crossobj=EEMAC07.cross.1, 
                           method="QK", 
                           provide.K=FALSE,
                           covariates=pca.geno$scores, 
                           trait="IK", 
                           threshold="Li&Ji", p=0.05,
                           out.file="GWAS Q + K model"))$selected

# Naive
(naive.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                             method="naive",
                             provide.K=FALSE,
                             covariates=FALSE, 
                             trait="IK", 
                             threshold="Li&Ji",
                             p=0.05, 
                             out.file="GWAS naive model"))$selected

########FR########
plotPheno (EEMAC07.cross.1, pheno.col="FR")

#Mixed model: Eigenanalysis (PCA as random component)
(pcaR.GWAS <- gwas.analysis (crossobj= EEMAC07.cross.1,
                             method="eigenstrat",
                             provide.K=FALSE, 
                             covariates=pca.geno$scores, 
                             trait="FR",
                             threshold="Li&Ji", 
                             p=0.05, 
                             out.file="GWAS PCA as Random model"))$selected

#'#Mixed model: Kinship model
(k.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                         method="kinship",
                         provide.K=FALSE, covariates=FALSE, 
                         trait="FR",
                         threshold="Li&Ji", p=0.05,
                         out.file =" GWAS K as Random model "))$selected

#'#Mixed model: Q+K
(qk.GWAS <- gwas.analysis (crossobj=EEMAC07.cross.1, 
                           method="QK", 
                           provide.K=FALSE,
                           covariates=pca.geno$scores, 
                           trait="FR", 
                           threshold="Li&Ji", p=0.05,
                           out.file="GWAS Q + K model"))$selected

# Naive
(naive.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                             method="naive",
                             provide.K=FALSE,
                             covariates=FALSE, 
                             trait="FR", 
                             threshold="Li&Ji",
                             p=0.05, 
                             out.file="GWAS naive model"))$selected

########VIS#######
plotPheno (EEMAC07.cross.1, pheno.col="VIS")



#Mixed model: Eigenanalysis (PCA as random component)
(pcaR.GWAS <- gwas.analysis (crossobj= EEMAC07.cross.1,
                             method="eigenstrat",
                             provide.K=FALSE, 
                             covariates=pca.geno$scores, 
                             trait="VIS",
                             threshold="Li&Ji", 
                             p=0.05, 
                             out.file="GWAS PCA as Random model"))$selected

#'#Mixed model: Kinship model
(k.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                         method="kinship",
                         provide.K=FALSE, covariates=FALSE, 
                         trait="VIS",
                         threshold="Li&Ji", p=0.05,
                         out.file =" GWAS K as Random model "))$selected

#'#Mixed model: Q+K
(qk.GWAS <- gwas.analysis (crossobj=EEMAC07.cross.1, 
                           method="QK", 
                           provide.K=FALSE,
                           covariates=pca.geno$scores, 
                           trait="VIS", 
                           threshold="Li&Ji", p=0.05,
                           out.file="GWAS Q + K model"))$selected

# Naive
(naive.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                             method="naive",
                             provide.K=FALSE,
                             covariates=FALSE, 
                             trait="VIS", 
                             threshold="Li&Ji",
                             p=0.05, 
                             out.file="GWAS naive model"))$selected

#########PD######
plotPheno (EEMAC07.cross.1, pheno.col="PD")



#Mixed model: Eigenanalysis (PCA as random component)
(pcaR.GWAS <- gwas.analysis (crossobj= EEMAC07.cross.1,
                             method="eigenstrat",
                             provide.K=FALSE, 
                             covariates=pca.geno$scores, 
                             trait="PD",
                             threshold="Li&Ji", 
                             p=0.05, 
                             out.file="GWAS PCA as Random model"))$selected

#'#Mixed model: Kinship model
(k.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                         method="kinship",
                         provide.K=FALSE, covariates=FALSE, 
                         trait="PD",
                         threshold="Li&Ji", p=0.05,
                         out.file =" GWAS K as Random model "))$selected

#'#Mixed model: Q+K
(qk.GWAS <- gwas.analysis (crossobj=EEMAC07.cross.1, 
                           method="QK", 
                           provide.K=FALSE,
                           covariates=pca.geno$scores, 
                           trait="PD", 
                           threshold="Li&Ji", p=0.05,
                           out.file="GWAS Q + K model"))$selected

# Naive
(naive.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                             method="naive",
                             provide.K=FALSE,
                             covariates=FALSE, 
                             trait="PD", 
                             threshold="Li&Ji",
                             p=0.05, 
                             out.file="GWAS naive model"))$selected

########BG#####
plotPheno (EEMAC07.cross.1, pheno.col="BG")



#Mixed model: Eigenanalysis (PCA as random component)
(pcaR.GWAS <- gwas.analysis (crossobj= EEMAC07.cross.1,
                             method="eigenstrat",
                             provide.K=FALSE, 
                             covariates=pca.geno$scores, 
                             trait="BG",
                             threshold="Li&Ji", 
                             p=0.05, 
                             out.file="GWAS PCA as Random model"))$selected

#'#Mixed model: Kinship model
(k.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                         method="kinship",
                         provide.K=FALSE, covariates=FALSE, 
                         trait="BG",
                         threshold="Li&Ji", p=0.05,
                         out.file =" GWAS K as Random model "))$selected

#'#Mixed model: Q+K
(qk.GWAS <- gwas.analysis (crossobj=EEMAC07.cross.1, 
                           method="QK", 
                           provide.K=FALSE,
                           covariates=pca.geno$scores, 
                           trait="BG", 
                           threshold="Li&Ji", p=0.05,
                           out.file="GWAS Q + K model"))$selected

# Naive
(naive.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                             method="naive",
                             provide.K=FALSE,
                             covariates=FALSE, 
                             trait="BG", 
                             threshold="Li&Ji",
                             p=0.05, 
                             out.file="GWAS naive model"))$selected

##########NO ATL#######
plotPheno (EEMAC07.cross.1, pheno.col="ATL")



#Mixed model: Eigenanalysis (PCA as random component)
(pcaR.GWAS <- gwas.analysis (crossobj= EEMAC07.cross.1,
                             method="eigenstrat",
                             provide.K=FALSE, 
                             covariates=pca.geno$scores, 
                             trait="ATL",
                             threshold="Li&Ji", 
                             p=0.05, 
                             out.file="GWAS PCA as Random model"))$selected

#'#Mixed model: Kinship model
(k.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                         method="kinship",
                         provide.K=FALSE, covariates=FALSE, 
                         trait="ATL",
                         threshold="Li&Ji", p=0.05,
                         out.file =" GWAS K as Random model "))$selected

#'#Mixed model: Q+K
(qk.GWAS <- gwas.analysis (crossobj=EEMAC07.cross.1, 
                           method="QK", 
                           provide.K=FALSE,
                           covariates=pca.geno$scores, 
                           trait="ATL", 
                           threshold="Li&Ji", p=0.05,
                           out.file="GWAS Q + K model"))$selected

# Naive
(naive.GWAS <- gwas.analysis(crossobj=EEMAC07.cross.1, 
                             method="naive",
                             provide.K=FALSE,
                             covariates=FALSE, 
                             trait="ATL", 
                             threshold="Li&Ji",
                             p=0.05, 
                             out.file="GWAS naive model"))$selected

