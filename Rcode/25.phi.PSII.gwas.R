###################################################################
#  codigo para hacer GWAS en la subpoblacion de indica            #
# los trait son de fotosintesis                                   #
# Gaston Quero - Sebastian Simondi                                #
# 3/4/2019                                                        #
###################################################################

# fijar Directorio
getwd()
setwd("E:/Work/gquero/gwas.photo.arroz/phiPSII")


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
library("qtl")
library(stringr)
library(data.table)
library(svMisc)

# cargo los datos

gwas.photo.indica <- read.cross (format="csv",
                                dir="./Data/rawdata", 
                                file= "gwas.photo.indica.1a.csv",
                                na.strings="NA", 
                                genotypes=c("0","1"),
                                crosstype="dh",
                                #alleles=c("0","1"),
                                estimate.map=FALSE,
                                convertXdata=TRUE, error.prob=0.0001)


## Esto hay que hacerlo con una codigo antes del r/qtl
# MAF se sacan los marcadores con un frecuencia menor al 5% y mayor al 95

total <- nind (gwas.photo.indica)
gt <- geno.table (gwas.photo.indica)
missing <- gt[,2]
num.obs <- c(total - missing)
obs.alelo <- gt[, c(3,4)]
Geno.freq <- as.vector(obs.alelo/num.obs)

SNPs.MAF <-  Geno.freq[,2] < 0.05 
rownames(Geno.freq [SNPs.MAF,])
names.marker <- c (rownames(Geno.freq [SNPs.MAF,]))
length(names.marker)

gwas.photo.indica.1 <- drop.markers (gwas.photo.indica, names.marker)

#PCA
pca.phiPSII <- pca.analysis (crossobj=gwas.photo.indica.1, p.val=0.05)

str()

plotPheno(gwas.photo.indica.1, pheno.col="em.phiII")

# Exportar los datos 
write.cross (gwas.photo.indica.1, format=c("csv"),
             filestem="./Data/procdata/gwas.photo.indica.1")

#Mixed model: Eigenanalysis (PCA as random component)
(pcaR.GWAS <- gwas.analysis (crossobj= gwas.photo.indica.1,
                             method="eigenstrat",
                             provide.K=FALSE, 
                             covariates=pca.phiPSII$scores, 
                             trait="em.phiII",
                             threshold="Li&Ji", 
                             p=0.05, 
                             out.file="GWAS PCA as Random model"))$selected

#'#Mixed model: Kinship model
(k.GWAS <- gwas.analysis(crossobj=gwas.photo.indica.1, 
                         method="kinship",
                         provide.K=FALSE, covariates=FALSE, 
                         trait="em.phiII",
                         threshold="Li&Ji", p=0.05,
                         out.file =" GWAS K as Random model "))$selected

#'#Mixed model: Q+K
(qk.GWAS <- gwas.analysis (crossobj=gwas.photo.indica.1, 
                           method="QK", 
                           provide.K=FALSE,
                           covariates=pca.phiPSII$scores, 
                           trait="em.phiII", 
                           threshold="Li&Ji", p=0.05,
                           out.file="GWAS Q + K model"))$selected

# Naive
(naive.GWAS <- gwas.analysis(crossobj=gwas.photo.indica.1, 
                             method="naive",
                             provide.K=FALSE,
                             covariates=FALSE, 
                             trait="em.phiII", 
                             threshold="Li&Ji",
                             p=0.05, 
                             out.file="GWAS naive model"))$selected



chr.1 <- linkdis.plots (crossobj=gwas.photo.indica.1,
                         heterozygotes=FALSE, chr=1)

chr.2 <- linkdis.plots (crossobj=gwas.photo.indica.1,
                        heterozygotes=FALSE, chr=2)

chr.3 <- linkdis.plots (crossobj=gwas.photo.indica.1,
                        heterozygotes=FALSE, chr=3)


chr.4 <- linkdis.plots (crossobj=gwas.photo.indica.1,
                        heterozygotes=FALSE, chr=4)


chr.5 <- linkdis.plots (crossobj=gwas.photo.indica.1,
                        heterozygotes=FALSE, chr=5)

chr.6 <- linkdis.plots (crossobj=gwas.photo.indica.1,
                        heterozygotes=FALSE, chr=6)

chr.7 <- linkdis.plots (crossobj=gwas.photo.indica.1,
                        heterozygotes=FALSE, chr=7)

chr.8 <- linkdis.plots (crossobj=gwas.photo.indica.1,
                        heterozygotes=FALSE, chr=8)

chr.9 <- linkdis.plots (crossobj=gwas.photo.indica.1,
                        heterozygotes=FALSE, chr=9)

chr.10 <- linkdis.plots (crossobj=gwas.photo.indica.1,
                        heterozygotes=FALSE, chr=10)
chr.11 <- linkdis.plots (crossobj=gwas.photo.indica.1,
                         heterozygotes=FALSE, chr=11)

chr.12 <- linkdis.plots (crossobj=gwas.photo.indica.1,
                         heterozygotes=FALSE, chr=12)
