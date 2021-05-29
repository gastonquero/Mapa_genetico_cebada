
getwd()

setwd("C:/Users/Usuario/OneDrive/Documentos/Paper_Cebada")


# paquetes
library (dplyr)
library (xtable)
library (tidyverse)
library (emmeans)
library (qtl)
library (qtl2)
library (stringr)
library (data.table)
library (svMisc)
library (ggpubr)
library (ggsci)
library (FactoMineR)
library (factoextra)
library (corrplot)
library (onemap)
library (broman) 
library (ASMap)
library(lmem.qtler)

## vamos a cargar la matriz de bill con el mapa de bill
set.seed(5)

# cross Ill50K.2
bill.a_Mapeo.qtl <- read.cross (format="tidy",
                                dir="./Data/procdata",
                                genfile ="geno.bill.a.csv", 
                                mapfile ="pmap_bill.a.csv", 
                                phefile ="pheno_nico.2.csv" , 
                                na.strings="-",
                                genotypes=c("AA","BB"), 
                                alleles=c("A","B"),
                                estimate.map=FALSE, 
                                convertXdata=TRUE, error.prob=0.0001,
                                map.function="kosambi",
                                #F.gen=6, 
                                crosstype ="riself")

plotMap(bill.a_Mapeo.qtl)

bill.a_Mapeo.qtl$pheno
ci.e.a.gd

## QTL_SMA
QTL.result.ci.e.a.gd <- qtl.analysis (crossobj=bill.a_Mapeo.qtl,
                                      step=0, method='SIM',
                                      trait="ci.e.a.gd", 
                                      threshold="Li&Ji", 
                                      distance=30,  
                                      cofactors=NULL, 
                                      window.size=30)

### QTL_SIM
QTL.result.ci.e.a.gd  <- qtl.analysis ( crossobj=bill.a_Mapeo.qtl, 
                                        step=5,method='SIM',
                                        trait="ci.e.a.gd", 
                                        threshold="Li&Ji",
                                        distance=30,
                                        cofactors=NULL, 
                                        window.size=30)
### QTL CIM
cofactors <- as.vector (QTL.result.ci.e.a.gd$selected$marker)
QTL.result.ci.e.a.gd <- qtl.analysis ( crossobj=bill.a_Mapeo.qtl,
                                       step=5, method='CIM', 
                                       trait="ci.e.a.gd", 
                                       threshold="Li&Ji", 
                                       distance=30,
                                       cofactors=cofactors, 
                                       window.size=30)

####
# cross Ill50K.2
Mapeo.qtl.4 <- read.cross (format="tidy",
                                dir="./Data/procdata",
                                genfile ="geno.bill.4.csv", 
                                mapfile ="bill_Mapeo.qtl.4.a_map.csv", 
                                phefile ="pheno_nico.2.csv" , 
                                na.strings="-",
                                genotypes=c("AA","BB"), 
                                alleles=c("A","B"),
                                estimate.map=FALSE, 
                                convertXdata=TRUE, error.prob=0.0001,
                                map.function="kosambi",
                                #F.gen=6, 
                                crosstype ="riself")

Mapeo.qtl.4 <- jittermap(Mapeo.qtl.4)


## QTL_SMA
QTL.result.ci.e.a.gd.4 <- qtl.analysis (crossobj=Mapeo.qtl.4,
                                      step=0, method='SIM',
                                      trait="ci.e.a.gd", 
                                      threshold="Li&Ji", 
                                      distance=30,  
                                      cofactors=NULL, 
                                      window.size=30)

### QTL_SIM
QTL.result.ci.e.a.gd.4  <- qtl.analysis ( crossobj=Mapeo.qtl.4, 
                                        step=5,method='SIM',
                                        trait="ci.e.a.gd", 
                                        threshold="Li&Ji",
                                        distance=30,
                                        cofactors=NULL, 
                                        window.size=30)
### QTL CIM
cofactors <- as.vector (QTL.result.ci.e.a.gd.4$selected$marker)
QTL.result.ci.e.a.gd.4 <- qtl.analysis ( crossobj=Mapeo.qtl.4,
                                       step=5, method='CIM', 
                                       trait="ci.e.a.gd", 
                                       threshold="Li&Ji", 
                                       distance=30,
                                       cofactors=cofactors, 
                                       window.size=30)





bill.a_Mapeo.qtl$pheno
