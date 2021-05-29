##############################################################################
# analis de GWAS tesis 
# datos FPTA 
# Gaston Quero - Lorena Cammarota
# 21-6-2019 
##############################################
getwd()
setwd("C:/Users/Usuario/Dropbox/Tesis_Lorena")


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
library(ggpubr)
library("ggsci")
library("FactoMineR")
library("factoextra")
library("corrplot")




#### cargar los datos segun R/qtl
dt.EEMAC07<- read.cross(format="csv",
                           dir="./Data/rawdata", file= "cross_EEMAC_07.csv",
                           na.strings="NA", 
                           genotypes=c("0","1"),
                           #alleles=c("0","1"),
                           estimate.map=FALSE, convertXdata=TRUE, error.prob=0.0001)


plotMissing(dt.EEMAC07)

crossobj = dt.EEMAC07
I.threshold = 0.1
I.quant = FALSE
p.val = 0.01
na.cutoff = 0.1

#mq_comparegenotypes_plot <- function(crossobj) {
  jittermap(crossobj)
  par(mfrow = c(1, 1))
  output <- comparegeno (crossobj)
  n.ind <- nind (crossobj)
  image(1:n.ind, 1:n.ind, output,
        col = pal2,
        main = "Pairwise comparation of genotypes")
  box()
  
  
  image(1:n.ind, 1:n.ind, output,
        col = gray((0:99) / 99), breaks = seq(0, 1, len = 101),
        main = "Pairwise comparation of genotypes")
  
  
  palx <- col="red"
  #devtools::install_github("aljrico/gameofthrones")
  #install.packages("gameofthrones")
  library(gameofthrones)
  pal <- got(10, alpha = 1, option = "Targaryen2", direction = -1)
  pal1 <- got(100, option = "Stark")
  pal12 <- got(100, option = "Tully")
  pal3 <- got(250, option = "Lannister", direction = -1)
  pal4 <- got(100, option = "Martell", direction = 1)
  pal2 <- got(100, option = "Daenerys", direction = 1)
  
#}





## Esto hay que hacerlo con una codigo antes del r/qtl
# MAF se sacan los marcadores con un frecuencia menor al 5% y mayor al 95

total <- nind (dt.EEMAC07)
gt <- geno.table (dt.EEMAC07)
missing <- gt[,2]
num.obs <- c(total - missing)
obs.alelo <- gt[, c(3,4)]
Geno.freq <- as.vector(obs.alelo/num.obs)

SNPs.MAF <-  Geno.freq[,2] < 0.05 
rownames(Geno.freq [SNPs.MAF,])
names.marker <- c (rownames(Geno.freq [SNPs.MAF,]))
length(names.marker)

dt.EEMAC07.1 <- drop.markers (dt.EEMAC07, names.marker)

summary (dt.EEMAC07.1)

plotMissing (dt.EEMAC07.1 )

geno.image (dt.EEMAC07.1 , reorder=FALSE, main="Genotype data_FPTA",
            alternate.chrid=FALSE)


#PCA
pca.EX <- pca.analysis (crossobj=dt.EEMAC07.1, p.val=0.05)
class(dt.EEMAC07.1)
str()

plotPheno (dt.EEMAC07.1, pheno.col="EX")


mq.g.diagnostics (dt.EEMAC07.1, I.threshold = 0.1,
                             I.quant = FALSE, p.val = 0.01, na.cutoff = 0.1) 


### 



PCA.GENO.1 <- read_delim (file="./Data/rawdata/PCA.GENO.txt", 
                                        delim = "\t", na = "NA")

PCA.GENO.1b <- type_convert (PCA.GENO.1, 
                                          col_types =c(genotypo = col_character(),
                                                       .default = col_integer()))



class (PCA.EEMAC07.1)

res.PCA.GENO <- PCA (PCA.GENO.1b, scale.unit = TRUE, ncp = 5, 
                  ind.sup = NULL, 
                  quanti.sup = 0, quali.sup = 1, row.w = NULL, 
                  col.w = NULL, graph = TRUE, axes = c(1,2))


eig.val.1 <- get_eigenvalue (res.PCA.GENO )
eig.val.1

fviz_eig (res.PCA.GENO , addlabels = TRUE, ylim = c(0, 15))

##### individuos
ind <- get_pca_ind(res.PCA.GENO)
#ind
fviz_pca_ind(res.PCA.GENO)
#res.pca.1$quali
fviz_pca_ind (res.PCA.GENO,geom = c("point"),
              palette = "jco", repel = TRUE)

# Clustering, auto nb of clusters:
hc.geno <- HCPC (res.PCA.GENO, nb.clust=-1)

#hc$data.clust

clust.pca.geno <- hc.geno$data.clust 
class (clust.pca.geno)
str(clust.pca.geno)
clust.pca.geno.1 <- clust.pca.geno %>%
                    dplyr::select(c(genotipo, clust))



write_delim (clust.pca.geno.1, "./Data/procdata/PCA.GENO.txt", delim = ",", 
             na = "NA", append = FALSE, quote_escape = "double")

#View (clust.pca.pheno)

fviz_pca_ind (res.PCA.GENO,axes = c(1, 2),
              geom.ind = "point", # show points only (nbut not "text")
              pointsize = 2, #Size and shape of plot elements
              col.ind = clust.pca.geno$clust, # color by groupspointshape = 21
              pointshape = 16,
              palette = c("orange", "royalblue", "firebrick1"),
              #palette = "uchicago",
              addEllipses = TRUE,# Concentration ellipses
              legend.title = "Groups")




#Mixed model: Eigenanalysis (PCA as random component)
(pcaR.GWAS <- gwas.analysis (crossobj= dt.EEMAC07.1,
                             method="eigenstrat",
                             provide.K=FALSE, 
                             covariates=pca.EX$scores, 
                             trait="EX",
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

