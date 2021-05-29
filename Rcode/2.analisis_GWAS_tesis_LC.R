#########################################################################################
#  Codigo para el analisis de los datos de la tesis de Lorena Camarota                 ##
#                                                                                      ##
# Los datos fueron enviados                                                            ##
# de:	Lorena Cammarota <lorena.cammarota@gmail.com>                                    ##
# para:	Gastón Quero <gastonquero@gmail.com>                                           ##  
# fecha:	24 oct. 2020 9:29                                                            ##
# asunto:	archivos                                                                     ##
#                                                                                      ##
# Gaston Quero                                                                         ##
# 26/10/2020                                                                           ##
#########################################################################################

getwd ()

setwd ("C:/Users/Usuario/Dropbox/Tesis_Lorena")


# cargar paquetes
library("qtl")
qtlversion()
library("FactoMineR")
library("dplyr")
library("car")
library (lme4)
library (emmeans)
library ("nlmrt")
library ("easynls")
library (tidyverse)
library ("lattice")
library ("latticeExtra")
library (multcompView)
library("ggridges")
library("viridis")
library("lmerTest")
library("lubridate")
library (ggcorrplot)
library (sjPlot)
library (sjlabelled)
library (sjmisc)
library (ggpubr)
library("ggsci")
library("factoextra")
library("corrplot")


#### cargo los datos de la matriz de genotipos 

cebada_geno <- read_delim ("./Data/rawdata/mapa_training_population.txt" , 
                         delim = "\t", na = "NA")


cebada_geno <- cebada_geno %>%
               dplyr::mutate(Chr= fct_recode  (Chr, "1" = "1H"))%>%
               dplyr::mutate(Chr= fct_recode  (Chr, "2" = "2H"))%>%
               dplyr::mutate(Chr= fct_recode  (Chr, "3" = "3H"))%>%
               dplyr::mutate(Chr= fct_recode  (Chr, "4" = "4H"))%>%
               dplyr::mutate(Chr= fct_recode  (Chr, "5" = "5H"))%>%
               dplyr::mutate(Chr= fct_recode  (Chr, "6" = "6H"))%>%
               dplyr::mutate(Chr= fct_recode  (Chr, "7" = "7H"))


# saco las columnas que no me sirve y renombro los marcadores
str(cebada_geno )

##
cebada_geno.1 <- cebada_geno %>%
                 dplyr::select (-c(Index, Address)) %>%
                 dplyr::mutate (id.1 = str_c("S", Chr))%>%
                 dplyr::mutate (marker = str_c(id.1, pos, sep="_"))%>%
                 dplyr::select (-id.1)%>%
                 dplyr::select (marker, everything())

cebada_geno.2 <- cebada_geno.1 %>%
                 dplyr::select (-rs)
#dt=cebada_geno.2

run_format_gwas.cross <- function (dt= NULL){
  

list.mrk <- unique (dt$marker)
  
  dt.1 <- bind_cols(lapply (list.mrk , function (filt.mrk) {
    #filt.mrk = "M_2_130417669"
    
    print (filt.mrk)
    mrk.x <- dt %>%
             dplyr::filter (marker == filt.mrk) 
    
    mrk.x1 <- mrk.x %>%
              tidyr::pivot_longer(-c(marker, Chr,pos), names_to = "genotypes", values_to = "X1")
    
    #mrk.x1$X1 <- as.factor(mrk.x1$X1)
  # unique(mrk.x1$X1)
   
   mrk.x1a <- mrk.x1 %>%
             dplyr::select ( -c(marker,Chr,pos))
   
    mrk.x2 <- mrk.x1a %>%
               dplyr::mutate( X1 = str_replace (X1, "AA", "1"))%>%
               dplyr::mutate( X1 = str_replace (X1, "BB", "2"))%>%
               dplyr::mutate( X1 = str_replace (X1, "[^12]", "-"))
 
    names (mrk.x2) <- c("genotypes", filt.mrk)
  
    mrk.x3 <- mrk.x2 %>%
              dplyr::select(-genotypes)
 
  }))
  
genos.id <- dt [1,]
  
genos.idx <- dt [1,] %>%
            tidyr::pivot_longer(-c(marker, Chr,pos), names_to = "genotypes", values_to = "X1")%>%
            dplyr::select (genotypes)
  
dt.2 <- bind_cols ( genos.idx,dt.1) %>%
        dplyr::select (genotypes, everything())

write_delim(dt.2, path="./Data/procdata/cebada_geno_1.txt", delim = "\t", na = "-", 
            append = FALSE, col_names = TRUE, quote_escape = "double")

map.1 <- dt %>%
         dplyr::select (marker, Chr,pos)


write_delim(map.1, path="./Data/procdata/cebada_map_1.txt", delim = "\t", na = "-", 
            append = FALSE, col_names = FALSE, quote_escape = "double")

return(dt.2)
}# aca termmina run_format_gwas.cross 
  
  
cebada_geno.3 <- run_format_gwas.cross (dt= cebada_geno.2 )


#### 
EEMAC_pheno <- read_delim (file = "./Data/rawdata/EEMAC_pheno.txt" ,
                            col_names = TRUE, delim = "\t", na = "NA")



EEMAC_pheno <- EEMAC_pheno%>%
               dplyr::arrange( genotypes )

head( EEMAC_pheno )
## Geno data
G.EEMAC.data <- read.table ("./Data/procdata/cebada_geno_1.txt",
                            header = T, sep = "\t",
                            dec = ".", na.strings = "-")

setdiff (EEMAC_pheno$genotypes,G.EEMAC.data$genotypes )

setdiff (G.EEMAC.data$genotypes,EEMAC_pheno$genotypes)

G.EEMAC.data.1 <- G.EEMAC.data %>%
                  dplyr::filter (genotypes !="Inno2_25")%>%
                  dplyr::filter (genotypes !="CLE268")%>%
                  dplyr::arrange( genotypes )

setdiff (EEMAC_pheno$genotypes,G.EEMAC.data.1$genotypes)

setdiff (G.EEMAC.data.1$genotypes,EEMAC_pheno$genotypes)


G.EEMAC.data.1$genotypes == EEMAC_pheno$genotypes

G.EEMAC.data.1 <- G.EEMAC.data.1 %>%
                  dplyr::mutate (genotype = paste ("G",G.EEMAC.data.1$genotypes,sep="_"))%>%
                  dplyr::select (-genotypes) %>%
                  dplyr::select (genotype, everything())%>%
                  dplyr::arrange (genotype)

EEMAC_pheno <- EEMAC_pheno%>%
               dplyr::mutate (genotype = paste ("G",EEMAC_pheno$genotypes,sep="_"))%>%
               dplyr::select (-genotypes)%>%
               dplyr::select (genotype, everything())%>%
               dplyr::arrange (genotype)

setdiff (EEMAC_pheno$genotype,G.EEMAC.data.1$genotype)

setdiff (G.EEMAC.data.1$genotype,EEMAC_pheno$genotype)

G.EEMAC.data.1$genotype == EEMAC_pheno$genotype


map.cebada.data <- read.table ("./Data/procdata/cebada_map_1.txt",
                             header = F, sep = "\t",
                             dec = ".", na.strings = "NA")


setdiff (EEMAC_pheno$genotype,G.EEMAC.data.1$genotype)

setdiff (G.EEMAC.data.1$genotype,EEMAC_pheno$genotype)




P.data   <- EEMAC_pheno
G.data   <- G.EEMAC.data.1
map.data <- map.cebada.data

dt = G.data 
      
list.geno <- unique (G.data$genotype)

dt.genotype.missing <- bind_rows (lapply (list.geno  , function (filt.geno) {
 # filt.geno = "G_Inno1_111"
  print (filt.geno)
  geno.x <- dt %>%
           dplyr::filter (genotype ==  filt.geno) 
  
  geno.x1 <- geno.x %>%
             tidyr::pivot_longer(-genotype, names_to = "marcador", values_to = "alelo")
  
  geno.x1$alelo <- as.character(geno.x1$alelo)
  
              
    Nax <- geno.x1 %>%
           group_by (genotype)%>%
           count (alelo) %>%
           dplyr::mutate (total.mrk = sum (n))%>%
           dplyr::mutate (porc = (n * 100)/total.mrk)%>%
           dplyr::ungroup()
    
    
    if (nrow (Nax) != 3){
      xx <- unique (Nax$total.mrk)
      gg <- tibble (genotype = filt.geno, alelo="x", n=0, total.mrk =xx , porc=0)
      
      Nax <- bind_rows (Nax, gg )
      return (Nax)
      }
    return (Nax)
  }))



N1 <- dt.genotype.missing %>%
      dplyr::filter ( alelo  == 1)

N2 <- dt.genotype.missing %>%
  dplyr::filter ( alelo  == 2)


Nx <- dt.genotype.missing %>%
  dplyr::filter ( alelo  == "x")

Nas <- dt.genotype.missing %>%
  dplyr::filter ( is.na (alelo))%>%
  dplyr::arrange (desc(n))%>%
  dplyr::filter (porc > 40)

id <- Nas$genotype

G.data <- G.data%>%
          dplyr::filter (genotype != id [1])%>%
          dplyr::filter (genotype != id [2])

P.data <- P.data%>%
          dplyr::filter (genotype != id [1])%>%
          dplyr::filter (genotype != id [2])
  
##### loop para datos faltantes
na_count <- apply (G.data, 1, function(x) sum(is.na(x)))
porcentaje <- 51 # menos del % de missing data

hp1 <- (porcentaje * (ncol(G.data)-1))/ 100
hp2 <- NULL
hp3 <- NULL # son los que quedan 
length(hp3)
#i=1
for(i in 1:length(na_count)){
  if (na_count [i] < hp1 ){
    hp2 <- cbind(hp2, na_count[i])
    hp3 <- cbind (hp3, i)
  }
}



G.data  <- G.data[hp3, ]
dim (G.data)
str(G.data)

G.data [is.na(G.data)] <- "-"
G.data[G.data == 1] <- "1"
G.data[G.data == 2] <- "2"

str(map.data)

#Cargar los datos lmem.gwaser 
cross.data.cebada.EEMAC <- gwas.cross (P.data, G.data , map.data,
                                 cross='gwas', heterozygotes=FALSE)


total <- nind (cross.data.cebada.EEMAC)
gt <- geno.table (cross.data.cebada.EEMAC)
missing <- gt[,2]
num.obs <- c(total - missing)
obs.alelo <- gt[, c(3,4)]
Geno.freq <- as.vector(obs.alelo/num.obs)

SNPs.MAF <-  Geno.freq[,2] < 0.1 
rownames(Geno.freq [SNPs.MAF,])
names.marker <- c (rownames(Geno.freq [SNPs.MAF,]))
length(names.marker)

cross.data.cebada.EEMAC.1 <- drop.markers (cross.data.cebada.EEMAC, names.marker)

summary (cross.data.cebada.EEMAC.1)

plotMissing (cross.data.cebada.EEMAC.1 )

geno.image (cross.data.cebada.EEMAC.1  , reorder=FALSE, main="cross.data.cebada.EEMAC.1 ",
            alternate.chrid=FALSE)


pca.EEMAC.1 <- pca.analysis (crossobj=cross.data.cebada.EEMAC.1, p.val=0.05)







G <- G.data %>%
  inner_join (P.data, by="genotype") %>%
  dplyr::select (c(genotype,sitio,EX,NS,BG,PT, everything()))


G.PCA.1 <- G  %>%
  dplyr::select (-c(sitio))

write_delim(G.PCA.1, path= "./Data/procdata/G.PCA.1.txt", 
            delim = ",", na = "-")



G.PCA.2 <- read_delim (file = "./Data/procdata/G.PCA.1.txt" ,
                       col_names = TRUE, delim = ",", na = "-")

G.PCA.2 [1:5, 1:5]

res.pca.1 <- PCA (G.PCA.2, scale.unit = TRUE, ncp = 10, 
                  ind.sup = NULL, 
                  quanti.sup = c(2,3,4,5), quali.sup = 1, row.w = NULL, 
                  col.w = NULL, graph = TRUE, axes = c(1,2))
print (res.pca.1)


fviz_eig (res.pca.1, addlabels = TRUE)


##### individuos
ind <- get_pca_ind(res.pca.1)
ind
fviz_pca_ind(res.pca.1)
res.pca.1$quali

#View (clust.pca.pheno)

fviz_pca_ind (res.pca.1,
              geom.ind = "point", # show points only (nbut not "text")
              #col.ind = G.PCA.2$ciclo.x, # color by groups
              #palette = c("orange", "navyblue", "red"),
              #addEllipses = TRUE, # Concentration ellipses
              #legend.title = "Groups"
              )

fviz_pca_ind (res.pca.1,
              geom.ind = "point", # show points only (nbut not "text")
              col.ind = G.PCA.2$grupo, # color by groups
              palette = c("darkgreen", "navyblue", "red"),
              #addEllipses = TRUE, # Concentration ellipses
              legend.title = "Groups")


#fake.1 <- subset(cross.data.cebada.EEMAC, chr=c("1", "2", "3","4", "5", "6", "7"))

#cross.chr.2 <- subset(cross.data.cebada.EEMAC, chr="2")
#cross.chr.1 <- subset(cross.data.cebada.EEMAC, chr="1")
#cross.chr.3 <- subset(cross.data.cebada.EEMAC, chr=c("1","3"))

#geno.image (cross.chr.2  , reorder=FALSE, main="cross.data.cebada.EEMAC.1 ",
            #alternate.chrid=FALSE)



fake.1 <- subset(cross.data.cebada.EEMAC.1, chr=c("1"))
fake.2 <- subset(cross.data.cebada.EEMAC.1, chr="2")
fake.3 <- subset(cross.data.cebada.EEMAC.1, chr="3")
fake.4 <- subset(cross.data.cebada.EEMAC.1, chr="4")


pca.EEMAC.1 <- pca.analysis (crossobj=fake.1 , p.val=0.05)

pca.EEMAC.2 <- pca.analysis (crossobj=fake.2 , p.val=0.05)
pca.EEMAC.3 <- pca.analysis (crossobj=fake.3 , p.val=0.05)
pca.EEMAC.4 <- pca.analysis (crossobj=fake.4 , p.val=0.05)

ind$coord
### IS.rel.z2
# Mixed model: Eigenanalysis (PCA as random component)
(pca.cebada.EEMAC.EX <- gwas.analysis (crossobj=cross.data.cebada.EEMAC, method="eigenstrat",
                                      provide.K=FALSE, covariates=ind$coord,
                                      trait="EX", threshold="Li&Ji",
                                      p=0.05,  out.file="GWAS PCA as Random model"))$selected








#PCA

# Chr 1 
write.cross (cross.data.z14_15, format="csv",
             filestem="./Data/procdata/cross.data.z14_15.chr1", chr=1, digits=NULL)




pca.EEMAC <- pca.analysis (crossobj=cross.data.cebada.EEMAC, p.val=0.05)


