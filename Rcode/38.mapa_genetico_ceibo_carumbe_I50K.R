###############################################################
#  codigo para la generacion del mapa de genetico            ##
#  de una poblacion RIL F6 de Ceibo x carumbe                ##
#  Gaston Quero - Nicolas Mastandrea                         ##
# 20/02/2020                                                 ##
###############################################################

getwd ()

setwd ("R:/Mapa_genetico_cebada")

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

#library(lmem.qtler)

## Resumen
# Construcción mapa de cebada 
# Se usaron los datos geneticos de Cebada x Carumbe
# datos enviados por Ariel Castro                                   
# fecha:	28 de junio de 2018, 15:26
# Se incorpora el mapa de Nicolas Mastandrea 
# enviados el 20/02/2020

### Pre_mapa
# Identificacion de ID
Ids_ceibo_carumbe <- read_delim (file="./Data/rawdata/ids_rils_ceibo_carumbe.txt", 
                     delim = "\t", na = "NA")


Ids_mkrs_eurekaSSR <- read_delim (file="./Data/rawdata/ids_markers_eurekaSSR.txt", 
                                 delim = "\t", na = "NA")


Ids_mkrs_pia <- read_delim (file="./Data/rawdata/ids_markers_pia.txt", 
                                  delim = "\t", na = "NA")

##### 
I50K <- read_delim (file="./Data/rawdata/Ill50K.txt", 
                      delim = "\t", na = "NA")

### generar el geno ####
genos.I50K <- I50K  %>%
              dplyr::select (-c (Index, Name, Chr, Position,...7,...16,...29,...37,...42,...76,...80))

list.name <- (colnames(genos.I50K))

geno.I50k.1 <- bind_cols (lapply(list.name , function (filt.name ){
  
  #filt.name = "Mo1-001.GType"
  id.x <- genos.I50K %>%
          dplyr::select (all_of(filt.name ))
  
  id.y <- Ids_ceibo_carumbe %>%
          dplyr::filter (cod_1 == filt.name )
  
  colnames (id.x) <-  id.y$IND
  
  id.x
}))

mrk.I50K <- I50K %>%
            dplyr::select (Name)

geno.I50k.1 <- bind_cols (mrk.I50K ,geno.I50k.1)

write_delim (geno.I50k.1, path ="./Data/procdata/geno_I50K.1.csv",
             delim = ",", na = "-")

## ahora genero pmap
pmap_I50K <- I50K  %>%
             dplyr::select (Name,Chr,Position)

write_delim (pmap_I50K, file ="./Data/procdata/pmap_I50K.1.csv",
             delim = ",", na = "-")


### generar la pheno
colnames (geno.I50k.1)
id <- colnames (geno.I50k.1) [2:81] 
random <- rnorm(80, 0, 1)
length(random)
length(id)
pheno_I50K.1 <- rbind(id, random)
pheno_I50K.1 <- as.data.frame (pheno_I50K.1)
colnames (pheno_I50K.1) <- id
pheno_I50K.1  <- pheno_I50K.1 [-1,]

pheno_I50K.1 <- pheno_I50K.1 %>%
                dplyr::mutate (X="trait")%>%
                dplyr::select (X,everything())

write_delim (pheno_I50K.1, path ="./Data/procdata/pheno_I50K.1.csv",
             delim = ",", na = "-")

##### 
## Datos
# Se cargan los datos originales en el formato de *rqtl*
# cross I50K.2
cross_I50K.2 <- read.cross (format="tidy",
                                dir="./Data/procdata",
                                genfile ="geno_I50K.1.csv", 
                                mapfile ="pmap_I50K.1.csv" , 
                                phefile ="pheno_I50K.1.csv" , 
                                na.strings="-",
                                alleles=c("A","B"),
                                estimate.map=FALSE, 
                                convertXdata=TRUE, error.prob=0.0001,
                                map.function="kosambi",
                                #F.gen=6, 
                                crosstype ="riself")



## ## Análisis
### Curado de Datos

plotMissing (cross_I50K.2)

#Se eliminan los marcadores con mas del 50 % de datos faltantes
n.missing <- nmissing (cross_I50K.2, what="mar")[(nmissing(cross_I50K.2, 
                                                               what="mar"))/sum(summary(cross_I50K.2)$n.ind) > 0.51]

#Por datos faltantes se eliminados *
length((names.marker <- c (names(n.missing))))

#Los marcadores eliminados fueron
(names.marker <- c (names(n.missing)))

# aca se eliminan
cross_I50K.3  <- drop.markers (cross_I50K.2, names.marker)

#Se eliminan los individuos con mas de 70 % no genotipado
cross_I50K.4  <- subset (cross_I50K.3  , 
                               ind=(ntyped(cross_I50K.3 ) 
                                    >((totmar(cross_I50K.3 ) * 70)/100)))

indiv <- subset (cross_I50K.3 , 
                 ind=(ntyped(cross_I50K.3 )  
                      < ((totmar(cross_I50K.3 ) * 51)/100)))

#Los individuos eliminados fueron
#(indiv$pheno$id)

#Esta es la matriz despues de filtrar por datos faltantes
plotMissing (cross_I50K.4 )

#### Individuos duplicados

#Resulta útil comparar los genotipos entre todos
#los pares de individuos, 
#para revelar pares con genotipos inusualmente similares,
#que pueden indicar duplicaciones de muestra o gemelos monocigóticos. 
#En cualquier caso, querremos eliminar a un individuo de cada par. 

dat.cg <- comparegeno (cross_I50K.4)
hist(dat.cg[lower.tri(dat.cg)], breaks=seq(0, 1, len=101), xlab="matching genotypes", main= NULL)
rug(dat.cg[lower.tri(dat.cg)])

#Por ahora no voy a sacar ningun individuo.

#### Distorsion de segregación
gt <- geno.table (cross_I50K.4)
gt [gt$P.value < 0.05/totmar(cross_I50K.4),]
total <- nind (cross_I50K.4)

gt <- geno.table (cross_I50K.4)
missing <- gt[,2]
num.obs <- c(total - missing)
obs.alelo <- gt[, c(2,3,4)]

Geno.freq <- as.vector(obs.alelo/num.obs)

gfrq <- Geno.freq %>%
  mutate (mkrs = rownames(Geno.freq)) %>%
  select(mkrs, everything()) 
imkrs <- unique (rownames(Geno.freq))

freq.geno <- lapply (imkrs, function (filtro) { 
  GF <- gfrq %>% 
    dplyr::filter (mkrs == filtro) %>%
    tidyr::gather ("AA","BB", key="alelle",value ="frq")
  return (GF)
})

df.geno.freq <- as_tibble(do.call (rbind, freq.geno))

ggbarplot (df.geno.freq, "mkrs", "frq",
           fill = "alelle", color = "alelle", 
           palette =c("orangered1","royalblue4"),
           label = FALSE, lab.col = "white", lab.pos = "out") +
  geom_hline(yintercept = 0.5 ,  colour = "black", linetype = "dashed") +
  rremove("x.text")

#Se sacaron los siguientes marcadores 

(todrop <- rownames (gt[gt$P.value < 1e-10,]))

# Aca se sacan
cross_I50K.5 <- drop.markers (cross_I50K.4, todrop)

#### Distorsion de segregación
gt <- geno.table (cross_I50K.5)
gt [gt$P.value < 0.05/totmar(cross_I50K.5),]
total <- nind (cross_I50K.5)

gt <- geno.table (cross_I50K.5)
missing <- gt[,2]
num.obs <- c(total - missing)
obs.alelo <- gt[, c(2,3,4)]

Geno.freq <- as.vector(obs.alelo/num.obs)

gfrq <- Geno.freq %>%
  mutate (mkrs = rownames(Geno.freq)) %>%
  select(mkrs, everything()) 
imkrs <- unique (rownames(Geno.freq))

freq.geno <- lapply (imkrs, function (filtro) { 
  GF <- gfrq %>% 
    dplyr::filter (mkrs == filtro) %>%
    tidyr::gather ("AA","BB", key="alelle",value ="frq")
  return (GF)
})

df.geno.freq <- as_tibble(do.call (rbind, freq.geno))

ggbarplot (df.geno.freq, "mkrs", "frq",
           fill = "alelle", color = "alelle", 
           palette =c("orangered1","royalblue4"),
           label = FALSE, lab.col = "white", lab.pos = "out") +
  geom_hline(yintercept = 0.5 ,  colour = "black", linetype = "dashed") +
  rremove("x.text")

#Ahora quedan 

plotMap (cross_I50K.5)


# Realizo un analisis de calidad de la nueva matriz uso el *ASMap*


#################
cross_I50K.6 <- mstmap (cross_I50K.5, 
                              id = "id",bychr = TRUE, 
                              anchor = TRUE,
                              dist.fun = "kosambi", trace = FALSE)


#sg.d <- statGen (cross_I50K.6,id = "id",  bychr = TRUE, 
 #                stat.type = c("xo","dxo","miss"))

#Estos individuos son clones
gc.6 <- genClones (cross_I50K.6,id = "id", tol = 0.95)
gc.6$cgd

#Se fijan los clones

cg5737 <- gc.6$cgd

cross_I50K.7 <- fixClones (cross_I50K.6, 
                           cg5737, id = "id", 
                                 consensus = TRUE)

#Ahora quedan 
(nind (cross_I50K.7))

#Sin un mapa construido es imposible determinar el origen de la segregación.
#distorsión. Sin embargo, el uso ciego de marcadores distorsionados también puede
#generar problemas en la creacion del mapa.
#Puede ser más sensato dejar a un lado los marcadores distorsionados y
#construir el mapa con marcadores menos problemáticos. 
#Una vez que se construye el mapa se pueden introducir marcadores más problemáticos
#para determinar si tienen una utilidad o efecto nocivo en el mapa.

mm.e.AA <- statMark (cross_I50K.7, stat.type = "marker")$marker$AA

cross_I50K.8 <- drop.markers(cross_I50K.7, c(markernames(cross_I50K.7)[mm.e.AA > 0.98],
                                                         markernames(cross_I50K.7)[mm.e.AA < 0.2]))

cross_I50K.8 <- pullCross(cross_I50K.8, type = "missing", pars = list(miss.thresh = 0.1))
cross_I50K.8 <- pullCross(cross_I50K.8, type = "seg.distortion", pars = list(seg.thresh = "bonf"))

sum(ncol(cross_I50K.8$missing$data), ncol(cross_I50K.8$seg.dist$data), ncol(cross_I50K.8$co.located$data))

summary (cross_I50K.8  )

plotMissing(cross_I50K.8)

cross_I50K.9 <- mstmap (cross_I50K.8,id = "id",
                        bychr = FALSE, trace = TRUE, dist.fun = "kosambi", p.value = 1e-06)

chrlen (cross_I50K.9)
heatMap(cross_I50K.9, lmax = 50)

plotMap (cross_I50K.9)

class (cross_I50K.9 )

write.cross (cross_I50K.9, format= "tidy",
            filestem="./Data/procdata/cross_I50K.9", digits=4)

##### ACA termina  la generacion del mapa de ilumina #####################