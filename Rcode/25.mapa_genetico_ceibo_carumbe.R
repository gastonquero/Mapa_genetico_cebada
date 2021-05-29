###############################################################
#  codigo para la generacion del mapa de genetico            ##
#  de una poblacion RIL F6 de Ceibo x carumbe                ##
#  Gaston Quero - Nicolas Mastandrea                         ##
# 20/02/2020                                                 ##
###############################################################

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


geno_eurekaSSR <- read_delim (file="./Data/rawdata/geno_eurekaSSR.txt", 
                            delim = "\t", na = "NA")


III50K <- read_delim (file="./Data/rawdata/Ill50K.txt", 
                              delim = "\t", na = "NA")

head (III50K)
###### limpiar la matriz III50K

III50K.1 <- III50K %>%
            dplyr::select (c(Name, Chr,  Position),starts_with("Mo1"))


## cambiar los IDs de la matriz III50K

id_genotypes_III50K.1 <- colnames (III50K.1) [4:83]
class (id_genotypes_III50K.1)
length(id_genotypes_III50K.1)

Ids_ceibo_carumbe.1 <- Ids_ceibo_carumbe %>%
                       dplyr::filter (cod_1 %in%id_genotypes_III50K.1)

id_genotypes_III50K.1 == Ids_ceibo_carumbe.1$cod_1
colnames (III50K.1) [4:83] == Ids_ceibo_carumbe.1$cod_1

colnames (III50K.1) [4:83] <- Ids_ceibo_carumbe.1$IND


# verificar la cantidad de genotipos entre III50K y eureka 

colnames (geno_eurekaSSR)[2:91]

geno_eurekaSSR.1 <- geno_eurekaSSR %>%
                    dplyr::select(c(cod_1, colnames (III50K.1)[4:83]))
                    
colnames (geno_eurekaSSR.1) [2:81] == colnames (III50K.1) [4:83]

dim(geno_eurekaSSR.1)



### agregar a III50K los marcadores de Eureka

geno_eurekaSSR.2 <- geno_eurekaSSR.1 %>%
                    inner_join (Ids_mkrs_eurekaSSR, by="cod_1")
 
head(geno_eurekaSSR.2)
geno_eurekaSSR.3 <- geno_eurekaSSR.2 %>%
                    dplyr::select (-c(cod_1, marker_1))%>%
                    dplyr::select (marker_2, everything())%>%
                    dplyr::rename (Name =marker_2) %>%
                    dplyr::mutate (Chr="XH")%>%
                    dplyr::mutate  (Position = 1000) %>%
                    dplyr::select  (Name, Chr,Position,everything() )
tail (geno_eurekaSSR.3)
dim(geno_eurekaSSR.3) 
dim(III50K.1)
1684 + 148

III50K.2 <- rbind (III50K.1,geno_eurekaSSR.3)
dim(III50K.2)
 
## ahora genero pmap

pmap_III50K.2 <- III50K.2 %>%
                 dplyr::select (Name,Chr,Position)
                

write_delim (pmap_III50K.2, path ="./Data/procdata/pmap_III50K.2.csv",
             delim = ",", na = "-")

### generar el geno

geno_III50K.2 <- III50K.2 %>%
                 dplyr::select (-c(Chr,Position))


write_delim (geno_III50K.2, path ="./Data/procdata/geno_III50K.2.csv",
             delim = ",", na = "-")

### generar la pheno
colnames (III50K.2)
id <- colnames (III50K.2) [4:83]
random <- rnorm(80, 0, 1)
length(random)
length(id)
pheno_III50K.2 <- rbind(id, random)
pheno_III50K.2 <- as.data.frame (pheno_III50K.2)
colnames (pheno_III50K.2) <- id
pheno_III50K.2  <- pheno_III50K.2 [-1,]

pheno_III50K.2 <- pheno_III50K.2 %>%
                  dplyr::mutate (X="trait")%>%
                  dplyr::select (X,everything())


write_delim (pheno_III50K.2, path ="./Data/procdata/pheno_III50K.2.csv",
             delim = ",", na = "-")

## Datos
# Se cargan los datos originales en el formato de *rqtl*
# cross Ill50K.2
Ill50K_Mapeo.qtl <- read.cross (format="tidy",
                                dir="./Data/procdata",
                                genfile ="geno_III50K.2.csv", 
                                mapfile ="pmap_III50K.2.csv" , 
                                phefile ="pheno_III50K.2.csv" , 
                                na.strings="-",
                                alleles=c("A","B"),
                                estimate.map=FALSE, 
                                convertXdata=TRUE, error.prob=0.0001,
                                map.function="kosambi",
                                #F.gen=6, 
                                crosstype ="riself")

Ill50K_Mapeo.qtl <- jittermap(Ill50K_Mapeo.qtl)
plotMap(Ill50K_Mapeo.qtl)
summary (Ill50K_Mapeo.qtl)

## ## Análisis
### Curado de Datos

#### Datos faltantes 
# Datos faltantes de la matriz original

plotMissing (Ill50K_Mapeo.qtl)

#Se eliminan los marcadores con mas del 50 % de datos faltantes
n.missing <- nmissing (Ill50K_Mapeo.qtl, what="mar")[(nmissing(Ill50K_Mapeo.qtl, 
             what="mar"))/sum(summary(Ill50K_Mapeo.qtl)$n.ind) > 0.51]

#Por datos faltantes se eliminados *
length((names.marker <- c (names(n.missing))))

#Los marcadores eliminados fueron
(names.marker <- c (names(n.missing)))

# aca se eliminan
Ill50K_Mapeo.qtl.a <- drop.markers (Ill50K_Mapeo.qtl, names.marker)

#Se eliminan los individuos con mas de 70% no genotipado
Ill50K_Mapeo.qtl.b  <- subset (Ill50K_Mapeo.qtl.a , 
                               ind=(ntyped(Ill50K_Mapeo.qtl.a) 
                                    >((totmar(Ill50K_Mapeo.qtl.a) * 70)/100)))

indiv <- subset (Ill50K_Mapeo.qtl.a, 
                 ind=(ntyped(Ill50K_Mapeo.qtl.a)  
                      < ((totmar(Ill50K_Mapeo.qtl.a) * 70)/100)))

#Los individuos eliminados fueron
(indiv$pheno$id)

#Esta es la matriz despues de filtrar por datos faltantes
plotMissing (Ill50K_Mapeo.qtl.b)

#### Individuos duplicados

#Resulta útil comparar los genotipos entre todos
#los pares de individuos, 
#para revelar pares con genotipos inusualmente similares,
#que pueden indicar duplicaciones de muestra o gemelos monocigóticos. 
#En cualquier caso, querremos eliminar a un individuo de cada par. 

dat.cg <- comparegeno (Ill50K_Mapeo.qtl.b)
hist(dat.cg[lower.tri(dat.cg)], breaks=seq(0, 1, len=101), xlab="matching genotypes", main= NULL)
rug(dat.cg[lower.tri(dat.cg)])

#Por ahora no voy a sacar ningun individuo.

#### Distorsion de segregación

gt <- geno.table (Ill50K_Mapeo.qtl.b)
gt [gt$P.value < 0.05/totmar(Ill50K_Mapeo.qtl.b),]
total <- nind (Ill50K_Mapeo.qtl.b)

gt <- geno.table(Ill50K_Mapeo.qtl.b)
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
Ill50K_Mapeo.qtl.c <- drop.markers (Ill50K_Mapeo.qtl.b, todrop)

#Ahora quedan 
totmar(Ill50K_Mapeo.qtl.c)

# Realizo un analisis de calidad de la nueva matriz uso el *ASMap*


#################
Ill50K_Mapeo.qtl.d <- mstmap (Ill50K_Mapeo.qtl.c, 
                              id = "id",bychr = TRUE, 
                              anchor = TRUE,
                              dist.fun = "kosambi", trace = FALSE)


sg.d <- statGen (Ill50K_Mapeo.qtl.d,id = "id",  bychr = TRUE, 
               stat.type = c("xo","dxo","miss"))

#Estos individuos son clones
gc.d <- genClones (Ill50K_Mapeo.qtl.d,id = "id", tol = 0.95)
gc.d$cgd

#Se fijan los clones

cgd <- gc.d$cgd

Ill50K_Mapeo.qtl.e <- fixClones (Ill50K_Mapeo.qtl.d, 
                                 cgd, id = "id", 
                                 consensus = TRUE)

#Ahora quedan 
(nind (Ill50K_Mapeo.qtl.e))

#Sin un mapa construido es imposible determinar el origen de la segregación.
#distorsión. Sin embargo, el uso ciego de marcadores distorsionados también puede
#generar problemas en la creacion del mapa.
#Puede ser más sensato dejar a un lado los marcadores distorsionados y
#construir el mapa con marcadores menos problemáticos. 
#Una vez que se construye el mapa se pueden introducir marcadores más problemáticos
#para determinar si tienen una utilidad o efecto nocivo en el mapa.

mm.e.AA <- statMark(Ill50K_Mapeo.qtl.e, stat.type = "marker")$marker$AA

Ill50K_Mapeo.qtl.f <- drop.markers(Ill50K_Mapeo.qtl.e, c(markernames(Ill50K_Mapeo.qtl.e)[mm.e.AA > 0.98],
                                                         markernames(Ill50K_Mapeo.qtl.e)[mm.e.AA < 0.2]))

Ill50K_Mapeo.qtl.g <- pullCross (Ill50K_Mapeo.qtl.f, 
                                 type = "seg.distortion", 
                                 pars = list(seg.thresh = "bonf"))


#Ahora quedan 
(totmar(Ill50K_Mapeo.qtl.g))

#Esta es la poblacion con la cual hago el primer mapa.

### Grupos de ligamiento
# Se construyen los grupos de ligamiento con  *ASMap*
  
Ill50K_Mapeo.qtl.h <- mstmap (Ill50K_Mapeo.qtl.g, 
                              id="id",
                              bychr = FALSE, 
                              trace = TRUE, 
                              dist.fun ="kosambi")


plotMap(Ill50K_Mapeo.qtl.h)
#heatMap (Ill50K_Mapeo.qtl.h, markDiagonal = TRUE, 
 #        lmax = 50, main=NULL )


### Asignación de los grupos de ligamiento

#Se cargan los datos de Bill Thomas de un mapa de referencia de 
#**Golden Promise x Morex **
 # Los datos fueron enviados por : 
#  de:	Maximiliano Verocai <maxiverocaipay@gmail.com>
#  para:	gastonquero@gmail.com
#fecha:	18 abr. 2019 10:21
#asunto:	Fwd: Posiciones de marcadores

GPxM <- 
  read.table (file="./Data/rawdata/GPxM_Bill_Thomas.csv",
              header = TRUE, 
              sep = ",",
              quote = "\"'", 
              dec = ".",
              na.strings = "-")
str (GPxM)

# Reemplazar los caracteres de los ID
#mr <- GPxM$snp
#mr  <- data.frame(sapply(mr , function(x) {
  #gsub("-", "_", x)
#}))
#xx <- mr [,1]
#GPxM <- GPxM %>%
 #       dplyr::mutate (snp = xx)

GPxM$snp <- as.character (GPxM$snp)
dim (GPxM)
head(GPxM)

#Este es el mapa que se hizo con la población  
# **Ill50K_Mapeo.qtl.h**
  
GPxM$snp
map.1 <- pull.map (Ill50K_Mapeo.qtl.h, as.table=TRUE)

map.1.dt <- map.1 %>%
            dplyr::mutate(mrk_names = rownames(map.1)) %>%
            dplyr::select ("mrk_names", "chr" ,"pos")

unique(map.1.dt$chr)
# Se filtran en el mapa de Bill Thomas los marcadores que coinciden con los nuestros.
 # Nosotros tenemos 
(length (unique (map.1.dt$mrk_names)))
#y 
# Bill Thomas 
(length (unique (GPxM$snp)))


mrk.map.1 <- unique (map.1.dt$mrk_names)

str(GPxM)
str (map.1.dt)


GPxM_1 <- GPxM %>%
          dplyr::filter (snp %in% mrk.map.1) %>%
          dplyr::arrange (genetic_map_chr)


#Coinciden 
(length (GPxM_1$snp))

#Esos marcadores los vamos a usar para como *tags* 
# de los grupos de ligamientos generados

#Cromosoma **1H**

GPxM_1_1H <-  GPxM_1 %>%
              dplyr::filter (genetic_map_chr =="1H")

mrk.GPxM_1.1H <- unique (GPxM_1_1H$snp)

(map.1.dt_1H <- map.1.dt %>%
                dplyr::filter  (mrk_names %in% mrk.GPxM_1.1H) %>%
                dplyr::arrange (chr))


#El **cromosoma 1H** se asocia con los grupos de ligamiento 
 (unique (map.1.dt_1H$chr))


#Cromosoma  **2H**

GPxM_1_2H <-  GPxM_1 %>%
               dplyr::filter (genetic_map_chr =="2H")

mrk.GPxM_1.2H <- unique (GPxM_1_2H$snp)

map.1.dt_2H <- map.1.dt %>%
               dplyr::filter  (mrk_names %in% mrk.GPxM_1.2H) %>%
               dplyr::arrange (chr)

# El **cromosoma 2H** se asocia con los grupos de ligamiento 
(unique (map.1.dt_2H$chr))

# Cromosoma  **3H**

GPxM_1_3H <-  GPxM_1 %>%
              dplyr::filter (genetic_map_chr =="3H")

mrk.GPxM_1.3H <- unique (GPxM_1_3H$snp)

(map.1.dt_3H <- map.1.dt %>%
    dplyr::filter  (mrk_names %in% mrk.GPxM_1.3H) %>%
    dplyr::arrange (chr))

#El **el cromosoma 3H** se asocia con los grupos de ligamiento 
(unique (map.1.dt_3H$chr))
  
#Cromosoma  **4H**

GPxM_1_4H <-  GPxM_1 %>%
  dplyr::filter (genetic_map_chr =="4H")

mrk.GPxM_1.4H <- unique (GPxM_1_4H$snp)

(map.1.dt_4H <- map.1.dt %>%
    dplyr::filter  (mrk_names %in% mrk.GPxM_1.4H) %>%
    dplyr::arrange (chr))


# El **cromosoma 4H** se asocia con los grupos de ligamiento  
 (unique (map.1.dt_4H$chr))
  
individuos
summary (Ill50K_Mapeo.qtl.e)

#Cromosoma  **5H**

GPxM_1_5H <-  GPxM_1 %>%
               dplyr::filter (genetic_map_chr =="5H")

mrk.GPxM_1.5H <- unique (GPxM_1_5H$snp)

(map.1.dt_5H <- map.1.dt %>%
    dplyr::filter  (mrk_names %in% mrk.GPxM_1.5H) %>%
    dplyr::arrange (chr))

# El **cromosoma 5H** se asocia con los grupos de ligamiento 
(unique (map.1.dt_5H$chr))


# Cromosoma  **6H**

GPxM_1_6H <-  GPxM_1 %>%
              dplyr::filter (genetic_map_chr =="6H")

mrk.GPxM_1.6H <- unique (GPxM_1_6H$snp)

map.1.dt_6H <- map.1.dt %>%
               dplyr::filter  (mrk_names %in% mrk.GPxM_1.6H) %>%
               dplyr::arrange (chr)

# El **cromosoma 6H** se asocia con los grupos de ligamiento 
(unique (map.1.dt_6H$chr))

# Cromosoma  **7H** 

GPxM_1_7H <-  GPxM_1 %>%
              dplyr::filter (genetic_map_chr =="7H")

mrk.GPxM_1.7H <- unique (GPxM_1_7H$snp)

(map.1.dt_7H <- map.1.dt %>%
    dplyr::filter  (mrk_names %in% mrk.GPxM_1.7H) %>%
    dplyr::arrange (chr))

# El **cromosoma 7H** se asocia con los grupos de ligamiento 
(unique (map.1.dt_7H$chr))


  
# Los unicos grupos que coinciden son el **L.3** que se tiene *tags* en  **1H** y **6H**, 
# y el **L.13** con **3H** y **6H** 
  
 
# Solo marcador **JHI_Hv50k_2016_23594** del grupo **L.3** esta en el **1H** mientras que 60 marcadores del **L.3**  estan en el **6H** 
  
# Solo marcador **JHI_Hv50k_2016_420337** del grupo **L.13** esta en el **6H**
  #mientras que 41 marcadores del **L.13**  estan en el **3H** 
  
#Los marcadores antes mencionados (**JHI_Hv50k_2016_23594** , **JHI_Hv50k_2016_420337**) se quitan del analisis

Ill50K_Mapeo.qtl.i  <- drop.markers (Ill50K_Mapeo.qtl.h,
                                     c("JHI-Hv50k-2016-23594","JHI-Hv50k-2016-420337"))


# El grupo **L.13** se vincula al **3H**  y el **L.3** al **6H** 

### Ordenamiento de marcadores dentro de los grupos de ligamiento
class (map.1)
unique (map.1$chr)
no_asignados <- c("L.2","L.8","L.18","L.21",
                     "L.22","L.24","L.25","L.26", "L.27", "L.28", "L.29")


LG_no_asignados <- map.1.dt %>%
                   dplyr::filter(chr %in% no_asignados)%>%
                   dplyr::mutate (marker_2 = mrk_names)


head (LG_no_asignados)

LG_no_asignados.1 <- LG_no_asignados%>%
                     inner_join (Ids_mkrs_eurekaSSR, by="marker_2") %>%
                     dplyr::rename (marker =marker_1)
                  


LG_no_asignados.2 <- LG_no_asignados.1%>%
                      inner_join (Ids_mkrs_pia, by="marker")

### genero el pmap de nuevo ##

map.2 <- pull.map (Ill50K_Mapeo.qtl.i, as.table=TRUE)

map.2.dt <- map.2 %>%
            dplyr::mutate(mrk_names = rownames(map.2)) %>%
            dplyr::select ("mrk_names", "chr" ,"pos")

str(map.2.dt)

map.2.dt$chr <- as.character (map.2.dt$chr)

# 1H
map.2.L1 <-  map.2.dt %>%
             dplyr::filter (chr=="L.1")

map.2.L4 <-  map.2.dt %>%
             dplyr::filter (chr=="L.4")

map.2.L6 <-  map.2.dt %>%
             dplyr::filter (chr=="L.6")

map.2.L7 <-  map.2.dt %>%
             dplyr::filter (chr=="L.7")

map.2.L9 <-  map.2.dt %>%
             dplyr::filter (chr=="L.9")

map.2.1H <- rbind (map.2.L1,
                   map.2.L4,
                   map.2.L6, 
                   map.2.L7,
                   map.2.L9)

## 2H
map.2.L10 <-  map.2.dt %>%
              dplyr::filter (chr=="L.10")

map.2.L11 <-  map.2.dt %>%
              dplyr::filter (chr=="L.11")

map.2.2H <- rbind (map.2.L10,
                   map.2.L11)

# 3H
map.2.L13 <-  map.2.dt %>%
              dplyr::filter (chr=="L.13")

map.2.L16 <-  map.2.dt %>%
              dplyr::filter (chr=="L.16")

map.2.L27 <-  map.2.dt %>%
              dplyr::filter (chr=="L.27")


map.2.3H <- rbind (map.2.L13,
                   map.2.L16,
                   map.2.L27)

# 4H
map.2.L12 <-  map.2.dt %>%
             dplyr::filter (chr=="L.12")

map.2.L17 <-  map.2.dt %>%
  dplyr::filter (chr=="L.17")

map.2.L19 <-  map.2.dt %>%
  dplyr::filter (chr=="L.19")

map.2.L26 <-  map.2.dt %>%
              dplyr::filter (chr=="L.26")

map.2.4H <- rbind (map.2.L12,
                   map.2.L17,
                   map.2.L19, 
                   map.2.L26)

# 5H
map.2.L5 <-  map.2.dt %>%
             dplyr::filter (chr=="L.5")

map.2.L15 <-  map.2.dt %>%
  dplyr::filter (chr=="L.15")

map.2.L20 <-  map.2.dt %>%
  dplyr::filter (chr=="L.20")

map.2.L22 <-  map.2.dt %>%
  dplyr::filter (chr=="L.22")

map.2.L29 <-  map.2.dt %>%
              dplyr::filter (chr=="L.29")

map.2.5H <- rbind (map.2.L5,
                   map.2.L15,
                   map.2.L20, 
                   map.2.L22,
                   map.2.L29)

# 6H
map.2.L3 <-  map.2.dt %>%
  dplyr::filter (chr=="L.3")

map.2.L25 <-  map.2.dt %>%
  dplyr::filter (chr=="L.25")

map.2.6H <- rbind (map.2.L3,
                   map.2.L25)

# 7H
map.2.L14 <-  map.2.dt %>%
  dplyr::filter (chr=="L.14")

map.2.L23 <-  map.2.dt %>%
  dplyr::filter (chr=="L.23")

map.2.L28 <-  map.2.dt %>%
  dplyr::filter (chr=="L.28")

map.2.7H <- rbind (map.2.L14,
                   map.2.L23,
                   map.2.L28)


####### asigno el H

map.2.1H <- map.2.1H %>%
            dplyr::mutate (chrom = "1H")

map.2.2H <- map.2.2H %>%
            dplyr::mutate (chrom = "2H")

map.2.3H <- map.2.3H %>%
  dplyr::mutate (chrom = "3H")

map.2.4H <- map.2.4H %>%
  dplyr::mutate (chrom = "4H")

map.2.5H <- map.2.5H %>%
  dplyr::mutate (chrom = "5H")

map.2.6H <- map.2.6H %>%
  dplyr::mutate (chrom = "6H")

map.2.7H <- map.2.7H %>%
  dplyr::mutate (chrom = "7H")

map.3 <- rbind (map.2.1H,
                map.2.2H,
                map.2.3H,
                map.2.4H,
                map.2.5H,
                map.2.6H,
                map.2.7H)
head (map.3)

map.3 <- map.3 %>%
          dplyr::select (mrk_names,chrom,pos)

#### 
pmap_III50K.3 <- map.3 
head(pmap_III50K.3)
length (pmap_III50K.3$mrk_names)

write_delim (pmap_III50K.3, path ="./Data/procdata/pmap_III50K.3.csv",
             delim = ",", na = "-")


### generar el geno
III50K.3 <- III50K.2 %>%
            dplyr::rename (mrk_names= Name)


geno_III50K.3 <- III50K.3 %>%
                 dplyr::filter (mrk_names %in%pmap_III50K.3$mrk_names)

geno_III50K.3 <-geno_III50K.3%>%
                 dplyr::select (-c(Chr,Position))


head(geno_III50K.3)
dim(geno_III50K.3)
(Ill50K_Mapeo.qtl.i$pheno$id)
genotypes <- as.character (unique (Ill50K_Mapeo.qtl.i$pheno$id))

class(genotypes)

genotypes [28] <- "IND37"
 xx <- c("mrk_names", genotypes )

geno_III50K.3 <- geno_III50K.3 %>%
                 dplyr::select (xx)


write_delim (geno_III50K.3, path ="./Data/procdata/geno_III50K.3.csv",
             delim = ",", na = "-")

### generar el pheno 
### generar la pheno
colnames (geno_III50K.3)

id <- colnames (geno_III50K.3) [2:76]
random <- rnorm(75, 0, 1)
length(random)
length(id)
pheno_III50K.3 <- rbind(id, random)
pheno_III50K.3 <- as.data.frame (pheno_III50K.3)
colnames (pheno_III50K.3) <- id
pheno_III50K.3  <- pheno_III50K.3 [-1,]

pheno_III50K.3 <- pheno_III50K.3 %>%
  dplyr::mutate (X="trait")%>%
  dplyr::select (X,everything())


write_delim (pheno_III50K.3, path ="./Data/procdata/pheno_III50K.3.csv",
             delim = ",", na = "-")

############

Ill50K_Mapeo.qtl.3 <- read.cross (format="tidy",
                                dir="./Data/procdata",
                                genfile ="geno_III50K.3.csv", 
                                mapfile ="pmap_III50K.3.csv" , 
                                phefile ="pheno_III50K.3.csv" , 
                                na.strings="-",
                                alleles=c("A","B"),
                                estimate.map=TRUE, 
                                convertXdata=TRUE, error.prob=0.0001,
                                map.function="kosambi",
                                #F.gen=6, 
                                crosstype ="riself")

Ill50K_Mapeo.qtl.3 <- jittermap(Ill50K_Mapeo.qtl.3)

plotMap(Ill50K_Mapeo.qtl.3)
summary (Ill50K_Mapeo.qtl.3)
plotMissing (Ill50K_Mapeo.qtl.3, main="")

cg.3 <- comparegeno(Ill50K_Mapeo.qtl.3)
hist(cg.3[lower.tri(cg.3)], breaks=seq(0, 1, len=101), xlab="No. matching genotypes")
rug(cg.3[lower.tri(cg.3)])

wh.3 <- which(cg.3 > 0.9, arr=TRUE)
wh.3 <- wh.3[wh.3[,1] < wh.3[,2],]
wh.3

####
#### Distorsion de segregación


gt <- geno.table (Ill50K_Mapeo.qtl.3)
gt [gt$P.value < 0.05/totmar(Ill50K_Mapeo.qtl.3),]
total <- nind (Ill50K_Mapeo.qtl.3)

gt <- geno.table(Ill50K_Mapeo.qtl.3)
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

gt.3 <- geno.table (Ill50K_Mapeo.qtl.3)
gt.3[gt.3$P.value < 0.05/totmar(Ill50K_Mapeo.qtl.3),]

todrop <- rownames(gt.3[gt.3$P.value < 1e-10,])


g3 <- pull.geno(Ill50K_Mapeo.qtl.3)
gfreq <- apply(g3, 1, function(a) table(factor(a, levels=1:2)))
gfreq <- t(t(gfreq) / colSums(gfreq))
par(mfrow=c(1,2), las=1)
for(i in 1:2)
  plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "BB")[i],
       ylim=c(0,1))


Ill50K_Mapeo.qtl.3 <- est.rf (Ill50K_Mapeo.qtl.3)


###################################################
### code chunk number 25: checkAlleles
###################################################
checkAlleles (Ill50K_Mapeo.qtl.3, threshold=5)

Ill50K_Mapeo.qtl.3  <- drop.markers(Ill50K_Mapeo.qtl.3 ,
                                    c("JHI-Hv50k-2016-51452"," JHI-Hv50k-2016-2369"))
     
Ill50K_Mapeo.qtl.3 <- est.rf (Ill50K_Mapeo.qtl.3)


checkAlleles (Ill50K_Mapeo.qtl.3, threshold=5)


Ill50K_Mapeo.qtl.3  <- drop.markers(Ill50K_Mapeo.qtl.3 ,
                                    c("JHI-Hv50k-2016-23694"))

checkAlleles (Ill50K_Mapeo.qtl.3, threshold=5)



Ill50K_Mapeo.qtl.4 <- mstmap (Ill50K_Mapeo.qtl.3, 
                              id = "id",bychr = TRUE, 
                              anchor = TRUE,
                              dist.fun = "kosambi", trace = FALSE)

#Ill50K_Mapeo.qtl.5 <- mstmap(Ill50K_Mapeo.qtl.3, 
                             #bychr = TRUE, id = "id",
                             #dist.fun = "kosambi", 
                             #anchor = TRUE,trace = TRUE)

nmar(Ill50K_Mapeo.qtl.4)

Ill50K_Mapeo.qtl.5 <- mergeCross(Ill50K_Mapeo.qtl.4, 
                                 merge = list(`1H` = c("1H.1","1H.2","1H.3", "1H.4", "1H.5"),
                                              `2H` = c("2H.1","2H.2","2H.8"),
                                              `3H` = c("3H.1","3H.11","3H.2", "3H.3"),
                                              `4H` = c("4H.1","4H.2","4H.3", "4H.4", "4H.6", "4H.7"),
                                              `5H` = c("5H.1","5H.2","5H.3"),
                                              `6H` = c("6H.1","6H.2","6H.7"),
                                              `7H` = c("7H.1","7H.11","7H.16","7H.2", "7H.22", "7H.25", 
                                                       "7H.27", "7H.31", "7H.41", "7H.5", "7H.6")))
plotMap( Ill50K_Mapeo.qtl.5, chr=c("1H", "2H", 
                                   "3H", "4H",
                                   "5H", "6H",
                                   "7H") )


Ill50K_Mapeo.qtl.6 <- subset( Ill50K_Mapeo.qtl.5, chr=c("1H", "2H", 
                                                           "3H", "4H",
                                                           "5H", "6H",
                                                           "7H")) 
plotMap( Ill50K_Mapeo.qtl.6)
heatMap(Ill50K_Mapeo.qtl.6, lmax = 70)

Ill50K_Mapeo.qtl.7 <- mstmap(Ill50K_Mapeo.qtl.6,id = "id",
                             bychr = FALSE, trace = TRUE, 
                             dist.fun = "kosambi")



plottMap ( Ill50K_Mapeo.qtl.7)


mapx <- est.map( Ill50K_Mapeo.qtl.7, map.function = "kosambi")

plotMap (Ill50K_Mapeo.qtl.7, mapx)


Ill50K_Mapeo.qtl.7 <- replace.map(Ill50K_Mapeo.qtl.7,mapx)

plotMap (Ill50K_Mapeo.qtl.7)


write.cross (Ill50K_Mapeo.qtl.7, format="tidy",
             filestem="./Data/procdata/Ill50K_Mapeo.qtl.7")


### comparar los phenos 

pheno_nico_1 <- read_delim (file="./Data/rawdata/pheno_nico_1a.txt", 
                            delim = "\t", na = "NA")


pheno_Ill50K_Mapeo.qtl.7 <- read_delim (file="./Data/procdata/Ill50K_Mapeo.qtl.7_phe.csv", 
                              delim = ",", na = "NA")


colnames (pheno_Ill50K_Mapeo.qtl.7)

pheno_Ill50K_Mapeo.qtl.7 <- pheno_Ill50K_Mapeo.qtl.7 %>%
                            dplyr::select (-c("IND3", "IND99"))

colnames (pheno_nico_1)


pheno_nico_2 <- pheno_nico_1 %>%
                dplyr::select (colnames (pheno_Ill50K_Mapeo.qtl.7))


write_delim (pheno_nico_2, path ="./Data/procdata/pheno_nico_2.csv",
             delim = ",", na = "-")


geno_Ill50K_Mapeo.qtl.7 <- read_delim (file="./Data/procdata/Ill50K_Mapeo.qtl.7_gen.csv", 
                                        delim = ",", na = "NA")

geno_Ill50K_Mapeo.qtl.8 <- geno_Ill50K_Mapeo.qtl.7 %>%
                           dplyr::select (-c("IND3", "IND99"))


write_delim (geno_Ill50K_Mapeo.qtl.8, path ="./Data/procdata/geno_Ill50K_Mapeo.qtl.8.csv",
             delim = ",", na = "-")



Ill50K_Mapeo.qtl.8 <- read.cross (format="tidy",
                                  dir="./Data/procdata",
                                  genfile ="geno_Ill50K_Mapeo.qtl.8.csv", 
                                  mapfile ="Ill50K_Mapeo.qtl.8_map.csv" , 
                                  phefile ="pheno_nico_2.csv" , 
                                  na.strings="-",
                                  genotypes=c("AA","BB"),
                                  alleles=c("A","B"),
                                  estimate.map=TRUE, 
                                  convertXdata=TRUE, 
                                  error.prob=0.0001,
                                  map.function="kosambi",
                                  #F.gen=6, 
                                  crosstype ="riself")

plotMap(Ill50K_Mapeo.qtl.8)
geno.image(Ill50K_Mapeo.qtl.8)
plot(Ill50K_Mapeo.qtl.8)
## Pheno Quality
pq.diagnostics (crossobj=Ill50K_Mapeo.qtl.8)

## Marker Quality
mq.diagnostics (crossobj=Ill50K_Mapeo.qtl.8,I.threshold=0.1,
                p.val=0.01,na.cutoff=0.1)

Ill50K_Mapeo.qtl.8$pheno
### QTL_SMA
QTL.result.V1 <- qtl.analysis (crossobj=Ill50K_Mapeo.qtl.8,
                                     step=0, method='SIM',  
                                     trait="v1", threshold="Li&Ji",
                               distance=30,  
                                      cofactors=NULL, window.size=30)

### QTL_SIM
QTL.result.V1 <- qtl.analysis ( crossobj=Ill50K_Mapeo.qtl.8, step=5,
                             method='SIM',trait="v1", threshold="Li&Ji",
                             distance=30,cofactors=NULL, window.size=30)


### QTL CIM
cofactors <- as.vector (QTL.result.V1$selected$marker)
QTL.result.v1 <- qtl.analysis ( crossobj=Ill50K_Mapeo.qtl.8, 
                                step=5,
                             method='CIM', trait="v1", threshold="Li&Ji", 
                             distance=30, cofactors=cofactors, window.size=30)
