#######################################################
# Construccion del Mapa de cebada                     #
# gaston quero                                        #
#18.12.2020                                           #
#######################################################

getwd()
setwd ("C:/Users/Usuario/OneDrive/Documentos/Mapa_cebada")


### Onemap
install.packages("devtools")
library(devtools)

install_github("augusto-garcia/onemap")

setRepositories(ind = 1:2)
install.packages("onemap", dependencies=TRUE)


system.file(package="onemap")

# paquetes
library (onemap)
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


#### cargo los datos

## esta es la referencia anterior de la posicion 
pia.mcd <- read_delim (file="./Data/rawdata/posiciones.pia.mcd.txt", 
                              delim = "\t", na = "-")

pia.mcd <- pia.mcd %>%
           dplyr::mutate (marker = str_c ("M", marker))

### cargo los datos iniciales 
pre.onemap.pia <- read_delim (file="./Data/rawdata/pre_onemap_todo2017Pia.txt", 
                         delim = "\t", na = "-")


### genero los id de las marcadores 
colnames (pre.onemap.pia) [-1]

IDs.genotypes.pia  <- str_c ( "ID", colnames (pre.onemap.pia) [-1])


### es la informacion de los estados alelicos
pre.one.1 <- pre.onemap.pia %>%
             dplyr::select (-ID)

### reeplazo los ID de los genotipoas 
colnames (pre.one.1) <- IDs.genotypes.pia 

# reemplazo los alelos
pre.one.1 [pre.one.1 =="A"] <- "a" 
pre.one.1 [pre.one.1 =="B"] <- "b" 

# estos son las etiquetas de los marcadores        
mrks <- pre.onemap.pia %>%
        dplyr::select (ID) %>%
        dplyr::mutate (ID = str_c ("*M", ID)) %>%
        dplyr::mutate (sgre = "A.B")

#uno las matrices
pre.one.2 <- bind_cols (mrks, pre.one.1 )

#exporto la matriz

write_delim (pre.one.2 , file ="./Data/procdata/pre.one.pia2.txt",
             delim = " ", na = "-")

#armo los metados
info.cross.pia <- "data type ri self" 

n.genotypes <- length (IDs.genotypes.pia)
n.mrk <- nrow (mrks)

info.cross.num.pia <-c (n.genotypes,n.mrk,0,0,0)

## hay que escribir los metadatos en el txt

### reingreso los datos de PIA
onemap_pia_riself <- read_onemap (dir="./Data/procdata", 
                                 inputfile = "pre.one.pia2.raw" )
str(onemap_pia_riself )
plot(onemap_pia_riself)


### Test de segregacion
ri_test <- test_segregation (onemap_pia_riself)

print (ri_test)

Bonferroni_alpha (ri_test)

plot (ri_test)

select_segreg (ri_test, distorted = TRUE)
print ( ri_test )

### los marcadores que no muestran distorsion de segregacion
no_dist <- select_segreg(ri_test, distorted = FALSE, numbers = TRUE) #to show the markers numbers without segregation distortion
no_dist


### los marcadores que muestran distorsion de segregacion
dist <- select_segreg(ri_test, distorted = TRUE, numbers = TRUE) #to show the markers numbers without segregation distortion
dist


# Estimating two-point recombination fractions
twopts_ri.LOD3 <- rf_2pts (onemap_pia_riself)

# Using only recombinations informations.
mark_all_ri.LOD3  <- make_seq (twopts_ri.LOD3, "all")

# Forming the groups
LGs_ri.Pia.LOD3  <- group (mark_all_ri.LOD3)

LGs_ri.Pia.LOD3

##### ESTA ES LA FUNCION ###########

run_onemap <-  function (dt= NULL, LOD=NULL,map.fun =NULL, data.base= NULL ){
  
  lod=LOD
  map.f =map.fun 
  db =data.base
  set_map_fun (type =  map.f)
  nlg <- dt$n.groups
  list.nlg <- seq (1:nlg)
  
  lapply (list.nlg, function (filt.lg){
    #filt.lg = 3
    LGx <- make_seq (dt, filt.lg)
    LGx.ord <- order_seq (input.seq = LGx, n.init = 5,
                          subset.search = "twopt",
                          twopt.alg = "rcd", THRES = 3)
    print( LGx.ord)
    
    LGx.ord.final <- make_seq ( LGx.ord, "force")
    print( LGx.ord.final)
    
    file.id <- str_c (db, "LG", filt.lg,"LOD",lod, map.f, sep="_")
    
    write_map ( map.list = LGx.ord.final, 
                file.out= str_c ("./Data/procdata/", file.id,".txt"))
    
    #LG.map <- LGx.ord.final
    
    #return (LG.map)
  })
  
  
  
}

#######################################

#### mapa pia LOD3 _ haldane

run_onemap (dt= LGs_ri.Pia.LOD3, LOD=3 ,map.fun ="haldane", data.base= "PIA" )

# reingreso los datos 
files.PIA.haldane.LOD3 <- dir(str_c("./Data/procdata"), pattern = "*_LOD_3_haldane.txt")

grupos.PIA.haldane.LOD3 <- bind_rows (lapply(files.PIA.haldane.LOD3 , function (filt.LG) {
  #filt.LG = "PIA_LG_1_LOD_3_haldane.txt" 
  dt <- read_delim (file = str_c("./Data/procdata/", filt.LG) , 
                    delim =" ", quote = "\"",
                    escape_backslash = FALSE,
                    escape_double = TRUE, 
                    col_names = FALSE, 
                    col_types = NULL,
                    locale = default_locale(), 
                    na = "NA")
  dt <- dt %>%
    dplyr::mutate (id.file=filt.LG)%>%
    dplyr::mutate (x=filt.LG) %>%
    tidyr::separate(x , c(NA, NA, "LG",NA ,NA, NA ), sep ="_")%>%
    dplyr::rename (marker=X2)%>%
    dplyr::rename (pos = X3)%>%
    dplyr::select ( LG,marker,pos, id.file)
  
}))

#### comparo los grupos de ligamiento 

unique (pia.mcd$LG)
dt.ref = pia.mcd 
dt.com = grupos.PIA.haldane.LOD3 

run_compare <- function (dt.ref= NULL, dt.com=NULL){
  
  list.ref <- unique (dt.ref$LG)
  df.x <- bind_rows (lapply ( list.ref, function (filt.LG){
    #filt.LG= "1Ha"
    
    dt.1 <- dt.ref %>%
      dplyr::filter (LG ==filt.LG )
    
    dt.2 <- dt.1 %>%
              inner_join( dt.com, by="marker" )
    
    print(dt.2)
    
  }))
return (df.x )
}

PIA.LOD.3 <- run_compare (dt.ref =pia.mcd, dt.com =grupos.PIA.haldane.LOD3 )

list.LG <- unique (PIA.LOD.3$LG.y)


run_plot <- function (dt = NULL){
  
  
  
  
}



PIA.LOD.mcd <- PIA.LOD.3 %>%
             dplyr::select (marker, LG.x,  pos.x ) %>%
             dplyr::rename (mar = marker) %>%
             dplyr::rename (chr = LG.x) %>%
             dplyr::rename (pos = pos.x) %>%
             dplyr::mutate (data="Pia")

PIA.LOD.onemap <- PIA.LOD.3 %>%
                  dplyr::select (marker, LG.y,  pos.y ) %>%
                  dplyr::rename (mar = marker) %>%
                  dplyr::rename (chr = LG.y) %>%
                  dplyr::rename (pos = pos.y)%>%
                  dplyr::mutate (data="onemap") 

XX <- bind_rows (PIA.LOD.mcd, PIA.LOD.onemap )

LG.1H <- PIA.LOD.3 %>%
         dplyr::filter (LG.x == "1Ha"  )

PIA.LOD.mcd <- LG.1H %>%
  dplyr::select (marker, LG.x,  pos.x ) %>%
  dplyr::rename (mar = marker) %>%
  dplyr::rename (chr = LG.x) %>%
  dplyr::rename (pos = pos.x) %>%
  dplyr::mutate (data="Pia")

PIA.LOD.onemap <- LG.1H %>%
  dplyr::select (marker, LG.y,  pos.y ) %>%
  dplyr::rename (mar = marker) %>%
  dplyr::rename (chr = LG.y) %>%
  dplyr::rename (pos = pos.y)%>%
  dplyr::mutate (data="onemap") 




ggpaired (XX, 
          x = "chr", y = "pos",
          id = "mar",
          #color = "genotype", 
          line.color = "black", 
          line.size = 0.5,
          width =0,
          point.size = 3,
          xlab = "data",
          ylab = "pos",
          label= "mar",
          repel = TRUE,
          #palette =  c("dodgerblue3", "firebrick2","cadetblue2", 
                       #"darkgoldenrod1")
          ) #+
  geom_hline (yintercept = 1,  linetype = 2)




ggpaired (LG.1H, cond1 = "LG.x", cond2 = "LG.y",
         fill = "marker", palette = "jco")

chr=c("1Ha" ,"1")
linkmap(XX ,chr=c("1Ha" ,"1"))

plotMap(PIA.LOD.mcd, PIA.LOD.onemap)

# Ordering markers within linkage groups
set_map_fun (type = "haldane")

### vamos por grupo 

LG_Pia.LOD3 <- make_seq (LGs_ri.Pia.LOD3, "all")

### grupo 1 ### 

LG1_Pia.LOD3 <- make_seq (LGs_ri.Pia.LOD3, 1)
LG1_Pia.LOD3$seq.num

LGs_ri.Pia.LOD3$marnames

LG1_Pia.LOD3_ord <- order_seq (input.seq = LG1_Pia.LOD3, n.init = 5,
                               subset.search = "twopt",
                               twopt.alg = "rcd", THRES = 3)

LG1_Pia.LOD3_force <- make_seq (LG1_Pia.LOD3_ord, "force")
class (LG1_Pia.LOD3_force )
marker_type (LG1_Pia.LOD3_force)

marnames <- colnames(get(LG1_Pia.LOD3_force$data.name, pos = 1)$geno)[LG1_Pia.LOD3_force$seq.num]

G1.pia.mcd <- pia.mcd %>%
              dplyr::filter (marker%in%marnames)





dt = LGs_ri.Pia.LOD3
LOD =3
map.fun = "haldane"
data.base= "pia" 











(LOD_sug <- suggest_lod (onemap_pia_riself))

Map.LG1.Pia.LOD3_force  <- map(LG1_Pia.LOD3_force)

str(LG1_Pia.LOD3_force)








map.list.haldane <- vector("list", LGs_ri.Pia.LOD3$n.groups)

for(i in 1:LGs_ri.Pia.LOD3$n.groups){
  ##create linkage group i
  LG.cur <- make_seq(LGs_ri.Pia.LOD3,i)
  ##ordering
  map.cur <- order_seq (LG.cur, subset.search = "sample")
  ##assign the map of the i-th group to the maps.list
  map.list.haldane [[i]] <- make_seq(map.cur, "force")
}


twopts_ri.LOD_sug <- rf_2pts (onemap_pia_riself, LOD =LOD_sug)





mark_no_dist_ri.LOD3 <- make_seq (twopts_ri.LOD3, no_dist)



mark_all_ri.LOD4.5     <- make_seq (twopts_ri.LOD_sug, "all")

mark_no_dist_ri.LOD4.5 <- make_seq (twopts_ri.LOD_sug, no_dist)

### 
LGs_ri.Pia.LOD3  <- group (mark_all_ri.LOD3)

LGs_ri.Pia.LOD3

LGs_ri.Pia.LOD3.no.dist  <- group (mark_no_dist_ri.LOD3)

LGs_ri.Pia.LOD3.no.dist 


LGs_ri.Pia.LOD4.5  <- group (mark_all_ri.LOD4.5, LOD = LOD_sug, max.rf = 0.5)

LGs_ri.Pia.LOD4.5

LGs_ri.Pia.LOD4.5.no.dist  <- group (mark_no_dist_ri.LOD4.5, LOD = LOD_sug, max.rf = 0.5)

LGs_ri.Pia.LOD4.5.no.dist 

#Ordering markers within linkage groups

set_map_fun(type = "haldane")

### vamos por grupo 
LG1_Pia.LOD3 <- make_seq (LGs_ri.Pia.LOD3, 1)
LG1_Pia.LOD3

LG1_ser_Pia.LOD3 <- seriation (LG1_Pia.LOD3)
LG1_rcd_Pia.LOD3 <- rcd (LG1_Pia.LOD3)
LG1_rec_Pia.LOD3 <- record (LG1_Pia.LOD3)
LG1_ug_Pia.LOD3  <- ug (LG1_Pia.LOD3)

LG1_Pia.LOD3_ord <- order_seq (input.seq = LG1_Pia.LOD3, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3)

LG1_Pia.LOD3_ord 

LG1_Pia.LOD3_safe <- make_seq (LG1_Pia.LOD3_ord, "safe")

LG1_Pia.LOD3_force <- make_seq (LG1_Pia.LOD3_ord, "force")


LG1_Pia.LOD3_ord.x <- order_seq(input.seq = LG1_Pia.LOD3, n.init = 5,
                        subset.search = "twopt",
                        twopt.alg = "rcd", THRES = 3,
                        touchdown = TRUE)


LG1_Pia.LOD3_force.x <- make_seq (LG1_Pia.LOD3_ord.x, "force")

ripple_seq (LG1_Pia.LOD3_force.x , ws = 5, LOD = 3)

LG1_Pia.LOD3_force.x$seq.num

rf_graph_table (LG1_Pia.LOD3_force.x)

######### 

draw_map2 (LG1_Pia.LOD3_force.x, main = "Only linkage information", 
          group.names = "LG1", output = "tiff")



######

## Not run: 
data(mapmaker_example_f2)
twopt <-rf_2pts(mapmaker_example_f2)
lg <-group(make_seq(twopt, "all"))

##"pre-allocate" an empty list of length lg$n.groups (3, in this case)
maps.list<-vector("list", lg$n.groups)

for(i in 1:lg$n.groups){
  ##create linkage group i
  LG.cur <- make_seq(lg,i)
  ##ordering
  map.cur<-order_seq(LG.cur, subset.search = "sample")
  ##assign the map of the i-th group to the maps.list
  maps.list[[i]]<-make_seq(map.cur, "force")
}

##write maps.list to "mapmaker_example_f2.map" file
write_map(maps.list, "mapmaker_example_f2.map")

##Using R/qtl
##you must install the package  'qtl'
##install.packages("qtl")

require(qtl)
file<-paste(system.file("example",package="onemap"),"mapmaker_example_f2.raw", sep="/")
dat1 <- read.cross("mm", file=file, mapfile="mapmaker_example_f2.map")
newmap <- est.map(dat1, tol=1e-6, map.function="kosambi")

(logliks <- sapply(newmap, attr, "loglik"))
plot.map(dat1, newmap)



##Using R/qtl to generate QTL Cartographer input files (.map and .cro)
write.cross(dat1, format="qtlcart", filestem="mapmaker_example_f2")


## End(Not run)



