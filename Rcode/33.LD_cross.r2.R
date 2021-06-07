####################################################################
# Codigo para el LD y LD decay de diferentes poblaciones de mapeo  #
#                                                                  #
# Gaston Quero - Sebastian Simomndi                                #
#                                                                  #
# 29/5/2021                                                        #
#                                                                  #
####################################################################

getwd ()
setwd ("R:/Mapa_genetico_cebada")


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
library(LDheatmap)
library(genetics)
library(assertive)

######### 
# primero cargar la informacion que se tenga del mapa.

mapa.cebada.bill  <- read_delim (file = "./Data/rawdata/8.Bill_thomas_table_3.txt" ,
                      col_names = TRUE, delim = "\t", na = "NA")

## La informacion de la posicion fisica debe estar en Mb 


mapa.cebada.bill <- mapa.cebada.bill %>%
                       dplyr::mutate (Mb = bp * 1e-6)
   



## Se cargan los datos de la poblacion de mapeo con la cual se va a hacer la caida del LD.

EEMAC.cross <- read.cross (format="csv",
                           dir="./Data/rawdata", file= "1.cross_EEMAC.1.csv",
                           na.strings="NA", 
                           genotypes=c("0","1"),
                           #alleles=c("0","1"),
                           estimate.map=FALSE,
                           crosstype ="dh",	
                           convertXdata=TRUE, error.prob=0.0001)

summary (EEMAC.cross )

### se extrae el Mapa de la poblacion de mapeo.

mapa.EEMAC.cross <- pull.map (EEMAC.cross, as.table=TRUE)

row.names (mapa.EEMAC.cross)

mapa.EEMAC.cross <- mapa.EEMAC.cross  %>%
                    dplyr::mutate (snp = row.names (mapa.EEMAC.cross)) 


#row.names (mapa.EEMAC.cross) <- NULL

### Veo cuantos de los marcadores de la poblacion de mapeo estan en la info del mapa general


mapa.EEMAC.cross.bill <- mapa.EEMAC.cross %>%
                         left_join(mapa.cebada.bill, by="snp")
                         

## Extraer los maracadores que tienen las posiciones en bp y en cM

mapa.EEMAC.cross.bill.Mb <- mapa.EEMAC.cross.bill  %>%
                            dplyr::filter (!is.na(Mb))

pos.unique <- unique (mapa.EEMAC.cross.bill.Mb$pos)

mapa.EEMAC.cross.cM <- mapa.EEMAC.cross %>%
                       dplyr::filter (pos %in% pos.unique)


#### hago una lista de los marcadores que debo sacar del mapa de la poblacion de mapeo

head ( mapa.EEMAC.cross.bill.Mb)

list.mark <- unique (mapa.EEMAC.cross.bill.Mb$snp)

drop.list.mark <- mapa.EEMAC.cross %>%
                  dplyr::filter (! snp %in%  list.mark)

head (drop.list.mark)

### saco los marcadores de la poblacion de mapeo

EEMAC.cross.1 <- drop.markers (EEMAC.cross, c(drop.list.mark$snp))

summary (EEMAC.cross.1)

write.cross (EEMAC.cross.1, format= "tidy",
            filestem= "./Data/procdata/2.cross_EEMAC.1")

### cargo de nuevo la matriz pero en formato tidy 

map.cross.2 <-  read_delim ( file ="./Data/procdata/2.cross_EEMAC.1_map.csv" , delim= ",",
                             quote = "\"", escape_backslash = FALSE,escape_double = TRUE,
                             col_names = TRUE, col_types = NULL, na = "NA")
head(map.cross.2)

map.cross.2 <- map.cross.2 %>%
               dplyr::rename (snp =  X1 )


mapa.EEMAC.cross.bill.Mb.select <- mapa.EEMAC.cross.bill.Mb %>%
                                   dplyr::select (c ("chr", "snp", "Mb"))

head (map.cross.2)

map.cross.2.Mb <- map.cross.2 %>%
                  dplyr::inner_join (mapa.EEMAC.cross.bill.Mb.select, by="snp") %>%
                  dplyr::select (c ("snp","chr.x", "Mb")) %>%
                  dplyr::rename (chr = chr.x) %>%
                  dplyr::rename (pos = Mb)

map.cross.2.Mb  <- write_delim (map.cross.2.Mb, 
                               file ="./Data/procdata/2.cross_EEMAC.1_map.Mb.csv", delim = ",")

#### Antes de cargar de nuevo borrar la etiqueta de la columna "snp"

EEMAC.cross.1.Mb <- read.cross (format="tidy",
                               dir="./Data/procdata",
                               genfile = "2.cross_EEMAC.1_gen.csv",
                               mapfile = "2.cross_EEMAC.1_map.Mb.csv", 
                               phefile =  "2.cross_EEMAC.1_phe.csv",
                               na.strings="-", 
                               genotypes=c("AA","BB"),
                               crosstype ="dh",
                               estimate.map=FALSE, error.prob=0.0001)
  
summary (EEMAC.cross.1.Mb)

### hay que mejor esto
plotMap (EEMAC.cross.1.Mb, horizontal=FALSE, show.marker.names=FALSE)

## Poblacion de mapeo (GWAS) soja

soja.cross <- read.cross(format="csv",
                                 dir="./Data/rawdata", file= "9.cross.data.z14_15.csv",
                                 na.strings="-", 
                                 genotypes=c("AA","BB"),
                                 #alleles=c("0","1"),
                                 estimate.map=FALSE,
                                 crosstype ='dh',
                                 convertXdata=TRUE, 
                                 error.prob=0.0001)


# Nota: se debe unificar los formatos a tibble todo 


## Hay que adaptar una funcion para graficar los mapas. 

## info de la poblacion
summary (soja.cross)

plotMap (soja.cross, horizontal=FALSE, show.marker.names=FALSE)


### los argumentos de la funcion para correr la funcion por dentro

#data.cross = soja.cross 
#heterozygotes = FALSE
#data.tibble = soja.cross.tibble
#distance.unit = "Mb"
#id.cross = "GWAS.soja"

############# CEBADA 
##  funcion jittermap para 
# evitar la misma posicion de los marcadores.


#EEMAC.cross.1 <- jittermap (EEMAC.cross.1, amount=1e-6)



data.cross = EEMAC.cross.1.Mb
heterozygotes = FALSE
distance.unit = "Mb"
id.cross = "EEMAC.cross.1.Mb"



# 1 bp => 1e-6 Mb

#filt.chr = 1
#data.tibble [1:2, 3610:3618]

## esta funcion calcula el r2 y grafica el heatmap.

run_plot_heatmap_LD <- function (id.cross = NULL , 
                                 data.cross = NULL,  
                                 heterozygotes = FALSE,
                                 distance.unit = NULL) {

# Nota: se debe crear la estructura del directorio

print (str_c("Se encontraron " , nchr (data.cross), " grupos de ligamientos")) # verifico el numero de cromosomas

if   ( distance.unit != "cM" & distance.unit != "Mb") {
  
  stop (str_c ("Debe definir una unidad de distancia valida"))
  
  
}

if (is.null (distance.unit)) {
  
  stop (str_c ("Falta definir unidad de distancia"))
  
}   

if   ( distance.unit == "Mb"  ) {
  
  print (str_c ("El mapa es un mapa fisico"))
  
}

if   ( distance.unit == "cM") {
  print (str_c ("El mapa es un mapa genetico"))

}

#  numero maximo de grupos de ligamientos

nchr.max <- nchr (data.cross)  

list.chr <- 1:nchr.max

### este es temporal para corre el lapply

#filt.chr = 1
###########  

lapply(list.chr, function (filt.chr){
  
  print (str_c("Estimando LD LG= " ,filt.chr))
  
## hago un subset de cada cromosoma
  
  crossobj <- subset(data.cross, chr=filt.chr)
  
  chr <- filt.chr
  
#  Esta seccion es para formatear la matriz para para LDheatmap.
#  Nota: ver si con los nuevos formatos cambia o no
  
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
### Esta es la funcion que calcula el LD, verificar la documentacion
##
##############

#######################3
# Nota : hay que agregar barra de progreso  
  
    ld <- LD (genos)       
  
map.cross <- pull.map ( crossobj, as.table = TRUE)
   
# plot LD heatmap
 
   if   ( distance.unit == "Mb" ) {
    
    plot.hm <- LDheatmap ( ld$"R^2", 
                           genetic.distances= map.cross$pos, 
                           distances="physical",
                           title=str_c(id.cross, "_Pairwise LD_LG.",filt.chr),
                           color = colorRampPalette(c("red4", "red","orangered", "orange","yellow1",  "blue4"))(60))
    
    }
   
   if   ( distance.unit == "cM") {
     
     plot.hm <- LDheatmap ( ld$"R^2", genetic.distances= map.cross$pos, 
                            distances="genetic",
                            title=str_c(id.cross, "Pairwise LD_LG.", filt.chr),
                            color = colorRampPalette(c("red4", "red","orangered", "orange","yellow1",  "blue4"))(60))
   }
   
### esta es la matriz que sale de la funcion LD
   
  LD.cross.matrix <- plot.hm$LDmatrix

  write.table (LD.cross.matrix , 
               file = str_c("./Data/procdata/",id.cross, "_LD.cross.matrix_", filt.chr,".txt"),
               append = FALSE, quote = TRUE, sep = ",",
               eol = "\n", na = "NA", dec = ".", row.names = TRUE,
               col.names = TRUE)
  
### genero el df de los marcadores y sus distancias
  
 cross.tibble_dist <- map.cross %>%
                            dplyr::mutate ( mrks  = rownames(map.cross)) %>%
                            dplyr::select ( mrks, pos)
  

 cross.tibble_dist <- as_tibble (cross.tibble_dist)

dt <-  cross.tibble_dist 
list.pos <- (unique (dt$pos))
list.mrks <- (unique (dt$mrks))

#### aca empieza la distancia 

print (str_c("Estimando diff.dist LG= " ,filt.chr))


#start.time <- Sys.time()
dt.diff.dist <- bind_rows (lapply (list.pos, function (filtro.x1) { 
  
   #filtro.x1 =  0.277230 
   #print (filtro.x1)
  
    dt.x1 <- dt %>% 
             dplyr::filter (pos == filtro.x1)
  
    dt.x2 <- dt 
    
    dt.z <- data.frame (dt.x1,  dt.x2) 
    
    dt.z <- dt.z  %>%
            dplyr::mutate (diff.dist = abs (pos - pos.1)) %>%
            dplyr::select ("mrks", "mrks.1", "diff.dist" ) %>%
            dplyr::rename (mrk.1 = mrks) %>% 
            dplyr::rename (mrk.2 = mrks.1) 
  
}))

#end.time <- Sys.time()
#time.taken <- end.time - start.time
#time.taken

#start.time <- Sys.time()

df.LD.decay <- bind_rows ( lapply (list.mrks, function (filt.mrk) {
  
  #filt.mrk= "JHI-Hv50k-2016-270"

  #print (filt.mrk)
  
  dt.diff.dist.mrk <-  dt.diff.dist %>% 
                       dplyr::filter (mrk.1 == filt.mrk)  # las diff de distancias de un marcador contra el resto
  
  colnames (LD.cross.matrix ) <- list.mrks
  rownames (LD.cross.matrix)  <- list.mrks
 
  LD.cross <-  as_tibble (LD.cross.matrix) # convierto en tibble la matriz de R2 que calcula ld
 
  LD.cross <- LD.cross %>%
              dplyr::mutate (mrk.id = list.mrks)%>% ## agrego una columna
              dplyr::select (mrk.id, everything())
  
  LD.mrk.1 <- LD.cross %>%                          # me quedo con el marcador para el que calcule las distancias
              dplyr::filter (mrk.id== filt.mrk)     # aca estan los r2 de ese marcador y el resto
  
  id.gather <- colnames (LD.mrk.1) [-1]
  
  LD.cross.2 <- LD.mrk.1 %>%           # traspongo el tibble de R2
                gather (all_of (id.gather), key="mrk.2" , value = "R2") %>%
                dplyr::rename (mrk.1 = mrk.id) 
  
  LD.dist.cross.3 <- LD.cross.2 %>% # aca uno los datos de distancia y R2
                    dplyr::inner_join ( dt.diff.dist.mrk ,  LD.cross.2, by= c("mrk.1", "mrk.2"))
  
  return (LD.dist.cross.3)
  
}))

#end.time <- Sys.time()
#time.taken <- end.time - start.time
#time.taken


df.LD.decay <- df.LD.decay  %>%
               dplyr::mutate (LG = str_c ("lg.", filt.chr)) %>%
               dplyr::arrange (diff.dist) 

### verificar estas medidas

write_delim (df.LD.decay  , file =str_c("./Data/procdata/", id.cross, "LD.decay.LG", filt.chr,".txt"),
             delim = ",", na = "NA")

if   (distance.unit == "Mb" ) {
  
  df.LD.decay.NA <- df.LD.decay %>%
                    dplyr::filter (!is.na (R2))
    
  x2 <- max (df.LD.decay$diff.dist, na.rm = TRUE)
  
if (x2 > 100) {

 plot.LD.decay <-  ggscatter (df.LD.decay.NA, x = "diff.dist", y = "R2",
                   title =str_c(id.cross,".LD.decay.LG_", filt.chr)) +
                   scale_x_continuous (name="(Distance (Mb)",
                        breaks =seq(0,x2,100),
                        #labels=NULL,
                        limits=c(0, x2))
}

  if (x2 < 100) {
    
    plot.LD.decay <-  ggscatter (df.LD.decay.NA, x = "diff.dist", y = "R2",
                                 title =str_c(id.cross,".LD.decay.LG_", filt.chr)) +
      scale_x_continuous (name="(Distance (Mb)",
                          breaks =seq(0,x2,10),
                          #labels=NULL,
                          limits=c(0, x2))
  }
  
  plot.LD.decay  %>% 
    ggexport(filename = str_c("./Figures/ plot.LD.decay.",chr,".png"))
  
  print (plot.LD.decay)

}

#if   ( distance.unit == "cM" ) {
  
  #png (filename = str_c("./Figures/",id.cross,".LD.decay.LG_", filt.chr,".png"),
  #     width = 480, height = 480, units = "px", pointsize = 12,
   #    bg = "white", res = NA)
  
#  plot (x = df.LD.decay$diff.dist , y = df.LD.decay$R2, 
 #       main=str_c(id.cross,".LD.decay.LG_", filt.chr),
#        pch = 20, 
 #       type ="n",
  #      xaxt="none",
   #     yaxt="none",
    #    axes = F,
     #   xlim = c(0, max (df.LD.decay$diff.dist)), 
      #  ylim = c(0, max (df.LD.decay$R2, na.rm = TRUE)),
       # ylab = expression(LD ~ (r^2)),
        #xlab = expression(Distance ~ (cM))) 
  
 # axis(side = 2, las = 1)
  #x2 <- max (df.LD.decay$diff.dist, na.rm = TRUE)
  #axis (side=1,at=seq(0,x2,1),las = 1)
  
  
#  points (df.LD.decay$diff.dist, df.LD.decay$R2, 
          #pch = 20, cex=1.5, col="gray28") 
 # box()
  #dev.off()
  
  #plot (x = df.LD.decay$diff.dist, y = df.LD.decay$R2, 
   #     main=str_c(id.cross,".LD.decay.LG_", filt.chr),
        #pch = 20, 
        #type ="n",
        #xaxt="none",
        #yaxt="none",
        #axes = F,
        #xlim = c(0, max (df.LD.decay$diff.dist)), 
        #ylim = c(0, max (df.LD.decay$R2, na.rm = TRUE)),
        #ylab = expression(LD ~ (r^2)),
        #xlab = expression(Distance ~ (cM))) 
  
#  axis(side = 2, las = 1)
  #x2 <- max (df.LD.decay$diff.dist, na.rm = TRUE)
  #axis (side=1,at=seq(0,x2,1),las = 1)
  
  
  #points (df.LD.decay$diff.dist,df.LD.decay$R2, 
          #pch = 20, cex=1.5, col="gray28") 
  #box()
#}#

})

}

start.time <- Sys.time()

LD.EEMAC.cross.1.Mb <- run_plot_heatmap_LD (id.cross ="EEMAC.cross.1.Mb" , 
                                            data.cross = EEMAC.cross.1.Mb,  
                                            heterozygotes = FALSE,
                                            distance.unit = "Mb")
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

# Cargo las matrices de datos ############

#start_time <- Sys.time()

# Time difference of 1.000327 min
## esto lo tenes setear
index.chr <- 1:7

# se cargan las matrices generadas antes
# se genear una lista con las tibbles que exporta la funcion run_plot_heatmap_LD

GENO.crom <- lapply (index.chr, function (filtro) {
  
  G <- read_delim (file=paste("./Data/procdata/EEMAC.cross.1.bpLD.decay.LG",filtro,".txt",
                              sep=""), delim = ",",
                   na = "NA", quote = "\"",col_names = TRUE)
  
  
  G <- G %>%
    dplyr::mutate (chrom = filtro)
  
  #G <- G %>%
  #    dplyr::mutate (chr = str_c("chr_",filtro))
  print (G)
  return (G)
})

df.geno.crom <- do.call (rbind, GENO.crom)

head (df.geno.crom)
summary (df.geno.crom)


unique ( df.geno.crom$chrom )
###############################
#### funciones para analisis de LD
max (df.geno.crom$diff.dist)


#### pruebo solo el LG.1

#df.geno.crom.2 <- df.geno.crom %>%
 #                 dplyr::filter (chrom == 2)

#max (df.geno.crom.2$diff.dist)

#########
# Argumentos anteriores

#id.cross = "EEMAC.cross.1"
#data= df.geno.crom
#ini = 1

run_table_quantiles <- function (id.cross=NULL,data= NULL, ini = 1, l1= 0.25, l2=0.5, l3=0.75, seq1= NULL, distance.unit = NULL) {
  
  if   ( distance.unit != "cM" &  distance.unit != "Mb") {
    
    stop (str_c ("Debe definir una unidad de distancia valida"))
    
    
  }
  
  if (is.null (distance.unit)) {
    
    stop (str_c ("Falta definir unidad de distancia"))
    
  }
  
  if   ( distance.unit == "Mb") {
    
    print (str_c ("El mapa es un mapa fisico"))
    
  }
  
  if   ( distance.unit == "cM") {
    print (str_c ("El mapa es un mapa genetico"))
    
  }
  
  
  list.chrom <- unique (data$chrom)
  #filt.crom =2  ###
  
  df.qq <- bind_rows (lapply (list.chrom, function (filt.crom){
    
    print(filt.crom)
    
    dat.1 <- data %>%
      dplyr::filter (chrom == filt.crom) %>%
      dplyr:::filter (R2 != "NA")
    
    list.dist <- seq1
    
    # filt.dist = 1 ######!!!!!!!
    
    df.hist.plot <- bind_rows (lapply (list.dist, function (filt.dist){
      
      print (filt.dist)
      
      ############### hay que revisar estos numeros ##############
      # if ( distance.unit == "Mb" ) {
      
      x.dist.Mb <- dat.1 %>%
        dplyr::filter (diff.dist <= filt.dist * 1e5)
      
      bin <- (filt.dist* 1e5)/1e6
      
      x.dist.Mb.1 <- x.dist.Mb %>%
        dplyr::mutate (bin = str_c (bin, "Mb"))
      
      #return (x.dist.Mb.1)
      
      #}
      
      #########################################################################
      
      ############### hay que revisar estos numeros ##############
      #if   ( distance.unit == "cM" ) {
      
      # x.dist.cM <- dat.1 %>%
      #  dplyr::filter ( diff.dist <= filt.dist)
      
      #  bin <- filt.dist
      
      # x.dist.cM <-x.dist.cM %>%
      #            dplyr::mutate (bin = str_c (bin, "cM"))
      #}
      
      ############### hay que revisar estos numeros ##############
      
    }))
    
    # ggv <-  ggviolin(df.violin.plot, x = "bin", y = "R2", title = filt.crom,
    #          color = "black", fill = "azure")
    
    #print (ggv)
    
    
    ggh <- gghistogram (df.hist.plot, y = "..density..", x = "R2",
                        facet.by = "bin",title =filt.crom,
                        #fill = "lightgray",
                        fill = "bin", palette =   "RdBu",
                        add = "median", rug = TRUE)
    print (ggh)
    
    ggh  %>%
      ggexport(filename = str_c("./Figures/ggh_",id.cross, "_", filt.crom,".png"))
    
    
    qqunif <- ggplot(df.hist.plot, aes(R2, colour =  bin , fill = bin), size= 2) + stat_ecdf() +
      theme_bw()+
      stat_function(fun=punif,args=list(0,1))+
      scale_color_manual(values=c("green", "black", "blue", "red"))+
      labs(title=str_c("ECDF and theoretical CDF","_", filt.crom)) +
      labs(y = "Theoretical Quantiles", x = "Sample Quantiles")
    
    print (qqunif)
    
    qqunif %>%
      ggexport(filename = str_c("./Figures/qqunif_",id.cross, "_", filt.crom,".png"))
    
    # if   ( distance.unit == "cM" ) {
    
    # Nota: aca los estoy estoy tomando cada 0.5 cM hay que ver si esto reemplaza a 0.1 Mb
    
    #list.seqCM <- seq (from = 0.5, to = 20, by = 0.5)
    #filt.distCM <- 0.5
    #df.QQ.seq.cM <- bind_rows (lapply (list.seqCM, function (filt.distCM){
    
    # x.ini.cM <- dat.1 %>%
    #  dplyr::filter (diff.dist <= ini * filt.distCM)
    
    #QQ <- quantile (x.ini.cM$R2)
    
    #l1 <- l1
    #l2 <- l2
    #l3 <- l3
    #  q1 <- quantile (QQ, l1)
    # q2 <- quantile (QQ, l2)
    #  q3 <- quantile (QQ, l3)
    
    # XX <- data.frame( HLE=round (q1 [[1]],2) ,  H1= round (q2 [[1]],2), HLD = round (q3 [[1]],2), chrom = filt.crom)
    
    #dt.QQ.seq.cM <- XX %>%
    # dplyr::mutate (inter.cM = filt.distCM ) %>%
    #dplyr::select (chrom, inter.cM, HLE,   H1,  HLD )
    
    #}))
    
    
    #df.QQ.seq.cM.long <- df.QQ.seq.cM %>%
    # pivot_longer( ! c(chrom,inter.cM), names_to = "cat.LD", values_to = "unqq.R2")
    
    
    #qq.scatt <- ggscatter (df.QQ.seq.cM.long, x = "inter.cM", y = "unqq.R2",
    #                      title = str_c("Caida de qqR2 en funcion del bin", id.cross, "_", filt.crom),
    #                     color = "cat.LD", shape = "cat.LD",
    #                    palette = c("navyblue", "gray48", "darkorange"),
    #                   add = "loess", rug= TRUE)
    
    #print (qq.scatt)
    
    #qq.scatt  %>%
    # ggexport(filename = str_c("./Figures/qq.scatt_",id.cross, "_", filt.crom,".png"))
    
    #  }
    
  }))
  return (df.qq)
}

#l2=0.5
#l3=0.75
#seq1 = c(1, 10, 50,100)
#distance.unit = "Mb"


start.time <- Sys.time()

EEMAC.cross.1.qq <- run_table_quantiles  ( id.cross = "EEMAC.cross.1",
                       data= df.geno.crom,
                       ini = 1,
                       l1= 0.25,
                       l2=0.5,
                       l3=0.75,
                       seq1 = c(1, 10, 50, 100),
                       distance.unit = "Mb")

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



#data = df.geno.crom 
#chr= 1
#ini=10
#l1= 0.25
#l2=0.5
#l3=0.75
#keep.Mb = 0.1
#prob.HLMK = 51
#prob.HUMK = 50


run_freq_decay <- function (data, chr=NULL, ini = NULL, l1= 0.25, l2=0.5, l3=0.75, keep.Mb =NULL,
                            prob.HLMK = NULL, 
                            prob.HUMK = NULL){
  
  # el primer index p
  #index.chrom <- unique (data$chr)
  
  # control de la clases de los objetos
  assert_is_data.frame (data)
  assert_is_numeric (chr)
  #assert_is_numeric(ini)
  assert_is_numeric(l1)
  assert_is_numeric(l2)
  assert_is_numeric(l3)
  
  
 # if (any (is_non_positive(ini), na.rm = TRUE)) {
    # Throw an error
   # stop ("ini contains non-positive values, so no puede caminar.")
  #}
  
  if (any (is_non_positive(c(l1, l2, l3)), na.rm = TRUE)) {
    # Throw an error
    stop ("limit contains non-positive values, so no se puede calcular.")
  }
  
  
  if (l1 > 1) {
    # Throw an error
    stop ("l1 contains valor mayor a 1, so no se puede calcular.")
  }
  
  
  if (l2 > 1) {
    # Throw an error
    stop ("l2 contains valor mayor a 1, so no se puede calcular.")
  }
  
  if (l3 > 1) {
    # Throw an error
    stop ("l2 contains valor mayor a 1, so no se puede calcular.")
  }
  
  
  chr <- chr
  print (chr)
  ini <- ini
  
  x1 <- data %>%
        dplyr::filter (chrom == chr)
  
  #summary (x1)
  
  x1.na <- x1 %>%
           dplyr:::filter (R2 != "NA")
  
  
  ### secuencia para el lapply
  
  #max_delta_bp <- max(x1.na$diff.dist)
  max_delta_Mb <- max(x1.na$diff.dist)
  max_delta_Mb.r <- round(max_delta_Mb + 1,0)
  #index.Mb <- seq (1, (round(max_delta_Mb + 1,0)), 1)
  
  index.Mb <- seq (keep.Mb, max_delta_Mb.r + keep.Mb, keep.Mb)
  

  
  dt.plot.LD <- bind_rows ( lapply (index.Mb, function (filt.Mb) { 
 
  ###### verificar esto con sebas     
    x.ini.Mb <- x1.na %>%
                dplyr::filter (diff.dist <= 1 * 1e5)  ### aca me quedo con la primer 0.1 Mb siempre
    
    QQ <- quantile (x.ini.Mb$R2)
    
    l1 <- l1
    l2 <- l2
    l3 <- l3
    q1 <- quantile (QQ, l1)
    q2 <- quantile (QQ, l2)
    q3 <- quantile (QQ, l3)
    
    
    
    ###################
    XX <- data.frame( HUMk=q1 [[1]] ,  H1= q2 [[1]], HLMk =q3 [[1]], chrom = chr)
    ################### 
    
    
    print (str_c (filt.Mb, "Mb_chr_", chr))
    
    
    x2 <- x1.na %>%
          dplyr::filter (diff.dist <= filt.Mb) #este tiene que ser el argumento que cambiar
                                                          # con algun criterio segun la logitud del cromosoma
    #dplyr::mutate (Mb = "0.1")
    
    num.total <- nrow (x2) # este hay que usar despues.
    
    
    # pp <- # el summary de x2$R2
    
    p1 <- x2 %>%
          dplyr::filter (R2 > q3[[1]]) %>%
          dplyr::mutate (prob="HLMk") #%>% ### este tiene que cambia 
    
    np1 <- nrow (p1)/num.total
    
    p2 <- x2 %>%
      dplyr::filter (R2 > q2[[1]]) %>% ### este tiene que cambia 
      dplyr::filter (R2 <= q3[[1]]) %>%
      dplyr::mutate (prob= "LMk")### este tiene que cambia 
    
    np2 <- nrow (p2)/num.total
    
    p3 <- x2 %>%
      dplyr::filter (R2 >  q1[[1]]) %>% ### este tiene que cambia 
      dplyr::filter (R2 <= q2[[1]]) %>%
      dplyr::mutate (prob= "UMk") ### este tiene que cambia 
    
    
    np3 <- nrow (p3)/num.total
    
    p4 <- x2 %>%
      dplyr::filter (R2 <= q1[[1]]) %>%
      dplyr::mutate (prob= "HUMk") ### este tiene que cambia 
    
    np4 <- nrow (p4)/num.total
    
    #XX <- rbind (p1, p2,p3,p4) ## este para que esta??
    
    XX.1 <- as.numeric (rbind (np1, np2, np3, np4))
    
    XX.1.1 <- c("HLMk","LMk","UMk","HUMk")
    
    XX.2 <- data.frame (ratio = XX.1, class = XX.1.1, Mb=filt.Mb)
    
    XX.3 <- data.frame (XX.2, num.total = num.total, chrom= chr)
    
    #XX.2$clase <- factor(df2$Genotype, levels = c("Genotype 2", "Genotype 3", "Genotype 1")).
    
    # Convert the cyl variable to a factor
    XX.3$class <- as.factor(XX.3$class)
    XX.3$class <- factor(XX.3$class, levels = c("HUMk", "UMk", "LMk", "HLMk"))
    
   return (XX.3)
    
  }))
  

  df.plot.LD <- as_tibble (dt.plot.LD)
  
  df.plot.LD.W <- df.plot.LD %>%
                  dplyr::select (-num.total)
  
  
  df.plot.LD.W <- df.plot.LD %>%
    dplyr::select (-num.total)%>%
    pivot_wider(names_from = class, values_from =ratio) %>%
    dplyr::select (chrom, Mb, everything())
  
  
  write_delim (df.plot.LD.W, file=str_c("./Data/procdata/df.plot.LD.W.", chr,".txt"), 
               delim = ",", na = "NA") 
  
  
  
  df.plot.LD <- df.plot.LD %>%
                dplyr::mutate (num.cat = num.total * ratio) %>%
                dplyr::select (ratio, class, Mb, num.cat,num.total, chrom )
  
  #df.plot.LD <- df.plot.LD %>%
  #dplyr::mutate(Chrom=chrom)
  
  write_csv (df.plot.LD, file= str_c("./Data/procdata/df.freq.LD_", chr,".csv"), 
             na = "NA", append = FALSE)
  
  
  df.plot.1 <- df.plot.LD %>%
    dplyr::select (c(Mb, ratio, class))
  
  
  bplt <- ggbarplot (df.plot.1 , "Mb",  "ratio",
                     title = str_c("LD.decay by interval chr_", chr),
                     fill = "class", 
                     border ="white",
                     #sort.val = "desc",
                     x.text.angle = 90,
                     xlab = str_c ("(bin ", keep.Mb ,"Mb)"),
                     color = "class",
                     palette =c("navyblue","royalblue3","orange", "red4"),
                     label = FALSE, lab.col = "white", lab.pos = "in")
  
  bplt1 <- bplt + 
    rremove ("x.text")  
  
  bplt1 %>% 
    ggexport(filename = str_c("./Figures/plot.freq.decay.",chr,".png"))
  
  
  print (bplt1)
  
  
  #### df para el el grafico de numero de parwise
  clases <- c("HUMk", "HLMk")
  
  
  df.num.HX <-  df.plot.LD %>%
                dplyr::filter (class %in% clases) %>%
                dplyr::select (Mb, class, num.cat, chrom )
  
  df.num.HX$class <- as.character(df.num.HX$class)
  
  df.num.total <-  df.plot.LD %>%
                 #dplyr::filter (class == "HUMk") %>%
                 dplyr::select (Mb, class, num.total, chrom )%>%
                 dplyr::rename (num.cat = num.total) %>%
                 dplyr::distinct_at (vars(Mb,num.cat, chrom)) %>%
                 dplyr::mutate (class =  "Ntotal")%>%
                 dplyr::select (Mb, class, num.cat, chrom)
  
  
  
  df.num <-  bind_rows (df.num.HX, df.num.total) %>%
             dplyr::arrange (Mb)
  
  df.num$class <- as.factor (df.num$class)
  
  scatter.num <- ggscatter (df.num, x = "Mb", y = "num.cat",
                            color="class", 
                            palette = c ("red", "navyblue", "gray48"),
                            title = str_c("num.parwise by interval_", chr),
                            xlab = str_c ("(bin ", keep.Mb ,"Mb)"),
                            ylab = "num.parwise") 
  
  
  scatter.num %>%
    ggexport(filename = str_c("./Figures/plot.scatter.num",chr,".png"))
  
  df.num.HLMk <- df.num %>%
                 dplyr::filter (class == "HLMk")
  
  
  scatter.num.HLMk <- ggscatter (df.num.HLMk, x = "Mb", y = "num.cat",
                                color="red", 
                                title = str_c("num.parwise by interval_", chr),
                                xlab = str_c ("(bin ", keep.Mb ,"Mb)"),
                                ylab = "num.parwise")
  
  scatter.num.HLMk %>%
    ggexport(filename = str_c("./Figures/plot.scatter.num.HLMk",chr,".png"))
  
  print (scatter.num)
  print (scatter.num.HLMk)
  
  ## datos para el grafico de proporciones 
  
  df.plot.LD.HUMk <- df.plot.LD %>%
                     dplyr::filter (class == "HUMk")
  

  df.plot.LD.HLMk <- df.plot.LD %>%
                     dplyr::filter (class == "HLMk")

### revisar el ratio si esta bien o hay otra mejor forma  
  
  df.plot.prob.HLMk_HUMk <-  df.plot.LD.HLMk %>%
                             dplyr::mutate (ratio.HLHU = num.cat/df.plot.LD.HUMk$num.cat) %>%
                             dplyr::select (ratio.HLHU, Mb,  chrom)
                             
  df.plot.prob.HLMk_HUMk$ratio.HLHU [which(!is.finite( df.plot.prob.HLMk_HUMk$ratio.HLHU))] <- NA
  
  
  
  df.plot.LD.HLMk_HUMk <- bind_rows (df.plot.LD.HLMk, df.plot.LD.HUMk) %>%
                          dplyr::arrange (Mb)
  
  df.plot.LD.HLMk_HUMk <- df.plot.LD.HLMk_HUMk %>%
                          dplyr::select ( ratio, class, Mb,chrom )
  

  umbral.H <- df.plot.LD.W %>%
              dplyr::filter (HUMk < prob.HUMK/100 & HLMk >= prob.HLMK/100)%>%
              dplyr::filter (Mb == max(Mb))
  
  scatter.prop.total <- ggscatter (df.plot.LD.HLMk_HUMk,  x = "Mb", y = "ratio",
                                   ylim=c(0,1),
                                   title = str_c("proption by interval_", chr),
                                   color = "class",
                                   palette = c( "navyblue",  "red"),
                                   xlab = str_c ("(bin ", keep.Mb ,"Mb)")) +
    geom_hline(yintercept = umbral.H$HUMk, color = "navyblue") +
    geom_hline(yintercept = umbral.H$HLMk, color ="red") +
    geom_vline(xintercept = umbral.H$Mb, color ="black") 
  
  scatter.prop.total %>%
    ggexport(filename = str_c("./Figures/scatter.prop.total",chr,".png"))
  
  print (scatter.prop.total)
  
  
  umbral.HH <- df.plot.prob.HLMk_HUMk %>%
               dplyr::filter (Mb == umbral.H$Mb )
  

  
  scatter.propHLMk_HUMk <- ggscatter (df.plot.prob.HLMk_HUMk, x = "Mb", y = "ratio.HLHU",
                                    title = str_c("proportion HLMk/HUMk_", chr),
                                    #ylim=c(0,1),
                                    color = "black",
                                    xlab =  str_c ("(bin ", keep.Mb ,"Mb)")) +
    geom_vline(xintercept = umbral.H$Mb, color ="black") +
    geom_hline(yintercept = umbral.HH$ratio.HLHU, color ="black") 
  
  
  scatter.propHLMk_HUMk %>%
    ggexport(filename = str_c("./Figures/scatter.propHLMk_HUMk",chr,".png"))
  
  print (scatter.propHLMk_HUMk)

  
  umbral.H <- umbral.H %>%
              dplyr::mutate (HLMk_HUMk = umbral.HH$ratio.HLHU)
  

  write_csv (umbral.H, file= str_c("./Data/procdata/umbral.HLMk_HUMk_", chr,".csv"), 
             na = "NA", append = FALSE)
  
  return (df.plot.LD)
}






head(df.geno.crom)
#unique (df.geno.crom$chr)
## Argumentos




# Voy a correr todos los cromosomas 
start.time <- Sys.time()

index.chrom <- unique (df.geno.crom$chrom)

dt.all.chrom <- lapply (index.chrom, function (filt.chr) {
  
  run_freq_decay (df.geno.crom , 
                  chr=filt.chr, 
                  ini = NULL, l1= 0.25, l2=0.5, l3=0.75, keep.Mb =0.1,
                  prob.HLMK = 51, 
                  prob.HUMK = 50)
  
})

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



#names (dt.all.chrom) <- index.chrom

df.all.chrom <- do.call (rbind, dt.all.chrom)

write_delim (df.all.chrom, path = "./Data/procdata/freq_LD.txt", 
             delim = ",", na = "NA")

#### cargo los umbrales

index.chr <- 1:20

# se cargan las matrices generadas antes 
# se genear una lista con 20 tibbles
dt.unbrales.crom <- lapply (index.chr, function (filtro) {
  
  H <- read_delim (file=str_c("./Data/procdata/umbral.H_chr_", filtro,".csv"), delim = ",", 
                   na = "NA", quote = "\"",col_names = TRUE)
  #G <- G %>%
  #    dplyr::mutate (chr = str_c("chr_",filtro))
  print (H)
  return (H)
})

df.umbrales.crom <- do.call (rbind, dt.unbrales.crom)








end_time <- Sys.time()

end_time - start_time


index.chr <- 1:20

# se cargan las matrices generadas antes 
# se genear una lista con 20 tibbles
df.nt.bin_chr <- lapply (index.chr, function (filtro) { 
  
  G.1 <- read_delim (file=paste("./Data/procdata/df.nt.bin_chr_",filtro,".csv",
                                sep=""), delim = ",", 
                     na = "NA", quote = "\"",col_names = TRUE)
  G.1 <- G.1 %>%
    dplyr::mutate (chr = str_c("chr_",filtro))
  print (G.1)
  return (G.1)
})

df.nt.bin.crom <- do.call (rbind, df.nt.bin_chr)

unique (df.nt.bin.crom$chr)

chr20 <- df.nt.bin.crom %>%
  dplyr::filter (chr == "chr_20")


plot (chr20$Num.total ~ chr20$Mb, pch=16, col="red")



## filtro los HLD  ####

run_HLD_decay <- function (data) {
  
  index.Chrom <- unique(data$Chrom)
  
  dt.HLD.all.chrom <- lapply (index.Chrom, function (filtro) {
    print (filtro)
    
    dt.HLD <-  df.all.chrom %>%
      dplyr::filter (Chrom == filtro )%>%
      dplyr::filter (class == "HLD")
    
    As1 <- 0.16 
    ts1 <- 0.1
    As2 <- 0.8
    ts2 <- 10
    ys0 <- 0.01
    
    eq.HLD.decay <- nls (ratio ~ A1*exp (-Mb/t1) +  A2*exp (-Mb/t2) + y0,               # esta es la estimacion de los parametros por Gauss-Newton
                         start = list (A1 = As1,
                                       t1= ts1 ,
                                       A2 = As2,
                                       t2 = ts2,
                                       y0 = ys0 ),
                         trace = FALSE , 
                         data= dt.HLD, 
                         nls.control(warnOnly=TRUE))
    
    dt_param <- cbind (chromosome = filtro, coef (eq.HLD.decay))
    
    A1.eq <- coef (eq.HLD.decay)[1] 
    t1.eq <- coef (eq.HLD.decay)[2]    
    A2.eq <- coef (eq.HLD.decay)[3]
    t2.eq <- coef (eq.HLD.decay)[4]
    y0.eq <- coef (eq.HLD.decay)[5]
    
    pG <- t1.eq * t2.eq / (t2.eq - t1.eq)* log(A1.eq/A2.eq)
    
    dt_param <- data.frame (chromosome = filtro, 
                            A1= round (A1.eq, 3),   
                            t1= round (t1.eq, 3), 
                            A2= round (A2.eq, 3) , 
                            t2= round (t2.eq, 3),   
                            y0= round (y0.eq, 3), 
                            pinf= round (pG, 3))
  })
  
  # hay que agregar el pearson 
  # agregar el R2 del desequilibrio
  
  df.plot.LD <- as_tibble (do.call (rbind, dt.HLD.all.chrom))
  
  
  return (df.plot.LD)
  
}

run_plot_HLD_decay <- function (data) {
  
  index.Chrom <- unique(data$Chrom)
  
  dt.HLD.all.chrom <- lapply (index.Chrom, function (filtro) {
    print (filtro)
    
    dt.HLD <-  df.all.chrom %>%
      dplyr::filter (Chrom == filtro )%>%
      dplyr::filter (class == "HLD")
    
    As1 <- 0.16 
    ts1 <- 0.1
    As2 <- 0.8
    ts2 <- 10
    ys0 <- 0.01
    
    eq.HLD.decay <- nls (ratio ~ A1*exp (-Mb/t1) +  A2*exp (-Mb/t2) + y0,               # esta es la estimacion de los parametros por Gauss-Newton
                         start = list (A1 = As1,
                                       t1= ts1 ,
                                       A2 = As2,
                                       t2 = ts2,
                                       y0 = ys0 ),
                         trace = FALSE , 
                         data= dt.HLD, 
                         nls.control(warnOnly=TRUE))
    
    ## mostrar el coeficiente de correlacion de pearson
    
    dt_param <- cbind (chromosome = filtro, coef (eq.HLD.decay))
    
    A1.eq <- coef (eq.HLD.decay)[1] 
    t1.eq <- coef (eq.HLD.decay)[2]    
    A2.eq <- coef (eq.HLD.decay)[3]
    t2.eq <- coef (eq.HLD.decay)[4]
    y0.eq <- coef (eq.HLD.decay)[5]
    
    pG <- t1.eq * t2.eq / (t2.eq - t1.eq)* log(A1.eq/A2.eq)
    
    
    # el K es la cantidad de veces que aporta una
    # con respecto a la otra
    
    # pG.K <- t1.eq * t2.eq / (t2.eq - t1.eq)* log(A1.eq/(k * A2.eq))
    
    # pG.1K <- t1.eq * t2.eq / (t2.eq - t1.eq)* log(A1.eq/(1/k * A2.eq))
    
    
    py <- A1.eq*exp (-pG/t1.eq) +  y0.eq
    
    
    plot (ratio ~ Mb,
          ylim = c(0,max(dt.HLD$ratio + 0.05)),
          #xlim= c(0,max(dt.HLD$cantidad + 0.5)),
          pch=16,
          col="gray",
          cex=1,
          axes=FALSE,
          xlab="", 
          ylab="", 
          type="n",
          main = str_c("HLD_decay_", filtro),
          data=dt.HLD)
    
    axis(2, ylim=c(0,max(dt.HLD$ratio + 0.05)),col="black",las=1)  ## las=1 makes horizontal labels
    mtext("ratio.HLD",side=2,line=2.5)
    axis(1)
    mtext("delta.Mb",side=1,line=2.5)
    box()
    
    points (ratio ~ Mb,
            pch=16,
            col="gray48",
            cex=1.3,
            data=dt.HLD)
    
    curve  (A1.eq*exp (-x/t1.eq) +  A2.eq*exp (-x/t2.eq) + y0.eq, 
            from=0, to = max(dt.HLD$Mb) + 1,
            xlab="Mb", ylab="f(x)" ,lwd=2,
            add = TRUE, col="black",lty= 4)
    
    
    segments(x0=pG, y0=0, x1 = pG, y1 = py,
             col = "black", lty = 1, lwd = 2)
    
    segments(x0=-3, y0=py, x1 = pG, y1 = py,
             col = "black", lty = 1, lwd = 2)
    
    abline (h=0 , col="black", lwd=1)
    
    #abline (v=pG + 1, col="black", lwd=2)
    
    curve  (A2.eq*exp (-x/t2.eq) + y0.eq, 
            from=0, to = max(dt.HLD$Mb) + 1,
            xlab="Mb", ylab="f(x)" ,lwd=2,
            add = TRUE, col="navyblue",lty= 1)
    
    curve  (A1.eq*exp (-x/t1.eq) + y0.eq, 
            from=0, to =  max(dt.HLD$Mb),
            xlab="Mb", ylab="f(x)" ,lwd=2,
            add = TRUE, col="darkorange",lty= 1)
  })
  
}


table_HLD <- run_HLD_decay (data=df.all.chrom)

write_delim (table_HLD, path="./Data/procdata/table_HLD.txt",
             delim = ",", na = "NA")


plot_HLD <- run_plot_HLD_decay (data=df.all.chrom)

# punto medio de transicion 

end_time <- Sys.time()

end_time - start_time



plot (cantidad ~ Mb, 
      pch=16,
      col="gray",
      data= HLD.chr.3)


As1 <- 0.16 
ts1 <- 0.1
As2 <- 0.8
ts2 <- 10
ys0 <- 0.01

eq.LD.decay <- nls (cantidad ~ A1*exp (-Mb/t1) +  A2*exp (-Mb/t2) + y0,               # esta es la estimacion de los parametros por Gauss-Newton
                    start = list (A1 = As1,
                                  t1= ts1 ,
                                  A2 = As2,
                                  t2 = ts2,
                                  y0 = ys0 ),
                    trace = FALSE , 
                    data= HLD.chr.1, nls.control(warnOnly=TRUE))



# From previous step
groom_model <- function(model) {
  list(
    model = glance(model),
    coefficients = tidy(model),
    observations = augment(model)
  )
}

# Call groom_model on model, assigning to 3 variables
c(mdl, cff, obs)%<-%groom_model (model)

# See these individual variables
mdl; cff; obs
mdl [1]





dt <- HLD.chr.4
As1 <- 0.16 
ts1 <- 0.1
As2 <- 0.8
ts2 <- 10
ys0 <- 0.01

eq.LD.decay <- nls (cantidad ~ A1*exp (-Mb/t1) +  A2*exp (-Mb/t2) + y0,               # esta es la estimacion de los parametros por Gauss-Newton
                    start = list (A1 = As1,
                                  t1= ts1 ,
                                  A2 = As2,
                                  t2 = ts2,
                                  y0 = ys0 ),
                    trace = FALSE , 
                    data= dt, nls.control(warnOnly=TRUE))

A1.eq <- coef (eq.LD.decay)[1]      # este es el A1 de Gauss-Newton
t1.eq <- coef (eq.LD.decay)[2]    
A2.eq <- coef (eq.LD.decay)[3]
t2.eq <- coef (eq.LD.decay)[4]
y0.eq <- coef (eq.LD.decay)[5]


pG <- t1.eq * t2.eq / (t2.eq - t1.eq)* log(A1.eq/A2.eq)

pGG1 <- A1.eq*exp (-(pG+x) /t1.eq) 
pGG2 <-  A2.eq*exp (-pG /t2.eq) 




plot (cantidad ~ Mb, 
      pch=16,
      cex=1,
      col="gray48",
      data= dt)


curve  (A1.eq*exp (-x/t1.eq) +  A2.eq*exp (-x/t2.eq) + y0.eq, 
        from=0, to = max(dt$Mb) + 1,
        xlab="Mb", ylab="f(x)" ,lwd=2,
        add = TRUE, col="blue",lty= 1)

abline (v=pG + 0.1, col="red", lwd=2)

abline (v=pG + 1, col="black", lwd=2)

curve  (A2.eq*exp (-x/t2.eq) + y0.eq, 
        from=0, to = max(dt$Mb) + 1,
        xlab="Mb", ylab="f(x)" ,lwd=2,
        add = TRUE, col="darkgreen",lty= 1)

curve  (A1.eq*exp (-x/t1.eq) + y0.eq, 
        from=0, to = pG,
        xlab="Mb", ylab="f(x)" ,lwd=2,
        add = TRUE, col="darkorange",lty= 1)

plot (cantidad ~ Mb, 
      xlim = c(0,1.5),
      pch=16,
      cex=1.5,
      col="gray48",
      data= dt)


curve  (A1.eq*exp (-(pG+x) /t1.eq), 
        from=0, to = 1,
        xlab="Mb", ylab="f(x)" ,lwd=2,
        add = TRUE, col="black",lty= 1)

curve  (A2.eq*exp (-(pG+x) /t2.eq) , 
        from=0, to = 1,
        xlab="Mb", ylab="f(x)" ,lwd=2,
        add = TRUE, col="red",lty= 1)






curve  (B.eq.07 *(1- exp (- k.eq.07  * x)), 
        from=0, to = max(dt$time) + 1, xlab="time", ylab="ET" ,lwd=1,
        add = TRUE, col="red",lty= 2)









###################################### ACA TERMINA EL CODIGO 33 ################################
Fun.sebin <- function (data =NULL, r2=0, b.i=0, step=1, uni ="Mb"){
  
  #list.crom <- index.chr
  data <- as.data.frame(data)
  list.crom <- str_c("chr_",index.chr)
  unique (data$chr) == list.crom
  r2 <- r2
  b.i <- b.i
  step <- step
  uni <- uni
  #head(data)
  
  if (uni == "Mb") {
    data <- data %>%
      dplyr::mutate (delta.Mb = delta.bp/1e6) %>%
      dplyr::filter (R2 > r2)
  }
  
  if (uni == "Kb") {
    data <- data %>%
      dplyr::mutate (delta.Kb = delta.bp/1e3)%>%
      dplyr::filter (R2 > r2)
  }
  
  dt.bin.1 <- lapply (list.crom, function (filtro){
    
    dt.1 <- data %>%
      dplyr::filter (chr == filtro)
    
    if (uni == "Mb") {
      xbin <- dt.1 %>%
        dplyr::select(delta.Mb)
      
    }
    
    if (uni == "Kb") {
      xbin <- dt.1 %>%
        dplyr::select(delta.Kb)
    }
    
    #############
    xf <- (round (max(xbin),0)-1)
    
    lisbin <- (seq (b.i ,xf, step))
    
    dt.bin.2 <- lapply (lisbin, function (x) { 
      
      #print(filtro)
      if (uni == "Mb"){
        if ( x == b.i ) {
          dt.2 <-  dt.1  %>%
            dplyr::filter ( delta.Mb > x  )%>%
            dplyr::filter ( delta.Mb <= x + step )%>%
            dplyr::mutate (bin = str_c(b.i + step,uni))
          #print(dt.2)
          return (dt.2)
        }
        
        if ( x > b.i ) { 
          if ( x < xf ) {
            dt.2  <-  dt.1  %>%
              dplyr::filter ( delta.Mb > x  )%>%
              dplyr::filter ( delta.Mb <= x + step )%>%
              dplyr::mutate (bin = str_c(x + step,uni))
            #print(dt.2 ) 
            return (dt.2 )
          }
        } 
      }
      
      
      if (uni == "Kb"){
        if ( x == b.i ) {
          dt.2 <-  dt.1  %>%
            dplyr::filter ( delta.Kb > x  )%>%
            dplyr::filter ( delta.Kb <= x + step )%>%
            dplyr::mutate (bin = str_c(b.i +  step, uni))
          return (dt.2)
        }
        
        if ( x > b.i ) { 
          if ( x < xf ) {
            dt.2  <-  dt.1  %>%
              dplyr::filter (delta.Kb > x  )%>%
              dplyr::filter (delta.Kb <= x + step )%>%
              dplyr::mutate (bin = str_c(x + step,uni))
            return (dt.2 )
          }
        } 
      }
      
      
    })
  })
  
  #return (dt.bin.1)
  # names (dt.bin.1) <-  list.crom
  #print (dt.bin.1)
  #return (data)
} ### aca termina sebin


start.time <- Sys.time()
X <- Fun.sebin (data =df.geno.crom, r2=0, b.i=0, step=1, uni ="Mb")
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


list.crom <- str_c("chr_",index.chr)

lapply (index.chr , function (x){
  plotX <- X [[x]]
  plotX1 <- do.call (rbind, plotX)
  
  plotX2 <- plotX1%>%
    dplyr::arrange (delta.bp)
  
  ggboxplot (plotX2, x="bin", y = "R2", xlab = str_c("bin", "(Mb)") ,
             title= str_c ("chr_", index.chr), 
             font.label = list(size =1, color = "black"),
             width = 0.8)
  
  #print (plotX2)
})








# cargo la matrices con el R2
#
# esto lo tenes setear
index.chr <- 1:20


soja.crom.r2 <- lapply (index.chr, function (filt.chr) { 

  Gr2 <- read_delim (file=str_c("./Data/procdata/LD.soja.cross.matrix.",filt.chr,".txt"),
                               delim = ",",  na = "NA", quote = "\"",col_names = TRUE)
  Gr2 <- Gr2 %>%
    dplyr::mutate (chr = str_c("chr_",filt.chr))
  print (Gr2)
  return (Gr2)
})



#### 

### 
soja.cross.tibble <- read_delim (file="./Data/procdata/cross.data.z14_15.csv",
                                         delim=",", quote = "\"", 
                                         na = "-", col_names = TRUE)


soja.cross.chr1.tibble_snp <- soja.cross.tibble %>%
                        dplyr::select (starts_with("S01"))


#soja.cross.chr1.tibble <- soja.cross.chr1.tibble %>%
 #                          dplyr::select(-S)


#soja.cross.chr1.tibble_snp  <- soja.cross.chr1.tibble %>%
 #                                 dplyr::select (starts_with("S"))



soja.cross.chr1.tibble_snp.1 <- soja.cross.chr1.tibble_snp [-c(1,2),]


soja.cross.chr1.tibble_bp <- soja.cross.chr1.tibble_snp [2,]

mkr.chr1 <- colnames (soja.cross.chr1.tibble_snp)

x.mk <- data.frame (mrks = mkr.chr1)

pos.chr1    <- soja.cross.chr1.tibble_snp [2,]
pos.chr1.1  <- data.frame (pos = pos.chr1[1,])

pos.chr1.2  <- as.numeric(t(pos.chr1.1))

soja.cross.chr1.tibble_bp.1 <- x.mk %>%
                              dplyr::mutate (pos=pos.chr1.2)

### 
dt <- soja.cross.chr1.tibble_bp.1 
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
    
    colnames (LD.soja.cross.chr1.matrix) <- list.mrks
    rownames (LD.soja.cross.chr1.matrix) <- list.mrks
    LD_soja.cross.chr1 <-  as_tibble(LD.soja.cross.chr1.matrix)
    LD_soja.cross.chr1 <- LD_soja.cross.chr1 %>%
      dplyr::mutate (mrk.id =list.mrks)%>%
      dplyr::select (mrk.id, everything())
    
    id <- colnames (XX.x.bp.1)
    
    LD_soja.cross.chr1.1 <- LD_soja.cross.chr1%>%
      dplyr::filter (mrk.id==id)
    
    id.gather <- colnames (LD_soja.cross.chr1.1) [-1]
    
    LD_soja.cross.chr1.2 <- LD_soja.cross.chr1.1 %>%
                            gather (id.gather, key="mrk.2" , value = "R2") %>%
                            dplyr::rename (mrk.1 = mrk.id) 
    
    LD_soja.cross.chr1.2$mrk.1 <- as.character(LD_soja.cross.chr1.2$mrk.1)
    
    
    colnames (XX.x.bp.1) == unique ( LD_soja.cross.chr1.2$mrk.1)
    
    LD_soja.cross.chr1.3 <- LD_soja.cross.chr1.2 %>%
                            dplyr::mutate (delta.bp = delta.bp )
   
     print (LD_soja.cross.chr1.3)
    return ( LD_soja.cross.chr1.3)

})

plot.LD.decay.chr1 <- as.data.frame (do.call (rbind, dt.LD.decay))
unique(plot.LD.decay.chr1$mrk.1)


write_delim (plot.LD.decay.chr1 , path="./Data/procdata/plot.LD.decay.chr1.R2.txt",
             delim = ",", na = "NA")


head(plot.LD.decay.chr1)
plot.LD.decay.chr1a <- plot.LD.decay.chr1 %>%
                       dplyr::mutate (chr=1)


######


#### esta es la funcion anterior ############################

#dt.x.bp <- lapply (list.mrks, function (filtro) {
 # xx <- x.bp [[filtro]] 
  #dt.xx <- as.data.frame (do.call (rbind, xx ))
  #colnames (dt.xx) <- filtro
  #return (dt.xx)
  
#})

#datosx.bp <- as.data.frame (do.call (cbind, dt.x.bp))

#datosx.bp [[1]]

#View(datosx.bp)

#datosx.bp.1 <- datosx.bp %>%
 #              dplyr::mutate (mrk = list.mrks)%>%
  #             dplyr::select(mrk, everything())

#datosx.bp.1 [, 1:2]

#View (datosx.bp.1)

#MyHeatmap.Dprim <- LDheatmap (CEUSNP, CEUDist, LDmeasure="D'",
#title="Pairwise LD in D'", add.map=TRUE,
# SNP.name=c("rs2283092", "rs6979287"),
#color=grey.colors(20), name="myLDgrob",
#add.key=TRUE)


#Q <- LD.soja.cross.chr1.matrix 

#Q [is.na(Q)] <- 0

#QQ <- t (Q)

#Q.Dprim <- Q + QQ

#datos.Dprim <- as.data.frame (Q.Dprim)


#datos.Dprim <- datos.Dprim %>%
 # dplyr::mutate (mrk = rownames (datos.Dprim))%>%
  #dplyr::select (mrk, everything())

#View (datos.Dprim)

#sapply(list, function)

#dt.xy <- lapply (list.mrks, function (filtro) {
 # px <- datosx.bp.1 %>%
  #      dplyr::select (filtro)
  #py <- datos.Dprim %>%
   #     dplyr::select (filtro)
  #pxy <- cbind ( px,  py)
  #return (pxy)
#})


#m1 <- dt.xy [[1]]
#names(m1) [1] <- "dist.bp"
#names(m1) [2]  <- "dist.cor" 
#m1.1 <- m1 %>%
 #       filter (dist.bp > 0)

#plot (m1.1$dist.bp, m1.1$dist.cor)

#ggscatter(m1.1, x = "dist.bp", y = "dist.cor",
          #add = "loess",
          #conf.int = TRUE)

#m1.1 <- m1.1 %>%
 #       dplyr::mutate (chr=1) 

#write_delim (m1.1 , path="./Data/procdata/LD.decay.chr1.txt",
 #            delim = ",", na = "NA")


