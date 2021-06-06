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

######### 
# primero cargar la informacion que se tenga del mapa.

mapa.cebada.bill  <- read_delim (file = "./Data/rawdata/8.Bill_thomas_table_3.txt" ,
                      col_names = TRUE, delim = "\t", na = "NA")

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

mapa.EEMAC.cross.bill.bp <- mapa.EEMAC.cross.bill  %>%
                            dplyr::filter (!is.na(bp))

pos.unique <- unique (mapa.EEMAC.cross.bill.bp$pos)

mapa.EEMAC.cross.cM <- mapa.EEMAC.cross %>%
                       dplyr::filter (pos %in% pos.unique)


#### hago una lista de los marcadores que debo sacar del mapa de la poblacion de mapeo

head ( mapa.EEMAC.cross.bill.bp)

list.mark <- unique (mapa.EEMAC.cross.bill.bp$snp)

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


mapa.EEMAC.cross.bill.bp.select <- mapa.EEMAC.cross.bill.bp %>%
                                   dplyr::select (c ("chr", "snp", "bp"))

head (map.cross.2.bp  )

map.cross.2.bp <- map.cross.2 %>%
                  dplyr::inner_join (mapa.EEMAC.cross.bill.bp.select, by="snp") %>%
                  dplyr::select (c ("snp","chr.x", "bp")) %>%
                  dplyr::rename (chr = chr.x) %>%
                  dplyr::rename (pos = bp)

map.cross.2.bp  <- write_delim (map.cross.2.bp, 
                               file ="./Data/procdata/2.cross_EEMAC.1_map.bp.csv", delim = ",")

#### Antes de cargar de nuevo borrar la etiqueta de la columna "snp"

EEMAC.cross.1.bp <- read.cross (format="tidy",
                               dir="./Data/procdata",
                               genfile = "2.cross_EEMAC.1_gen.csv",
                               mapfile = "2.cross_EEMAC.1_map.bp.csv", 
                               phefile =  "2.cross_EEMAC.1_phe.csv",
                               na.strings="-", 
                               genotypes=c("AA","BB"),
                               crosstype ="dh",
                               estimate.map=FALSE, error.prob=0.0001)
  
summary (EEMAC.cross.1.bp)





plotMap (EEMAC.cross.1.bp, horizontal=FALSE, show.marker.names=FALSE)

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



data.cross = EEMAC.cross.1.bp
heterozygotes = FALSE
distance.unit = "bp"
id.cross = "EEMAC.cross.1.bp"



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

if   ( distance.unit != "cM" &  distance.unit != "bp" &  distance.unit != "kb" & distance.unit != "Mb") {
  
  stop (str_c ("Debe definir una unidad de distancia valida"))
  
  
}

if (is.null (distance.unit)) {
  
  stop (str_c ("Falta definir unidad de distancia"))
  
}   

if   (distance.unit == "bp" | distance.unit == "Kb" | distance.unit == "Mb"  ) {
  
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
# Nota: las figuras se podrian guardar solas como png.
 
   if   ( distance.unit == "bp" | distance.unit == "Kb" | distance.unit == "Mb" ) {
    
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



dt.diff.dist <- bind_rows (lapply (list.pos, function (filtro.x1) { 
  
   #filtro.x1 =  277230
   #print (filtro.x1)
  
 dt.dist.2 <- bind_rows ( lapply (list.pos, function (filtro.x2) {
    
    #filtro.x2 =  2415604 
    #print (filtro.x2)
    
    dt.x1 <- dt %>% 
             dplyr::filter (pos == filtro.x1)
  
    dt.x2 <- dt %>% 
             dplyr::filter (pos == filtro.x2)
    
    dt.z <- data.frame (dt.x1,  dt.x2) 
    
    dt.z <- dt.z  %>%
            dplyr::mutate (diff.dist = abs (pos - pos.1)) %>%
            dplyr::select ("mrks", "mrks.1", "diff.dist" ) %>%
            dplyr::rename (mrk.1 = mrks) %>% 
            dplyr::rename (mrk.2 = mrks.1) 
    
  #print (dt.z)
    
  #return (dt.z)
      
  }))
  
}))

start.time <- Sys.time()

df.LD.decay <- bind_rows ( lapply (list.mrks, function (filt.mrk) {
  
  #filt.mrk= "JHI-Hv50k-2016-270"

  #print (filt.mrk)
  
  dt.diff.dist.mrk <-  dt.diff.dist %>% 
                       dplyr::filter (mrk.1 == filt.mrk) 
  
  colnames (LD.cross.matrix ) <- list.mrks
  rownames (LD.cross.matrix)  <- list.mrks
 
  LD.cross <-  as_tibble (LD.cross.matrix)
  LD.cross <- LD.cross %>%
                   dplyr::mutate (mrk.id = list.mrks)%>%
                   dplyr::select (mrk.id, everything())
  
  LD.mrk.1 <- LD.cross %>%
              dplyr::filter (mrk.id== filt.mrk)
  
  id.gather <- colnames (LD.mrk.1) [-1]
  
  LD.cross.2 <- LD.mrk.1 %>%
                gather (all_of (id.gather), key="mrk.2" , value = "R2") %>%
                dplyr::rename (mrk.1 = mrk.id) 
  
  LD.dist.cross.3 <- LD.cross.2 %>%
                    dplyr::inner_join ( dt.diff.dist.mrk ,  LD.cross.2, by= c("mrk.1", "mrk.2"))
  
  #return (LD.dist.cross.3)
  
}))




df.LD.decay <- df.LD.decay  %>%
               dplyr::mutate (LG = str_c ("lg.", filt.chr)) %>%
               dplyr::arrange (diff.dist) 

### verificar estas medidas

if   ( distance.unit == "bp" ) {
  
  df.LD.decay <- df.LD.decay  %>%
    dplyr::mutate (diff.dist = (diff.dist*1e-6))
  
}


if   ( distance.unit == "Kb" ) {
  
  df.LD.decay <- df.LD.decay  %>%
                 dplyr::mutate (diff.dist = diff.dist/1e6)

}

write_delim (df.LD.decay  , file =str_c("./Data/procdata/", id.cross, "LD.decay.LG", filt.chr,".txt"),
             delim = ",", na = "NA")

if   ( distance.unit == "bp" | distance.unit == "Kb" | distance.unit == "bp" ) {

  png (filename = str_c("./Figures/",id.cross,".LD.decay.LG_", filt.chr,".png"),
    width = 480, height = 480, units = "px", pointsize = 12,
    bg = "white", res = NA)

plot (x = df.LD.decay$diff.dist , y = df.LD.decay$R2, 
      main=str_c(id.cross,".LD.decay.LG_", filt.chr),
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (df.LD.decay$diff.dist)), 
      ylim = c(0, max (df.LD.decay$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 

axis(side = 2, las = 1)
x2 <- max (df.LD.decay$diff.dist, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (df.LD.decay$diff.dist,df.LD.decay$R2, 
        pch = 20, cex=1.5, col="gray28") 
box()
dev.off()

plot (x = df.LD.decay$diff.dist , y = df.LD.decay$R2, 
      main=str_c(id.cross,".LD.decay.LG_", filt.chr),
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (df.LD.decay$diff.dist)), 
      ylim = c(0, max (df.LD.decay$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 

axis(side = 2, las = 1)
x2 <- max (df.LD.decay$diff.dist, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)

points (df.LD.decay$diff.dist, df.LD.decay$R2, 
        pch = 20, cex=1.5, col="gray28") 
box()
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


# Time difference of 1.000327 mins
start.time <- Sys.time()

LD.EEMAC.cross.1.bp <- run_plot_heatmap_LD (id.cross ="EEMAC.cross.1.bp" , 
                                            data.cross = EEMAC.cross.1.bp,  
                                            heterozygotes = FALSE,
                                            distance.unit = "bp")
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

max (df.geno.crom.2$diff.dist)

#########
# Argumentos anteriores

id.cross = "EEMAC.cross.1"
data= df.geno.crom
ini = 1
l1= 0.25
l2=0.5
l3=0.75
seq1 = c(1, 10, 50,100,)
distance.unit = "Mb"


run_table_quantiles <- function (id.cross=NULL,data= NULL, ini = 1, l1= 0.25, l2=0.5, l3=0.75, seq1= NULL, distance.unit = NULL) {
  
  if   ( distance.unit != "cM" &  distance.unit != "Mb" & distance.unit != "Kb" & distance.unit != "bp" ) {
    
    stop (str_c ("Debe definir una unidad de distancia valida"))
    
    
  }
  
  if (is.null (distance.unit)) {
    
    stop (str_c ("Falta definir unidad de distancia"))
    
  }
  
  if   ( distance.unit == "Mb" | distance.unit == "Kb" | distance.unit == "bp") {
    
    print (str_c ("El mapa es un mapa fisico"))
    
  }
  
  if   ( distance.unit == "cM") {
    print (str_c ("El mapa es un mapa genetico"))
    
  }
  
  
  list.chrom <- unique (data$chrom)
  filt.crom =2  ###
  
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
     # if   ( distance.unit == "Mb" ) {
        
        x.dist.Mb <- dat.1 %>%
                     dplyr::filter (diff.dist <= filt.dist * 1e5)
        
        bin <- (filt.dist* 1e5)/1e6
        
        x.dist.Mb.1 <- x.dist.Mb %>%
                       dplyr::mutate (bin = str_c (bin, "Mb"))
        
        print (x.dist.Mb.1)
        
        
        
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
    
    
    ggh <- gghistogram (df.hist.plot, y = "..density..", x = "R2", facet.by = "bin",title =filt.crom,
                       #fill = "lightgray",
                       fill = "bin", palette =   "RdBu",
                       add = "median", rug = TRUE)
    print (ggh)
    
    ggh  %>%
      ggexport(filename = str_c("./Figures/ggh_",id.cross, "_", filt.crom,".png"))
    
    
    qqunif <- ggplot(df.hist.plot, aes(R2, colour =  bin , fill = bin), size=1.2) + stat_ecdf() +
      theme_bw()+
      stat_function(fun=punif,args=list(0,1))+
      scale_color_manual(values=c("green", "black", "blue", "red"))+
      labs(title=str_c("ECDF and theoretical CDF","_", filt.crom)) +
      labs(y = "Theoretical Quantiles", x = "Sample Quantiles")
    
    print (qqunif)
    
    qqunif %>%
      ggexport(filename = str_c("./Figures/qqunif_",id.cross, "_", filt.crom,".png"))
    
    if   ( distance.unit == "Kb" ) {
      x.ini.Mb <- dat.1 %>%
        dplyr::filter (diff.dist <= ini * 1e5)
      
      QQ <- quantile (x.ini.Mb$R2)
      
      l1 <- l1
      l2 <- l2
      l3 <- l3
      q1 <- quantile (QQ, l1)
      q2 <- quantile (QQ, l2)
      q3 <- quantile (QQ, l3)
      
      XX <- data.frame( HLE=round (q1 [[1]],2) ,  H1= round (q2 [[1]],2), HLD = round (q3 [[1]],2), chrom = filt.crom)
    }
    
    
    
    if   ( distance.unit == "bp" ) {
      x.ini.Mb <- dat.1 %>%
                  dplyr::filter (diff.dist <= ini * 1e-6)
      
      QQ <- quantile (x.ini.Mb$R2)
      
      l1 <- l1
      l2 <- l2
      l3 <- l3
      q1 <- quantile (QQ, l1)
      q2 <- quantile (QQ, l2)
      q3 <- quantile (QQ, l3)
      
      XX <- data.frame( HLE=round (q1 [[1]],2) ,  H1= round (q2 [[1]],2), HLD = round (q3 [[1]],2), chrom = filt.crom)
    }
    
    
    
    if   ( distance.unit == "cM" ) {
      
      # Nota: aca los estoy estoy tomando cada 0.5 cM hay que ver si esto reemplaza a 0.1 Mb
      
      list.seqCM <- seq (from = 0.5, to = 20, by = 0.5)
      #filt.distCM <- 0.5
      df.QQ.seq.cM <- bind_rows (lapply (list.seqCM, function (filt.distCM){
        
        x.ini.cM <- dat.1 %>%
          dplyr::filter (diff.dist <= ini * filt.distCM)
        
        
        
        QQ <- quantile (x.ini.cM$R2)
        
        l1 <- l1
        l2 <- l2
        l3 <- l3
        q1 <- quantile (QQ, l1)
        q2 <- quantile (QQ, l2)
        q3 <- quantile (QQ, l3)
        
        XX <- data.frame( HLE=round (q1 [[1]],2) ,  H1= round (q2 [[1]],2), HLD = round (q3 [[1]],2), chrom = filt.crom)
        
        dt.QQ.seq.cM <- XX %>%
          dplyr::mutate (inter.cM = filt.distCM ) %>%
          dplyr::select (chrom, inter.cM, HLE,   H1,  HLD )
        
      }))
      
      
      df.QQ.seq.cM.long <- df.QQ.seq.cM %>%
        pivot_longer( ! c(chrom,inter.cM), names_to = "cat.LD", values_to = "unqq.R2")
      
      
      qq.scatt <- ggscatter (df.QQ.seq.cM.long, x = "inter.cM", y = "unqq.R2",
                             title = str_c("Caida de qqR2 en funcion del bin", id.cross, "_", filt.crom),
                             color = "cat.LD", shape = "cat.LD",
                             palette = c("navyblue", "gray48", "darkorange"),
                             add = "loess", rug= TRUE)
      
      print (qq.scatt)
      
      qq.scatt  %>%
        ggexport(filename = str_c("./Figures/qq.scatt_",id.cross, "_", filt.crom,".png"))
      
    }
    
  }))
  return (df.qq)
}


start.time <- Sys.time()

EEMAC.cross.1.qq <- run_table_quantiles  ( id.cross = "EEMAC.cross.1",
                       data= df.geno.crom,
                       ini = 1,
                       l1= 0.25,
                       l2=0.5,
                       l3=0.75,
                       seq1 = c(1, 5, 10, 50, 100, 500),
                       distance.unit = "bp")

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken



















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


