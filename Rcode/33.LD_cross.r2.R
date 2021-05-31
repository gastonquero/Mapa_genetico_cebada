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

data.cross = soja.cross 
heterozygotes = FALSE
data.tibble = soja.cross.tibble
distance.unit = "Mb"
id.cross = "GWAS.soja"

############# CEBADA 
## Si el mapa esta en cM se debe correr la funcion jittermap para 
# evitar la misma posicion de los marcadores.

EEMAC.cross <- read.cross (format="csv",
                           dir="./Data/rawdata", file= "1.cross_EEMAC.1.csv",
                           na.strings="NA", 
                           genotypes=c("0","1"),
                           #alleles=c("0","1"),
                           estimate.map=FALSE, 
                           convertXdata=TRUE, error.prob=0.0001)

EEMAC.cross <- jittermap (EEMAC.cross, amount=1e-6)

#data.cross = EEMAC.cross
#heterozygotes = FALSE
#distance.unit = "cM"
#id.cross = "EEMAC.1"
#distance.unit = "cM"

#filt.chr = 1
#data.tibble [1:2, 3610:3618]

## esta funcion calcula el r2 y grafica el heatmap.

run_plot_heatmap_LD <- function (id.cross = NULL , 
                                 data.cross = NULL,  
                                 heterozygotes = FALSE,
                                 distance.unit = NULL) {

# Nota: se debe crear la estructura del directorio

print (str_c("Se encontraron " , nchr (data.cross), " grupos de ligamientos")) # verifico el numero de cromosomas

if   ( distance.unit != "cM" &  distance.unit != "Mb") {
  
  stop (str_c ("Debe definir una unidad de distancia valida"))
  
  
}

if (is.null (distance.unit)) {
  
  stop (str_c ("Falta definir unidad de distancia"))
  
}   

if   ( distance.unit == "Mb" | distance.unit == "Kb" ) {
  
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
#  Nota: ver si con los nuevos formaatos cambia o no
  
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
 
   if   ( distance.unit == "Mb" | distance.unit == "Kb" ) {
    
    plot.hm <- LDheatmap ( ld$"R^2", 
                           genetic.distances= map.cross$pos, 
                           distances="physical",
                           title=str_c(id.cross, "Pairwise LD_LG.",filt.chr),
                           color = colorRampPalette(c("red4", "red","orangered", "orange","yellow1",  "blue4"))(60))
  }
   
   if   ( distance.unit == "cM") {
     
     plot.hm <- LDheatmap ( ld$"R^2", genetic.distances= map.cross$pos, 
                            distances="genetic",
                            title=str_c(id.cross, "Pairwise LD_LG.",filt.chr),
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
  
   #filtro.x1 = 153541 
   #print (filtro.x1)
  
 dt.dist.2 <- bind_rows ( lapply (list.pos, function (filtro.x2) {
    
    #filtro.x2 = 293684
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
  
df.LD.decay <- bind_rows ( lapply (list.mrks, function (filt.mrk) {
  
  #filtro = "S01_153541"

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

if   ( distance.unit == "Kb" ) {
  
  df.LD.decay <- df.LD.decay  %>%
                 dplyr::mutate (diff.dist = diff.dist/1e6)

}

write_delim (df.LD.decay  , file =str_c("./Data/procdata/", id.cross, "LD.decay.LG", filt.chr,".txt"),
             delim = ",", na = "NA")

if   ( distance.unit == "Mb" | distance.unit == "Kb" ) {

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

if   ( distance.unit == "cM" ) {
  
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
        xlab = expression(Distance ~ (cM))) 
  
  axis(side = 2, las = 1)
  x2 <- max (df.LD.decay$diff.dist, na.rm = TRUE)
  axis (side=1,at=seq(0,x2,1),las = 1)
  
  
  points (df.LD.decay$diff.dist, df.LD.decay$R2, 
          pch = 20, cex=1.5, col="gray28") 
  box()
  dev.off()
  
  plot (x = df.LD.decay$diff.dist, y = df.LD.decay$R2, 
        main=str_c(id.cross,".LD.decay.LG_", filt.chr),
        pch = 20, 
        type ="n",
        xaxt="none",
        yaxt="none",
        axes = F,
        xlim = c(0, max (df.LD.decay$diff.dist)), 
        ylim = c(0, max (df.LD.decay$R2, na.rm = TRUE)),
        ylab = expression(LD ~ (r^2)),
        xlab = expression(Distance ~ (cM))) 
  
  axis(side = 2, las = 1)
  x2 <- max (df.LD.decay$diff.dist, na.rm = TRUE)
  axis (side=1,at=seq(0,x2,1),las = 1)
  
  
  points (df.LD.decay$diff.dist,df.LD.decay$R2, 
          pch = 20, cex=1.5, col="gray28") 
  box()
}

})

}






start.time <- Sys.time()
X <- run_plot_heatmap_LD ( id.cross = "EEMAC.cross.1" , 
                           data.cross = EEMAC.cross,  
                           heterozygotes = FALSE,
                           distance.unit = "cM" )
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken









# Time difference of 1.000327 mins

run_plot_heatmap_LD (data.cross  = soja.cross, 
                     data.tibble = soja.cross.tibble,  
                     heterozygotes = FALSE)



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


