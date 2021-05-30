##########################################################
# Codigo para el LD y LD decay de la ploblacion de GWAS_Soja
#
# Gaston Quero - Sebastian Simomndi 
# 3/3/2020
###############################################
getwd ()
setwd ("C:/Users/Usuario/OneDrive/Documentos/Paper_Crop_Science")


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


## reingreso la poblacion solo de cada cromosoma
soja.cross <- read.cross(format="csv",
                                 dir="./Data/procdata", file= "cross.data.z14_15.csv",
                                 na.strings="-", 
                                 genotypes=c("AA","BB"),
                                 #alleles=c("0","1"),
                                 estimate.map=FALSE,
                                 crosstype ='dh',
                                 convertXdata=TRUE, 
                                 error.prob=0.0001)

soja.cross.tibble <- read_delim (file="./Data/procdata/cross.data.z14_15.csv",
                                 delim=",", quote = "\"", 
                                 na = "-", col_names = TRUE)

plotMap (soja.cross, horizontal=TRUE, show.marker.names=FALSE)

data.cross = soja.cross 
heterozygotes = FALSE
data.tibble = soja.cross.tibble

filt.chr = 1
data.tibble [1:2, 3610:3618]
## esta funcion calcula el r2 y grafica el heatmap.
run_plot_heatmap_LD <- function(data.cross = NULL, data.tibble =NULL,  heterozygotes = FALSE){
  
nchr.max <- nchr (data.cross)
  
list.chr <- 1:nchr.max
filt.chr = 1


lapply(list.chr, function (filt.chr){
  
  print (filt.chr)
  crossobj <- subset(data.cross, chr=filt.chr)
  chr=filt.chr
  
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
  ld <- LD(genos)
  
  # plot LD heatmap
  plot.hm <- LDheatmap(ld$"R^2", title=str_c("Pairwise LD_Chr.",filt.chr),
                       color = colorRampPalette(c("red4", "red","orangered", "orange","yellow1",  "blue4"))(60))


  LD.soja.cross.matrix <- plot.hm$LDmatrix

  
  write.table (LD.soja.cross.matrix , 
               file = str_c("./Data/procdata/LD.soja.cross.matrix.",filt.chr,".txt"),
               append = FALSE, quote = TRUE, sep = ",",
               eol = "\n", na = "NA", dec = ".", row.names = TRUE,
               col.names = TRUE)
  
  if(filt.chr < 10 ) {
   
     soja.cross.tibble_snp <- data.tibble  %>%
                              dplyr::select (starts_with(str_c("S0", filt.chr)))
  }
  
  if(filt.chr >= 10 ) {
    
soja.cross.tibble_snp <- data.tibble  %>%
      dplyr::select (starts_with(str_c("S", filt.chr)))
 
     }

#soja.cross.chr.tibble_snp.1 <- soja.cross.tibble_snp [-c(1,2),]
  
  
  soja.cross.tibble_bp <- soja.cross.tibble_snp [2,]
  
  mkr.chr <- colnames (soja.cross.tibble_snp)
  
  x.mk <- data.frame (mrks = mkr.chr)
  
  pos.chr    <- soja.cross.tibble_snp [2,]
  pos.chr.1  <- data.frame (pos = pos.chr[1,])
  
pos.chr.2  <- as.numeric(t(pos.chr.1))

soja.cross.tibble_bp.1 <- x.mk %>%
                              dplyr::mutate (pos=pos.chr.2)

### 
dt <- soja.cross.tibble_bp.1 
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
  
names(x.bp) <- list.mrks
XX.x.bp <- as.data.frame (do.call (cbind, x.bp))
  
dt.LD.decay <- lapply (list.mrks, function (filtro) {
  
  XX.x.bp.1 <- XX.x.bp %>% 
               dplyr::select (filtro) 
  delta.bp <- unlist (XX.x.bp.1 , use.names=FALSE)
  
  colnames (LD.soja.cross.matrix ) <- list.mrks
  rownames (LD.soja.cross.matrix) <- list.mrks
  LD_soja.cross <-  as_tibble(LD.soja.cross.matrix)
  LD_soja.cross <- LD_soja.cross %>%
                   dplyr::mutate (mrk.id =list.mrks)%>%
                   dplyr::select (mrk.id, everything())
  
  id <- colnames (XX.x.bp.1)
  
  LD_soja.cross.1 <- LD_soja.cross%>%
                    dplyr::filter (mrk.id==id)
  
  id.gather <- colnames (LD_soja.cross.1) [-1]
  
  LD_soja.cross.2 <- LD_soja.cross.1 %>%
                     gather (id.gather, key="mrk.2" , value = "R2") %>%
                     dplyr::rename (mrk.1 = mrk.id) 
  
  LD_soja.cross.2$mrk.1 <- as.character(LD_soja.cross.2$mrk.1)
  
  
  colnames (XX.x.bp.1) == unique ( LD_soja.cross.2$mrk.1)
  
  LD_soja.cross.3 <- LD_soja.cross.2 %>%
                     dplyr::mutate (delta.bp = delta.bp )
  
  print (LD_soja.cross.3)
  
  

  return ( LD_soja.cross.3)
  
})

df.LD.decay <- as.data.frame (do.call (rbind, dt.LD.decay))

df.LD.decay <- df.LD.decay  %>%
               dplyr::mutate (chrom = str_c ("chr_", filt.chr))%>%
               dplyr::arrange(delta.bp)%>%
               dplyr::mutate (delta.bp.1 = delta.bp/1e6)

write_delim (df.LD.decay  , path=str_c("./Data/procdata/plot.LD.decay.chr", filt.chr,".txt"),
             delim = ",", na = "NA")

png (filename = str_c("./Figures/plot.r2.decay.chr_", filt.chr,".png"),
    width = 480, height = 480, units = "px", pointsize = 12,
    bg = "white", res = NA)

plot (x = df.LD.decay$delta.bp.1 , y = df.LD.decay$R2, 
      main=str_c("LD.decay.", filt.chr),
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (df.LD.decay$delta.bp.1)), 
      ylim = c(0, max (df.LD.decay$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 

axis(side = 2, las = 1)
x2 <- max (df.LD.decay$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (df.LD.decay$delta.bp.1,df.LD.decay$R2, 
        pch = 20, cex=1.5, col="gray28") 
box()
dev.off()

plot (x = df.LD.decay$delta.bp.1 , y = df.LD.decay$R2, 
      main=str_c("LD.decay.", filt.chr),
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (df.LD.decay$delta.bp.1)), 
      ylim = c(0, max (df.LD.decay$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 

axis(side = 2, las = 1)
x2 <- max (df.LD.decay$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (df.LD.decay$delta.bp.1,df.LD.decay$R2, 
        pch = 20, cex=1.5, col="gray28") 
box()

})
}


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


