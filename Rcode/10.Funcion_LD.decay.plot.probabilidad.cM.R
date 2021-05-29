################################################
# este es el codigo para graficar la caida del LD
#
# Gaston Quero - Sebastian simondi
# 3/3/2020
#################################################

getwd ()
setwd ("C:/Users/Usuario/OneDrive/Documentos/Paper_Crop_Science")

## cargar los paquetes
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
library("assertive")

## Cargo las matrices de datos ############

start_time <- Sys.time()

# Time difference of 1.000327 min
## esto lo tenes setear
index.chr <- 1:20

# se cargan las matrices generadas antes 
# se genear una lista con 20 tibbles
GENO.crom <- lapply (index.chr, function (filtro) { 
  
  G <- read_delim (file=paste("./Data/procdata/plot.LD.decay.chr",filtro,".txt",
                              sep=""), delim = ",", 
                              na = "NA", quote = "\"",col_names = TRUE)
  G <- G %>%
       dplyr::mutate (chr = str_c("chr_",filtro))
  print (G)
  return (G)
})

df.geno.crom <- do.call (rbind, GENO.crom)
head (df.geno.crom )
#head(df.geno.crom)
#unique (df.geno.crom$chr)
## Argumentos


data=
  
chrom=NULL
ini = 1
l1= 0.25
l2=0.5
l3=0.75

 filtro = "chr_1"
   
data = df.geno.crom 
chrom= filtro
ini=1
l1= 0.25
l2=0.5
l3=0.75

run_LD_decay <- function (data, chrom=NULL, ini = 1, l1= 0.25, l2=0.5, l3=0.75){

# el primer index p
#index.chrom <- unique (data$chr)
  
# control de la clases de los objetos
assert_is_data.frame (data)
assert_is_character (chrom)
assert_is_numeric(ini)
assert_is_numeric(l1)
assert_is_numeric(l2)
assert_is_numeric(l3)
  
  
  if (any (is_non_positive(ini), na.rm = TRUE)) {
    # Throw an error
    stop ("ini contains non-positive values, so no puede caminar.")
  }
  
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
  
  
chrom <- chrom

ini <- ini

x1 <- data %>%
      dplyr::filter (chr ==  chrom )

#summary (x1)

x1.na <- x1 %>%
         dplyr:::filter (R2 != "NA")

#summary (x1.na)
 


### secuencia para el lapply

max_delta_bp <- max(x1.na$delta.bp)
max_delta_Mb <- max_delta_bp /1e6

index.Mb <- seq (1, (round(max_delta_Mb + 1,0)*10), 1)

#filt.Mb <- 590

dt.plot.LD <- lapply (index.Mb, function (filt.Mb) { 

x.ini.Mb <- x1.na %>%
            dplyr::filter (delta.bp <= ini * 1e5)
  
QQ <- quantile (x.ini.Mb$R2)

l1 <- l1
l2 <- l2
l3 <- l3
q1 <- quantile (QQ, l1)
q2 <- quantile (QQ, l2)
q3 <- quantile (QQ, l3)

  
x2 <- x1.na %>%
      dplyr::filter (delta.bp <= filt.Mb * 1e5) #%>%
      #dplyr::mutate (Mb = "0.1")

num.total <- nrow (x2) # este hay que usar despues.


# pp <- # el summary de x2$R2

p1 <- x2 %>%
      dplyr::filter (R2 > q3[[1]]) %>%
      dplyr::mutate (prob="HLD") #%>% ### este tiene que cambia 

np1 <- nrow (p1)/num.total

p2 <- x2 %>%
      dplyr::filter (R2 > q2[[1]]) %>% ### este tiene que cambia 
      dplyr::filter (R2 <= q3[[1]]) %>%
      dplyr::mutate (prob= "LD")### este tiene que cambia 

np2 <- nrow (p2)/num.total

p3 <- x2 %>%
      dplyr::filter (R2 >  q1[[1]]) %>% ### este tiene que cambia 
      dplyr::filter (R2 <= q2[[1]]) %>%
      dplyr::mutate (prob= "LE") ### este tiene que cambia 


np3 <- nrow (p3)/num.total

p4 <- x2 %>%
      dplyr::filter (R2 <= q1[[1]]) %>%
      dplyr::mutate (prob= "HLE") ### este tiene que cambia 

np4 <- nrow (p4)/num.total

#XX <- rbind (p1, p2,p3,p4) ## este para que esta??

XX.1 <- as.numeric (rbind (np1, np2, np3, np4))

XX.1.1 <- c("HLD","LD","LE","HLE")

XX.2 <-  data.frame (ratio = XX.1, class = XX.1.1, Mb=0.1*filt.Mb)

#XX.2$clase <- factor(df2$Genotype, levels = c("Genotype 2", "Genotype 3", "Genotype 1")).

# Convert the cyl variable to a factor
XX.2$class <- as.factor(XX.2$class )
XX.2$class <- factor(XX.2$class, levels = c("HLE", "LE", "LD", "HLD"))

print (XX.2)


})


df.plot.LD <- as_tibble (do.call (rbind, dt.plot.LD))

df.plot.LD <- df.plot.LD %>%
              dplyr::mutate(Chrom=chrom)

write_csv (df.plot.LD, path= str_c("./Data/procdata/df.plot.LD_", chrom,".csv"), 
           na = "NA", append = FALSE)


bplt <- ggbarplot (df.plot.LD , "Mb",  "ratio",
           title = str_c("LD.decay by interval_", chrom),
           fill = "class", 
           border ="white",
           #sort.val = "desc",
           x.text.angle = 90,
           xlab = "bin (0.1 Mb)",
           color = "class", palette =c("navyblue","royalblue3","orange", "red4"),
           label = FALSE, lab.col = "white", lab.pos = "in")

bplt1 <- bplt + 
         rremove ("x.text") 

print (bplt1)

return (df.plot.LD)

}

# Voy a correr todos los cromosomas 

index.chrom <- unique (df.geno.crom$chrom)


index.chrom <- c("chr_1", "chr_18")

dt.all.chrom <- lapply (index.chrom, function (filtro) {
  
  run_LD_decay (data = df.geno.crom , chrom= filtro, ini=1,l1= 0.25, l2=0.5, l3=0.75)
  
})

#names (dt.all.chrom) <- index.chrom

df.all.chrom <- do.call (rbind, dt.all.chrom)

write_delim (df.all.chrom, path = "./Data/procdata/freq_LD.txt", 
             delim = ",", na = "NA")


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




LD.df.geno.chr3 <- run_LD_decay (data = df.geno.crom , chrom="chr_3", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr4 <- run_LD_decay (data = df.geno.crom , chrom="chr_4", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr5 <- run_LD_decay (data = df.geno.crom , chrom="chr_5", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr6 <- run_LD_decay (data = df.geno.crom , chrom="chr_6", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr7 <- run_LD_decay (data = df.geno.crom , chrom="chr_7", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr8 <- run_LD_decay (data = df.geno.crom , chrom="chr_8", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr9 <- run_LD_decay (data = df.geno.crom , chrom="chr_9", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr10 <- run_LD_decay (data = df.geno.crom , chrom="chr_10", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr11 <- run_LD_decay (data = df.geno.crom , chrom="chr_11", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr12 <- run_LD_decay (data = df.geno.crom , chrom="chr_12", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr13 <- run_LD_decay (data = df.geno.crom , chrom="chr_13", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr14 <- run_LD_decay (data = df.geno.crom , chrom="chr_14", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr15 <- run_LD_decay (data = df.geno.crom , chrom="chr_15", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr16 <- run_LD_decay (data = df.geno.crom , chrom="chr_16", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr17 <- run_LD_decay (data = df.geno.crom , chrom="chr_17", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr18 <- run_LD_decay (data = df.geno.crom , chrom="chr_18", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr19 <- run_LD_decay (data = df.geno.crom , chrom="chr_19", ini=1,l1= 0.25, l2=0.5, l3=0.75)

LD.df.geno.chr20 <- run_LD_decay (data = df.geno.crom , chrom="chr_20", ini=1,l1= 0.25, l2=0.5, l3=0.75)












malos <- df.plot.LD %>%
         dplyr::filter (clase =="A")



buenos <- df.plot.LD %>%
          dplyr::filter (clase =="D") %>%
          dplyr::mutate (nada = cantidad /(1-cantidad))

head (buenos)
dim(buenos)
indice [1]
indice [2]
head (buenos)

buenos$cantidad [1] + D.1 [1,2] *0.1

buenos$cantidad [2]

buenos$cantidad [2] 

D.1 [1,2]

 D <- (buenos$cantidad [indice + 1]- buenos$cantidad [indice])/0.1
 class (D)
 #colnames (D ) <- "derivada"
 
 D.1 <- as.data.frame(D)
 colnames (D.1 ) <- "derivada"
 
 
 

 
 D.1 <- D.1 %>%
        dplyr::mutate (index =indice )  %>%
        dplyr::select (index, derivada)
 
 class (D.1)
 dim(D.1)
 View (D)
 dim( deriv)
 
 deriv <- cbind (indice, D)

head ( deriv, 20 )



 plot (D ~ indice, data=deriv  ,type="p",
       pch=16, col="Blue")
 
 
xmalos <- malos %>%
          dplyr::filter (cantidad >=0.7) 

xm <- min(xmalos$Mb)

head(malos)

plot (cantidad ~ Mb, data=malos ,
      pch=16, col="#004616")


plot (nada ~ Mb, 
      ylim=c(0,1),
      data=buenos ,
      pch=16, col="orange")


points (cantidad ~ Mb, data=buenos ,
        pch=16, col="#004616")


plot ((cantidad /(1-cantidad) ~ Mb, 
      ylim=c(0,1),
      data=buenos ,
      pch=16, col="orange")




abline (v=0, lwd=1, lty=1, col= "black")

abline (h=0.7, lwd=2, lty=1, col= "#004616")
abline (v=xm, lwd=2, lty=1, col= "#004616")


ggdotplot (malos, "Mb",  "cantidad",
           fill = "clase", 
           border ="white",
           color = "clase", palette =c("#004616"),
           label = FALSE, lab.col = "white", lab.pos = "in")


x3 <- x2 %>%
      dplyr::filter (R2 >= 0.8354771 ) ## este  valor debe cambian

num.ex <- nrow (x3)

p <- num.ex/(num.total *0.25)

x4 <- cbind (filtro, p, num.total)

return (x4)

}) 
#} ## aca termina Fun.plot.LD













#Fun.probabilidad <-  function { 
  
 #dt.X <- lapply (indice , function (filtro) { 

    
          x2 <- x1.na %>%
                dplyr::filter (delta.bp <= 1 * 1e5)
          
          summary (x2)
          
          num.total <- nrow (x2)
          
          x3 <-x2 %>%
               dplyr::filter (R2 >= 0.8354771 ) ## este  valor debe cambian
          
          num.ex <- nrow (x3)
          
          p <- num.ex/(num.total *0.25)
          
          x4 <- cbind (filtro, p, num.total)
          
          return (x4)
  
  #})
  
 print(dt.X)
 
 df.X <- as.data.frame (do.call (rbind,dt.X))
 
 
 plot (x = df.X$filtro , y = df.X$p, 
       #main=str_c("LD.decay.", filtro),
       pch = 20, 
       type ="n",
       xaxt="none",
       yaxt="none",
       axes = F,
       xlim = c(0, max (df.X$filtro)), 
       ylim = c(0, 1)
       #ylab = expression(LD ~ (r^2)),
       #xlab = expression(Distance ~ (Mb))
       ) 
 box()
 axis(side = 2, las = 1)
 x2 <- max (df.X$filtro, na.rm = TRUE)
 axis (side=1,at=seq(0,x2,1),las = 1)
 tail
 points (df.X$filtro, df.X$p, 
         pch = 20, cex=1, col="gray69")
 
 px <- df.X %>%
       dplyr::filter (p >= 0.4)
 
 max (px$filtro)
 
 
 abline (h=0.75, lwd=2, lty=3)
 
 abline (h=0.5, lwd=2, lty=1, col="black")
 
 abline (h=0.43, lwd=2, lty=2, col="red")
 
 abline (v= max (px$filtro), lwd=2, lty=2)
 
 abline (v= max (px$filtro), lwd=2, lty=1,col="blue" )
 
 dim (dt.X)
 
 
 
 
 # }# aca termina Fun.probabilidad 


summary (x1.na )


M1.1 <- x1.na %>%
      dplyr::filter (mrk.1 == "S01_56670080")

M2.1 <- x1.na %>%
      dplyr::filter (mrk.2 == "S01_56670080")


summary (x1.na)


########## 1 ##############
x1.na.1Mb <-  x1.na %>%
              dplyr::filter (delta.bp <= 1e6)

summary (x1.na.1Mb)

M1.2 <- x1.na.1Mb  %>%
        dplyr::filter (mrk.1 == "S01_56670080")

M2.2 <- x1.na.1Mb  %>%
        dplyr::filter (mrk.2 == "S01_56670080")

View (M2.2)

############# 2 ############

x1.na.100Kb <-  x1.na.1Mb %>%
                dplyr::filter (delta.bp <= 1e5)

summary (x1.na.100Kb)

M1.3 <-x1.na.100Kb   %>%
     dplyr::filter (mrk.1 == "S01_56670080")

M2.3 <- x1.na.100Kb  %>%
      dplyr::filter (mrk.2 == "S01_56670080")



########   3  #####
x1.na.1Mb.r2.100kb <- x1.na.1Mb %>%
                      dplyr::filter (R2 >= 0.1707791)


summary (x1.na.1Mb.r2.100kb)

########### 4 ############### 

x1.na.r2.100kb <- x1.na %>%
                  dplyr::filter (R2 >= 0.1707791)


summary (x1.na.r2.100kb)


View(x1.1Mb.r2.100kb)

M1 <-x1.1Mb.r2.100kb  %>%
  dplyr::filter (mrk.1 == "S01_56670080")

M2 <-x1.1Mb.r2.100kb %>%
  dplyr::filter (mrk.2 == "S01_56670080")


summary (x1.1Mb.r2.100kb)


######## 4 ##########3
x1.r2 <-  x1 %>%
          dplyr::filter (R2>=0.1708)

summary (x1.r2)

M1 <-x1.r2   %>%
        dplyr::filter (mrk.1 == "S01_56670080")

M2 <- x1.r2 %>%
      dplyr::filter (mrk.2 == "S01_56670080")

# 1Mb = 2116
# 
                dplyr::filter (delta.bp >= 0) %>%
                dplyr::filter (R2 >= 0.1706)

summary (x1.1Mb.b.r2)

x1.1Mb.a.r2 <-  x1 %>%
              dplyr::filter (delta.bp <= 1e6) %>%
              dplyr::filter (delta.bp >= 1e5) %>%
             dplyr::filter (R2 >= 0.1706)



summary (x1.1Mb.a.r2)
unique(x1.1Mb$mrk.1)





x1.1Mb.r2 <-  x1 %>%
              dplyr::filter (delta.bp <= 1e6)%>%
              dplyr::filter (R2 >= 0.171)

summary (x1.1Mb.r2)
unique(x1.1Mb$mrk.1)


#unique (data$chr) == list.crom
Fun.bin.Mb <- function (data =NULL, r2=0, b.i=0, step=1, uni ="Mb"){
  
 index.chr <- 1:20
 data <- as.data.frame(data)
 list.crom <- str_c("chr_",index.chr)
 unique (data$chr) == list.crom
 r2 <- r2
 b.i <- b.i
 step <- step
 uni <- uni
 
   data <- data %>%
           dplyr::mutate (delta.Mb = delta.bp/1e6) %>%
           dplyr::filter (R2 > r2)
 
 dt.bin.1 <- lapply (list.crom, function (filtro){
   
   dt.1 <- data %>%
           dplyr::filter (chr == filtro)
   
     xbin <- dt.1 %>%
             dplyr::select(delta.Mb)
     
   xf <- (round (max(xbin),0)-1)

   lisbin <- (seq (b.i ,xf, step))

   dt.bin.2 <- lapply (lisbin, function (x) { 

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

     
   })
   
 })
 #df.bin.1 <- do.call (rbind, dt.bin.1)
 #return ( df.bin.1 )
} ### aca termina Fun.bin.Mb
 
start.time <- Sys.time()

crom.bin  <- Fun.bin.Mb (data =df.geno.crom, r2=0, b.i=0, step=1, uni ="Mb")

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken


#View(mkr_100)

list.crom <- str_c("chr_",index.chr)

 dt.crom.bin <- lapply (index.chr , function (x){
 
   binX <- crom.bin  [[x]]
   binX1 <- do.call (rbind, binX)
})

df.crom.bin <- do.call (rbind,  dt.crom.bin)

list.crom <- str_c("chr_",index.chr)

head(df.crom.bin)

cHR.1.1Mb.r2 <- df.crom.bin %>%
             dplyr::filter(chr=="chr_1")%>%
             dplyr::filter(bin=="1Mb")%>%
             dplyr::filter (R2>= 0.171)


dim(cHR.1.1Mb)

summary (cHR.1.1Mb)
unique(cHR.1.1Mb$mrk.1)

cHR.1.1Mb.r2

dim(cHR.1.1Mb.r2)

summary (cHR.1.1Mb.r2)
unique(cHR.1.1Mb.r2$mrk.1)




#head (df.crom.bin)
lapply (list.crom  , function (filtro){
   
   plotX <- df.crom.bin %>%
            dplyr::filter (chr==filtro)
  
   ggboxplot (plotX , x = "bin", y = "R2",title =filtro,  width = 0.8)
   
 })


Fun.bin.Kb <- function (data =NULL, r2=0, b.i=0){
  
  index.chr <- 1:20
  data <- as.data.frame(data)
  list.crom <- str_c("chr_",index.chr)
r2 <- r2
b.i <- b.i
step <- 100
uni <- "Kb"

data <- data %>%
      dplyr::mutate (delta.Kb = delta.bp/1e3) %>%
      dplyr::filter (R2 > r2)


dt.bin.Kb.1 <- lapply (list.crom, function (filtro){

  
  dt.kb.1 <- data %>%
             dplyr::filter (chr == filtro)%>%
             dplyr::filter (bin == "1Mb")
  
  xbin.kb <- dt.kb.1 %>%
             dplyr::select (delta.Kb)
  
  xf.kb <- (round (max(xbin.kb),0))
  
  lisbin.kb <- (seq (b.i ,xf.kb, step))
  
dt.bin.Kb.2 <- lapply (lisbin.kb, function (x) { 
  
  if ( x == b.i ) {
    dt.kb.2 <-   dt.kb.1  %>%
      dplyr::filter ( delta.Kb > x  )%>%
      dplyr::filter ( delta.Kb <= x + step )%>%
      dplyr::mutate (bin = str_c(b.i + step,uni))
    #print(dt.2)
    return (dt.kb.2)
  }
  
  if ( x > b.i ) { 
    if ( x < xf.kb ) {
      dt.kb.2  <-  dt.kb.1  %>%
        dplyr::filter ( delta.Kb > x  )%>%
        dplyr::filter ( delta.Kb <= x + step )%>%
        dplyr::mutate (bin = str_c(x + step,uni))
      #print(dt.2 ) 
      return (dt.kb.2 )
    }
  } 
  
})


})

}


start.time <- Sys.time()

crom.bin.kb  <- Fun.bin.Kb (data =df.crom.bin, r2=0, b.i=0)

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#list. <- str_c("chr_",index.chr)

dt.crom.bin.kb <- lapply (index.chr , function (x){
  
  binX <- crom.bin.kb  [[x]]
  binX1 <- do.call (rbind, binX)
})

df.crom.bin.kb <- do.call (rbind,  dt.crom.bin.kb)

list.crom <- str_c("chr_",index.chr)
#head (df.crom.bin)
lapply (list.crom  , function (filtro){
  
  plotX <- df.crom.bin.kb %>%
    dplyr::filter (chr==filtro)
  
  ggboxplot (plotX , x = "bin", y = "R2",title =filtro,  width = 0.8)
  
})

#data <- df.crom.bin.kb
Fun.r2 <- function (data=NULL, r2=0, b.i=0) {
  
  index.chr <- 1:20
  data <- as.data.frame(data)
  list.crom <- str_c("chr_",index.chr)
  
  lapply (list.crom  , function (filtro){
    
  dt.x.1 <- data  %>%
            dplyr::filter (chr == filtro)%>%
            dplyr::filter (bin == "100Kb")
  
  r2.bin <-  dt.x.1$R2

  q50.r2.bin <-  quantile (r2.bin, .50, na.rm = TRUE)
  q25.r2.bin <-  quantile (r2.bin, .25, na.rm = TRUE)
  q75.r2.bin <-  quantile  (r2.bin, .75, na.rm = TRUE)


  
  h.R2  <- hist (r2.bin , breaks ="Sturges", plot = F)
  h.R2$density <- h.R2$counts/sum (h.R2$counts)*100
  
  
  plot (h.R2, freq=FALSE, xlab="r2",main = "",
        ylab="Percent of total", col=NULL,

        #,xlim = c(400,1000)
         xaxt='l')
  
  par(new=TRUE)
  
  hist (r2.bin, main=filtro, 
        freq=F, xlab="", ylab = "",
        axes=FALSE,col="gray85")
  
  lines (density (r2.bin, na.rm = T), lwd=2,
         col="navyblue", lty=1)
  
  dens.r2.bin <- density (r2.bin, na.rm = T)
  x1.r2.bin <- min (which (dens.r2.bin$x >= 0))
  x2.r2.bin <- max (which (dens.r2.bin$x <= q25.r2.bin))
  
  with (dens.r2.bin, polygon(x=c(x[c(x1.r2.bin,x1.r2.bin:x2.r2.bin,x2.r2.bin)]),
                         y= c(0, y[x1.r2.bin:x2.r2.bin], 0),border ='navyblue',
                         col=scales::alpha( 'navyblue',.5)))
  
  x3.r2.bin  <- min (which (dens.r2.bin$x >= q75.r2.bin))
  x4.r2.bin  <- max (which (dens.r2.bin$x > 0))
  with (dens.r2.bin, polygon(x=c(x[c(x3.r2.bin,x3.r2.bin:x4.r2.bin,x4.r2.bin)]),
                         y= c(0, y[x3.r2.bin:x4.r2.bin], 0),border ='navyblue',
                         col=scales::alpha( 'navyblue',.5)))
  
  abline (v=q50.r2.bin, col="navyblue", lwd=1.5, lty =2)
  
  box()
  
}) ## aca termina el lapply del plot 
  
  
  dt.q  <-  lapply (list.crom  , function (filtro){
    
    x1 <- data  %>%
              dplyr::filter (chr ==  filtro)%>%
              dplyr::filter (bin == "100Kb")
 
    r2.bin <-  x1$R2
    x2 <- x1 [1,]%>%
          dplyr::select (chr) %>%
          dplyr::mutate (Q1.r2 = round (quantile (r2.bin, .25, na.rm = TRUE),3 ))%>%
          dplyr::mutate (Q2.r2 = round (quantile (r2.bin, .50, na.rm = TRUE),3 ))%>%
          dplyr::mutate (Q3.r2 = round (quantile (r2.bin, .75, na.rm = TRUE),3 ))%>%
          dplyr::mutate (IQR.r2 = Q3.r2 - Q1.r2)
    
  }) ## aca termina el lapply de los quartiles
  
  df.q <- do.call (rbind,dt.q)
  
}

Tr.r2 <- Fun.r2 (data=df.crom.bin.kb, r2=0, b.i=0)

#q.data <- XX
#data <- df.geno.crom
#head(data)
Fun.Mb <- function(data = NULL, q.data = NULL) {
  index.chr <- 1:20
  data <- as.data.frame(data)
  list.crom <- str_c("chr_",index.chr)
  
  dt.mb  <-  lapply (list.crom  , function (filtro){
    
    Q1.xq <-  q.data %>%
      dplyr::filter (chr == filtro )%>%
      dplyr::select (Q1.r2)
    
    Q1.xq [[1]]
    
    Q1.data <- data %>%
      dplyr::filter (chr == filtro )%>%
      dplyr::filter (R2 >=   Q1.xq [[1]]) %>%
      dplyr::mutate (delta.Mb= delta.bp/1e6) %>%
      dplyr::arrange (delta.Mb)
    
    dMb <-  Q1.data$delta.Mb
    
    
  
    x3 <- Q1.data [1,]%>%
      dplyr::select (chr) %>%
      dplyr::mutate (Q1.Mb = round (quantile ( dMb, .25, na.rm = TRUE),3 ))%>%
      dplyr::mutate (Q2.Mb = round (quantile ( dMb, .50, na.rm = TRUE),3 ))%>%
      dplyr::mutate (Q3.Mb = round (quantile ( dMb, .75, na.rm = TRUE),3 ))%>%
      dplyr::mutate (IQR.Mb = Q3.Mb - Q1.Mb)
    
  }) ## aca termina el lapply de los quartiles
  
  
  df.Mb <- do.call (rbind, dt.mb)
  
}
  
  
  
Tr.Mb <- Fun.Mb ( data=df.geno.crom, q.data = Tr.r2 )
  
Fun.plot.LD.decay <- function (data=NULL, Mb=NULL, r2=NULL) {
  
  index.chr <- 1:20
   data <- as.data.frame(data)
   list.crom <- str_c("chr_",index.chr)
   lapply (list.crom  , function (filtro){ 

   ZZ <- data %>%
     dplyr::filter (chr==filtro)%>%
     dplyr::mutate (delta.Mb = delta.bp/1e6)
   
   r2.1 <- r2 %>%
           dplyr::filter (chr==filtro)
   
  
   Mb.1 <- Mb %>%
           dplyr::filter (chr==filtro)
   
 
   Q1.r2 <-   r2.1$Q1.r2 [[1]]
   Q2.r2 <-   r2.1$Q2.r2 [[1]]
   Q3.r2 <-   r2.1$Q3.r2 [[1]]
   
   
   Q1.Mb <-  Mb.1$Q1.Mb [[1]]
   Q2.Mb <-  Mb.1$Q2.Mb [[1]]
   Q3.Mb <-  Mb.1$Q3.Mb [[1]]
   
   plot (x = ZZ$delta.Mb , y = ZZ$R2, 
         main=str_c("LD.decay.", filtro),
         pch = 20, 
         type ="n",
         xaxt="none",
         yaxt="none",
         axes = F,
         xlim = c(0, max (ZZ$delta.Mb)), 
         ylim = c(0, 1),
         ylab = expression(LD ~ (r^2)),
         xlab = expression(Distance ~ (Mb))) 
   box()
   axis(side = 2, las = 1)
   x2 <- max (ZZ$delta.Mb, na.rm = TRUE)
   axis (side=1,at=seq(0,x2,1),las = 1)
   
   points (ZZ$delta.Mb, ZZ$R2, 
           pch = 20, cex=1, col="gray69")
   
   abline (h=Q1.r2, lwd=2, lty=3)
   #abline (h=Q2.r2, lwd=2, lty=3)
   #abline (h=Q3.r2, lwd=2, lty=3)
   
   
   
   abline (v=Q1.Mb, lwd=2, lty=1, col="red")
   abline (v=Q2.Mb, lwd=2, lty=1, col="black")
   abline (v=Q3.Mb, lwd=2, lty=1, col="blue")
   
   points (Q1.Mb, Q1.r2 , 
                  pch = 20, cex=1.5, col="red")
   
   points (Q2.Mb, Q1.r2 , 
           pch = 20, cex=1.5, col="black")
   
   points (Q3.Mb, Q1.r2 , 
           pch = 20, cex=1.5, col="blue")
   
   #points (Q2.Mb, Q1.r2 , 
    #       pch = 10, cex=1, col="black")
   
   
   #segments(x0=Q1.Mb, y0=Q1.r2, x1 = Q2.Mb, y1 = Q2.r2,
    #        col ="black", lty =1, lwd = 2)
   
  # segments(x0=Q2.Mb, y0=Q2.r2, x1 = Q3.Mb, y1 = Q3.r2,
   #         col ="black", lty =1, lwd = 2)
   
   
   #points (Q1.Mb, Q1.r2 , 
   #        pch = 16, cex=1, col="black")
   
   
   })
   
} 


Fun.plot.LD.decay ( data=df.geno.crom, Mb=Tr.Mb, r2=Tr.r2 )


sebin <- Tr.Mb %>%
         inner_join(Tr.r2, by="chr")


write_delim(sebin, path ="./Data/procdata/sebin.txt" , 
            delim = ", ", na = "NA")



#############################3

#######

plot.LD.decay.chr1.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr1.R2.txt",
                                      delim=",", quote = "\"", 
                                      na = "NA", col_names = TRUE)

head (plot.LD.decay.chr1.R2)


plot.LD.decay.chr1.R2 <- plot.LD.decay.chr1.R2%>%
                         dplyr::arrange(delta.bp)%>%
                         dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr1.R2$delta.bp.1 , y = plot.LD.decay.chr1.R2$R2, 
      main="LD.decay.chr1",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr1.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr1.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr1.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)
  

points (plot.LD.decay.chr1.R2$delta.bp.1, plot.LD.decay.chr1.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

##### voy por los boxplots

min (plot.LD.decay.chr1.R2$delta.bp.1)
max (plot.LD.decay.chr1.R2$delta.bp.1)
Fbin <-  function (dt = NULL,r2=NULL, b.i=0, step=1) {
  
b.i <- b.i
step <- step
r2 <- r2
dt <- dt %>%
       dplyr::filter (R2 > r2)
 
xbin <- dt$delta.bp.1

xf <- (round (max(xbin),0)-1)

lisbin <- (seq (b.i ,xf, step))

class (lisbin[1])

# seleccion 
dt.bin <- lapply (lisbin, function (filtro) { 
  #print(filtro)
  if ( filtro == b.i ) {
    datos <- dt %>%
      dplyr::filter ( delta.bp.1 > filtro  )%>%
      dplyr::filter ( delta.bp.1 <= filtro + step )%>%
      dplyr::mutate (bin = str_c(b.i +  step, "Mb"))
    print(datos)
    return (datos)
  }
  
  if ( filtro > b.i ) { 
    if ( filtro < xf ) {
      datos <- dt %>%
        dplyr::filter ( delta.bp.1 > filtro  )%>%
        dplyr::filter ( delta.bp.1 <= filtro + step )%>%
        dplyr::mutate (bin = str_c(filtro + step, "Mb"))
      print(datos) 
      return (datos)
    }
  
  } 
  
  #if ( filtro == xf ) {
   # print("termine")
  #} 

})
  
df.bin <- do.call(rbind,dt.bin )
}
  
ccc <- na.omit (plot.LD.decay.chr1.R2)
  

x1Mb <- plot.LD.decay.chr1.R2 %>%
        dplyr::filter (delta.bp.1 <= 1)%>%
        dplyr::mutate (delta.bp.1 = delta.bp.1*10)


trols <- Fbin (dt=x1Mb,
                 r2=0,
                 b.i=0,
                 step=1)

ggboxplot (trols , x = "bin", y = "R2", width = 0.8)



b1 <- trols %>%
      filter (bin=="1Mb")

summary(b1)


###

plot.LD.decay.chr2.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr2.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)

head (plot.LD.decay.chr2.R2)


plot.LD.decay.chr2.R2 <- plot.LD.decay.chr2.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr2.R2$delta.bp.1 , y = plot.LD.decay.chr2.R2$R2, 
      main="LD.decay.chr2",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr2.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr2.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr2.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr2.R2$delta.bp.1, plot.LD.decay.chr2.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

##### voy por los boxplots

min (plot.LD.decay.chr2.R2$delta.bp.1)
max (plot.LD.decay.chr2.R2$delta.bp.1)
Fbin <-  function (dt = NULL,r2=NULL, b.i=0, step=1) {
  
  b.i <- b.i
  step <- step
  r2 <- r2
  dt <- dt %>%
    dplyr::filter (R2 > r2)
  
  xbin <- dt$delta.bp.1
  
  xf <- (round (max(xbin),0)-1)
  
  lisbin <- (seq (b.i ,xf, step))
  
  class (lisbin[1])
  
  # seleccion 
  dt.bin <- lapply (lisbin, function (filtro) { 
    #print(filtro)
    if ( filtro == b.i ) {
      datos <- dt %>%
        dplyr::filter ( delta.bp.1 > filtro  )%>%
        dplyr::filter ( delta.bp.1 <= filtro + step )%>%
        dplyr::mutate (bin = str_c(b.i +  step, "Mb"))
      print(datos)
      return (datos)
    }
    
    if ( filtro > b.i ) { 
      if ( filtro < xf ) {
        datos <- dt %>%
          dplyr::filter ( delta.bp.1 > filtro  )%>%
          dplyr::filter ( delta.bp.1 <= filtro + step )%>%
          dplyr::mutate (bin = str_c(filtro + step, "Mb"))
        print(datos) 
        return (datos)
      }
      
    } 
    
    #if ( filtro == xf ) {
    # print("termine")
    #} 
    
  })
  
  df.bin <- do.call(rbind,dt.bin )
}

ccc <- na.omit (plot.LD.decay.chr2.R2)


trols.2 <- Fbin (dt=plot.LD.decay.chr2.R2,
                 r2=0,
                 b.i=0,
                 step=1)
ggboxplot (trols.2 , x = "bin", y = "R2", width = 0.8)


x2.1Mb <- plot.LD.decay.chr2.R2 %>%
          dplyr::filter (delta.bp.1 <= 1)%>%
          dplyr::mutate (delta.bp.1 = delta.bp.1*10)


summary (x2.1Mb)


trols.2.1 <- Fbin (dt=x2.1Mb,
                 r2=0,
                 b.i=0,
                 step=1)

ggboxplot (trols.2.1 , x = "bin", y = "R2", width = 0.8)


b1.2 <- trols.2.1 %>%
         filter (bin=="1Mb")

summary(b1.2)
b1.2.1 <- b1.2 %>%
          filter (R2 > 0.24)

summary(b1.2.1)


trols.2.2 <- Fbin (dt=plot.LD.decay.chr2.R2,
               r2=0.24,
               b.i=0,
               step=1)

ggboxplot (trols.2.2 , x = "bin", y = "R2", 
           width = 0.8)


plot.LD.decay.chr4.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr4.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)

head (plot.LD.decay.chr4.R2)


plot.LD.decay.chr4.R2 <- plot.LD.decay.chr4.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr4.R2$delta.bp.1 , y = plot.LD.decay.chr4.R2$R2, 
      main="LD.decay.chr4",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr4.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr4.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr4.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr4.R2$delta.bp.1, plot.LD.decay.chr4.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

##### voy por los boxplots

min (plot.LD.decay.chr4.R2$delta.bp.1)
max (plot.LD.decay.chr4.R2$delta.bp.1)

ccc <- na.omit (plot.LD.decay.chr4.R2)


trols.4 <- Fbin (dt=plot.LD.decay.chr4.R2,
                 r2=0,
                 b.i=0,
                 step=1)
ggboxplot (trols.4 , x = "bin", y = "R2", width = 0.8)


x4.1Mb <- plot.LD.decay.chr4.R2 %>%
  dplyr::filter (delta.bp.1 <= 1)%>%
  dplyr::mutate (delta.bp.1 = delta.bp.1*10)


summary (x4.1Mb)


trols.4.1 <- Fbin (dt=x4.1Mb,
                   r2=0,
                   b.i=0,
                   step=1)

ggboxplot (trols.4.1 , x = "bin", y = "R2", width = 0.8)


b1.4 <- trols.4.1 %>%
  filter (bin=="1Mb")

summary(b1.4)

b1.4.1 <- b1.4 %>%
          filter (R2 > 0.41)

summary(b1.4.1)


trols.4.2 <- Fbin (dt=plot.LD.decay.chr4.R2,
                   r2=0.41,
                   b.i=0,
                   step=1)

ggboxplot (trols.4.2 , x = "bin", y = "R2", 
           width = 0.8)




b1.2.1 <- trols.2.1 %>%
  filter (bin=="1Mb")






plot (trols.0$delta.bp.1, trols.$R2, 
      ylim=c(0,1))


, 
           width = 0.8)


b100K<-trols%>%
   dplyr::filter (bin=="1Mb")
summary (b100K)


plot (x = x1Mb$delta.bp.1 , y = x1Mb$R2, 
      main="LD.decay.chr1",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (x1Mb$delta.bp.1)), 
      ylim = c(0, max (x1Mb$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (x1Mb$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,10,1),las = 1)


points (x1Mb$delta.bp.1, x1Mb$R2, 
        pch = 20, cex=1.5, col="gray28")


summary (x1Mb)



unique(trols$bin)

names (trols)
put <- c("1Mb","2Mb","3Mb","4Mb","5Mb","6Mb","7Mb","8Mb")
put <- c("0.5Mb","1Mb","1.5Mb","2Mb","2.5Mb","3Mb","3.5Mb","4Mb","4.5Mb", "5Mb")
bin.1.8 <- trols %>%
         dplyr::filter (bin %in% put)  



bin.2 <- trols %>%
         dplyr::filter (bin=="2Mb")

summary (bin.1)
summary (bin.2)



ggboxplot (bin.1.8 , x = "bin", y = "R2", width = 0.8)



boxplot (bin.1.$R2)



boxplot.stats(bin.2$R2, coef = 1.5, do.conf = TRUE, do.out = TRUE)$out

bin.2 <- plot.LD.decay.chr1.R2 %>%
         dplyr::filter ( x < delta.bp.1)%>%
         dplyr::filter ( delta.bp.1 <= x+1)%>%
         dplyr::mutate (bin = str_c(x+1, "Mb"))





  
  
  plot.LD.decay.chr1.R2.bin <- rbind (bin.1, bin.2)
  
  ggboxplot (plot.LD.decay.chr1.R2.bin , x = "bin", y = "R2", width = 0.8)

######### chr 2 #######

plot.LD.decay.chr2.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr2.R2.txt",
                                    delim=",", quote = "\"", 
                                    na = "NA", col_names = TRUE)

head (plot.LD.decay.chr2.R2)


plot.LD.decay.chr2.R2 <- plot.LD.decay.chr2.R2%>%
                        dplyr::arrange(delta.bp)%>%
                        dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr2.R2$delta.bp.1 , y = plot.LD.decay.chr2.R2$R2, 
      main="LD.decay.chr2",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr2.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr2.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr2.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr2.R2$delta.bp.1, plot.LD.decay.chr2.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

######### ch3 #######

plot.LD.decay.chr3.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr3.R2.txt",
                                    delim=",", quote = "\"", 
                                    na = "NA", col_names = TRUE)


plot.LD.decay.chr3.R2 <- plot.LD.decay.chr3.R2%>%
                        dplyr::arrange(delta.bp)%>%
                       dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr3.R2$delta.bp.1 , y = plot.LD.decay.chr3.R2$R2, 
      main="LD.decay.chr3",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr3.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr3.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr3.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr3.R2$delta.bp.1, plot.LD.decay.chr3.R2$R2, 
        pch = 20, cex=1.5, col="gray28")


######### chr4 #######

plot.LD.decay.chr4.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr4.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr4.R2 <- plot.LD.decay.chr4.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr4.R2$delta.bp.1 , y = plot.LD.decay.chr4.R2$R2, 
      main="LD.decay.chr4",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr4.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr4.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr4.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr4.R2$delta.bp.1, plot.LD.decay.chr4.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

######### chr5 #######

plot.LD.decay.chr5.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr5.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr5.R2 <- plot.LD.decay.chr5.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr5.R2$delta.bp.1 , y = plot.LD.decay.chr5.R2$R2, 
      main="LD.decay.chr5",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr5.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr5.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr5.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr5.R2$delta.bp.1, plot.LD.decay.chr5.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

###
######### chr6 #######

plot.LD.decay.chr6.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr6.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr6.R2 <- plot.LD.decay.chr6.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr6.R2$delta.bp.1 , y = plot.LD.decay.chr6.R2$R2, 
      main="LD.decay.chr6",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr6.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr6.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr6.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr6.R2$delta.bp.1, plot.LD.decay.chr6.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

########### chr7 #######

plot.LD.decay.chr7.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr7.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr7.R2 <- plot.LD.decay.chr7.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr7.R2$delta.bp.1 , y = plot.LD.decay.chr7.R2$R2, 
      main="LD.decay.chr7",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr7.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr7.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr7.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr7.R2$delta.bp.1, plot.LD.decay.chr7.R2$R2, 
        pch = 20, cex=1.5, col="gray28")
##
######### chr8 #######

plot.LD.decay.chr8.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr8.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr8.R2 <- plot.LD.decay.chr8.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr8.R2$delta.bp.1 , y = plot.LD.decay.chr8.R2$R2, 
      main="LD.decay.chr8",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr8.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr8.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr8.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr8.R2$delta.bp.1, plot.LD.decay.chr8.R2$R2, 
        pch = 20, cex=1.5, col="gray28")
##
######### chr9 #######

plot.LD.decay.chr9.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr9.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr9.R2 <- plot.LD.decay.chr9.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr9.R2$delta.bp.1 , y = plot.LD.decay.chr9.R2$R2, 
      main="LD.decay.chr9",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr9.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr9.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr9.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr9.R2$delta.bp.1, plot.LD.decay.chr9.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

##
######### chr10 #######

plot.LD.decay.chr10.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr10.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr10.R2 <- plot.LD.decay.chr10.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr10.R2$delta.bp.1 , y = plot.LD.decay.chr10.R2$R2, 
      main="LD.decay.chr10",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr10.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr10.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr10.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr10.R2$delta.bp.1, plot.LD.decay.chr10.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

##
######### chr11 #######

plot.LD.decay.chr11.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr11.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr11.R2 <- plot.LD.decay.chr11.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr11.R2$delta.bp.1 , y = plot.LD.decay.chr11.R2$R2, 
      main="LD.decay.chr11",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr11.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr11.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr11.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr11.R2$delta.bp.1, plot.LD.decay.chr11.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

#
######### chr12 #######

plot.LD.decay.chr12.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr12.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr12.R2 <- plot.LD.decay.chr12.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr12.R2$delta.bp.1 , y = plot.LD.decay.chr12.R2$R2, 
      main="LD.decay.chr12",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr12.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr12.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr12.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr12.R2$delta.bp.1, plot.LD.decay.chr12.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

#
######### chr13 #######

plot.LD.decay.chr13.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr13.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr13.R2 <- plot.LD.decay.chr13.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr13.R2$delta.bp.1 , y = plot.LD.decay.chr13.R2$R2, 
      main="LD.decay.chr13",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr13.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr13.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr13.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr13.R2$delta.bp.1, plot.LD.decay.chr13.R2$R2, 
        pch = 20, cex=1.5, col="gray28")
#
######### chr14 #######

plot.LD.decay.chr14.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr14.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr14.R2 <- plot.LD.decay.chr14.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr14.R2$delta.bp.1 , y = plot.LD.decay.chr14.R2$R2, 
      main="LD.decay.chr14",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr14.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr14.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr14.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr14.R2$delta.bp.1, plot.LD.decay.chr14.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

##
######### chr15 #######

plot.LD.decay.chr15.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr15.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr15.R2 <- plot.LD.decay.chr15.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr15.R2$delta.bp.1 , y = plot.LD.decay.chr15.R2$R2, 
      main="LD.decay.chr15",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr15.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr15.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr15.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr15.R2$delta.bp.1, plot.LD.decay.chr15.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

#
######### chr16 #######

plot.LD.decay.chr16.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr16.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr16.R2 <- plot.LD.decay.chr16.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr16.R2$delta.bp.1 , y = plot.LD.decay.chr16.R2$R2, 
      main="LD.decay.chr16",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr16.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr16.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr16.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr16.R2$delta.bp.1, plot.LD.decay.chr16.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

#
######### chr17 #######

plot.LD.decay.chr17.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr17.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr17.R2 <- plot.LD.decay.chr17.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr17.R2$delta.bp.1 , y = plot.LD.decay.chr17.R2$R2, 
      main="LD.decay.chr17",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr17.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr17.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr17.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr17.R2$delta.bp.1, plot.LD.decay.chr17.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

#
######### chr18 #######

plot.LD.decay.chr18.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr18.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr18.R2 <- plot.LD.decay.chr18.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr18.R2$delta.bp.1 , y = plot.LD.decay.chr18.R2$R2, 
      main="LD.decay.chr18",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr18.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr18.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr18.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr18.R2$delta.bp.1, plot.LD.decay.chr18.R2$R2, 
        pch = 20, cex=1.5, col="gray28")

#
########## chr19 #######

plot.LD.decay.chr19.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr19.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr19.R2 <- plot.LD.decay.chr19.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr19.R2$delta.bp.1 , y = plot.LD.decay.chr19.R2$R2, 
      main="LD.decay.chr19",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr19.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr19.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr19.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr19.R2$delta.bp.1, plot.LD.decay.chr19.R2$R2, 
        pch = 20, cex=1.5, col="gray28")
#
######### chr20 #######

plot.LD.decay.chr20.R2 <- read_delim (file="./Data/procdata/plot.LD.decay.chr20.R2.txt",
                                     delim=",", quote = "\"", 
                                     na = "NA", col_names = TRUE)


plot.LD.decay.chr20.R2 <- plot.LD.decay.chr20.R2%>%
  dplyr::arrange(delta.bp)%>%
  dplyr::mutate (delta.bp.1 = delta.bp/1e6)


plot (x = plot.LD.decay.chr20.R2$delta.bp.1 , y = plot.LD.decay.chr20.R2$R2, 
      main="LD.decay.chr20",
      pch = 20, 
      type ="n",
      xaxt="none",
      yaxt="none",
      axes = F,
      xlim = c(0, max (plot.LD.decay.chr20.R2$delta.bp.1)), 
      ylim = c(0, max (plot.LD.decay.chr20.R2$R2, na.rm = TRUE)),
      ylab = expression(LD ~ (r^2)),
      xlab = expression(Distance ~ (Mb))) 
box()
axis(side = 2, las = 1)
x2 <- max (plot.LD.decay.chr20.R2$delta.bp.1, na.rm = TRUE)
axis (side=1,at=seq(0,x2,10),las = 1)


points (plot.LD.decay.chr20.R2$delta.bp.1, plot.LD.decay.chr20.R2$R2, 
        pch = 20, cex=1.5, col="gray28")