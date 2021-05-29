##############################################################################
# analis de fenotipo tesis 
# datos FPTA 
# Gaston Quero - Lorena Cammarota
# 6-7-2020
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
library(lmerTest)
library(jtools)
library(tidytext)
library(stringi)
library(stringr)
library(psych)
library(rstatix)
library(datarium)
library(lubridate)
library(devtools)
library(arm)
#### Cargar los datos 

EEMAC  <- read_delim (file = "./Data/rawdata/EEMAC.txt" ,
                      col_names = TRUE, delim = "\t", na = "NA")


#rm(list=ls())

# EEMAC2018A=read.table("./LORE1.txt", header = TRUE, sep = "\t", dec = ".") #los datos faltantes tienen que estar vacios
# colnames(x)
# colnames(x)<-c("parcela","block") rename columns
#head(EEMAC2018A)
#str(EEMAC2018A)  #LAS VARIABLES TIENEN QUE SER NUMERIC
#attach(EEMAC2018A)
options(digits = 2)

#table(entry); table(check)
#En caso de precisar cambiar DataType
#EEMAC2018A$Rinde_corr_12=as.numeric(Rinde_corr_12)

EEMAC$BLOCK <- as.factor(EEMAC$BLOCK)
EEMAC$CHECK <- as.factor(EEMAC$CHECK)
EEMAC$GENOTIPOS <- as.factor(EEMAC$GENOTIPOS)
EEMAC$LINEA <- as.factor(EEMAC$LINEA )
EEMAC$PLOT <- as.factor(EEMAC$PLOT)

str (EEMAC)


unique (EEMAC$GENOTIPOS )
length (unique (EEMAC$LINEA))


## analisis exploratorio ###
# nitrogeno soluble

ggdotplot (EEMAC, x = "LINEA", y = "NS",
           #fill = "LINEA",
           size=0.8,
           add = "mean_sd")

LL99.NS <- EEMAC %>%
        dplyr::filter (GENOTIPOS == "LL99" ) %>%
        dplyr::select (LINEA, PLOT,  BLOCK, CHECK, GENOTIPOS,NS)

unique (LL99.NS$CHECK)

ggdotplot (LL99, x = "LINEA", y = "NS",
           fill = "LINEA",
           size=0.8,
           add = "mean_sd")





ggboxplot(LL99, x = "LINEA", y = "NS",
          fill =  "LINEA", 
          add = "jitter")


out.data <- LL99.NS %>%
  group_by (LINEA) %>%
  identify_outliers("NS") 


#### modelos
##3

#medias.a <- LL99 %>%
            #group_by (LINEA)%>%
            #mean (NS, na.rm = TRUE)


EEMAC.NS.mod1 <- lm (NS ~ LINEA + BLOCK, data = EEMAC)

anova (EEMAC.NS.mod1)
em.NS.EEMAC <-  emmeans(EEMAC.NS.mod1,  ~ LINEA)

em.NS.EEMAC.sum <- summary(em.NS.EEMAC)





ggdotplot (em.NS.EEMAC.sum, x = "LINEA", y = "emmean",
           #fill = "LINEA",
           size=0.8,
           #add = "mean_sd"
           )
plot ( NS ~ LINEA, data =EEMAC)

hist (em.NS.EEMAC.sum$emmean)


points (emmean ~ LINEA, data =em.NS.EEMAC.sum,
        pch=16, cex=1.2, col="red")


class (x )
summary (EEMAC.NS.mod1 )

## modelo mixto

EEMAC.NS.mod2 <- lmer (NS ~ CHECK + BLOCK + (1|CHECK:GENOTIPOS), data = EEMAC)

blup_NS <- ranef (EEMAC.NS.mod2, condVar = TRUE)[[1]]

blup_NS <- blup_NS %>%
           dplyr::mutate ( id = row.names(blup_NS))%>%
           dplyr::select (id, everything())


SE_blup_NS <- se.ranef (EEMAC.NS.mod2)[[1]] 



blup_NS_EEMAC <- data.frame (blup_NS,
                              SE_blup_NS)


colnames(blup_NS_EEMAC) <- c("id", "ranef_NS","SE_NS")
                             
                         

imod2 <- fixef(EEMAC.NS.mod2)[1]


blup_NS_EEMAC1 <- blup_NS_EEMAC %>%
                 dplyr::mutate ( blup_NS =  ranef_NS + imod2 ) %>%
                 separate (id, into =c ("A", "B"), sep=":")
        
blup_NS_EEMAC1 <-blup_NS_EEMAC1 %>%
                 dplyr::filter (B != "LL99")

blup_NS_EEMAC2 <- blup_NS_EEMAC1 [,-1]


em.NS.EEMAC.mod.2 <-  emmeans(EEMAC.NS.mod2 ,  ~CHECK)

blue.NS <- summary (em.NS.EEMAC.mod.2)






coef(EEMAC.NS.mod2)[1]

blup_NS_EEMAC$Tobs <- abs(blup_NS_EEMAC$blup_NS)/blup_NS_EEMAC$SE_blup_NS


blup_Rinde_EEMAC18A <- blup_Rinde_EEMAC18A %>%
  dplyr::mutate (Blup)



anova (EEMAC.NS.mod2)


ranef(EEMAC.NS.mod2)

#EEMAC2018A <-x[c(1751:1927),c(1:35)] #creo un objeto cortando de la fila 1 a la 177, por 17 columnas
#str(EEMAC2018A)

#Rinde_corr_12
Rinde_EEMAC18A <- lmer(Rinde_corr_12 ~ check + (1|block) + (1|check:geno), data = EEMAC2018A)


#check=efecto fijo  #Block=random #Geno=random anidado dentro de check

             
blup_Rinde_EEMAC18A <- ranef (Rinde_EEMAC18A, condVar = TRUE)[[1]]


SE_Rinde_EEMAC18A <- se.ranef (Rinde_EEMAC18A)[[1]] 



blup_Rinde_EEMAC18A <- data.frame (blup_Rinde_EEMAC18A,
                                   SE_Rinde_EEMAC18A, 
                                   rep (fixef(Rinde_EEMAC18A)[1],153)) 





colnames(blup_Rinde_EEMAC18A) <- c("blup_Rinde_EEMAC18A","SE_Rinde_EEMAC18A", 
                                   "interceptrindeEEMAC18A")


coef(Rinde_EEMAC18A)[1]
blup_Rinde_EEMAC18A$Tobs <- abs(blup_Rinde_EEMAC18A$blup_Rinde_EEMAC18A)/blup_Rinde_EEMAC18A$SE_Rinde_EEMAC18A



blup_Rinde_EEMAC18A$pvalor <- pt (blup_Rinde_EEMAC18A$Tobs, df.residual(Rinde_EEMAC18A), lower.tail=F)


blup_Rinde_EEMAC18A$blup_RindeLINEA_EEMAC18A <- abs(blup_Rinde_EEMAC18A$blup_Rinde_EEMAC18A + blup_Rinde_EEMAC18A$interceptrindeEEMAC18A)





#Porcentaje 1a+2a
Porc1a2aEEMAC18A <- lmer(Porc_1a_2a ~ check + (1 | block) + (1 | check:geno), data = EEMAC2018A) 


blupsPorc1a2aEEMAC18A <- ranef(Porc1a2aEEMAC18A, condVar = TRUE)[[1]]
SE_Porc1a2aEEMAC18A <- se.ranef(Porc1a2aEEMAC18A)[[1]] 
blups_Porc1a2aEEMAC18A <- data.frame(blupsPorc1a2aEEMAC18A, SE_Porc1a2aEEMAC18A, rep(fixef(Porc1a2aEEMAC18A)[1],153)) 
colnames(blups_Porc1a2aEEMAC18A) <- c("blups_Porc1a2aEEMAC18A","SE_Porc1a2aEEMAC18A", "intercept1a2aEEMAC18A")
coef(Porc1a2aEEMAC18A)[1]
blups_Porc1a2aEEMAC18A$Tobs <- abs(blups_Porc1a2aEEMAC18A$blups_Porc1a2aEEMAC18A)/blups_Porc1a2aEEMAC18A$SE_Porc1a2aEEMAC18A
blups_Porc1a2aEEMAC18A$pvalor <- pt(blups_Porc1a2aEEMAC18A$Tobs, df.residual(Porc1a2aEEMAC18A), lower.tail=F)
blups_Porc1a2aEEMAC18A$Blup_LINEAporc1a2aEEMAC18A <- abs(blups_Porc1a2aEEMAC18A$blups_Porc1a2aEEMAC18A + blups_Porc1a2aEEMAC18A$intercept1a2aEEMAC18A)

#TGW_corr
TGW_corr_EEMAC2018A <- lmer(TGW_corr ~ check + (1 | block) + (1 | check:geno), data = EEMAC2018A) 
blupsTGW_corr_EEMAC2018A <- ranef(TGW_corr_EEMAC2018A, condVar = TRUE)[[1]]
SE_TGW_corr_EEMAC2018A <- se.ranef(TGW_corr_EEMAC2018A)[[1]] 
blupsTGW_corr_EEMAC2018A <- data.frame(blupsTGW_corr_EEMAC2018A, SE_TGW_corr_EEMAC2018A, rep(fixef(TGW_corr_EEMAC2018A)[1],153)) 
colnames(blupsTGW_corr_EEMAC2018A) <- c("blupsTGW_corr_EEMAC2018A","SE_TGW_corr_EEMAC2018A", "interceptTGWEEMAC18A")
coef(TGW_corr_EEMAC2018A)[1]
blupsTGW_corr_EEMAC2018A$Tobs <- abs(blupsTGW_corr_EEMAC2018A$blupsTGW_corr_EEMAC2018A)/blupsTGW_corr_EEMAC2018A$SE_TGW_corr_EEMAC2018A
blupsTGW_corr_EEMAC2018A$pvalor <- pt(blupsTGW_corr_EEMAC2018A$Tobs, df.residual(TGW_corr_EEMAC2018A), lower.tail=F)
blupsTGW_corr_EEMAC2018A$Blup_LINEATGWEEMAC18A <- abs(blupsTGW_corr_EEMAC2018A$blupsTGW_corr_EEMAC2018A + blupsTGW_corr_EEMAC2018A$interceptTGWEEMAC18A)

#GR_m2 
GR_m2_EEMAC18A <- lmer(GR_m2 ~ check + (1 | block) + (1 | check:geno), data = EEMAC2018A) 
blupGR_m2_EEMAC18A <- ranef(GR_m2_EEMAC18A, condVar = TRUE)[[1]]
SE_GR_m2_EEMAC18A <- se.ranef(GR_m2_EEMAC18A)[[1]] 
blupGR_m2_EEMAC18A <- data.frame(blupGR_m2_EEMAC18A, SE_GR_m2_EEMAC18A, rep(fixef(GR_m2_EEMAC18A)[1],153)) 
colnames(blupGR_m2_EEMAC18A) <- c("blupGR_m2_EEMAC18A","SE_GR_m2_EEMAC18A", "interceptGR_m2_EEMAC18A")
coef(GR_m2_EEMAC18A)[1]
blupGR_m2_EEMAC18A$Tobs <- abs(blupGR_m2_EEMAC18A$blupGR_m2_EEMAC18A)/blupGR_m2_EEMAC18A$SE_GR_m2_EEMAC18A
blupGR_m2_EEMAC18A$pvalor <- pt(blupGR_m2_EEMAC18A$Tobs, df.residual(GR_m2_EEMAC18A), lower.tail=F)
blupGR_m2_EEMAC18A$blupGR_LINEAGRm2EEMAC18A <- abs(blupGR_m2_EEMAC18A$blupGR_m2_EEMAC18A + blupGR_m2_EEMAC18A$interceptGR_m2_EEMAC18A)

#PI_Yield
PIYield_EEMAC18A <- lmer(PIYield ~ check + (1 | block) + (1 | check:geno), data = EEMAC2018A) 
blupPIYield_EEMAC18A <- ranef(PIYield_EEMAC18A, condVar = TRUE)[[1]]
SE_PIYield_EEMAC18A <- se.ranef(PIYield_EEMAC18A)[[1]] 
blupPIYield_EEMAC18A <- data.frame(blupPIYield_EEMAC18A, SE_PIYield_EEMAC18A, rep(fixef(PIYield_EEMAC18A)[1],153)) 
colnames(blupPIYield_EEMAC18A) <- c("blupPIYield_EEMAC18A","SE_PIYield_EEMAC18A", "interceptPIYield_EEMAC18A")
coef(PIYield_EEMAC18A)[1]
blupPIYield_EEMAC18A$Tobs <- abs(blupPIYield_EEMAC18A$blupPIYield_EEMAC18A)/blupPIYield_EEMAC18A$SE_PIYield_EEMAC18A
blupPIYield_EEMAC18A$pvalor <- pt(blupPIYield_EEMAC18A$Tobs, df.residual(PIYield_EEMAC18A), lower.tail=F)
blupPIYield_EEMAC18A$Blup_LINEAPIYield_EEMAC18A <- abs(blupPIYield_EEMAC18A$blupPIYield_EEMAC18A + blupPIYield_EEMAC18A$interceptPIYield_EEMAC18A)
##
#ESPIGA M2
espigasm2EEMAC18A <- lmer(espiga.m2 ~ check + (1 | block) + (1 | check:geno), data = EEMAC2018A) 
blupespigasm2EEMAC18A <- ranef(espigasm2EEMAC18A, condVar = TRUE)[[1]]
SE_espigasm2EEMAC18A <- se.ranef(espigasm2EEMAC18A)[[1]] 
blupespigasm2EEMAC18A <- data.frame(blupespigasm2EEMAC18A, SE_espigasm2EEMAC18A, rep(fixef(espigasm2EEMAC18A)[1],153)) 
colnames(blupespigasm2EEMAC18A) <- c("blupespigasm2EEMAC18A","SE_espigasm2EEMAC18A", "interceptespigasm2EEMAC18A")
coef(espigasm2EEMAC18A)[1]
blupespigasm2EEMAC18A$Tobs <- abs(blupespigasm2EEMAC18A$blupespigasm2EEMAC18A)/blupespigasm2EEMAC18A$SE_espigasm2EEMAC18A
blupespigasm2EEMAC18A$pvalor <- pt(blupespigasm2EEMAC18A$Tobs, df.residual(espigasm2EEMAC18A), lower.tail=F)
blupespigasm2EEMAC18A$Blup_LINEAespigasm2EEMAC18A <- abs(blupespigasm2EEMAC18A$blupespigasm2EEMAC18A + blupespigasm2EEMAC18A$interceptespigasm2EEMAC18A)

#GRANOS POR ESPIGA
granosESPIGA_EEMAC18A <- lmer(granos.espiga ~ check + (1 | block) + (1 | check:geno), data = EEMAC2018A) 
blupgranosESPIGA_EEMAC18A <- ranef(granosESPIGA_EEMAC18A, condVar = TRUE)[[1]]
SE_granosESPIGA_EEMAC18A <- se.ranef(granosESPIGA_EEMAC18A)[[1]] 
blupgranosESPIGA_EEMAC18A <- data.frame(blupgranosESPIGA_EEMAC18A, SE_granosESPIGA_EEMAC18A, rep(fixef(granosESPIGA_EEMAC18A)[1],153)) 
colnames(blupgranosESPIGA_EEMAC18A) <- c("blupgranosESPIGA_EEMAC18A","SE_granosESPIGA_EEMAC18A", "interceptespigasm2EEMAC18A")
coef(granosESPIGA_EEMAC18A)[1]
blupgranosESPIGA_EEMAC18A$Tobs <- abs(blupgranosESPIGA_EEMAC18A$blupgranosESPIGA_EEMAC18A)/blupgranosESPIGA_EEMAC18A$SE_granosESPIGA_EEMAC18A
blupgranosESPIGA_EEMAC18A$pvalor <- pt(blupgranosESPIGA_EEMAC18A$Tobs, df.residual(granosESPIGA_EEMAC18A), lower.tail=F)
blupgranosESPIGA_EEMAC18A$Blup_LINEAgranosESPIGAEEMAC18A <- abs(blupgranosESPIGA_EEMAC18A$blupgranosESPIGA_EEMAC18A + blupgranosESPIGA_EEMAC18A$interceptespigasm2EEMAC18A)

#CADA OBJETO QUE FUI CREANDO TIENE 6 COLUMNAS(BLUP, SE, INTERCEP...LA ULTIMA=BLUP DE LA LINEA)
#EL ARCHIVO ESTA MEDIO ENTREVERADO (LO S?)..
#A CONTINUACION LO QUE HAGO, ES UNIR TODOS LOS OBJETOS QUE FUI CREANDO 


BlupsAgroEEMAC18A <- cbind(blup_Rinde_EEMAC18A,blups_Porc1a2aEEMAC18A, blupsTGW_corr_EEMAC2018A, blupGR_m2_EEMAC18A, blupPIYield_EEMAC18A, 
                           blupespigasm2EEMAC18A, blupgranosESPIGA_EEMAC18A)#rbind y cbind une filas y columnas

#Y QUEDARME CON LA COLUMNA 6,12,18,24 Y ASI..QUE ES EL BLUP PARA CADA VARIABLE
BlupsAgroEEMAC18A <- BlupsAgroEEMAC18A[,-c(1:5, 7:11, 13:17, 19:23, 25:29, 31:35, 37:41)]
##BlupsCombinados$blupZ20 <- NULL  # para borrar una columna utilizando su nombre
write.table(BlupsAgroEEMAC18A, file= 'blupsAgroEEMAC2018A.txt')
#GUARDO UN ARCHIVO CON LOS BLUPS DEFINITIVOS QUE ES LA QUE VOY A CARGAR PARA HACER GWAS
