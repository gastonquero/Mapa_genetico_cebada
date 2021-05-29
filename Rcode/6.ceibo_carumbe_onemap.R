#####################################################################
# Datos para la construccion del mapa geneticos de Cebada x Carumbe #
# datos enviados por Ariel Castro                                   #
# de:	vontruch@fagro.edu.uy
# para:	gastonquero@gmail.com
# Cc:	lviega@fagro.edu.uy
# fecha:	28 de junio de 2018, 15:26
# asunto:	Marcadores Ceibo carumbe
# enviado por:	fagro.edu.uy
# Gaston Quero
# 19/7/2018 
############################################################


getwd()
setwd("E:/Proyecto_Cebada")

("devtools")
library(devtools)
install_github("augusto-garcia/onemap")
library(onemap)
library(raster)
# plot raster without axes and box

onemap_example_f2 <- read_onemap(dir="./genetic_map/Data/rawdata", 
                                 inputfile = "onemap_example_f2.raw")


onemap_ceibo.carumbe_RIL <- read_onemap(dir="./genetic_map/Data/rawdata", 
                                 inputfile = "ceibo_carumbe_onemap_2.txt")


class(onemap_ceibo.carumbe_RIL)

plot(onemap_example_f2)
plot(onemap_ceibo.carumbe_RIL)

# Segregation tests
f2_test <- test_segregation(onemap_example_f2)
print(f2_test)
Bonferroni_alpha(f2_test)

# Segregation tests
DM.SO_f2_test <- test_segregation (onemap_DM.SO_f2)
print(f2_test)
Bonferroni_alpha(f2_test)
plot(f2_test)

#### Esta es la matriz con la cual trabajo 
onemap_DM.SO_f2.0.8 <- read_onemap (dir="./Data/procdata", 
                               inputfile = "DM.SO.onemap_0.8.txt")

onemap_DM.SO_f2.0.8.1 <- read_onemap (dir="./Data/procdata", 
                                   inputfile = "DM.SO.onemap_0.8_chr_pos.txt")



DM.SO_f2_test.0.8 <- test_segregation (onemap_DM.SO_f2.0.8)
DM.SO_f2_test.0.8.1 <- test_segregation (onemap_DM.SO_f2.0.8.1)

print(DM.SO_f2_test.0.8)
print(DM.SO_f2_test.0.8.1)

# la función Bonferroni_alpha muestra el valor alfa que se debe considerar 
# para este número de loci si se aplica la 
# corrección de Bonferroni con alfa global de 0,05

Bonferroni_alpha(DM.SO_f2_test.0.8)

Bonferroni_alpha(DM.SO_f2_test.0.8.1)

# ver qué marcadores están distorsionados bajo el criterio de Bonferroni

plot (DM.SO_f2_test.0.8)


# El gráfico se explica por sí mismo: 
# los valores de p se transformaron utilizando -log10 (valores p)
# para una mejor visualización. 
# Una línea vertical muestra el umbral para las pruebas si se aplica
# la corrección de Bonferroni. 
# Se identifican pruebas significativas y no significativas.

# Por favor, recuerde que la corrección de Bonferroni 
# es conservadora, y también que el descarte de datos de marcador 
# puede no ser un buen enfoque para su análisis. 
# Este gráfico es solo para sugerir un criterio, 
# así que úsalo con precaución

# list of markers with non-distorted segregation
lnds <- select_segreg (DM.SO_f2_test.0.8 , distorted = FALSE, numbers = FALSE) 
length(lnd)

lds <- select_segreg (DM.SO_f2_test.0.8 , distorted = TRUE, numbers = TRUE) 
length (lds)

6000 - 1897

# si max.rf <0.5 podemos afirmar que los marcadores están vinculados.
# La puntuación LOD es la estadística utilizada para evaluar 
# la importancia de la prueba para max.rf = 0,50.

twopts_DM.SO_f2 <- rf_2pts (onemap_DM.SO_f2.0.8)



print(twopts_DM.SO_f2)

print(twopts_DM.SO_f2.1)

###############
(LOD_sug_DM.SO_f2 <- suggest_lod(onemap_DM.SO_f2.0.8))

(LOD_sug_DM.SO_f2.1 <- suggest_lod(onemap_DM.SO_f2.0.8.1))

twopts_DM.SO_f2a <- rf_2pts (onemap_DM.SO_f2.0.8, LOD =  6.701726 )


print(twopts_DM.SO_f2a)

print(twopts_DM.SO_f2.1a)

# Se va a usar dos estrategias 

# A : Usando solo la frecuencia de  recombinacion 
## A.1 Assigning markers to linkage groups


mark_all_DM.SO_f2.1 <- make_seq (twopts_DM.SO_f2.1, "1")
class(mark_all_DM.SO_f2)

mark_all_DM.SO_f2a <- make_seq (twopts_DM.SO_f2a, "all")
class(mark_all_DM.SO_f2a)
## A.2 Forming the groups

LGs_DM.SO_f2 <- group(mark_all_DM.SO_f2)

LGs_DM.SO_f2a <- group(mark_all_DM.SO_f2a)

LGs_DM.SO_f2

LGs_DM.SO_f2a

# B : Usando la frecuencia de  recombinacion y la referencia del genoma
# B.1 
twopts_DM.SO_f2.1 <- rf_2pts (onemap_DM.SO_f2.0.8.1)
twopts_DM.SO_f2.1a <- rf_2pts (onemap_DM.SO_f2.0.8.1, LOD =  6.701726)

#### chr. 1

     CHR1       <- make_seq (twopts_DM.SO_f2.1a, "1")
     CHR1_ord   <- order_seq (CHR1)
     CHR1_frame    <- make_seq(CHR1_ord, "force")
     CHR1_test_map <- map(CHR1_frame)

     
     
svg (filename="./Figures/Figures_GENO/DM.SO.chr1.svg", 
          width=7, 
          height=5, 
          pointsize=12)

rf_graph_table(CHR1_test_map, main = "chr1", inter = FALSE)

dev.off()




#### chr. 2

     CHR2 <- make_seq (twopts_DM.SO_f2.1a, "2")
     CHR2_ord    <- order_seq (CHR2)
     CHR2_frame    <- make_seq(CHR2_ord, "force")
     CHR2_test_map <- map(CHR2_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr2.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR2_test_map, main = "chr2", inter = FALSE)

     dev.off()  
     
 #### chr. 3

     CHR3 <- make_seq (twopts_DM.SO_f2.1a, "3")
     CHR3_ord    <- order_seq (CHR3)
     CHR3_frame    <- make_seq(CHR3_ord, "force")
     CHR3_test_map <- map(CHR3_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr3.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR3_test_map, main = "chr3", inter = FALSE)
     
     dev.off()     

#### chr.4

     CHR4 <- make_seq (twopts_DM.SO_f2.1a, "4")
     CHR4_ord    <- order_seq (CHR4)
     CHR4_frame    <- make_seq(CHR4_ord, "force")
     CHR4_test_map <- map(CHR4_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr4.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR4_test_map, main = "chr4", inter = FALSE)
     
     dev.off()     

#### chr. 5
     CHR5 <- make_seq (twopts_DM.SO_f2.1a, "5")
     CHR5_ord    <- order_seq (CHR5)
     CHR5_frame    <- make_seq(CHR5_ord, "force")
     CHR5_test_map <- map(CHR5_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr5.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR5_test_map, main = "chr5", inter = FALSE)
     
     dev.off()    

#### chr. 6
     CHR6 <- make_seq (twopts_DM.SO_f2.1a, "6")
     CHR6_ord    <- order_seq (CHR6)
     CHR6_frame    <- make_seq(CHR6_ord, "force")
     CHR6_test_map <- map(CHR6_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr6.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR6_test_map, main = "chr6", inter = FALSE)
     
     dev.off()    
     
#### chr. 7
     CHR7 <- make_seq (twopts_DM.SO_f2.1a, "7")
     CHR7_ord    <- order_seq (CHR7)
     CHR7_frame    <- make_seq(CHR7_ord, "force")
     CHR7_test_map <- map(CHR7_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr7.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR7_test_map, main = "chr7", inter = FALSE)
     
     dev.off()  

#### chr. 8
     CHR8 <- make_seq (twopts_DM.SO_f2.1a, "8")
     CHR8_ord    <- order_seq (CHR8)
     CHR8_frame    <- make_seq(CHR8_ord, "force")
     CHR8_test_map <- map(CHR8_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr8.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR8_test_map, main = "chr8", inter = FALSE)
     
     dev.off()     
     
#### chr. 9
     CHR9 <- make_seq (twopts_DM.SO_f2.1a, c9)
     print(CHR9)
     CHR9_ord    <- order_seq (CHR9)
     CHR9_frame    <- make_seq(CHR9_ord, "force")
     CHR9_test_map <- map(CHR9_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr9.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR9_test_map, main = "chr9", inter = FALSE)
     
     dev.off() 

#### chr. 10
     CHR10 <- make_seq (twopts_DM.SO_f2.1a, "10")
     CHR10_ord    <- order_seq (CHR10)
     CHR10_frame    <- make_seq(CHR10_ord, "force")
     CHR10_test_map <- map(CHR10_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr10.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR10_test_map, main = "chr10", inter = FALSE)
     
     dev.off()      

 #### chr. 11
     CHR11 <- make_seq (twopts_DM.SO_f2.1a, "11")
     CHR11_ord    <- order_seq (CHR11)
     CHR11_frame    <- make_seq(CHR11_ord, "force")
     CHR11_test_map <- map(CHR11_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr11.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR11_test_map, main = "chr11", inter = FALSE)
     
     dev.off()     
 
#### chr. 12
     CHR12 <- make_seq (twopts_DM.SO_f2.1a, "12")
     CHR12_ord    <- order_seq (CHR12)
     CHR12_frame    <- make_seq(CHR12_ord, "force")
     CHR12_test_map <- map(CHR12_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr12.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR12_test_map, main = "chr12", inter = FALSE)
     
     dev.off()     
 
#### chr. 13
     CHR13 <- make_seq (twopts_DM.SO_f2.1a, "13")
     CHR13_ord    <- order_seq (CHR13)
     CHR13_frame    <- make_seq(CHR13_ord, "force")
     CHR13_test_map <- map(CHR13_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr13.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR13_test_map, main = "chr13", inter = FALSE)
     
     dev.off()     
     
#### chr. 14
     CHR14 <- make_seq (twopts_DM.SO_f2.1a, "14")
     CHR14_ord    <- order_seq (CHR14)
     CHR14_frame    <- make_seq(CHR14_ord, "force")
     CHR14_test_map <- map(CHR14_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr14.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR14_test_map, main = "chr14", inter = FALSE)
     
     dev.off() 
     
#### chr. 15
     CHR15 <- make_seq (twopts_DM.SO_f2.1a, "15")
     CHR15_ord    <- order_seq (CHR15)
     CHR15_frame    <- make_seq(CHR15_ord, "force")
     CHR15_test_map <- map(CHR15_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr15.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR15_test_map, main = "chr15", inter = FALSE)
     
     dev.off()     
     
#### chr. 16
     CHR16 <- make_seq (twopts_DM.SO_f2.1a, "16")
     CHR16_ord    <- order_seq (CHR16)
     CHR16_frame    <- make_seq(CHR16_ord, "force")
     CHR16_test_map <- map(CHR16_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr16.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR16_test_map, main = "chr16", inter = FALSE)
     
     dev.off()     

#### chr. 17
     CHR17 <- make_seq (twopts_DM.SO_f2.1a, "17")
     CHR17_ord    <- order_seq (CHR17)
     CHR17_frame    <- make_seq(CHR17_ord, "force")
     CHR17_test_map <- map(CHR17_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr17.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR17_test_map, main = "chr17", inter = FALSE)
     
     dev.off()     
  
 #### chr. 18
     CHR18 <- make_seq (twopts_DM.SO_f2.1a, "18")
     CHR18_ord    <- order_seq (CHR18)
     CHR18_frame    <- make_seq(CHR18_ord, "force")
     CHR18_test_map <- map(CHR18_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr18.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR18_test_map, main = "chr18", inter = FALSE)
     
     dev.off()     
  
#### chr. 19
     CHR19 <- make_seq (twopts_DM.SO_f2.1a, "19")
     CHR19_ord    <- order_seq (CHR19)
     CHR19_frame    <- make_seq(CHR19_ord, "force")
     CHR19_test_map <- map(CHR19_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr19.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR19_test_map, main = "chr19", inter = FALSE)
     
     dev.off()     
     
#### chr. 20
     CHR20 <- make_seq (twopts_DM.SO_f2.1a, "20")
     CHR20_ord    <- order_seq (CHR20)
     CHR20_frame    <- make_seq(CHR20_ord, "force")
     CHR20_test_map <- map(CHR20_frame)
     
     svg (filename="./Figures/Figures_GENO/DM.SO.chr20.svg", 
          width=7, 
          height=5, 
          pointsize=12)
     
     rf_graph_table(CHR20_test_map, main = "chr20", inter = FALSE)
     
     dev.off()     
     
########### riple
CHR1_final <- CHR1_test_map
CHR2_final <- CHR2_test_map
CHR3_final <- CHR3_test_map
CHR4_final <- CHR4_test_map
CHR5_final <- CHR5_test_map
CHR6_final <- CHR6_test_map
CHR7_final <- CHR7_test_map
CHR8_final <- CHR8_test_map
CHR9_final <- CHR9_test_map
CHR10_final <- CHR10_test_map
CHR11_final <- CHR11_test_map
CHR12_final <- CHR12_test_map
CHR13_final <- CHR13_test_map
CHR14_final <- CHR14_test_map
CHR15_final <- CHR15_test_map
CHR16_final <- CHR16_test_map
CHR17_final <- CHR17_test_map
CHR18_final <- CHR18_test_map
CHR19_final <- CHR19_test_map
CHR20_final <- CHR20_test_map


ripple_seq(CHR1_final)
ripple_seq(CHR2_final)
ripple_seq(CHR3_final)
ripple_seq(CHR4_final)
ripple_seq(CHR5_final)
ripple_seq(CHR6_final)
ripple_seq(CHR7_final)
ripple_seq(CHR8_final)
ripple_seq(CHR9_final)
ripple_seq(CHR10_final)
ripple_seq(CHR11_final)
ripple_seq(CHR12_final)
ripple_seq(CHR13_final)
ripple_seq(CHR14_final)
ripple_seq(CHR15_final)
ripple_seq(CHR16_final)
ripple_seq(CHR17_final)
ripple_seq(CHR18_final)
ripple_seq(CHR19_final)
ripple_seq(CHR20_final)

maps_list <- list(CHR1_final
                  ,CHR2_final
                  ,CHR3_final
                  ,CHR4_final
                  ,CHR5_final
                  ,CHR6_final
                  ,CHR7_final
                  ,CHR8_final
                  ,CHR9_final
                  ,CHR10_final
                  ,CHR11_final
                  ,CHR12_final
                  ,CHR13_final
                  ,CHR14_final
                  ,CHR15_final
                  ,CHR16_final
                  ,CHR17_final
                  ,CHR18_final
                  ,CHR19_final
                  ,CHR20_final)
                  
draw_map(maps_list, names = FALSE, grid = FALSE, cex.mrk = 0.7) 
  

#### ##### arbitrario 
c1 <- 1:59
c2 <- 60:206
c3 <- 207:315
c4 <- 316:383
c5 <- 384:464
c6 <- 465:576
c7 <- 577:722
c8 <- 723:841
c9 <- 842:911
c10 <- 912:995
c11 <- 996:1055
c12 <- 1056:1120
c13 <- 1121:1276
c14 <- 1277:1342
c15 <- 1343:1423 
c16 <- 1424:1497
c17 <- 1498:1585
c18 <- 1586:1664
c19 <- 1665:1799
c20 <- 1800:1897



CHR1.fisc     <- make_seq (twopts_DM.SO_f2.1a, c1)
CHR1_fisc_map <- map(CHR1.fisc)

CHR2.fisc     <- make_seq (twopts_DM.SO_f2.1a, c2)
CHR2_fisc_map <- map(CHR2.fisc)

CHR3.fisc     <- make_seq (twopts_DM.SO_f2.1a, c3)
CHR3_fisc_map <- map(CHR3.fisc)

CHR4.fisc     <- make_seq (twopts_DM.SO_f2.1a, c4)
CHR4_fisc_map <- map(CHR4.fisc)

CHR5.fisc     <- make_seq (twopts_DM.SO_f2.1a, c5)
CHR5_fisc_map <- map(CHR5.fisc)

CHR6.fisc     <- make_seq (twopts_DM.SO_f2.1a, c6)
CHR6_fisc_map <- map(CHR6.fisc)

CHR7.fisc     <- make_seq (twopts_DM.SO_f2.1a, c7)
CHR7_fisc_map <- map(CHR7.fisc)

CHR8.fisc     <- make_seq (twopts_DM.SO_f2.1a, c8)
CHR8_fisc_map <- map(CHR8.fisc)

CHR9.fisc     <- make_seq (twopts_DM.SO_f2.1a, c9)
CHR9_fisc_map <- map(CHR9.fisc)

CHR10.fisc     <- make_seq (twopts_DM.SO_f2.1a, c10)
CHR10_fisc_map <- map(CHR10.fisc)

CHR11.fisc     <- make_seq (twopts_DM.SO_f2.1a, c11)
CHR11_fisc_map <- map(CHR11.fisc)

CHR12.fisc     <- make_seq (twopts_DM.SO_f2.1a, c12)
CHR12_fisc_map <- map(CHR12.fisc)

CHR13.fisc     <- make_seq (twopts_DM.SO_f2.1a, c13)
CHR13_fisc_map <- map(CHR13.fisc)

CHR14.fisc     <- make_seq (twopts_DM.SO_f2.1a, c14)
CHR14_fisc_map <- map(CHR14.fisc)

CHR15.fisc     <- make_seq (twopts_DM.SO_f2.1a, c15)
CHR15_fisc_map <- map(CHR15.fisc)

CHR16.fisc     <- make_seq (twopts_DM.SO_f2.1a, c16)
CHR16_fisc_map <- map(CHR16.fisc)

CHR17.fisc     <- make_seq (twopts_DM.SO_f2.1a, c17)
CHR17_fisc_map <- map(CHR17.fisc)

CHR18.fisc     <- make_seq (twopts_DM.SO_f2.1a, c18)
CHR18_fisc_map <- map(CHR18.fisc)

CHR19.fisc     <- make_seq (twopts_DM.SO_f2.1a, c19)
CHR19_fisc_map <- map(CHR19.fisc)
                     
CHR20.fisc     <- make_seq (twopts_DM.SO_f2.1a, c20)
CHR20_fisc_map <- map(CHR20.fisc)


maps_fisc_list <- list(CHR1_fisc_map
                  ,CHR2_fisc_map
                  ,CHR3_fisc_map
                  ,CHR4_fisc_map
                  ,CHR5_fisc_map
                  ,CHR6_fisc_map
                  ,CHR7_fisc_map
                  ,CHR8_fisc_map
                  ,CHR9_fisc_map
                  ,CHR10_fisc_map
                  ,CHR11_fisc_map
                  ,CHR12_fisc_map
                  ,CHR13_fisc_map
                  ,CHR14_fisc_map
                  ,CHR15_fisc_map
                  ,CHR16_fisc_map
                  ,CHR17_fisc_map
                  ,CHR18_fisc_map
                  ,CHR19_fisc_map
                  ,CHR20_fisc_map)

draw_map(maps_fisc_list, names = FALSE, grid = FALSE, cex.mrk = 0.7) 

write_map (maps_fisc_list, "./Data/procdata/maps_fisc_list.map")



file <- paste(,"maps_fisc_list.raw", sep="/")
dat1 <- read.cross("mm", file=file, mapfile="mapmaker_example_f2.map")
newmap <- est.map(dat1, tol=1e-6, map.function="kosambi")


rf_graph_table(CHR19_test_map, main = "chr19", inter = FALSE)
rf_graph_table(CHR20_test_map, main = "chr20", inter = FALSE)

  LG1_f2_final, LG2_f2_final, LG3_f2_final)
     
#### Todos los cromosomas juntos #######

     CHR_mks_DM.SO_f2.1a <- group_seq (twopts_DM.SO_f2.1a,
                                    seqs = list(CHR1=CHR1,
                                                CHR2=CHR2,
                                                CHR3=CHR3,
                                                CHR4=CHR4,
                                                CHR5=CHR5,
                                                CHR6=CHR6,
                                                CHR7=CHR7,
                                                CHR8=CHR8,
                                                CHR9=CHR9,
                                                CHR10=CHR10,
                                                CHR11=CHR11,
                                                CHR12=CHR12,
                                                CHR13=CHR13,
                                                CHR14=CHR14,
                                                CHR15=CHR15,
                                                CHR16=CHR16,
                                                CHR17=CHR17,
                                                CHR18=CHR18,
                                                CHR19=CHR19,
                                                CHR20=CHR20), 
                                    rm.repeated=FALSE)
CHR_mks_DM.SO_f2.1a
CHR_mks_DM.SO_f2.1a$repeated   
CHR_mks_DM.SO_f2.1a$sequences$CHR1

 twopts_DM.SO_f2.1 <- rf_2pts (onemap_DM.SO_f2.0.8.1)
     
     
mark_all_DM.SO_f2a <- make_seq (twopts_DM.SO_f2a, "all")
class(mark_all_DM.SO_f2a)
## A.2 Forming the groups

LGs_DM.SO_f2 <- group(mark_all_DM.SO_f2)

LGs_DM.SO_f2a <- group(mark_all_DM.SO_f2a)

LGs_DM.SO_f2

LGs_DM.SO_f2a