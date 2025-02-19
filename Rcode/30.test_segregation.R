#######################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: test_segregation.R                                            ##
## Contains: test_segregation_of_a_marker,                             ##
## test_segregation, plot.onemap_segreg_test,                          ##
## print.onemap_segreg_test, Bonferroni_alpha, select_segreg           ##
##                                                                     ##
## Written by Antonio Augusto Franco Garcia with minor modifications   ##
## by Marcelo Mollinari and Cristiane Taniguti                         ##
## copyright (c) 2015 Antonio Augusto Franco Garcia                    ##
##                                                                     ##
## First version: 2015/04/18                                           ##
## Last update: 2017/11/08                                             ##
## License: GNU General Public License version 3 or later              ##
##                                                                     ##
#######################################################################

globalVariables(c("Marker", "p.value"))

##' test_segregation_of_a_marker
##'
##' Applies the chi-square test to check if markers are following the
##' expected segregation pattern, i. e., 1:1:1:1 (A), 1:2:1 (B), 3:1 (C) and 1:1 (D)
##' according to OneMap's notation. It does not use Yate's correction.
##'
##' First, the function selects the correct segregation pattern, then it
##' defines the H0 hypothesis, and then tests it, together with percentage of
##' missing data.
##' 
##' @importFrom stats chisq.test cor.test dist na.omit qchisq rf setNames
##'
##' @param x an object of class \code{onemap}, with data and additional information.
##' @param marker the marker which will be tested for its segregation.
##'
##' @return a list with the H0 hypothesis being tested, the chi-square statistics,
##' the associated p-values, and the \% of individuals genotyped.
##'
##' ##'@examples
##' data(onemap_example_bc) # Loads a fake backcross dataset installed with onemap
##' test_segregation_of_a_marker(onemap_example_bc,1)
##'
##' data(onemap_example_out) # Loads a fake outcross dataset installed with onemap
##' test_segregation_of_a_marker(onemap_example_out,1)
##' @export
test_segregation_of_a_marker <- function(x, marker) {
  options(warn=-1) #Supress warnings from function chisq.test - when cell have < 5 elements
  ## Segregation pattern for each marker type
  p.a <- rep(1/4, 4); p.b <- c(1/4, 1/2, 1/4); p.c <- c(3/4, 1/4); p.d <- rep(1/2, 2)
  ## Counting each category, removing missing data (coded as 0)
  count <- table(x$geno[,marker], exclude=0)
  ## Do the chisq test, using the appropriate expected segregation
  ## grepl() allows finding the marker type (it has the letter in the argument)
  ## Some markers may not have data for some classes, so fill them with 0
  if (grepl("A.H.B",x$segr.type[marker])) {
    if (is.element(1,x$geno[,marker])) c1 <- count[names(count)==1] else c1 <- 0
    if (is.element(2,x$geno[,marker])) c2 <- count[names(count)==2] else c2 <- 0
    if (is.element(3,x$geno[,marker])) c3 <- count[names(count)==3] else c3 <- 0
    qui <- chisq.test(c(c1,c2,c3), p=p.b, correct = FALSE)
    H0 <- "1:2:1"
  }
  else if (grepl("C.A",x$segr.type[marker])) {
    if (is.element(1,x$geno[,marker])) c1 <- count[names(count)==1] else c1 <- 0
    if (is.element(5,x$geno[,marker])) c2 <- count[names(count)==5] else c2 <- 0
    qui <- chisq.test(as.vector(c(c1,c2)), p=rev(p.c), correct = FALSE)
    H0 <- "3:1"
  }
  else if (grepl("D.B",x$segr.type[marker])) {
    if (is.element(3,x$geno[,marker])) c1 <- count[names(count)==3] else c1 <- 0
    if (is.element(4,x$geno[,marker])) c2 <- count[names(count)==4] else c2 <- 0
    qui <- chisq.test(as.vector(c(c1,c2)), p=rev(p.c), correct = FALSE)
    H0 <- "3:1"
  }
  else if (grepl("A.H",x$segr.type[marker])) {
    if (is.element(1,x$geno[,marker])) c1 <- count[names(count)==1] else c1 <- 0
    if (is.element(2,x$geno[,marker])) c2 <- count[names(count)==2] else c2 <- 0
    qui <- chisq.test(as.vector(c(c1,c2)), p=p.d, correct = FALSE)
    H0 <- "1:1"
  }
  else if (grepl("A.B",x$segr.type[marker])) {
    if (is.element(1,x$geno[,marker])) c1 <- count[names(count)==1] else c1 <- 0
    if (is.element(3,x$geno[,marker])) c2 <- count[names(count)==3] else c2 <- 0
    qui <- chisq.test(as.vector(c(c1,c2)), p=p.d, correct = FALSE)
    H0 <- "1:1"
  }
  else if (grepl("A",x$segr.type[marker]) & class(x)[2]=="outcross") {
    if (is.element(1,x$geno[,marker])) c1 <- count[names(count)==1] else c1 <- 0
    if (is.element(2,x$geno[,marker])) c2 <- count[names(count)==2] else c2 <- 0
    if (is.element(3,x$geno[,marker])) c3 <- count[names(count)==3] else c3 <- 0
    if (is.element(4,x$geno[,marker])) c4 <- count[names(count)==4] else c4 <- 0
    qui <- chisq.test(as.vector(c(c1,c2,c3,c4)), p=p.a, correct = FALSE)
    H0 <- "1:1:1:1"
  }
  else if (grepl("B",x$segr.type[marker]) & class(x)[2]=="outcross") {
    if (is.element(1,x$geno[,marker])) c1 <- count[names(count)==1] else c1 <- 0
    if (is.element(2,x$geno[,marker])) c2 <- count[names(count)==2] else c2 <- 0
    if (is.element(3,x$geno[,marker])) c3 <- count[names(count)==3] else c3 <- 0
    qui <- chisq.test(as.vector(c(c1,c2,c3)), p=p.b, correct = FALSE)
    H0 <- "1:2:1"
  }
  else if (grepl("C",x$segr.type[marker]) & class(x)[2]=="outcross") {
    if (is.element(1,x$geno[,marker])) c1 <- count[names(count)==1] else c1 <- 0
    if (is.element(2,x$geno[,marker])) c2 <- count[names(count)==2] else c2 <- 0
    qui <- chisq.test(as.vector(c(c1,c2)), p=p.c, correct = FALSE)
    H0 <- "3:1"
  }
  else if (grepl("D",x$segr.type[marker]) & class(x)[2]=="outcross") {
    if (is.element(1,x$geno[,marker])) c1 <- count[names(count)==1] else c1 <- 0
    if (is.element(2,x$geno[,marker])) c2 <- count[names(count)==2] else c2 <- 0
    qui <- chisq.test(as.vector(c(c1,c2)), p=p.d, correct = FALSE)
    H0 <- "1:1"
  }
  #impossible to test: dominant and co-dominant mixed in the same marker
  #however, it will not work with "qui <- NA"; using NULL instead
  else if (grepl("M.X",x$segr.type[marker])) {
    qui <- NULL
    qui$statistic <- NA
    qui$p.value <- NA
    H0 <- "Check data!"
  }
  
  return(list(Hypothesis=H0, qui.quad=qui$statistic, p.val=qui$p.value,
              perc.genot=100*(sum(table(x$geno[,marker], exclude=0))/x$n.ind)))
}
##'

##' test_segregation
##'
##' Using OneMap internal function test_segregation_of_a_marker(),
##' performs the Chi-square test to check if all markers in a dataset are following
##' the expected segregation pattern, i. e., 1:1:1:1 (A), 1:2:1 (B), 3:1 (C) and 1:1 (D)
##' according to OneMap's notation.
##'
##' First, it identifies the correct segregation pattern and corresponding H0 hypothesis,
##' and then tests it.
##'
##' @param x an object of class \code{onemap}, with data and additional information.
##'
##' @return an object of class onemap_segreg_test, which is a list with marker name,
##' H0 hypothesis being tested, the chi-square statistics, the associated p-values
##' and the \% of individuals genotyped. To see the object, it is necessary to print
##' it.
##'
##' @examples
##' data(onemap_example_out) # Loads a fake outcross dataset installed with onemap
##' Chi <- test_segregation(onemap_example_out) # Performs the chi-square test for all markers
##' print(Chi) # Shows the results
##'
##' @export
test_segregation <- function(x) {
  if (is(x,"onemap")) {
    y <- list(Marker=dimnames(x$geno)[[2]],
              Results.of.tests=sapply(1:x$n.mar, function(onemap.object, marker)
                test_segregation_of_a_marker(onemap.object, marker),
                onemap.object=x))
    # sapply iterates from 1 to x$n.mar; x is fixed (onemap object with data)
    class(y) <- c("onemap_segreg_test")
    invisible(y) #returns y without showing it
  }
  else stop("This is not an OneMap object with raw data")
}
##'

##' Show the results of segregation tests
##'
##' It shows the results of Chisquare tests performed for all markers in a onemap object
##' of cross type outcross, backcross, F2 intercross or recombinant inbred lines.
##'
##' @param x an object of class onemap_segreg_test
##'
##' @param ... currently ignored
##'
##' @return a dataframe with marker name, H0 hypothesis, chi-square statistics,
##' p-values, and % of individuals genotyped.
##'
##' @examples
##' data(onemap_example_out) # Loads a fake outcross dataset installed with onemap
##' Chi <- test_segregation(onemap_example_out) # Performs the chi-square test for all markers
##' print(Chi) # Shows the results
##'
##' @method print onemap_segreg_test
##' @export
print.onemap_segreg_test <- function(x,...) {
  Z <- data.frame(Marker=x$Marker,
                  H0=unlist(x$Results.of.tests[1,]),
                  Chi.square=unlist(x$Results.of.tests[2,]),
                  p.value=unlist(x$Results.of.tests[3,]),
                  Perc.genot=round(unlist(x$Results.of.tests[4,]),2))
  colnames(Z) <- c("Marker","H0","Chi-square","p-value","% genot.")
  return(Z)
}
##'

##' Plot p-values for chi-square tests of expected segregation
##'
##' Draw a graphic showing the p-values (re-scaled to -log10(p-values)) associated with the
##' chi-square tests for the expected segregation patterns for all markers in a dataset.
##' It includes a vertical line showing the threshold for declaring statistical significance
##' if Bonferroni's correction is considered, as well as the percentage of markers that
##' will be discarded if this criterion is used.
##'
##' @param x an object of class onemap_segreg_test (produced by onemap's function
##' test_segregation()), i. e., after performing segregation tests
##' @param order a variable to define if p-values will be ordered in the plot
##'
##' @param ... currently ignored
##'
##' @return a ggplot graphic
##'
##' @import ggplot2
##'
##' @examples
##'
##' data(onemap_example_bc) # load OneMap's fake dataset for a backcross population
##' BC.seg <- test_segregation(onemap_example_bc) # Applies chi-square tests
##' print(BC.seg) # Shows the results
##' plot(BC.seg) # Plot the graph, ordering the p-values
##' plot(BC.seg, order=FALSE) # Plot the graph showing the results keeping the order in the dataset
##' # You can store the graphic in an object, then save it.
##' # For details, see the help of ggplot2's function ggsave()
##' # g <- plot(BC.seg)
##' # ggplot2::ggsave("SegregationTests.jpg", g, width=7, height=5, dpi=600)
##'
##' data(onemap_example_out) # load OneMap's fake dataset for an outcrossing population
##' Out.seg <- test_segregation(onemap_example_out) # Applies chi-square tests
##' print(Out.seg) # Shows the results
##' plot(Out.seg) # Plot the graph, ordering the p-values
##' plot(Out.seg, order=FALSE) # Plot the graph showing the results keeping the order in the dataset
##' # You can store the graphic in an object, then save it.
##' # For details, see the help of ggplot2's function ggsave()
##' g <- plot(Out.seg)
##' ggplot2::ggsave("SegregationTests.jpg", g, width=7, height=5, dpi=600)
##' @method plot onemap_segreg_test
##' @export
plot.onemap_segreg_test <- function(x, order=TRUE,...) {
  # Create a data frame
  Z <- data.frame(Marker=x$Marker,
                  X.square=unlist(x$Results.of.tests[2,]),
                  p.value=unlist(x$Results.of.tests[3,]))
  Bonf <- -log10(.05/nrow(Z)) #Bonferroni's threshold'
  Z$signif <- factor(ifelse(-log10(Z$p.value)<Bonf,"non sign.","sign."))
  Z$order <- 1:nrow(Z)
  # % of distorted
  perc <- 100*(1-(table(Z$signif)[1]/nrow(Z)))
  # Keeping markers in their original order (not alphanumeric), or by p-values (default)
  if (order!=TRUE) Z$Marker <- factor(Z$Marker, levels = Z$Marker[order(Z$order)])
  else Z$Marker <- factor(Z$Marker, levels = Z$Marker[order(Z$p.value, decreasing=TRUE)])
  # Plotting
  g <- ggplot(data=Z, aes(x=Marker, y=-log10(p.value)))
  g <- g + ylab(expression(-log[10](p-value)))
  g <- g + geom_point(aes(color=signif), stat="identity",size=2.5)
  g <- g + scale_colour_manual(name=paste("Bonferroni\n","(",round(perc,0),"% distorted)",sep=""),
                               values = c("#46ACC8","#B40F20"))
  g <- g + geom_hline(yintercept = Bonf, colour="#E58601", linetype = "longdash")
  g <- g + coord_flip()
  if (nrow(Z)>30) g <- g + theme(axis.text.y = element_blank())
  g
}
##'


##' Calculates individual significance level to be used to achieve a global alpha (with Bonferroni)
##'
##' It shows the alpha value to be used in each chi-square segregation test, in order to achieve
##' a given global type I error. To do so, it uses Bonferroni's criteria.
##'
##'@importFrom methods is
##'
##' @param x an object of class onemap_segreg_test
##' @param global.alpha the global alpha that
##'
##' @return the alpha value for each test (numeric)
##'
##' @examples
##' data(onemap_example_bc) # Loads a fake backcross dataset installed with onemap
##' Chi <- test_segregation(onemap_example_bc) # Performs the chi-square test for all markers
##' print(Chi) # Shows the results of the Chi-square tests
##' Bonferroni_alpha (Chi) # Shows the individual alpha level to be used
##'
##' @export
Bonferroni_alpha <- function(x, global.alpha=0.05) {
  if (!is(x,"onemap_segreg_test")) stop("This is not an object of class onemap_segreg_test")
  alpha.Bonf <- global.alpha/length(x$Marker)
  return(alpha.Bonf)
}
##'

##' Show markers with/without segregation distortion
##'
##' A function to shows which marker have segregation distortion if Bonferroni's correction is
##' applied for the Chi-square tests of mendelian segregation.
##'
##' @param x an object of class onemap_segreg_test
##' @param distorted a TRUE/FALSE variable to show distorted or non-distorted markers
##' @param numbers a TRUE/FALSE variable to show the numbers or the names of the markers
##'
##' @return a vector with marker names or numbers, according to the option for "distorted" and "numbers"
##'
##' @examples
##' # Loads a fake backcross dataset installed with onemap
##' data(onemap_example_bc)
##' # Performs the chi-square test for all markers
##' Chi <- test_segregation(onemap_example_bc)
##' # To show non-distorted markers
##' select_segreg(Chi)
##' # To show markers with segregation distortion
##' select_segreg(Chi, distorted=TRUE)
##' # To show the numbers of the markers with segregation distortion
##' select_segreg(Chi, distorted=TRUE, numbers=TRUE)
##'
##' @export
select_segreg <- function(x, distorted=FALSE, numbers=FALSE) {
  if (!is(x,"onemap_segreg_test")) stop("This is not an object of class onemap_segreg_test")
  Z <- data.frame(Marker=x$Marker,
                  p.value=unlist(x$Results.of.tests[3,]))
  if (distorted==FALSE) Z <- subset(Z, p.value>=Bonferroni_alpha(x))
  else Z <- subset(Z, p.value<Bonferroni_alpha(x))
  if (numbers==TRUE) return(which(x$Marker %in% as.vector(Z[,1])))
  else return(as.vector(Z[,1]))
}
##'
