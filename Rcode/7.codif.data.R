codif_data
function (geno.in, segr.type.in, cross = c("outcross", "f2", 
                                           "backcross", "riself", "risib")) 
{
  cross <- match.arg(cross)
  geno.out <- matrix(NA, nrow(geno.in), ncol(geno.in))
  segr.type.out <- rep(NA, length(segr.type.in))
  geno.out[is.na(geno.in)] <- 0
  for (i in 1:length(segr.type.in)) {
    if (cross == "outcross") {
      switch(EXPR = segr.type.in[i], A.1 = {
        geno.out[which(geno.in[, i] == "ac"), i] <- 1
        geno.out[which(geno.in[, i] == "ad"), i] <- 2
        geno.out[which(geno.in[, i] == "bc"), i] <- 3
        geno.out[which(geno.in[, i] == "bd"), i] <- 4
        segr.type.out[i] <- 1
      }, A.2 = {
        geno.out[which(geno.in[, i] == "a"), i] <- 1
        geno.out[which(geno.in[, i] == "ac"), i] <- 2
        geno.out[which(geno.in[, i] == "ba"), i] <- 3
        geno.out[which(geno.in[, i] == "bc"), i] <- 4
        segr.type.out[i] <- 1
      }, A.3 = {
        geno.out[which(geno.in[, i] == "ac"), i] <- 1
        geno.out[which(geno.in[, i] == "a"), i] <- 2
        geno.out[which(geno.in[, i] == "bc"), i] <- 3
        geno.out[which(geno.in[, i] == "b"), i] <- 4
        segr.type.out[i] <- 1
      }, A.4 = {
        geno.out[which(geno.in[, i] == "ab"), i] <- 1
        geno.out[which(geno.in[, i] == "a"), i] <- 2
        geno.out[which(geno.in[, i] == "b"), i] <- 3
        geno.out[which(geno.in[, i] == "o"), i] <- 4
        segr.type.out[i] <- 1
      }, B1.5 = {
        geno.out[which(geno.in[, i] == "a"), i] <- 1
        geno.out[which(geno.in[, i] == "ab"), i] <- 2
        geno.out[which(geno.in[, i] == "b"), i] <- 3
        segr.type.out[i] <- 2
      }, B2.6 = {
        geno.out[which(geno.in[, i] == "a"), i] <- 1
        geno.out[which(geno.in[, i] == "ab"), i] <- 2
        geno.out[which(geno.in[, i] == "b"), i] <- 3
        segr.type.out[i] <- 3
      }, B3.7 = {
        geno.out[which(geno.in[, i] == "a"), i] <- 1
        geno.out[which(geno.in[, i] == "ab"), i] <- 2
        geno.out[which(geno.in[, i] == "b"), i] <- 3
        segr.type.out[i] <- 4
      }, C.8 = {
        geno.out[which(geno.in[, i] == "a"), i] <- 1
        geno.out[which(geno.in[, i] == "o"), i] <- 2
        segr.type.out[i] <- 5
      }, D1.9 = {
        geno.out[which(geno.in[, i] == "ac"), i] <- 1
        geno.out[which(geno.in[, i] == "bc"), i] <- 2
        segr.type.out[i] <- 6
      }, D1.10 = {
        geno.out[which(geno.in[, i] == "a"), i] <- 1
        geno.out[which(geno.in[, i] == "ab"), i] <- 2
        segr.type.out[i] <- 6
      }, D1.11 = {
        geno.out[which(geno.in[, i] == "a"), i] <- 1
        geno.out[which(geno.in[, i] == "b"), i] <- 2
        segr.type.out[i] <- 6
      }, D1.12 = {
        geno.out[which(geno.in[, i] == "ab"), i] <- 1
        geno.out[which(geno.in[, i] == "a"), i] <- 2
        segr.type.out[i] <- 6
      }, D1.13 = {
        geno.out[which(geno.in[, i] == "a"), i] <- 1
        geno.out[which(geno.in[, i] == "o"), i] <- 2
        segr.type.out[i] <- 6
      }, D2.14 = {
        geno.out[which(geno.in[, i] == "ac"), i] <- 1
        geno.out[which(geno.in[, i] == "bc"), i] <- 2
        segr.type.out[i] <- 7
      }, D2.15 = {
        geno.out[which(geno.in[, i] == "a"), i] <- 1
        geno.out[which(geno.in[, i] == "ab"), i] <- 2
        segr.type.out[i] <- 7
      }, D2.16 = {
        geno.out[which(geno.in[, i] == "a"), i] <- 1
        geno.out[which(geno.in[, i] == "b"), i] <- 2
        segr.type.out[i] <- 7
      }, D2.17 = {
        geno.out[which(geno.in[, i] == "ab"), i] <- 1
        geno.out[which(geno.in[, i] == "a"), i] <- 2
        segr.type.out[i] <- 7
      }, D2.18 = {
        geno.out[which(geno.in[, i] == "a"), i] <- 1
        geno.out[which(geno.in[, i] == "o"), i] <- 2
        segr.type.out[i] <- 7
      })
    }
    else if (cross == "f2") {
      switch(EXPR = segr.type.in[i], A.H.B = {
        geno.out[which(geno.in[, i] == "a"), i] <- 1
        geno.out[which(geno.in[, i] == "ab"), i] <- 2
        geno.out[which(geno.in[, i] == "b"), i] <- 3
        segr.type.out[i] <- 1
      }, D.B = {
        geno.out[which(geno.in[, i] == "b"), i] <- 3
        geno.out[which(geno.in[, i] == "d"), i] <- 4
        segr.type.out[i] <- 2
      }, C.A = {
        geno.out[which(geno.in[, i] == "a"), i] <- 1
        geno.out[which(geno.in[, i] == "c"), i] <- 5
        segr.type.out[i] <- 3
      })
    }
    else if (cross == "backcross") {
      geno.out[which(geno.in[, i] == "a"), i] <- 1
      geno.out[which(geno.in[, i] == "ab"), i] <- 2
      segr.type.out[i] <- NA
    }
    else if (cross == "riself" || cross == "risib") {
      geno.out[which(geno.in[, i] == "a"), i] <- 1
      geno.out[which(geno.in[, i] == "b"), i] <- 3
      segr.type.out[i] <- NA
    }
    if (any(is.na(geno.out[, i]))) 
      stop(paste("Invalid marker codification. Please check data for marker", 
                 colnames(geno.in)[i]), ".", sep = "")
  }
  dimnames(geno.out) <- dimnames(geno.in)
  return(list(geno.out, segr.type.out))
}
<bytecode: 0x0000000007f581d8>
  <environment: namespace:onemap>