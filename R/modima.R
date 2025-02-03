MODIMAstat <- function(exposure, mediator, response){
  EM = bcdcor(exposure, mediator)
  MRE = pdcor(mediator, response, exposure)
  MODIMAstat <- EM*MRE
  attr(MODIMAstat, "names") <- "MODIMA stat"
  return(MODIMAstat)
}

subset.distance <- function(d, sub) as.dist(as.matrix(d)[sub, sub])

permuteDist <- function(d){
  d=as.matrix(d)
  n=nrow(d)
  p=sample(n)
  as.dist(d[p,p])
}

modima <- function(exposure, mediator, response, nrep=999){
  method <- "WARNING:Specify the number of replicates nrep>0 to perform the test"
  if (!is.null(nrep)){
    nrep <- floor(nrep)
    if (nrep <1)
      nrep <- 0
    if (nrep > 0)
      method <- paste("MODIMA: a Method for Multivariate Omnibus Distance Mediation Analysis")
  }
  else {
    nrep <- 0
  }
  modima_stat <- MODIMAstat(exposure, mediator, response)
  p1 <- energy::dcor.test(exposure, mediator, R = nrep)$p.value
  p2 <- energy::dcor.test(mediator, response, R = nrep)$p.value
  p3 <- spdcov.test(mediator, response, exposure, R = nrep)$p.value
  p4 <- spdcov.test(response, mediator, exposure, R = nrep)$p.value
  p5 <- energy::pdcov.test(mediator, response, exposure, R = nrep)$p.value
  bcdcorEM <- energy::bcdcor(exposure, mediator)
  bcdcorER <- energy::bcdcor(exposure, response)
  bcdcorMR <- energy::bcdcor(mediator, response)
  pdcorMRE <- energy::pdcor(mediator, response, exposure)
  attr(bcdcorEM, "names") <- paste0("exposure", "\u2013", "mediator", " bcdcor")
  attr(bcdcorER, "names") <- paste0("exposure", "\u2013", "response", " bcdcor")
  attr(bcdcorMR, "names") <- paste0("mediator", "\u2013", "response", " bcdcor")
  attr(pdcorMRE, "names") <- paste0("mediator", "\u2013", "response", "\u2013", "exposure", " pdcor")
  e <- list(method = method,
            data.name = base::paste("exposure = ", deparse(substitute(exposure)),
                                    "\n       mediator = ", deparse(substitute(mediator)),
                                    "\n       response = ", deparse(substitute(response)),
                                    "\nnumber of permutations=", nrep+1,
                                    "\nsample estimates are \n \t-bias-corrected distance correlation (energy::bcdcor) of indicated pairs and \n \t-partial distance correlation (energy::pdcor) of indicated triple"
                                    #deparse(substitute(exposure)), "and", deparse(substitute(mediator)), "removing", deparse(substitute(response))
                                    ),
            statistic = modima_stat, 
            p.value = max(p1, p2, p3, p4, p5),
            estimates = c(bcdcorEM, bcdcorER, bcdcorMR, pdcorMRE)
            #other options for getAnywhere(htest): parameter, alternative, null.value, conf.int
            )
  class(e) <- "htest"
  return(e)
  }

spdcov.test <- function (x, y, z, R){
  if (!(class(x) == "dist")) 
    x <- dist(x)
  if (!(class(y) == "dist")) 
    y <- dist(y)
  if (!(class(z) == "dist")) 
    z <- dist(z)
  Dx <- as.matrix(x)
  Dy <- as.matrix(y)
  Dz <- as.matrix(z)
  n <- nrow(Dx)
  Pxz <- energy:::projection(Dx, Dz)  ##we use ::: since these objects are not exported, otherwise we get this error: Error in { : task 1 failed - "'projection' is not an exported object from 'namespace:energy'"
  Py <- energy:::U_center(Dy)
  
  teststat <- n * energy:::U_product(Pxz, Py)
  den <- sqrt(energy:::U_product(Pxz, Pxz) * energy:::U_product(Py, Py))
  if (den > 0) {
    estimate <- teststat/(n * den)
  }
  else estimate <- 0
  bootfn <- function(Pxz, i, Py) {
    energy:::U_product(Pxz[i, i], Py)
  }
  reps <- replicate(R, expr = {
    i <- sample(1:n)
    bootfn(Pxz, i, Py = Py)
  })
  replicates <- n * reps
  pval <- (1 + sum(replicates >= teststat))/(1 + R)
  dataname <- paste("replicates ", R, sep = "")
  names(estimate) <- "spdcor"
  names(teststat) <- "n V^*"
  e <- list(call = match.call(), method = paste("semi-pdcov test", 
                                                sep = ""), statistic = teststat, estimate = estimate, 
            p.value = pval, n = n, replicates = replicates, data.name = dataname)
  class(e) <- "htest"
  return(e)
}
