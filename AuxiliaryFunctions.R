propAbundResamp <- function(CommMat, Index){
  Species <- colnames(CommMat)
  Species <- sample(Species,
                    length(Species),
                    prob = CommMat[Index + 1, ] + 1) # Start with a random species!
  SpeciesResampled <- rep(NA_character_, length(Species))
  A <- CommMat[Index, ] + 1
  Target <- length(Species)
  Counter <- 1
  while(Counter <= Target){
    DrawSpecies <- sample(Species, 1, prob = A)
    SpeciesResampled[Counter] <- DrawSpecies
    W <- which(Species == DrawSpecies)
    Species <- Species[-W]
    A <- A[-W]
    Counter <- Counter + 1
  }
  return(SpeciesResampled)
}


sesMpdOneComm <- function(Samp, Dis, NullModel, AbundanceWeighted, Runs, 
                          CommMat = NULL, Index = NULL){
  Obs <- mpd(Samp, Dis, AbundanceWeighted) # mntd
  Shuf <- rep(0, Runs)
  if(NullModel == "phylo.labels"){
    for(i in 1:Runs){
      S <- sample(colnames(Dis))
      DisTmp <- Dis
      colnames(DisTmp) <- S
      rownames(DisTmp) <- S
      Shuf[i] <- mpd(Samp, DisTmp, AbundanceWeighted) # mntd
    }
  }
  # Does not differ much from phylo.labels:
  if(NullModel == "comm.labels"){
    for(i in 1:Runs){
      S <- sample(colnames(Samp))
      colnames(Samp) <- S
      Shuf[i] <- mpd(Samp, DisTmp, AbundanceWeighted) # mntd 
    }
  }
  # Label reshuffling proportional to abundancies from the previous community
  if(NullModel == "phylo.labels.propabun"){
    for(i in 1:Runs){
      S <- propAbundResamp(CommMat, Index)
      DisTmp <- Dis
      colnames(DisTmp) <- S
      rownames(DisTmp) <- S
      Shuf[i] <- mpd(Samp, DisTmp, AbundanceWeighted) # mntd
    }
  }
  Z <- (Obs - mean(Shuf)) / sd(Shuf)
  return(Z)
}


sesMPD <- function (Comm, CophTree, Runs, i) {
  SesMpd <- NA
  if (sum(Comm[i, -c(1:3)]) != 0) {
    TmpCophTree <- CophTree - 2*Comm$Age[i]/1000
    diag(TmpCophTree) <- 0
    SesMpd <- sesMpdOneComm(Samp = Comm[i, -c(1:3)],
                            Dis = TmpCophTree, 
                            NullModel = "phylo.labels.propabun",
                            AbundanceWeighted = TRUE,
                            Runs = Runs,
                            CommMat = Comm[, -c(1:3)],
                            Index = i - 1)
  }
  return(SesMpd)
}


getSesMpd <- function (Comm, CophTree, Runs, Ncores) {
  registerDoParallel(Ncores)
  SesMpd <- foreach(iter = 2:nrow(Comm),
                    .combine = c,
                    .packages = c("picante"),
                    .inorder = TRUE) %dopar%  sesMPD(Comm, CophTree, Runs, iter)
  stopImplicitCluster()
  SesMpd <- c(NA, SesMpd)
  return(SesMpd)
}



getSesFric <- function (Comm, TraitsScaled, Runs, Ncores) {
  
  sesFricOneComm <- function(i, Comm, TraitsScaled) {
    S <- propAbundResamp(CommMat = Comm, Index = i - 1)
    TraitsScaledTmp <- TraitsScaled[S, ]
    rownames(TraitsScaledTmp) <- colnames(Comm)
    TraitsDistTmp <- daisy(TraitsScaledTmp, metric = "gower")
    FdTmp <- dbFD(x = TraitsDistTmp,
                  a = rbind(Comm[i, ], rep(400, 17)),
                  m = 2,
                  w.abun = TRUE,
                  ord = "podani",
                  corr = "cailliez",
                  calc.FRic = TRUE,
                  stand.FRic = FALSE,
                  calc.CWM = FALSE,
                  calc.FDiv = FALSE,
                  messages = FALSE)
    return(FdTmp$FRic[1])
  }
  
  TraitsDist <- daisy(TraitsScaled, metric = "gower")
  SesFric <- rep(NA_real_, nrow(Comm))
  registerDoParallel(Ncores)
  for (i in 2:nrow(Comm)) {
    FdObs <- dbFD(x = TraitsDist,
                  a = rbind(Comm[i, ], rep(400, 17)),
                  m = 2,
                  w.abun = TRUE,
                  ord = "podani",
                  corr = "cailliez",
                  calc.FRic = TRUE,
                  stand.FRic = FALSE,
                  calc.CWM = FALSE,
                  calc.FDiv = FALSE,
                  messages = FALSE)
    FdNull <- foreach(iter = 1:Runs,
                      .combine = c,
                      .packages = c("FD", "cluster"),
                      .inorder = FALSE) %dopar%  sesFricOneComm(i, Comm, TraitsScaled)
    SesFric[i] <- (FdObs$FRic[1] - mean(FdNull, na.rm = TRUE)) / sd(FdNull, na.rm = TRUE)
  }
  stopImplicitCluster()
  return(SesFric)
}


stochasticRarefaction <- function (EntireComm, Size = 1, Reps = 100,
                                   Planktonics = NULL, EndThrough = FALSE, Ncores = 1) {
  StochCommList <- vector(mode = "list", length = Reps)
  
  stochasticComm <- function (EntireComm, Size = 1) {
    CommRan <- as.data.frame(matrix(0, nrow = nrow(EntireComm), ncol = ncol(EntireComm)))
    colnames(CommRan) <- colnames(EntireComm)
    for (i in 1:nrow(EntireComm)) {
      Obs <- EntireComm[i, ]
      Present <- Obs > 0
      Species <- names(Obs)[Present]
      Obs <- Obs[Present]
      R <- rep(Species, times = Obs)
      SamInd <- sample(R, size = Size, replace = TRUE)
      Tab <- table(SamInd)
      CommRan[i, names(Tab)] <- c(Tab)
    }
    return(CommRan)
  }
  
  parStochComm <- function (EntireComm, Size) {
    StochComm <- stochasticComm (EntireComm, Size)
    if (!is.null(Planktonics)) {
      StochComm <- StochComm[, which(colnames(StochComm) %in% Planktonics)]
    }
    if (EndThrough) {
      N <- nrow(StochComm)
      if (any(StochComm$Cyclotella_bifacialis > 0)) {
        Tmp <- rep(0, N)
        Tmp[min(which(StochComm$Cyclotella_bifacialis > 0)):N] <- 1
        StochComm$Cyclotella_bifacialis <- ifelse(Tmp >= StochComm$Cyclotella_bifacialis,
                                                  Tmp,
                                                  StochComm$Cyclotella_bifacialis)
      }
      if (any(StochComm$Cyclotella_fottii > 0)) {
        Tmp <- rep(0, N)
        Tmp[min(which(StochComm$Cyclotella_fottii > 0)):N] <- 1
        StochComm$Cyclotella_fottii <- ifelse(Tmp >= StochComm$Cyclotella_fottii,
                                              Tmp,
                                              StochComm$Cyclotella_fottii)
      }
      if (any(StochComm$Lindavia_thienemannii > 0)) {
        Tmp <- rep(0, N)
        Tmp[min(which(StochComm$Lindavia_thienemannii > 0)):N] <- 1
        StochComm$Lindavia_thienemannii <- ifelse(Tmp >= StochComm$Lindavia_thienemannii,
                                                  Tmp,
                                                  StochComm$Lindavia_thienemannii)
      }
      # Endemic species but today absent i.e. endemic fossils
      # Those should include an abundance of 1 from first until last appearance datum
      EndAbs <- c("Cribrionella_ohridana", "Cyclotella_cavitata", "Cyclotella_sollevata",
                  "Pantocsekiella_multiocellata", "Pantocsekiella_preocellata")
      for(y in 1:length(EndAbs)){ # Ugly loop
        Tmp <- StochComm[, EndAbs[y]]
        if (sum(Tmp > 0) > 1) {
          R <- range(which(Tmp > 0))
          Tmp[ which(Tmp[R[1]:R[2]]==0) + R[1]-1 ] <- 1
          StochComm[, EndAbs[y]] <- Tmp
        }
      }
    }
    return(StochComm)
  }
  
  registerDoParallel(Ncores)
  StochCommList <- foreach(iter = 1:Reps,
                           .inorder = FALSE) %dopar%  parStochComm(EntireComm, Size)
  stopImplicitCluster()
  SrList <- lapply(StochCommList, function(x) rowSums(ifelse(x > 0, 1, 0)))
  SrMat <- do.call("cbind", SrList)
  SrMean <- rowMeans(SrMat)
  SrCI <- apply(SrMat, 1, function(x) quantile(x, c(0.025, 0.975)))
  SrRarefied <- rbind(SrMean, SrCI)
  rownames(SrRarefied) <- c("Sr", "LwrCI", "UprCI")
  Out <- list()
  Out[[1]] <- SrRarefied
  Out[[2]] <- StochCommList
  return(Out)
}


plotLR04 <- function (LR04, Age) {
  for(i in 1:nrow(LR04)){
    if( -LR04[i,2] > min(Age) & (LR04[i,1] %% 2 == 1) ){
      rect(xleft = -LR04[i,2], ybottom = par("usr")[3],
           xright = -LR04[i-1,2], ytop = par("usr")[4],
           border = NA, col = adjustcolor("grey", alpha.f = 0.1))
    }
  }
}

# Fork of pgirmess::correlog, which fails with high temporal resolution due to matrix size
correlogFork <- function (time, z, method = "Moran", nbclass = NULL, breaks = NULL, ...) {
  coords <- cbind(0, abs(time))
  if (is.null(breaks)) {
    if (is.null(nbclass)) {
      l <- as.numeric(nrow(coords))
      n <- (l * l - l) / 2
      nbclass <- ceiling(log2(n) + 1)
    }
    etendue <- c( min(abs(diff(coords[, 2]))), abs(diff(coords[c(1, nrow(coords)), 2]))  )
    breaks1 <- seq(etendue[1], etendue[2], l = nbclass + 1)
    breaks2 <- breaks1 + 1e-06
    breaks <- cbind(breaks1[1:length(breaks1) - 1], breaks2[2:length(breaks2)])
    breaks[1, 1] <- breaks[1, 1] - 1e-06
  }
  # nbreaks <- nrow(breaks)
  nbclass <- nrow(breaks)
  lst.nb1 <- rep(list(NA), nbclass)
  lst.z1 <- rep(list(NA), nbclass)
  for (i in 1:nbclass) {
    lst.z1[[i]] <- z
    lst.nb1[[i]] <- dnearneigh(coords, breaks[i, 1], breaks[i, 2])
    zero <- which(card(lst.nb1[[i]]) == 0)
    if (length(zero) > 0) {
      lst.nb1[[i]] <- dnearneigh(coords[-zero, ], breaks[i, 1], breaks[i, 2])
      lst.z1[[i]] <- z[-zero]
    }
  }
  lst.res1 <- rep(list(NA), nbclass)
  for (i in 1:nbclass) {
    xt <- switch(pmatch(method, c("Moran", "Geary"), nomatch = 3),
                 try(moran.test(lst.z1[[i]], nb2listw(lst.nb1[[i]], style = "W"), ...), silent = TRUE),
                 try(geary.test(lst.z1[[i]], nb2listw(lst.nb1[[i]], style = "W"), ...), silent = TRUE),
                 stop("Method must be 'Moran' or 'Geary'"))
    if (inherits(xt, "try-error")) {
      stop("Bad selection of class breaks, try another one...")
    }
    else {
      x <- xt$estimate[1]
      p <- xt$p.value
      N <- sum(card(lst.nb1[[i]]))
    }
    lst.res1[[i]] <- c(x = x, p = p, N = N)
  }
  meth <- names(xt[[3]][1])
  mat <- matrix(unlist(lst.res1), ncol = 3, byrow = TRUE)
  res <- cbind(dist.class = rowMeans(breaks),
               coef = mat[, 1], p.value = mat[, 2], n = mat[, 3])
  attributes(res) <- c(attributes(res), list(Method = meth))
  class(res) <- c("correlog", "matrix")
  res
}


getWeight <- function(Time) {
  TimeDiff <- diff(Time)
  Weight <- rep(NA_real_, length(Time))
  for (i in 2:length(Time)) {
    Weight[i] <- mean(TimeDiff[(i-1):i])
  }
  Weight[is.na(Weight)] <- mean(Weight, na.rm = TRUE)
  Weight <- Weight / mean(Weight)
  return(Weight)
}


setInitsPollen <- function(Seed, Length, K) {
  InitList <- vector(mode = "list", length = Length)
  for (i in 1:Length) {
    set.seed(Seed[i])
    Inter <- runif(n = 1, min = 0, max = 0.5)
    InitList[[i]] <- list(
      "Intercept" = Inter,
      "sdgp_1" = array(runif(n = 1, min = 4, max = 6)),
      "lscale_1" = matrix(runif(n = 1, min = 0.0001, max = 0.01), 1, 1),
      "zgp_1" = array(runif(n = K, min = -1, max = 1)),
      "b_Intercept" = Inter
    )
  }
  return(InitList)
}


bind <- function(...) cbind(...) # brms uses cbind internally


setInitsIsotope <- function(Seed, Length, K) {
  InitList <- vector(mode = "list", length = Length)
  for (i in 1:Length) {
    set.seed(Seed[i])
    Inter <- runif(n = 1, min = -6, max = -4)
    InitList[[i]] <- list(
      "Intercept" = Inter,
      "sdgp_1" = array(runif(n = 1, min = 0, max = 2)),
      "lscale_1" = matrix(runif(n = 1, min = 0.0001, max = 0.01), 1, 1),
      "zgp_1" = array(runif(n = K, min = -1, max = 1)),
      "sigma" = runif(n = 1, min = 0, max = 1),
      "b_Intercept" = Inter
    )
  }
  return(InitList)
}


setInitsTic <- function(Seed, Length, K) {
  InitList <- vector(mode = "list", length = Length)
  for (i in 1:Length) {
    set.seed(Seed[i])
    InitList[[i]] <- list(
      "Intercept" = runif(n = 1, min = -5, max = -3),
      "sdgp_1" = array(runif(n = 1, min = 0.0001, max = 0.1)),
      "lscale_1" = matrix(runif(n = 1, min = 0.0001, max = 0.01), 1, 1),
      "zgp_1" = array(runif(n = K, min = -1, max = 1)),
      "phi" = runif(n = 1, min = 0.0001, max = 0.1),
      "b_Intercept" = runif(n = 1, min = -5, max = -3)
    )
  }
  return(InitList)
}


setInitsToc <- function(Seed, Length, K) {
  InitList <- vector(mode = "list", length = Length)
  for (i in 1:Length) {
    set.seed(Seed[i])
    Inter <- runif(n = 1, min = 0.01, max = 2)
    InitList[[i]] <- list(
      "Intercept" = Inter,
      "sdgp_1" = array(runif(n = 1, min = 0.0001, max = 0.1)),
      "lscale_1" = matrix(runif(n = 1, min = 0.0001, max = 0.01), 1, 1),
      "zgp_1" = array(runif(n = K, min = -1, max = 1)),
      "phi" = runif(n = 1, min = 0.0001, max = 0.1),
      "b_Intercept" = Inter
    )
  }
  return(InitList)
}


setInitsK <- function(Seed, Length, K) {
  InitList <- vector(mode = "list", length = Length)
  for (i in 1:Length) {
    set.seed(Seed[i])
    Inter <- runif(n = 1, min = 5, max = 6)
    InitList[[i]] <- list(
      "Intercept" = Inter,
      "sdgp_1" = array(runif(n = 1, min = 1.5, max = 2)),
      "lscale_1" = matrix(runif(n = 1, min = 0.0001, max = 0.001), 1, 1),
      "zgp_1" = array(runif(n = K, min = -1, max = 1)),
      "sigma" = runif(n = 1, min = 0.3, max = 0.5),
      "b_Intercept" = Inter
    )
  }
  return(InitList)
}


setInitsQuartz <- function(Seed, Length, K) {
  InitList <- vector(mode = "list", length = Length)
  for (i in 1:Length) {
    set.seed(Seed[i])
    Inter <- runif(n = 1, min = -2, max = -1)
    InitList[[i]] <- list(
      "Intercept" = Inter,
      "sdgp_1" = array(runif(n = 1, min = 0.4, max = 0.6)),
      "lscale_1" = matrix(runif(n = 1, min = 0.0001, max = 0.01), 1, 1),
      "zgp_1" = array(runif(n = K, min = -1, max = 1)),
      "phi" = runif(n = 1, min = 130, max = 150),
      "b_Intercept" = Inter
    )
  }
  return(InitList)
}


setInitsRegr <- function (Seed, Length) {
  InitList <- vector(mode = "list", length = Length)
  for (i in 1:Length) {
    set.seed(Seed[i])
    Inter <- runif(n = 1, min = -6, max = -4)
    InitList[[i]] <- list(
      "b_b0" = array(runif(n = 1, min = -0.69, max = 0.71)),
      "b_b1" = array(runif(n = 1, min = 1.29, max = 1.31)),
      "b_b2" = array(runif(n = 1, min = -1.21, max = -1.19)),
      "b_b3" = array(runif(n = 1, min = 1.89, max = 1.91)),
      "b_b4" = array(runif(n = 1, min = -1.71, max = -1.69)),
      "b_b5" = array(runif(n = 1, min = 0.59, max = 0.61)),
      "b_b6" = array(runif(n = 1, min = -0.81, max = 0.79)),
      "b_bAge" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bGrain" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bIns" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bIso" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bK" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bMed" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bOaks" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bToc" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_S1" = array(runif(n = 1, min = -4.51, max = -4.49)),
      "b_S2" = array(runif(n = 1, min = 0.64, max = 0.66)),
      "b_S3" = array(runif(n = 1, min = 1.69, max = 1.71)),
      "b_S4" = array(runif(n = 1, min = 0.16, max = 0.17)),
      "b_S5" = array(runif(n = 1, min = 1.16, max = 1.17)),
      "b_S6" = array(runif(n = 1,  min = 0.63, max = 0.65)),
      "sigma" = runif(n = 1, min = 0.1, max = 0.3)
    )
  }
  return(InitList)
}
