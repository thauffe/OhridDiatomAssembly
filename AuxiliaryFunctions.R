# Function for sampling a species proportional to abundancies
propAbundResamp <- function(CommMat, Index){
  Species <- colnames(CommMat) # Get all species
  # Shuffle order of species according to their abundancies in the previous time slice
  Species <- sample(Species, # Species names to be shuffled
                    length(Species), # Number of draws equal to the number of species
                    # Abundancies from the previous time slice (i.e. index + 1)
                    # Add + 1 to the abundancies for not having a probability of 0 to draw a species
                    prob = CommMat[Index + 1, ] + 1) 
  SpeciesResampled <- rep(NA_character_, length(Species)) # Initialize empty vector
  A <- CommMat[Index, ] + 1 # Abundancies from the previous time slice
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
  Z <- (Obs - mean(Shuf)) / sd(Shuf) # Calculated standardized effect size
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
                  # We need two communities to use this function (i.e. a matrix)
                  # A vector or 1-dimensional matrix for our focal community causes an error
                  # (probably because of some matrix operations)
                  # Hence we create a matrix of our focal community in the first row and a fake community in the second.
                  # We only extract the functional richness of the former and dump the latter.
                  a = rbind(Comm[i, ], # Focal community
                            rep(400, 17)), # Fake community
                  m = 2,
                  w.abun = TRUE, # Weight by species' abundancies
                  ord = "podani",
                  corr = "cailliez", # Make non-euclidean distances euclidean using the cailliez correction
                  calc.FRic = TRUE, # Calculate functional richness
                  # Do not calculate other indices of functional diversity to safe some time
                  stand.FRic = FALSE,
                  calc.CWM = FALSE,
                  calc.FDiv = FALSE,
                  messages = FALSE) # Print no messages or progress on screen
    return(FdTmp$FRic[1])
  }
  # Pairwise distance matrix between species according to their differences in traits.
  # We use the gower distance because many traits are categorical/discrete traits (e.g. samll, medium, large)
  TraitsDist <- daisy(TraitsScaled, metric = "gower")
  # Initialize empty vector to store Functional richness values
  SesFric <- rep(NA_real_, nrow(Comm))
  registerDoParallel(Ncores)
  # For each community starting with the 2nd oldest 
  # (the oldest does not have a community before to do the resampling of species)
  for (i in 2:nrow(Comm)) { 
    # Observed Functional richness of the focal community (i.e. at one moment in time)
    FdObs <- dbFD(x = TraitsDist,
                  # We need two communities to use this function (i.e. a matrix)
                  # A vector or 1-dimensional matrix for our focal community causes an error
                  # (probably because of some matrix operations)
                  # Hence we create a matrix of our focal community in the first row and a fake community in the second.
                  # We only extract the functional richness of the former and dump the latter.
                  a = rbind(Comm[i, ], # Focal community
                            rep(400, 17)), # Fake community
                  m = 2,
                  w.abun = TRUE, # Weight by species' abundancies
                  ord = "podani",
                  corr = "cailliez", # Make non-euclidean distances euclidean using the cailliez correction
                  calc.FRic = TRUE, # Calculate functional richness
                  # Do not calculate other indices of functional diversity to safe some time
                  stand.FRic = FALSE, 
                  calc.CWM = FALSE,
                  calc.FDiv = FALSE,
                  messages = FALSE) # Print no messages or progress on screen
    # Get Runs-times (e.g. 999x) a functional richness for a community assembled according to the Null model
    # We use parallel computation on Ncores to speed this up
    FdNull <- foreach(iter = 1:Runs, 
                      .combine = c, # combine results in vector
                      .packages = c("FD", "cluster"), # Packages needed
                      # There is no specific order of the Null communities
                      .inorder = FALSE) %dopar%  sesFricOneComm(i, Comm, TraitsScaled) # Get functional richness for the i-th sample in our community matrix
    # Calculate standardized effect size of functional richness
    # (Observed FR - mean FR of the Null model) / standard deviation of the Null model
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

# Function to plot glacial/interglacial cycles as grey rectangles
plotLR04 <- function (LR04, Age) {
  # Do for all rows in the table with ages when cylce begin and end: 
  for(i in 1:nrow(LR04)){ 
    # If end if earlier than min age of our data (1.36 million years ago) and if it is a glacial period:
    if( -LR04[i,2] > min(Age) & (LR04[i,1] %% 2 == 1) ){
      # Draw rectangle
      rect(xleft = -LR04[i,2],  # Left of the rectangle equal to the begin of the cycle
           ybottom = par("usr")[3], # Bottom of the rectangle equal to the height of the plot
           xright = -LR04[i-1,2], # Right of the rectangle equal to the end of the cycle
           ytop = par("usr")[4], # Top of the rectangle equal to the height of the plot
           border = NA, # No outline of the rectangle
           col = adjustcolor("grey", alpha.f = 0.1)) # Grey rectangle with 10% opacity
    }
  }
}

# Fork of pgirmess::correlog, which is made for 2D spatial analyses.
# Time is like a 1D space but the function fails because it calculate pairwise distance between all samples
# to get the maximum distance between two samples.
# For time, getting the maximum distance is much easier because 
# it is the difference between the first and the last sample.
correlogFork <- function (time, z, method = "Moran", nbclass = NULL, breaks = NULL, ...) {
  coords <- cbind(0, abs(time))
  if (is.null(breaks)) {
    if (is.null(nbclass)) {
      l <- as.numeric(nrow(coords))
      n <- (l * l - l) / 2
      nbclass <- ceiling(log2(n) + 1)
    }
    etendue <- c( min(abs(diff(coords[, 2]))), 
                  # Here is the distance between the first and the last sample
                  abs(diff(coords[c(1, nrow(coords)), 2]))  )
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

# Get weights for Gaussian process smoothing and the regression analysis of MPD 
getWeight <- function(Time) {
  TimeDiff <- diff(Time) # Time difference between samples
  Weight <- rep(NA_real_, length(Time)) # Initialize empty vector
  for (i in 2:length(Time)) {
    # Weight for a focal sample is from half the distance towards the previous sample
    # until half the distance until the next sample
    Weight[i] <- mean(TimeDiff[(i-1):i])
  }
  # First and last sample do not have a weight because there is no younger or earlier sample. 
  # Thus, they are NA. We simply replace the NAs by the mean weight of the other 379 samples.
  Weight[is.na(Weight)] <- mean(Weight, na.rm = TRUE)
  # Until here 'weights' are a time difference. We divided them by the mean weights.
  # Thus, if all samples would have the same distance between them, all weights would be exactly 1.
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

# According to the brms manual, brms uses cbind internally
# and we need an alias for the multivariate grainsize smoothening
bind <- function(...) cbind(...)


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


# Function to get the rolling mean and standard error of Potassium 
getRollingMeanSE <- function (x) {
  M <- mean(x, na.rm = TRUE)
  S <- sd(x, na.rm = TRUE) / sqrt(length(na.omit(x)))
  return( c(M, S))
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

# Function to get a cumulative difference (Similar to cumsum).
# We need this because ages of shifts are parameterized in the function of earliest shift minus 2nd minus 3rd etc.
cumDiff <- function(x) {
  cumsum(x * c(1, rep(-1, length(x) - 1)))
}


setInitsRegr <- function (Seed, Length) {
  InitList <- vector(mode = "list", length = Length)
  for (i in 1:Length)  {
    set.seed(Seed[i])
    # Ages for shifts are difficult to initialize because we cannot just draw ages from a uniform distribution from present until origin of Lake Ohrid
    # because ages of the shifts are constrained to be consecutive/serial (i.e. ordered in time).
    # i.e. youngest shift is parameterized as brms::nlf(o1 ~ inv_logit(S6 - S5 - S4 - S3 - S2 - S1) * 3.242474)
    # Therefore, we first define sequence of 1000 ages (i.e. a uniform distribution) 
    # between the origin of Lake Ohrid and the midpoint between the earliest shift and 
    # 2nd oldest shift (i.e. 520 kilo years ago) identified in the change point analysis.
    # The sequence will be logit-transformed because brms uses an inverse-logit transformation to model the shifts.
    # (the inv_logit is like a step function setting returning 0 for all ages before the shift and 1 for after.
    # 0 and 1 will be multiplied with the effect of the shift \beta_{S_{k}})
    PotS6 <- -log((1/seq(520/1364.7, 1355/1364.7, length.out = 1000)) - 1) 
    # Uniform sequence from midpoint between shifts 4/5 to the midpoint between shifts 5/6
    # e.g. (Age shift 4 + Age shift 5) / 2 = 292
    PotS5 <- -log((1/seq(292/1364.7, 520/1364.7, length.out = 1000)) - 1)
    PotS4 <- -log((1/seq(155/1364.7, 292/1364.7, length.out = 1000)) - 1)
    PotS3 <- -log((1/seq(87/1364.7, 155/1364.7, length.out = 1000)) - 1)
    PotS2 <- -log((1/seq(25/1364.7, 87/1364.7, length.out = 1000)) - 1)
    PotS1 <- -log((1/seq(10/1364.7, 25/1364.7, length.out = 1000)) - 1)
    # Sample one oft ages from the sequence of 1000
    RanS6 <- sample(PotS6, 1)
    RanS5 <- sample(PotS5, 1)
    RanS4 <- sample(PotS4, 1)
    RanS3 <- sample(PotS3, 1)
    RanS2 <- sample(PotS2, 1)
    RanS1 <- sample(PotS1, 1)
    # Thus, we need to take the age of the earliest shift minus the age of the next shift 
    NewRanS5 <- RanS6 - RanS5
    NewRanS4 <- RanS5 - RanS4
    NewRanS3 <- RanS4 - RanS3
    NewRanS2 <- RanS3 - RanS2
    NewRanS1 <- RanS2 - RanS1
    # Check whether initial ages are ordered in time
    # brms:::inv_logit(cumDiff(c(RanS6, NewRanS5, NewRanS4, NewRanS3, NewRanS2, NewRanS1)))*1364.7
    # One list for each chain with initial values for all model parameters
    InitList[[i]] <- list(
      # Effect of the environmental and climatic predictors,
      # uniform distribution specifying a relatively small initial effect [-0.1, 0.1]
      "b_bGrain" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bIns" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bK" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bOaks" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bToc" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bIso" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bMed" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bAge" = array(runif(n = 1, min = -0.1, max = 0.1)),
      # Effect of the shift, uniform distribution specifying a relatively small effect [-0.1, 0.1]
      "b_bS1" = array(runif(n = 1, min = -0.1, max = 0.1)), 
      "b_bS2" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bS3" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bS4" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bS5" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bS6" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_bS7" = array(runif(n = 1, min = -0.1, max = 0.1)),
      "b_S1" = array(NewRanS1),
      "b_S2" = array(NewRanS2),
      "b_S3" = array(NewRanS3),
      "b_S4" = array(NewRanS4),
      "b_S5" = array(NewRanS5),
      "b_S6" = array(RanS6),
      # Initial value for the model error epsilon, must be positive
      "sigma" = runif(n = 1, min = 0.1, max = 0.3)
    )
  }
  return(InitList)
}

# Get conditional effects of one predictors while the remaining predictors are set to their mean
# This function is forked from brms because the regression analysis of MPD contains two time the predictor age.
# First: 'Age' used for the shifts, where it is in the range [0, 3.242474] because it needs to be positive. 
# (3.24 is scale(Age) + min(scale(Age)) * -1)
# Second: 'AgeScaled', which is the same than 'Age' but centered to 0 and to variance of 1. 
# That way it is treated the same way as all paleoenvironmental and climatic predictors.
# brms:::conditional_effects will set on of these two age predictors to its mean 
# and the resulting plot will not show the response of MPD to the combined effect of shifts and the linear age relationship
# (either the linear effect or the shifts only)
# We modified the brms function and added the argument 'same_as_effect' where we can specify that AgeScaled is the same than age.
# Then the function varies both from their minimum to maximum value to get the combined effect.
conditional_effects_forked <- function(x, effects = NULL, same_as_effects = NULL,
                                       conditions = NULL, 
                                       int_conditions = NULL, re_formula = NA, 
                                       prob = 0.95, robust = TRUE, 
                                       method = "posterior_epred",
                                       spaghetti = FALSE, surface = FALSE,
                                       categorical = FALSE, ordinal = FALSE,
                                       transform = NULL, resolution = 100, 
                                       select_points = 0, too_far = 0,
                                       probs = NULL, ...) {
  probs <- brms:::validate_ci_bounds(prob, probs = probs)
  method <- brms:::validate_pp_method(method)
  spaghetti <- brms:::as_one_logical(spaghetti)
  surface <- brms:::as_one_logical(surface)
  categorical <- brms:::as_one_logical(categorical)
  ordinal <- brms:::as_one_logical(ordinal)
  brms:::contains_draws(x)
  x <- restructure(x)
  new_formula <- brms:::update_re_terms(x$formula, re_formula = re_formula)
  bterms <- brmsterms(new_formula)
  
  if (!is.null(transform) && method != "posterior_predict") {
    stop2("'transform' is only allowed if 'method = posterior_predict'.")
  }
  if (ordinal) {
    warning2("Argument 'ordinal' is deprecated. ", 
             "Please use 'categorical' instead.")
  }
  rsv_vars <- brms:::rsv_vars(bterms)
  use_def_effects <- is.null(effects)
  if (use_def_effects) {
    effects <- brms:::get_all_effects(bterms, rsv_vars = rsv_vars)
  } else {
    # allow to define interactions in any order
    effects <- strsplit(as.character(effects), split = ":")
    if (any(unique(unlist(effects)) %in% rsv_vars)) {
      stop2("Variables ", collapse_comma(rsv_vars),
            " should not be used as effects for this model")
    }
    if (any(lengths(effects) > 2L)) {
      stop2("To display interactions of order higher than 2 ",
            "please use the 'conditions' argument.")
    }
    all_effects <- brms:::get_all_effects(
      bterms, rsv_vars = rsv_vars, comb_all = TRUE
    )
    ae_coll <- all_effects[lengths(all_effects) == 1L]
    ae_coll <- brms:::ulapply(ae_coll, paste, collapse = ":")
    matches <- match(lapply(all_effects, sort), lapply(effects, sort), 0L)
    if (sum(matches) > 0 && sum(matches > 0) < length(effects)) {
      invalid <- effects[setdiff(seq_along(effects), sort(matches))]  
      invalid <- brms:::ulapply(invalid, paste, collapse = ":")
      warning2(
        "Some specified effects are invalid for this model: ",
        collapse_comma(invalid), "\nValid effects are ", 
        "(combinations of): ", collapse_comma(ae_coll)
      )
    }
    effects <- unique(effects[sort(matches)])
    if (!length(effects)) {
      stop2(
        "All specified effects are invalid for this model.\n", 
        "Valid effects are (combinations of): ", 
        collapse_comma(ae_coll)
      )
    }
  }
  if (categorical || ordinal) {
    int_effs <- lengths(effects) == 2L
    if (any(int_effs)) {
      effects <- effects[!int_effs]
      warning2(
        "Interactions cannot be plotted directly if 'categorical' ", 
        "is TRUE. Please use argument 'conditions' instead."
      )
    }
  }
  if (!length(effects)) {
    stop2("No valid effects detected.")
  }
  mf <- model.frame(x)
  conditions <- brms:::prepare_conditions(
    x, conditions = conditions, effects = effects, 
    re_formula = re_formula, rsv_vars = rsv_vars
  )
  int_conditions <- lapply(int_conditions, 
                           function(x) if (is.numeric(x)) sort(x, TRUE) else x)
  int_vars <- brms:::get_int_vars(bterms)
  group_vars <- brms:::get_group_vars(bterms)
  out <- list()
  for (i in seq_along(effects)) {
    eff <- effects[[i]]
    cond_data <- brms:::prepare_cond_data(
      mf[, eff, drop = FALSE], conditions = conditions, 
      int_conditions = int_conditions, int_vars = int_vars,
      group_vars = group_vars, surface = surface, 
      resolution = resolution, reorder = use_def_effects
    )
    # Add the scaled time to the conditional data
    if (!is.null(same_as_effects)) {
      cond_data[, same_as_effects] <- seq(max(mf[, same_as_effects]),
                                          min(mf[, same_as_effects]),
                                          length.out = resolution) 
    }
    if (surface && length(eff) == 2L && too_far > 0) {
      # exclude prediction grid points too far from data
      ex_too_far <- mgcv::exclude.too.far(
        g1 = cond_data[[eff[1]]], 
        g2 = cond_data[[eff[2]]], 
        d1 = mf[, eff[1]],
        d2 = mf[, eff[2]],
        dist = too_far)
      cond_data <- cond_data[!ex_too_far, ]  
    }
    out[[i]] <- conditional_effects(
      bterms, fit = x, cond_data = cond_data, method = method, 
      surface = surface, spaghetti = spaghetti, categorical = categorical, 
      ordinal = ordinal, re_formula = re_formula, transform = transform, 
      conditions = conditions, int_conditions = int_conditions, 
      select_points = select_points, probs = probs, robust = robust,
      ...
    )[[1]]
  }
  names(out) <- names(effects)
  structure(out, class = "brms_conditional_effects")
}


# Function for caterpillar plot
caterpillar <- function (Draws, Prob = 0.95, LwdProb = 1,
                         Xlab = "", Ylab = "", Yticks = "",
                         Xlim = NULL,
                         ColLines = "black", 
                         AddDensity = TRUE, ColDensity = "grey",
                         AddVertZero = TRUE, Xaxt = NULL) {
  NcolDraws <- ncol(Draws) # How many columns (i.e. "predictors") in Draws table?
  Prob <- sort(Prob, decreasing = TRUE) # Sort probabilities descending to have lower Prob on top of higher Prob
  LwdProb <- sort(LwdProb, decreasing = FALSE) # Lower Prob thinner line
  # Get credible interval for all predictors
  # Number of list-entries equals how many probabilities we specify
  CredInt <- list()
  for (i in 1:length(Prob)) {
    # One matrix with credible interval in columns and lower/upper values in rows
    CredInt[[i]] <- matrix(NA_real_, nrow = 2, ncol = NcolDraws)
    for (y in 1:ncol(Draws)) {
      P <- (1 - Prob[i]) / 2 # Lower and upper percentile
      CredInt[[i]][, y] <- quantile(Draws[, y], c(P, 1 - P)) # Get 95% range of the predictors
    }
  }
  # Mean of the parameter (i.e. of the columns in the Draws table)
  Means <- colMeans(Draws)
  # Get limits for the x-axis
  RangeDensity <- c(0, 0) # Set range of the density distribution to zeros in case we do not use densities
  AddY <- 0 # Some padding for the y-axis of the plot
  AddTicksY <- 0 # Shift of tick-labels on the y-axis
  if (AddDensity) { # Range of the density distribution
    AddY <- 1
    AddTicksY <- 0.5
    # Get density distribution for all predictors
    DensList <- sapply(1:NcolDraws, function(x) density(Draws[, x],
                                                        from = min(Draws[, x]),
                                                        to = max(Draws[, x])),
                       simplify = FALSE)
    RangeDensity <- range(unlist(lapply(DensList, function(x) x$x)))
  }
  RangeCredInt <- range(c(unlist(CredInt), RangeDensity)) # Combine ranges of credible intervals and densities
  if (!is.null(Xlim)) { # Allow manual setting of limits to the x-axis
    RangeCredInt <- Xlim
  }
  # Plot 
  if (length(ColLines) == 1 && NcolDraws > 1) { # Same color for all lines?
    ColLines <- rep(ColLines, NcolDraws) # Repeat color for all predictors
  }
  if (AddDensity && length(ColDensity) == 1 && NcolDraws > 1) { # Same color for all density polygones?
    ColDensity <- rep(ColDensity, NcolDraws) # Repeat color for all predictors
  }
  if ( length(Yticks) == 1) { # Are there are no custom labels for the y-axis?
    Yticks <- colnames(Draws) # Names of the predictors according to the Draws table
  }
  plot(1, 1, type = "n", xlim = c(RangeCredInt[1], RangeCredInt[2]), ylim = c(1, NcolDraws + AddY),
       xlab = Xlab, ylab = Ylab, yaxt = "n", xaxt = Xaxt) # Create empty canvas
  axis(side = 2, at = 1:NcolDraws + AddTicksY, labels = Yticks) # Add labels at y-axis ticks
  for (i in 1:NcolDraws) { # Do for all columns in Draws table:
    if (AddDensity) { # If we want to plot the density distribution:
      X <- DensList[[i]]$x # Get all values for the x-axis
      Y <- DensList[[i]]$y # Get all values for the y-axis
      Y <- (Y / max(Y)) * 0.8 # Scale height of y to 80% row-height in the plot
      # Polygon showing the density distribution of the sampled parameter
      polygon(x = c(X, rev(X)),
              y = i + c(Y, rep(0, length(Y))), 
              col = ColDensity[i], # Color of the density distribution
              border = NA) # No outline for the polygon
    }
    for (p in 1:length(Prob)) {
      # Lines with different thickness for the different Probs
      lines(x = CredInt[[p]][, i], # Take Credible interval of the i-th predictor for the p-th probability 
            y = rep(i, 2), # Position on y-axis (1 for the first predicto, 2 for the 2nd etc.)
            lwd = LwdProb[p], # Thickness of the lines
            col = ColLines[i]) # Color of the lines
    }
    # Add a dot at the mean of the predictor parameter
    points(Means[i], i, pch = 19, col = ColLines[i], cex = max(LwdProb) * 0.75)
  }
  if (AddVertZero) {
    abline(v = 0, lty = 2) # Dashed line at 0 indicating no effect
  }
}


# Function to plot posterior versus prior draws
plotPriorVerusPosterior <- function (Prior, Post,
                                 DensFrom, DensTo, DensN,
                                 Xlim, Ylim, 
                                 Xlab, 
                                 ColPrior, YTicks, PriorTickLabels, PriorMulti,
                                 Xticks, XticksLabels,
                                 LabAxis4X, LabAxis4Y,
                                 TitleBarY, TitleX, TitleY, Title) {
  DenPrior <- density(Prior, from = DensFrom, to = DensTo, n = DensN) # Get prior density
  DenPost <- density(Post) # Get posterior density
  plot(DenPost$x, DenPost$y, type = "l", ylim = Ylim, xlim = Xlim,
       xlab = Xlab, ylab = "Posterior density", xaxt = "n", yaxt = "n") # Plot posterior density
  axis(side = 2, at = YTicks) # Add ticks to the left axis
  axis(side = 4, at = YTicks, label = PriorTickLabels,
       col.ticks = ColPrior, col.axis = ColPrior) # Add ticks to the right axis
  axis(side = 1, at = Xticks, labels = XticksLabels) # Add tick-labels to the x-axis
  text(x = par("usr")[2] + LabAxis4X, y = LabAxis4Y, labels = "Prior density", col = ColPrior,
       xpd = TRUE, srt = -90) # Add label to the axis on the right
  lines(DenPrior$x, DenPrior$y * PriorMulti, col = ColPrior) # Draw line for prior density
  rect(xleft = par("usr")[1], ybottom = TitleBarY, xright = par("usr")[2], ytop = par("usr")[4],
       col = "grey70", border = "black") # Draw grey rectangle as background for the plot title
  text(x = TitleX, y = TitleY, labels = Title, adj = c(0.5, 0.5)) # Plot title
}



# Function to check if moved periods are plus/minus 25% of their original duration.
movedPeriodsAccepted <- function (DataRegr, PerMoved) {
  OK <- TRUE # Will get false if duration criteria is not fulfilled
  Per <- levels(DataRegr$Shift) # Which periods do we have?
  for (i in 1:length(unique(Per))) { # Do for all periods
    OrigDuration <- diff(range(DataRegr[DataRegr$Shift == as.numeric(Per[i]), "Age"])) # Original duration
    MovedDuration <- diff(range(PerMoved[PerMoved$Shift == as.numeric(Per[i]), "Age"])) # Duration of moved periods
    # Moved periods are plus/minus 25% of their original duration?
    if ( (MovedDuration * 1.25) < OrigDuration || (MovedDuration * 0.75) > OrigDuration ) {
      OK <- FALSE
    }
  }
  return(OK)
}


movePeriods <- function (DataRegr) {
  PerMoved <- DataRegr
  Per <- levels(DataRegr$Shift) # Which periods do we have?
  MovesOK <- FALSE
  while (!MovesOK) {
    MovedPeriods <- rep(NA, nrow(DataRegr))
    # Do for all periods
    # We start with the most recent period because it is the shortest one
    Start <- nrow(DataRegr)
    for (i in length(unique(Per)):1) { 
      W <- which(DataRegr$Shift == as.numeric(Per[i]))
      MaxChange <- round(length(W) * 0.25)
      Change <- sample(-MaxChange:MaxChange, 1)
      End <- min(W) + Change
      # Force earliest period to date back until 1.36 Ma
      if (i == 1) {
        End <- 1
      }
      MovedPeriods[Start:End] <- as.numeric(Per[i]) 
      Start <- End - 1
    }
    PerMoved$Shift <- as.factor(MovedPeriods)
    if (length(levels(PerMoved$Shift)) == 7) { # Do we have 7 periods / 6 shifts ?
      # Check if moved periods are plus/minus 25% of their original duration
      MovesOK <- movedPeriodsAccepted(DataRegr, PerMoved)
    }

  }
  return(PerMoved)
}


