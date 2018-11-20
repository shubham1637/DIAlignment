library(devtools)
devtools::install_github("Roestlab/DIAlign")
library(DIAlign)

##########
OuterProdMeanNormAll6Func <- function (data, pep, runA, runB) 
{
  # pep = "55711_TFISPIK/2"
  # runA = "Plasma5"; runB = "Plasma18";
  # data <- StrepChromsPlasma
  num_of_frag <- length(data[[runA]][[pep]])
  num_of_samplesA <- length(data[[runA]][[pep]][[1]][, 1])
  num_of_samplesB <- length(data[[runB]][[pep]][[1]][, 1])
  MeanNormA <- sapply(data[[runA]][[pep]], function(x) sum(x[, 
                                                             2])/num_of_samplesA)
  MeanNormA <- mean(MeanNormA)
  MeanNormB <- sapply(data[[runB]][[pep]], function(x) sum(x[, 
                                                             2])/num_of_samplesB)
  MeanNormB <- mean(MeanNormB)
  outerProdList <- list()
  for (i in 1:num_of_frag) {
    print(i)
    NormIntensityA <- data[[runA]][[pep]][[i]][, 2]/MeanNormA
    NormIntensityB <- data[[runB]][[pep]][[i]][, 2]/MeanNormB
    outerProdList[[i]] <- outer(NormIntensityA, NormIntensityB)
  }
  return(outerProdList)
}

# when merging single run with already merged smth: C + ABD
OuterProdMeanNormAll6Func_forMerged <- function (data, dataMerged, pep, runA, runB) 
  {
    num_of_frag <- length(data[[runA]][[pep]])
    num_of_samplesA <- length(data[[runA]][[pep]][[1]][, 1])
    num_of_samplesB <- length(dataMerged[[runB]][[pep]][[1]][, 1])
    MeanNormA <- sapply(data[[runA]][[pep]], function(x) sum(x[, 
                                                               2])/num_of_samplesA)
    MeanNormA <- mean(MeanNormA)
    MeanNormB <- sapply(dataMerged[[runB]][[pep]], function(x) sum(x[, 
                                                               2])/num_of_samplesB)
    MeanNormB <- mean(MeanNormB)
    outerProdList <- list()
    for (i in 1:num_of_frag) {
      NormIntensityA <- data[[runA]][[pep]][[i]][, 2]/MeanNormA
      NormIntensityB <- dataMerged[[runB]][[pep]][[i]][, 2]/MeanNormB
      outerProdList[[i]] <- outer(NormIntensityA, NormIntensityB)
    }
    return(outerProdList)
}

getSimilarityMatrix_forMerged <- function (data, dataMerged, pep, runA, runB, type = c("dotProductMasked", 
                                                                           "dotProduct"), dotProdThresh = 0.96, cosAngleThresh = 0.3) 
{
  type <- match.arg(type)
  switch(type, dotProduct = {
    OuterProdNormAll6 <- OuterProdMeanNormAll6Func(data, 
                                                   pep, runA, runB)
    s <- add(OuterProdNormAll6)
  }, dotProductMasked = {
    OuterProdNormAll6 <- OuterProdMeanNormAll6Func_forMerged(data, dataMerged, 
                                                   pep, runA, runB)
    s1 <- add(OuterProdNormAll6)
    OuterProdL2NormAll <- OuterProdL2NormAllFunc_forMerged(data,dataMerged, pep, 
                                                 runA, runB)
    s2 <- cos(2 * acos(pmin(add(OuterProdL2NormAll), 1)))
    MASK <- (s1 > quantile(s1, dotProdThresh))
    AngleGreat <- (((1 * MASK) * s2) + (1 - MASK)) > cosAngleThresh
    s <- s1 * (1 * AngleGreat)
  })
  return(s)
}

OuterProdL2NormAllFunc_forMerged <- function(data,dataMerged, pep, runA, runB){
  num_of_frag <- length(dataMerged[[runB]][[pep]])
  L2NormA <- sapply(data[[runA]][[pep]], function(x) x[,2])
  L2NormA <- sqrt(rowSums(L2NormA^2))
  L2NormB <- sapply(dataMerged[[runB]][[pep]], function(x) x[,2])
  L2NormB <- sqrt(rowSums(L2NormB^2))
  outerProdList <- list()
  for (i in 1:num_of_frag){
    NormIntensityA <- data[[runA]][[pep]][[i]][,2]/L2NormA
    NormIntensityA[is.nan(NormIntensityA)] <-0
    NormIntensityB <- dataMerged[[runB]][[pep]][[i]][,2]/L2NormB
    NormIntensityB[is.nan(NormIntensityB)] <-0
    outerProdList[[i]] <- outer(NormIntensityA, NormIntensityB)
  }
  return(outerProdList) }





obj2  <- getAlignObj(s, go, FreeEndGaps=F)
obj3 <- getAffineAlignObj(s, go, ge, FreeEndGaps= F)
# MaxScores  <- getAlignment(obj3)$TrB[,"score"]
for (matrix in getAlignment(obj3)) {
  MaxScores <- matrix[, "score"]
}

### Function TerminalNAnumber, gives as an output vector of 2 values: number of NAs in the beginning 
### and number of NAs in the end of the input vector
TerminalNAnumber <- function (vector) {
  NAstart = 0 
  NAend = 0
  count = 0
  for (i in 1:length(vector)) {
    if (NAstart < count) {
      break
    }
    if (is.na(vector[i])) {
      NAstart <- NAstart + 1
    }
    count <- count + 1
  }
  count1 = 0
  for (i in length(vector):1) {
    if (NAend < count1) {
      break
    }
    if (is.na(vector[i])) {
        NAend <- NAend + 1
        }
    count1 <- count1 + 1
   }
    
  return(c(NAstart, NAend))
}

### Function SumNearMaxScores gives as an output vector of values: background sums (4x4) for each max score
SumNearMaxScores <- function(alignObj) {
    MaxScores <- getAlignment(alignObj)[[1]][,"score"]
    RowIndices <- getAlignment(alignObj)[[1]][,"indexA_aligned"]
    ColIndices <- getAlignment(alignObj)[[1]][,"indexB_aligned"]
  
  # alignObj@M - score matrix
  NumRow <- nrow(alignObj@M)
  NumCol <- ncol(alignObj@M)

  # first maximum score - M[n,n], calculate the background sum of the 4x4 square with this score in the lowest right corner
  MaxIndices <- c(NumRow, NumCol)
  SumNearMax <- c(sum(alignObj@M[(NumRow - 3):NumRow, (NumCol - 3):NumCol]))
  # check if there are shifts in the alignment (NAs in the beginning/end) and do not consider these terminal NAs further
  # get start and end rows (after/before terminal NAs)
  if (TerminalNAnumber(RowIndices)[1] > TerminalNAnumber(ColIndices)[1]) {
    Start <- 1 + TerminalNAnumber(RowIndices)[1]
  }
  else {
    Start <- 1 + TerminalNAnumber(ColIndices)[1]
  }
  if (TerminalNAnumber(RowIndices)[2] > TerminalNAnumber(ColIndices)[2]) {
    End <- length(RowIndices) - TerminalNAnumber(RowIndices)[2]
  }
  else {
    End <- length(ColIndices) - TerminalNAnumber(ColIndices)[2]
  }
  
  if (TerminalNAnumber(RowIndices[1:End])[2] > 0) {
    End <- End - TerminalNAnumber(RowIndices[1:End])[2]
  }
  if (TerminalNAnumber(ColIndices[1:End])[2] > 0) {
    End <- End - TerminalNAnumber(ColIndices[1:End])[2]
  }
  
  print(c("End ", End,"Start ", Start))
  
  if (!is.na(RowIndices[End])) {
    LastRowNotNA <- RowIndices[End]
  }
  else {
    LastRowNotNA <- RowIndices[End-1]
  }
  if (!is.na(ColIndices[End])) {
    LastColNotNA <- ColIndices[End]
  }
  else {
    LastColNotNA <- ColIndices[End-1]
  }
  LastColNotNA <- ColIndices[End]
  
  
  # firstly, substitute NAs in the middle with previous indices
  for (i in (End - 1):Start) {
    # start with (End -1) because we've already calculated sum for the first maximum score M[n,n]
    if (is.na(RowIndices[i])) {
      if (!is.na(RowIndices[i-1])) {
        RowIndices[i] <- RowIndices[i-1]
        LastRowNotNA <-  RowIndices[i-1]
      }
      else {
        RowIndices[i] <- LastRowNotNA
      }
    }
    if (is.na(ColIndices[i])) {
      if (!is.na(ColIndices[i-1])) {
        ColIndices[i] <-  ColIndices[i-1]
        LastColNotNA <-  ColIndices[i-1]
      }
      else {
        ColIndices[i] <- LastColNotNA
      }
    }
      
    if (alignObj@M[RowIndices[i] + 1,ColIndices[i] + 1] == (MaxScores[i])) {
      MaxIndices <- c(RowIndices[i] + 1,ColIndices[i] + 1)
      if (((MaxIndices[1] - 3) > 0) & ((MaxIndices[2] - 3) > 0)) {
      SumNearMax <- c(SumNearMax, sum(alignObj@M[(MaxIndices[1] - 3): MaxIndices[1], (MaxIndices[2] - 3): MaxIndices[2]]))
      }
    }
    #else {
      #print(c(alignObj@M[RowIndices[i] + 1,ColIndices[i] + 1] , MaxScores[i]))
      #print(c(RowIndices[i] + 1,ColIndices[i] + 1))
   # }
  }

  return(SumNearMax)
}

SumNearMaxScores(Alignobj)
getAlignment(Alignobj)

############################################################################
######################## PLOTS CHROM ALIGNMENT #############################
############################################################################
########### function plotSingleAlignedChrom - plots aligned chromatograms
plotChromatogram <- function(data, runname, peptide, ObservedRT, printTitle = TRUE){
  # run <- filenames[runname]
  df <- do.call("cbind", data[[runname]][[peptide]])
  df <- df[,!duplicated(colnames(df))]
  df <- melt(df, id.vars="time", value.name = "Intensity")
  g <- ggplot(df, aes(time, Intensity, col=variable)) + geom_line(show.legend = FALSE) + theme_bw()
  if(printTitle) g <- g + ggtitle(paste0(runname, ", ",peptide)) + theme(plot.title = element_text(size = 9, hjust = 0.5))
  g <- g + geom_vline(xintercept=ObservedRT[peptide, runname], lty="dotted", size = 0.4)
  return(g)
}

plotSingleAlignedChrom <- function(data, peptide, runname, filenames, idx, t, printTitle = TRUE){
  # data <- StrepChroms
  # run <- runname
  # peptide = "9091_NDYGNTTLALR/2"
  # idx <- AlignedIndices[[1]][,"indexA_aligned"]
  # run <- filenames[runname]
  intensity <- list()
  for(k in 1:length(data[[runname]][[peptide]])){
    mutateInt <- data[[runname]][[peptide]][[k]][idx, 2]
    mutateInt <- na.locf(na.locf(mutateInt, na.rm = FALSE),fromLast = TRUE)
    intensity[[k]] <- mutateInt
  }
  df <- do.call("cbind", intensity)
  Index <- 1:nrow(df)
  df <- cbind(Index, as.data.frame(df))
  df <- melt(df, id.vars="Index", value.name = "Intensity")
  
  g <- ggplot(df, aes(Index, Intensity, col=variable)) + geom_line(show.legend = FALSE) + theme_bw()
  if(printTitle) g <- g + ggtitle(paste0(run, ", ",peptide)) + theme(plot.title = element_text(hjust = 0.5)) 
  return(g)
}

plotAlignedChroms <- function(data, pair, peptide, ObservedRT, #AlignErrorinSec,
                              AlignedIndex, tA, tB, FourOrTwo = TRUE){
  runA <- strsplit(pair, split = "_")[[1]][1]
  runB <- strsplit(pair, split = "_")[[1]][2]
  pTL <- plotChromatogram(data, runA, peptide, ObservedRT, FALSE)
  pBL <- plotChromatogram(data, runB, peptide, ObservedRT, FALSE)
  #pBL <- pBL + geom_vline(xintercept=ObservedRT[peptide, runB]+AlignErrorinSec[peptide, pair], lty="dashed", size = 0.4, color = "blue")
  if(FourOrTwo){
    pTR <- plotSingleAlignedChrom(data, peptide, runA, filenames, AlignedIndex[[1]][, "indexA_aligned"], tA, FALSE) + geom_vline(xintercept=which.min(abs(tA - ObservedRT[peptide, runA])), lty="dotted", size = 0.4)
    pBR <- plotSingleAlignedChrom(data, peptide, runB, filenames, AlignedIndex[[1]][, "indexB_aligned"], tB, FALSE) + geom_vline(xintercept=which.min(abs(tB - ObservedRT[peptide, runB])), lty="dotted", size = 0.4)
    pBR <- pBR + geom_vline(xintercept=which.min(abs(tA - ObservedRT[peptide, runA])), lty="dashed", size = 0.4, color = "red")
    p <- grid.arrange(pTL, pTR, pBL, pBR, nrow=2, ncol=2, top = paste0(pair, ", ",peptide))
  }
  else {
    pBL <- pBL + geom_vline(xintercept=tB[which.min(abs(tA - ObservedRT[peptide, runA]))], lty="dashed", size = 0.4, color = "red")
    grid.arrange(pTL, pBL, nrow=2, ncol=1, top = paste0(pair, ", ",peptide))
  }
}

# This function plots chromatogram after merging runs (doesn't consider actual retention time, uses only indices)
plotChromatogramMerged <- function(data, runname, peptide, #type = c("Mean", "Median", "WeightedMean"),
                                   printTitle = TRUE){
  # type <- match.arg(type)
  # switch(type, Mean = {
  for (i in 1:length(data[[runname]][[peptide]])) {
     # data[[runname]][[peptide]][[i]] <- data[[runname]][[peptide]][[i]][, c(1:2)]
      names(data[[runname]][[peptide]][[i]]) <- c("time", paste(i,peptide,sep = "_"))
   } 
  #   } , Median = {
  #   for (i in 1:length(data[[runname]][[peptide]])) {
  #     data[[runname]][[peptide]][[i]] <- data[[runname]][[peptide]][[i]][, c(1,3)]
  #     names(data[[runname]][[peptide]][[i]]) <- c("time", paste(i,peptide,sep = "_")) 
  #   }
  #   } ,
  #  WeightedMean = {
  #   for (i in 1:length(data[[runname]][[peptide]])) {
  #     data[[runname]][[peptide]][[i]] <- data[[runname]][[peptide]][[i]][, c(1,4)]
  #     names(data[[runname]][[peptide]][[i]]) <- c("time", paste(i,peptide,sep = "_"))
  #   }
  # } 
  # )
  
  df <- do.call("cbind", data[[runname]][[peptide]])
  df <- df[,!duplicated(colnames(df))]
  df <- melt(df, id.vars="time", value.name = "Intensity")
  g <- ggplot(df, aes(time, Intensity, col=variable)) + geom_line(show.legend = FALSE) + theme_bw()
  if(printTitle) g <- g + ggtitle(paste0(runname, ", ", peptide)) + #,subtitle = type) 
                                   theme(plot.title = element_text(size = 9,hjust = 0.5)) #, plot.subtitle = element_text(size=8, hjust=0.5, face="italic"))
  # g <- g + geom_vline(xintercept=ObservedRT[peptide, runname], lty="dotted", size = 0.4)
  return(g)
}
  
for (i in 1:20) {
  print(sum(MaxScoresAllRandom[[i]], na.rm = T))
}
############################################################################
#################### END OF PLOTS CHROM ALIGNMENT ##########################
############################################################################

Score4runsMatrices <- list()
for (peptide in peptides) {
  Score4runsMatrices[[peptide]] <- matrix(nrow = length(runs), ncol = length(runs))
  colnames(Score4runsMatrices[[peptide]]) <- runs
  rownames(Score4runsMatrices[[peptide]]) <- runs
}


# Fit a local alignment between chromatogram groups
pair_names <- vector(); runs <- names(StrepChroms)
for (i in 1:(length(runs)-1)){
  for (j in (i+1): length(runs)){
    pair_names <- c(paste(runs[i], runs[j], sep = "_"), pair_names)
  }}

AlignedChromsPlots <- list()
for (peptide in peptides) {
  AlignedChromsPlots[[peptide]] <- list()
}

for(peptide in peptides){
  pdf(paste("AlignedChroms",gsub("/","-",peptide),"4runs.pdf", sep = "_"))
  for(pair in pair_names){
    print(pair)
    run_pair <- strsplit(pair, split = "_")[[1]]
    gapQuantile <- 0.5; goFactor <- 1/8; geFactor <- 40
    simMeasure <- "dotProductMasked"
    #run_pair <- c("run1", "run2")
    Err <- matrix(NA, nrow = length(peptides), ncol = 1)
    rownames(Err) <- peptides
    ListOfSimMatrices <- list()
    ListOfMaxScores <- list()
    ListOfBackgroundSum <- list()
    ListOfAlignedIndices <- list()
    NumberOfGaps <- c()
    NumberOfGapsList <- list()
    AverageMaxScore <- list()
    MedianOfM <- list()
  
    s <- getSimilarityMatrix(StrepChroms, peptide, run_pair[1], run_pair[2], type = simMeasure)
    gapPenalty <- getGapPenalty(s, gapQuantile, type = simMeasure)
    Alignobj <- getAffineAlignObj(s, go = gapPenalty*goFactor, ge = gapPenalty*geFactor)
    ListOfSimMatrices[[paste0("", peptide)]] <- Alignobj@M
    AlignedIndices <- getAlignment(Alignobj)
    
    NumberOfGaps <-c(NumberOfGaps, sum(is.na(AlignedIndices[[1]][TerminalNAnumber(AlignedIndices[[1]][,1])[1]:TerminalNAnumber(AlignedIndices[[1]][,1])[2],1]), 
                                       is.na(AlignedIndices[[1]][TerminalNAnumber(AlignedIndices[[1]][,1])[1]:TerminalNAnumber(AlignedIndices[[1]][,1])[2],2])))
    
    NumberOfGapsList[[peptide]] <- sum(is.na(AlignedIndices[[1]][TerminalNAnumber(AlignedIndices[[1]][,1])[1]:TerminalNAnumber(AlignedIndices[[1]][,1])[2],1]), 
                                       is.na(AlignedIndices[[1]][TerminalNAnumber(AlignedIndices[[1]][,1])[1]:TerminalNAnumber(AlignedIndices[[1]][,1])[2],2]))
  
    ListOfAlignedIndices[[paste0("", peptide)]] <- AlignedIndices
    AverageMaxScore[[peptide]] <- AlignedIndices[[1]][length(AlignedIndices[[1]][,3]), 3] / length(AlignedIndices[[1]][ ,3])
    MedianOfM[[peptide]] <- median(Alignobj@M)
    ListOfMaxScores[[paste0("", peptide)]] <- AlignedIndices[[1]][,"score"]
    ListOfBackgroundSum[[paste0("", peptide)]] <- SumNearMaxScores(Alignobj)
    tA <- StrepChroms[[run_pair[1]]][[peptide]][[1]][["time"]]
    tB <- StrepChroms[[run_pair[2]]][[peptide]][[1]][["time"]]
    tA.aligned <- mapIdxToTime(tA, AlignedIndices[[1]][,"indexA_aligned"])
    tB.aligned <- mapIdxToTime(tB, AlignedIndices[[1]][,"indexB_aligned"])
    predictTime <- tB.aligned[which.min(abs(tA.aligned - StrepAnnot[peptide, run_pair[1]]))]
    deltaT <- predictTime - StrepAnnot[peptide, run_pair[2]]
    Err[peptide, 1] <- deltaT
    
    plot(plotAlignedChroms(StrepChroms, pair, peptide, StrepAnnot, #AlignErrorinSec,
                                                              AlignedIndices, tA, tB, FourOrTwo = TRUE))
               
  }
  dev.off() 
  
  
  # from lists of cumulative scores create lists of separate max and background scores
  # and list of dataframes for each peptide (to create a plot)
  BackgroundSum <- list()
  MaxScoresAll <- list()
  Ratio <- list()
  plots <- list()
  plotsRatio <- list()
  plotStats <- list()
  ScoreDataFrames <- list()
  StatsDataFrame <- list()
  peptidenames <- names(ListOfBackgroundSum)
  for (peptide in peptidenames) {
    BackgroundSum[[paste0("", peptide)]] <- c(rep(NA, length(ListOfBackgroundSum[[peptide]])))
    for (i in length(BackgroundSum[[paste0("", peptide)]]):2) {
      BackgroundSum[[paste0("", peptide)]][i] <- ListOfBackgroundSum[[paste0("", peptide)]][i-1] -
                                                      ListOfBackgroundSum[[paste0("", peptide)]][i]
      if (is.na(BackgroundSum[[paste0("", peptide)]][i]) | (BackgroundSum[[paste0("", peptide)]][i] %in% Inf)) {
        BackgroundSum[[paste0("", peptide)]] <- BackgroundSum[[paste0("", peptide)]][-i]
      }
    }
    MaxScoresAll[[paste0("", peptide)]] <- c(rep(NA,length(ListOfMaxScores[[peptide]])))
    for (i in (length(MaxScoresAll[[paste0("", peptide)]]) -1):1) {
      MaxScoresAll[[paste0("", peptide)]][i] <- ListOfMaxScores[[peptide]][i+1] - ListOfMaxScores[[peptide]][i]
      if (is.na(MaxScoresAll[[paste0("", peptide)]][i]) | (MaxScoresAll[[paste0("", peptide)]][i] %in% Inf)) {
        MaxScoresAll[[paste0("", peptide)]] <- MaxScoresAll[[paste0("", peptide)]][-i]
      }
    }
    Ratio[[peptide]] <- MaxScoresAll[[peptide]][length(BackgroundSum[[peptide]]):1] / BackgroundSum[[peptide]]
    
    m <- ggplot(as.data.frame(Ratio[[peptide]]), aes(x = as.numeric(Ratio[[peptide]]))) 
    m <- m + geom_density() 
    p <- ggplot_build(m)
    y_max <- max(p$data[[1]]$y)
    AUC[[peptide]] <- auc(p$data[[1]]$x,p$data[[1]]$y)
    
    ScoreDataFrames[[peptide]] <- data.frame(rbind(cbind(BackgroundSum[[peptide]], c(rep("Background Sum", 
                                                                                         length(BackgroundSum[[peptide]])))),
                                                   cbind(MaxScoresAll[[peptide]], c(rep("Maximum Score",
                                                                                         length(MaxScoresAll[[peptide]]))))) )
    names(ScoreDataFrames[[peptide]]) <- c("Score", "Type")
    
    StatsDataFrame[[peptide]] <- data.frame(rbind(c(NumberOfGapsList[[peptide]], "Number of gaps"),
                                                  c(AverageMaxScore[[peptide]], "Average maximum score"),
                                                  c(MedianOfM[[peptide]], "Median of the score matrics"),
                                                  c(median(Ratio[[peptide]], na.rm = T), "Median of a MaxScore/Background ratio")))
    
    Score4runsMatrices[[peptide]][run_pair[1],run_pair[2]] <- (1 + NumberOfGapsList[[peptide]])/AUC[[peptide]]
    Score4runsMatrices[[peptide]][run_pair[2],run_pair[1]] <- (1 + NumberOfGapsList[[peptide]])/AUC[[peptide]]
  }
  
}


pdf(paste("Trees_4runs_GapsAUC.pdf"))
for(peptide in peptidenames) {
  DistMat = as.dist(Score4runsMatrices[[peptide]], diag = TRUE)
  Clustered <- hclust(DistMat, method = "single")
  my_tree <- as.phylo(Clustered) 
  plot(my_tree, main = peptide)
}
dev.off()

## save some plots 
require(MESS)
for (peptide in peptidenames) {
  plots[[peptide]] <- ggplot(ScoreDataFrames[[peptide]], aes(x = as.numeric(Score), colour = Type, fill = Type)) +
    geom_density(alpha = 0.3) +
    labs(x = "Score") +
    theme(legend.title = element_blank())
  
  m <- ggplot(as.data.frame(Ratio[[peptide]]), aes(x = as.numeric(Ratio[[peptide]]))) 
  m <- m + geom_density() 
  p <- ggplot_build(m)
  y_max <- max(p$data[[1]]$y)
  AUC[[peptide]] <- auc(p$data[[1]]$x,p$data[[1]]$y)
  
  plotsRatio[[peptide]] <- ggplot(as.data.frame(Ratio[[peptide]]), aes(x = as.numeric(Ratio[[peptide]]))) +   
    geom_density(alpha = 0.3, color = "darkblue", fill = "darkblue") +
    labs(x = "Maximum score at each step /
         Background sum (4x4)") +
    theme(axis.title.y = element_text(size = 7),
          axis.title.x = element_text(size = 6)) +
    annotate("text", label = paste("AUC = ", round(AUC[[peptide]],2)), x = max(Ratio[[peptide]], na.rm = T) - 
               max(Ratio[[peptide]], na.rm = T)/10, y = y_max - y_max/10, hjust = 1) +
    annotate("text", label = paste("mean = ", round(mean(Ratio[[peptide]],na.rm = T),2)),
             x = max(Ratio[[peptide]], na.rm = T) - max(Ratio[[peptide]], na.rm = T)/10, y = y_max - 2*y_max/10,hjust = 1) +
    annotate("text", label = paste("sd = ", round(sd(Ratio[[peptide]], na.rm = T),2)),
             x = max(Ratio[[peptide]], na.rm = T) - max(Ratio[[peptide]], na.rm = T)/10, y = y_max - 3*y_max/10,hjust = 1)
  
  plotStats[[peptide]] <-   ggplot(subset(StatsDataFrame[[peptide]], !(StatsDataFrame[[peptide]]$V2 %in% c("Median of a MaxScore/Background ratio"))),
           aes(x = V2, y = as.numeric(as.character(score)))) +
            geom_bar(stat = "identity", fill = "darkblue") +
      labs( x = "",
            y = "") +
      theme(axis.text.y = element_text(size = 7),
            axis.text.x = element_text(size = 6)) 
    
  
}

require(lattice)
statplots <- list()
for(peptide in peptidenames) {
  statplots[[peptide]] <- grid.arrange(plotsRatio[[peptide]], plotStats[[peptide]], ncol = 2)
}

pdf(paste("ChromAndDistribAndStats_run1-2_4x4background.pdf"))
for(peptide in peptidenames) {
  print(grid.arrange(plotChromatogram(StrepChroms, "run1", peptide,StrepAnnot, printTitle =TRUE ),
                     plotChromatogram(StrepChroms, "run2", peptide,StrepAnnot, printTitle =TRUE ),
                     statplots[[peptide]],
                     ncol =1))
}
dev.off()

# modify a few functions in order to run alignment of random peptides
### 1: add pepA and pepB instead of pep (everywhere)
OuterProdMeanNormAll6Func_random <- function(data, pepA, pepB, runA, runB) 
{
  num_of_frag <- length(data[[runA]][[pepA]])
  num_of_samplesA <- length(data[[runA]][[pepA]][[1]][, 1])
  num_of_samplesB <- length(data[[runB]][[pepB]][[1]][, 1])
  MeanNormA <- sapply(data[[runA]][[pepA]], function(x) sum(x[, 
                                                             2])/num_of_samplesA)
  MeanNormA <- mean(MeanNormA)
  MeanNormB <- sapply(data[[runB]][[pepB]], function(x) sum(x[, 
                                                             2])/num_of_samplesB)
  MeanNormB <- mean(MeanNormB)
  outerProdList <- list()
  for (i in 1:num_of_frag) {
    NormIntensityA <- data[[runA]][[pepA]][[i]][, 2]/MeanNormA
    NormIntensityB <- data[[runB]][[pepB]][[i]][, 2]/MeanNormB
    outerProdList[[i]] <- outer(NormIntensityA, NormIntensityB)
  }
  return(outerProdList)
}

### 2: add pepA and pepB instead of pep (everywhere)
OuterProdL2NormAllFunc_random <- function(data, pepA,pepB, runA, runB){
  num_of_frag <- length(data[[runA]][[pepA]])
  L2NormA <- sapply(data[[runA]][[pepA]], function(x) x[,2])
  L2NormA <- sqrt(rowSums(L2NormA^2))
  L2NormB <- sapply(data[[runB]][[pepB]], function(x) x[,2])
  L2NormB <- sqrt(rowSums(L2NormB^2))
  outerProdList <- list()
  for (i in 1:num_of_frag){
    NormIntensityA <- data[[runA]][[pepA]][[i]][,2]/L2NormA
    NormIntensityA[is.nan(NormIntensityA)] <-0
    NormIntensityB <- data[[runB]][[pepB]][[i]][,2]/L2NormB
    NormIntensityB[is.nan(NormIntensityB)] <-0
    outerProdList[[i]] <- outer(NormIntensityA, NormIntensityB)
  }
  return(outerProdList) }

### 3: add pepA and pepB instead of pep (everywhere)
library(FuzzyStatTra)
library(Omisc)
getSimilarityMatrixForRandom <- function (data, pepA, pepB, runA, runB, type = c("dotProductMasked", 
                                                                          "dotProduct", "cosineAngle", "cosine2Angle", "euclideanDist", 
                                                                          "covariance", "correlation"), dotProdThresh = 0.96, cosAngleThresh = 0.3) 
{
  type <- match.arg(type)
  switch(type, dotProduct = {
    OuterProdNormAll6 <- OuterProdMeanNormAll6Func(data, 
                                                   pep, runA, runB)
    s <- add(OuterProdNormAll6)
  }, cosineAngle = {
    OuterProdL2NormAll <- OuterProdL2NormAllFunc(data, pep, 
                                                 runA, runB)
    s <- add(OuterProdL2NormAll)
  }, cosine2Angle = {
    OuterProdL2NormAll <- OuterProdL2NormAllFunc(data, pep, 
                                                 runA, runB)
    s <- cos(2 * acos(pmin(add(OuterProdL2NormAll), 1)))
  }, dotProductMasked = {
    OuterProdNormAll6 <- OuterProdMeanNormAll6Func_random(data, 
                                                   pepA, pepB, runA, runB) #added pepA, pepB instead of pep and '_random' to the
                                                                           #name of function
    s1 <- add(OuterProdNormAll6)
    OuterProdL2NormAll <- OuterProdL2NormAllFunc_random(data, pepA,pepB, #added pepA, pepB instead of pep and '_random' to the
                                                        runA, runB)      #name of function
                                                 
    s2 <- cos(2 * acos(pmin(add(OuterProdL2NormAll), 1)))
    MASK <- (s1 > quantile(s1, dotProdThresh))
    AngleGreat <- (((1 * MASK) * s2) + (1 - MASK)) > cosAngleThresh
    s <- s1 * (1 * AngleGreat)
  }, euclideanDist = {
    OuterProdEucl <- OuterProdEuclFunc(data, pep, runA, runB)
    s <- 1/(1 + sqrt(add(OuterProdEucl)))
  }, covariance = {
    s <- OuterProdCovFunc(data, pep, runA, runB)
  }, correlation = {
    s <- OuterProdCorFunc(data, pep, runA, runB)
    s[is.na(s)] <- 0
  })
  return(s)
}

## Then repeat the local alignment and building plots with 2 random peptides from different runs. 
gapQuantile <- 0.5; goFactor <- 1/8; geFactor <- 40
simMeasure <- "dotProductMasked"
run_pair <- c("run2", "run3")
Err <- matrix(NA, nrow = length(PeptidesA), ncol = 1)
rownames(Err) <- PeptidesA
ListOfSimMatricesRandom <- list()
ListOfMaxScoresRandom <- list()
ListOfBackgroundSumRandom <- list()
ListOfAlignedIndicesRandom <- list()
NumberOfGapsRandom <- c()
NumberOfGapsRandomList <- list()
AverageMaxScoreRandom <- list()
MedianOfMRandom <- list()
PeptidesA <- peptides[8:12]
PeptidesB <- peptides[13:17]
for (peptideA in PeptidesA){
  for (peptideB in PeptidesB) {
    s <- getSimilarityMatrixForRandom(StrepChroms, peptideA, peptideB, run_pair[1], run_pair[2], type = simMeasure)
    gapPenalty <- getGapPenalty(s, gapQuantile, type = simMeasure)
    Alignobj <- getAffineAlignObj(s, go = gapPenalty*goFactor, ge = gapPenalty*geFactor)
    ListOfSimMatricesRandom[[paste0("", peptideA," - ",peptideB)]] <- Alignobj@M
    AlignedIndices <- getAlignment(Alignobj)
    ListOfAlignedIndicesRandom[[paste0("", peptideA," - ",peptideB)]] <- AlignedIndices
    ListOfMaxScoresRandom[[paste0("", peptideA," - ",peptideB)]] <- AlignedIndices[[1]][,"score"]
    ListOfBackgroundSumRandom[[paste0("", peptideA," - ",peptideB)]] <- SumNearMaxScores(Alignobj)
    
    NumberOfGapsRandom <-c(NumberOfGapsRandom, sum(is.na(AlignedIndices[[1]][TerminalNAnumber(AlignedIndices[[1]][,1])[1]:TerminalNAnumber(AlignedIndices[[1]][,1])[2],1]), 
                                                   is.na(AlignedIndices[[1]][TerminalNAnumber(AlignedIndices[[1]][,1])[1]:TerminalNAnumber(AlignedIndices[[1]][,1])[2],2]),
                                                   is.na(AlignedIndices[[1]][TerminalNAnumber(AlignedIndices[[1]][,2])[1]:TerminalNAnumber(AlignedIndices[[1]][,2])[2],1]), 
                                                   is.na(AlignedIndices[[1]][TerminalNAnumber(AlignedIndices[[1]][,2])[1]:TerminalNAnumber(AlignedIndices[[1]][,2])[2],2])))
    
    
    NumberOfGapsRandomList[[paste0("", peptideA," - ",peptideB)]] <- sum(is.na(AlignedIndices[[1]][TerminalNAnumber(AlignedIndices[[1]][,1])[1]:TerminalNAnumber(AlignedIndices[[1]][,1])[2],1]), 
                                                                         is.na(AlignedIndices[[1]][TerminalNAnumber(AlignedIndices[[1]][,1])[1]:TerminalNAnumber(AlignedIndices[[1]][,1])[2],2]),
                                                                         is.na(AlignedIndices[[1]][TerminalNAnumber(AlignedIndices[[1]][,2])[1]:TerminalNAnumber(AlignedIndices[[1]][,2])[2],1]), 
                                                                         is.na(AlignedIndices[[1]][TerminalNAnumber(AlignedIndices[[1]][,2])[1]:TerminalNAnumber(AlignedIndices[[1]][,2])[2],2]))

    AverageMaxScoreRandom[[paste0("", peptideA," - ",peptideB)]] <-
      AlignedIndices[[1]][length(AlignedIndices[[1]][,3]), 3] / length(AlignedIndices[[1]][ ,3])
    MedianOfMRandom[[paste0("", peptideA," - ",peptideB)]] <- median(Alignobj@M)
    tA <- StrepChroms[[run_pair[1]]][[peptideA]][[1]][["time"]]
    tB <- StrepChroms[[run_pair[2]]][[peptideB]][[1]][["time"]]
    tA.aligned <- mapIdxToTime(tA, AlignedIndices[[1]][,"indexA_aligned"])
    tB.aligned <- mapIdxToTime(tB, AlignedIndices[[1]][,"indexB_aligned"])
    predictTime <- tB.aligned[which.min(abs(tA.aligned - StrepAnnot[peptideA, run_pair[1]]))]
    deltaT <- predictTime - StrepAnnot[peptideB, run_pair[2]]
   # Err[peptideA, 1] <- deltaT
  }
}

##### Same procedures for Random peptides
# from lists of cumulative scores create lists of separate max and background scores
# and list of dataframes for each peptide (to create a plot)
BackgroundSumRandom <- list()
MaxScoresAllRandom <- list()
plotsRandom <- list()
ScoreDataFramesRandom <- list()
StatsDataFrameRandom <- list()
RatioRandom <- list()
plotsRatioRandom <- list()
plotStatsRandom <- list()
peptidenamesRandom <- names(ListOfBackgroundSumRandom)
for (peptide in peptidenamesRandom) {
  #pdf(paste("AlignedChroms",gsub("/","-",peptide),"runsRandom.pdf", sep = "_"))
  BackgroundSumRandom[[paste0("", peptide)]] <- c(rep(NA, length(ListOfBackgroundSumRandom[[peptide]])))
  for (i in length(BackgroundSumRandom[[paste0("", peptide)]]):2) {
    BackgroundSumRandom[[paste0("", peptide)]][i] <- ListOfBackgroundSumRandom[[paste0("", peptide)]][i-1] -
      ListOfBackgroundSumRandom[[paste0("", peptide)]][i]
    if (is.na(BackgroundSumRandom[[paste0("", peptide)]][i]) | (BackgroundSumRandom[[paste0("", peptide)]][i] %in% Inf)) {
      BackgroundSumRandom[[paste0("", peptide)]] <- BackgroundSumRandom[[paste0("", peptide)]][-i]
    }
  }
  MaxScoresAllRandom[[paste0("", peptide)]] <- c(rep(NA,length(ListOfMaxScoresRandom[[peptide]])))
  for (i in (length(MaxScoresAllRandom[[paste0("", peptide)]]) -1):1) {
    MaxScoresAllRandom[[paste0("", peptide)]][i] <- ListOfMaxScoresRandom[[peptide]][i+1] - ListOfMaxScoresRandom[[peptide]][i]
    if (is.na(MaxScoresAllRandom[[paste0("", peptide)]][i]) | (MaxScoresAllRandom[[paste0("", peptide)]][i] %in% Inf)) {
      MaxScoresAllRandom[[paste0("", peptide)]] <- MaxScoresAllRandom[[paste0("", peptide)]][-i]
    }
  }
  RatioRandom[[peptide]] <- MaxScoresAllRandom[[peptide]][length(BackgroundSumRandom[[peptide]]):1] / BackgroundSumRandom[[peptide]]
  
  m <- ggplot(as.data.frame(RatioRandom[[peptide]]), aes(x = as.numeric(RatioRandom[[peptide]]))) 
  m <- m + geom_density() 
  p <- ggplot_build(m)
  y_max <- max(p$data[[1]]$y)
  AUC[[peptide]] <- auc(p$data[[1]]$x,p$data[[1]]$y)
  print(AUC[[peptide]])


  ScoreDataFramesRandom[[peptide]] <- data.frame(rbind(cbind(BackgroundSumRandom[[peptide]], c(rep("Background Sum", 
                                                                                       length(BackgroundSumRandom[[peptide]])))),
                                                 cbind(MaxScoresAllRandom[[peptide]], c(rep("Maximum Score",
                                                                                      length(MaxScoresAllRandom[[peptide]]))))) )
  names(ScoreDataFramesRandom[[peptide]]) <- c("Score", "Type")
  
  StatsDataFrameRandom[[peptide]] <- data.frame(rbind(c(NumberOfGapsRandomList[[peptide]], "Number of gaps"),
                                                c(AverageMaxScoreRandom[[peptide]], "Average maximum score"),
                                                c(MedianOfMRandom[[peptide]], "Median of the score matrics"),
                                                c(median(RatioRandom[[peptide]], na.rm = T), "Median of a MaxScore/Background ratio")))
  names(StatsDataFrameRandom[[peptide]]) <- c("score", "V2")
  
  # store some plots
  plotsRatioRandom[[peptide]] <- ggplot(as.data.frame(RatioRandom[[peptide]]), aes(x = as.numeric(RatioRandom[[peptide]]))) +   
    geom_density(alpha = 0.3, color = "darkblue", fill = "darkblue") +
    labs(x = "Maximum score at each step /
         Background sum (4x4)") +
    theme(axis.title.y = element_text(size = 7),
          axis.title.x = element_text(size = 6)) +
    annotate("text", label = paste("AUC = ", round(AUC[[peptide]],2)), x = max(RatioRandom[[peptide]], na.rm = T) - 
               max(RatioRandom[[peptide]], na.rm = T)/10, y_max - y_max/10, hjust = 1) +
    annotate("text", label = paste("mean = ", round(mean(RatioRandom[[peptide]],na.rm = T),2)),
             x = max(RatioRandom[[peptide]], na.rm = T) - max(RatioRandom[[peptide]], na.rm = T)/10, y = y_max - 2*y_max/10,hjust = 1) +
    annotate("text", label = paste("median = ", round(median(RatioRandom[[peptide]], na.rm = T),2)),
             x = max(RatioRandom[[peptide]], na.rm = T) - max(RatioRandom[[peptide]], na.rm = T)/10, y = y_max - 3*y_max/10,hjust = 1) +
    annotate("text", label = paste("sd = ", round(sd(RatioRandom[[peptide]], na.rm = T),2)),
             x = max(RatioRandom[[peptide]], na.rm = T) - max(RatioRandom[[peptide]], na.rm = T)/10, y = y_max - 4*y_max/10,hjust = 1)
  
  plotStatsRandom[[peptide]] <-   ggplot(subset(StatsDataFrameRandom[[peptide]], !(StatsDataFrameRandom[[peptide]]$V2 %in% c("Median of a MaxScore/Background ratio"))),
                                         aes(x = V2, y = as.numeric(as.character(score)))) +
    geom_bar(stat = "identity", fill = "darkblue") +
    labs( x = "",
          y = "") +
    theme(axis.text.y = element_text(size = 7),
          axis.text.x = element_text(size = 6)) 
}
  
  # statplotsPlasma[[peptide]][[pair]] <- grid.arrange(plotsRatioPlasma[[peptide]][[pair]], plotStatsPlasma[[peptide]][[pair]], ncol = 2)



statplotsRandom <- list()
for(peptide in peptidenamesRandom) {
  statplotsRandom[[peptide]] <- grid.arrange(plotsRatioRandom[[peptide]], plotStatsRandom[[peptide]], ncol = 2)
}

pdf(paste("ChromAndDistribAndStats_run2-3_Random_4x4background"))
for(peptideA in PeptidesA) {
  for (peptideB in PeptidesB) {
    print(grid.arrange(plotChromatogram(StrepChroms, "run2", peptideA, StrepAnnot, printTitle =TRUE ),
                     plotChromatogram(StrepChroms, "run3", peptideB,StrepAnnot, printTitle =TRUE ),
                     statplotsRandom[[paste0("", peptideA," - ",peptideB)]],
                     ncol = 1))
  }
}
dev.off()


# Save plots to pdf
pdf(paste("ChromAndDistribAndAligned_run1-2_RandomPeptides_4x4background"))
for(peptideA in PeptidesA) {
  for (peptideB in PeptidesB) {
    print(grid.arrange(plotChromatogram(StrepChroms, "run1", peptideA, StrepAnnot, printTitle =TRUE ), 
                       plotSingleAlignedChrom(StrepChroms, peptideA, ),
                       plotChromatogram(StrepChroms, "run2", peptideB,StrepAnnot, printTitle =TRUE ),
                       plotsRandom[[paste0("", peptideA," - ",peptideB)]],
                       ncol = 2))
  }
}
dev.off()



# Fit a global alignment function between runs
run_pair <- c("run1", "run2")
loess.fit <- getLOESSfit(run_pair, peptides, oswOutStrep, 0.15)
StrepAnnot <- as.data.frame(StrepAnnot) # output of openSWATH (extracted Features - RT)
predict.run2 <- predict(loess.fit, data.frame(RUN1 = StrepAnnot[, run_pair[1]]))
Err <- predict.run2 - StrepAnnot[,run_pair[2]]


RSE4runsMatrix <- data.frame()
names <- c()
for (i in 1:length(globalStrep["RSE",]) ) {
  names <- unique(c(names, strsplit(rownames(as.data.frame(globalStrep["RSE",])), '_')[[i]]))
}

RSE4runsMatrix <- matrix(nrow = length(names), ncol = length(names))
colnames(RSE4runsMatrix) <- names
rownames(RSE4runsMatrix) <- names

for (i in 1:length(globalStrep["RSE",])) {
  RSE4runsMatrix[strsplit(strsplit(rownames(as.data.frame(globalStrep["RSE",])), '_')[[i]], " ")[[1]], 
                 strsplit(strsplit(rownames(as.data.frame(globalStrep["RSE",])), '_')[[i]], " ")[[2]] ] <- globalStrep["RSE",i]
  
  RSE4runsMatrix[strsplit(strsplit(rownames(as.data.frame(globalStrep["RSE",])), '_')[[i]], " ")[[2]], 
                 strsplit(strsplit(rownames(as.data.frame(globalStrep["RSE",])), '_')[[i]], " ")[[1]] ] <- globalStrep["RSE",i]
}

DistMat = as.dist(RSE4runsMatrix, diag = TRUE)
Clustered <- hclust(DistMat, method = "single")
library(ape)
class(Clustered) # must be hclust class
my_tree <- as.phylo(Clustered) 
plot(my_tree)


###################################################
########## PLASMA DATASET STARTS HERE #############
###################################################

filenames <- list.files("BloodPlasmaChroms/", pattern="*.rds")
filenames <- paste("BloodPlasmaChroms", filenames, sep="/")
count <- -1
StrepChromsPlasma <- list()
for (file in filenames) {
  count <- count + 1 
  Name <- paste("Plasma",count, sep = "") 
  StrepChromsPlasma[[Name]] <- readRDS(file)
}

filenames <- list.files(pattern="*.rds")
count <- -1
StrepChromsPlasma2 <- list()
for (file in filenames) {
  count <- count + 1 
  Name <- paste("Plasma",count, sep = "") 
  StrepChromsPlasma2[[Name]] <- readRDS(file)
}

filenames <- list.files("SwathOutputPlasma/", pattern="*.csv")
filenames <- paste("SwathOutputPlasma", filenames, sep="/")
count <- -1
oswOutStrepPlasma <- list()
for (file in filenames) {
  count <- count + 1 
  Name <- paste("Plasma",count, sep = "") 
  oswOutStrepPlasma[[Name]] <- read.csv(file, sep = '\t')
  oswOutStrepPlasma[[Name]] <- subset(oswOutStrepPlasma[[Name]], select =  c("transition_group_id", "RT", "m_score"))
}


runsPlasma <- names(StrepChromsPlasma)
peptidesPlasma <- c()
for (i in 1:length(StrepChromsPlasma$Plasma1)) {
  peptidesPlasma <- c(peptidesPlasma, names(StrepChromsPlasma$Plasma1[i])) 
}

AlignedChromsPlots <- list()
for (run in runsPlasma) {
  for (peptide in peptidesPlasma) {
    AlignedChromsPlots[[run]][[peptide]] <- NA
  }
}

library(readxl)
StrepAnnotPlasma <- read_excel("StrepAnnot_forPlasmaSet.xlsx")
names(StrepAnnotPlasma) <- gsub(' \\(sec\\)','', names(StrepAnnotPlasma))
names(StrepAnnotPlasma) <- gsub('run','Plasma', names(StrepAnnotPlasma))
rownames <- StrepAnnotPlasma$peptide_group_label
StrepAnnotPlasma <- as.data.frame(StrepAnnotPlasma[,2:25])
rownames(StrepAnnotPlasma) <- rownames

# Fit a local alignment between chromatogram groups
pair_names <- vector(); runs <- names(StrepChromsPlasma)
for (i in 1:(length(runsPlasma)-1)){
  for (j in (i+1): length(runsPlasma)){
    pair_names <- c(paste(runsPlasma[i], runsPlasma[j], sep = "_"), pair_names)
  }}

Score24runsPlasmaMatrices <- list()
for (peptide in peptidesPlasma) {
  Score24runsPlasmaMatrices[[peptide]] <- matrix(nrow = length(runsPlasma), ncol = length(runsPlasma))
  colnames(Score24runsPlasmaMatrices[[peptide]]) <- runsPlasma
  rownames(Score24runsPlasmaMatrices[[peptide]]) <- runsPlasma
}

# Alignment, calculation of total number of gaps, background sums, maxximum scores (local), ratio, 
AlignedChromsPlots <- list()
BackgroundSumPlasma <- list()
MaxScoresAllPlasma <- list()
RatioPlasma <- list()
plotsPlasma <- list()
plotsRatioPlasma <- list()
plotStatsPlasma <- list()
ScoreDataFramesPlasma <- list()
StatsDataFramePlasma <- list()
AUC <- list()
for(peptide in peptidesPlasma[1]){
  pdf(paste("AlignedChroms",gsub("/","-",peptide),"runsPlasma.pdf", sep = "_"))
  for(pair in pair_names){
    print(pair)
    run_pair <- strsplit(pair, split = "_")[[1]]
    gapQuantile <- 0.5; goFactor <- 1/8; geFactor <- 40
    simMeasure <- "dotProductMasked"
    Err <- matrix(NA, nrow = length(peptidesPlasma[1]), ncol = 1)
    rownames(Err) <- peptidesPlasma[1]
    ListOfSimMatricesPlasma <- list()
    ListOfMaxScoresPlasma <- list()
    ListOfBackgroundSumPlasma <- list()
    ListOfAlignedIndicesPlasma <- list()
    NumberOfGapsPlasma <- c()
    NumberOfGapsListPlasma <- list()
    AverageMaxScorePlasma <- list()
    MedianOfMPlasma <- list()

    s <- getSimilarityMatrix(StrepChromsPlasma, peptide, run_pair[1], run_pair[2], type = simMeasure)
    gapPenalty <- getGapPenalty(s, gapQuantile, type = simMeasure)
    if (round(gapPenalty,2) == 0) {
      while (gapPenalty < 0.01) {
        gapQuantile <- gapQuantile + 0.1
        gapPenalty <- getGapPenalty(s, gapQuantile, type = simMeasure)
      }
    }
    print(gapPenalty)
    
    Alignobj <- getAffineAlignObj(s, go = gapPenalty*goFactor, ge = gapPenalty*geFactor)
    ListOfSimMatricesPlasma[[peptide]][[pair]] <- Alignobj@M
    AlignedIndices[[peptide]][[pair]] <- getAlignment(Alignobj)
    
    NumberOfGapsPlasma <-c(NumberOfGapsPlasma, sum(is.na(AlignedIndices[[peptide]][[pair]][[1]][TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,1])[1]:TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,1])[2],1]), 
                                       is.na(AlignedIndices[[peptide]][[pair]][[1]][TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,1])[1]:TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,1])[2],2]),
                                       is.na(AlignedIndices[[peptide]][[pair]][[1]][TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,2])[1]:TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,2])[2],1]), 
                                       is.na(AlignedIndices[[peptide]][[pair]][[1]][TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,2])[1]:TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,2])[2],2])))
    
    NumberOfGapsListPlasma[[peptide]][[pair]] <- sum(is.na(AlignedIndices[[peptide]][[pair]][[1]][TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,1])[1]:TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,1])[2],1]), 
                                                     is.na(AlignedIndices[[peptide]][[pair]][[1]][TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,1])[1]:TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,1])[2],2]),
                                                     is.na(AlignedIndices[[peptide]][[pair]][[1]][TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,2])[1]:TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,2])[2],1]), 
                                                     is.na(AlignedIndices[[peptide]][[pair]][[1]][TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,2])[1]:TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,2])[2],2]))

    ListOfAlignedIndicesPlasma[[peptide]][[pair]] <- AlignedIndices
    AverageMaxScorePlasma[[peptide]][[pair]] <- AlignedIndices[[peptide]][[pair]][[1]][length(AlignedIndices[[peptide]][[pair]][[1]][,3]), 3] /
                                                           length(AlignedIndices[[peptide]][[pair]][[1]][ ,3])
    MedianOfMPlasma[[peptide]][[pair]] <- median(Alignobj@M)
    ListOfMaxScoresPlasma[[peptide]][[pair]] <- AlignedIndices[[peptide]][[pair]][[1]][,"score"]
    ListOfBackgroundSumPlasma[[peptide]][[pair]] <- SumNearMaxScores(Alignobj)
    tA <- StrepChromsPlasma[[run_pair[1]]][[peptide]][[1]][["time"]]
    tB <- StrepChromsPlasma[[run_pair[2]]][[peptide]][[1]][["time"]]
    tA.aligned <- mapIdxToTime(tA, AlignedIndices[[peptide]][[pair]][[1]][,"indexA_aligned"])
    tB.aligned <- mapIdxToTime(tB, AlignedIndices[[peptide]][[pair]][[1]][,"indexB_aligned"])
    predictTime <- tB.aligned[which.min(abs(tA.aligned - StrepAnnotPlasma[peptide, run_pair[1]]))]
    deltaT <- predictTime - StrepAnnotPlasma[peptide, run_pair[2]]
    Err[peptide, 1] <- deltaT
    
    # plotAlignedChroms(StrepChromsPlasma, pair, peptide, StrepAnnotPlasma, #AlignErrorinSec,
    #                                                                      AlignedIndices, tA, tB, FourOrTwo = TRUE)
  
  
    # from lists of cumulative scores create lists of separate max and background scores
    # and list of dataframes for each peptide (to create a plot)

    BackgroundSumPlasma[[peptide]][[pair]] <- c(rep(NA, length(ListOfBackgroundSumPlasma[[peptide]][[pair]])))
    for (i in length(BackgroundSumPlasma[[peptide]][[pair]]):2) {
      BackgroundSumPlasma[[peptide]][[pair]][i] <- ListOfBackgroundSumPlasma[[peptide]][[pair]][i-1] -
        ListOfBackgroundSumPlasma[[peptide]][[pair]][i]
      if (is.na(BackgroundSumPlasma[[peptide]][[pair]][i]) | (BackgroundSumPlasma[[peptide]][[pair]][i] %in% Inf)) {
        BackgroundSumPlasma[[peptide]][[pair]] <- BackgroundSumPlasma[[peptide]][[pair]][-i]
      }
    }
    MaxScoresAllPlasma[[peptide]][[pair]] <- c(rep(NA,length(ListOfMaxScoresPlasma[[peptide]][[pair]])))
    for (i in (length(MaxScoresAllPlasma[[peptide]][[pair]]) -1):1) {
      MaxScoresAllPlasma[[peptide]][[pair]][i] <- ListOfMaxScoresPlasma[[peptide]][[pair]][i+1] - ListOfMaxScoresPlasma[[peptide]][[pair]][i]
      if (is.na(MaxScoresAllPlasma[[peptide]][[pair]][i]) | (MaxScoresAllPlasma[[paste0("", peptide)]][[pair]][i] %in% Inf)) {
        MaxScoresAllPlasma[[peptide]][[pair]] <- MaxScoresAllPlasma[[peptide]][[pair]][-i]
      }
    }
    RatioPlasma[[peptide]][[pair]] <- MaxScoresAllPlasma[[peptide]][[pair]][length(BackgroundSumPlasma[[peptide]][[pair]]):1] / BackgroundSumPlasma[[peptide]][[pair]]
    
    #### AUC calculation 
    m <- ggplot(as.data.frame(RatioPlasma[[peptide]][[pair]]), aes(x = as.numeric(RatioPlasma[[peptide]][[pair]]))) 
    m <- m + geom_density() 
    p <- ggplot_build(m)
    y_max <- max(p$data[[1]]$y)
    AUC[[peptide]][[pair]] <- auc(p$data[[1]]$x,p$data[[1]]$y)
    print(AUC[[peptide]][[pair]])
    
    ScoreDataFramesPlasma[[peptide]][[pair]] <- data.frame(rbind(cbind(BackgroundSumPlasma[[peptide]], c(rep("Background Sum", 
                                                                                         length(BackgroundSumPlasma[[peptide]][[pair]])))),
                                                   cbind(MaxScoresAllPlasma[[peptide]][[pair]], c(rep("Maximum Score",
                                                                                        length(MaxScoresAllPlasma[[peptide]][[pair]]))))) )
    names(ScoreDataFramesPlasma[[peptide]][[pair]]) <- c("Score", "Type")
    
    StatsDataFramePlasma[[peptide]][[pair]] <- data.frame(rbind(c(NumberOfGapsListPlasma[[peptide]][[pair]], "Number of gaps"),
                                                  c(AverageMaxScorePlasma[[peptide]][[pair]], "Average maximum score"),
                                                  c(MedianOfMPlasma[[peptide]][[pair]], "Median of the score matrics"),
                                                  c(median(RatioPlasma[[peptide]][[pair]], na.rm = T), "Median of a MaxScore/Background ratio")))
    
    # Fill in score matrix with values for subsequent h-clustering
    Score24runsPlasmaMatrices[[peptide]][run_pair[1],run_pair[2]] <- (1 + NumberOfGapsListPlasma[[peptide]][[pair]])/AUC[[peptide]][[pair]]
    Score24runsPlasmaMatrices[[peptide]][run_pair[2],run_pair[1]] <- (1 + NumberOfGapsListPlasma[[peptide]][[pair]])/AUC[[peptide]][[pair]]
    if (is.na(Score24runsPlasmaMatrices[[peptide]][run_pair[1],run_pair[2]])) {
      Score24runsPlasmaMatrices[[peptide]][run_pair[1],run_pair[2]] <- median(Score24runsPlasmaMatrices[[peptide]][run_pair[1],],na.rm = T)
      Score24runsPlasmaMatrices[[peptide]][run_pair[2],run_pair[1]] <- median(Score24runsPlasmaMatrices[[peptide]][run_pair[1],],na.rm = T)
    }
    
    # store some plots
    plotsRatioPlasma[[peptide]][[pair]] <- ggplot(as.data.frame(RatioPlasma[[peptide]][[pair]]), aes(x = as.numeric(RatioPlasma[[peptide]][[pair]]))) +   
      geom_density(alpha = 0.3, color = "darkblue", fill = "darkblue") +
      labs(x = "Maximum score at each step /
           Background sum (4x4)") +
      theme(axis.title.y = element_text(size = 7),
            axis.title.x = element_text(size = 6)) +
      annotate("text", label = paste("AUC = ", round(AUC[[peptide]][[pair]],2)), x = max(RatioPlasma[[peptide]][[pair]], na.rm = T) - 
                 max(RatioPlasma[[peptide]][[pair]], na.rm = T)/10, y_max - y_max/10, hjust = 1) +
      annotate("text", label = paste("mean = ", round(mean(RatioPlasma[[peptide]][[pair]],na.rm = T),2)),
               x = max(RatioPlasma[[peptide]][[pair]], na.rm = T) - max(RatioPlasma[[peptide]][[pair]], na.rm = T)/10, y = y_max - 2*y_max/10,hjust = 1) +
      annotate("text", label = paste("median = ", round(median(RatioPlasma[[peptide]][[pair]], na.rm = T),2)),
               x = max(RatioPlasma[[peptide]][[pair]], na.rm = T) - max(RatioPlasma[[peptide]][[pair]], na.rm = T)/10, y = y_max - 3*y_max/10,hjust = 1) +
      annotate("text", label = paste("sd = ", round(sd(RatioPlasma[[peptide]][[pair]], na.rm = T),2)),
               x = max(RatioPlasma[[peptide]][[pair]], na.rm = T) - max(RatioPlasma[[peptide]][[pair]], na.rm = T)/10, y = y_max - 4*y_max/10,hjust = 1)
    
    plotStatsPlasma[[peptide]][[pair]] <- ggplot(subset(StatsDataFramePlasma[[peptide]][[pair]], !(StatsDataFramePlasma[[peptide]][[pair]]$X2 %in%
                                                                                                     c("Median of a MaxScore/Background ratio"))),
                                           aes(x = X2, y = as.numeric(as.character(X1)))) +
      geom_bar(stat = "identity", fill = "darkblue") +
      labs( x = "",
            y = "") +
      theme(axis.text.y = element_text(size = 7),
            axis.text.x = element_text(size = 6)) 
    
   # statplotsPlasma[[peptide]][[pair]] <- grid.arrange(plotsRatioPlasma[[peptide]][[pair]], plotStatsPlasma[[peptide]][[pair]], ncol = 2)
    
    # save plots to pdf
    print(grid.arrange(plotAlignedChroms(StrepChromsPlasma, pair, peptide, StrepAnnotPlasma, #AlignErrorinSec,
                                         AlignedIndices[[peptide]][[pair]], tA, tB, FourOrTwo = TRUE),
                grid.arrange(plotsRatioPlasma[[peptide]][[pair]], plotStatsPlasma[[peptide]][[pair]], ncol = 2),
                       ncol =1))
    
    }
  dev.off()
}

# Clustering using the score matrices produced as a result of the chosen distance (1+NGaps)/AUC
Clusters1Peptide <- list()
DistMatrices <- list()
pdf(paste("Trees_24runsPlasma_10peptides_upgma_GapsAUC.pdf"))
for(peptide in peptidesPlasma) {
  print(peptide)
  DistMatrices[[peptide]] = as.dist(Score24runsPlasmaMatrices[[peptide]], diag = TRUE)
  Clustered <- hclust(DistMatrices[[peptide]], method = "average")
  my_tree <- as.phylo(Clustered) 
  plot(my_tree,main = peptide,label.offset = 0.01)
  
  # save tree in the newick format
  #write.tree(phy=my_tree, file= paste(gsub("/","-",peptide),"_plasma_tree.newick", sep = ""))
}
dev.off()

# Weights obtained from the tree
Weights <- list()
weights <- read.table("weights_116795_TFISPIK-2.output", sep = "\t", header = F)
rownames(weights) <- weights[,1]
weights[,1] <- NULL
Weights[[peptide]] <- weights

peptide <- "116795_TFISPIK/2"

### UPGMA algorithm to save the sequence of merging runs within the tree
### Gives as an output list like: (A, AB, CD, ABCD, ABCDE, KF, KFX, ABCDEKFX)

ScoreMatrix <- as.matrix(Score24runsPlasmaMatrices[[1]])
for (i in 1:length(ScoreMatrix[1,])) {ScoreMatrix[i,i] <- NA}

# first, find the cell of a matrix with lowest value
# Locate the smallest cell in the table
lowest_cell <- function(ScoreMatrix) {
  # Set default to infinity
  min_cell <- Inf
  x <- -1
  y <- -1

  # Go through every cell, looking for the lowest
  for (i in 1:(length(ScoreMatrix[1,]))) {
    for (j in 1:(length(ScoreMatrix[i,]))) {
      if (!is.na(ScoreMatrix[i,j])) {
        if (ScoreMatrix[i,j] < min_cell) {
         min_cell <- ScoreMatrix[i,j]
         x <- i
         y <- j
      }}}}
  return(c(x,y))
}

# join_labels: combine two merging labels (names of runs) and add to list of labels to keep
# the sequence of merge events. Like: c(A, AB, CD, ABCD, ABCDE, KF, KFX, ABCDEKFX)
join_labels <- function(ScoreMatrix, labels, a, b) {
  # Swap if the indices are not ordered
  if (b < a) {
    temp <- a
    a <- b
    b <- temp
  }
  # Join the labels (merging event AB + CDE --> ABCDE added to list)
  labels <- c(labels, paste(dimnames(ScoreMatrix)[[1]][a], 
                            dimnames(ScoreMatrix)[[2]][b], sep = "_"))
  return(labels)
}

# join_table: Join the entries of a table on the cell (a, b) by averaging their data entries
join_table <- function(ScoreMatrix, a, b) {
  # Swap if the indices are not ordered
  if (b < a) {
    temp <- a
    a <- b
    b <- temp
  }
  # For the lower index, reconstruct the entire row (A, i), where i < A
  row <- c()
  for (i in 1:a) {
    row <- c(row, ((ScoreMatrix[a,i] + ScoreMatrix[b,i])/2))
  }
  ScoreMatrix[a, c(1:a)] <- row
 
  # Then, reconstruct the entire column (i, A), where i > A
  # Note: Since the matrix is lower triangular, row b only contains values for indices < b
  for (i in (a+1):b) {
    ScoreMatrix[a,i] <- (ScoreMatrix[a,i] + ScoreMatrix[b,i])/2
  }
  #   We get the rest of the values from row i
  if ((b+1) < length(ScoreMatrix[1,])) {
    for (i in (b+1):length(ScoreMatrix[1,])) {
      ScoreMatrix[i,a] = (ScoreMatrix[i,a] + ScoreMatrix[i,b])/2
    # Remove the (now redundant) second index column entry
      ScoreMatrix[i,b] <- NA
  }
  }
  
  # Remove the (now redundant) second index row and column and rename the first index row and column
  LastMergeName <- paste(names(ScoreMatrix[b,])[a],names(ScoreMatrix[a,])[b], sep = "_")
  rownames(ScoreMatrix)[a] <- LastMergeName
  colnames(ScoreMatrix)[a] <- LastMergeNameт
  ScoreMatrix <- ScoreMatrix[-b,-b] 
  return(ScoreMatrix)
}

UPGMA <- function(ScoreMatrix, labels) {
  for (i in 1:(length(ScoreMatrix[1,]) - 1)) {
  # Locate lowest cell in the ScoreMatrix
    x <- lowest_cell(ScoreMatrix)[1]
    y <- lowest_cell(ScoreMatrix)[2]
    # Update the labels 
    labels <- join_labels(ScoreMatrix,labels, x, y)
     # Join the ScoreMatrix on the cell co-ordinates
    ScoreMatrix <- join_table(ScoreMatrix, x, y)
  }
  
    # Return the final label
  return(labels)
}

labels <- c()
MergingSequenceLabels <- UPGMA(ScoreMatrix, labels)


Intensity <- list()
Intensity_mean <- list()
Intensity_median <- list()
Intensity_weightedmean <- list()
mergedStrepChromsPlasma <- list()

########################################################
########### MAIN FUNCTION FOR MERGING ##################
########################################################
AlignAndMergeRuns <- function(data, mergedData, run1, run2, peptide, newLengthofData, allNAs = TRUE, randomNAs = FALSE,
                              typeOfRuns = c("single", "singleAndMerged","merged"), resample = TRUE) {
  # if single run + merged, always assign single run to run1 !!!
  
  pair <- paste(run1, run2, sep = "_")
  run_pair <- c(run1, run2)
  gapQuantile <- 0.5; goFactor <- 1/8; geFactor <- 40
  simMeasure <- "dotProductMasked"
  
  typeOfRuns <- match.arg(typeOfRuns)
  switch(typeOfRuns, single = {
    data1 <- data
    data2 <- data
    s <- getSimilarityMatrix(data, peptide, run_pair[1], run_pair[2], type = simMeasure)
  },singleAndMerged = {
    data1 <- data
    data2 <- mergedData
    s <- getSimilarityMatrix_forMerged(data, mergedData, peptide, run_pair[1], run_pair[2], type = simMeasure)
  }, merged = {
    data1 <- mergedData
    data2 <- mergedData
    s <- getSimilarityMatrix(mergedData, peptide, run_pair[1], run_pair[2], type = simMeasure)
  })
  
  print(run_pair[1])
  print(run_pair[2])

  gapPenalty <- getGapPenalty(s, gapQuantile, type = simMeasure)
  if (round(gapPenalty,2) == 0) {
    while (gapPenalty < 0.01) {
      gapQuantile <- gapQuantile + 0.1
      gapPenalty <- getGapPenalty(s, gapQuantile, type = simMeasure)
    }
  }
  Alignobj <- getAffineAlignObj(s, go = gapPenalty*goFactor, ge = gapPenalty*geFactor)
  AlignedIndices[[peptide]][[pair]] <- getAlignment(Alignobj)
  timepointsNumber <- length(AlignedIndices[[peptide]][[pair]][[1]][,1])

  Intensity[[peptide]][[pair]] <- list()
  
  for (i in 1:length(data1[[run_pair[1]]][[peptide]])) {
    IndicesToBeSkipped <- c()
    Intensity[[peptide]][[pair]][[i]] <- c(rep(NA, 220))
  
    for (index in 1:timepointsNumber) {
      ### calculate new intensities according to the relabeling model considering all NAs
      # the index corresponding to the timepoint[j] in StrepChromsPlasma is AlignedIndices[[peptide]][[pair]][[1]][timepoint[j],1]
      if (allNAs) {
        if (is.na(AlignedIndices[[peptide]][[pair]][[1]][index,1][[1]])) {
          if (index <= TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,1])[1])  {
            Intensity[[peptide]][[pair]][[i]][index] <- mean(c(data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,1][[1]], 2],
                                                               data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,2][[1]], 2]))
          }
          else if (index > (timepointsNumber - TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,1])[2])) {
            Intensity[[peptide]][[pair]][[i]][index] <- mean(c(data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][timepointsNumber - 
                                                                                 TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][, 1])[2], 1], 2],
                                                               data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index, 2][[1]], 2]))
            
          }
          else {
            Intensity[[peptide]][[pair]][[i]][index] <- mean(c(mean(c(data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index - 1,1][[1]], 2],
                                                                      data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index + 1,1][[1]], 2])),
                                                               data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,2][[1]], 2] ))
          }
        }

        else if (!is.na(AlignedIndices[[peptide]][[pair]][[1]][index, 1][[1]]))  {
          if (!is.na(AlignedIndices[[peptide]][[pair]][[1]][index, 2][[1]])) {
            Intensity[[peptide]][[pair]][[i]][index] <- mean(c(data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,1][[1]], 2],
                                                               data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,2][[1]], 2]))
          }
          if (is.na(AlignedIndices[[peptide]][[pair]][[1]][index, 2])) {
            if (index <= TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,2])[1])  {
              Intensity[[peptide]][[pair]][[i]][index] <- mean(c(data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,1][[1]], 2],
                                                                 data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,2][[1]], 2]))
            }
            else if (index > (timepointsNumber - TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,2])[2])) {
              Intensity[[peptide]][[pair]][[i]][index] <- mean(c(data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,1][[1]], 2],
                                                                 data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index, 2][[1]] -
                                                                                                                    TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][, 2])[2], 2]))
              
            }
            
            else {
              Intensity[[peptide]][[pair]][[i]][index] <- mean(c(mean(c(data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index - 1,2][[1]], 2],
                                                                        data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index + 1,2][[1]], 2])),
                                                                 data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,1][[1]], 2] ))
            }
          }
        }
      }
      ### calculate new intensities according to the relabeling model skipping all intermediate NAs (non-terminal)
      if (!allNAs) {
        if (is.na(AlignedIndices[[peptide]][[pair]][[1]][index,1][[1]])) {
          if (index <= TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,1])[1])  {
            Intensity[[peptide]][[pair]][[i]][index] <- mean(c(data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,1][[1]], 2],
                                                               data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,2][[1]], 2]))
          }
          else if (index > (timepointsNumber - TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,1])[2])) {
            Intensity[[peptide]][[pair]][[i]][index] <- mean(c(data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][timepointsNumber - 
                                                                                                                                             TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][, 1])[2], 1], 2],
                                                               data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index, 2][[1]], 2]))
            
          }
          else {
            IndicesToBeSkipped <- c(IndicesToBeSkipped, index)
          }
        }
        
        else if (!is.na(AlignedIndices[[peptide]][[pair]][[1]][index, 1][[1]]))  {
          if (!is.na(AlignedIndices[[peptide]][[pair]][[1]][index, 2][[1]])) {
            Intensity[[peptide]][[pair]][[i]][index] <- mean(c(data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,1][[1]], 2],
                                                               data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,2][[1]], 2]))
          }
          if (is.na(AlignedIndices[[peptide]][[pair]][[1]][index, 2])) {
            if (index <= TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,2])[1])  {
              Intensity[[peptide]][[pair]][[i]][index] <- mean(c(data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,1][[1]], 2],
                                                                 data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,2][[1]], 2]))
            }
            else if (index > (timepointsNumber - TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][,2])[2])) {
              Intensity[[peptide]][[pair]][[i]][index] <- mean(c(data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index,1][[1]], 2],
                                                                 data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][index, 2][[1]] -
                                                                                                        TerminalNAnumber(AlignedIndices[[peptide]][[pair]][[1]][, 2])[2], 2]))
              
            }
            else {
              IndicesToBeSkipped <- c(IndicesToBeSkipped, index)
               }
          }
        }
      }
      ### skipping randomly NAs
      if (randomNAs == TRUE) {
        LeftNAsIndices <- sample(IndicesToBeSkipped, size = round(length(IndicesToBeSkipped)/2))
        for (idx in LeftNAsIndices) {
          if (is.na(AlignedIndices[[peptide]][[pair]][[1]][idx, 2])) {
            Intensity[[peptide]][[pair]][[i]][idx] <- mean(c(mean(c(data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][idx - 1,2][[1]], 2],
                                                                      data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][idx + 1,2][[1]], 2])),
                                                               data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][idx,1][[1]], 2] ))
    
          }
          if (is.na(AlignedIndices[[peptide]][[pair]][[1]][idx, 1]))   {
            Intensity[[peptide]][[pair]][[i]][idx] <- mean(c(mean(c(data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][idx - 1,1][[1]], 2],
                                                                      data1[[run_pair[1]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][idx + 1,1][[1]], 2])),
                                                               data2[[run_pair[2]]][[peptide]][[i]][AlignedIndices[[peptide]][[pair]][[1]][idx,2][[1]], 2] ))
          }
        }
      }
    }
    
    Intensity[[peptide]][[pair]][[i]] <- Intensity[[peptide]][[pair]][[i]][!is.na(Intensity[[peptide]][[pair]][[i]])]
   
    mergedData[[pair]][[peptide]][[i]] <- data.frame(cbind(c(1:length(Intensity[[peptide]][[pair]][[i]][1])), Intensity[[peptide]][[pair]][[i]]))
  
    # resample back to newLengthofData points:
    if (resample) {
      interpolation <- approx(mergedData[[pair]][[peptide]][[i]][,2],
                            method = "linear",
                            n = newLengthofData)
      mergedData[[pair]][[peptide]][[i]] <- data.frame(cbind(c(1:newLengthofData), interpolation$y))
      print(i)
      
      # optional step, just didn't manage to save mergedData in other way soooo...
      # mergedStrepChromsPlasma[[pair]][[peptide]][[i]] <- mergedData[[pair]][[peptide]][[i]]
        
    }
  }

  if (typeOfRuns == "single") {
    BeforeMergeSingle <- grid.arrange(plotChromatogram(data1, run_pair[1],peptide,StrepAnnotPlasma,TRUE),
                                    plotChromatogram(data2,run_pair[2],peptide,StrepAnnotPlasma, TRUE), ncol = 2)
    plot <- grid.arrange(BeforeMergeSingle, 
                 plotChromatogramMerged(mergedData,pair,peptide, TRUE))
  }
  if (typeOfRuns == "singleAndMerged") {
    BeforeMergeBoth <- grid.arrange(plotChromatogramMerged(data2, run_pair[2], peptide, TRUE) ,
                                  plotChromatogram(data1,run_pair[1],peptide,StrepAnnotPlasma, TRUE), ncol = 2)
    plot <- grid.arrange(BeforeMergeBoth, 
                 plotChromatogramMerged(mergedData,pair,peptide, TRUE))
  }
  if (typeOfRuns == "merged") {
    BeforeMergeMerged <- grid.arrange(plotChromatogramMerged(data1, run_pair[1], peptide, TRUE) ,
                                    plotChromatogramMerged(data2, run_pair[2], peptide, TRUE), ncol = 2)
    plot <- grid.arrange(BeforeMergeMerged, 
                 plotChromatogramMerged(mergedData,pair,peptide, TRUE))
  }
  
  if (allNAs) {
    saveRDS(mergedData, file = paste(paste("mergedStrepChroms_with", pair, gsub("/","-", peptide), sep = "_"), "rds", sep = "."))
  }
  else if ((!allNAs) & (randomNAs == FALSE)) {
    saveRDS(mergedData, file = paste(paste("mergedStrepChromsSkipNAs_with", pair, gsub("/","-", peptide), sep = "_"), "rds", sep = "."))
  }
  else if ((!allNAs) & (randomNAs == TRUE)) {
    saveRDS(mergedData, file = paste(paste("mergedStrepChromsRandomNAs_with", pair, gsub("/","-", peptide), sep = "_"), "rds", sep = "."))
  }
  return(plot)
}

### First, Plasma5-18 and 20 --> 20-5-18 + 7 --> 7-20-5-18 + 0-19 --> 0-19-7-20-5-18 + 13-15 --> 0-19-7-20-5-18-13-15
### Secondly, plasma-17-6, 10-2 --> 17-6-10-2, then 8 and 16-12 --> 8-16-12 and then join: 17-6-10-2-8-16-12
### Thirdly, 0-19-7-20-5-18-13-15 + 17-6-10-2-8-16-12 --->  0-19-7-20-5-18-13-15-17-6-10-2-8-16-12
### Fourth, add 9-4
### Then, plasma22-3 and 11-14 + 23--> 23-11-14 + 22-3 --> 23-11-14-22-3
### Pre-last : 9-4-0-19-7-20-5-18-13-15-17-6-10-2-8-16-12 + 23-11-14-22-3
### Last: + 1-21

### Perform Alignment and merging in order according to MergingSequenceLabels
for (i in 1:length(MergingSequenceLabels)) {
  if (length(strsplit(MergingSequenceLabels[i], "_")) == 2) {
    AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaRandomNAs, strsplit(MergingSequenceLabels[i])[1],
                      strsplit(MergingSequenceLabels[i])[2], allNAs = FALSE, randomNAs = TRUE, 
                      peptide, newLengthofData = 195, typeOfRuns = "single", resample = TRUE)
    AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, strsplit(MergingSequenceLabels[i])[1],
                      strsplit(MergingSequenceLabels[i])[2], allNAs = FALSE, randomNAs = FALSE, 
                      peptide, newLengthofData = 195, typeOfRuns = "single", resample = TRUE)
    AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaAllNAs, strsplit(MergingSequenceLabels[i])[1],
                      strsplit(MergingSequenceLabels[i])[2], allNAs = TRUE, randomNAs = FALSE, 
                      peptide, newLengthofData = 195, typeOfRuns = "single", resample = TRUE)
    
    
  }
  else if  {
    
  }
}

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaRandomNAs, "Plasma5", "Plasma18", allNAs = FALSE, randomNAs = TRUE, 
                  peptide, newLengthofData = 195, typeOfRuns = "single", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma5_Plasma18_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma5_Plasma18_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma5_Plasma18_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaRandomNAs, "Plasma20","Plasma5_Plasma18",allNAs = FALSE, randomNAs = T,
                  peptide, newLengthofData = 195, typeOfRuns = "singleAndMerged", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma20_Plasma5_Plasma18_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma20_Plasma5_Plasma18_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma20_Plasma5_Plasma18_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma7","Plasma20_Plasma5_Plasma18",allNAs = FALSE, randomNAs = T,
                  peptide, newLengthofData = 195, typeOfRuns = "singleAndMerged", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma7_Plasma20_Plasma5_Plasma18_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma7_Plasma20_Plasma5_Plasma18_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma7_Plasma20_Plasma5_Plasma18_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaRandomNAs, "Plasma0", "Plasma19",allNAs = FALSE, randomNAs = T,
                  peptide, newLengthofData = 195, typeOfRuns = "single", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma0_Plasma19_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma0_Plasma19_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma0_Plasma19_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma0_Plasma19","Plasma7_Plasma20_Plasma5_Plasma18",allNAs = FALSE, randomNAs = F,
                  peptide, newLengthofData = 195, typeOfRuns = "merged", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma15", "Plasma13",allNAs = FALSE, randomNAs = F,
                  peptide, newLengthofData = 195, typeOfRuns = "single", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma15_Plasma13_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma15_Plasma13_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma15_Plasma13_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs,"Plasma15_Plasma13", "Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18",
                  allNAs = FALSE, randomNAs = F,
                  peptide, newLengthofData = 195, typeOfRuns = "merged", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma15_Plasma13_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma15_Plasma13_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma15_Plasma13_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma17", "Plasma6",allNAs = FALSE, randomNAs = F,
                  peptide, newLengthofData = 195, typeOfRuns = "single", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma17_Plasma6_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma17_Plasma6_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma17_Plasma6_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasma, "Plasma10", "Plasma2",allNAs = TRUE, randomNAs = F,
                  peptide, newLengthofData = 195, typeOfRuns = "single", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma10_Plasma2_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma10_Plasma2_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma10_Plasma2_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs,"Plasma10_Plasma2", "Plasma17_Plasma6",allNAs = FALSE, randomNAs = F,
                  peptide, newLengthofData = 195, typeOfRuns = "merged", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma10_Plasma2_Plasma17_Plasma6_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma10_Plasma2_Plasma17_Plasma6_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma10_Plasma2_Plasma17_Plasma6_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma16", "Plasma12",allNAs = FALSE, randomNAs = F,
                  peptide, newLengthofData = 195, typeOfRuns = "single", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma16_Plasma12_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma16_Plasma12_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma16_Plasma12_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma8","Plasma16_Plasma12",allNAs = FALSE, randomNAs = F,
                  peptide, newLengthofData = 195, typeOfRuns = "singleAndMerged", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma8_Plasma16_Plasma12_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma8_Plasma16_Plasma12_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma8_Plasma16_Plasma12_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs,"Plasma8_Plasma16_Plasma12", "Plasma10_Plasma2_Plasma17_Plasma6", 
                  allNAs = FALSE, randomNAs = F,
                  peptide, newLengthofData = 195, typeOfRuns = "merged", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18",
                  "Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6",peptide, newLengthofData = 195,
                  allNAs = FALSE, randomNAs = F, typeOfRuns = "merged", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma4", "Plasma9",allNAs = FALSE, randomNAs = F,
                  peptide, newLengthofData = 195, typeOfRuns = "single", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma4_Plasma9_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma4_Plasma9_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma4_Plasma9_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6",
                  "Plasma4_Plasma9",peptide, newLengthofData = 195, allNAs = FALSE, randomNAs = F,
                  typeOfRuns = "merged", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_Plasma4_Plasma9_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_Plasma4_Plasma9_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_Plasma4_Plasma9_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma23","Plasma11_Plasma14",allNAs = FALSE, randomNAs = F,
                  peptide, newLengthofData = 195, typeOfRuns = "singleAndMerged", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma23_Plasma11_Plasma14_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma23_Plasma11_Plasma14_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma23_Plasma11_Plasma14_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma3", "Plasma22",allNAs = FALSE, randomNAs = F,
                  peptide, newLengthofData = 195, typeOfRuns = "single", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma3_Plasma22_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma3_Plasma22_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma3_Plasma22_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma3_Plasma22","Plasma23_Plasma11_Plasma14",
                  allNAs = FALSE, randomNAs = F,
                  peptide, newLengthofData = 195, typeOfRuns = "merged", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma3_Plasma22_Plasma23_Plasma11_Plasma14_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma3_Plasma22_Plasma23_Plasma11_Plasma14_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma3_Plasma22_Plasma23_Plasma11_Plasma14_116795_TFISPIK-2.rds")
 
AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma3_Plasma22_Plasma23_Plasma11_Plasma14",
                  "Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_Plasma4_Plasma9",
                  allNAs = FALSE, randomNAs = F, peptide, newLengthofData = 195, typeOfRuns = "merged", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma3_Plasma22_Plasma23_Plasma11_Plasma14_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_Plasma4_Plasma9_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma3_Plasma22_Plasma23_Plasma11_Plasma14_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_Plasma4_Plasma9_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma3_Plasma22_Plasma23_Plasma11_Plasma14_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_Plasma4_Plasma9_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma1", "Plasma21",allNAs = FALSE, randomNAs = F,
                  peptide, newLengthofData = 195, typeOfRuns = "single", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma1_Plasma21_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma1_Plasma21_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma1_Plasma21_116795_TFISPIK-2.rds")

AlignAndMergeRuns(StrepChromsPlasma, mergedStrepChromsPlasmaSkipNAs, "Plasma1_Plasma21",
                  "Plasma3_Plasma22_Plasma23_Plasma11_Plasma14_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_Plasma4_Plasma9",
                  allNAs = FALSE, randomNAs = F, peptide, newLengthofData = 195, typeOfRuns = "merged", resample = TRUE)

mergedStrepChromsPlasma <- readRDS("mergedStrepChroms_with_Plasma1_Plasma21_Plasma3_Plasma22_Plasma23_Plasma11_Plasma14_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_Plasma4_Plasma9_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaSkipNAs <- readRDS("mergedStrepChromsSkipNAs_with_Plasma1_Plasma21_Plasma3_Plasma22_Plasma23_Plasma11_Plasma14_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_Plasma4_Plasma9_116795_TFISPIK-2.rds")
mergedStrepChromsPlasmaRandomNAs <- readRDS("mergedStrepChromsRandomNAs_with_Plasma1_Plasma21_Plasma3_Plasma22_Plasma23_Plasma11_Plasma14_Plasma0_Plasma19_Plasma7_Plasma20_Plasma5_Plasma18_Plasma8_Plasma16_Plasma12_Plasma10_Plasma2_Plasma17_Plasma6_Plasma4_Plasma9_116795_TFISPIK-2.rds")



######## Get RSE and create a distance table
samples4gradient <- 100; RSEdistFactor <- 3.5; hardConstrain <- FALSE
pair_names <- vector(); runs <- names(StrepChromsPlasma)
for (i in 1:(length(runs)-1)){
  for (j in (i+1): length(runs)){
    pair_names <- c(paste(runs[i], runs[j], sep = "_"), pair_names)
  }}
globalStrep <- matrix(NA, nrow = 1, ncol = length(pair_names))
colnames(globalStrep) <- pair_names
rownames(globalStrep) <- c("RSE")

for(pair in pair_names){
  run_pair <- strsplit(pair, split = "_")[[1]]
  Loess.fit <- getLOESSfit(run_pair, peptidesPlasma[1:50], oswOutStrepPlasma, 0.1)
  globalStrep["RSE", pair] <- Loess.fit$s
}

RSErunsPlasmaMatrix <- matrix(nrow = length(runsPlasma), ncol = length(runsPlasma))
colnames(RSErunsPlasmaMatrix) <- runsPlasma
rownames(RSErunsPlasmaMatrix) <- runsPlasma

for (i in 1:length(globalStrep["RSE",])) {
  RSErunsPlasmaMatrix[strsplit(strsplit(rownames(as.data.frame(globalStrep["RSE",])), '_')[[i]], " ")[[1]], 
                 strsplit(strsplit(rownames(as.data.frame(globalStrep["RSE",])), '_')[[i]], " ")[[2]] ] <- globalStrep["RSE",i]
  
  RSErunsPlasmaMatrix[strsplit(strsplit(rownames(as.data.frame(globalStrep["RSE",])), '_')[[i]], " ")[[2]], 
                 strsplit(strsplit(rownames(as.data.frame(globalStrep["RSE",])), '_')[[i]], " ")[[1]] ] <- globalStrep["RSE",i]
}

## Build a global tree for 24 runs (plasma) based on LOESS prediction
DistMatPlasma = as.dist(RSErunsPlasmaMatrix, diag = TRUE)
ClusteredPlasma <- hclust(DistMatPlasma, method = "average")
library(ape)
class(ClusteredPlasma) # must be hclust class
my_tree <- as.phylo(ClusteredPlasma) 
plot(my_tree)


