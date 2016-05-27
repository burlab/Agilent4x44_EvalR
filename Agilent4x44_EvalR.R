#####################################################################################################
#####################################################################################################
#                                                                                                   #
# Scripts for Evaluating Agilent4x44 TWO COLOR Microarray Data                                      #
# -------------------------------------------------------------                                     #
#                                                                                                   #
# Version 19 (29.04.2015)                                                                           #                                                                                                 #
#                                                                                                   #
#                                                                                                   #
#                                                                                                   #
# R scripts for: - Import of target matrix (table containing sample information)                    #
#                - Import of raw data (Data obtained from Agilent Feature extraction)               #
#                - Background substraction and normalization of raw data                            #
#                - QC plots                                                                         #  
#                - Subsetting of data based on parameters defined in the target matrix              #
#                - Replicate probe averaging, removal of control feature data                       #
#                - Adding gene symbols to probe IDs                                                 #
#                - LIMMA analyses (including definition of linear models, contrast matrices etc)    #
#                - Display and plotting of LIMMA results (heat maps, Venn diagrams etc)             #
#                - Spacial plots of microarrays to identify spacial artefacts                       #
#                - More                                                                             #
#                                                                                                   #
#                 Using R and R/Bioconductor packages                                               #
#                                                                                                   #
#                                                                                                   #
# © Bo Burla, University Hospital Zürich, Switzerland    (bo.burla@gmail.com)                       #
#                                                                                                   #      
# 2013-2015                                                                                         #
#                                                                                                   #
#####################################################################################################
#####################################################################################################


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Initialization of libraries and functions (load before start)
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



#====================================================================================================
# Libraries
#====================================================================================================

# Load Libaries
# R packages
library("gplots")
library("convert")

# R/Bioconductor packages (see 'http://www.bioconductor.org/install/' on how to install)
library(limma)
library(HsAgilentDesign026652.db)
#library(hgug4112a.db)

library(marray)
library(annotate)
library(arrayQualityMetrics)
library(Biobase)
library(genefilter)

library(topGO)
library(GO.db)
library(RankProd)

#library(vsn)



rm(list=ls())    # delete all defined variables, functions etc from the memory!
old.par <- par(no.readonly = TRUE)   #


#####################################################################################################
#                                                                                                   #
# Global functions for scripts "Evaluating Agilent4x44 Microarray Data"                             #
#                                                                                                   #
#####################################################################################################

#====================================================================================================
# Plot MA plots of all arrays
#     TODO: Plotting of plot titles ("main") is not working
#====================================================================================================

plotMA_AllArrays <- function(data, arrayCount) {
  o.par <- par(no.readonly = TRUE)
  on.exit(par(o.par))                 #reset graphic parameters when leaving function
  
  par(mfrow=c(ceiling(arrayCount/4),4), # plot 4 charts per row
      oma = c(1, 0, 1 , 0),  # ? Check
      mar = c(0, 0, 0, 0),   # ? Check
      mgp = c(1, 0, 0))      # ? Check
  
  
  for (i in 1:(arrayCount)){
    limma::plotMA(data, array=i, legend=FALSE,cex=0.2, tck=0.03)
    cat(i) 
  }
  return(cat("Done!"))
}

#====================================================================================================
# Plot R vs G plots of all arrays 
#     TODO: Plotting of plot titles ("main") is not working
#====================================================================================================

plotRG_AllArrays <- function(data, arrayCount) {
  o.par <- par(no.readonly = TRUE)
  on.exit(par(o.par))                 #reset graphic parameters when leaving function
  par(mfrow=c(ceiling(arrayCount/3),3),
      oma = c(1, 0, 1 , 0),  # ? Check
      mar = c(0, 0, 0, 0),   # ? Check
      mgp = c(1, 0, 0))      # ? Check
  spotColors <- attr(data$genes$Status, "col")
  spotTypeCol=spotColors[ordered(data$genes$Status, levels= attr(data$genes$Status, "values"))]
  for (i in 1:arrayCount){
    plot(log2(data$G[,i]),log2(data$R[,i]), main = colnames(data[,i]), col=spotTypeCol)
    cat(i)
  }
  par(old.par)
}

#====================================================================================================
# Spacial Plots of Arrays
#
# Agilent Feature Extraction software doesn't print rows for blank spots and the like, so the block has
# missing rows, and this causes a problem because the limma functions assume complete blocks
# https://stat.ethz.ch/pipermail/bioconductor/2007-January/015532.html
#====================================================================================================

allArrayplot <- function(data, arrayCount, channel) {
  o.par <- par(no.readonly = TRUE)
  on.exit(par(o.par))                 #reset graphic parameters when leaving function
  
  par(mfrow=c(ceiling(arrayCount/4),4),
      oma = c(2, 2, 1, 1),  # rows of text at the outer margins
      mar = c(1, 1, 0, 0),  # rows to separate plots
      mgp = c(1, 0, 0))     # axis and tick labels distance
  
  data$genes$Block <- 1
  names(data$genes)[2] <- "Column"
  data$printer <- getLayout(data$genes)
  
  # Fill missing rows
  r <- data$genes$Row
  c <- data$genes$Col
  nr <- max(r)
  nc <- max(c)
  y <- rep(NA,nr*nc)
  i <- (r-1)*nc+c
  
  
  for (z in 1:(arrayCount)){
    cat(z)
    if (channel=="R"){
      y[i] <- log2(data$R[,z])
    }
    else if (channel=="G"){
      y[i] <- log2(data$G[,z])
    }
    else if (channel=="RG"){
      y[i] <- log2(data$R[,z])/log2(data$G[,z])
    }
    else if (channel=="E"){
      y[i] <- log2(data$E[,z])
    }
    
    imageplot(y,data$printer, low="white", high="red", ncolors=64, zerocenter=F)  #zlim could be used to define range c(min, max)
  }
  par(old.par)
  return(NULL)
}

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Select and import data
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# ADJUST PATH of Agilent raw files
setwd("H:/USZ/Projects/Microarrays HFE/HFE_Microarrays_2013/Series_2_2012_TwoColor")  # Path of data

#====================================================================================================
# Define samples to be imported (target matrix)
#====================================================================================================


# Read target matrix (contains array names, file names, additional data/parameters for analysis)
# Adjust path if necessary
targetsALL=readTargets("../TargetMatrix/TargetMatrix_HFE_Microarray_SubjectData_V10_CompleteDataset.txt")

# Define which array to import: This can be done here, or later just before the "Between-Array Normalization"
#   to save time by avoid re-importing all data again when testing different sample subsets

targetsALL = targetsALL[targetsALL$FileName_Series2 != "",] # Select only Array Files existing for Series 2
targetsSelected = targetsALL
# Select  arrays according to their array QC results (0=bad, 1=questionable/ok, 2= good)
targetsSelected = targetsSelected[targetsSelected$ArrayQC_Series2 > 0,] 
targetsSelected = targetsSelected[targetsSelected$Priority == 1,]      # OPTIONAL: Select only Samples with Priority-Level 1
#targetsSelected = targetsSelected[!(targetsSelected$RNAID == "11"),]   # Exclude a specififc array

print(targetsSelected)

#====================================================================================================
# Rename sample names with more informative names (ID_Genotype_Age_Sex_Ferritin_Severity) 
#====================================================================================================

newArrayNames =paste(targetsSelected$RNAID,"_", 
                targetsSelected$Genotype,"_", 
                targetsSelected$Age,
                targetsSelected$Sex,
                "_Fer",targetsSelected$Ferritin,
                "-", targetsSelected$SeverityHH, 
                sep="")
targetsSelected$SampleDescription = newArrayNames

#====================================================================================================
# Importing selected data from array data files (files from Agilent Future Extraction Software)
#====================================================================================================


RGraw <- read.maimages(targetsSelected$FileName_Series2, 
                       source="agilent", 
                       annotation = c("Row", "Col","FeatureNum", "SystematicName", "ControlType","ProbeName"), 
                       names = newArrayNames)
#
# or import PROCESSED VALUES (rProcessedSignal and gProcessedSignal) of files from Agilent Future Extraction Software 
# RGb <- read.maimages(targetsALL, source="agilent", columns= list( R = "rProcessedSignal", G = "gProcessedSignal") )
RGraw$targets <- targetsSelected

#Add spot type definitions (Gene, Spike-in, Corner, Ctrl etc = QC spots of arrays) to the values 
spottypes <- readSpotTypes("SpotTypes.txt")
RGraw$genes$Status <- controlStatus(spottypes, RGraw)

# Set a global variable storing number of arrays imported
arrayCount <- nrow(RGraw$targets)

#====================================================================================================
# QC plots of important raw data
#====================================================================================================

# Plots spacial plots of all arrays
allArrayplot(RGraw, arrayCount, "R")    # R channel
allArrayplot(RGraw, arrayCount, "G")    # G channel
allArrayplot(RGraw, arrayCount, "RG")   # R/G ratio


## Plots R vs. G raw intensity value plots all arrays. Legend to spot colors see spot type definition before
#plotRG_AllArrays(RGraw, arrayCount)

# Plot background signals of green and red channels
boxplot(data.frame(log2(RGraw$Gb)),main="Green background", ylab="Log Intensity") 
boxplot(data.frame(log2(RGraw$Rb)),main="Red background", ylab="Log Intensity")

# Plot  signals of green and red channels
boxplot(data.frame(log2(RGraw$G)),main="Green Channel", ylab="Log Intensity") 
boxplot(data.frame(log2(RGraw$R)),main="Red Channel", ylab="Log Intensity")

# Plot distributions of intensities of raw signals 
limma::plotDensities(RGraw)
limma::plotDensities(RGraw, group = RGraw$target$Arrayslide_Series2, 
                     col = c("orange", "red", "green", "chartreuse4", "skyblue2", "blue"), log=TRUE)
plot(density(log(RGraw$G)))
     
# MA plots of raw data from all imported arrays
plotMA_AllArrays(RGraw, arrayCount)
     
# MA plots of raw data from one array
limma::plotMA(RGraw, array=3, main = "MA plot of raw signal")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Array data pre-processing:  - Background correction, normalization, replicate averaging
#                             - Quality Control (QC) plots
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#====================================================================================================
# Background correction 
#====================================================================================================
RGb <- backgroundCorrect(RGraw, method="normexp", offset=50)   #offset = 50 to avoid negative values

# QC Plots of background substrated data
# ======================================
# Plot distributions of intensities of background corrected signals 
limma::plotDensities(RGb) ## plotDensities(RGb) only displays red channel...BUG
limma::plotDensities(RGb, group = RGb$target$RNAID, col=c(rep("blue",20), "red", rep("green", 3)), log=TRUE)
limma::plotDensities(RGb, group = RGb$target$RNAID, log=TRUE)
limma::plotDensities(RGb, group = RGb$target$RNAID, col=rainbow(23), log=TRUE)

# MA plots of all imported and background-corrected arrays
limma::plotMA(RGb, array=1, main="MA plot of Array 1", ylim=c(-4, 4))  #MA plot of 1 array
plotMA_AllArrays(RGb, arrayCount)  #MA plot of all arrays (load function in the header of this file)

#==================================================================================
# Intra-array normalization by LOESS etc
#====================================================================================================

MAnormwithin <- normalizeWithinArrays(RGb, method="loess")

## TODO: Normalize after downweight and upweight specific controls or features 
# w <- modifyWeights(array(1,dim(RG)), RG$genes$Status, c("Spikein-E1A","NegCtrl", "Spikein-ETG", "Spikein-ERCC", "Spikein-DCP"), c(0,2,2,2,2))
# MA <- normalizeWithinArrays(RGb, weights=w, method="loess")
# QC Plots of background substrated data

# ======================================
limma::plotMA(MAnormwithin, array=3, main = "MA plot of BKG substracted and WITHINARRAY-NORMALISED signals")
plotMA_AllArrays(MAnormwithin, arrayCount)  # MA plots of all background-corrected and normalized arrays

limma::plotDensities(MAnormwithin)
limma::plotDensities(MAnormwithin, group = MAnormwithin$target$Arrayslide_Series2, 
                     col = c("orange", "red", "green", "chartreuse4", "skyblue2", "blue"), log=TRUE)
limma::plotDensities(MAnormwithin, group = MAnormwithin$target$RNAID, col=rainbow(23), log=TRUE)
limma::plotDensities(MAnormwithin, group = MAnormwithin$target$RNAID, log=TRUE)

# Create arrayQualityMetrics report -> will be created in the subfoder "QualityMetrics_Reports"
arrayQualityMetrics(expressionset = MAnormwithin,outdir="QualityMetrics_Reports_MAnormwithin",force=TRUE)

#====================================================================================================
# Normalization by Variance Stabilization  (vsn package) 
#
# Alternative method for normalization of Arrays: Instead of background correction, intra and inter-array
# normalization
#
# ATTENTION: This script is may be incorrect
#====================================================================================================

# MAvsn <- normalizeVSN(RGraw)
# meanSdPlot(justvsn(RGraw))
# 
# meanSdPlot(MAvsn)
# plotMA_AllArrays(MAvsn, arrayCount)  # MA plots of all background-corrected and normalized arrays


#====================================================================================================
# Subsetting of data according to parameters defined in the target matrix 
#  --> Can (should) be done before Between-Array Normalization ?
#====================================================================================================

MAnws <- MAnormwithin[,]
# MAnws <- MAnormwithin[, MAnormwithin$targets$Priority =="1" ] # Select only samples with priority level = 1
# MAnws <- MAnormwithin[, MAnormwithin$targets$Sex =="M" ]

#====================================================================================================
# Between-Array  normalization
# ----------------------------
# (TOCHECK: Really required or may it even be detrimental to the data?)
#====================================================================================================

MAnbs <- normalizeBetweenArrays(MAnws, method="Aquantile")  # or Gquantile

# QC plot
plotMA_AllArrays(MAnbs, arrayCount)  # MA plots of all background-corrected and normalized arrays
limma::plotDensities(MAnbs, group = RGraw$target$Arrayslide_Series2, 
                     col = c("orange", "red", "green", "chartreuse4", "skyblue2", "blue"), log=TRUE)
limma::plotDensities(MAnbs, group = RGraw$target$RNAID, col=rainbow(23), log=TRUE)

# Create arrayQualityMetrics report -> will be created in the subfoder "QualityMetrics_Reports"
arrayQualityMetrics(expressionset = MAnbs,outdir="QualityMetrics_Reports_MAnbs",force=TRUE)

#====================================================================================================
# Replicate probe averaging, remove array control features, adding gene symbols 
#====================================================================================================

# Adding gene symbols to the probes for better identification of genes (Proben names are array-specific)
rownames(MAnbs$M) = MAnbs$genes$ProbeName
MAnbs$genes$Symbol <- getSYMBOL(MAnbs$genes$ProbeName, "HsAgilentDesign026652.db")
MAnbs$genes$Symbol[MAnbs$genes$Symbol==""] = "NoName"    # Replace empty names with "NoName"

# Remove array control (QC) features (so that they are not considered as measured value in sunsequent analyses)
MAnsample = MAnbs[MAnbs$gene$ControlType == 0,]     # Sample: ControlType = 0         

# Average replicate probes within arrays  (there are multiple probes with the same probe ID accros the array)
MAn = avereps(MAnsample, ID=MAnsample$genes$ProbeName)          
#MAn = avereps(MAn, ID=MAn$genes$SystematicName)   # TOCHECK: Why does this not work ????????????????????????



#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Statistics/Analyses to test for differentially expressed (DE) genes
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#====================================================================================================
# Retrieve experimental LEVELS and values defined in the targets matrix 
#
# Has to be adapted according to the target matrix
#====================================================================================================

fGenotype <- factor(MAn$targets$Genotype)                # FACTORS 
fSex <- factor(MAn$targets$Sex)
vAge<- MAn$targets$Age
vFerritin <- MAn$targets$Ferritin

fSlide <- factor(MAn$targets$Arrayslide_Series2)
fBatch <- factor(MAn$targets$Batch_Series2)
fSeverityHH <- factor(MAn$targets$SeverityHH)
fFerritinLevel <- factor(MAn$targets$Ferritin_Level)
fTransfSatLevel <- factor(MAn$targets$TSAT_Level)
fSolTransfLevel <- factor(MAn$targets$sTfR_Level)

  
#====================================================================================================
# Define Design/Contrast MATRICES 
#====================================================================================================

designALL <- model.matrix(~0+fGenotype)
contrastsALL<- makeContrasts(fGenotypehfe-fGenotypectrl, levels=designALL)

designALL <- model.matrix(~0+fGenotype+fSex)
contrastsALL<- makeContrasts(fGenotypehfe-fGenotypectrl, levels=designALL)

designALL <- model.matrix(~0+fGenotype+fSex+fSlide)
contrastsALL<- makeContrasts(fGenotypehfe-fGenotypectrl, levels=designALL)

designALL <- model.matrix(~0+fGenotype+fSlide)
contrastsALL<- makeContrasts(fGenotypehfe-fGenotypectrl, levels=designALL)

designALL <- model.matrix(~0+fGenotype+fSex+fBatch)
contrastsALL<- makeContrasts(fGenotypehfe-fGenotypectrl, levels=designALL)

designALL <- model.matrix(~0+fGenotype+fSex+fBatch+fSlide)
contrastsALL<- makeContrasts(fGenotypehfe-fGenotypectrl, levels=designALL)

# designALL <- model.matrix(~0+vAge+fSex)
# TOCHECK: How to make contast matrix in this case....?   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#   or 
# designALL <- model.matrix(~0+fGenotype)
# contrastsALL<- makeContrasts(fGenotypehfe-fGenotypectrl, levels=designALL)
#   or 
designALL <- model.matrix(~0+fSex)
contrastsALL<- makeContrasts(fSexF-fSexM, levels=designALL)

designALL <- model.matrix(~0+fSex+fSlide)
contrastsALL<- makeContrasts(fSexF-fSexM, levels=designALL)

designALL <- model.matrix(~0+fSex+fBatch)
contrastsALL<- makeContrasts(fSexF-fSexM, levels=designALL)

# designALL <- model.matrix(~0+fGenotype*fSex)
# colnames(designALL) = c("ctrl", "hfe", "sexM", "hfeSex")
# contrastsALL<- makeContrasts(ctrl-hfe, (ctrl-hfeSex) - (hfe-hfeSex), levels=designALL)
# colnames(designALL) <- c("A", "B", "F", "FF")

# designALL <- model.matrix(~0+fGenotype+fFerritinLevel+fSex)
# contrastsALL<- makeContrasts(fGenotypehfe-fGenotypectrl, levels=designALL)

# designALL <- model.matrix(~0+fFerritinLevel+fSex)
# contrastsALL<- makeContrasts(fGenotypehfe-fGenotypectrl, levels=designALL)

#====================================================================================================
# Determine and plot Array Weights: Array weights (arayw) can be used later in Limma analyses
#====================================================================================================

arrayw <- arrayWeights(MAn, design=designALL)
barplot(cex.names=0.8,
        arrayw, 
        xlab="", ylab="Weight", col="white", las=2, 
        main="~0+fGenotype+fSlide", 
        names.arg=MAn$targets$SampleName)
abline(h=1, lwd=1, lty=2)


#====================================================================================================
# Linear model fit
#====================================================================================================

fitA   <- lmFit(object=MAn,design=designALL)
fitAcontrast <- contrasts.fit(fitA, contrastsALL)
fitAe <- eBayes(fitAcontrast)

topTableFull <- topTable(fitAe, coef=1, adjust="BH",sort.by = "logFC", confint=TRUE, number=Inf)
#write.table(topTableFull,"clipboard-50000",sep="\t", row.names=F)  # -> can be directly pasted into an EXCEL sheet

# Display distribution of p-values (histogram)
hist(x=fitAe$p.value, freq=TRUE, breaks=20, main="Histogram of p values")

# Create a vulcano plot
volcanoplot(fitAe,coef=1,highlight=5, cex = 0.7)

#====================================================================================================
# Decide tests
#====================================================================================================

results <- decideTests(fitAe)
summary(results)
vennDiagram(results)

#====================================================================================================
# Display Results (Toptable) obtained from LIMMA analyses 
#====================================================================================================


# Copy top 100 features (sorted by abs(FC)) to clipboard, only display selected columns (fields) of the TopTable
topTableSubset <- topTableFull[1:100,-c(1:5, 7)]
write.table(topTableSubset,"clipboard",sep="\t", row.names=F)
print(topTableSubset)

# Display top 100 features (features) with adjusted p-values < 0.05 and fold-change of 1.5
topTableSubset <- topTableFull[(topTableFull$adj.P.Val < 0.05) & (abs(topTableFull$logFC) > log2(1.5)),][,-c(1:5, 7)]
write.table(topTableSubset,"clipboard-10000",sep="\t", row.names=F)
nrow(topTableSubset)-1

# Display top 100 features (features) with p-values < 0.05 and fold-change of 1.5
topTableSubset <- topTableFull[(topTableFull$P.Value < 0.05) & (abs(topTableFull$logFC) > log2(1.5)),][,-c(1:5, 7)]
write.table(topTableSubset,"clipboard-10000",sep="\t", row.names=F)
nrow(topTableSubset)-1

print(topTableSubset)

# Display values of a specific gene, e.g. "OASL"
#topTableSubset <- topTableFull[(topTableFull$Symbol %in% "OASL"),][,-c(1:5, 7)]
#print(topTableSubset)


#====================================================================================================
# Display Heatmaps from Subset of Toptable (topTableOutput)
#====================================================================================================


# Subsets the MAn object with the features (genes) in listed in the topTableSubject object
MAnSubset <- MAn[match(topTableSubset$ProbeName,MAn$genes$ProbeName),]
probeNames <- rownames(MAnSubset)

color.map <- function(Genotype) { 
    if (Genotype=="ctrl") "#0000FF" else "#FF0000" 
}  

# Create vector with colors to indicate control and treatment groups
treatmentcolors <- unlist(lapply(MAn$targets$Genotype, color.map)) # colors for the colorbar

par(old.par)
par(oma=c(5,1,1,1))
heatmap.2(2^MAnSubset$M, 
          col=redgreen(75), 
          scale="row", 
          ColSideColors=treatmentcolors,
          key=TRUE, 
          symkey=FALSE, 
          density.info="none", 
          trace="none", 
          cexRow=1, 
          Colv = TRUE, 
          Rowv = TRUE,
          labRow = MAn$genes[match(probeNames,MAn$genes$ProbeName),]$Symbol,    #Rename probe names with gene symbols
          )

#====================================================================================================
# "Manual" search for Differentially Expressed Genes
#
# t-Tests between top 3 highest values from patients vs. all controls for all genes
#====================================================================================================

allControls=MAn[,grep("ctrl", colnames(MAn$M))]$M    #AllControls
topHighThreePatients=t(apply(MAn[,grep("hfe", colnames(MAn$M))]$M,1,function (x) sort(x,decreasing = TRUE)[1:3]))
topLowThreePatients=t(apply(MAn[,grep("hfe", colnames(MAn$M))]$M,1,function (x) sort(x,decreasing = FALSE)[1:3]))

fdataUp= cbind(allControls,topHighThreePatients)
fdataDown= cbind(allControls,topLowThreePatients)

# Test for top 3 up-regulated genes in patients
#tResults = apply(fdata, 1, function (tdat) {t.test(tdat[c(1:9)], tdat[10:12])$p.value})
tResults = rowttests(fdataUp, factor(c(rep("C",9), rep("P",3))))
tResults$Symbol = MAn$genes$Symbol
tResults$ProbeID = MAn$genes$ProbeName
tResults$logFC = apply(fdataUp, 1, function (tdat) {mean(tdat[c(10:12)])}) - apply(fdataUp, 1, function (tdat) {mean(tdat[c(1:9)])})
tResults$adj.p.value = p.adjust(tResults$p.value, method="fdr")

# Test for top 3 down-regulated genes in patients
tResults = rowttests(fdataDown, factor(c(rep("C",9), rep("P",3))))
tResults$Symbol = MAn$genes$Symbol
tResults$ProbeID = MAn$genes$ProbeName
tResults$logFC = apply(fdataDown, 1, function (tdat) {mean(tdat[c(10:12)])}) - apply(fdataDown, 1, function (tdat) {mean(tdat[c(1:9)])})
tResults$adj.p.value = p.adjust(tResults$p.value, method="fdr")


tResultsSortP= tResults[order(tResults$p.value),]
topTableSubsetP <- tResultsSortP[(tResultsSortP$p.value < 0.05) & (abs(tResultsSortP$logFC) > log2(2)),]
write.table(topTableSubsetP,"clipboard-10000",sep="\t", row.names=F)

tResultsSortAdjP= tResults[order(tResults$adj.p.value),]
topTableSubsetAdjP <- tResultsSortAdjP[(tResultsSortAdjP$adj.p.value < 0.05) & (abs(tResultsSortAdjP$logFC) > log2(2)),]
write.table(topTableSubsetAdjP,"clipboard-10000",sep="\t", row.names=F)


#====================================================================================================
# Conduct non-parametric Mann-Whitney test between "Ctrl" and "hfe" groups for each gene
#====================================================================================================

MAn$NPtest$p.values<- sapply(1:nrow(MAn$M), function(i) wilcox.test(MAn$M[i,grep("hfe", colnames(MAn$M))], 
                                                     MAn$M[i,grep("ctrl", colnames(MAn$M))])$p.value)
MAn$NPtest$adj.p.values <- p.adjust(MAn$NPtest$p.values,method="BH")

names(MAn$NPtest$p.values) =  MAn$genes$Symbol
print(sort(MAn$NPtest$p.values, decreasing = FALSE)[1:100])

names(MAn$NPtest$adj.p.values) = MAn$genes$Symbol
print(sort(MAn$NPtest$adj.p.values, decreasing = FALSE)[1:100])

#====================================================================================================
# Find genes with highest variances or differences in min and max expression values across all arrays
#====================================================================================================

diff = apply(MAn$M, 1, sd) / apply(MAn$M, 1, mean)
names(diff)=MAn$genes$Symbol
diffsort= sort(diff, decreasing = TRUE)
diffsort = diffsort[1:150]
print(diffsort)



#====================================================================================================
# Plot expression values of a single gene  (Dot plot and BOX PLOT)
#====================================================================================================

genes <- c("TACSTD2")
maintext = genes

# Obtain expression values of "ctrl" and "hf" samples of a selected gene
eC=2^MAn[match(genes, MAn$genes$Symbol),grep("ctrl", colnames(MAn$M))]$M
eP=2^MAn[match(genes, MAn$genes$Symbol),grep("hfe", colnames(MAn$M))]$M

# probes = c("A_24_P854913")
# maintext = probes
# eC=2^MAn[match(probes, MAn$genes$ProbeName),grep("ctrl", colnames(MAn$M))]$M
# eP=2^MAn[match(probes, MAn$genes$ProbeName),grep("hfe", colnames(MAn$M))]$M

# Dot plot of the expression values of the with two groups
  
dotplot.2 <- function(valuesA, valuesB, textA, textB, genename) { 
    plot(pch=19, col=rainbow(9), cex = 2.5,jitter(rep(1, times=length(valuesA),
                    factor=3, amount=21)), 
                    eC, 
                    axes=FALSE, 
                    frame.plot=TRUE , 
                    ylab="Expression level",
                    main = maintext,
                    xlab="", 
                    xlim=c(0.5,2.5),            
                    ylim=c(min(eC,eP,0),max(eC,eP)))     # or set min = 0...
    
    points(pch=19, cex=2.5,col=rainbow(11), jitter(rep(2, times=length(valuesB),factor=1, amount=11)), valuesB)
    axis(2)
    mtext(textA, side = 1, at=1, padj=1)
    mtext(textB, side = 1, at=2, padj=1)
    legend("topleft",pch = c(rep(19, 6)), cex=0.9,col=c(rainbow(9), rainbow(11)),
           legend = c(colnames(MAn[,grep("ctrl", colnames(MAn$M))]$M),colnames(MAn[,grep("hfe", colnames(MAn$M))]$M)))
}


dotplot.2(eC, eP, "Control", "hfe -/-", genes)

    
# Boxplot of the expression values of the with two groups
boxplot(list(eC, eP), outline=TRUE, names=c("Ctrl", "hfe"))


#============================================================================================================================
# topGO:  Gene enrichment analyses
#
#============================================================================================================================

allGeneNames <- MAn$genes$ProbeName
topGeneNames <- topTableSubset$ProbeName
geneList <- factor(as.integer (allGeneNames %in% topGeneNames))
names(geneList) <- allGeneNames

GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList, annot = annFUN.db, affyLib = "hgug4112a.db")
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultFisherTable = GenTable(GOdata, classicFisher = resultFisher, topNodes = 30)
write.table(resultFisherTable,"clipboard-10000",sep="\t", row.names=F)

resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
GenTable(GOdata, classicKS = resultKS, topNodes = 30)
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
GenTable(GOdata, elimKS = resultKS.elim, topNodes = 30)


allRes <- GenTable(GOdata, classicFisher = resultFisher, classicKS = resultKS, elimKS = resultKS.elim, 
                   orderBy = "elimKS", ranksOf = "classicFisher", topNodes = 30)


pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]
gstat <- termStat(GOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
colMap <- function(x) {
  .col <- rep(rev(heat.colors(length(unique(x)))), time = table(x))
  return(.col[match(1:length(x), order(x))])
}
gCol <- colMap(gstat$Significant)
library(Rgraphviz)
plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim",pch = 19, cex = gSize, col = gCol)

showSigOfNodes(GOdata, score(resultKS.elim), firstSigNodes = 3, useInfo = 'all')

#============================================================================================================================
# RankProd:  Rank product non-parametric method (http://www.bioconductor.org/packages/release/bioc/html/RankProd.html) 
#
# PROBLEM: Groups (control vs hfe) are not balanced in sex ratio, therefore many sex-specific genes appear in the list (e.g XIST) 
#
#============================================================================================================================

clGenotypes <- ifelse(targetsSelected$Genotype == "ctrl", 0, 1)     #  Define a vector with groups (only 2 groups permitted)

origin <- rep(x=1, times = ncol(MAn$M))                    #  Define origin/batches of data (here assuming all data from one batch)
# origin <- targetsSelected$Slide
# origin <- targetsSelected$Batch
                   
RPoutGenotype <- RPadvance(MAn$M,clGenotypes,origin,num.perm=100, logged=TRUE,rand=123)
plotRP(RPoutGenotype, cutoff=0.05)
topGenotypeDEgenes <- topGene(RPoutGenotype,cutoff=0.05,method="pfp",logged=TRUE,logbase=2,gene.names=make.unique(MAn$genes$Symbol)) 



