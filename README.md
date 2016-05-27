# Agilent4x44_EvalR
### Agilent Human 4x44K Microarrays Data Analysis Workflow in R


Collection of R scripts to import, pre-process and statistically evaluate Agilent Human 4x44 Two Color microarray data  

  - Import of target matrix (table containing sample information)
  - Import of raw data (Data obtained from Agilent Feature extraction)
  - Background substraction and normalization of raw data
  - QC plots
  - Subsetting of data based on parameters defined in the target matrix
  - Replicate probe averaging, removal of control feature data
  - Adding gene symbols to probe IDs
  - LIMMA analyses (including definition of linear models, contrast matrices etc)
  - Display and plotting of LIMMA results (heat maps, Venn diagrams etc)
  - Spacial plots of microarrays to identify spacial artefacts
  - More

Uses R and R/Bioconductor packages  

© Bo Burla, 2015
