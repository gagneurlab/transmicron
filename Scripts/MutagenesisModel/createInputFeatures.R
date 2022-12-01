##################################################
# description
##################################################

# this script creates GRanges object of the chromatin and genic features used in the mutagenesis model
# depending on user preferences, the table contains 1) only the pre-defined features, 2) only custom, user-supplied features,
# or 3) a combination of both.

################################################
# load packages
################################################
suppressPackageStartupMessages({
library(GenomicFeatures)
library(data.table)
library(regioneR)
library(GenomicRanges)
})
#load genome
#Tx <- makeTxDbFromGFF("/s/project/transposon/Output/PublishedPipeline/Input/Annotations/Genes/mm10.refGene.gtf.gz", format = "gtf")
Tx <- makeTxDbFromGFF(snakemake@params[["REFSEQgtf"]], format = "gtf")


################################################
# create pre-defined features
################################################

#prepare pre-defined features
if(snakemake@wildcards[["mutaFeatures"]]!="OnlyCustomFeatures"){
  
  #read pre-supplied chromatin information
  DNAse_mESC <- fread(snakemake@params[["DNAse_mESC"]])
  ATAC_mESC <- fread(snakemake@params[["ATAC_mESC"]])
  
  #conver to GRanges object
  DNAse_mESC <-toGRanges(DNAse_mESC)
  ATAC_mESC <- toGRanges(ATAC_mESC)

  #genic features
  genes <- sort(genes(Tx))
  TSS <- resize(transcripts(Tx), width=1, fix='start')
  TTS <- resize(transcripts(Tx), width=1, fix='end')
  
  
  FeatureList <- list(DNAse_mESC,
                     ATAC_mESC,
                     genes,
                     TSS,
                     TTS)
  
  names(FeatureList) <- c("DNAse_mESC",
                         "ATAC_mESC",
                         "genes",
                         "TSS",
                         "TTS")
} 

################################################
# create user-defined features
################################################

#if custom features are used, read BED-files and covert to GRanges
if(snakemake@wildcards[["mutaFeatures"]]=="OnlyCustomFeatures"){
FeatureList <- lapply(snakemake@params[["customFeatures"]], fread) #read in list of bed files supplied by used
FeatureList <- lapply(FeatureList, toGRanges) #convert them to GRanges objects
}

################################################
# combined custom and pre-defined feautures 
################################################

#if user wants to add their features to ours, read their BEDs and append to our table
if(snakemake@wildcards[["mutaFeatures"]]=="AddCustomFeatures"){
  
    CustomFeatures <- lapply(snakemake@params[["customFeatures"]], fread) #read in list of bed files supplied by used
    CustomFeatures <- lapply(CustomFeatures, toGRanges) #convert them to GRanges objects
    
    #combine the custom features with the pre-supplied features
    FeatureList <- c(FeatureList,CustomFeatures)
}

################################################
# save output
################################################ 

saveRDS(FeatureList,snakemake@output[["featureList"]])

