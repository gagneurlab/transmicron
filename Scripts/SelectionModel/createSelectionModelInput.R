##################################################
# description
##################################################

#this script creates the input for the selection model

##################################################
# #load packages
##################################################
#load packages
suppressPackageStartupMessages({
library(data.table)
library(GenomicRanges)
library(dplyr)
})
##################################################
# #read input
##################################################

#GRanges object of all insertions in the user supplied BED-file
InsertionLocations <- readRDS(snakemake@input[["Granges_insertions"]])

#read matrix with insertion rates and chromosome for each genetic element
annotationMatrix <- readRDS(snakemake@input[["annotationMatrix"]]) #this is computed separately / precomputed

#Granges object of the target annotation (e.g. all genes)
annotationGRanges <- readRDS(snakemake@input[["AnnotationGRanges"]])

##################################################
#get observed insertions per sample
##################################################

#function to count insertions per sample, per gene 
createInsertionCountMatrix <- function(InsertionLocations,annotationGRanges){
  
  #get vector of unique sample names
  samples<-unique(as.character(as.data.frame(InsertionLocations)$TumorID))
  matrixReal<-matrix(nrow=length(samples), ncol=length(annotationGRanges)) #prepare empty result matrix

  # UNIQUE counts by genetic element, by mouse
    for(i in 1:length(samples)){ #iterate over sample names
      
      #step 1: subset for insertions of each sample
      chosen<-which(as.data.frame(InsertionLocations)$TumorID==samples[i])
      grangessample<-InsertionLocations[chosen]
      
      #step 2: overlap insertions with all genes and store results in matrix
      overlaps<-countOverlaps( annotationGRanges,grangessample, ignore.strand=TRUE)
      matrixReal[i,]<-as.numeric(overlaps)
    }
  
    colnames(matrixReal)<-as.character(names(overlaps)) #names of columns: name of gene / feature 
    matrixReal<-as.data.table(cbind(samples,matrixReal)) #assign column with sample names
    count_matrix_molten<-as.data.frame(melt(matrixReal, id.vars = "samples")) #melt matrix
    colnames(count_matrix_molten)<-c("samples", "feature_name", "insertion_count_ij")
    return(count_matrix_molten)
}

##################################################
#function to join samplewise, gene-wise insertion counts with gene-wise insertion rates and chromosomes
##################################################

createModelInput <- function(countMatrix,annotationMatrix){
  
  # join the two matrices
  SelectionModelInput <- left_join(countMatrix,annotationMatrix)
  
  #format the output
  SelectionModelInput$insertion_count_ij<-as.numeric(SelectionModelInput$insertion_count_ij)
  SelectionModelInput$sum_class_probabilities <-as.numeric(SelectionModelInput$sum_class_probabilities)
  SelectionModelInput$samples <-as.character(SelectionModelInput$samples)
  SelectionModelInput$chromosome <-as.character(SelectionModelInput$chromosome)
  SelectionModelInput$feature_name<-as.character(SelectionModelInput$feature_name)
  
  return(SelectionModelInput)
}

countMatrix <- createInsertionCountMatrix(InsertionLocations,annotationGRanges) #create count matrix for genes
ModelInput <-createModelInput(countMatrix,annotationMatrix)

##########################
#save output
##########################  

saveRDS(ModelInput, snakemake@output[["SelectionModelInput"]])


