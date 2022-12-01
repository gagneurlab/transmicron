##################################################
# description
##################################################

# this script reads the BED-file containing the loci of transposon insertions (supplied by user)
# it conducts certain quality checks and converts the insertions to a GRanges object

################################################
#load packages
################################################
suppressPackageStartupMessages({
library(GenomicRanges)
library(data.table)
library(regioneR)
library(BSgenome.Mmusculus.UCSC.mm10)
})
################################################
#read input
################################################

#insertionBED <- read.table("/s/project/transposon/Output/PublishedPipeline/Input/BEDInsertionTesting/cuSCC.BED")

insertionBED <- fread(snakemake@params[["insertionFile"]], colClasses = 'character')

#conversion to GRanges

Granges_insertions <- toGRanges(insertionBED)
names(Granges_insertions) <- 1:length(Granges_insertions)


################################################
# quality check of user-supplied BED-file with insertions
################################################

#check if width is 2 bp (needed later for sequence retrieval)
if(unique(width(Granges_insertions))!=2){
  stop("Error: Please make sure your BED file specifies insertion sites with start and end position (-1, +1)")  
}

#check if tumor IDs are provided
if(length(Granges_insertions$TumorID)!=length(Granges_insertions)){
  stop("Error: Please make sure your you supply a TumorID / Sample ID for each insertion in your BED file. Please name the column 'TumorID'.")  
}

#check if orientations are provided
if(length(Granges_insertions$orientation)!=length(Granges_insertions)){
  stop("Error: Please make sure your you supply an orientation for each insertion (either '-', '+', or '*'). Please name the column 'Orientation'.")  
}

#check if only "-" and "+" are in orientation column
if(FALSE%in%(unique(Granges_insertions$orientation)%in%c("-","+"))){
  stop("Warning: You supplied values other than '-' and '+' in the orientation column of your BED-file") 
}

################################################
#for easier handling, insertion sites are converted to one bp locations
################################################

#assign strand to GRanges object
strand(Granges_insertions) <- Granges_insertions$orientation
Granges_insertions$orientation <- NULL


get_oneBP_insertion_site<-function(Granges_insertions){

  strand<-as.vector(strand(Granges_insertions))
  index_minus<-which(strand=="-")
  index_plus<-which(strand=="+")

  end(Granges_insertions[index_plus]) <- end(Granges_insertions[index_plus]) -1
  start(Granges_insertions[index_minus]) <- start(Granges_insertions[index_minus]) +1

  return(Granges_insertions)
}

Granges_insertions <- get_oneBP_insertion_site(Granges_insertions) 

#remove insertions non-standard chromosomes, if there are any
nonStandardCHR <- NULL
nonStandardCHR <- which(as.vector(seqnames(Granges_insertions))%in%
                          (seqnames(Mmusculus)[1:21])==FALSE)
if(length(nonStandardCHR)!=0){
  Granges_insertions<- Granges_insertions[-nonStandardCHR]
}


################################################
#save output
################################################
saveRDS(Granges_insertions, snakemake@output[["Granges_insertions"]])

