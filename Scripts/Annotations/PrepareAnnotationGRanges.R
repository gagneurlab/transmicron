#########################
# description
#########################

# this script prepares prepares a GRanges object of the target annoations

#########################
# load packages
#########################

library(regioneR)
library(data.table)
library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)
library(org.Mm.eg.db)
library(GenomicFeatures)

######################
# 1. custom annotation
######################

if(snakemake@wildcards["annotation"]=="custom"){
  AnnotationGRanges <- fread(snakemake@params["CustomAnnotation"])
  
  #convert input BED to GRanges object
  AnnotationGRanges <- toGRanges(AnnotationGRanges)
  
  #assign feature name:
  
  #if users supply column "name", use this column
  if(length(AnnotationGRanges$name)==length(AnnotationGRanges)){
    names(AnnotationGRanges) <- AnnotationGRanges$name
  } else{ #if no name column supplied: chr_start_end
    names(AnnotationGRanges) <- paste(as.vector(seqnames(AnnotationGRanges)),
                                      start(AnnotationGRanges),
                                      end(AnnotationGRanges), sep="_")
  }
}

##############################
# 2. Prepared annotation1: mm 10 REFseq genes (no promoters)
##############################
if(snakemake@wildcards["annotation"]=="genes"){
  #create TXDB object from pre-supplied gtf-file for refseq genes 
  #Tx <- makeTxDbFromGFF("/s/project/transposon/Output/PublishedPipeline/Input/Annotations/mm10.refGene.gtf.gz", format = "gtf")
  Tx <- makeTxDbFromGFF(snakemake@params[["REFSEQgtf"]], format = "gtf")
  genes <- genes(Tx)
  AnnotationGRanges<-split(genes, genes$gene_id )
}

##############################
# 3. Prepared annotation2: genomic bins from 10 - 150kb  (mm 10 genome)
##############################

if(grepl("kb",unlist(snakemake@wildcards["annotation"]) , fixed = TRUE)==TRUE){
  kb<-as.numeric(unlist(strsplit(unlist(snakemake@wildcards["annotation"]), "kb")))
  
  #function to bin the genome in bins of 10kb * X
  tile_genome <- function(kb){
    tiles <- tileGenome(seqinfo(Mmusculus), tilewidth=kb*10000,cut.last.tile.in.chrom = TRUE)
    nonstandard_chromosome <- which(as.vector(seqnames(tiles))%in%seqnames(Mmusculus)[1:21]==FALSE)
    tiles <- tiles[-nonstandard_chromosome]  
    return(tiles)
  }
  
  AnnotationGRanges <- tile_genome(kb)
  
  #create names of the bins: e.g.: chr1_10000 (chromosome, _, start of bin)
  names(AnnotationGRanges) <- paste(as.vector(seqnames(AnnotationGRanges)), #chromosome
                                    as.vector(start(AnnotationGRanges)), sep="_") #start of bin
}

##############################
# save output
##############################


saveRDS(AnnotationGRanges, snakemake@output[["AnnotationGRanges"]])