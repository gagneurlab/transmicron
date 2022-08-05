#########################
# description
#########################

# this script prepares the annotation defined by users 

#########################
# load packages
#########################

library(regioneR)
library(data.table)
library(BSgenome.Mmusculus.UCSC.mm10)
library(dplyr)
library(org.Mm.eg.db)

#########################
# define functions
#########################


#function sum insertion rates for all single TTAA / TA sites within the genetic feature
overlapforeachgene<-function(x,one_bp_ranges_negatives){
  
  temp<-as.data.frame(findOverlaps(x,one_bp_ranges_negatives, ignore.strand=TRUE))
  rows<-as.data.frame(temp)$subjectHits
  temp$probs<-unlist(InsertionRatesGenome[rows,])
  aggregate<-aggregate(. ~ queryHits  , data=temp, sum)
  sumvector<-rep(0, length(x))
  sumvector[aggregate$queryHits]<-aggregate$probs
  return(sumvector)
}

#function to create the matrix

prepareAnnotation <- function(AnnotationGRanges, one_bp_ranges_negatives, InsertionRatesGenome){ #@ata: genetic element is a GrangesList object, each element one place in the list, name of the feature = names of the list elements
  #sum the single insertion rate in every gene
  sum_class_probabilites<-overlapforeachgene(AnnotationGRanges, one_bp_ranges_negatives)
  sum_class_probabilites<-as.data.frame(cbind(names(AnnotationGRanges), sum_class_probabilites))
  colnames(sum_class_probabilites)<-c("feature_name", "sum_class_probabilities")
  
  #prepate matrix with: 1) name of the gene, 2) chromosome, 3) gene-wise insertion rate
  chromosomes<-data.table(feature_name=as.character(unlist(names(AnnotationGRanges))), #names of each gene
                          chromosome=as.character(as.data.frame(AnnotationGRanges)$seqnames)) #chromosome of each gene
  
  preparedInputMatrix <- left_join(chromosomes, sum_class_probabilites)
  return(preparedInputMatrix)
}


######################
#read input
######################

#matrix with insertion rates
#InsertionRatesGenome <- fread("/s/project/transposon/Output/Snakemake/DLBCL_PB/sequence_DNAsemESCSRX1452763_AtacseqmESC_distancetoGenerelative_distancettsrelative_distancetssrelative/InsertionRatesGenome.csv")
InsertionRatesGenome <- fread(snakemake@input[["InsertionRatesGenome"]])

#granges object with the positions of all TTAA / TA sites in the genome
#one_bp_ranges_negatives <- readRDS("/s/project/transposon/Output/Snakemake/DLBCL_PB/one_bp_ranges_negatives.RData")
one_bp_ranges_negatives <- readRDS(snakemake@input[["one_bp_ranges_negatives"]])

AnnotationGRanges <- readRDS(snakemake@input[["AnnotationGRanges"]])


######################
# create matrix
######################

AnnotationMatrix <- prepareAnnotation(AnnotationGRanges, 
                                      one_bp_ranges_negatives=one_bp_ranges_negatives,
                                     InsertionRatesGenome=whole_genome_predictions)
######################
# save output
######################
#save.image("annotationMatrix.RData")
saveRDS(AnnotationMatrix, snakemake@output[["AnnotationMatrix"]])


# 
# 
# gtf <- fread("/s/project/transposon/Output/PublishedPipeline/Input/mm10.refGene.gtf.gz")
# TXFromGTF <- makeTxDbFromGFF("/s/project/transposon/Output/PublishedPipeline/Input/mm10.refGene.gtf.gz", format = "gtf")
# genesfromgtf <- genes(TXFromGTF)
# genesfromgtf$gene_id
# AnnotationGRanges<-split(genesfromgtf, genesfromgtf$gene_id )
# plyr::count(as.vector(seqnames(genesfromgtf)))
# genesfromgtf
