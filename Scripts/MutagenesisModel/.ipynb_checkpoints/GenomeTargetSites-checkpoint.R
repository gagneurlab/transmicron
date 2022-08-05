##################################################
# description
##################################################

# this script creates a table of all TTAA / TA sites (depending on the transposon system) in the genome
# the table contains feature information used by the mutagenesis model (sequence, chromatin etc.)

##################################################
# load packages
##################################################

library(GenomicRanges)
library(GenomicFeatures)
library(tidymodels)
library(caret)
library(dplyr)
library(regioneR)
library(data.table)
library(BSgenome.Mmusculus.UCSC.mm10)
library(rhdf5)
library(stringr)

#prepare HDF5 session; might fail otherwise
h5closeAll()
h5disableFileLocking()

##################################################
# assign input and parameters
##################################################

#source functions
#source("/data/nasif12/home_if12/bredthau/workspace/transposon/publishedPipeline/MutagenesisModel/define_functions.R"
source(snakemake@params[["InputFunctions"]])

#assign genome
Tx <- makeTxDbFromGFF(snakemake@params[["REFSEQgtf"]], format = "gtf")
genomemouse <- BSgenome.Mmusculus.UCSC.mm10

#assign transposon system
transposon_system <- snakemake@wildcards["transposonSystem"]    

if(transposon_system=="PB"){
  matching_pattern <- "TTAA"
} 
if(transposon_system=="SB"){
  matching_pattern <- "TA"
} 

#determines the size of the sequence context around the insertion size used for mutagenesis correction
window_size<- 10

#read in a a list of features (granges objects)
FeatureList <- readRDS(snakemake@input[["featureList"]])


##################################################
# get TTAA / TA sites in the genome and create window around them
##################################################

#create object containing all TTAA / TA sites in the genome
granges_pattern <- vmatchPattern(matching_pattern,genomemouse)

#remove non-standard chromosomes
granges_pattern<- granges_pattern[-which(as.vector(seqnames(granges_pattern))%in%
                                           (seqnames(genomemouse)[1:21])==FALSE)]


#reduce to one base pair insertion sites(TA motif -> T only)
one_bp_ranges_genome<- get_insertion_site_fromPattern(granges_pattern, matching_pattern) # function is sourced from snakemake@params["InputFunctions"]

#################################################
# create sequence windows around TTAA / TA sites
#################################################

#retrieve a window of -/+ 10bp around the TTAA / TA sites
window_ranges_genome <- retrieve_sequence_window (one_bp_ranges_genome, window_size)  # function is sourced from snakemake@params["InputFunctions"]


#remove locations where sequence context cannot be properly retrieved (e.g. site too close to chromosome boundary)
trimmed <- trim(window_ranges_genome)
removed <- which(window_ranges_genome%in%trimmed==FALSE)  
if(length(removed)!=0){
  window_ranges_genome <- window_ranges_genome[-removed]
  one_bp_ranges_genome <- one_bp_ranges_genome[-removed]
}

wrong_width<-which(width(window_ranges_genome)!=window_size*2+1)
if(length(wrong_width)!=0){
    window_ranges_genome<-window_ranges_genome[-wrong_width]
    one_bp_ranges_genome<-one_bp_ranges_genome[-wrong_width]
}

##################################################
# create feature matrix (negatives, all TTAA / TA sites in the genome)
##################################################

#build the input matrix for the mutagenesis model (iteratively over chromosomes to reduce memory requirements)
for (i in 1:21){
  
  #assign chromosome name i
  seqname <- seqnames(genomemouse)[i]
  
  #get TTAA / TA sites on chromsome i
  one_bp_ranges_chromosome <- one_bp_ranges_genome[which(seqnames(one_bp_ranges_genome)==seqname)]
  window_ranges_chromosome <- window_ranges_genome[which(seqnames(window_ranges_genome)==seqname)] #windows around TTAA / TA sites on chromosome i
  
  #create feature matrices
  combined_input_chromosome <- create_features(one_bp_ranges_chromosome, window_ranges_chromosome, FeatureList) #function from script sourced above

  #remove columns with only 1 level (e.g. TTAA, TA motif around the insertion site); one-hot-enconding fails otherwise
  oneBaseOnly <- which(apply(combined_input_chromosome,2,uniqueN)==1)
  if (length(oneBaseOnly)!=0){
    combined_input_chromosome <- combined_input_chromosome[,-oneBaseOnly]
  }
  
  #remove sequences with non TCGA base (N)
  remove_N <- unique(which(combined_input_chromosome=="N", arr.ind = TRUE)[,1])
  if(length(remove_N)!=0){
    combined_input_chromosome <- combined_input_chromosome[-remove_N,]
    remove <- which(one_bp_ranges_genome %in% one_bp_ranges_chromosome[remove_N])
    one_bp_ranges_chromosome <- one_bp_ranges_chromosome[-remove_N,]
    one_bp_ranges_genome <- one_bp_ranges_genome[-remove] #remove these from the object containing all TTAA / TA sites in the genome (which is used later)
    window_ranges_genome <- window_ranges_genome[-remove]
  }
  
  ###### one hot encode the non-numeric variables  
  # for the first chromosome, create and safe a receipe to use for other inputs (to ensure consistent one-hot-encoding)
  if(i ==1){ 
    formula <- recipe(~ .,
                data = combined_input_chromosome)
    
    receipe <- formula %>%
      step_dummy(all_nominal_predictors()) %>%
      prep(training = combined_input_chromosome)
  
    saveRDS(receipe, snakemake@output[["one_hot_encoder"]])
  }
  
  #use the receipe to create one-hot-encoded matrix
  combined_input_chromosome <- bake(receipe, new_data = combined_input_chromosome)
  
  #remove NAs
  contains_NA<-unique(which(is.na(combined_input_chromosome), arr.ind=TRUE)[,1])
  if(length(contains_NA)!=0){
    combined_input_chromosome <-combined_input_chromosome[-contains_NA,]
    remove <- which(one_bp_ranges_genome %in% one_bp_ranges_chromosome[contains_NA])
    one_bp_ranges_genome <- one_bp_ranges_genome[-remove]#remove these from the object containing all TTAA / TA sites in the genome (which is used later)
    window_ranges_genome <- window_ranges_genome[-remove]
  }

  #save colnames to write to hdf5 dataset later
  colnames <- colnames(combined_input_chromosome)
  #transpose the output; we read it into python for prediction. doing that automaitcally transposes the matrix
  combined_input_chromosome <- t(as.matrix(combined_input_chromosome))


  #save chromosomewise output to hdf5-file
  try(h5createFile(snakemake@output[["step1_input_negatives_whole_genome"]]))
  try(h5delete(snakemake@output[["step1_input_negatives_whole_genome"]], paste("chromsome",str_pad(i, 2, pad = "0"),sep="_")))
  try(h5createDataset(snakemake@output[["step1_input_negatives_whole_genome"]], paste("chromsome",str_pad(i, 2, pad = "0"),sep="_"),  level=9,chunk = c(dim(combined_input_chromosome)[1],1000000)))
  h5write(combined_input_chromosome, snakemake@output[["step1_input_negatives_whole_genome"]], paste("chromsome",str_pad(i, 2, pad = "0"),sep="_"))

  #write colnames as attribute
  file=H5Fopen(snakemake@output[["step1_input_negatives_whole_genome"]])
  did <- H5Dopen(file,paste("chromsome",str_pad(i, 2, pad = "0"),sep="_"))
  h5writeAttribute(did, attr=colnames(paste("chromsome",str_pad(i, 2, pad = "0"),sep="_")),name="colnames") 
}

##################################################
# create random subset of whole genome matrix to limit memory requirements in future calcucations
##################################################

#pick 2 million random TTAA / TA sites
random_subset <- sample(1:length(one_bp_ranges_genome),2000000) 
random_subset_onebp <- one_bp_ranges_genome[random_subset,]
random_subset_onewindow <- window_ranges_genome[random_subset,]

#create feature table for the random subset
combined_input_chromosome <- create_features(random_subset_onebp, random_subset_onewindow, FeatureList)

#clean results
remove_N <- NULL
try(remove_N <- which(combined_input_chromosome=="N", arr.ind = TRUE)[,1])
if(length(remove_N)!=0){
  combined_input_chromosome <- combined_input_chromosome[-remove_N,]
  random_subset_onebp <- random_subset_onebp[-remove_N]
  random_subset_onewindow <- random_subset_onewindow[-removeN]
}

#remove columns with only 1 level (applies to TTAA, TA around the insertion site)
combined_input_chromosome <- combined_input_chromosome[,-which(apply(combined_input_chromosome,2,uniqueN)==1)]

#one hot encode the non-numeric variables  
step1_input_negatives_controls <- bake(receipe, new_data = combined_input_chromosome)

#remove NAs
contains_NA<-which(is.na(step1_input_negatives_controls), arr.ind=TRUE)[,1]
if(length(contains_NA)!=0){
  step1_input_negatives_controls <-step1_input_negatives_controls[-contains_NA,]
  random_subset_onebp <- random_subset_onebp[-contains_NA]
  random_subset_onewindow <- random_subset_onewindow[-contains_NA]
  combined_input_chromosome <- combined_input_chromosome[-contains_NA,]
}

##################################################
# save output
##################################################

saveRDS(random_subset_onebp, snakemake@output[["granges_negatives_controls"]])
saveRDS(step1_input_negatives_controls, snakemake@output[["step1_input_negatives_controls"]])
saveRDS(random_subset_onebp, snakemake@output[["granges_negatives_controls"]])
saveRDS(one_bp_ranges_genome, snakemake@output[["one_bp_ranges_negatives"]])

