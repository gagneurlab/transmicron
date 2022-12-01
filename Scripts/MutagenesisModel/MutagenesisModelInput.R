##################################################
# description
##################################################

# this script creates 1) a table containing feature information (sequence context, chromatin etc.)
# for all insertion site in the user-supplied dataset (positives)
# 2) this table is combined with a feature table of a random subset of TTAA / TA sites (negatives).
# the combined table will be used to train the mutagenesis model in the next step

##################################################
# load packages
##################################################
suppressPackageStartupMessages({
library(GenomicFeatures)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(BSgenome.Mmusculus.UCSC.mm10)
library(stringr)
library(data.table)
library(regioneR)
library(tidymodels)
library(caret)
library(dplyr)
})
########################################
#read input
########################################

#source functions used to create the feature table
#source("/data/nasif12/home_if12/bredthau/workspace/transposon/Snakemake/Scripts/define_functions.R")
source(snakemake@params[["InputFunctions"]])

#assign genome
Tx <- makeTxDbFromGFF(snakemake@params[["REFSEQgtf"]], format = "gtf")
genomemouse <- BSgenome.Mmusculus.UCSC.mm10

#determines the size of the sequence context around the insertion size used for mutagenesis correction
window_size <- 10

#granges object with location of insertions
Granges_insertions <- readRDS(snakemake@input[["Granges_insertions"]])

#we have to use the same one-hot-encoder for the positives (insertion sites) that was previously used for the negatives (random TTAA / TA sites)
receipe <- readRDS(snakemake@input[["one_hot_encoder"]])

#read in a list of granges objects (locations of TSS, chromatin peaks etc.) to calculate feature vectors
featureList <- readRDS(snakemake@input[["featureList"]])


random_subset_negatives <- readRDS(snakemake@input[["step1_input_negatives_controls"]])
granges_subset_negatives <- readRDS(snakemake@input[["granges_negatives_controls"]])
Granges_insertions <- readRDS(snakemake@input[["Granges_insertions"]])


########################################
# clean the insertion list
########################################

#retrieve the sequence window needed for mutagenesis correction around the insertion site
window_ranges_positives <- retrieve_sequence_window(Granges_insertions, window_size)

#remove locations where sequence context cannot be properly retrieved
trimmed <- trim(window_ranges_positives)
removed <- which(window_ranges_positives%in%trimmed==FALSE)  
if(length(removed)!=0){
  window_ranges_positives <- window_ranges_positives[-removed]
  Granges_insertions <- Granges_insertions[-removed]
}

wrong_width<-which(width(window_ranges_positives)!=window_size*2+1)
if(length(wrong_width)!=0){
  window_ranges_positives<-window_ranges_positives[-wrong_width]
  Granges_insertions<-Granges_insertions[-wrong_width]
}

########################################
# build feature matrix (positives = insertion sites)
########################################


combined_input_chromosome <- create_features(one_bp_granges_object= Granges_insertions , window_ranges_positives, featureList)

#one hot encode the non-numeric variables  
step1_input_positives <- bake(receipe, new_data = combined_input_chromosome)

#remove NAs
contains_NA <-which(is.na(step1_input_positives), arr.ind=TRUE)[,1]
if(length(contains_NA)!=0){
  step1_input_positives <-step1_input_positives[-contains_NA,]
}

#for our model to be mathematically correct, the mutagenesis model has to be trained on unique insertion locations 
# -> multiple insertions at the same locus in different samples are filtered out
duplicates <- which(duplicated(Granges_insertions)==TRUE)
if(length(duplicates)!=0){
  step1_input_positives <- step1_input_positives[-duplicates,]
}

######################################
# remove insertion sites from the list of tandom TTAA / TA sites (negatives)
######################################

# we use a random subset of all TTAA / TA sites as negative controls in step 1
# here,  we exclude loci which correspond to insertion sites from the list of negatives


overlap_with_insertions <- which(granges_subset_negatives%in%Granges_insertions)
granges_subset_negatives <- granges_subset_negatives[-overlap_with_insertions]
random_subset_negatives <- random_subset_negatives[-overlap_with_insertions,]

######################################
# create combined data-frame used for training the mutagenesis model
######################################

#add vector containing the class for training
class_vector_positives <- rep(1, nrow(step1_input_positives)) # 1 = positives
step1_input_positives <- cbind(class_vector_positives,step1_input_positives)

class_vector_negatives <- rep(0, nrow(random_subset_negatives)) # 0 = negatives
random_subset_negatives <- cbind(class_vector_negatives, random_subset_negatives)

#randomly pick negatives for their number to equal the insertion sites
negatives_equal_number <- random_subset_negatives[sample(1:nrow(random_subset_negatives), nrow(step1_input_positives)),]
colnames(negatives_equal_number)[1] <- colnames(step1_input_positives)[1]

#combine positives and negatives into a single data frame
combined_data <- as.data.frame(rbind(step1_input_positives, negatives_equal_number))

#create training and testing data    
partition <- createDataPartition(as.factor(combined_data$class_vector_positives) , p=0.1)
combined_testing_data<- combined_data[partition$Resample1,]
combined_training_data <- combined_data[-partition$Resample1,]

########################################
# save output
########################################

saveRDS(combined_input_chromosome,snakemake@output[["step1_input_positives_RAW"]])
saveRDS(step1_input_positives,snakemake@output[["step1_input_positives"]])
saveRDS(combined_training_data,snakemake@output[["MutagenesisTrainingData"]])
saveRDS(combined_testing_data,snakemake@output[["MutagenesisTestData"]])





