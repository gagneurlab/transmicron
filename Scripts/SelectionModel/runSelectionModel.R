##################################################
# description
##################################################

# this script runs the selection model (step2) and creates the output table

#####################################################
# load packages
#####################################################
suppressPackageStartupMessages({
library(data.table)
library(doParallel)
library(MASS)
library(dplyr)
})
#####################################################
#read input data
#####################################################

#SelectionModelInput50kb <- readRDS("/s/project/transposon/Output/PublishedPipeline/ModelInput/ModelInput10kb.RData")
#SelectionModelInputTargetAnnotation <- readRDS("/s/project/transposon/Output/PublishedPipeline/ModelInput/ModelInputGenes.RData")
SelectionModelInput50kb<-readRDS(snakemake@input[["SelectionModelInput50kb"]])
SelectionModelInputTargetAnnotation<-readRDS(snakemake@input[["SelectionModelInputTargetAnnotation"]])

multest_correction <- snakemake@params[["multest_correction"]]

#####################################################
# fit model per genomic region (50kb)
# -> equation 4 manuscript
#####################################################

# clean up the input tables: if there are 0 TA / TTAA sites in a gene / bin, the model does not fit
# pvalue = NA will be assigned to these genes / bins

remove_no_TTAA<- function(inputTable){
  remove_no_TTAA <- which(inputTable$sum_class_probabilities==0)
  if(length(remove_no_TTAA)!=0){
    inputTable<-as.data.frame(as.data.table(inputTable[-remove_no_TTAA,]))
  }
  return(inputTable)
}

SelectionModelInput50kbClean <- remove_no_TTAA(SelectionModelInput50kb)
SelectionModelInputTargetAnnotationClean <- remove_no_TTAA(SelectionModelInputTargetAnnotation)


####
#fit sample-coefficient (iteratively)
####

#for parallelization
# @ata: you can choose whatever parallelization works best
# maxcores <-detectCores()
# used_cores <- min(5,maxcores)
#cl <- parallel::makeCluster(used_cores,outfile="")
#doParallel::registerDoParallel(cl)

sample_coefs <-foreach(i=1:uniqueN(SelectionModelInput50kbClean$samples), .combine = "rbind")%do%{
  #subset for single sample
  sample_name <- unique(SelectionModelInput50kbClean$samples)[i]
  SelectionModelInput50kbClean$sample <- SelectionModelInput50kbClean$samples==sample_name
  #fit model
  model_poisson_sample<-glm.fit(SelectionModelInput50kbClean$sample,SelectionModelInput50kbClean$insertion_count_ij, family = poisson())
  #collect results
  sample_coefs_i <- c(sample_name,coef(model_poisson_sample))
  sample_coefs_i
}

#stopCluster(cl)

#format output
sample_coefs <- as.data.frame(sample_coefs)
colnames(sample_coefs) <- c("sample_names", "coefficient")
SelectionModelInput50kbClean$sample_coefficient <- as.numeric(sample_coefs$coefficient[match(SelectionModelInput50kbClean$samples, sample_coefs$sample_name)])
SelectionModelInput50kbClean$sample <- NULL #removes temportary sample logical

####
# fit complete model, incl. chromosome coefficient (see manuscript, equation 4)
####

model_global_complete <- glm(insertion_count_ij ~ 0+ offset(log(sum_class_probabilities))+ sample_coefficient + chromosome,
                             data =SelectionModelInput50kbClean,
                             maxit=10000,
                             family = "poisson"(link="log"))

#unused: fit the glm in one step; this is equivalent to the iteraitve process above, but to slow and memory intensive
# model_complete <- (glm(insertion_count_ij ~ offset(log(sum_class_probabilities))+samples + chromosome,
#                        data =SelectionModelInput50kbClean,maxit=10000, family = "poisson"(link="log")))


######################################
# gene-wise GLM-based one-sided Wald-test (manuscript equation 5)
######################################

#generate matrix with actual and fitted values
create_result_matrix<- function (sample_coefs,model_complete, input_data){ # function to create matrix with real and fitted values, x = model-object, y= input data frame used to fit the model
  
  #look-up sample coefficients for the target annotation
  input_data$sample_coefficient <- 
      SelectionModelInput50kbClean$sample_coefficient[
        match(input_data$samples,
              SelectionModelInput50kbClean$samples)]
  
  #predict expected insertions for every gene (or other target element defined by user):
  #above, we fitted equation 4 on 50kb bins; here, we have to create scores for genes / target elements
  
  fitted_complete <- predict(model_global_complete, input_data,type="response")
  
  result_matrix<-data.table(feature_name=input_data$feature_name,
                            actual=input_data$insertion_count_ij,
                            fitted=as.numeric(fitted_complete))
  result_matrix$actual<-as.numeric(as.character(result_matrix$actual))
  result_matrix$fitted<-as.numeric(as.character(result_matrix$fitted))
  return(result_matrix)
}

input_featurewise_theta<-create_result_matrix(sample_coefs,
                                              model_global_complete,
                                              SelectionModelInputTargetAnnotationClean)

#aggregate input table to remove genes with 0 insertions (model fails otherwise)
input_featurewise_thetaaggr <- aggregate(.~feature_name ,input_featurewise_theta, sum)
GeneWithoutInsertions <- unique(input_featurewise_thetaaggr$feature_name[which(input_featurewise_thetaaggr$actual == 0)])
modelInput <- input_featurewise_theta[-which(input_featurewise_theta$feature_name%in%GeneWithoutInsertions),]


#use for parallelization
# @ata: you can choose whatever parallelization works best
maxcores <-detectCores()
used_cores <- min(21,maxcores)
cl <- parallel::makeCluster(used_cores,outfile='log.txt')
doParallel::registerDoParallel(cl)

#fit actual vs. fitted for single genes
summary_frame <-foreach(i=1:uniqueN(modelInput$feature_name), .combine = "rbind")%dopar%{
  
  #subset for single genes / features
  data = modelInput[which(modelInput$feature_name==unique(modelInput$feature_name)[i]),]
  
  #fit model 
  model_temp <- glm(actual~offset(log(fitted)),
                    data=data,
                    maxit=10000,family="poisson")
  
  #calculate pvalues of intercept (cannot take it from summary directly, we need one-sided pvalues)
  p_value <- pnorm(summary.glm(model_temp)$coefficients[1,3], lower.tail = FALSE)
  intercept <- coefficients(summary.glm(model_temp))[1]
  
  c(unique(modelInput$feature_name)[i], p_value, intercept)
}

stopCluster(cl)


######################################################################
# generate and format output 
######################################################################

summary_frame <- as.data.frame(summary_frame)
rownames(summary_frame) <- 1:nrow(summary_frame)
names(summary_frame) <- c("feature_name", "unadjusted_p_value", "intercept")
summary_frame$unadjusted_p_value <- as.numeric(summary_frame$unadjusted_p_value)
summary_frame$intercept <- as.numeric(summary_frame$intercept)

#generate result frame
results <- aggregate(insertion_count_ij ~ feature_name, #aggregate non-filtered model input matrix (to retain all genes, including without TTAA / TA sites and 0 insertions)
                     SelectionModelInputTargetAnnotation,
                     FUN=sum)

results$expectedInsertions <- input_featurewise_thetaaggr$fitted[match(results$feature_name, input_featurewise_thetaaggr$feature_name)]

results <- left_join(results, summary_frame)


#adjust p-value for multiple testing
results$adjusted_pvalues_featurewise <- p.adjust(results$unadjusted_p_value, method=as.character(multest_correction))

#add expected insertions 


#add chromosomes
results$chromosome <- SelectionModelInputTargetAnnotation[match(results$feature_name,SelectionModelInputTargetAnnotation$feature_name),]$chromosome

#add list of samples in which the target genes / bins are hit
results$samples <- NULL
SelectionModelInputTargetAnnotation$samples <- as.character(SelectionModelInputTargetAnnotation$samples)
exclude0_insertions <- which(SelectionModelInputTargetAnnotation$insertion_count_ij==0)
SelectionModelInputTargetAnnotation_aggregated_nonzero <- aggregate(.~feature_name, SelectionModelInputTargetAnnotation[-exclude0_insertions,], FUN=list)
samples <- SelectionModelInputTargetAnnotation_aggregated_nonzero[,c(1,2)]
results <- left_join(results, samples, by="feature_name")
results$samples <- gsub('"', "",results$samples)
results$samples <- gsub("c(","",results$samples, fixed = T)

#get sample count for each gene
count_by_sample <- plyr::count(SelectionModelInputTargetAnnotation[-exclude0_insertions,]$feature_name)
colnames(count_by_sample) <-c("feature_name", "tumor_count")
results <- left_join(results, count_by_sample, by="feature_name")

#get pattern count (how many TTAA / TA sites per gene or bin)



###############################################
# sort and filter the output
###############################################

#by default, we use this rule to determine the threshold of insertions
# -> all genes below the threshold should be assigned NA p-value
N_tumors <- uniqueN(SelectionModelInputTargetAnnotation$samples)

if(floor(0.05*uniqueN(N_tumors))>3){
  insertion_cutoff <-floor(0.05*uniqueN(N_tumors))
} else{
  insertion_cutoff <-3
}

# @ata: TBD should we let users decide what to filter for? if so, user input should override default
#insertion_cutoff <- as.numeric(snakemake@params[["insertion_cutoff"]]) 

#remove genes below a pre-specified insertion cutoff
remove_insertion_cutoff <- results[which((results$actual<insertion_cutoff)==TRUE),]$feature_name

#set p-values to NA
remove_names <- unique(remove_insertion_cutoff)
if(length(remove_names)!=0){
  results[which(results$feature_name%in%remove_names),]$unadjusted_p_value <- NA
  results[which(results$feature_name%in%remove_names),]$adjusted_pvalues_featurewise <- NA
}

# sort results by p-value
results <- results[order(results$adjusted_pvalues_featurewise),]

fwrite(results, snakemake@output[["results"]])
results <- fread("/s/project/transposon/Output/PublishedPipeline/DLBCLPB/PB/results/results_genes_NewInsRates_precompFeatures.csv")
# # 
#  catalogue <- readRDS("/s/project/transposon/Output/Snakemake/benchmarks/common/COSMICcomplete.RData")
# 
#  test <- readRDS("/s/project/transposon/Output/PublishedPipeline/DLBCL_PB_PB/MutagenesisTestData.RData")
#  library(stringr)
#   results$OncoVar <- str_to_title(results$feature_name)%in%str_to_title(catalogue)
#   results<- results[-which(results$adjusted_pvalues_featurewise>0.45),]
#   results<- results[-which(is.na(results$adjusted_pvalues_featurewise)),]
#   sum(results$OncoVar)/length(results$OncoVar)
# resultsControls <- readRDS("/s/project/transposon/Output/Snakemake/DLBCL_PB/sequence_DNAsemESCSRX1452763_AtacseqmESC_distancetoGenerelative_distancettsrelative_distancetssrelative/unique/genes/keep/poisson/predict/elementwise/sorted_result_poisson_withlog_aggregate_bonferroni.RData")
# resultsControls$rank <- rank(resultsControls$adjusted_pvalues_featurewise)
# resultsControls$ranknew <- results$rank[match(resultsControls$feature_name, results$feature_name)]
# resultsControls$rank50 <- results50kb$rank[match(resultsControls$feature_name, results50kb$feature_name)]

# 
# library(ggplot2)
# ggplot(data=resultsControls, aes(x=rank, y=ranknew))+
#   geom_point()+
#   scale_x_log10()+
#   scale_y_log10()
# #############################
# #done
# #############################
# 
# 
# # inputQQ <-results[-which(str_to_title(results$feature_name)%in%str_to_title(catalogue)),]
# # library(qqplotr)
# # df_mutag <-data.table(inputQQ$unadjusted_p_value,group ="Our Model Mutagenesis")
# # qqplot <- ggplot(df_mutag, aes(sample=as.numeric(V1)))+
# #   geom_qq(distribution = stats::qunif, aes(color=group)) + 
# #   #geom_qq_band(distribution = "unif", aes(color=group)) + 
# #   geom_abline()
# #   #ggplot_theme
# # qqplot




