##################################################
# description
##################################################

#this script defines functions that will be used to create the feature input tables 
# (for training the mutagenesis model)


##################################################
# sequence window management
##################################################

#this function is used to extend one bp insertion sites by a user determined window (for mutagenesis correction of sequence context)
retrieve_sequence_window <-function(one_bp_object, window_size){
  start(one_bp_object) <- start(one_bp_object) - window_size
  end(one_bp_object) <- end(one_bp_object) + window_size
  return(one_bp_object)
}

#this function is used to convert a granges object containing TTAA / TA sites to a granges object with single bp insertion loci
get_insertion_site_fromPattern <-function(granges_pattern, matching_pattern){
  
  
  strand<-as.vector(strand(granges_pattern))
  index_minus<-which(strand=="-")
  index_plus<-which(strand=="+")
  
  if(matching_pattern == "TTAA"){
    
    end(granges_pattern[index_plus]) <- end(granges_pattern[index_plus]) -2
    start(granges_pattern[index_minus]) <- start(granges_pattern[index_minus]) +2
    
    start(granges_pattern[index_plus]) <- start(granges_pattern[index_plus]) +1
    end(granges_pattern[index_minus]) <- end(granges_pattern[index_minus]) -1
  }
  
  else if(matching_pattern == "TA"){
    end(granges_pattern[index_plus]) <- end(granges_pattern[index_plus]) -1
    start(granges_pattern[index_minus]) <- start(granges_pattern[index_minus]) +1
  }
  
  return(granges_pattern)
}

##################################################
# create feature table (can be used ireatively over chromosomes)
##################################################

#this function is used to create the input tables containing the features for every insertions
create_features <- function(one_bp_granges_object, window_ranges_object, list_inputs){
  
  
  ######################################
  #sequence
  ######################################
  
  sequence <- as.matrix(DNAStringSet(getSeq(genomemouse,window_ranges_object))) #retrieve the sequence
  colnames(sequence) <- paste(rep("sequence", ncol(sequence)), 1:ncol(sequence), sep="_")
  
  ######################################
  #other features
  ######################################
  
  #this function is used to calculate the strand-specific distance of insertion sites / TTAA/TA sites to the features 
  calculateDistanceWithSign <- function(feature, one_bp_granges_object){
    distanceToNearestNeighbor <- as.data.frame(distanceToNearest(one_bp_granges_object, feature, ignore.strand=TRUE)) #find nearest neighbor
    downstream <- which(start(one_bp_granges_object[distanceToNearestNeighbor$queryHits])<start(feature[distanceToNearestNeighbor$subjectHits])) #check if insertion site is to the left of nearest neighbour
    distanceToNearestNeighbor$distance[downstream] <- distanceToNearestNeighbor$distance[downstream]*-1 #if so: make distance negative
    
    if(length(which(strand(feature[distanceToNearestNeighbor$subjectHits])=="-"))){ #if nearest neighbor is assigned to minus strand, reverse sign
      distanceToNearestNeighbor$distance[which(strand(feature[distanceToNearestNeighbor$subjectHits])=="-")]<- distanceToNearestNeighbor$distance[which(strand(feature[distanceToNearestNeighbor$subjectHits])=="-")]*-1
    }
    return(distanceToNearestNeighbor$distance)
  }
  
  #calculate distance between features and TTAA / TA / insertion sites
  other_features <- lapply(list_inputs, calculateDistanceWithSign, one_bp_granges_object=one_bp_granges_object )
  #convert to data frame
  other_features <- as.data.frame(bind_cols(other_features))

  #combine sequence matric and other features
  completeInputData <- cbind(sequence, other_features)
  
  return(completeInputData)
}

