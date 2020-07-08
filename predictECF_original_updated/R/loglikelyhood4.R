#taking the log-likely-hood function outside of the predict functions! I think this will be good for consistency and 
#will let me write the function better.
#version 3.2 tries to be fastes by using unique 
#version 4.0 can accept weights as input. 
#4.0 also removes option for loglikelyhood - we havent used it in like literaly years...

loglikelyhood <- function(actualresidue, protein, dna, row_weights = NA){
  
  #first make sure that the number of residues = the size of the proteinmatrix
  if (length(actualresidue)!=dim(as.matrix(protein))[2]){return(c(NA,NA,NA,NA))}
  
  #if there are no weights, all the weights are one
  if (is.na(row_weights[1])){row_weights <- rep(1, length(dna))}
  
  #this piece of code gets all of the things in protein which are the same as the thing of interest
  #it starts with everything in play and then narrows down only the things that match the sequence
  prots_in_play <- which(as.matrix(protein)[,1]==actualresidue[1])
  
  if (length(actualresidue) >1){
    for (i in 2:length(actualresidue)){
      
      prots_in_play <- intersect(prots_in_play, which(as.matrix(protein)[,i]==actualresidue[i]))
      
    }  
  }
  
  #if the combination doesnt exist, then what?
  if (sum(row_weights[prots_in_play]) <= 100 & length(actualresidue)>1){
    
    return(loglikelyhood_separate(actualresidue, protein, dna, row_weights))
    
  }
  
  #for log probability only we will optimize separately from the others
  #since this is all there is to this, and since this the the one I really use
  #adding acgt so theyll always be in the table output
  #ok, so we cant do the table here anymore, but we can still force acgt
  dna_of_intrest <- dna[prots_in_play]
  row_weights_of_int <- row_weights[prots_in_play]
  table_dna <- rep(NA,4)
  table_dna[1] <- sum(row_weights_of_int[which(dna_of_intrest == "a")])
  table_dna[2] <- sum(row_weights_of_int[which(dna_of_intrest == "c")])
  table_dna[3] <- sum(row_weights_of_int[which(dna_of_intrest == "g")])
  table_dna[4] <- sum(row_weights_of_int[which(dna_of_intrest == "t")])
  
  result <- log(table_dna/sum(table_dna),2)
  
  for (i in 1:4){
    if (is.na(result[i])){result[i] <- log(0.01 ,2)}  
    if (result[i] < -10000){result[i] <- log(0.01 ,2)}
  }
  
  return(result)
  
}

#############
#############
#this is the function that deals with separate ones - ie when pprot = 0 

loglikelyhood_separate <- function(actualresidue, protein, dna, row_weights = NA){
  
  temp_results <- NULL
  quantity <- c(0,0)
  
  for (i in 1:length(actualresidue)){
    
    res <- loglikelyhood(actualresidue[i], protein[,i], dna, row_weights)
    temp_results <- cbind(temp_results, res)
    
  }
  
  return(log(rowMeans(2^(temp_results)),2))

}
  