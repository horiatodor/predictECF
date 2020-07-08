#
#
predict10w <- function(sequence, dnamatrix, proteinmatrix, row_weights = NA){
  
  motif <- matrix(0, 4,7)   
  rownames(motif) <- c('a','c','g','t')
  sequence <- tolower(sequence)
  
  #-13
  #this is a function of d84
  #if (sequence[53]=='d'){
  motif[,1] <- loglikelyhood(sequence[53],proteinmatrix[,53], dnamatrix[,9],row_weights)
#  motif[,1] <- loglikelyhood(c(sequence[53],sequence[56]),
#                             cbind(proteinmatrix[,53],proteinmatrix[,56]), dnamatrix[,9],row_weights)
  
#  if (sequence[53] == "n"){
#    of_int <- which(proteinmatrix[,53] == "n")  
#    motif[,1] <- loglikelyhood(sequence[52],proteinmatrix[of_int,52], dnamatrix[of_int,9],typeofmetric)
#  }
  
  #-12
  #this is a fucntion of d84
  motif[,2] <- loglikelyhood(sequence[53],proteinmatrix[,53], dnamatrix[,10],row_weights)
  
#  if (sequence[53] == "n"){
#    of_int <- which(proteinmatrix[,53] == "n")  
#    motif[,2] <- loglikelyhood(sequence[52],proteinmatrix[of_int,52], dnamatrix[of_int,10],typeofmetric)
#  }
  
  #-11
  #this is a fucntion of k56 and asn80 and maybe arg77
  motif[,3] <- loglikelyhood(c(sequence[12],sequence[48]),
                             cbind(proteinmatrix[,12],proteinmatrix[,48]), dnamatrix[,11],row_weights)

#  if (sequence[12] %in% c("r","k")){
#    motif[,3] <- loglikelyhood(sequence[12],proteinmatrix[,12], dnamatrix[,11],row_weights)
#  }
  
  #-10
  #this is a fucntion of r65, d67,s68
  motif[,4] <- loglikelyhood(c(sequence[22],sequence[26]), 
                             cbind(proteinmatrix[,22],proteinmatrix[,26]),
                             dnamatrix[,12],row_weights)
  
  #-9
  #this is a fucntion of 
  motif[,5] <- loglikelyhood(sequence[40],proteinmatrix[,40], dnamatrix[,13],row_weights)

  #-8
  #this is a fucntion of ser68 (27)
  motif[,6] <- loglikelyhood(sequence[45],proteinmatrix[,45], dnamatrix[,14],row_weights)

  #-7
  #this is a fucntion of tyr75
  motif[,7] <- loglikelyhood(sequence[41],proteinmatrix[,41], dnamatrix[,15],row_weights)
  
  return(motif)
 
}
