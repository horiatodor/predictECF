#
predict35w <- function(sequence, dnamatrix, proteinmatrix, row_weights=NA){

 motif <- matrix(0, 4,7)   
 rownames(motif) <- c('a','c','g','t')
 sequence <- tolower(sequence)
 
 #-36
 #this is a function of R176 
 motif[,1] <- loglikelyhood(sequence[43],proteinmatrix[,43], dnamatrix[,6], row_weights)
 
 #-35
 #this is a function of R176 
 motif[,2] <- loglikelyhood(sequence[46],proteinmatrix[,46], dnamatrix[,7], row_weights)

 #-34
 #this is a fucntion of S172 and R173 
 motif[,3] <- loglikelyhood(sequence[42:43],proteinmatrix[,42:43], dnamatrix[,8],row_weights)

 #-33
 motif[,4] <- loglikelyhood(sequence[c(38,42)], proteinmatrix[,c(38,42)], dnamatrix[,9],row_weights)

  #-32
 #this is a fucntion of s175 and and r176
 motif[,5] <- loglikelyhood(sequence[c(45,41)], proteinmatrix[,c(45,41)], dnamatrix[,10],row_weights)

 #-31
 #this is a fucntion of r171
 motif[,6] <- loglikelyhood(sequence[41],proteinmatrix[,41], dnamatrix[,11],row_weights)

 #-30
 #this is a fucntion of r149
 motif[,7] <- loglikelyhood(sequence[41],proteinmatrix[,41], dnamatrix[,12],row_weights)
 
 return(motif)
 
}

