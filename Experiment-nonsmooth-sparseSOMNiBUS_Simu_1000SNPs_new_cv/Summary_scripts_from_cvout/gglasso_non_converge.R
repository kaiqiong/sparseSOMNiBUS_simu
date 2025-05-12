
M = 50
bestLams <- matrix(NA, nrow = M, ncol = 2)

for( i in 1:M){
  
    try({
    load(paste0("Res", i, ".RData"))
    

    
    
    if(length(cvout$bestLambda1)>2){
      print(i)
    }
   # bestLams[i,] <- c(cvout$bestLambda1, cvout$bestLambda2)
    
    })  
    
}
