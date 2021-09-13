
setwd("")


# Define a function for generations
# assume two allels S/P, p/q denotes the frequency of the allele, respectively
# a,b,c denote the additional benefit or deleterious effect of the genotype
# z denotes the strength of frequency-dependent selection
# gSS/gPP/gSP denotes the genotype frequency of SS/PP/SP, respectively


gen <- function(p,q,gens){
  for(i in 1:gens){
    p.prime <- p
    q.prime <- q
    
    pr <- q
    
    random.matrix <- matrix(runif(initial.popsize,0,1), nrow= initial.popsize, ncol=1) 
    allele.m <-matrix(0, nrow=initial.popsize, ncol=1)
    for(i in 1:initial.popsize){
      if (random.matrix[i,] >= pr){allele.m[i,]<-"p"}
      else {allele.m[i,]<-"q"}
    }
    psum<-sum(allele.m=='p')
    qsum<-sum(allele.m=='q')
    
    p <- psum/initial.popsize
    q <- qsum/initial.popsize
    
    SS <- p^2
    SP <- 2*p*q
    PP <- q^2
    
    
    wPP <-a*(1-z*PP)
    
    wSS <-b*(1-z*SS-z*SP)
    wSP <-c*(1-z*SP-z*SS)
    
    
    wmean <- wSS*SS+wSP*SP+wPP*PP
    
    if(q*p == 0){
      drift = drift + 1 
      freq <- gen(p.prime,q.prime,1)
      p <- freq[1]
      q <- freq[2]
      
    }
    else{
      p <- (wSS*SS+wSP*SP*0.5)/wmean
      q <- (wPP*PP+wSP*SP*0.5)/wmean
      gSS <- p^2*wSS/wmean
      gPP <- q^2*wPP/wmean
      gSP <- 2*p*q*wSP/wmean
    }
  }
  
  return(c(p,q,drift,gSS,gPP,gSP))
}

#starting drift is 0

for(i in seq(1,1,1)){
    runs <-
    drift <- 0
    test <- matrix(0, nrow=6, ncol=runs)
  
    z <- 
    a <- 
    b <- 
    c <-
    
    
    initial.pfreq<- 
    initial.qfreq <- 1-initial.pfreq
  initial.popsize <- 
    
    generations <- 		#Number of gens to cycle through for allele frequencies
    
    for(j in seq(1,runs,1)){
      test[,j] <- gen(initial.pfreq, initial.qfreq,generations)
    }
  test <- as.data.frame(test)
  rownames(test) <- c("p","q","drift","gSS","gPP","gSP")
  write.table(test,file="",col.names = TRUE, append=TRUE, sep = ",")
}
