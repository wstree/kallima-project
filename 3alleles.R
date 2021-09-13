
setwd("")

#the structure of this script is the same as 2alleles.R
#the parameters are adopted in the paper

runs <-10
drift <- 0

gen <- function(p,q,r,gens){
  for(i in 1:gens){
    p.prime <- p
    q.prime <- q
    r.prime <- r
    
    pr <- q+r
    qr <- r
    
    random.matrix <- matrix(runif(initial.popsize,0,1), nrow= initial.popsize, ncol=1) 
    allele.m <-matrix(0, nrow=initial.popsize, ncol=1)
    for(i in 1:initial.popsize){
      if (random.matrix[i,] >= pr){allele.m[i,]<-"p"}
      else if(random.matrix[i,] >= qr & random.matrix[i,] < pr){allele.m[i,]<-"q"}
      else {allele.m[i,]<-"r"}
    }
    psum<-sum(allele.m=='p')
    qsum<-sum(allele.m=='q')
    rsum<-sum(allele.m=='r')
    
    p <- psum/initial.popsize
    q <- qsum/initial.popsize
    r <- rsum/initial.popsize
    
    SS <- p^2
    SP <- 2*p*q
    PP <- q^2
    
    SR <- 2*p*r
    PR <- 2*q*r
    RR <- r^2
    
    wPP <- aPP*(1-z*PP)
    
    wSS <- aSS*(1-z*SS-z*SP)
    wSP <- aSP*(1-z*SP-z*SS)
    
    wSR <- aSR*(1-z*SR)
    wPR <- aPR*(1-z*PR)
    
    wRR <- 0
    
    
    wmean <- wSS*SS+wSP*SP+wPP*PP+wSR*SR+wPR*PR+wRR*RR
    
    if(r*q*p == 0){
      drift = drift + 1 
      freq <- gen(p.prime,q.prime,r.prime,1)
      p <- freq[1]
      q <- freq[2]
      r <- freq[3]
      
    }
    else{
      p <- (wSS*SS+wSP*SP*0.5+wSR*SR*0.5)/wmean
      q <- (wPP*PP+wSP*SP*0.5+wPR*PR*0.5)/wmean
      r <- (wRR*RR+wPR*PR*0.5+wSR*SR*0.5)/wmean
      
    }
    
  }
  gSS <- p^2*wSS/wmean
  gPP <- q^2*wPP/wmean
  gSP <- 2*p*q*wSP/wmean
  gRR <- r^2*wRR/wmean
  gRP <- 2*r*q*wPR/wmean
  gSR <- 2*p*r*wSR/wmean  
  return(c(p,q,r,drift,gSS,gPP,gSP,gRR,gRP,gSR))
}


for(i in seq(0.6,1,0.02)){
  test <- matrix(0, nrow=10, ncol=runs)
  
  z <- 0.2
  aPR <- i
  aSP <- 1
  aPP <- 1
  
  aSS <- 1
  aSP <- 1
  
  aSR <- 1
  
  initial.pfreq<- 0.2799
  initial.qfreq <- 0.72
  initial.rfreq<- 0.0001
  initial.popsize <- 10000
  
  generations <- 400	
  
  for(j in seq(1,runs,1)){
    test[,j] <- gen(initial.pfreq, initial.qfreq,initial.rfreq,generations)
  }
  test <- as.data.frame(test)
  rownames(test) <- c("p","q","r","drift","gSS","gPP","gPS","gRR","gPR","gSR")
  write.table(test,file="",col.names = TRUE, append=TRUE, sep = ",")
}
