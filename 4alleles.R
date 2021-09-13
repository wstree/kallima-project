
setwd("")

#the structure of this script is the same as 2alleles.R
#the parameters are adopted in the paper

runs <-10 
drift <- 0

gen <- function(p,q,r,m,gens){
  for(i in 1:gens){
    p.prime <- p
    q.prime <- q
    r.prime <- r
    m.prime <- m
    
    
    pr <- q+r+m
    qr <- r+m
    rr <- m
    
    
    random.matrix <- matrix(runif(initial.popsize,0,1), nrow= initial.popsize, ncol=1) 
    allele.m <-matrix(0, nrow=initial.popsize, ncol=1)
    for(i in 1:initial.popsize){
      if (random.matrix[i,] >= pr){allele.m[i,]<-"p"}
      else if(random.matrix[i,] >= qr & random.matrix[i,] < pr){allele.m[i,]<-"q"}
      else if(random.matrix[i,] >= rr & random.matrix[i,] < qr){allele.m[i,]<-"r"}
      else {allele.m[i,]<-"m"}
    }
    psum<-sum(allele.m=='p')
    qsum<-sum(allele.m=='q')
    rsum<-sum(allele.m=='r')
    msum<-sum(allele.m=='m')
    
    p <- psum/initial.popsize
    q <- qsum/initial.popsize
    r <- rsum/initial.popsize
    m <- msum/initial.popsize
    
    SS <- p^2
    SP <- 2*p*q
    PP <- q^2
    MS <- 2*p*m
    MM <- m^2
    MP <- 2*q*m
    
    SR <- 2*p*r
    PR <- 2*q*r
    MR <- 2*m*r
    RR <- r^2
    
    wPP <- aPP*(1-z*PP)
    
    wSS <- aSS*(1-z*SS-z*SP)
    wSP <- aSP*(1-z*SP-z*SS)
    
    wMP <- aMP*(1-z*MP)
    wMS <- aMS*(1-z*MS)
    
    wSR <- aSR*(1-z*SR)
    wPR <- aPR*(1-z*PR)
    
    wMR <- aMR*(1-z*MR)
    
    wRR <- 0
    wMM <- 0
    
    wmean <- wSS*SS+wSP*SP+wPP*PP+wMP*MP+wMS*MS+wSR*SR+wPR*PR+wMR*MR+wRR*RR
    
    if(r*q*m*p == 0){
      drift = drift + 1 
      freq <- gen(p.prime,q.prime,r.prime,m.prime,1)
      p <- freq[1]
      q <- freq[2]
      r <- freq[3]
      m <- freq[4]
      
    }
    else{
      p <- (wSS*SS+wSP*SP*0.5+wMS*MS*0.5+wSR*SR*0.5)/wmean
      q <- (wPP*PP+wSP*SP*0.5+wMP*MP*0.5+wPR*PR*0.5)/wmean
      r <- (wRR*RR+wPR*PR*0.5+wSR*SR*0.5+wMR*MR*0.5)/wmean
      m <- (wMS*MS*0.5+wMP*MP*0.5+wMR*MR*0.5)/wmean
    }
  }
  gSS <- p^2*wSS/wmean
  gPP <- q^2*wPP/wmean
  gSP <- 2*p*q*wSP/wmean
  
  gRR <- r^2*wRR/wmean
  gPR <- 2*r*q*wPR/wmean
  gSR <- 2*p*r*wSR/wmean 
  
  gMP <- MP*wMP/wmean
  gMS <- MS*wMS/wmean
  gMR <- MR*wMR/wmean
  gMM <- MM*wMM/wmean
  
  return(c(p,q,r,m,drift,gSS,gPP,gRR,gMM,gSP,gPR,gSR,gMP,gMS,gMR))
}

#result <- matrix(0, nrow = 2, ncol = 1)
# Phase 1 S and P
# Fitness of homo equals fitness of hetero
for(i in seq(0.8,1,0.02)){
  test <- matrix(0, nrow=15, ncol=runs)
  
  z <- 0.2
  aPP <- 1
  
  aSS <- 1
  aSP <- 0.98
  
  aMP <- 1*i
  aMS <- 1.1
  
  aSR <- 1
  aPR <- 0.9
  aRR <- 1
  
  aMR <- 1*i
  
  
  initial.pfreq<- 0.29
  initial.qfreq <- 0.6899
  initial.rfreq<- 0.02
  initial.mfreq<- 0.0001
  
  initial.popsize <- 10000
  
  generations <- 400		
  
  for(j in seq(1,runs,1)){
    test[,j] <- gen(initial.pfreq, initial.qfreq,initial.rfreq,initial.mfreq,generations)
  }
  test <- as.data.frame(test)
  rownames(test) <- c("p","q","r","m","drift","gSS","gPP","gRR","gMM","gPS","gPR","gSR","gPM","gSM","gRM")
  write.table(test,file="./4allele.csv",col.names = TRUE, append=TRUE, sep = ",")
}
