
setwd("")
#the structure of this script is the same as 2alleles.R
#the parameters are adopted in the paper

runs <-10
drift <- 0

gen <- function(p,q,r,m,h,gens){
  for(i in 1:gens){
    p.prime <- p
    q.prime <- q
    r.prime <- r
    m.prime <- m
    h.prime <- h
    
    pr <- q+r+m+h 
    qr <- r+m+h
    rr <- m+h
    mr <- h
    
    random.matrix <- matrix(runif(initial.popsize,0,1), nrow= initial.popsize, ncol=1) 
    allele.m <-matrix(0, nrow=initial.popsize, ncol=1)
    for(i in 1:initial.popsize){
      if (random.matrix[i,] >= pr){allele.m[i,]<-"p"}
      else if(random.matrix[i,] >= qr & random.matrix[i,] < pr){allele.m[i,]<-"q"}
      else if(random.matrix[i,] >= rr & random.matrix[i,] < qr){allele.m[i,]<-"r"}
      else if(random.matrix[i,] >= mr & random.matrix[i,] < rr){allele.m[i,]<-"m"}
      else {allele.m[i,]<-"h"}
    }
    psum<-sum(allele.m=='p')
    qsum<-sum(allele.m=='q')
    rsum<-sum(allele.m=='r')
    msum<-sum(allele.m=='m')
    hsum<-sum(allele.m=='h')
    
    p <- psum/initial.popsize
    q <- qsum/initial.popsize
    r <- rsum/initial.popsize
    m <- msum/initial.popsize
    h <- hsum/initial.popsize
    
    YY <- p^2
    YI <- 2*p*q
    II <- q^2
    MY <- 2*p*m
    MM <- m^2
    MI <- 2*q*m
    
    YR <- 2*p*r
    IR <- 2*q*r
    MR <- 2*m*r
    RR <- r^2
    HY <- 2*p*h
    HI <- 2*q*h
    HR <- 2*h*r
    HM <- 2*h*m
    HH <- h^2
    
    wII <- aII*(1-z*II)
    
    wYY <- aYY*(1-z*YY-z*HY-z*YI)
    wHY <- aHY*(1-z*HY-z*YI-z*YY)
    wYI <- aYI*(1-z*YI-z*YY-z*HY)
    
    wHI <- aHI*(1-z*HI-z*HH)
    wHH <- aHH*(1-z*HH-z*HI)
    
    wMI <- aM*aMI*(1-z*MI)
    wMY <- aMY*(1-z*MY)
    
    wYR <- aYR*(1-z*YR)
    wIR <- aIR*(1-z*IR)
    wHR <- aHR*(1-z*HR)
    
    wHM <- aM*aHM*(1-z*HM)
    wMR <- aM*aMR*(1-z*MR)
    
    wRR <- 0
    
    wMM <- 0
    
    
    wmean <- wYY*YY+wYI*YI+wII*II+wMI*MI+wMY*MY+wHY*HY+wHI*HI+wHH*HH+wYR*YR+wIR*IR+wHR*HR+wMR*MR+wRR*RR+wHM*HM
    
    if(r*q*m*h*p == 0){
      drift = drift + 1 
      freq <- gen(p.prime,q.prime,r.prime,m.prime,h.prime,1)
      p <- freq[1]
      q <- freq[2]
      r <- freq[3]
      m <- freq[4]
      h <- freq[5]
      
    }
    else{
      p <- (wYY*YY+wYI*YI*0.5+wMY*MY*0.5+wYR*YR*0.5+wHY*HY*0.5)/wmean
      q <- (wII*II+wYI*YI*0.5+wMI*MI*0.5+wIR*IR*0.5+wHI*HI*0.5)/wmean
      r <- (wRR*RR+wIR*IR*0.5+wYR*YR*0.5+wHR*HR*0.5+wMR*MR*0.5)/wmean
      m <- (wMY*MY*0.5+wMI*MI*0.5+wMR*MR*0.5+wHM*HM*0.5)/wmean
      h <- (wHH*HH+wHI*HI*0.5+wHY*HY*0.5+wHR*HR*0.5+wHM*HM*0.5)/wmean
    }
  }

  gYY<-YY*wYY/wmean
  gII<-II*wII/wmean
  gHH<-HH*wHH/wmean
  
  gHY<-HY*wHY/wmean
  gHI<-HI*wHI/wmean
  gYI<-YI*wYI/wmean
  
  gMM<-MM*wMM/wmean
  gRR<-RR*wRR/wmean
  
  gMR<-MR*wMR/wmean
  gHR<-HR*wHR/wmean
  gIR<-IR*wIR/wmean
  gYR<-YR*wYR/wmean
  
  gHM<-HM*wHM/wmean
  gMI<-MI*wMI/wmean
  gMY<-MY*wMY/wmean
  
  return(c(p,q,r,m,h,drift,gII,gHH,gYY,gRR,gMM,gHI,gYI,gIR,gMI,gHY,gHR,gHM,gYR,gMY,gMR))
}


for(i in seq(1,1,1)){
  test <- matrix(0, nrow=21, ncol=runs)
  
  z <- 0.2
  
  aII <- 1
  aYY <- 1.04
  aHH <- 0.96
  
  aHY <- 1.045
  aYI <- 0.98
  aHI <- 1
  
  aMI <- 1
  aMY <- 1.1
  aHM <- 1
  
  aYR <- 1
  aIR <- 0.9
  aHR <- 1.04
  aMR <- 1 
  
  aM <- 0.9
  
  initial.pfreq<- 0.2695
  initial.qfreq <- 0.7016
  initial.rfreq<- 0.0238
  initial.mfreq<- 0.0049
  initial.hfreq<- 0.0001
  initial.popsize <- 10000
  
  generations <- 400		
  
  for(j in seq(1,runs,1)){
    test[,j] <- gen(initial.pfreq, initial.qfreq,initial.rfreq,initial.mfreq,initial.hfreq,generations)
  }

  test <- as.data.frame(test)
  rownames(test) <- c("s","p","r","m","h","drift","gPP","gVV","gSS","gRR","gMM","gPV","gSP","gPR","gPM","gSV","gVR","gVM","gSR","gSM","gMR")
  write.table(test,file="./5allele.csv",col.names = TRUE, append=TRUE, sep = ",")
}
par(mfrow=c(1,2),mar = c(4,4,1,1))
box<-t(test)
box<-subset(box,select = -c(drift))
box<-box[,6:20]


boxplot(box,ylim=c(0,0.35))
real<-read.csv("../real_for_plot.csv")
library(forcats)
real$genotype<-fct_inorder(real$genotype)
plot(real,type="dot",ylim=c(0,0.35))
