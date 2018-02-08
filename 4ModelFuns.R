## TODO
## remove most sanity checks so can be sourced
## think LHS
## source("3ModelDefn.R")

## read in data
PZ <- parse.parmtable(data = read.csv('data/PTBHHCT.csv')) #no test
## PZ <- parse.parmtable(data = read.csv('data/PTBHHCT.csv'), #test!
##                       outfile='test/zzout.csv',testdir = 'test')

names(PZ)


## =========== function definitions ===============

## ## what function are needed?
## clx <- showParmz(kexp)$calcs
## clx                                     #avoid the 1- or the 1
## clx[!grepl('1',clx)]

## testing

## == co-prevalence (empirical)
coprev <- function(a){
  tmp <- PZ$coprev04$r(length(a))
  tmp[a>=5] <- PZ$coprev514$r(sum(a>=5))
  tmp
}
coprev(1:10)

## == case detection
CDR <- function(mn,ab){
  ## mn <- mn/runif(length(mn),min=0.5,max=1) #CDR adjustment
  mn <- mn*(1/0.5 - (1/0.5-1)*runif(length(mn))) #CDR adjustment 2
  mn <- pmin(mn,1)
  a <- mn*ab
  b <- (1-mn)*ab
  rbeta(n=length(mn),shape1 = a,shape2 = b)
  ## 0.4
}
## TODO check and change

## == CFR on tx
CFRtxY <- function(a){#,hiv,art){        #NB optimized for clarity not speed
  tmp <- PZ$ontxY$r(length(a))
  tmp[a>=5] <- PZ$ontxO$r(sum(a>=5))  #NB this could be achieved in  the tree model
  ## hivartOR
  tmp
}
CFRtxY(1:10)                            #test
## TODO hivartOR

## == CFR off tx
CFRtxN <- function(a){#,hiv,art){
  tmp <- PZ$notxY$r(length(a))
  ## notxHY
  ## notxHAY
  tmp[a>=5] <- PZ$notxO$r(sum(a>=5))
  ## notxHO
  ## notxHAO
  tmp
}
CFRtxN(1:10)                            #test
## TODO hiv

## == LTBI infection probability
#NB this is LTBI given not active: it is taken to be max(0,LTBI-coprev)
ltbi.prev <- function(a,coprev){
  tmp <- PZ$LTBI04$r(length(a))
  tmp[a>=5] <- PZ$LTBI514$r(sum(a>=5))
  pmax(0,tmp - coprev)
}
ltbi.prev(1:10,0.1)


## ===  progression function (finer grained)
progprob <- function(a,bcgcov,lat){
  n <- length(a)
  pp <- pd <- rep(0,n)
  ## age categories
  A1 <- which(a<1); n1 <- length(A1)
  A2 <- which(a>=1 & a<2); n2 <- length(A2)
  A3 <- which(a>=2 & a<5); n3 <- length(A3)
  A4 <- which(a>=5 & a<10); n4 <- length(A4)
  A5 <- which(a>=10 & a<15); n5 <- length(A5)
  ## progn to TB by age
  pp[A1] <- PZ$pp1$r(n1)
  pp[A2] <- PZ$pp2$r(n2)
  pp[A3] <- PZ$pp3$r(n3)
  pp[A4] <- PZ$pp4$r(n4)
  pp[A5] <- PZ$pp5$r(n5)
  ## frac of TB DTB by age
  pd[A1] <- PZ$pd1$r(n1)
  pd[A2] <- PZ$pd2$r(n2)
  pd[A3] <- PZ$pd3$r(n3)
  pd[A4] <- PZ$pd4$r(n4)
  pd[A5] <- PZ$pd5$r(n5)
  ## BCG
  dBCG <- PZ$dBCG$r(n)                  #protn against DTB
  pBCG <- PZ$pBCG$r(n)                  #frac of above to PTB
  vBCG <- PZ$vBCG                       #frac protn lots to equator
  f <- (1-vBCG*(1-abs(lat)/90))         #frac protn at latitude
  ## combine
  ## tb <- pp*(1-(1-dBCG)*pBCG*f)
  ptb <- pp*(1-pd)*(1-bcgcov*(1-dBCG)*pBCG*f)  #pulmonary
  dtb <- pp*pd*(1-bcgcov*(1-dBCG)*f)           #non pulmonary
  ptb + dtb                             #all TB
}
## summary(progprob(runif(1e3,0,1),0.99,90))

## average for u5s
avu5progprob <- function(a1,a2,a3,a4,a5,lat){
  zs <- rep(0,length(a1))
  (progprob(zs+0.5,a1*1e-2,lat)+progprob(zs+1.5,a2*1e-2,lat)+
   progprob(zs+2.5,a3*1e-2,lat)+progprob(zs+3.5,a4*1e-2,lat)+
   progprob(zs+4.5,a5*1e-2,lat))/5
}

## average for o5s
avo5progprob <- function(a6,a7,a8,a9,a10,
                         a11,a12,a13,a14,a15,
                         lat){
  zs <- rep(5,length(a6))
  (progprob(zs+0.5,a6*1e-2,lat)+progprob(zs+1.5,a7*1e-2,lat)+
   progprob(zs+2.5,a8*1e-2,lat)+progprob(zs+3.5,a9*1e-2,lat)+
   progprob(zs+4.5,a10*1e-2,lat)+progprob(zs+5.5,a11*1e-2,lat)+
   progprob(zs+6.5,a12*1e-2,lat)+progprob(zs+7.5,a13*1e-2,lat)+
   progprob(zs+8.5,a14*1e-2,lat)+progprob(zs+9.5,a15*1e-2,lat)
   )/10
}


## === disaggregating by TST status? may not be used??
RRtst <- function(a){
  PZ$RRtst10$r(length(a))
}
## RRtst(rep(1,10))

## === IPT efficacy
IPTrr <- function(a){
  PZ$iptRR$r(length(a))
}
## IPTrr(1:10)





## ------ additionally ----
## names(PZ)
