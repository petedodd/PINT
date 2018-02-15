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
oddit <- function(x) x/(1-x)
ioddit <- function(x) x/(1+x)
logit <- function(x) log(oddit(x))
ilogit <- function(x) ioddit(exp(x))

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
CFRtxY <- function(a,hiv=0,art=0){#NB optimized for clarity not speed
  if(length(a)>1 & length(hiv)==1) hiv <- rep(hiv,length(a))
  if(length(a)>1 & length(art)==1) art <- rep(art,length(a))
  tmp <- PZ$ontxY$r(length(a))
  tmp[a>=5] <- PZ$ontxO$r(sum(a>=5))  #NB this could be achieved in  the tree model
  ## hivartOR
  Z <- PZ$hivartOR$r(length(a))
  hor <- rep(1,length(a))
  tmp <- logit(tmp)                     #transformt
  tmp[hiv>0] <- tmp[hiv>0]+Z[hiv>0,1]
  tmp[art>0] <- tmp[art>0]+Z[art>0,2]
  tmp <- ilogit(tmp)                    #inverse transform
  tmp
}
CFRtxY(1:10)                            #test
summary(CFRtxY(1:1e3))
summary(CFRtxY(1:1e3,hiv=1))
summary(CFRtxY(1:1e3,hiv=1,art=1))


## == CFR off tx
CFRtxN <- function(a,hiv=0,art=0){
  if(length(a)>1 & length(hiv)==1) hiv <- rep(hiv,length(a))
  if(length(a)>1 & length(art)==1) art <- rep(art,length(a))
  tmp <- PZ$notxY$r(length(a))          #default a<5 and hiv=art=0
  tmp[a<5 & hiv>0 & art==0] <- PZ$notxHY$r(sum(a<5 & hiv>0 & art==0)) #u5,HIV+,ART-
  tmp[a<5 & hiv>0 & art>0] <- PZ$notxHAY$r(sum(a<5 & hiv>0 & art>0)) #u5,HIV+,ART+
  tmp[a>=5] <- PZ$notxO$r(sum(a>=5))    #o5, HIV-ve
  tmp[a>=5 & hiv>0 & art==0] <- PZ$notxHO$r(sum(a>=5 & hiv>0 & art==0)) #o5,HIV+,ART-
  tmp[a>=5 & hiv>0 & art>0] <- PZ$notxHAO$r(sum(a>=5 & hiv>0 & art>0)) #o5,HIV+,ART+
  tmp
}
CFRtxN(1:10)                            #test
summary(CFRtxN(1:1e3))
summary(CFRtxN(1:1e3,hiv=1))
summary(CFRtxN(1:1e3,hiv=1,art=1))

## == LTBI infection probability
#NB this is LTBI given not active: it is taken to be max(0,LTBI-coprev)
ltbi.prev <- function(a,coprev){
  tmp <- PZ$LTBI04$r(length(a))
  tmp[a>=5] <- PZ$LTBI514$r(sum(a>=5))
  tmp
  ## pmax(0,tmp - coprev) # already taken into account with decision tree
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
avu5progprob <- function(a1,a2,a3,a4,a5,lat,hiv=0,art=0){
  zs <- rep(0,length(a1))
  if(length(a1)>1 & length(hiv)==1) hiv <- rep(hiv,length(a1))
  if(length(a1)>1 & length(art)==1) art <- rep(art,length(a1))
  ans <- (progprob(zs+0.5,a1*1e-2,lat)+progprob(zs+1.5,a2*1e-2,lat)+
          progprob(zs+2.5,a3*1e-2,lat)+progprob(zs+3.5,a4*1e-2,lat)+
          progprob(zs+4.5,a5*1e-2,lat))/5
  if(any(hiv>0)){ #treat as IRR for escape
    hr <- PZ$hivpi$r(sum(hiv>0))
    ans[hiv>0] <- 1-(1-ans[hiv>0])^hr
  }
  if(any(art>0)){ #treat as IRR for escape
    hr <- PZ$artp$r(sum(art>0))
    ans[art>0] <- 1-(1-ans[art>0])^hr
  }
  ans
}
ate <- rep(.9,1e3)                      #test
summary(avu5progprob(ate,ate,ate,ate,ate,lat=ate))
summary(avu5progprob(ate,ate,ate,ate,ate,lat=ate,hiv=1))
summary(avu5progprob(ate,ate,ate,ate,ate,lat=ate,hiv=1,art=1))

## average for o5s
avo5progprob <- function(a6,a7,a8,a9,a10,
                         a11,a12,a13,a14,a15,
                         lat,hiv=0,art=0){
  zs <- rep(5,length(a6))
  if(length(a6)>1 & length(hiv)==1) hiv <- rep(hiv,length(a6))
  if(length(a6)>1 & length(art)==1) art <- rep(art,length(a6))
  ans <- (progprob(zs+0.5,a6*1e-2,lat)+progprob(zs+1.5,a7*1e-2,lat)+
          progprob(zs+2.5,a8*1e-2,lat)+progprob(zs+3.5,a9*1e-2,lat)+
          progprob(zs+4.5,a10*1e-2,lat)+progprob(zs+5.5,a11*1e-2,lat)+
          progprob(zs+6.5,a12*1e-2,lat)+progprob(zs+7.5,a13*1e-2,lat)+
          progprob(zs+8.5,a14*1e-2,lat)+progprob(zs+9.5,a15*1e-2,lat)
  )/10
  if(any(hiv>0)){ #treat as IRR for escape
    hr <- PZ$hivpi$r(sum(hiv>0))
    ans[hiv>0] <- 1-(1-ans[hiv>0])^hr
  }
  if(any(art>0)){ #treat as IRR for escape
    hr <- PZ$artp$r(sum(art>0))
    ans[art>0] <- 1-(1-ans[art>0])^hr
  }
  ans
}
summary(avo5progprob(ate,ate,ate,ate,ate,ate,ate,ate,ate,ate,lat=ate))
summary(avo5progprob(ate,ate,ate,ate,ate,ate,ate,ate,ate,ate,lat=ate,hiv=1))
summary(avo5progprob(ate,ate,ate,ate,ate,ate,ate,ate,ate,ate,lat=ate,hiv=1,art=1))



## === disaggregating by TST status - used for TST-driven interventions
RRtst <- function(a){
  PZ$RRtst10$r(length(a))
}
## RRtst(rep(1,10))

## === IPT efficacy
IPTrr <- function(a){
  PZ$iptRR$r(length(a))
}
## IPTrr(1:10)

