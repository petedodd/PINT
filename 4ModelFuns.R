## TODO
## CDR BCG coverage functions etc
## remove most sanity checks so can be sourced
## think LHS

## read in data
PZ <- parse.parmtable(data = read.csv('data/PTBHHCT.csv')) #no test
## PZ <- parse.parmtable(data = read.csv('data/PTBHHCT.csv'), #test!
##                       outfile='test/zzout.csv',testdir = 'test')

names(PZ)


## =========== function definitions ===============

## what function are needed?
clx <- showParmz(kexp)$calcs
clx                                     #avoid the 1- or the 1
clx[!grepl('1',clx)]

## testing
##TODO define functions

## == co-prevalence (empirical)
coprev <- function(a){
  tmp <- PZ$coprev04$r(length(a))
  tmp[a>=5] <- PZ$coprev514$r(sum(a>=5))
  tmp
}
coprev(1:10)

## == case detection
CDR <- function(a){
  ## TODO from notification/estimates
  0.4
}

## == CFR on tx
CFRtxY <- function(a){#,hiv,art){        #NB optimized for clarity not speed
  tmp <- PZ$ontxY$r(length(a))
  tmp[a>=5] <- PZ$ontxO$r(sum(a>=5))  #NB this could be achieved in  the tree model
  ## hivartOR
  tmp
}
CFRtxY(1:10)                            #test

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

summary(progprob(runif(1e3,0,1),0.99,90))

## === disaggregating by TST status? may not be used??
RRtst <- function(a){
  PZ$RRtst10$r(length(a))
}

## === IPT efficacy
IPTrr <- function(a){
  PZ$iptRR$r(length(a))
}
IPTrr(1:10)





## ------ additionally ----
names(PZ)

## ## ================= tree testing =========


## tdf <- data.table(a=runif(1e4,0,15),iso3='LSO',lat=0,bcgcov=0.8)

## ## compute other data
## tdf[,CDR:=CDR(a)]
## tdf[,CFRtxY:=CFRtxY(a)]
## tdf[,CFRtxN:=CFRtxN(a)]
## tdf[,coprev:=coprev(a)]
## tdf[,IPTrr:=IPTrr(a)]
## tdf[,ltbi.prev:=ltbi.prev(a,coprev)]
## tdf[,pprogn:=progprob(a,bcgcov,lat)]
## tdf[,rrtst:=RRtst(a)]
## tdf[,inc:=ltbi.prev * pprogn]
## tdf[,progn.LP.PTn:=inc]                     #change for RR!
## tdf[,progn.LN.PTn:=0]                       #change for RR!
## tdf[,progn.LP.PTp:=progn.LP.PTn*IPTrr]
## tdf[,progn.LN.PTp:=progn.LN.PTn*IPTrr]
## tdf[,PTcov.N:=0]
## tdf[,PTcov.P:=0]
## tdf[,LE:=LE(a)]

## ## age categories
## tdf[,acs:=cut(a,breaks = 0:15,include.lowest = TRUE,right=FALSE,ordered_result = TRUE)]
## tdf[,ac:=cut(a,breaks=c(0,5,15),include.lowest = TRUE,right=FALSE,ordered_result = TRUE)]
## tdf[,acb:=cut(a,breaks=c(0,1,2,5,10,15),include.lowest = TRUE,right=FALSE,ordered_result = TRUE)]
## levels(tdf$acb)


## ## just looking at outcomes as simpler
## G <- makeTfuns(outcomes,unique(outcomes$fieldsAll))
## getAQ(outcomes,'check')

## summary(G$checkfun(tdf))
## print(outcomes,'check','p','LE')


## summary(G$de(tdf))
## summary(G$death(tdf))
## summary(G$deathfun(tdf))
## summary(G$LEfun(tdf))
## summary(G$incidencefun(tdf))


## summary(G$LEfun(tdf))
## getAQ(outcomes,'LE')

## ## print(outcomes,'death','treatments','incidence','LE','p')
## ## print(kexp,'p')
## ## print(kexp,'incidence','LE')
## ## print(kexp,'death','treatments')
## ## print(kexp,'prevalent','LTBI')



## F <- makeTfuns(kexp,unique(kexp$fieldsAll))
## getAQ(kexp,'LE')
## summary(F$checkfun(tdf))                         #!!




## summary(F$prevalentfun(tdf))
## summary(F$incidencefun(tdf))
## summary(F$deathfun(tdf))
## summary(F$LEfun(tdf))


## tdf$e.prevalent <- F$prevalentfun(tdf)
## tdf$e.incidence <- F$incidencefun(tdf)
## tdf$e.deaths <- F$deathfun(tdf)
## tdf$e.LE <- F$LEfun(tdf)



## ## print(kexp,'LTBI','p')


## ## print(kexp,'incidence','LE')
## ## print(kexp,'death','treatments')
## ## print(kexp,'prevalent','LTBI')

## plotter(kexp, varz=c('name','LE'), edgelabel = TRUE)
## plotter(kexp, varz=c('name','incidence'), edgelabel = TRUE)



## ## comparison by category
## tdf[,.(ep=mean(e.prevalent),
##        ei=mean(e.incidence),
##        ed=mean(e.deaths))]

## tdf[,.(ep=mean(e.prevalent),
##        ei=mean(e.incidence),
##        ed=mean(e.deaths)),
##     by=ac]

## ## incidence by TST status?


