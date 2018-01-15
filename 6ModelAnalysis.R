## TODO
rm(list=ls())
load('data/DL.Rdata')
source('3ModelDefn.R')
source('4ModelFuns.R')


DL[iso3=='ABW']
DL[iso3=='ZWE']

## simplifying for u5 only
DL[,LE:=(`[0,1)`+4*`[1,5)`)/5]
DL[,`[0,1)`:=NULL]; DL[,`[1,5)`:=NULL]; DL[,`[5,10)`:=NULL]; DL[,`[10,15)`:=NULL]; 

## country level version
DLC <- unique(DL[,.(iso3,g_whoregion,LAT,a1,a2,a3,a4,a5,cdr04,cdr04ab,LE)])

DLC                                     #error

## NB this is value hh stuff ( acat adult, but we can aggregate these)
DLK <- DL[,.(iso3,acat,sex,HHu5mu HHu5logsd)]      #will probably need to invlude HIV stuff here



## ages to harmonize or split for!
## LE
## bcgcov
## progression

## tasks:
## DF with numbers found?
## split numbers by age
## merge in PSA for parameters?
nrep <- 10                              #NB 100 gives 60Mb
DLL <- DL[rep(1:nrow(DL),nrep)]
DLL[,repn:=rep(1:nrep,each=nrow(DL))]

## TODO need to advance for on off interventions

## size means probably need to do as acat


## ================= tree testing =========


tdf <- data.table(a=runif(1e4,0,15),iso3='LSO',lat=0,bcgcov=0.8)

## compute other data
tdf[,CDR:=CDR(a)]
tdf[,CFRtxY:=CFRtxY(a)]
tdf[,CFRtxN:=CFRtxN(a)]
tdf[,coprev:=coprev(a)]
tdf[,IPTrr:=IPTrr(a)]
tdf[,ltbi.prev:=ltbi.prev(a,coprev)]
tdf[,pprogn:=progprob(a,bcgcov,lat)]
tdf[,rrtst:=RRtst(a)]
tdf[,inc:=ltbi.prev * pprogn]
tdf[,progn.LP.PTn:=inc]                     #change for RR!
tdf[,progn.LN.PTn:=0]                       #change for RR!
tdf[,progn.LP.PTp:=progn.LP.PTn*IPTrr]
tdf[,progn.LN.PTp:=progn.LN.PTn*IPTrr]
tdf[,PTcov.N:=0]
tdf[,PTcov.P:=0]
## tdf[,LE:=LE(a)]

## age categories
tdf[,acs:=cut(a,breaks = 0:15,include.lowest = TRUE,right=FALSE,ordered_result = TRUE)]
tdf[,ac:=cut(a,breaks=c(0,5,15),include.lowest = TRUE,right=FALSE,ordered_result = TRUE)]
tdf[,acb:=cut(a,breaks=c(0,1,2,5,10,15),include.lowest = TRUE,right=FALSE,ordered_result = TRUE)]
levels(tdf$acb)


## just looking at outcomes as simpler
G <- makeTfuns(outcomes,unique(outcomes$fieldsAll))
getAQ(outcomes,'check')

summary(G$checkfun(tdf))
print(outcomes,'check','p','LE')


summary(G$de(tdf))
summary(G$death(tdf))
summary(G$deathfun(tdf))
summary(G$LEfun(tdf))
summary(G$incidencefun(tdf))


summary(G$LEfun(tdf))
getAQ(outcomes,'LE')

## print(outcomes,'death','treatments','incidence','LE','p')
## print(kexp,'p')
## print(kexp,'incidence','LE')
## print(kexp,'death','treatments')
## print(kexp,'prevalent','LTBI')



F <- makeTfuns(kexp,unique(kexp$fieldsAll))
getAQ(kexp,'LE')
summary(F$checkfun(tdf))                         #!!




summary(F$prevalentfun(tdf))
summary(F$incidencefun(tdf))
summary(F$deathfun(tdf))
summary(F$LEfun(tdf))


tdf$e.prevalent <- F$prevalentfun(tdf)
tdf$e.incidence <- F$incidencefun(tdf)
tdf$e.deaths <- F$deathfun(tdf)
tdf$e.LE <- F$LEfun(tdf)



## print(kexp,'LTBI','p')


## print(kexp,'incidence','LE')
## print(kexp,'death','treatments')
## print(kexp,'prevalent','LTBI')

plotter(kexp, varz=c('name','LE'), edgelabel = TRUE)
plotter(kexp, varz=c('name','incidence'), edgelabel = TRUE)



## comparison by category
tdf[,.(ep=mean(e.prevalent),
       ei=mean(e.incidence),
       ed=mean(e.deaths))]

tdf[,.(ep=mean(e.prevalent),
       ei=mean(e.incidence),
       ed=mean(e.deaths)),
    by=ac]

## incidence by TST status?


