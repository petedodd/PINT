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
DLC

## NB this is value hh stuff ( acat adult, but we can aggregate these)
DLK <- DL[,.(iso3,acat,sex,value,HHu5mu,HHu5logsd)]      #will probably need to invlude HIV stuff here


nrep <- 1e3
DLKL <- DLK[rep(1:nrow(DLK),nrep),]
DLKL[,repn:=rep(1:nrep,each=nrow(DLK))]

## u5 contacts found
DLKL[,phh := rlnorm(n = nrow(DLKL), meanlog=HHu5mu, sdlog=HHu5logsd)]
## print(DLKL[phh>1e2,],n=Inf)                            #TODO problem for rich countries!
DLKL[phh>5,phh:=runif(sum(phh>5),max=5)]               #TODO temporary!


DLKL[,u5hhc := value * phh]

DLKL[,summary(phh)]
DLKL[,qplot(phh)]           #

## aggregates
chhc <- DLKL[,.(u5hhc=sum(u5hhc),notes=sum(value),phh=mean(phh)),
             by=.(repn,iso3)] #country aggregates

ghhc <- DLKL[,.(u5hhc=sum(u5hhc),notes=sum(value),phh=mean(phh)),
             by=.(repn)] #global aggregates

summary(ghhc)
ghhc[,mean(u5hhc)*1e-6]                 #around 3
chhc[iso3=='ZAF',summary(u5hhc)]/DLK[iso3=='ZAF',sum(value)]

chhc[iso3=='ZAF',qplot(u5hhc)]
chhc[iso3=='AFG',qplot(u5hhc)]

chhc <- chhc[,.(u5hhc=mean(u5hhc),u5hhc.sd=sd(u5hhc),
                u5hhc.l=mean(log(u5hhc)),u5hhc.sdl=sd(log(u5hhc))),by=iso3]

save(chhc,file='data/chhc.Rdata')

## merge into parent data table
DLC <- merge(DLC,chhc,by='iso3')


## LE
## bcgcov
## progression

## TODO need to advance for on off interventions


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
## ## tdf[,LE:=LE(a)]

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


