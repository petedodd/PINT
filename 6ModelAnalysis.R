## TODO
## HIV functions
## CDR upped to reflect HHCT households
## TODO courses treats screens
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

## merge to make parent data table for PSA
DLC <- merge(DLC,chhc,by='iso3')

## build PSA dt
nrep <- 1e2
PSA <- DLC[rep(1:nrow(DLC),nrep),]
PSA[,repn:=rep(1:nrep,each=nrow(DLC))]


## TODO need to advance for on off interventions

## compute variables for PSA data table
azu5 <- rep(1,nrow(PSA))                #dummy ages for functions defined like that
PSA[,CDR:=CDR(cdr04,cdr04ab)]
PSA[,CFRtxY:=CFRtxY(azu5)]
PSA[,CFRtxN:=CFRtxN(azu5)]
PSA[,coprev:=coprev(azu5)]
PSA[,IPTrr:=IPTrr(azu5)]
PSA[,ltbi.prev:=ltbi.prev(azu5,coprev)]
PSA[,pprogn:=avu5progprob(a1,a2,a3,a4,a5,LAT)]
PSA[,rrtst:=RRtst(azu5)]
PSA[,inc:=ltbi.prev * pprogn]
PSA[,progn.LP.PTn:=inc]                     
PSA[,progn.LN.PTn:=0]                       
PSA[,progn.LP.PTp:=progn.LP.PTn*IPTrr]
PSA[,progn.LN.PTp:=progn.LN.PTn*IPTrr]
PSA[,PTcov.N:=0]
PSA[,PTcov.P:=0]

## intervention set
## TODO NB CDR also to change!
npsa <- nrow(PSA)
PSA <- PSA[rep(1:npsa,2),]
PSA[,intervention:=c(rep('basecase',npsa),rep('full',npsa))]
PSA[intervention=='full',PTcov.N:=1]
PSA[intervention=='full',PTcov.P:=1]
PSA[intervention=='full',CDR:=1]




## ------ tree model calculations


F <- makeTfuns(kexp,unique(kexp$fieldsAll))
getAQ(kexp,'LE')
summary(F$checkfun(PSA))                         #!!

## summary(F$prevalentfun(PSA))
## summary(F$incidencefun(PSA))
## summary(F$deathfun(PSA))
## summary(F$LEfun(PSA))


PSA$e.prevalent <- F$prevalentfun(PSA)*PSA$u5hhc
PSA$e.incidence <- F$incidencefun(PSA)*PSA$u5hhc
PSA$e.deaths <- F$deathfun(PSA)*PSA$u5hhc
PSA$e.LE <- F$LEfun(PSA)*PSA$u5hhc

## print(kexp,'incidence','LE')
## print(kexp,'death','treatments')
## print(kexp,'prevalent','LTBI')

plotter(kexp, varz=c('name','LE'), edgelabel = TRUE)
plotter(kexp, varz=c('name','incidence'), edgelabel = TRUE)



## comparison by category
PSAR <- PSA[,.(ep=sum(e.prevalent),
               ei=sum(e.incidence),
               ed=sum(e.deaths),
               el=sum(e.LE)),by=.(repn,intervention)]

(PSARg <- PSAR[,.(ep=mean(ep),ei=mean(ei),ed=mean(ed),el=mean(el)*1e-6),by=intervention])
## 50:50 cop:inc and ~110K d

(PSARg[1,ed] - PSARg[2,ed])             #will be a bit lower with CDR correction
## TODO why is ei about the same?? 
## TODO IPT RR wrong way? check

