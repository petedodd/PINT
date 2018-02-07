## TODO
## HIV functions
## TODO courses treats screens
## see end file TODO 
rm(list=ls())
source('3ModelDefn.R')
source('4ModelFuns.R')


## ======= children 0-4
load('data/DLC.Rdata')                  #parent data

## build PSA data frame
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
names(F)
## getAQ(kexp,'LE')
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
print(kexp,'incidence')


## comparison by category
PSAR <- PSA[,.(ep=sum(e.prevalent),
               ei=sum(e.incidence),
               ed=sum(e.deaths),
               el=sum(e.LE)),by=.(repn,intervention)]

(PSARg <- PSAR[,.(ep=mean(ep),ei=mean(ei),
                  ed=mean(ed),ed.sd=sd(ed),
                  el=mean(el)*1e-6),by=intervention])
## 50:50 cop:inc and ~110K d

(PSARg[1,ed] - PSARg[2,ed])             #will be a bit lower with CDR correction
(PSARg[1,el] - PSARg[2,el])             #will be a bit lower with CDR correction
## TODO IPT read more
## TODO check NAs

PSARg

summary(PSA)

## TODO
## u5 intervention and outputs
## number of tracing and courses IPT
## o5 intervention and outputs
## HIV functions
## document assumptions esp CDR
## CY compare
