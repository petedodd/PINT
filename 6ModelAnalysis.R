## TODO
## HIV functions
## TODO courses treats screens
## see end file TODO 
rm(list=ls())
source('3ModelDefn.R')    #load the decision tree model logic/definition
source('4ModelFuns.R')    #function defns for decision tree (and distribution parameters)

nrep <- 1e2                             #number of reps for PSA

## ======= children 0-4
load('data/DLC.Rdata')                  #parent data for children 0-4
PSA <- DLC[rep(1:nrow(DLC),nrep),]      #build PSA data frame
PSA[,repn:=rep(1:nrep,each=nrow(DLC))]

## compute variables for PSA data table
azu5 <- rep(1,nrow(PSA))                #dummy ages for functions defined like that
PSA[,CDR:=CDR(cdr04,cdr04ab)]           #CDR
PSA[,CFRtxY:=CFRtxY(azu5)]              #CDR on ATT
PSA[,CFRtxN:=CFRtxN(azu5)]              #CDR not on ATT
PSA[,coprev:=coprev(azu5)]              #coprevalent TB
PSA[,IPTrr:=IPTrr(azu5)]                #IPT RR for incident TB
PSA[,ltbi.prev:=ltbi.prev(azu5,coprev)] #LTBI prevalence
PSA[,pprogn:=avu5progprob(a1,a2,a3,a4,a5,LAT)] #progression probability (averaged 0-4)
PSA[,rrtst:=RRtst(azu5)]                #RR for incidence if TST+ve
PSA[,inc:=ltbi.prev * pprogn]           #TB incidence, total
PSA[,progn.LP.PTn:=inc]                 #TB incidence in LTBI +ve PT-ve
PSA[,progn.LN.PTn:=0]                   #TB incidence in LTBI -ve PT-ve
PSA[,progn.LP.PTp:=progn.LP.PTn*IPTrr]  #TB incidence in LTBI +ve PT+ve
PSA[,progn.LN.PTp:=progn.LN.PTn*IPTrr]  #TB incidence in LTBI -ve PT+ve
PSA[,PTcov.N:=0]                        #coverage of PT in LTBI -ve
PSA[,PTcov.P:=0]                        #coverage of PT in LTBI +ve

## intervention set
npsa <- nrow(PSA)
PSA <- PSA[rep(1:npsa,3),]              #replicates by intervention
PSA[,intervention:=c(rep('No intervention',npsa),
                     rep('Under 5 & HIV+ve',npsa),
                     rep('Under 5 & HIV+ve & LTBI+',npsa))]
PSA[intervention!='No intervention',PTcov.N:=1]
PSA[intervention!='No intervention',PTcov.P:=1]
PSA[intervention!='No intervention',CDR:=1]
PSA[,acat:="[0,5)"]                     #age group

## ======= children 5-14
load('data/DLO.Rdata')                  #parent data for children 5-14
PSO <- DLO[rep(1:nrow(DLO),nrep),]      #build PSA data frame
PSO[,repn:=rep(1:nrep,each=nrow(DLO))]

## compute variables for PSA data table
azu5 <- rep(10,nrow(PSO))                #dummy ages for functions defined like that
PSO[,CDR:=CDR(cdr514,cdr514ab)]         #CDR
PSO[,CFRtxY:=CFRtxY(azu5)]              #CDR on ATT
PSO[,CFRtxN:=CFRtxN(azu5)]              #CDR not on ATT
PSO[,coprev:=coprev(azu5)]              #coprevalent TB
PSO[,IPTrr:=IPTrr(azu5)]                #IPT RR for incident TB
PSO[,ltbi.prev:=ltbi.prev(azu5,coprev)] #LTBI prevalence
PSO[,pprogn:=avo5progprob(a6,a7,a8,a9,a10,
                          a11,a12,a13,a14,a15,
                          LAT)] #progression probability (averaged 0-4)
PSO[,rrtst:=RRtst(azu5)]                #RR for incidence if TST+ve
PSO[,inc:=ltbi.prev * pprogn]           #TB incidence, total
PSO[,progn.LP.PTn:=inc]                 #TB incidence in LTBI +ve PT-ve

PSO[,progn.LN.PTn:=0]                   #TB incidence in LTBI -ve PT-ve
PSO[,progn.LP.PTp:=progn.LP.PTn*IPTrr]  #TB incidence in LTBI +ve PT+ve
PSO[,progn.LN.PTp:=progn.LN.PTn*IPTrr]  #TB incidence in LTBI -ve PT+ve
PSO[,PTcov.N:=0]                        #coverage of PT in LTBI -ve
PSO[,PTcov.P:=0]                        #coverage of PT in LTBI +ve

## intervention set
PSO <- PSO[rep(1:npsa,3),]              #replicates by intervention
PSO[,intervention:=c(rep('No intervention',npsa),
                     rep('Under 5 & HIV+ve',npsa),
                     rep('Under 5 & HIV+ve & LTBI+',npsa))]
PSO[intervention!='No intervention',PTcov.N:=1]
PSO[intervention!='No intervention',PTcov.P:=1]
PSO[intervention!='No intervention',CDR:=1]

## other changes TODO

PSO[,acat:="[5,15)"]                    #age group

## ==== join
## ditch BCG coverage by age now
PSA[,a1:=NULL];PSA[,a2:=NULL];PSA[,a3:=NULL];PSA[,a4:=NULL];PSA[,a5:=NULL];
PSO[,a6:=NULL];PSO[,a7:=NULL];PSO[,a8:=NULL];PSO[,a9:=NULL];PSO[,a10:=NULL];
PSO[,a11:=NULL];PSO[,a12:=NULL];PSO[,a13:=NULL];PSO[,a14:=NULL];PSO[,a15:=NULL];
PSA[,cdr04:=NULL]; PSA[,cdr04ab:=NULL]; PSO[,cdr514:=NULL]; PSO[,cdr514ab:=NULL];
PSA[,u5hhc.l:=NULL]; PSA[,u5hhc.sdl:=NULL]; PSO[,o5hhc.l:=NULL]; PSO[,o5hhc.sdl:=NULL];
names(PSA)[5:6] <- names(PSO)[5:6] <- c('hhc','hhc.sd')

PSA <- rbind(PSA,PSO)
rm(PSO)
PSA

## ========================================
## ------ tree model calculations
F <- makeTfuns(kexp,unique(kexp$fieldsAll))
names(F)
## getAQ(kexp,'LE')
summary(F$checkfun(PSA))                         #!!

## summary(F$prevalentfun(PSA))
## summary(F$incidencefun(PSA))
## summary(F$deathfun(PSA))
## summary(F$LEfun(PSA))

## TODO include hhc variance

PSA$e.prevalent <- F$prevalentfun(PSA)*PSA$hhc
PSA$e.incidence <- F$incidencefun(PSA)*PSA$hhc
PSA$e.deaths <- F$deathfun(PSA)*PSA$hhc
PSA$e.LE <- F$LEfun(PSA)*PSA$hhc
PSA$e.txs <- F$treatmentsfun(PSA)*PSA$hhc

## print(kexp,'incidence','LE')
## print(kexp,'death','treatments')
## print(kexp,'prevalent','LTBI')

## plotter(kexp, varz=c('name','LE'), edgelabel = TRUE)
## plotter(kexp, varz=c('name','incidence'), edgelabel = TRUE)
print(kexp,'incidence')



## comparison by category
PSAR <- PSA[,.(ep=sum(e.prevalent),
               ei=sum(e.incidence),
               ed=sum(e.deaths),
               et=sum(e.txs),
               el=sum(e.LE)),by=.(repn,intervention)]


pp <- function(x,sf=3,ns=0) format(signif(round(x),sf), nsmall=ns, big.mark=",")

(PSARg <- PSAR[,.(ep=pp(mean(ep)),ei=pp(mean(ei)),
                  ed=pp(mean(ed)),ed.sd=pp(sd(ed)),
                  et=pp(mean(et)),et.sd=pp(sd(et)),
                  el=pp(mean(el)*1e-6)),by=intervention])
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
