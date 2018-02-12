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
PSA[,progn.LP.PTn:=inc*rrtst/(1+rrtst)] #TB incidence in LTBI +ve PT-ve
PSA[,progn.LN.PTn:=inc*1/(1+rrtst)]     #TB incidence in LTBI -ve PT-ve
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
PSA[intervention!='No intervention',CDR:=1] #screening
PSA[intervention!='No intervention',PTcov.N:=1] #IPT
PSA[intervention!='No intervention',PTcov.P:=1] #IPT
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
PSO[,progn.LP.PTn:=inc*rrtst/(1+rrtst)] #TB incidence in LTBI +ve PT-ve
PSO[,progn.LN.PTn:=inc*1/(1+rrtst)]     #TB incidence in LTBI -ve PT-ve
PSO[,progn.LP.PTp:=progn.LP.PTn*IPTrr]  #TB incidence in LTBI +ve PT+ve
PSO[,progn.LN.PTp:=progn.LN.PTn*IPTrr]  #TB incidence in LTBI -ve PT+ve
PSO[,PTcov.N:=0]                        #coverage of PT in LTBI -ve
PSO[,PTcov.P:=0]                        #coverage of PT in LTBI +ve

## intervention set
PSO <- PSO[rep(1:npsa,3),]              #replicates by intervention
PSO[,intervention:=c(rep('No intervention',npsa),
                     rep('Under 5 & HIV+ve',npsa),
                     rep('Under 5 & HIV+ve & LTBI+',npsa))]
PSO[intervention!='No intervention',CDR:=1] #screening
PSO[intervention=='Under 5 & HIV+ve & LTBI+',PTcov.P:=1] #PT for TST+


## other changes once HIV included TODO

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

## ==========  HIV/ART
load('data/DL.Rdata')
unique(DL[,.(iso3,hivprop,artprop)])


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
PSA$e.IPT <- F$IPTfun(PSA)*PSA$hhc

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
               ept=sum(e.IPT),
               el=sum(e.LE)),by=.(repn,intervention)]


pp <- function(x,sf=3,ns=0) format(signif(round(x),sf), nsmall=ns, big.mark=",")
 
(PSARg <- PSAR[,.(ep=pp(mean(ep)),ei=pp(mean(ei)),
                  ed=pp(mean(ed)),ed.sd=pp(sd(ed)),
                  et=pp(mean(et)),et.sd=pp(sd(et)),
                  ept=pp(mean(ept)),ept.sd=pp(sd(ept)),
                  el=pp(mean(el)*1e-6)),by=intervention])




## (PSARg[1,ed] - PSARg[2,ed])             #will be a bit lower with CDR correction
## (PSARg[1,el] - PSARg[2,el])             #will be a bit lower with CDR correction
## ## TODO IPT read more
## ## TODO check NAs
## PSARg
## summary(PSA)

## ==== full results table w/o uncertainty for now ===
## global
PSAG <- PSA[,.(ei=sum(e.incidence),
               ed=sum(e.deaths),
               et=sum(e.txs),
               ept=sum(e.IPT),
               el=sum(e.LE)),by=.(repn,intervention)]
## regional
PSAR <- PSA[,.(ei=sum(e.incidence),
               ed=sum(e.deaths),
               et=sum(e.txs),
               ept=sum(e.IPT),
               el=sum(e.LE)),by=.(repn,intervention,g_whoregion)]
## age
PSAA <- PSA[,.(ei=sum(e.incidence),
               ed=sum(e.deaths),
               et=sum(e.txs),
               ept=sum(e.IPT),
               el=sum(e.LE)),by=.(repn,intervention,acat)]


## == diffs
## global
tmp1 <- PSAG[,.SD[intervention=='Under 5 & HIV+ve'] - .SD[intervention=='No intervention'],.SDcols=3:7,by=repn]
tmp1[,intervention:='B-A']
tmp <- PSAG[,.SD[intervention=='Under 5 & HIV+ve & LTBI+'] - .SD[intervention=='No intervention'],.SDcols=3:7,by=repn]
tmp[,intervention:='C-A']
PSAG <- rbind(PSAG,tmp1,tmp)

## regional
tmp1 <- PSAR[,.SD[intervention=='Under 5 & HIV+ve'] - .SD[intervention=='No intervention'],.SDcols=4:8,by=.(repn,g_whoregion)]
tmp1[,intervention:='B-A']
tmp <- PSAR[,.SD[intervention=='Under 5 & HIV+ve & LTBI+'] - .SD[intervention=='No intervention'],.SDcols=4:8,by=.(repn,g_whoregion)]
tmp[,intervention:='C-A']
PSAR <- rbind(PSAR,tmp1,tmp)

## age
tmp1 <- PSAA[,.SD[intervention=='Under 5 & HIV+ve'] - .SD[intervention=='No intervention'],.SDcols=4:8,by=.(repn,acat)]
tmp1[,intervention:='B-A']
tmp <- PSAA[,.SD[intervention=='Under 5 & HIV+ve & LTBI+'] - .SD[intervention=='No intervention'],.SDcols=4:8,by=.(repn,acat)]
tmp[,intervention:='C-A']
PSAA <- rbind(PSAA,tmp1,tmp)

## == summary
PSAGm <- PSAG[,.(incidence=mean(ei),deaths=mean(ed),ATT=mean(et),IPT=mean(ept),kLY=mean(el)*1e-3),by=.(intervention)]
PSARm <- PSAR[,.(incidence=mean(ei),deaths=mean(ed),ATT=mean(et),IPT=mean(ept),kLY=mean(el)*1e-3),by=.(intervention,g_whoregion)]
PSAAm <- PSAA[,.(incidence=mean(ei),deaths=mean(ed),ATT=mean(et),IPT=mean(ept),kLY=mean(el)*1e-3),by=.(intervention,acat)]

## == visits, children
## - visits
load('data/HHV.Rdata')
HHV[,sum(visits)]
HHVR <- HHV[,.(visits=sum(visits)),by=g_whoregion]
## PSAGm[,visits:=HHV[,sum(visits)]];
## PSAGm[intervention %in% c('No intervention','B-A','C-A'),visits:=0]
## PSARm <- merge(PSARm,HHV[,.(visits=sum(visits)),by=g_whoregion],by='g_whoregion',all.x = TRUE)
## PSARm[intervention %in% c('No intervention','B-A','C-A'),visits:=0]

## - children TODO uncertainty
HCg <- PSA[,.(hhc=sum(hhc)),by=repn]
HCr <- PSA[,.(hhc=sum(hhc)),by=.(repn,g_whoregion)]
HCa <- PSA[,.(hhc=sum(hhc)),by=.(repn,acat)]
HCg <- HCg[,.(hhc=mean(hhc))]
HCr <- HCr[,.(hhc=mean(hhc)),by=g_whoregion]
HCa <- HCa[,.(hhc=mean(hhc)),by=acat]

## == gathering
intl <- c("No intervention","Under 5 & HIV+ve","Under 5 & HIV+ve & LTBI+","B-A","C-A")
varlv <- c('households','contacts','IPT','ATT','incidence','deaths','kLY')
nint <- length(intl)

## global
RTg <- melt(PSAGm,id='intervention')
RTx <- data.table(intervention=rep(intl,2),
                  variable=rep(c('households','contacts'),each=nint),
                  value=rep(0,2*nint))  #visits & contacts
RTx[grepl('5|-A',intervention) & variable=='households',value:=HHV[,sum(visits)]]
RTx[grepl('5|-A',intervention) & variable=='contacts',value:=HCg[,hhc]]
RTg <- rbind(RTg,RTx)

## regional
RTr <- melt(PSARm,id=c('intervention','g_whoregion'))

RTx <- data.table(intervention=rep(intl,2*6),
                  g_whoregion=rep(RTr[,unique(g_whoregion)],each=nint*2),
                  variable=rep(c('households','contacts'),each=nint),
                  value=rep(0,2*nint*6))  #visits & contacts
for(reg in RTx[,unique(g_whoregion)]) RTx[grepl('5|-A',intervention) & variable=='households' & g_whoregion==reg,value:=HHVR[g_whoregion==reg,sum(visits)]]
for(reg in RTx[,unique(g_whoregion)]) RTx[grepl('5|-A',intervention) & variable=='contacts' & g_whoregion==reg,value:=HCr[g_whoregion==reg,hhc]]
RTr <- rbind(RTr,RTx)
RTr <- dcast(RTr,intervention+variable ~ g_whoregion,value.var = 'value')

## age
RTa <- melt(PSAAm,id=c('intervention','acat'))
RTx <- data.table(intervention=rep(RTg[,unique(intervention)],2*2),
                  variable=rep(c('households','contacts'),each=nint),
                  acat=rep(RTa[,unique(acat)],each=nint*2),
                  value=rep(0,4*nint))  #visits & contacts
for(reg in RTx[,unique(acat)]) RTx[grepl('5|-A',intervention) & variable=='households' & acat==reg,value:=HHV[,sum(visits)]]
for(reg in RTx[,unique(acat)]) RTx[grepl('5|-A',intervention) & variable=='contacts' & acat==reg,value:=HCa[acat==reg,hhc]]
RTa <- rbind(RTa,RTx)
RTa <- dcast(RTa,intervention+variable ~ acat,value.var = 'value')

## merge
RTg$intervention <- factor(RTg$intervention,levels=intl,ordered=TRUE); RTg$variable <- factor(RTg$variable,levels=varlv,ordered=TRUE)
RTr$intervention <- factor(RTr$intervention,levels=intl,ordered=TRUE); RTr$variable <- factor(RTr$variable,levels=varlv,ordered=TRUE)
RTa$intervention <- factor(RTa$intervention,levels=intl,ordered=TRUE); RTa$variable <- factor(RTa$variable,levels=varlv,ordered=TRUE)
RTg <- RTg[order(intervention,variable),]; RTr <- RTr[order(intervention,variable),]; RTa <- RTa[order(intervention,variable),]
names(RTg)[3] <- 'Global'
RT <- cbind(RTg,RTr[,-c(1:2),with=FALSE],RTa[,-c(1:2),with=FALSE])

## NNx
x <- RT[intervention=='B-A',Global]
(ba <- c(hhn=-x[1]/x[6],hcn=-x[2]/x[6],ptn=-x[3]/x[6],txn=-x[4]/x[6]))
x <- RT[intervention=='C-A',Global]
(ca <- c(hhn=-x[1]/x[6],hcn=-x[2]/x[6],ptn=-x[3]/x[6],txn=-x[4]/x[6]))


## == formatting & output
library(officer)
library(magrittr)
library(flextable)

fwrite(RT,file='tables/RT.csv')

## formats
tmp <- RT[,lapply(.SD,pp),.SDcols=3:11]
tmp <- cbind(RT[,1:2,with=FALSE],tmp)
intz <- as.character(tmp$intervention)
intz[1:7] <- paste0('A: ',intz[1:7])
intz[1:7 + 7] <- paste0('B: ',intz[1:7 + 7])
intz[1:7 + 14] <- paste0('C: ',intz[1:7 + 14])
tmp$intervention <- factor(intz)

myft <- regulartable(tmp)
myft <- merge_v(myft, j = c("intervention", "variable") )
## myft <- autofit(myft)
myft

read_docx() %>%
  body_add_par(value = "Global, regional, and age outputs", style = "heading 1") %>%
  ## body_add_table(value = myft, style = "table_template" ) %>%
  body_add_flextable(value = myft) %>%
  body_end_section(continuous = FALSE, landscape = TRUE) %>% 
  print(target = "tables/RT.docx") %>% 
    invisible()

## all country output too

## ## writing out
## read_docx() %>%
##   body_add_par(value = "Global output table", style = "heading 1") %>%
##   body_add_table(value = PSARg, style = "table_template" ) %>%
##   print(target = "tables/globaloutput.docx") %>% 
##     invisible()

## write.csv(PSARg,file='tables/globaloutput.csv')

## TODO
## HIV functions
## document assumptions esp CDR
## CY compare
## uncertainty
## hh uncertainty
