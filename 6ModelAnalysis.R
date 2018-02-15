## TODO
## HIV functions
## TODO courses treats screens
## see end file TODO 
rm(list=ls())                           #clear
## next 4 for formatting & outputs
library(officer)
library(magrittr)
library(flextable)
pp <- function(x,sf=3,ns=0) format(signif(round(x),sf), nsmall=ns, big.mark=",",scientific = FALSE)
## load decision tree model and functions
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
PSA[,progn.LP.PTn:=pprogn*rrtst/(1+rrtst)] #TB incidence in LTBI +ve PT-ve
PSA[,progn.LN.PTn:=pprogn*1/(1+rrtst)]     #TB incidence in LTBI -ve PT-ve
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
PSO[,progn.LP.PTn:=pprogn*rrtst/(1+rrtst)] #TB incidence in LTBI +ve PT-ve
PSO[,progn.LN.PTn:=pprogn*1/(1+rrtst)]     #TB incidence in LTBI -ve PT-ve
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
PSA[,hiv:=0]
PSA[,art:=0]


## ==========  HIV/ART
load('data/DL.Rdata')
unique(DL[,.(iso3,hivprop,artprop)])

PSAH <- copy(PSA)
PSAH[,hiv:=1]
PSAHa <- copy(PSAH)
PSAHa[,art:=1]
PSAH <- rbind(PSAH,PSAHa)
rm(PSAHa)

## ========================================
## ------ tree model calculations
F <- makeTfuns(kexp,unique(kexp$fieldsAll))
names(F)
## getAQ(kexp,'LE')
summary(F$checkfun(PSA))                         #check

## TODO include hhc variance

PSA$e.prevalent <- F$prevalentfun(PSA)
PSA$e.incidence <- F$incidencefun(PSA)
PSA$e.deaths <- F$deathfun(PSA)
PSA$e.LE <- F$LEfun(PSA)
PSA$e.ATT <- F$treatmentsfun(PSA)
PSA$e.IPT <- F$IPTfun(PSA)
PSA$e.LTBI <- F$LTBIfun(PSA)
PSA$e.ATTprev <- F$ATTprevfun(PSA)

## multiply by number of children
PSA[,ehhc:=hhc]                         #TODO change to stochastic
ests <- grep('e\\.',names(PSA),value=TRUE)
nest <- length(ests)
PSA[,c(ests):=lapply(.SD,function(x) x*ehhc),.SDcols=ests] # needing to be: x HHC
## add this and visits to ests list to deal with generically
PSA[,e.hhc:=hhc]          #HH contacts join into things dealt with generically
load('data/HHV.Rdata')    #visit (notification) data
PSA <- merge(PSA,HHV[,.(iso3,e.households=visits)],all.x = TRUE) #merge visits in
PSA[acat==unique(acat)[1],e.households:=0]                       #avoid double counting
PSA[!grepl('5',intervention),e.households:=0]                    #no visits
ests <- grep('e\\.',names(PSA),value=TRUE)                       #regrab things to est
nest <- length(ests)


## ## TODO IPT read more

## ==== full results table ===
## global/regional/age
PSAG <- PSA[,lapply(.SD,sum),by=.(repn,intervention),.SDcols=ests]
PSAR <- PSA[,lapply(.SD,sum),by=.(repn,intervention,g_whoregion),.SDcols=ests]
PSAA <- PSA[,lapply(.SD,sum),by=.(repn,intervention,acat),.SDcols=ests]

## == diffs
## global
tmp1 <- PSAG[,.SD[intervention=='Under 5 & HIV+ve'] - .SD[intervention=='No intervention'],.SDcols=3:(2+nest),by=repn]
tmp1[,intervention:='B-A']
tmp <- PSAG[,.SD[intervention=='Under 5 & HIV+ve & LTBI+'] - .SD[intervention=='No intervention'],.SDcols=3:(2+nest),by=repn]
tmp[,intervention:='C-A']
PSAG <- rbind(PSAG,tmp1,tmp)

## regional
tmp1 <- PSAR[,.SD[intervention=='Under 5 & HIV+ve'] - .SD[intervention=='No intervention'],.SDcols=4:(3+nest),by=.(repn,g_whoregion)]
tmp1[,intervention:='B-A']
tmp <- PSAR[,.SD[intervention=='Under 5 & HIV+ve & LTBI+'] - .SD[intervention=='No intervention'],.SDcols=4:(3+nest),by=.(repn,g_whoregion)]
tmp[,intervention:='C-A']
PSAR <- rbind(PSAR,tmp1,tmp)

## age
tmp1 <- PSAA[,.SD[intervention=='Under 5 & HIV+ve'] - .SD[intervention=='No intervention'],.SDcols=4:(3+nest),by=.(repn,acat)]
tmp1[,intervention:='B-A']
tmp <- PSAA[,.SD[intervention=='Under 5 & HIV+ve & LTBI+'] - .SD[intervention=='No intervention'],.SDcols=4:(3+nest),by=.(repn,acat)]
tmp[,intervention:='C-A']
PSAA <- rbind(PSAA,tmp1,tmp)

## == summary
PSAGm <- PSAG[,lapply(.SD,mean),by=.(intervention),.SDcols=ests]
PSARm <- PSAR[,lapply(.SD,mean),by=.(intervention,g_whoregion),.SDcols=ests]
PSAAm <- PSAA[,lapply(.SD,mean),by=.(intervention,acat),.SDcols=ests]

## TODO uncertainty measures here

## == gathering
intl <- c("No intervention","Under 5 & HIV+ve","Under 5 & HIV+ve & LTBI+","B-A","C-A")
varlv <- ests[c(10,9,7,1,8,5,6,2,3,4)]  #reorder
nint <- length(intl)
hhrpl <- PSA[repn==1 & acat=="[5,15)" & intervention==unique(intervention)[2],sum(e.households)] #the non zero value

## global
RTg <- melt(PSAGm,id='intervention')
## regional
RTr <- melt(PSARm,id=c('intervention','g_whoregion'))
RTr <- dcast(RTr,intervention+variable ~ g_whoregion,value.var = 'value')
## age
RTa <- melt(PSAAm,id=c('intervention','acat'))
RTa <- dcast(RTa,intervention+variable ~ acat,value.var = 'value')

## merge
RTg$intervention <- factor(RTg$intervention,levels=intl,ordered=TRUE); RTg$variable <- factor(RTg$variable,levels=varlv,ordered=TRUE)
RTr$intervention <- factor(RTr$intervention,levels=intl,ordered=TRUE); RTr$variable <- factor(RTr$variable,levels=varlv,ordered=TRUE)
RTa$intervention <- factor(RTa$intervention,levels=intl,ordered=TRUE); RTa$variable <- factor(RTa$variable,levels=varlv,ordered=TRUE)
RTg <- RTg[order(intervention,variable),]; RTr <- RTr[order(intervention,variable),]; RTa <- RTa[order(intervention,variable),]
names(RTg)[3] <- 'Global'
RT <- cbind(RTg,RTr[,-c(1:2),with=FALSE],RTa[,-c(1:2),with=FALSE])

## intervention-related values for hhc & hh visits
RT[!grepl('5|-A',intervention) & variable=='e.hhc',(3:ncol(RT)):=0] #only for B & C
RT[grepl('5|-A',intervention) & variable=='e.households',(3:ncol(RT)):=hhrpl]
RT[!grepl('5|-A',intervention) & variable=='e.households',(3:ncol(RT)):=0]

## ## NNx TODO redo
## x <- RT[intervention=='B-A',Global]
## (ba <- c(hhn=-x[1]/x[6],hcn=-x[2]/x[6],ptn=-x[3]/x[6],txn=-x[4]/x[6]))
## x <- RT[intervention=='C-A',Global]
## (ca <- c(hhn=-x[1]/x[6],hcn=-x[2]/x[6],ptn=-x[3]/x[6],txn=-x[4]/x[6]))


## == formatting & output
## formats
whoz <- !RT$variable %in% c('e.ATTprev')
tmp <- RT[whoz,lapply(.SD,pp),.SDcols=3:ncol(RT)]
tmp <- cbind(RT[whoz,1:2,with=FALSE],tmp)
intz <- as.character(tmp$intervention)
uintz <- unique(intz)
intz[intz==uintz[1]] <- paste0('A: ',intz[intz==uintz[1]])
intz[intz==uintz[2]] <- paste0('B: ',intz[intz==uintz[2]])
intz[intz==uintz[3]] <- paste0('C: ',intz[intz==uintz[3]])
tmp$intervention <- factor(intz)
myft <- regulartable(tmp)
myft <- merge_v(myft, j = c("intervention", "variable") )
## myft <- autofit(myft)
myft


fwrite(RT,file='tables/RT.csv')

## simpler output for main article
RTS <- RT[,.(intervention,variable,Global,`[0,5)`,`[5,15)`)]
names(RTS)[2] <- 'quantity'
RTS <- melt(RTS,id.vars = c('intervention','quantity'))
RTS <- dcast(RTS,quantity + variable ~ intervention,value.var = 'value')
RTS$variable <- factor(RTS$variable,levels=c('[0,5)','[5,15)','Global'),ordered=TRUE)
RTS <- RTS[order(variable,quantity),]
tmp <- RTS[,lapply(.SD,pp),.SDcols=3:ncol(RTS)]
tmp <- cbind(RTS[,2:1,with=FALSE],tmp)
myft2 <- regulartable(tmp)
myft2 <- merge_v(myft2, j = c("variable","quantity") )
## myft <- autofit(myft)
myft2


read_docx() %>%
  body_add_par(value = "Global and age outputs", style = "heading 1") %>%
  ## body_add_table(value = myft, style = "table_template" ) %>%
  body_add_flextable(value = myft2) %>%
  body_end_section(continuous = FALSE, landscape = TRUE) %>% 
  print(target = "tables/RTS.docx") %>% 
    invisible()

## graph expts
names(RT)[2] <- 'quantity'
RTM <- melt(RT,id.vars = c('intervention','quantity'))

ggplot(data=RTM[!grepl('-A',intervention) &
                quantity %in% c('e.ATT','e.IPT','e.incidence','e.deaths','e.LE') &
                variable %in% c('AFR','EUR','WPR','SEA','WPO','AMR')],
       aes(x=variable,y=value,fill=intervention)) +
  geom_bar(stat='identity',position='dodge') + 
  coord_flip() + facet_wrap(~quantity,scales = 'free_x')

## ## or col by region:
## ggplot(data=RTM[!grepl('-A',intervention) &
##                 quantity %in% c('e.ATT','e.IPT','e.incidence','e.deaths','e.LE') &
##                 variable %in% c('AFR','EUR','WPR','SEA','WPO','AMR')],
##        aes(x=intervention,y=value,fill=variable)) +
##   geom_bar(stat='identity',position='dodge') + 
##   coord_flip() + facet_wrap(~quantity,scales = 'free_x')


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
## CHECK CDR for incident cases?
## consider different TST for HIV+
