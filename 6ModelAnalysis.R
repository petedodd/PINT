## This file is the modelling analysis. In addition to the libraries at the top of
## this file, the HEdtree and data.table libraries are needed (see the top of file 3).
## Once installed, these do not need separate loading as they are loaded when file 3 is
## sourced below. Outputs are sent to tables/
rm(list=ls())                           #clear
## next 4 for formatting & outputs
library(officer)
library(magrittr)
library(flextable)
pp <- function(x,sf=3,ns=0) format(signif(round(x),sf), nsmall=ns, big.mark=",",scientific = FALSE)
pps <- function(x,...) gsub("^\\s+|\\s+$", "", pp(x,...))
ppb <- function(x,y,...) paste0('(',pps(x,...),' \u2013 ',pps(y,...),')')
uq <- function(x)quantile(x,.75)
lq <- function(x)quantile(x,.25)

## load decision tree model and functions
source('3ModelDefn.R')    #load the decision tree model logic/definition
source('4ModelFuns.R')    #function defns for decision tree (and distribution parameters)

nrep <- 1e2                             #number of reps for PSA
## ======= children 0-4
load('data/DLC.Rdata')                  #parent data for children 0-4
PSA <- DLC[rep(1:nrow(DLC),nrep),]      #build PSA data frame
PSA[,repn:=rep(1:nrep,each=nrow(DLC))]
## ======= children 5-14
load('data/DLO.Rdata')                  #parent data for children 5-14
PSO <- DLO[rep(1:nrow(DLO),nrep),]      #build PSA data frame
PSO[,repn:=rep(1:nrep,each=nrow(DLO))]

## variables needed for PSA
PSA[,az:=rep(1,nrow(PSA))]               #dummy ages for functions defined like that
PSA[,cm:=cdr04]
PSA[,cab:=cdr04ab]
PSA[,acat:="[0,5)"]                     #age group
PSO[,az:=rep(10,nrow(PSO))]                #dummy ages for functions defined like that
PSO[,cm:=cdr514]
PSO[,cab:=cdr514ab]
PSO[,acat:="[5,15)"]                     #age group
## join these
names(PSA)
names(PSO)
PSA <- rbind(PSA,PSO,fill=TRUE)
rm(PSO)
## HIV copies of this data
PSAh <- copy(PSA)
PSAa <- copy(PSA)
PSA[,hiv:=0]; PSAh[,hiv:=1]; PSAa[,hiv:=1]
PSA[,art:=0]; PSAh[,art:=0]; PSAa[,art:=1]
PSA <- rbind(PSA,PSAh,PSAa)
rm(PSAh,PSAa)

## compute variables for PSA data table
PSA[,CDR:=CDR(cm,cab)]           #CDR
PSA[,coprev:=coprev(az)]                #coprevalent TB
PSA[,IPTrr:=IPTrr(az)]                  #IPT RR for incident TB
PSA[,ltbi.prev:=ltbi.prev(az,coprev)]   #LTBI prevalence
PSA[,rrtst:=RRtst(az)]                  #RR for incidence if TST+ve
PSA[,CFRtxY:=CFRtxY(az,hiv,art)]                #CFR on ATT
PSA[,CFRtxN:=CFRtxN(az,hiv,art)]                #CFR not on ATT
PSA[acat=="[0,5)",
    pprogn:=avu5progprob(a1,a2,a3,a4,a5,LAT,hiv,art)] #progression probability (averaged 0-4)
PSA[acat=="[5,15)",
    pprogn:=avo5progprob(a6,a7,a8,a9,a10,
                         a11,a12,a13,a14,a15,
                         LAT,hiv,art)] #progression probability (averaged 5-14)
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
PSA[acat=="[0,5)" & intervention!='No intervention',PTcov.N:=1] #IPT for all under 5s
PSA[acat=="[0,5)" & intervention!='No intervention',PTcov.P:=1] #IPT for all under 5s
PSA[acat=="[5,15)" & hiv==1 & intervention!='No intervention',PTcov.N:=1] #IPT o5s HIV+
PSA[acat=="[5,15)" & hiv==1 & intervention!='No intervention',PTcov.P:=1] #IPT o5s HIV+
PSA[acat=="[5,15)" & intervention=='Under 5 & HIV+ve & LTBI+',PTcov.P:=1] #IPT for o5s L+


## tidying
PSA[,c(paste0('a',1:15)):=NULL]
PSA[,c('cdr04','cdr04ab','cdr514','cdr514ab'):=NULL]
PSA[acat=="[0,5)",hhc:=u5hhc]; PSA[acat=="[0,5)",hhc.sd:=u5hhc.sd]
PSA[acat=="[5,15)",hhc:=o5hhc]; PSA[acat=="[5,15)",hhc.sd:=o5hhc.sd]
PSA[,c('u5hhc.l','u5hhc.sdl','o5hhc.l','o5hhc.sdl','u5hhc','u5hhc.sd','o5hhc','o5hhc.sd'):=NULL]
PSA[,c('az','cm','cab'):=NULL]




## ========================================
## ------ tree model calculations
F <- makeTfuns(kexp,unique(kexp$fieldsAll))
names(F)
## getAQ(kexp,'LE')
summary(F$checkfun(PSA))                         #check

PSA$e.prevalent <- F$prevalentfun(PSA)
PSA$e.incidence <- F$incidencefun(PSA)
PSA$e.deaths <- F$deathfun(PSA)
PSA$e.LE <- F$LEfun(PSA)
PSA$e.ATT <- F$treatmentsfun(PSA)
PSA$e.IPT <- F$IPTfun(PSA)
PSA$e.LTBI <- F$LTBIfun(PSA)
PSA$e.ATTprev <- F$ATTprevfun(PSA)


## multiply by number of children
## TODO include hhc variance
PSA[,ehhc:=hhc]                         #TODO change to stochastic
## HIV splits
load('data/DL.Rdata')
HIV <- unique(DL[,.(iso3,hivprop,artprop)])
PSA <- merge(PSA,HIV,by='iso3',all.x=TRUE)
PSA[hiv==0,ehhc:=ehhc*(1-hivprop)]
PSA[hiv==1 & art==0,ehhc:=ehhc*hivprop*(1-artprop)]
PSA[hiv==1 & art==1,ehhc:=ehhc*hivprop*artprop]
PSA[,c('hivprop','artprop'):=NULL]
ests <- grep('e\\.',names(PSA),value=TRUE)
nest <- length(ests)

PSA[,c(ests):=lapply(.SD,function(x) x*ehhc),.SDcols=ests] # needing to be: x HHC
## add this and visits to ests list to deal with generically
PSA[,e.hhc:=ehhc]          #HH contacts join into things dealt with generically
PSA[!grepl('5|-A',intervention),e.hhc:=0] #no intervention
load('data/HHV.Rdata')    #visit (notification) data
PSA <- merge(PSA,HHV[,.(iso3,e.households=visits)],all.x = TRUE) #merge visits in
PSA[acat==unique(acat)[1],e.households:=0]                       #avoid double counting
PSA[!grepl('5',intervention),e.households:=0]                    #no visits
PSA[acat=='[0,5)',e.households:=PSA[acat=='[5,15)',e.households]] #add for younger
PSA[,e.households:=e.households/6] # avoid x6 need x2 for age split
ests <- grep('e\\.',names(PSA),value=TRUE)                       #regrab things to est
nest <- length(ests)


## ## TODO IPT read more

## ==== full results table ===
## global/regional/age
PSAG <- PSA[,lapply(.SD,sum),by=.(repn,intervention),.SDcols=ests]
PSAR <- PSA[,lapply(.SD,sum),by=.(repn,intervention,g_whoregion),.SDcols=ests]
PSAA <- PSA[,lapply(.SD,sum),by=.(repn,intervention,acat),.SDcols=ests]
PSAA[,e.households:=e.households*2]     #same for each age cat really

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
## means
PSAGm <- PSAG[,lapply(.SD,mean),by=.(intervention),.SDcols=ests]
PSARm <- PSAR[,lapply(.SD,mean),by=.(intervention,g_whoregion),.SDcols=ests]
PSAAm <- PSAA[,lapply(.SD,mean),by=.(intervention,acat),.SDcols=ests]
PSACAm <- PSA[,lapply(.SD,mean),by=.(intervention,iso3,acat),.SDcols=ests]
PSACm <- PSA[,lapply(.SD,mean),by=.(intervention,iso3),.SDcols=ests]
## sds
PSAGs <- PSAG[,lapply(.SD,sd),by=.(intervention),.SDcols=ests]
PSARs <- PSAR[,lapply(.SD,sd),by=.(intervention,g_whoregion),.SDcols=ests]
PSAAs <- PSAA[,lapply(.SD,sd),by=.(intervention,acat),.SDcols=ests]
## IQRs
PSAGh <- PSAG[,lapply(.SD,uq),by=.(intervention),.SDcols=ests]
PSAAh <- PSAA[,lapply(.SD,uq),by=.(intervention,acat),.SDcols=ests]
PSAGl <- PSAG[,lapply(.SD,lq),by=.(intervention),.SDcols=ests]
PSAAl <- PSAA[,lapply(.SD,lq),by=.(intervention,acat),.SDcols=ests]
PSACAl <- PSA[,lapply(.SD,lq),by=.(intervention,iso3,acat),.SDcols=ests]
PSACl <- PSA[,lapply(.SD,lq),by=.(intervention,iso3),.SDcols=ests]
PSACAh <- PSA[,lapply(.SD,uq),by=.(intervention,iso3,acat),.SDcols=ests]
PSACh <- PSA[,lapply(.SD,uq),by=.(intervention,iso3),.SDcols=ests]

## == gathering
intl <- c("No intervention","Under 5 & HIV+ve","Under 5 & HIV+ve & LTBI+","B-A","C-A")
varlv <- ests[c(10,9,7,1,8,5,6,2,3,4)]  #reorder
nint <- length(intl)
hhrpl <- PSA[repn==1 & acat=="[5,15)" & hiv==0 & intervention==unique(intervention)[2],sum(e.households)] #the non zero value

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

## --- uncertainty
RTH <- melt(PSAGh,id='intervention'); names(RTH)[3] <- 'hi'
RTL <- melt(PSAGl,id='intervention'); names(RTL)[3] <- 'lo'
RTU <- merge(RTL,RTH); rm(RTH,RTL)
RTU[,u:=ppb(lo,hi,sf=4)]
RTU[,c('lo','hi'):=NULL]

## age
RTah <- melt(PSAAh,id=c('intervention','acat')); names(RTah)[4] <- 'hi'
RTal <- melt(PSAAl,id=c('intervention','acat')); names(RTal)[4] <- 'lo'
RTaU <- merge(RTal,RTah); rm(RTah,RTal)
RTaU[,u:=ppb(lo,hi,sf=4)]
RTaU <- dcast(RTaU,intervention+variable ~ acat,value.var = 'u')
names(RTaU)[3:4] <- c('ay','ao')

## sds for all
RTD <- melt(PSAGs,id='intervention'); names(RTD)[3] <- 'gsdev'
RTrD <- melt(PSARs,id=c('intervention','g_whoregion'))
RTrD <- dcast(RTrD,intervention + variable ~ g_whoregion,value='value')
names(RTrD)[3:ncol(RTrD)] <- paste0('s.',names(RTrD)[3:ncol(RTrD)])
RTaD <- melt(PSAAs,id=c('intervention','acat'))
RTaD <- dcast(RTaD,intervention + variable ~ acat,value='value')
names(RTaD)[3:4] <- c('sy','so')

## ## NNx TODO redo
## x <- RT[intervention=='B-A',Global]
## (ba <- c(hhn=-x[1]/x[6],hcn=-x[2]/x[6],ptn=-x[3]/x[6],txn=-x[4]/x[6]))
## x <- RT[intervention=='C-A',Global]
## (ca <- c(hhn=-x[1]/x[6],hcn=-x[2]/x[6],ptn=-x[3]/x[6],txn=-x[4]/x[6]))


## == formatting & output
## country-level output for supplementary
## by age too
names(PSACAh)[4:ncol(PSACAh)] <- paste0(names(PSACAh)[4:ncol(PSACAh)],'.hi')
names(PSACAl)[4:ncol(PSACAh)] <- paste0(names(PSACAl)[4:ncol(PSACAh)],'.lo')
PSACAm <- merge(PSACAm,PSACAh);PSACAm <- merge(PSACAm,PSACAl); rm(PSACAl,PSACAh)
PSACAm <- PSACAm[,c(names(PSACAm)[1:3],sort(names(PSACAm)[4:ncol(PSACAm)])),with=FALSE]
PSACAm[,c(grep('house',names(PSACAm),value=TRUE)):=NULL]
fwrite(PSACAm,file='tables/country_age_output.csv')
## just by country
names(PSACh)[3:ncol(PSACh)] <- paste0(names(PSACh)[3:ncol(PSACh)],'.hi')
names(PSACl)[3:ncol(PSACh)] <- paste0(names(PSACl)[3:ncol(PSACh)],'.lo')
PSACm <- merge(PSACm,PSACh);PSACm <- merge(PSACm,PSACl); rm(PSACl,PSACh)
PSACm[,c(grep('house',names(PSACm),value=TRUE)):=NULL]
PSACm <- merge(PSACm,HHV[,.(iso3,household.visits=visits)],by='iso3',all.x = TRUE)
PSACm[!grepl('5',intervention),household.visits:=0]                    #no visits
PSACm <- PSACm[,c(names(PSACm)[1:2],sort(names(PSACm)[3:ncol(PSACm)])),with=FALSE]
fwrite(PSACm,file='tables/country_output.csv')


## regional etc
RTL <- merge(RT,RTD)
RTL <- merge(RTL,RTaD)
RTL <- merge(RTL,RTrD)
whoz <- !RTL$variable %in% c('e.ATTprev')
RTL <- RTL[whoz]
RTL[,Global:=paste0(pps(Global,sf=4),'\n(',pps(gsdev,sf=4),')')]; RTL[,gsdev:=NULL]
regz <- PSA[,unique(g_whoregion)]; regzsd <- paste0('s.',regz)
for(i in seq_along(regz)){
  a <- RTL[,regz[i],with=FALSE]; b <- RTL[,regzsd[i],with=FALSE]
  RTL[,c(regz[i]):=paste0(pps(a,sf=4),'\n(',pps(b,sf=4),')')]
  RTL[,c(regzsd[i]):=NULL]
}
RTL[,`[0,5)`:=paste0(pps(`[0,5)`,sf=4),'\n(',pps(sy,sf=4),')')]; RTL[,sy:=NULL]
RTL[,`[5,15)`:=paste0(pps(`[5,15)`,sf=4),'\n(',pps(so,sf=4),')')]; RTL[,so:=NULL]
intz <- c(unique(PSA$intervention),'B-A','C-A')
RTL$intervention <- factor(RTL$intervention,levels=intz,ordered=TRUE)
RTL <- RTL[order(as.integer(intervention),variable),]
myft <- regulartable(RTL)
myft <- merge_v(myft, j = c("intervention", "variable") )
myft

## saving
save(RTL,file='tables/RTL.Rdata')
read_docx() %>%
  body_add_par(value = "Global, regional and age outputs", style = "heading 1") %>%
  body_add_flextable(value = myft2) %>%
  body_end_section(continuous = FALSE, landscape = TRUE) %>% 
  print(target = "tables/RTL.docx") %>% 
    invisible()


## simpler output for main article
RTS <- merge(RT,RTU,by=c('intervention','variable'))
RTS <- merge(RTS,RTaU,by=c('intervention','variable'))
RTS[,Global:=paste0(pps(Global,sf=4),'\n',u)]
RTS[,`[0,5)`:=paste0(pps(`[0,5)`,sf=4),'\n',ay)]
RTS[,`[5,15)`:=paste0(pps(`[5,15)`,sf=4),'\n',ao)]
RTS[,c('u','ay','ao'):=NULL]
names(RTS)[2] <- 'quantity'
RTS <- melt(RTS,id.vars = c('intervention','quantity'))
RTS <- dcast(RTS,quantity + variable ~ intervention,value.var = 'value')
RTS$variable <- factor(RTS$variable,levels=c('[0,5)','[5,15)','Global'),ordered=TRUE)
RTS <- RTS[order(variable,quantity),]
RTS <- RTS[,.(quantity,variable,`No intervention`,`Under 5 & HIV+ve`,
              `Under 5 & HIV+ve & LTBI+`,`B-A`,`C-A`)]
RTS <- RTS[!quantity %in% c('e.LTBI','e.ATTprev','')]
myft2 <- regulartable(RTS)
myft2 <- merge_v(myft2, j = c("variable","quantity") )
myft2

## saving
save(RTS,file='tables/RTS.Rdata')
read_docx() %>%
  body_add_par(value = "Global and age outputs", style = "heading 1") %>%
  ## body_add_table(value = myft, style = "table_template" ) %>%
  body_add_flextable(value = myft2) %>%
  body_end_section(continuous = FALSE, landscape = TRUE) %>% 
  print(target = "tables/RTS.docx") %>% 
    invisible()

## ## graph expts
## names(RT)[2] <- 'quantity'
## RTM <- melt(RT,id.vars = c('intervention','quantity'))

## ggplot(data=RTM[!grepl('-A',intervention) &
##                 quantity %in% c('e.ATT','e.IPT','e.incidence','e.deaths','e.LE') &
##                 variable %in% c('AFR','EUR','WPR','SEA','WPO','AMR')],
##        aes(x=variable,y=value,fill=intervention)) +
##   geom_bar(stat='identity',position='dodge') + 
##   coord_flip() + facet_wrap(~quantity,scales = 'free_x')

## ## or col by region:
## ggplot(data=RTM[!grepl('-A',intervention) &
##                 quantity %in% c('e.ATT','e.IPT','e.incidence','e.deaths','e.LE') &
##                 variable %in% c('AFR','EUR','WPR','SEA','WPO','AMR')],
##        aes(x=intervention,y=value,fill=variable)) +
##   geom_bar(stat='identity',position='dodge') + 
##   coord_flip() + facet_wrap(~quantity,scales = 'free_x')


## all country output too

## TODO
## document assumptions esp CDR
## CY compare
## hh uncertainty
## CHECK CDR for incident cases?
## consider different TST for HIV+
## consider different LE for HIV
