## DHS data analysis for survey-weighted mean (CI) number of children in household by:
##   country & year, & age category & sex of cohabitee
## Raw DHS input data for this analysis is not provided due to license restrictions.
##   (...it is also tens of gigabytes...)
## However it is freely available on application to: https://dhsprogram.com/data/
## the outputs of this analysis (which are included in this repo) are:
##   AM.Rdata (used in onward analyses)
##   graphs: India_04s.pdf, India_514s.pdf, AMplot1.pdf, AMplot2.pdf

rm(list=ls())
library(data.table)
library(haven)

dta2dt <- function(x) as.data.table(read_dta(x)) #utility function

## First lot of county data in STATA format
dhsfilelist <- c('Bangladesh/BDHR61DT/BDHR61FL.DTA',
                 'Brazil/BRHR31DT/BRHR31FL.DTA',
                 'Cambodia/KHHR61DT/KHHR61FL.DTA',
                 'Congo/CDHR50dt/cdhr50fl.dta',
                 'Ethiopia/ETHR61DT/ETHR61FL.DTA',
                 'India/IAPR52DT/IAPR52FL.dta',
                 'Indonesia/IDHR51dt/IDHR51FL.DTA',
                 'Kenya/KEHR52DT/KEHR52FL.DTA',
                 'Mozambique/MZHR61DT/MZHR61FL.DTA',
                 'Nigeria/NGHR61DT/NGHR61FL.DTA',
                 'Philippines/PHHR52DT/PHHR52FL.DTA',
                 'SouthAfrica/ZAHR31DT/ZAHR31FL.DTA',
                 'Tanzania/TZHR6ADT/TZHR6AFL.DTA',
                 'Uganda/UGHR6ADT/UGHR6AFL.DTA',
                 'Vietnam/VNHR52DT/VNHR52FL.DTA',
                 'Zimbabwe/ZWHR62DT/ZWHR62FL.DTA'
                 )

## see MAP file for mapping
varz <- c('hhid',  #id
          'hv000', #country code and phase
          'hv005', #weight
          'hv007', #year
          'hv012', #n usual hh member
          'hv013', #de facto hh members
          'hv014', #n usual 0-5 hh members
          'hv021', #cluster id
          'hv022', #strata var
          'hv025' #rural/urban
)
## 'hv104', #sex 'hv105'  #age

## make a list and save
all <- list()
for(fn in dhsfilelist){
  cat('working on ',fn,'...\n')
  P <- dta2dt(fn)                       #read in
  keep <- c(varz,grep('hv104',names(P),value=TRUE),grep('hv105',names(P),value=TRUE))
  if(P[,all(is.na(hv021))] & 'sh021' %in% names(P)) P[,hv021:=sh021] #for India
  P <- P[,keep,with=FALSE]                #restrict to variables
  all[[fn]] <- P
}

save(file='all.Rdata',all)


## Add second lot data from Courtney (in SPSS SAV format)
sav2dt <- function(x) as.data.table(read_spss(x)) #utility function

setwd('SAV')
dhsfilelist <- list.files(pattern='sav')

## make a list and save
all2 <- list()
for(fn in dhsfilelist){
  cat('working on ',fn,'...\n')
  P <- sav2dt(fn)                       #read in
  names(P) <- tolower(names(P))
  keep <- c(varz,grep('hv104',names(P),value=TRUE),grep('hv105',names(P),value=TRUE))
  ## if(P[,all(is.na(hv021))] & 'sh021' %in% names(P)) P[,hv021:=sh021] #for India
  P <- P[,keep,with=FALSE]                #restrict to variables
  all2[[fn]] <- P
}

setwd('..')

## tmp files
save(file='all2.Rdata',all2)
load('all.Rdata')
load('all2.Rdata')


## combine and tidy
ALL <- list()
cnt <- 0
for(i in 1:(length(all)+length(all2))){
  if(i<=length(all)){
    LV <- all[[i]]
  } else { LV <- all2[[i-length(all)]] }
  if(! LV[1,hv000] %in% names(ALL)){    #if not done already
    cat('working on',LV[1,hv000],'\n')
    cnt <- cnt+1
    tmp <- melt(data=LV,id.vars = varz)
    tmp[,N2:=length(variable),by=hhid]    #no records for hh
    tmp <- tmp[,age:=value[(N2[1]/2+1):N2[1]],by=hhid] #add age from bottom half
    tmp[,n:=1:N2[1],by=hhid]                           #record number in hh
    tmp[,keep:= (n<=N2[1]/2) & (!is.na(value) & !is.na(age)),by=hhid] #top half  & no-na
    ## restrict
    tmp <- tmp[keep==TRUE,.(srv=hv000,hhid=hhid,year=hv007,age=age,sex=value,
                            w=hv005,clid=hv021,strat=hv022,rurb=hv025,
                            hhnu=hv012,hhnd=hv013,hhk5u=hv014)]
    tmp <- tmp[order(hhid),]              #sort
    tmp[,iso2:=substr(srv,start=1,stop=2)] #ISO2 code
    if(tmp[1,nchar(year)]==2) tmp[,year:=paste0('19',year)] #correct year if necessary
    ALL[[tmp[1,srv]]] <- tmp                                #add
  }
}
print(cnt)
print(length(ALL))

save(ALL,file='ALLallRAW.Rdata')        #tmp file
load('ALLallRAW.Rdata')

## oops BRA is chr! & other issues to correct:
for(i in 1:length(ALL)){
  ALL[[i]][,srv:=as.factor(srv)]
  ALL[[i]][,hhid:=as.factor(hhid)]
  ALL[[i]][,year:=as.integer(year)]
  ALL[[i]][,age:=as.integer(age)]
  ALL[[i]][,sex:=as.factor(sex)]
  ALL[[i]][,w:=as.numeric(w)]
  ALL[[i]][,clid:=as.factor(clid)]
  ALL[[i]][,strat:=as.factor(strat)]
  ALL[[i]][,rurb:=as.factor(rurb)]
  ALL[[i]][,hhnu:=as.integer(hhnu)]
  ALL[[i]][,hhnd:=as.integer(hhnd)]
  ALL[[i]][,hhk5u:=as.integer(hhk5u)]
  ALL[[i]][,iso2:=as.factor(iso2)]
}


## collapse
ALL <- do.call('rbind',ALL)             #CD5 and CD6
ALL <- ALL[srv!='CD5',]                 #drop earlier survey
ALL[,length(unique(iso2))]              #68 surveys

## use Courtney's spreadsheet to correctly assign country ids and years
load('CY2.Rdata')
CY2 <- CY2[iso3!='PAK']
setkey(CY2,DHS)
ALL[,dhs:=substr(srv,start=1,stop=2)]
tmp <- CY2[as.character(ALL[,dhs])]
ALL[,iso3:=tmp$iso3]
ALL[,year:=tmp$year]
ALL[,dhs:=NULL]
load('isodict.Rdata')                   #ISO dictionary
ALL <- merge(ALL,ISO[,.(iso2,iso3,g_whoregion)],by='iso3') 
ALL[,length(unique(iso3))]              #68
ALL <- ALL[order(iso3,hhid),]
setkey(ALL,iso3)
save(ALL,file='ALLall2.Rdata')
load('ALLall2.Rdata')

## survey analysis accounting for sampling design
cnz <- ALL[,as.character(unique(iso3))]
library(survey)
options(survey.lonely.psu="adjust")

## add age categories
brks <- c(0,5,15,25,35,45,55,65,Inf)
ALL[,acat:=cut(age,breaks=brks,include.lowest = TRUE,right = FALSE)]

## add hh counts
ALL[,n04:=sum(age<5),by=.(iso3,hhid)]
ALL[,n514:=sum(age>4 & age<15),by=.(iso3,hhid)]

tmp <- ALL['IDN']
tmp[hhid==tmp[1,hhid]]              #check


ALL <- ALL[sex!='9',]                   #drop unknown sex
ALL$sex <- factor(ALL$sex)

## looped analysis for all countries, age groups, sex: WARNING takes a while!
## check: iso3 %in% c('ARM','AZE','BEN')
AM <- ALL[,{
  print(as.character(iso3[1]))
  print(as.character(acat[1]))
  print(as.character(sex[1]))
  options(survey.lonely.psu="adjust")
  strata <- strat
  if( any(is.na(strata)) ) strata <- rep(1,length(strata))
  DHSdesign <- svydesign(id = clid,strata=strata,weights = w/1000000,data=.SD)
  t1 <- svymean(~n04,DHSdesign)
  t2 <- svymean(~n514,DHSdesign)
  .SD[,.(n04_m=t1[1],
         n04_v=attributes(t1)$var[1],
         n514_m=t2[1],
         n514_v=attributes(t2)$var[1])
      ]
}, by=.(acat,sex,iso3,g_whoregion)]

save(AM,file='AM.Rdata')                #save out

## quick look
gp <- ggplot(AM,aes(x=acat,y=n04_m,col=iso3,group=iso3)) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=n04_m-1.96*sqrt(n04_v),ymax=n04_m+1.96*sqrt(n04_v)),width=0) + 
  facet_grid(g_whoregion~sex)

ggsave('graphs/0AMplot1.pdf',gp,height=14)                    #save out

gp <- ggplot(AM,aes(x=acat,y=n514_m,col=iso3,group=iso3)) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=n514_m-1.96*sqrt(n514_v),ymax=n514_m+1.96*sqrt(n514_v)),
                width=0) +
  facet_grid(g_whoregion~sex)

ggsave('graphs/0AMplot2.pdf',gp,height=14)                    #save out




## ===================================================

## load Indian DHS data set again
P <- dta2dt('IAPR52DT/IAPR52FL.dta')    # ~ 1.5GB

## see MAP file for mapping
varz <- c('hhid',
          'hv000',
          'hv102', #de jure household members, 1 if person and missing ow
          'hv005', #weight
          'hv011', #mum alive
          'hv012', #mum id
          'hv013', #dad alive
          'hv014', #dad id
          'hv021', #cluster id
          'hv022', #strata var
          'hv104', #region
          'hv025', #rural/urban
          'hv104', #sex
          'hv105', #age
          'hv226', #IAP
          'sh30'   #TB (india)
          )

## for india
P[,hv021:=sh021]

P <- P[,varz,with=FALSE]                #restrict to variables

## clean TB variable
P[,table(sh30)]
P[sh30==9,sh30:=0]
P[,TB:=factor(sh30)]
P[,sex:=c('Male','Female')[hv104]]      #recode sex
P[,region:=c('urban','rural')[hv025]]      #recode rural/urban

## add some more
P[,n:=sum(hv104>0),by=hhid]             #hh size
P[,n04:=sum(hv105<5),by=hhid]           #no kids <5 in hh
P[,n514:=sum(hv105>4 & hv105<15),by=hhid]
P[hhid==P[1,1,with=FALSE]]              #check

## patterns
## ------ 0-4 ------
ggplot(P[hv105>14],aes(x=hv105,y=n04,group=TB,col=TB)) +
  geom_smooth() +
  facet_grid(region ~ sex) +
  xlab('Age') + ylab('Number cohabiting children age 0-4') 
ggsave('graphs/0India_04s.pdf')

## ------ 5-14 ------
ggplot(P[hv105>14],aes(x=hv105,y=n514,group=TB,col=TB)) +
  geom_smooth() +
  facet_grid(region ~ sex) +
  xlab('Age') + ylab('No. cohabiting children age 5-14') 
ggsave('graphs/0India_514s.pdf')

