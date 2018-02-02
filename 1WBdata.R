## This file includes the data download and cleaning/preparation for:
## - world bank covariate data for household structure modelling
## - global health observatory life-table data

rm(list=ls())
## ===================== world bank covariate data
library(wbstats)

new_wb_cache <- wbcache()

IDS <- c("NY.GDP.PCAP.KD", "SP.DYN.LE00.IN", "SP.DYN.IMRT.IN", #GDP, LE, infant mort
         "SP.POP.TOTL","SP.POP.0014.TO",                       #total pop, <15 pop
         "SI.POV.GINI","SP.URB.TOTL.IN.ZS",                    #GINI & %urban
         "EN.POP.DNST","SP.DYN.TFRT.IN"                       #pop density,total fert
         )

wb_dat <- wb(indicator = IDS)
names(wb_dat)

wb_countries <- wbcountries() 
names(wb_countries)

wb_dat <- merge(wb_dat, y = wb_countries[c("iso2c", "region",'iso3c')], by = "iso2c", all.x = TRUE)

wb_dat <- subset(wb_dat, region != "Aggregates")

wb <- as.data.table(wb_dat)
wb[,iso3:=iso3c]
wb[,iso3c:=NULL]
setkey(wb,iso3)

wbmr <- wb[,.SD[date==max(date),],by=.(iso3,indicator)] #most recents
save(wbmr,file='data/wbmr.Rdata')

load('data/wbmr.Rdata')

load('data/CY2.Rdata')
wbcs <- merge(wb,CY2,by='iso3',all.x = TRUE)
wbcs <- wbcs[!is.na(DHS)]               #only the DHS countries
wbcs[,cy:=as.integer(date)]
wbcs[,dst:=abs(cy-year)]              #absolute difference from survey date
wbcs <- wbcs[,.SD[dst==min(dst),],by=iso3] #closest to surveys
wbcs[dst>0]
save(wbcs,file='data/wbcs.Rdata')


load('data/isodict.Rdata')
tmp <- ISO[,.(iso3,g_whoregion)]
mf <- model.frame(formula=~g_whoregion,data=tmp)
XX2 <- model.matrix(mf,data=tmp)


## cast for comparison
## for modelling
XX <- dcast(wbcs[,.(iso3,indicatorID,value)],iso3~indicatorID,value.var = 'value')
XX[,SI.POV.GINI:=NULL]
XX[,FRACK:=SP.POP.0014.TO/SP.POP.TOTL]
XX[,SP.POP.0014.TO:=NULL]
XX[,SP.POP.TOTL:=NULL]
head(XX)
dim(XX)
save(XX,file='data/XX.Rdata')

XX2 <- merge(ISO[,.(iso3,g_whoregion)],XX,by='iso3',all.x = FALSE,all.y=TRUE)
tmp <- XX2
tmp <- tmp[,2:ncol(tmp),with=FALSE]
mf <- model.frame(formula=~.,data=tmp)
tmp <- model.matrix(mf,data=tmp)
XX2 <- tmp
save(XX2,file='data/XX2.Rdata')


## for prediction
XP <- dcast(wbmr[,.(iso3,indicatorID,value)],iso3~indicatorID,value.var = 'value')
XP[,SI.POV.GINI:=NULL]
XP[,FRACK:=SP.POP.0014.TO/SP.POP.TOTL]
XP[,SP.POP.0014.TO:=NULL]
XP[,SP.POP.TOTL:=NULL]
XP
dim(XP)
XP <- XP[!is.na(EN.POP.DNST+NY.GDP.PCAP.KD+SP.DYN.IMRT.IN+SP.URB.TOTL.IN.ZS+FRACK+SP.DYN.LE00.IN+SP.DYN.TFRT.IN),]
save(XP,file='data/XP.Rdata')


XP2 <- merge(ISO[,.(iso3,g_whoregion)],XP,by='iso3',all.x = FALSE,all.y=TRUE)
tmp <- XP2
tmp <- tmp[,2:ncol(tmp),with=FALSE]
mf <- model.frame(formula=~.,data=tmp)
tmp <- model.matrix(mf,data=tmp)
XP2 <- tmp
save(XP2,file='data/XP2.Rdata')

## ============== life table data

library(rgho)
## nMx - age-specific death rate between ages x and x+n	LIFE_0000000029
## nqx - probability of dying between ages x and x+n	LIFE_0000000030
## lx - number of people left alive at age x	LIFE_0000000031
## ndx - number of people dying between ages x and x+n	LIFE_0000000032
## nLx - person-years lived between ages x and x+n	LIFE_0000000033
## Tx - person-years lived above age x	LIFE_0000000034
## ex - expectation of life at age x	LIFE_0000000035

## search_codes("life", dimension = "GHO")
## search_codes("LIFE_0000000035", dimension = "GHO")


result <- get_gho_data(dimension = "GHO",code = "LIFE_0000000035",
                         filter = list(
                           YEAR = "2015"
                         ))

## print(result)
## print(result, width = Inf)

res <- as.data.table(result)
res <- res[AGEGROUP %in% c("AGELT1","AGE1-4","AGE5-9","AGE10-14"),
           .(iso3=COUNTRY,AGEGROUP,REGION,sex=SEX,LE=Numeric)]

## same region codes
res[,g_whoregion:=REGION]
res[REGION=='SEAR',g_whoregion:='SEA']
res[,REGION:=NULL]

## new agegroup codes
res[,agegp:=AGEGROUP]
res[agegp=="AGELT1",agegp:='[0,1)']
res[agegp=="AGE1-4",agegp:='[1,5)']
res[agegp=="AGE5-9",agegp:='[5,10)']
res[agegp=="AGE10-14",agegp:='[10,15)']
res[,AGEGROUP:=NULL]
res$agegp <- factor(res$agegp,levels=c('[0,1)','[1,5)','[5,10)','[10,15)'),ordered=TRUE)

## ditch non-countries
res <- res[!is.na(iso3)]
res <- res[,.(LE=mean(LE)),by=.(iso3,g_whoregion,agegp)]
res <- res[order(iso3,agegp)]
setkey(res,iso3)

## check
res['AFG']
res['ZWE']
res['GBR']

## save out
LEA <- res
save(LEA,file='data/LEA.Rdata')


