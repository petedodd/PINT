## This file prepares the final data files needed in the modelling analysis in file 6
## i.e. in order to build the parent data.table for the PSA
## NB the first statement wipes all data so be careful when applied!

rm(list=ls())

## ================ initial data pooling ===========
## NB to save file space, the publicly available datasets in this section have not been included in the repo:
## the relevant variables are bundled in 'data/D.Rdata'; skip down to the corresponding load statement
## or just use the overall output 'data/DL.Rdata'

## --- WHO data from
## http://www.who.int/tb/country/data/download/en/
## <2018-01-11 Thu>

## - notifications (start)
N <- fread('/Users/pjd/Documents/WHO_TBreports/data2017/TB/TB_notifications_2018-01-11.csv')

## below introduces an analysis of the proportions pulmonary
## relevant variables
nmz <- grep("new_.{2}_.+",names(N),value=TRUE)  #keep right patterns
nmz <- nmz[!grepl("unk",nmz)] #drop unknown
nmz <- nmz[!grepl("fu|mu",nmz)]       #more unknowns dropped
nmz <- nmz[!grepl('plus|514|04|014',nmz)] #drop children & not-child catchall
nmz <- c('iso3','year',nmz)

## reduce to relevant data
NP <- N[,nmz,with=FALSE]
tmp <- NP[,.(total=rowSums(.SD)),.SDcols=3:ncol(NP),by=.(iso3,year)]
NP <- merge(NP,tmp,by=c('iso3','year'))
NP <- NP[!is.na(total)]
NP <- NP[total>200]                     #remove small-TB countries 

## reshape
NP[,total:=NULL]
NP <- melt(NP,id.vars = c('iso3','year'))

## re-parse variable as separate characteristics
NP[,tbtype:=ifelse(grepl("ep",variable),'EP','Pulm')]
NP[,sex:=ifelse(grepl("f",variable),'F','M')]
NP[,age:=gsub("[a-z]|_","",variable)]
NP[,age:=gsub("(\\d{2})(\\d*)","\\1-\\2",age,perl=TRUE)]

## aggregate over years
NPR <- NP[,.(value=sum(value)),by=.(iso3,tbtype,sex,age)]
NPR[,total:=sum(value),by=.(iso3,age,sex)]
NPR[,proportion:=value/total,by=.(iso3,tbtype,sex,age)]
NPR[,bad:=any(!is.finite(proportion)),by=iso3]
NPR <- NPR[bad==FALSE]                  #drop those with things amiss

## plotting
abslab <- function(x,...) {             #scale labeler
  abs(x)
}

## for all countries
GP <- ggplot(data=NPR,aes(x=age,y=proportion,fill=tbtype)) +
  coord_flip() +
  geom_bar(data=NPR[sex=='M',],stat='identity',aes(y=proportion))+
  geom_bar(data=NPR[sex=='F',],stat='identity',aes(y=-proportion))+
  ylab('New notifications') + xlab('Age group')+
  geom_hline(yintercept=0,col=1,lty=2) +
  annotate(geom='text',x='65-',y=0.5,label='M')+
  annotate(geom='text',x='65-',y=-0.5,label='F')+
  scale_y_continuous(labels = abslab) +
  ggtitle("New notifications EP vs pulmonary by age and sex") +
  theme(title = element_text(size=22)) +
  facet_wrap(~iso3)
ggsave(GP,file='graphs/5EPvP.pdf',w=12,h=12)
ggsave(GP,file='graphs/5EPvP.png',w=12,h=12)


## regional averages
NPR <- merge(NPR,unique(N[,.(iso3,g_whoregion)]),by='iso3',all.x=TRUE,all.y=FALSE)
NPRreg <- NPR[,.(proportion=mean(proportion)),by=.(g_whoregion,tbtype,sex,age)]

GP <- ggplot(data=NPRreg,aes(x=age,y=proportion,fill=tbtype)) +
  coord_flip() +
  geom_bar(data=NPRreg[sex=='M',],stat='identity',aes(y=proportion))+
  geom_bar(data=NPRreg[sex=='F',],stat='identity',aes(y=-proportion))+
  ylab('New notifications') + xlab('Age group')+
  geom_hline(yintercept=0,col=1,lty=2) +
  annotate(geom='text',x='65-',y=0.5,label='M')+
  annotate(geom='text',x='65-',y=-0.5,label='F')+
  scale_y_continuous(labels = abslab) +
  ggtitle("New notifications EP vs pulmonary by age and sex") +
  facet_wrap(~g_whoregion)
ggsave(GP,file='5EPvPreg.pdf',w=7,h=7)
ggsave(GP,file='5EPvPreg.png',w=7,h=7)


## use these regional averages to fill out to all countries
NPR2 <- merge(unique(N[,.(iso3,g_whoregion)]),
              NPRreg,by='g_whoregion',allow.cartesian = TRUE)#reg averges for all iso3
NPR2 <- NPR2[!(iso3 %in% NPR[,unique(iso3)])]                #restrict to those w/o data
NPR <- rbind(NPR[,.(iso3,sex,age,tbtype,proportion)],
             NPR2[,.(iso3,sex,age,tbtype,proportion)]) #add in w/o data iso3s as means
setkey(NPR,iso3)
save(NPR,file='data/NPR.Rdata')
## end of pulmonary analysis

N <- N[year==2016]                      #restrict to most recent for onward use
## - notifications (end)

## estimations
E <- fread('/Users/pjd/Documents/WHO_TBreports/data2017/TB/TB_burden_countries_2018-01-11.csv')
E <- E[year==2016]                      #restrict
## merge
D <- merge(N,E,by=c('country','iso2','iso3','iso_numeric','g_whoregion','year'))
W <- fread('data/WBIL.csv')

## --- LAT
load('/Users/pjd/Documents/WHO_TBreports/LAT.Rdata')
D <- merge(D,LAT[,c('iso3','LAT')],by='iso3',all.x=TRUE)
## deal with NAs as regional averages
latna <- which(is.na(D$LAT))
latnareg <- D$g_whoregion[latna]
mnlat <- rep(0,length(latnareg))
for(i in 1:length(mnlat))               #wasteful
  mnlat[i] <- D[g_whoregion==latnareg[i],mean(LAT,na.rm = TRUE)]
D[latna,LAT:=mnlat]


## --- WHO/UNICEF BCG
B <- fread('data/BCG_11_1_2018.csv',skip=1)
B <- B[,c(2,5:41),with=FALSE]
names(B) <- c('iso3',paste0('a',1:(ncol(B)-1)))
B <- B[,1:16,with=FALSE]

## NA handling
nacnt <- rowSums(is.na(as.matrix(B[,2:16,with=FALSE])))
comp <- nacnt>0 & nacnt<15
stna <- is.na(B$a1)
enna <- is.na(B$a15)
B[comp & enna]
B[comp & stna]

## for NA at the end, fill with last
enners <- which(comp & enna)
for(i in enners){
  tmp <- B[i,c(15-nacnt[i]+1),with=FALSE]
  rng <- (15-nacnt[i]+1):ncol(B)
  B[i,(rng):=tmp]
}
B[enners]

## for NA at start...
begs <- which(comp & stna)
for(i in begs){
  tmp <- B[i,c(nacnt[i]+2),with=FALSE]
  rng <- 2:(nacnt[i]+1)
  B[i,(rng):=tmp]
}
B[begs]

## merge
D <- merge(D,B,by='iso3',all.x = TRUE)

rmn <- function(x) round(mean(x,na.rm=TRUE))

K <- D[,lapply(.SD,rmn),.SDcols=(ncol(D)-15+1):ncol(D),by=g_whoregion]
setkey(K,g_whoregion)
nna <- c(is.na(D[,ncol(D),with=FALSE]))

for(i in nna){
  for(j in (ncol(D)-15+1):ncol(D)) D[i,(j):= K[ D[i,g_whoregion], (j-ncol(D)+ncol(K)), with=FALSE] ]
}
print(D[nna,(ncol(D)-15+1):ncol(D),with=FALSE],n=Inf)


## --- tidy notifications by age/sex
D[,totnotes:=c_newinc]
D[is.na(totnotes),totnotes:=0]

## males
D[,n_m_0_4:=newrel_m04]
D[,n_m_5_14:=newrel_m514]
D[,n_m_15_24:=newrel_m1524]
D[,n_m_25_34:=newrel_m2534]
D[,n_m_35_44:=newrel_m3544]
D[,n_m_45_54:=newrel_m4554]
D[,n_m_55_64:=newrel_m5564]
D[,n_m_65_Inf:=newrel_m65]

## females
D[,n_f_0_4:=newrel_f04]
D[,n_f_5_14:=newrel_f514]
D[,n_f_15_24:=newrel_f1524]
D[,n_f_25_34:=newrel_f2534]
D[,n_f_35_44:=newrel_f3544]
D[,n_f_45_54:=newrel_f4554]
D[,n_f_55_64:=newrel_f5564]
D[,n_f_65_Inf:=newrel_f65]


names(D)[(ncol(D)-16+1):ncol(D)]

## total from disagg
rnp <- D[,lapply(.SD,rmn),.SDcols=(ncol(D)-16+1):ncol(D),by=g_whoregion] #regional note pattern
setkey(rnp,g_whoregion)
thna <- rowSums(D[,(ncol(D)-16+1):ncol(D),with=FALSE])            #is there an NA?
tfd <- rowSums(D[,(ncol(D)-16+1):ncol(D),with=FALSE],na.rm=TRUE)  #totes w/o NA
tfd <- D$totnotes - tfd
tfd[tfd<0] <- 0

## - add in the remaining notifications following regional pattern
## build key
rtots <- rowSums(rnp[,2:ncol(rnp),with=FALSE]);names(rtots) <- rnp$g_whoregion
rnpat <- list()
for(reg in rnp$g_whoregion) rnpat[[reg]] <- rnp[reg,2:ncol(rnp),with=FALSE]/rtots[reg]

## loop through
for(j in (ncol(D)-16+1):ncol(D))        #set NAs to 0 here
  set(D,which(is.na(D[[j]])),j,0)

for(i in 1:nrow(D)){                    #allocate excess by regional pattern
  if(tfd[i]>0){
    addon <- round(tfd[i]*rnpat[[D[i,as.character(g_whoregion)]]]) #regional patterned addon
    ## print(addon)
    D[i,((ncol(D)-16+1):ncol(D)):=D[i,(ncol(D)-16+1):ncol(D),with=FALSE]+addon]
  }
}


## HIV/ART
D[,.(iso3,newrel_hivtest,newrel_hivpos,newrel_art)]
D[,.(iso3,newrel_hivpos/newrel_hivtest,newrel_art/newrel_hivpos)]
D[,hivprop:=newrel_hivpos/newrel_hivtest]
D[,artprop:=newrel_art/newrel_hivpos]
D[!is.finite(hivprop),artprop:=0]
D[!is.finite(artprop),artprop:=0]
D[!is.finite(hivprop),hivprop:=0]


## --- CDR calculations
## estimates
D[,.(iso3,e_inc_num_m014,e_inc_num_f014,n_m_0_4,n_m_5_14,n_f_0_4,n_f_5_14,
     e_inc_num_m014_hi-e_inc_num_m014_lo,e_inc_num_f014_hi-e_inc_num_f014_lo)]

## 04 CDR (split incidence evenly)
D[,cdr04:=(n_m_0_4+n_f_0_4)/((e_inc_num_m014 + e_inc_num_f014)/2)]
D[,cdr04ab:=(n_m_0_4+n_f_0_4)/((e_inc_num_m014 + e_inc_num_f014)/2)]
D[!is.finite(cdr04),cdr04:=0]
D[cdr04>1,cdr04:=1]
D[,cdr04ab:=((1-cdr04)/cdr04)/((e_inc_num_m014_hi-e_inc_num_m014_lo)/(3.92*e_inc_num_m014))^2-1]
D[!is.finite(cdr04ab) | cdr04ab<0, cdr04ab:=0] #NB CDR sampling needs to handle 0s


## 514 CDR (split incidence evenly)
D[,cdr514:=(n_m_5_14+n_f_5_14)/((e_inc_num_m014 + e_inc_num_f014)/2)]
D[,cdr514ab:=(n_m_5_14+n_f_5_14)/((e_inc_num_m014 + e_inc_num_f014)/2)]
D[!is.finite(cdr514),cdr514:=0]
D[cdr514>1,cdr514:=1]
D[,cdr514ab:=((1-cdr514)/cdr514)/((e_inc_num_m014_hi-e_inc_num_m014_lo)/(3.92*e_inc_num_m014))^2-1]
D[!is.finite(cdr514ab) | cdr514ab<0, cdr514ab:=0] 


## --- multiply by proportion pumonary
nmz <- grep('^n_',names(D),value=TRUE)  #notification columns
NPR <- dcast(NPR[tbtype=="Pulm"],iso3 ~ age + sex, value.var = 'proportion')
## harmonize names
nmzp <- names(NPR)[-1]
nmzp <- gsub("-","_",nmzp)
nmzp[grepl("M",nmzp)] <- paste0("pn_m_",nmzp[grepl("M",nmzp)])
nmzp[grepl("F",nmzp)] <- paste0("pn_f_",nmzp[grepl("F",nmzp)])
nmzp <- gsub("_M","",nmzp); nmzp <- gsub("_F","",nmzp)
nmzp <- gsub("65_","65_Inf",nmzp)
names(NPR)[2:ncol(NPR)] <- nmzp
D <- merge(D,NPR,by='iso3',all.x = TRUE,all.y=FALSE) #join in
D[,c(nmz):=lapply(.SD,as.numeric),.SDcols=nmz]       #make numeric from int
nmz <- nmz[!grepl('0_4',nmz)]; nmz <- nmz[!grepl('5_14',nmz)] #drop children
for(nm in nmz) # multiply by corresponding pulmonary factor
  D[,c(nm):=D[,nm,with=FALSE] * D[,paste0('p',nm),with=FALSE]]
D[,c(nmzp):=NULL]                       #ditch the pulmonary props


## --- restrict
D <- D[,.(iso3,country,g_whoregion,LAT,
          a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15, #BCG by age
          n_m_0_4,n_m_5_14,n_m_15_24,n_m_25_34,n_m_35_44,n_m_45_54,n_m_55_64,n_m_65_Inf,
          n_f_0_4,n_f_5_14,n_f_15_24,n_f_25_34,n_f_35_44,n_f_45_54,n_f_55_64,n_f_65_Inf,
          hivprop,artprop,cdr04,cdr04ab,cdr514,cdr514ab)]
D
save(D,file='data/D.Rdata')

## =========== merging and reshaping

rm(list=ls())
load('data/D.Rdata')                    #load parent data
load('data/U5.Rdata')                   #load HH size predictions
load('data/O5.Rdata')                   #load HH size predictions
load('data/LEA.Rdata')                  #load life-expectancy data

## reshape
DL <- melt(D,id.vars = c("iso3","country","g_whoregion","LAT",paste0("a",1:15),
                         "hivprop","artprop","cdr04","cdr04ab","cdr514","cdr514ab"))

## remap for consistency
U5$acat <- plyr::mapvalues(U5$acat,from=levels(U5$acat),
                           to=c(levels(U5$acat)[1:5],'[65,Inf)'))

O5$acat <- plyr::mapvalues(O5$acat,from=levels(O5$acat),
                           to=c(levels(O5$acat)[1:5],'[65,Inf)'))


## test <- head(DL$variable)               #for testing functions below

## extract sex from notification variables
getsex <- function(x){
  toupper(unlist(lapply(strsplit(as.character(x),"_"),function(x)x[[2]])))
}
## getsex(test)                            #test

## extract age category from notification variables
getacat <- function(x){
  bot <- unlist(lapply(strsplit(as.character(x),"_"),function(x)x[[3]]))
  top <- unlist(lapply(strsplit(as.character(x),"_"),function(x)x[[4]]))
  paste0('[',bot,',',top,')')
}
## getacat(test) #test

## tidy up
DL[,sex:=factor(getsex(variable))]      #add sex category
DL[,acat:=factor(getacat(variable))]    #add age category
DL <- DL[!acat %in%c('[0,4)','[5,14)')] #drop child notifications
DL$acat <- plyr::mapvalues(DL$acat,from=as.character(DL[,unique(acat)]),
                           to=levels(U5$acat))


## --- merge HH predictions
names(U5)[2:5] <- c('sex','acat','HHu5mu','HHu5logsd')
U5[,sex:=factor(c('M','F')[as.numeric(as.character(sex))])]
DL <- merge(DL,U5,by=c('iso3','acat','sex'),all.x=TRUE)

names(O5)[2:5] <- c('sex','acat','HHo5mu','HHo5logsd')
O5[,sex:=factor(c('M','F')[as.numeric(as.character(sex))])]
DL <- merge(DL,O5,by=c('iso3','acat','sex'),all.x=TRUE)

## U5[iso3=='ZWE']
## DL[iso3=='ZWE']

## regional average for NAs
wcols <- c('HHu5mu','HHu5logsd','HHo5mu','HHo5logsd')
DL[,(wcols):=lapply(wcols,function(x){
  x <- get(x)
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}),by=.(g_whoregion,acat,sex)]


## --- merge LEA
LEAW <- dcast(LEA[,.(iso3,agegp,LE)],iso3~agegp,value.var='LE')
DL <- merge(DL,LEAW,by='iso3',all.x=TRUE)

## fill in NAs with regional average
wcols <- names(DL)[(ncol(DL)-3):ncol(DL)]
DL[,(wcols):=lapply(wcols,function(x){
  x <- get(x)
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}),by=g_whoregion]

DL[iso3=='ABW']
DL[iso3=='ZWE']

save(DL,file='data/DL.Rdata')


## === making parent data frame for younger children and calculating HH contacts
load('data/DL.Rdata')

## simplifying for u5 only
DL[,LE:=(`[0,1)`+4*`[1,5)`)/5]
DL[,`[0,1)`:=NULL]; DL[,`[1,5)`:=NULL]; DL[,`[5,10)`:=NULL]; DL[,`[10,15)`:=NULL]; 

## country level version
DLC <- unique(DL[,.(iso3,g_whoregion,LAT,a1,a2,a3,a4,a5,cdr04,cdr04ab,LE)])
DLC

## NB this is value hh stuff ( acat adult, but we can aggregate these)
DLK <- DL[,.(iso3,acat,sex,value,HHu5mu,HHu5logsd)]      #will probably need to invlude HIV stuff here

## extend for calculating numbers mean numbers of HH contacts
nrep <- 1e3
DLKL <- DLK[rep(1:nrow(DLK),nrep),]
DLKL[,repn:=rep(1:nrep,each=nrow(DLK))]

## u5 contacts found
DLKL[,phh := rlnorm(n = nrow(DLKL), meanlog=HHu5mu, sdlog=HHu5logsd)]
DLKL[,u5hhc := value * phh]             #under 5 HH mean contacts
## DLKL[,summary(phh)]                     #check
## DLKL[,qplot(phh)]           #check

## aggregates
chhc <- DLKL[,.(u5hhc=sum(u5hhc),notes=sum(value),phh=mean(phh)),
             by=.(repn,iso3)] #country aggregates

## ## checking!
## CHK<- merge(chhc,AM[,.(ny=mean(n04_m),no=mean(n514_m)),by=iso3],all.x=FALSE,all.y=TRUE)
## CHK[,ny:=ny*notes]
## CHK
## CHK[,summary(u5hhc/ny)]
## qplot(data=CHK,y=u5hhc,x=ny) + geom_abline(intercept=0,slope=1,col=2)


ghhc <- DLKL[,.(u5hhc=sum(u5hhc),notes=sum(value),phh=mean(phh)),
             by=.(repn)] #global aggregates


summary(ghhc)
ghhc[,mean(u5hhc)*1e-6]                 #around 3 million u5 contacts
## sanity checks
chhc[iso3=='ZAF',summary(u5hhc)]/DLK[iso3=='ZAF',sum(value)] #check
chhc[iso3=='ZAF',qplot(u5hhc)]
chhc[iso3=='AFG',qplot(u5hhc)]

chhc <- chhc[,.(u5hhc=mean(u5hhc),u5hhc.sd=sd(u5hhc),
                u5hhc.l=mean(log(u5hhc)),u5hhc.sdl=sd(log(u5hhc))),by=iso3]


save(chhc,file='data/chhc.Rdata')

## merge to make parent data table for PSA
DLC <- merge(DLC,chhc,by='iso3')
save(DLC,file='data/DLC.Rdata')

## ====== same but for 5-15
load('data/DL.Rdata')

## simplifying for u5 only
DL[,LE:=(`[5,10)`+`[10,15)`)/2]
DL[,`[0,1)`:=NULL]; DL[,`[1,5)`:=NULL]; DL[,`[5,10)`:=NULL]; DL[,`[10,15)`:=NULL]; 

## country level version
DLO <- unique(DL[,.(iso3,g_whoregion,LAT,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,cdr514,cdr514ab,LE)])
DLO

## NB this is value hh stuff ( acat adult, but we can aggregate these)
DLK <- DL[,.(iso3,acat,sex,value,HHo5mu,HHo5logsd)]      #will probably need to invlude HIV stuff here

## extend for calculating numbers mean numbers of HH contacts
nrep <- 1e3
DLKL <- DLK[rep(1:nrow(DLK),nrep),]
DLKL[,repn:=rep(1:nrep,each=nrow(DLK))]

## u5 contacts found
DLKL[,phh := rlnorm(n = nrow(DLKL), meanlog=HHo5mu, sdlog=HHo5logsd)]
DLKL[,o5hhc := value * phh]             #under 5 HH mean contacts
## DLKL[,summary(phh)]                     #check
## DLKL[,qplot(phh)]           #check

## aggregates
ohhc <- DLKL[,.(o5hhc=sum(o5hhc),notes=sum(value),phh=mean(phh)),
             by=.(repn,iso3)] #country aggregates

ghhc <- DLKL[,.(o5hhc=sum(o5hhc),notes=sum(value),phh=mean(phh)),
             by=.(repn)] #global aggregates

summary(ghhc)
ghhc[,mean(o5hhc)*1e-6]                 #around 5 million u5 contacts
## sanity checks
ohhc[iso3=='ZAF',summary(o5hhc)]/DLK[iso3=='ZAF',sum(value)] #check
ohhc[iso3=='ZAF',qplot(o5hhc)]
ohhc[iso3=='AFG',qplot(o5hhc)]

## ## checking!
## CHK<- merge(ohhc,AM[,.(ny=mean(n04_m),no=mean(n514_m)),by=iso3],all.x=FALSE,all.y=TRUE)
## CHK[,no:=no*notes]
## CHK
## CHK[,summary(o5hhc/no)]
## qplot(data=CHK,y=o5hhc,x=no) + geom_abline(intercept=0,slope=1,col=2) 


ohhc <- ohhc[,.(o5hhc=mean(o5hhc),o5hhc.sd=sd(o5hhc),
                o5hhc.l=mean(log(o5hhc)),o5hhc.sdl=sd(log(o5hhc))),by=iso3]

save(ohhc,file='data/ohhc.Rdata')


DLO <- merge(DLO,ohhc,by='iso3')
save(DLO,file='data/DLO.Rdata')


## === household visits
load('data/DL.Rdata')
DL <- unique(DL[,.(iso3,g_whoregion,acat,sex,value)])
HHV <- DL[,.(visits=sum(value)),by=.(iso3,g_whoregion)]

save(HHV,file='data/HHV.Rdata')
