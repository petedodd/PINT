## TODO
## edit, neaten, check;
## O5s

rm(list=ls())
load('data/AM.Rdata')

## variance in log-transformed space:

## σ≈sinh(s/m)/2
## if
## (m−2s,m+2s)=m.(exp(2σ),exp(−2σ))

## mean exp(mu+sig^2/2)
## median exp(mu)
## var exp(2mu+sig^2)(exp(sig^2)-1) = m^2 (exp(sig^2)-1)
## var/m^2+1 = exp(sig^2)


NR <- length(unique(AM$iso3))
AS <- 2*length(unique(AM$acat))

AM[,y0v:=log(1+n04_v/n04_m^2)]                                    #log
AM[,y0m:=log(n04_m)]                                    #log
NP <- 1                                                 #initialize

AMW <- dcast(AM[!acat %in% c('[0,5)','[5,15)')],
             iso3 ~ acat+sex,value.var=c('y0m','y0v'))

## =================================================
## load up covariate data
load('data/XX.Rdata')
load('data/XX2.Rdata')
load('data/XP.Rdata')
load('data/XP2.Rdata')
load('data/isodict.Rdata')


## investigate transformation
summary(XP)
XP$EN.POP.DNST <- log(XP$EN.POP.DNST)
XP$NY.GDP.PCAP.KD <- log(XP$NY.GDP.PCAP.KD)
XP$SP.DYN.IMRT.IN <- log(XP$SP.DYN.IMRT.IN)
XP$SP.DYN.TFRT.IN <- log(XP$SP.DYN.TFRT.IN)
XP$SP.DYN.LE00.IN <- XP$SP.DYN.LE00.IN^2

XPM <- melt(XP,id='iso3')
ggplot(XPM,aes(value)) +
  facet_wrap(~variable,scales='free') + geom_histogram()

## transform:
for(nm in c('EN.POP.DNST','NY.GDP.PCAP.KD','SP.DYN.IMRT.IN','SP.DYN.TFRT.IN')){
  XX2[,nm] <- log(XX2[,nm])
  XP2[,nm] <- log(XP2[,nm])
}
XX2[,'SP.DYN.LE00.IN'] <- XX2[,'SP.DYN.LE00.IN']^2
XP2[,'SP.DYN.LE00.IN'] <- XP2[,'SP.DYN.LE00.IN']^2


## check covarates aligned with AMW
setdiff(XX[,as.character(iso3)],AMW[,as.character(iso3)])
setdiff(AMW[,as.character(iso3)],XX[,as.character(iso3)])

## order XXs to match AMW
der <- unlist(lapply(AMW[,as.character(iso3)],function(x)which(XX[,as.character(iso3)]==x)))

(bad <- which(XX[der,as.character(iso3)]!=AMW[,as.character(iso3)]))

XX <- XX[der,]
XX2 <- XX2[der,]


## =================================================
## utility functions

getacat <- function(x){
  x <- as.character(x)
  tmp <- strsplit(x,split = '_')
  unlist(lapply(tmp,function(X)X[2]))
}


getsex <- function(x){
  x <- as.character(x)
  tmp <- strsplit(x,split = '_')
  unlist(lapply(tmp,function(X)X[3]))
}

## =================================================
## final preparations

ZZ <- as.matrix(AMW[,.SD,.SDcols=names(AMW)[-1]])
EPS <- ZZ[,grepl('v',colnames(ZZ))]     #variances
ZZ <- ZZ[,grepl('m',colnames(ZZ))]      #means
K <- ncol(ZZ)
P <- ncol(XX2)
N <- nrow(ZZ)


## =================================================
## fitting

library(mvregerr)

BP <- matrix(5,nrow=P,ncol=K)           #5 us default

out <- mvregerrGS(ZZ,EPS,XX2, nchain = 5*4*50,
                init=list(Psi=diag(ncol(ZZ))*5e-2,nu=5,B=BP),
                every = 20,record = c('Y','beta'),XP=XX2)

tmp <- out$Y
tmpp <- out$YP                           #PREDICTION!!

y0 <- do.call('rbind',tmp)
colnames(y0) <- colnames(ZZ)
y0df <- as.data.frame(y0)
y0df$iso3 <- AMW[,as.character(iso3)]
y0df$rep <- rep(1:length(tmp),each=nrow(tmp[[1]]))
y0df <- as.data.table(y0df)
y0df <- melt(y0df,id=c('iso3','rep'))

y0p <- do.call('rbind',tmpp)
colnames(y0p) <- colnames(ZZ)
y0pdf <- as.data.frame(y0p)
y0pdf$iso3 <- AMW[,as.character(iso3)]
y0pdf$rep <- rep(1:length(tmp),each=nrow(tmp[[1]]))
y0pdf <- as.data.table(y0pdf)
y0pdf <- melt(y0pdf,id=c('iso3','rep'))


## add in acat & exp of value
y0df[,acat:=factor(getacat(variable))]
y0df[,sex:=factor(getsex(variable))]
y0df[,n04_m:=exp(value)]

y0pdf[,acat:=factor(getacat(variable))]
y0pdf[,sex:=factor(getsex(variable))]
y0pdf[,n04_m:=exp(value)]

## merge in g_whoregion
tmp <-  merge(y0df,ISO[,.(iso3,g_whoregion)],by='iso3',all.x=TRUE,all.y=FALSE)
y0df <- tmp

tmp <-  merge(y0pdf,ISO[,.(iso3,g_whoregion)],by='iso3',all.x=TRUE,all.y=FALSE)
y0pdf <- tmp



## ===== residuals/prediction errors
## ----- predictions
y1p <- y0pdf[,.(n04_m=mean(n04_m)),by=.(iso3,acat,sex,g_whoregion)]
y1p <- merge(y1p,AM[!acat %in% c('[0,5)','[5,15)'),.(datp=n04_m,acat,sex,iso3)],
            by=c('iso3','acat','sex'),all.x=TRUE)
rg <- c(0,2.5)


gp <- ggplot(y1p,aes(datp,n04_m,col=acat,shape=sex)) +
  geom_point() +
  geom_abline(intercept = 0,slope=1,lty=2) +
  coord_fixed(ratio=1) +
  xlim(rg) + ylim(rg) +
  facet_wrap(~g_whoregion) +
  xlab('Data') + ylab('Prediction')

ggsave('graphs/2Predictions.pdf',gp)
ggsave('graphs/2Predictions.png',gp)

## which countries have the prediction outliers
## y1p[abs(1-datp/n04_m)>.5,]
(bad <- y1p[abs(n04_m-datp)>.6,as.character(unique(iso3))])

gp <- ggplot(data=y1p[g_whoregion %in% c('AFR','EMR','EUR')],
       aes(x=acat,y=datp,group=iso3,col=iso3)) +
  facet_grid(g_whoregion~sex) +
  geom_point(data=y1p[iso3 %in% bad],aes(x=acat,y=datp),shape=2)  +
  geom_label(data=y1p[iso3 %in% bad],aes(x=acat,y=datp,label=iso3))+
  geom_point() + geom_line() +
  theme(legend.position="none",axis.text.x = element_text(angle=90))

ggsave('graphs/2BadPredictionsOutlie.pdf',gp)
ggsave('graphs/2BadPredictionsOutlie.png',gp)

## ----- residuals
y1 <- y0df[,.(n04_m=mean(n04_m)),by=.(iso3,acat,sex,g_whoregion)]
y1 <- merge(y1,AM[!acat %in% c('[0,5)','[5,15)'),.(datp=n04_m,acat,sex,iso3)],
            by=c('iso3','acat','sex'),all.x=TRUE)
rg <- c(0,2.5)

gp <- ggplot(y1,aes(datp,n04_m,col=acat,shape=sex)) +
  geom_point() +
  geom_abline(intercept = 0,slope=1,lty=2) +
  coord_fixed(ratio=1) +
  xlim(rg) + ylim(rg) +
  facet_wrap(~g_whoregion) +
  xlab('Data') + ylab('Fitted value')

ggsave('graphs/2Residuals.pdf',gp)
ggsave('graphs/2Residuals.png',gp)

## ========================
## prediction for all cns!

out <- mvregerrGS(ZZ,EPS,XX2, nchain = 5*4*50,
                init=list(Psi=diag(ncol(ZZ))*5e-2,nu=5,B=BP),
                every = 5*20,record = c('Y','beta'),XP=XP2)



tmpp <- out$YP                           #PREDICTION!!

y0 <- do.call('rbind',tmpp)
colnames(y0) <- colnames(ZZ)
y0df <- as.data.table(y0)
y0df[,iso3:=XP$iso3]
y0df[,rep:=rep <- rep(1:length(tmpp),each=nrow(tmpp[[1]]))]
## y0df[,]

y0df <- melt(y0df,id=c('iso3','rep'))

U5 <- y0df
U5[,sex:=factor(getsex(variable))]
U5[,acat:=factor(getacat(variable))]
U5[,variable:=NULL]
U5 <- U5[,.(value=mean(value),HHu5logsd=sd(value)),by=.(iso3,sex,acat)] ## mean
save(U5,file='data/U5.Rdata')


## check nothing ridiculous
load('data/isodict.Rdata')
U5 <- merge(U5,ISO[,.(iso3,g_whoregion)],by='iso3')

ggplot(data=U5[g_whoregion=='EUR'],
       aes(acat,exp(value),group=paste0(iso3,sex),col=iso3)) +
  geom_point() + geom_line()

ggplot(data=U5,
       aes(acat,exp(value),group=paste0(iso3,sex),col=iso3)) +
  geom_point() + geom_line() +
  facet_wrap(~g_whoregion) + theme(legend.position='none')

## TODO copy over o5s



