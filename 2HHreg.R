## TODO
## edit, neaten, check; change wdir to repo

## to look at doing regression in stan etc.
setwd('/Users/pjd/Documents/Rwork/hhstructure/hhint')
rm(list=ls())
load('AM.Rdata')

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
## do differently
NP <- 1



AMW <- dcast(AM[!acat %in% c('[0,5)','[5,15)')],
             iso3 ~ acat+sex,value.var=c('y0m','y0v'))

## corrections
AMW[iso3=='NCL',iso3:='NIC']
AMW[iso3=='CYM',iso3:='KGZ']
AMW[iso3=='GUM',iso3:='GTM']


## =================================================
## load up covariate data
load('XX.Rdata')
load('XX2.Rdata')
load('XP.Rdata')
load('XP2.Rdata')
load('~/Documents/WHO_TBreports/data/isodict.Rdata')

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
## EPS <- matrix(1e4,nrow=nrow(EPS),ncol=ncol(EPS))
ZZ <- ZZ[,grepl('m',colnames(ZZ))]      #means
## XX <- matrix(1,nrow=nrow(ZZ),ncol=1)
K <- ncol(ZZ)
P <- ncol(XX2)
N <- nrow(ZZ)

## check
(plt <- ggplot(data=AM[iso3=='ZWE' & !acat %in% c('[0,5)','[5,15)')],
              aes(x=acat,y=log(n04_m),group=iso3,col=iso3)) +
  geom_point() + geom_line() +
  geom_errorbar(aes(ymin=y0m-2*sqrt(y0v),ymax=y0m+2*sqrt(y0v)),width=0) + 
  facet_grid(sex~g_whoregion))


## MLE/MAP attempt for B
xi <- t(XX2) %*% XX2
## xi <- xi %x% diag(K)
xi <- diag(K) %x% xi                    #TODO THIS WAS WRONG!!
yy <- t(XX2) %*% ZZ  %*% diag(K)
vy <- c(yy)
bmap <- solve(xi,vy)
bmap <- matrix(bmap,nrow=P,ncol=K)
ypred <- XX2 %*% bmap
colnames(ypred) <- colnames(ZZ)
ypdf <- as.data.table(ypred)
ypdf$iso3 <- AMW[,iso3]
ypdf <- melt(ypdf,id='iso3')
ypdf[,sex:=getsex(variable)]
ypdf[,acat:=getacat(variable)]
ypdf[,n04_m:=exp(value)]

plt + geom_point(data=ypdf[iso3=='ZWE'],col=3,shape=2) 


## =================================================
## fitting

library(mvregerr)

## Gibbs sampling
## out <- mvregerrGS(ZZ,EPS,XX, nchain = 4*50,
##                 init=list(Psi=diag(ncol(ZZ)),nu=5),
##                 every = 20,record = c('Y','beta'))

BP <- matrix(5,nrow=P,ncol=K)           #5 us default
## smaller Psi? think about motivation -- seems sensitive to this

out <- mvregerrGS(ZZ,EPS,XX2, nchain = 5*4*50,
                init=list(Psi=diag(ncol(ZZ))*5e-2,nu=5,B=BP),
                every = 20,record = c('Y','beta'),XP=XX2)

## TODO prediction by who region

tmp <- out$Y
tmpp <- out$YP                           #PREDICTION!!
## for(i in setdiff(1:length(tmp),seq(from=50,to=length(tmp),by=4))) tmp[[i]] <- NULL #thin

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


## xx <- c(AMW[iso3=='ZWE'])
## xx$iso3 <- NULL
## xx <- unlist(xx)
## plot(xx[1:12])
## plot(xx[seq(from=1,by=2,len=6)])
## plot(out$Y[[50]][68,seq(from=1,by=2,len=6)])


## plot(out$beta[[50]][1,])

## TODO add acat and value alignment

## getacat(head(y0df$variable))
## getsex(head(y0df$variable))

## add in acat & exp of value
y0df[,acat:=factor(getacat(variable))]
y0df[,sex:=factor(getsex(variable))]
y0df[,n04_m:=exp(value)]

y0pdf[,acat:=factor(getacat(variable))]
y0pdf[,sex:=factor(getsex(variable))]
y0pdf[,n04_m:=exp(value)]

## merge in g_whoregion
## tmp <- AM[,.(iso3,g_whoregion)]
## setkey(tmp,iso3)
## tmp <- tmp[,.SD[1,],by=iso3]
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


ggplot(y1p,aes(datp,n04_m,col=acat,shape=sex)) +
  geom_point() +
  geom_abline(intercept = 0,slope=1,lty=2) +
  coord_fixed(ratio=1) +
  xlim(rg) + ylim(rg) +
  facet_wrap(~g_whoregion) +
  xlab('Data') + ylab('Prediction')

ggsave('ng/Predictions.pdf')

## which countries have the prediction outliers
## y1p[abs(1-datp/n04_m)>.5,]
(bad <- y1p[abs(n04_m-datp)>.6,as.character(unique(iso3))])


ggplot(data=y1p[iso3 %in% bad],
       aes(x=acat,y=datp,group=iso3,col=iso3)) +
  facet_grid(iso3~sex) +
  geom_point(aes(acat,n04_m),shape=2)  +
  geom_point() + geom_line() +
  theme(legend.position="none",axis.text.x = element_text(angle=90))
ggsave('ng/BadPredictions.pdf')

ggplot(data=y1p[g_whoregion %in% c('AFR','EMR','EUR')],
       aes(x=acat,y=datp,group=iso3,col=iso3)) +
  facet_grid(g_whoregion~sex) +
  geom_point(data=y1p[iso3 %in% bad],aes(x=acat,y=datp),shape=2)  +
  geom_label(data=y1p[iso3 %in% bad],aes(x=acat,y=datp,label=iso3))+
  geom_point() + geom_line() +
  theme(legend.position="none",axis.text.x = element_text(angle=90))
ggsave('ng/BadPredictionsOutlie.pdf')

## ----- residuals
y1 <- y0df[,.(n04_m=mean(n04_m)),by=.(iso3,acat,sex,g_whoregion)]
y1 <- merge(y1,AM[!acat %in% c('[0,5)','[5,15)'),.(datp=n04_m,acat,sex,iso3)],
            by=c('iso3','acat','sex'),all.x=TRUE)
rg <- c(0,2.5)

ggplot(y1,aes(datp,n04_m,col=acat,shape=sex)) +
  geom_point() +
  geom_abline(intercept = 0,slope=1,lty=2) +
  coord_fixed(ratio=1) +
  xlim(rg) + ylim(rg) +
  facet_wrap(~g_whoregion) +
  xlab('Data') + ylab('Fitted value')

ggsave('ng/Residuals.pdf')

## ========================
## ---- particular region

AMR <- AM[!acat %in% c('[0,5)','[5,15)') & g_whoregion=='EMR']
yr <- y0pdf[g_whoregion=='EMR']
## AMR <- AM[!acat %in% c('[0,5)','[5,15)') & iso3=='ZWE']
## yr <- y0df[iso3=='ZWE' & rep==1]
## yr <- y0df[iso3=='ZWE']


plt <- ggplot(data=AMR,
              aes(x=acat,y=log(n04_m),group=iso3,col=iso3)) +
  facet_grid(iso3~sex) +
  ## geom_point(data=yr,shape=2)  +
  geom_violin(data=yr,aes(x=acat,y=log(n04_m),group=factor(paste(acat,sex,iso3))))  +
  geom_point() + geom_line() +
  theme(legend.position="none",axis.text.x = element_text(angle=90))
  ## guides(color=FALSE)
  ## theme(legend.position="none")

plt

## ggsave(plt,filename = 'ngWX_AFR.pdf',width=7,height=30)
## ggsave(plt,filename = 'ngnoX_AFR.pdf',width=7,height=30)


## TODO
## prediction for all cns!
## longer chain / 
## other outputs?
## TODO over 5 too?

out <- mvregerrGS(ZZ,EPS,XX2, nchain = 5*4*50,
                init=list(Psi=diag(ncol(ZZ))*5e-2,nu=5,B=BP),
                every = 5*20,record = c('Y','beta'),XP=XP2)



tmpp <- out$YP                           #PREDICTION!!

y0 <- do.call('rbind',tmpp)
colnames(y0) <- colnames(ZZ)
y0df <- as.data.table(y0)
y0df[,iso3:=XP$iso3]
y0df[,rep:=rep <- rep(1:length(tmpp),each=nrow(tmpp[[1]]))]
y0df[,]

y0df <- melt(y0df,id=c('iso3','rep'))



