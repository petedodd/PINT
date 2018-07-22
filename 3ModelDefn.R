## This file defines the logic of the decision tree model and loads the relevant
## packages (see next para of comments). The functions the model neeeds to run
## are defined in file 4. This file can be sourced, although various test
## visualisations of the tree will print to console.
## To check the tree in a web browser use:
## plotter(kexp, varz=c('name','check'), edgelabel = TRUE)
## once constructed. 'check' is set to 1 for all leaves to check the associated
## function for the tree. Change this to one of the other variables computed for
## more general checking.

## Requirements: data.tree and HEdtree packages.
## install.packages('data.tree')
## should install data.tree.
## HEdtree can be installed by first installing devtools as above,
## and then calling:
## devtools::install_github('petedodd/HEdtree',dependencies=FALSE,build_vignettes=TRUE)

library(HEdtree)
library(data.tree)

## ==== tree for HHCT child
kexp <- Node$new('Child HH contact, age a, HIV/ART status')
atb <- kexp$AddChild('Prevalent TB disease')
natb <- kexp$AddChild('No prevalent TB disease') #misses tests

## kexp$Set(prevalent=0)
atb$prevalent <- 1
natb$prevalent <- 0

## not co-prevalent
ltbi <- natb$AddChild('LTBI (IGRA+ or TST+)')
noltbi <- natb$AddChild('no LTBI or TB disease')

## IPT
ltbiPTy <- ltbi$AddChild('LTBI, gets PT')
ltbiPTn <- ltbi$AddChild('LTBI, no PT')
noltbiPTy <- noltbi$AddChild('no LTBI or TB disease, gets PT') #FPs?
noltbiPTn <- noltbi$AddChild('no LTBI or TB disease, no PT')

## disease 
ltbiPTyT <- ltbiPTy$AddChild('TB disease <1y')
ltbiPTyN <- ltbiPTy$AddChild('no TB disease <1y')
ltbiPTnT <- ltbiPTn$AddChild('TB disease <1y')
ltbiPTnN <- ltbiPTn$AddChild('no TB disease <1y')
noltbiPTyT <- noltbiPTy$AddChild('TB disease <1y') 
noltbiPTyN <- noltbiPTy$AddChild('no TB disease <1y') 
noltbiPTnT <- noltbiPTn$AddChild('TB disease <1y')
noltbiPTnN <- noltbiPTn$AddChild('no TB disease <1y')

## ltbiPTyN$LE <- ltbiPTnN$LE <- noltbiPTyN$LE <-  noltbiPTnN <- 'LE'

## print(noltbi,'p','LE')

## parameters/functions here
atb$p <- 'coprev'                    #prevalence
natb$p <- '1-coprev'
ltbi$p <- 'ltbi.prev'                   #LTBI
noltbi$p <- '1-ltbi.prev'
ltbiPTy$p <- 'PTcov.P'                   #IPT
ltbiPTn$p <- '1-PTcov.P'
noltbiPTy$p <- 'PTcov.N'                 #redefine here
noltbiPTn$p <- '1-PTcov.N'
ltbiPTyT$p <- 'progn.LP.PTp'           #disease
ltbiPTyN$p <- '1-progn.LP.PTp'
ltbiPTnT$p <- 'progn.LP.PTn'
ltbiPTnN$p <- '1-progn.LP.PTn'
noltbiPTyT$p <- 'progn.LN.PTp'
noltbiPTyN$p <- '1-progn.LN.PTp'
noltbiPTnT$p <- 'progn.LN.PTn'
noltbiPTnN$p <- '1-progn.LN.PTn'

## print(kexp,'p')

## ===== disease outcomes
outcomes <- Node$new('TB disease outcomes')
octxY <- outcomes$AddChild('TB treatment')
octxN <- outcomes$AddChild('no TB treatment')
octxYd <- octxY$AddChild('dies')
octxYs <- octxY$AddChild('survives')
octxNd <- octxN$AddChild('dies')
octxNs <- octxN$AddChild('survives')

## --- parameters here
outcomes$p <- 1;  #needed as linking on to others
octxY$p <- 'CDR';
octxN$p <- '1-CDR';
octxYd$p <- 'CFRtxY';
octxNd$p <- 'CFRtxN';
octxYs$p <- '1-CFRtxY';
octxNs$p <- '1-CFRtxN';

## print(outcomes,'p')

## deaths
outcomes$Set(death=0)
octxYd$death <- octxNd$death <- 1

## tx
outcomes$Set(treatments=0)
octxY$treatments <- 1

## incidence
outcomes$Set(incidence=0)
outcomes$incidence <- 1

## life-expectancy
outcomes$Set(LE=0)
octxYs$LE <- octxNs$LE <- 'LE'

outcomes$Set(check=0)
outcomes$Set(check=1,filterFun=isLeaf)

print(outcomes,'death','treatments','incidence','LE')



## ===== combined tree
noltbiPTyT$AddChildNode(Clone(outcomes))
noltbiPTnT$AddChildNode(Clone(outcomes))
ltbiPTyT$AddChildNode(Clone(outcomes))
ltbiPTnT$AddChildNode(Clone(outcomes))
atb$AddChildNode(Clone(outcomes))


print(kexp,'p')
## plotter(kexp, varz=c('name'), edgelabel = FALSE)
## plotter(kexp, varz=c('name','check'), edgelabel = TRUE)

## fill in all blanks in indicator variables appropriately
## fill in incidence
kexp$Set(incidence=0,filterFun=function(x)is.null(x$incidence))
kexp$Set(incidence=0,filterFun=function(x)is.na(x$incidence))
kexp$`Prevalent TB disease`$`TB disease outcomes`$incidence <- 0 #coprev doesn't count

## fill in LE
kexp$Set(LE=0,filterFun=function(x)is.null(x$LE))
kexp$Set(LE=0,filterFun=function(x)is.na(x$LE))
kexp$Set(LE='LE',filterFun=function(x) x$name=='no TB disease <1y' & x$isLeaf)
## plotter(kexp,varz = c('name','incidence'))
print(kexp,'incidence','LE')

## fill in death
kexp$Set(death=0,filterFun=function(x)is.null(x$death))
kexp$Set(death=0,filterFun=function(x)is.na(x$death))
## fill in treatment
kexp$Set(treatments=0,filterFun=function(x)is.null(x$treatments))
kexp$Set(treatments=0,filterFun=function(x)is.na(x$treatments))

print(kexp,'death','treatments')

## fill in prevalent
kexp$Set(prevalent=0,filterFun=function(x)is.null(x$prevalent))
kexp$Set(prevalent=0,filterFun=function(x)is.na(x$prevalent))
## check
kexp$Set(check=0)
kexp$Set(check=1,filterFun=isLeaf)

## PT courses
kexp$Set(IPT=0)
kexp$Set(IPT=1,filterFun=function(x) grepl("gets PT",x$name) )
## print(kexp,'IPT')
## plotter(kexp,varz = c('name','IPT'))

print(kexp,'prevalent','check')
## plotter(kexp, varz=c('name'), edgelabel = FALSE) # inspect tree by plot

## LTBI
kexp$Set(LTBI=0)
kexp$Set(LTBI=1,filterFun=function(x) grepl("TST+",x$name) )
## print(kexp,'LTBI')

## ATTprev (to distinguish from treatments from detected incidence)
kexp$Set(ATTprev=0)
kexp$`Prevalent TB disease`$`TB disease outcomes`$`TB treatment`$Set(ATTprev=1,filterFun=function(x) !x$isLeaf)
## print(kexp,'ATTprev')

## ATTinc (to distinguish from treatments from detected incidence)
kexp$Set(ATTinc=0)
kexp$`No prevalent TB disease`$Set(ATTinc=1,filterFun=function(x) grepl('treatment',x$name) & !grepl('no',x$name))
## print(kexp,'ATTinc','ATTprev')

## cases both inc and prevalence
kexp$Set(cases=0)
kexp$Set(cases=1,filterFun=function(x) grepl('outcomes',x$name))
## print(kexp,'cases','prevalent','incidence')
## plotter(kexp,varz=c('name','cases','prevalent','incidence'))

## deaths by inc prev
kexp$Set(deathinc=0)
kexp$Set(deathprev=0)
kexp$`No prevalent TB disease`$Set(deathinc=1,filterFun=function(x) grepl('dies',x$name))
kexp$`Prevalent TB disease`$Set(deathprev=1,filterFun=function(x) grepl('dies',x$name))
## print(kexp,'deathinc','deathprev')

## different CDR for prevalent
kexp$`Prevalent TB disease`$`TB disease outcomes`$`TB treatment`$p <- "CDRp"
kexp$`Prevalent TB disease`$`TB disease outcomes`$`no TB treatment`$p <- "1-CDRp"

## ## checking
## print(kexp,'treatments','ATTprev')
## plotter(kexp,varz = c('name','ATTprev','treatments'))
plotter(kexp,varz = c('name','ATT'))
