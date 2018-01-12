## TODO
## still to write

load('data/D.Rdata')                    #load parent data

## reshape
DL <- melt(D,id.vars = c("iso3","country","g_whoregion","LAT",paste0("a",1:15),
                         "hivprop","artprop","cdr04","cdr04ab","cdr514","cdr514ab"))


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

## --- merge HH predictions
## TODO

