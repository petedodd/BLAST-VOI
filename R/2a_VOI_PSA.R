## Adapting to look at N-patch simulations
library(here)
library(data.table)
library(ggplot2)
library(epiR)
library(glue)

## === trying different data for
fn <- here("VOI/utils/pops.Rdata")
if(!file.exists(fn)){
  zones <- sf::read_sf(here("data/blantyre_7zone_update/blantyre_7zone_update.shp"))
  zones <- zones[order(zones$zone), ] # order by zone number
  pops <- zones$population # populations
  save(pops,file=fn)
} else {
  load(fn)
}

## === utils
source(here("VOI/utils/benefit.R")) #for benefit()

## for tables etc
gr <- function(m,l,h){
  m <- formatC(m,format="f",digits = 2)
  l <- formatC(l,format="f",digits = 2)
  h <- formatC(h,format="f",digits = 2)
  glue("{m} ({l} to {h})")
}


## ====================================================
## ** task 3 **: Analyse strategies & VOI with emulator


## these are from task 2

## patch veresion
np <- c(5,
        10,
        15,
        20,
        25,
        30,
        35,
        40,
        45,
        50)
RL <- PL <- list()
for(n in np){
  load(glue("~/Dropbox/Holocron/tmp/RP2/R{n}.Rdata"))
  load(glue("~/Dropbox/Holocron/tmp/RP2/P{n}.Rdata"))
  RL[[n]] <- R
  PL[[n]] <- P
}
R <- rbindlist(RL)
P <- rbindlist(PL)

R[,.N,by=.(iter,numpatches)]
P[,.N,by=.(iter,numpatches)]

all(R[,.N,by=.(iter,numpatches)]==P[,.N,by=.(iter,numpatches)])

RP <- R

## ## TODO BUG?
## RP <- merge(R[,.(slp,iter,numpatches)],
##             P[,.(real.prev.mean,iter,numpatches)],by=c("iter","numpatches")) #TODO do not understand
## ## RP <- merge(R,P,allow.cartesian = TRUE)


## ## ggplot(RP,aes(bet,slp))+geom_point()

## ## M <- lm(data=RP,slp ~ prev + ari0 + pf + alph + bet)
## ## qqnorm(M$residuals)

## ## RP[,pf:=NULL] #NOTE temp

## NOTE new conversion of slopes to have harmonized screenrate
## screenrate = 0.5 * screenrate0 / (num_patches/10)
## to get slope under screenrate0 -> x (num_patches/10) / 0.5 = num_patches/5
RP[,slp:=slp * numpatches / 5]


## ## trying to inspect slopes directly
## ## RP <- P #BUG
## RPS0 <- RP[,lapply(.SD,mean),.SDcols=4:ncol(RP),by=numpatches]
## RPS0

## ggplot(RPS0,aes(numpatches,slp)) +
##   geom_line() + geom_point() +
##   expand_limits(y=c(0,NA))

## ## NOTE we were saturating communities (staying too long)

## mean(RPS0$prev.CV/RPS0$real.prev.CV)
## plot(RPS0$prev.CV,RPS0$real.prev.CV); abline(a=0,b=1/2.3,col=2)

## RPS1 <- RP[,lapply(.SD,mean),.SDcols=4:ncol(RP),by=.(numpatches,variable)]
## RPS1

## ## TODO flat by zone? could up
## ggplot(RPS1,aes(variable,slp,group=numpatches))+
##   geom_line()+
##   facet_wrap(~numpatches) #,scales="free"

## compare with random, or population weighted order
## need to convert slopes from simulations into coverages

## slopes in simulations are: TBdeaths ~ t with screenrate = 1e4 per month
## to convert slopes in data to TBdeaths per person screened,
## need to multiply: popslope = simslope / 10^4



## difference in benefits (DbenefitFromSlps in utils/benefit.R)

## tests
true_slopes <- RP[iter==1 & numpatches==10,slp] / 1e4
pops <- rep(500e3/10,10)
DbenefitFromSlps(true_slopes,pops, #use pops as ranking
                  pops,0.5*sum(pops),verbose = TRUE,separate = TRUE)


## loop at half coverage:
B <- RP[,{
  DbenefitFromSlps(pmax(0,slp)/1e4,
                   ## rep(1,numpatches),
                   1:numpatches,
                   rep(500e3/numpatches,numpatches),
                   250e3,
                   separate = TRUE) #using com id as rank
},by=.(iter,numpatches)]

## ## slow way...
## B <- list()
## for(np in RP[,unique(numpatches)]){
##   cat("np=",np,"...\n")
##   for(k in 1:max(RP$iter)){
##     slopes <- RP[iter==k & numpatches==np,pmax(0,slp)/1e4]
##     ans <- DbenefitFromSlps(slopes,
##                             1:np,
##                             rep(500e3/np,np),
##                             250e3,
##                             verbose=FALSE,
##                             separate = TRUE) #using com id as rank
##     if(any(is.na(ans))) cat("...iter=",k,"\n")
##     tmp <- data.table(iter=k,numpatches=np)
##     tmp[,c("bnft_opt","bnft_sub"):=ans]
##     B[[paste(np,k)]] <- tmp
##   }
## }
## B <- rbindlist(B)
## Blong <- B

## check OK
## TEST <- merge(B,Blong,by=c("iter","numpatches"))
## TEST[,all(bnft_opt.x==bnft_opt.y)]
## TEST[,all(bnft_sub.x==bnft_sub.y)]


## np <- 50
## (test <- RP[iter==3 & numpatches==np][,pmax(0,slp)/1e4])
## pops <- rep(500e3/np,np)
## DbenefitFromSlps(test,
##                  1:np,
##                  pops,0.5*sum(pops),verbose = TRUE,separate = TRUE)


PB <- merge(P,B, by = c("iter","numpatches"))

save(PB,file="~/Dropbox/Holocron/tmp/RP2/PB.Rdata")

## RB <- merge(RP,B,by = c("iter","numpatches"))
## RB[is.na(bnft_opt),.(numpatches,slp)]

## PB[,c("bnft_opt",  "bnft_sub"):=NULL]
## RB[,c("bnft_opt",  "bnft_sub"):=NULL]

## TODO BUG still issue with 35?!
ggplot(PB,aes(numpatches,bnft_opt,group=numpatches)) + geom_boxplot()
## ggplot(PB,aes(numpatches,bnft_opt/bnft_sub,group=numpatches)) + geom_boxplot()


pairs(PB[,.(real.prev.mean, real.prev.CV, ari0, pf, alph, bet, RRbetaSD,  bnft_opt)],horInd = 8)

BM <- melt(PB[,.(real.prev.mean, real.prev.CV, ari0, pf, alph, bet, logRRbetaSD=log(RRbetaSD),
                 DB=bnft_opt,numpatches)],
           id=c("DB","numpatches")) #absolute

ggplot(BM[abs(DB)<1e3],aes(value,DB))+
  geom_point(shape=1)+
  facet_grid(numpatches~variable,scales="free") +
  ylab("TB deaths averted (optimal)")+
  theme_linedraw()

## ggsave(here("plots/VOI2_psa_total_multi.png"),w=10,h=10)


BM <- melt(PB[,.(real.Fprev.mean, real.Fprev.CV,real.foi.mean, real.foi.CV, alph,
                 DB=bnft_opt,numpatches)],
           id=c("DB","numpatches")) #absolute

ggplot(BM[abs(DB)<1e3],aes(value,DB))+
  geom_point(shape=1)+
  facet_grid(numpatches~variable,scales="free") +
  ylab("TB deaths averted (optimal)")+
  theme_linedraw()

ggsave(here("plots/VOI2_psa_total_multi2.png"),w=10,h=10)


BR <- PB[,.(real.prev.mean, real.prev.CV, ari0, pf, alph, bet,  logRRbetaSD=log(RRbetaSD),
            DB=bnft_opt-bnft_sub,numpatches)]
BM <- melt(BR,id=c("DB","numpatches")) #diff

ggplot(BM[abs(DB)<1e3],aes(value,DB))+
  geom_point(shape=1)+
  facet_grid(numpatches~variable,scales="free") +
  ylab("Difference in TB deaths averted (optimal-random)")+
  theme_linedraw()

## ggsave(here("plots/VOI2_EVPI.png"),w=10,h=7)



BM <- melt(PB[,.(real.Fprev.mean, real.Fprev.CV,real.foi.mean, real.foi.CV, alph,
                 DB=bnft_opt-bnft_sub,numpatches)],
           id=c("DB","numpatches")) #absolute

ggplot(BM[abs(DB)<1e3],aes(value,DB))+
  geom_point(shape=1)+
  facet_grid(numpatches~variable,scales="free") +
  ylab("Difference in TB deaths averted (optimal-random)")+
  theme_linedraw()

ggsave(here("plots/VOI2_EVPI2.png"),w=10,h=7)



lm(data=BR,DB ~.)

BRO <- PB[,.(real.Fprev.mean, real.Fprev.CV,real.foi.mean, real.foi.CV, alph,
            numpatches,
            DB=bnft_opt-bnft_sub)]

tab <- as.data.table(epi.prcc(BRO))
tab[,txt:=gr(est,lower,upper)]
tab

tabr <- tab[rev(order(abs(est)))][1:5][,.(variable=var,PRCC=txt)]
tabr
fwrite(tabr,file=here("plots/VOI2_psa_prccN.csv"))

BRO2 <- BRO[numpatches==50]
BRO2[,numpatches:=NULL]

tab <- as.data.table(epi.prcc(BRO2))
tab[,txt:=gr(est,lower,upper)]
tab

tabr <- tab[rev(order(abs(est)))][1:5][,.(variable=var,PRCC=txt)]
tabr
fwrite(tabr,file=here("plots/VOI2_psa_prccNR.csv"))


## unique(gsub("[[:punct:]]+|[[:digit:]]+", "", BLASTtbmod::get_cols))
## grep("TB_deaths",BLASTtbmod::get_cols,value=TRUE)



## Simulate across parameters exploring variation in:
##    Prevalence
##    Zone heterogeneity in prevalence
##    Transmission
##    Zone transmission heterogeneity & mixing
##    Progression
##    Zone progression heterogeneity
## Analyse & emulate slopes as function of inputs
## VOI from emulated slopes:
##    Perfect knowledge of relative prevalence vs transmission patterns,
##    within some strategy space

## Use the data with actual slopes to compute VOI
