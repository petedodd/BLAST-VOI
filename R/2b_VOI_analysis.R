library(here)
library(data.table)
library(ggplot2)
library(ggpubr)
library(epiR)
library(patchwork)

## load(file="~/Dropbox/Holocron/tmp/RP/PB.Rdata") #from Task3N
load(file="~/Dropbox/Holocron/tmp/RP2/PB.Rdata") #from Task3N aggregated version


## ============ individual figures

## ------------ fig 1a

## data:
BRO <- PB[,.(real.Fprev.mean, real.Fprev.CV,real.foi.mean, real.foi.CV, alph,
             numpatches,
             DB=bnft_opt-bnft_sub)]


tab <- as.data.table(epi.prcc(BRO))
tabr <- tab[rev(order(abs(est)))][1:5]
tabr$var <- factor(tabr$var,levels=rev(tabr$var),ordered=TRUE)
tabr[,var2:=fcase(var=="numpatches","Number of patches",
                  var=="real.Fprev.mean","TB prevalence mean",
                  var=="real.Fprev.CV","TB prevalence CV",
                  var=="real.foi.CV","FOI CV",
                  var=="real.foi.mean","FOI mean")]
tabr$var2 <- factor(tabr$var2,levels=rev(tabr$var2),ordered=TRUE)

## plot:
fig1a <- ggplot(tabr,aes(var2,y=est,ymin=ifelse(est>0,est,lower),ymax=ifelse(est>0,upper,est)))+
  geom_bar(stat="identity")+
  geom_errorbar(width=0,col=2,linewidth=1.5)+
  xlab("") + ylab("Partial rank correlation coefficient with VOI")+
  coord_flip()+
  theme_classic()+
  grids()

## fig1a  + theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust = 2))
## fig1a <- fig1a  + theme(axis.text.y = element_blank()) #maybe?
fig1a

## ------------ fig 1b
BN <- PB[,## !is.na(bnft_opt+bnft_sub), #BUG mode
         .(DB=mean(bnft_opt-bnft_sub),DB.sd=sd(bnft_opt-bnft_sub)),
         by=numpatches]

## ## BUG
## BN <- RP[,.(DB=mean(slp)),by=.(numpatches,iter)]
## BN <- BN[,.(DB=mean(DB),DB.sd=sd(DB)),by=numpatches]

fig1b <- ggplot(BN,aes(numpatches,DB))+
  geom_line()+geom_point()+
  theme_classic()+
  expand_limits(y=c(0,NA))+
  grids()+
  xlab("Number of zones")+
  ylab("Value of information")

fig1b


## ------------ fig 1c
## smooth
NB <- 10
PBR <- PB[numpatches>5 & real.Fprev.mean<0.02]
brks <- seq(from=0,to=2000,by=250)
PBR[,v1c:=cut(1e5*real.Fprev.mean,include.lowest = TRUE,ordered_result=TRUE,breaks=brks,dig.lab = 4)]
## PB[,v2c:=cut(bet,include.lowest = TRUE,ordered_result=TRUE,breaks=NB)]
PBR[,unique(v1c)]

XY <- PBR[,
         .(DB=mean(bnft_opt-bnft_sub)),
         by=.(v1c,numpatches)]

## TODO horrible - needs rethink
## plot
fig1c <- ggplot(XY,aes(v1c,numpatches,fill=DB))+
  geom_tile()+
  theme_minimal()+
  xlab("Prevalence per 100,000")+
  ylab("Number of patches")+
  scale_fill_viridis_c(option = "plasma")+ #plas,mag
  theme(legend.position = "top",
        legend.title = element_blank(),
        axis.text.x = element_text(angle=45,hjust = 1))

fig1c


## fig1c <- ggplot(PB,aes(alph,bet,col=bnft_opt-bnft_sub))+
##   geom_point()+
##   scale_y_log10()+
##   theme_minimal()+
##   theme(legend.position="none")


## ============ combined figure
## lhs <- ggarrange(plotlist=list(fig1a,fig1b),nrow=2,ncol=1,labels=LETTERS[1:2],align="h")
## rhs <- ggarrange(plotlist=list(fig1c),nrow=1,ncol=1,labels=LETTERS[3])

## Fig1 <- ggarrange(plotlist=list(lhs,rhs),nrow=1,ncol=2,widths=c(1,1.1))
## Fig1

## TODO align left axes etc etc
((fig1a/fig1b) | fig1c) + plot_layout(ncol=2,widths = c(1,1.5)) + plot_annotation(tag_levels = 'A')

ggsave(here("plots/VOI2_psa_fig.png"),w=10,h=6)

## TODO
## check units 1b
## variable names in figure 1a
## cuts in 1c//
## larger more aposite runs
