

rm(list=ls())
library(ggplot2)
library(ggpubr)
library(data.table)
library(here)
library(scales)
library(broom)
library(sensitivity)
library(patchwork)
set.seed(123)

# main data:
load(file=here('PL.Rdata'))
load(file=here('YAL.Rdata'))
load(file=here('YAL0.Rdata'))


# time series
YAL[id==1]
YAL[,t:=year + (Q-1)/4]; YAL0[,t:=year + (Q-1)/4]
YAL[1:20,diff(t)]
YAL[1:20,.(year,Q,t)]
kindex <- sum(YAL[id==1,t<100]) #ref index
YAL[(kindex-5):(kindex+5),.(year,Q,t)]
YAL[kindex,.(year,Q,t)]

YAL[,`Relative notifications`:=notif/notif[kindex],by=id]
YAL[,`Relative incidence`:=inc/inc[kindex],by=id]
YAL[,`Relative prevalence`:=prev/prev[kindex],by=id]



L <- 0.4; H <- 1.8
GP1 <- ggplot(YAL,aes(t-99.75,`Relative notifications`)) +
  stat_summary(geom="ribbon",fun.data=median_hilow,alpha=.3,
               fill="lightblue")+
  stat_summary(geom="line",fun=mean,col="red") +
  ylab('Relative quarterly notifications') +
  xlab('Time in years') +
  geom_vline(xintercept = 0,lty=2,col='darkgrey')+
  geom_vline(xintercept = 5,lty=2,col='darkgrey')+
  annotate("text",x=4,y=H+0.2,label='Notification peak',
           fontface = 'italic')+
  annotate("text",x=10,y=L,label='Notification trough',
           fontface = 'italic')+
  theme_classic() + ggpubr::grids()



GP2 <- ggplot(YAL,aes(t-99.75,`Relative incidence`)) +
  stat_summary(geom="ribbon",fun.data=median_hilow,alpha=.3,
               fill="lightblue")+
  stat_summary(geom="line",fun=mean,col="red") +
  geom_vline(xintercept = 0,lty=2,col='darkgrey')+
  geom_vline(xintercept = 5,lty=2,col='darkgrey')+
  ylab('Relative quarterly incidence') +
  xlab('Time in years') +
  ylim(c(0,1.2)) +
  annotate("text",x=9,y=0.75,label='Incidence trough',
           fontface = 'italic')+
  theme_classic() + ggpubr::grids()



GP3 <- ggplot(YAL,aes(t-99.75,`Relative prevalence`)) +
  stat_summary(geom="ribbon",fun.data=median_hilow,alpha=.3,
               fill="lightblue")+
  stat_summary(geom="line",fun=mean,col="red") +
  geom_vline(xintercept = 0,lty=2,col='darkgrey')+
  geom_vline(xintercept = 5,lty=2,col='darkgrey')+
  ylab('Relative prevalence') +
  xlab('Time in years') +
  ylim(c(0,1.2)) +
  annotate("text",x=9,y=0.3,
           label='Prevalence trough',fontface = 'italic')+
  theme_classic() + ggpubr::grids()



GPB <- ggarrange(GP1,GP2,GP3,ncol=1,align = 'h',labels = c('A','B','C'))
ggsave(GPB,file=here('ModelFigure1.pdf'),w=5,h=10)


# peaks and troughs

YAL[,tmax:=year[which.max(notif)],by=id] 
YAL[,nmax:=max(notif),by=id]             
YAL[,imin:=min(inc),by=id]               
YAL[,pmin:=min(prev),by=id]              
YAL[,rnmax:=max(`Relative notifications`),by=id] 
YAL[,rnmin:=min(`Relative notifications`),by=id] 
YAL[,rimin:=min(`Relative incidence`),by=id]     
YAL[,rpmin:=min(`Relative prevalence`),by=id]     

# rebounds
tmp <- YAL[t>=105 & t<=107,.(grown=`Relative notifications`/`Relative notifications`[1],
                             growi=`Relative incidence`/`Relative incidence`[1],
                             growp=`Relative prevalence`/`Relative prevalence`[1], 
                             t),by=id]
tmp[,t:=t-105]
tmpm <- melt(tmp,id.vars = c('id','t'))
tmp <- dcast(tmpm[t<=0.25], id+variable ~ t,value.var = 'value')
tmp[,f:=2-`0.25`]
tmp[,r:=-log(f)/0.25]        #1-exp(-rt) growth
tmp[,fg:=(`0.25`-`0`)/0.25]  #linear growth


tmp[,Td:=log(2)/r]
tmp[r>0,.(rlm=mean(fg),rld=median(fg), #lin growth
          rem=mean(r),red=median(r),   #exp growth
          Tdm=mean(Td),Tdd=median(Td), #doubling times
          lo=quantile(Td,0.025),hi=quantile(Td,0.975)),by=variable]

tmp <- dcast(tmp,id~variable,value.var = 'Td') #doubling times

# data for corr analysis
timax <- unique(YAL[,.(id,tmax,nmax,imin,rnmax,rnmin,rimin,pmin,rpmin)])
timax <- merge(timax,PL,by='id')


BL <- YAL[year==99 & Q==4,.(id,PNR=prev/notif,PR,CDR=notif/inc,prev=prev,
                            prev0=prev,inc0=inc,note0=notif)]

BL[,.(mean(inc0),quantile(inc0,c(.25,.5,.75)))] 

# baseline epi 
timax <- merge(timax,BL,by='id')
# merge in rebounds
timax <- merge(timax,tmp,by='id') 
# cumulative data
cumdat <- merge(YAL0[,.(cnotes0=sum(notif),cinc0=sum(inc)),by=id],
                YAL[,.(cnotes=sum(notif),cinc=sum(inc)),by=id],
                by='id')
cumdat[,c('rcnotes','rcinc'):=.(cnotes/cnotes0-1,cinc/cinc0-1)]
timax <- merge(timax,cumdat[,.(id,rcnotes,rcinc)],by='id')

names(timax)

summary(timax[,.(tmax,nmax,imin,rnmax,rimin,rnmin,pmin)])
timax[,qplot(tmax)] 


summary(timax[,.(rcinc,rcnotes)])
timaxr <- timax
names(timaxr)

(tmp1 <- pcc(timaxr[,.(rracf,hivPop,
                       CDR,PR,PNR,sut,prev)],timaxr[,rcinc])$PCC)
names(tmp1) <- 'Cumulative incidence'
(tmp2 <- pcc(timaxr[,.(rracf,hivPop,
                       CDR,PR,PNR,sut,prev)],timaxr[,rcnotes])$PCC)
names(tmp2) <- 'Cumulative notifications'
(tmp3 <- pcc(timaxr[,.(rracf,hivPop,
                       CDR,PR,PNR,sut,prev)],timaxr[,rnmax])$PCC)
names(tmp3) <- 'Notification peak'
(tmp4 <- pcc(timaxr[,.(rracf,hivPop,
                       CDR,PR,PNR,sut,prev)],timaxr[,rimin])$PCC)
names(tmp4) <- 'Incidence trough'
(tmp5 <- pcc(timaxr[,.(rracf,hivPop,
                       CDR,PR,PNR,sut,prev)],timaxr[,rnmin])$PCC)
names(tmp5) <- 'Notification trough'

## rebound
(tmp6 <- pcc(timaxr[is.finite(grown),.(rracf,hivPop,
                                       CDR,PR,PNR,sut,prev)],
             timaxr[is.finite(grown),grown])$PCC)
names(tmp6) <- 'Notification rebound timescale'
(tmp7 <- pcc(timaxr[,.(rracf,hivPop,
                       CDR,PR,PNR,sut,prev)],timaxr[,growi])$PCC)
names(tmp7) <- 'Incidence rebound timescale'

## prevalence trough
(tmp8 <- pcc(timaxr[,.(rracf,hivPop,
                       CDR,PR,PNR,sut,prev)],timaxr[,rpmin])$PCC)
names(tmp8) <- 'Prevalence trough'



SA3 <- cbind(tmp2,tmp1,tmp3,tmp5,tmp4,tmp6,tmp7,tmp8)
SA3$qty <- rownames(SA3)
names(SA3)
SA3M <- reshape2::melt(SA3,id='qty')

qty2 <- c('Screening HR','HIV prevalence',
          'CDR','Proportion recent','P:N ratio',
          'scale-up timescale','TB prevalence')

qtykey <- data.table(qty=c('rracf','hivPop',
                           'CDR','PR','PNR','sut','prev'),
                     qty2=qty2)


SA3M <- merge(SA3M,qtykey,by='qty')
SA3M <- as.data.table(SA3M)

smax <- SA3M[,.(value=max(abs(value)),value2=mean(abs(value))),by=qty2]
smax[order(value)]
smax[order(value2)]

SA3M$qty2 <- factor(SA3M$qty2,
                    levels = smax[order(value),qty2], #max abs val
                    ordered = TRUE)

vrs <- c('Incidence rebound timescale',
         'Prevalence trough',
         'Incidence trough',
         'Notification rebound timescale',
         'Notification peak',
         'Cumulative incidence',
         'Cumulative notifications',
         'Notification trough'
)
vrs <- vrs[c(5,3,2,8,4,1,6,7)] # order vars

SA3M$variable <- factor(SA3M$variable,
                        levels = vrs,
                        ordered = TRUE)



SA3P <- ggplot(SA3M,aes(variable,qty2,fill=value)) +
  geom_tile() +
  scale_fill_gradient2(midpoint=0,
                       low = "#0000FF",high ="#FF0000") +
  xlab('') + ylab('') +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=45,
                                   vjust=1,hjust=0,
                                   face='italic')) +
  theme(legend.position = 'bottom') +
  scale_x_discrete(position = "top") +
  labs(fill='Partial rank correlation coefficient')+
  theme(plot.background = element_rect(fill = 'lightgrey'),
        plot.margin = margin(.11,2,.11,.11,"cm"))
SA3P 

ggsave(SA3P,file=here('ModelFigure2c.pdf'),w=9,h=7)
ggsave(SA3P,file=here('ModelFigure2c.png'),w=9,h=7)


#====baseline epi=========

DT <-YAL[t==99.75,]

# calculate epi==========
DT[, IRR:=(`incall[2]`/`N[2]`)  /(`incall[1]`/`N[1]`)]
DT[, CDR:=notif/inc]
DT[, `P:N Ratio`:=prev/notif]

DT <- DT[, .(prev, inc, PR, `P:N Ratio`, CDR, IRR)]

DT <-cbind(DT, PL[,.(hivPop, sut, rracf)])

setnames(DT, c('prev','inc', 'PR','hivPop', 'sut', 'rracf'), 
         c('TB prevalence', 'TB incidence', 'Proportion recent',
           ' HIV prevalence', 'Scale-up timescale', 'Screening HR'))

DTL <- melt(DT)


# histogram of baseline epi===
bepi <-ggplot(DTL, aes(value)) + 
  geom_histogram() + facet_wrap(~variable, scale= 'free') + 
  ylab('Density') +xlab(NULL) +theme_bw()

ggsave(bepi, file= 'baselineEpi.pdf')




#=====relationship between inc,prev,notif=====

syear=100
eyear=105

inc   <- YAL[t>=syear & t<=eyear,.(imin=min(`Relative incidence`)),by=id]     
prev  <- YAL[t>=syear & t<=eyear,.(pmin=min(`Relative prevalence`)),by=id]   
notif <- YAL[t>=syear & t<=eyear,.(nmax=max(`Relative notifications`)),by=id] 

DTT <-cbind(inc, prev, notif)
DTT <-DTT[, .(imin, pmin, nmax)]


a<-ggplot(DTT, aes(imin, nmax,  group=1))+geom_point(alpha=0.4, size=1) + ylab('Relative peak notifications')+
  ggtitle('a')+ xlab('Minimum relative incidence') +theme_classic()

b<-ggplot(DTT, aes(pmin,nmax,  group=1))+geom_point(alpha=0.4, size=1)+ ylab('Relative peak notifications')+ 
  xlab('Minimum relative prevalence')+ggtitle('b') +theme_classic()

c<-ggplot(DTT, aes(imin, pmin, group=1))+geom_point(alpha=0.4, size=1) + ylab('Minimum relative prevalence')+
  ggtitle('c')+ xlab('Minimum relative incidence') +theme_classic() 



a+b+c


ggsave('inc_notif_prev_scatter.pdf', height = 3, width=8)


#==========

tmp <- timax[, .(rnmax, rnmin, rimin, rracf)]


d<-ggplot(tmp, aes(rracf, rnmax,  group=1))+geom_point(alpha=0.4, size=1) + ylab('Peak notifications(Relative)')+
  ggtitle('a')+ xlab('Screening HR') +theme_classic()

e<-ggplot(tmp, aes(rracf,rnmin,  group=1))+geom_point(alpha=0.4, size=1)+ ylab('Trough notifications (Relative)')+ 
  xlab('Screening HR')+ggtitle('b') +theme_classic()

f<-ggplot(tmp, aes(rracf, rimin, group=1))+geom_point(alpha=0.4, size=1) + ylab('Trough incidence (Relative)')+
  ggtitle('c')+  xlab('Screening HR') +theme_classic() 



d+e+f


ggsave('acf_HR_and_epi.pdf', height = 3, width=8)


