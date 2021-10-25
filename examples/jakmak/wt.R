# wt.R
radian
library(ggplot2)
library(tidyverse)
library(data.table)
source("../../R/readadmb.R")
source("../../R/helper.R")
getwd()

df<-read_dat("wt.dat")
names(df)

# get data (fishery)     
df_fsh <- data.frame(year=df$fshry_yrs,df$fishery)
names(df_fsh) <- c("year",1:12)
# Anomalies on the data
maxage=10 # For display only...
pivot_longer(df_fsh,cols = 2:13, names_to = "age", values_to = "wt") %>% group_by(age) %>% 
   mutate(age=as.numeric(age), mnwt=mean(wt)) %>% ungroup() %>% filter(age<=maxage) %>% mutate(anom=wt/mnwt-1,Anomaly=ifelse(abs(anom)>.5,NA,anom) ) %>%
   ggplot(aes(y=year,x=age,fill=Anomaly,label=round(wt,2))) + geom_tile() + 
   scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
   geom_text(size=3) + ylab("Year") + xlab("Age") + ggtitle("Jack mackerel data") +
   scale_y_reverse() + theme_minimal(base_size=18)


   scale_fill_gradient2(high = scales::muted("blue"), low = scales::muted("red"), na.value = "white",
   breaks=c(-0.2,0.2), labels=format(c(-0.2,0.2)), limits = c(-.5,.5) ) + 

#--------------------------
# Read in results/estimates
#--------------------------
   # show some examples
  df_res <- fn_get_pred(file='wt',source="Jack mackerel")
  pivot_longer(df_res, cols=2:13,names_to = "age", values_to = "wt") %>% group_by(age,source) %>% 
   mutate(age=as.numeric(age), mnwt=mean(wt)) %>% ungroup() %>% filter(year>=1975,age>=2,age<=maxage) %>% 
   mutate(anom=wt/mnwt-1,Anomaly=ifelse(abs(anom)>.5,NA,anom) ) %>%
   ggplot(aes(y=year,x=age,fill=Anomaly,label=round(wt,2))) + geom_tile() + 
     scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
     geom_text(size=3) + ylab("Year") + xlab("Age") + 
     scale_y_reverse() + theme_minimal(base_size=18) + ggtitle("Jack mackerel model")
     if (length(unique(dfin$source))>1) p1 <- p1 + facet_grid(.~source)
 }
 head(df_res)
fn_plot_anoms(df_res)
dfin <- df_res
names(dfin)[3:dim(dfin)[2]-1] 
names(dfin)[2]
  names(df)
  fn_plot_anoms(df_res)
  tt<- read_admb("wt")  
  tt$fit$names 
  wt_hist <- data.frame(matrix(tt$fit$est[grep("wt_hist",tt$fit$names )],ncol=13,byrow=TRUE),Year=df_rep$yr,"est")
  names(wt_hist) <- c(3:15,"Year","type")
  wt_std <- data.frame(matrix(tt$fit$std[grep("wt_hist",tt$fit$names )],ncol=13,byrow=TRUE),Year=df_rep$yr,"std")
  names(wt_std) <- c(3:15,"Year","type")
  wt_cv <- data.frame(wt_std[,1:13]/ wt_hist[,1:13],Year=df_rep$yr,"CV")
  names(wt_cv) <- c(3:15,"Year","type")
  wt_res<- rbind(wt_hist,wt_std,wt_cv)

pivot_longer(wt_res,cols = 1:13, names_to = "age", values_to = "CV") %>% 
   mutate(age=as.numeric(age)) %>% filter(age<=maxage,type=="CV",Year>1970) %>% 
   ggplot(aes(y=Year,x=age,fill=CV,label=round(CV,2))) + geom_tile() + 
   scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
   geom_text(size=3) + ylab("Year") + xlab("Age") + ggtitle("Weight-at-age CV") +
   scale_y_reverse() + theme_minimal(base_size=16) #+ facet_grid(.~source)

  df_rep <- read_rep(paste0("wt.rep"))
  df_pred<- data.frame(year=df_rep$yr,df_rep$wt_pre)
  names(df_pred) <- c("year",3:15)

  df_nosrv<-read_rep("wt_nosrv.rep")
  df_nosrv<- data.frame(year=df_nosrv$yr,df_nosrv$wt_pre)
  names(df_nosrv) <- c("year",3:15)

pivot_longer(df_pred,cols = 2:14, names_to = "age", values_to = "wt") %>% group_by(age) %>% 
   mutate(age=as.numeric(age), mnwt=mean(wt)) %>% ungroup() %>% filter(age<=maxage) %>% mutate(anom=wt/mnwt-1,Anomaly=ifelse(abs(anom)>.5,NA,anom) ) %>%
   ggplot(aes(y=year,x=age,fill=Anomaly,label=round(wt,2))) + geom_tile() + 
   scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
   geom_text(size=3) + ylab("Year") + xlab("Age") + 
   scale_y_reverse() + theme_minimal(base_size=18) #+ facet_grid(.~source)

df_all <- rbind(data.frame(df_fsh,source="fishery"),data.frame(df_srv,source="index"),data.frame(df_pred,source="model"))
df_all <- rbind(data.frame(df_fsh,source="fishery data"),data.frame(df_srv,source="survey data"),data.frame(df_pred,source="model with survey"),data.frame(df_nosrv,source="model no survey"))
names(df_all) <- c("year",3:15,"source")
tail(df_all)

# Functions
   file=c("wt");source=c("model");i=1
fn_get_pred <- function(file=c("wt"),source=c("model")){
    df_tmp <- data.frame()
    for (i in 1:length(file)){
      dd     <- read_rep(paste0(file[i],".rep"))
      df_tmp <- rbind(df_tmp, data.frame(year=dd$yr,dd$wt_pre,source=source[i]))
    }
    names(df_tmp) <- c("year",3:15,"source")
    return(df_tmp)
  }
fn_plot_anoms <- function(df, maxage=10){
  pivot_longer(df,cols = 2:14, names_to = "age", values_to = "wt") %>% group_by(age,source) %>% 
   mutate(age=as.numeric(age), mnwt=mean(wt)) %>% ungroup() %>% filter(year>=1982,age<=maxage) %>% mutate(anom=wt/mnwt-1,Anomaly=ifelse(abs(anom)>.5,NA,anom) ) %>%
   ggplot(aes(y=year,x=age,fill=Anomaly,label=round(wt,2))) + geom_tile() + 
   scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
   geom_text(size=3) + ylab("Year") + xlab("Age") + 
   scale_y_reverse() + theme_minimal(base_size=18) + facet_grid(.~source)
 }

pivot_longer(df_all,cols = 2:14, names_to = "age", values_to = "wt") %>% group_by(age,source) %>% 
   mutate(age=as.numeric(age), mnwt=mean(wt)) %>% ungroup() %>% filter(year>=1982,age<=maxage) %>% mutate(anom=wt/mnwt-1,Anomaly=ifelse(abs(anom)>.5,NA,anom) ) %>%
   ggplot(aes(y=year,x=age,fill=Anomaly,label=round(wt,2))) + geom_tile() + 
   scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
   geom_text(size=3) + ylab("Year") + xlab("Age") + 
   scale_y_reverse() + theme_minimal(base_size=18) + facet_grid(.~source)

#--Constant N-at-age to see weight only impact-----------
   names(df_all)
   unique(df_all$source)
   # Check age composition
   library(ggthemes)
pivot_longer(df_all,cols = 2:14, names_to = "age", values_to = "wt") %>% mutate(source=as.factor(source),N=exp(-.2*as.numeric(age))) %>% 
              mutate(source=factor(source, c("fishery data","survey data","model with survey","model no survey") ) ) %>%
              group_by(year,source) %>% summarise(Biomass=sum(wt*N)) %>% 
             # filter(source %in% c("fishery data","survey data","model with survey")) %>%
              ggplot(aes(x=year,y=Biomass,color=source)) + geom_line(size=2) + theme_few(base_size=18) + ggtitle("Constant numbers-at-age") + 
              expand_limits(y=0)
                #levels = source[c(1,4,2,3))])  #%>%

#--------------------------
# Prediction evaluations from means
#--------------------------
# Read in weights by age 3-15
stwt <- as.data.frame(read.table("~/OneDrive/ebswp/data/sampler/cases/ebswp/statwts.dat",as.is=TRUE))
stwt[11:13] <- 0
stwt
mnwt <- data.table(read.table("~/OneDrive/ebswp/data/sampler/cases/ebswp/wtage2016.dat",as.is=TRUE))
names(mnwt) <- 3:15;
mnwt$yr <- 1991:2015; 
mnwt
mnwt.m <- melt(data.table(mnwt), measure.vars = 1:13, variable.name = "age", value.name = "wt")
mnwt.m$age <- as.numeric(mnwt.m$age) +2
setkey(mnwt.m,yr,age) # NOTE the order of this matters...
mnwt.m$mygam <- as.numeric(stwt)
mnwt.m
reres[[45]]
length(reres)
repred <- data.table(subn=integer(),yr=double(),age=double(),pred=double(),obs=double(),pyr=integer(), ayr=double(),stat=double())[]
for (i in 1:length(reres))
{
  mysub <- ifelse(i<=15,"No survey",ifelse(i<=30,"With survey","Year effects only"))
  if(i!=15|i!=30|i!=45)
  {
    iyr <- reres[[i]]$cur_yr
    pred <- data.table(age=rep(3:15,2),yr=rep(c(iyr,iyr+1),each=13))
    dtmp <- data.table(yr=reres[[i]]$yr,reres[[i]]$wt_pre)
    dtmp
    names(dtmp)
    names(dtmp) <- c("yr",3:15)
    pred <- melt( dtmp ,measure.vars = 2:14, variable.name="age",value.name="pred")
    pred$age <- as.numeric(pred$age)+2
    setkey(pred,yr,age)
    pred <- pred[yr %between% c(iyr,iyr+1)] #,.(wt=mean(wt)),age][,wt]
    pred
    setkey(pred,yr,age)
    t     <- merge(pred,mnwt.m )[,.(pred=pred,obs=wt,stat=sum(mygam*(wt-pred)^2) ),.(yr,age)]
    t$ayr <- iyr; t$pyr <- t$yr - iyr; 
    if (i<=15) 
      t$subn <- mysub #"No survey" # i #paste0("Model_",1*integer(i/15),i)
    else
      t$subn <- mysub #"With survey"   # i #paste0("Model_",1*integer(i/15),i)
    repred   <- rbind(repred,t)
  }
}
repred$pyr <- as.factor(repred$pyr)
#repred$subn <- as.factor(rep(c("with","without"),each=15))
setkey(repred,age,pyr)
repred
mnpred <- data.table(subn=integer(),yr=double(),age=double(),pyr=integer(), ayr=double(),stat=double(),pred=double(),obs=double())[]
for (subn in c(1,3,5,10))
{
  for (i in 2001:2014){
    pred <- data.table(age=rep(3:15,2),yr=c(rep(i,13),rep(i+1,13)))
    #if (subn==0)
    #  pred$pred <- mnwt.m[yr == i,.(wt=mean(wt)),age][,wt]
    #else
    pred$pred <- mnwt.m[yr %between% c(i-subn,i-1),.(wt=mean(wt)),age][,wt]
    setkey(pred,yr,age)
    t     <- merge(pred,mnwt.m )[,.(pred=pred,obs=wt,stat=sum(mygam*(wt-pred)^2) ),.(yr,age)]
    t$ayr <- i; t$pyr <- t$yr - i; t$subn <- subn
    t
    mnpred
    mnpred   <- rbind(mnpred,t)
  }
}
mnpred$pyr <- as.factor(mnpred$pyr)
mnpred$subn <- as.factor(mnpred$subn)
mnpred
mnpred[,.(score=sum(stat)),.(subn,pyr)] %>% ggplot(aes(x=subn,y=score,fill=pyr)) + geom_bar(stat="identity",position="dodge") + mytheme #stat_identity())
setkey(mnpred,yr)
mnpred

#repred[subn=="With survey"&age<18&age>10] %>%  ggplot(aes(x=yr,y=pred,colour=as.factor(age))) + labs(y="Body weight (kg)",x="Year") + 
p <- repred[subn=="With survey"&age<8&age>3] %>%  ggplot(aes(x=yr,y=pred,colour=as.factor(age))) + labs(y="Body weight (kg)",x="Year") + 
                  mytheme + theme(panel.grid.major.x = element_line(colour="grey",linetype="dashed")) +
                  scale_x_continuous(breaks=2001:2015) + annotate("text", x=2001, y=0.49,colour="red",         label="Age 4",size=9) + annotate("text", x=2003, y=0.68,colour="limegreen", label="Age 5",size=9) + annotate("text", x=2002, y=0.82,colour="darkcyan",   label="Age 6",size=9) + annotate("text", x=2001, y=0.95,colour="purple",      label="Age 7",size=9) +
                  geom_point(aes(x=yr-.2*(as.numeric(pyr)-1.5),shape=pyr,size=1.2)) + 
                  geom_point(aes(x=yr,y=obs),shape=8,size=4) + guides(fill=FALSE,shape=FALSE,size=FALSE,colour=FALSE) #+
                  geom_line(data=mnpred[subn=="1"&age<8&age>3],aes(x=yr,y=pred,colour=as.factor(age)),size=1) 
print(p)
t <- repred[subn=="With survey"&age<8&age>3]
p <- ggplot(t, aes(x=yr,y=pred,colour=as.factor(age))) + labs(y="Body weight (kg)",x="Year") + 
                  mytheme + theme(panel.grid.major.x = element_line(colour="grey",linetype="dashed")) +
                  scale_x_continuous(breaks=2001:2015) + annotate("text", x=2001, y=0.49,colour="red", label="Age 4",size=9) + annotate("text", x=2003, y=0.68,colour="limegreen", label="Age 5",size=9) + annotate("text", x=2002, y=0.82,colour="darkcyan",   label="Age 6",size=9) + annotate("text", x=2001, y=0.95,colour="purple",      label="Age 7",size=9) +
                  geom_point(aes(x=yr-.2*(as.numeric(pyr)-1.5),shape=pyr,size=1.2)) + 
                  geom_point(aes(x=yr,y=obs),shape=8,size=4) + guides(fill=FALSE,shape=FALSE,size=FALSE,colour=FALSE) +
                  geom_line(data=mnpred[subn=="1"&age<8&age>3], aes(x=yr,y=pred,colour=as.factor(age)),size=1) + 
                  geom_point(data=mnpred[subn=="1"&age<8&age>3], aes(x=yr-.2*(as.numeric(pyr)-1.5),shape=pyr,size=1.2))  
                  #geom_line(data=mnpred[subn=="1"&age<18&age>10],aes(x=yr,y=pred,colour=as.factor(age),size=1.1)) 

print(p)

t <- rbind(mnpred[subn=="1"],repred[subn=="No survey"])
ggplot(t,aes(x=yr,y=stat,col=pyr,as.factor(age))) + geom_point() + mytheme #stat_identity())
ggplot(t[pyr=="0"],aes(x=yr,y=stat,col=subn)) + geom_point() +facet_grid(age~.) +  mytheme #stat_identity())

#--------------------------
# Prediction evaluations all models...
#--------------------------
tt <- rbind(mnpred[,.(score=sum(stat)),.(subn,pyr)],repred[,.(score=sum(stat)),.(subn,pyr)] )
setkey(tt,pyr,subn)
tt
ggplot(tt,aes(x=subn,y=score,fill=pyr)) +labs(x="Model",y="Weighted score",fill="Projection \n year") + geom_bar(stat="identity",position="dodge") + mytheme #stat_identity())
#--------------------------
stwt

tt
repred

stwt <- data.table(read.table("statwts.dat",as.is=TRUE))
mnwt <- data.table(read.table("~/OneDrive/ebswp/data/sampler/cases/ebswp/wtage2016.dat",as.is=TRUE))
names(mnwt) <- 3:15;mnwt$yr <- 1991:2015; 
stwt
mnwt
mnwt.m <- melt(data.table(mnwt), measure.vars = 1:13, variable.name = "age", value.name = "wt")
mnwt.m$age <- as.numeric(mnwt.m$age) +2
mnwt.m$mygam <- as.numeric(stwt)
setkey(mnwt.m,yr,age)
mnwt.m
str(mnwt.m)
i=2001
subn=1
mnpred[subn=="1"&age<8&age>3] %>%  ggplot(aes(x=yr,y=pred,colour=as.factor(age))) + geom_line(aes(size=1.0)) + geom_point(aes(shape=pyr,size=1.2)) + geom_point(aes(x=yr,y=obs),shape=3,size=4) +mytheme

mnpred[subn=="1"&age<8&age>3] %>%  ggplot(aes(x=yr,y=pred,colour=as.factor(age))) + geom_line() + labs(y="Body weight (kg)",x="Year") + 
                         annotate("text", x=2001, y=0.54,colour="red", label="Age 4",size=9) +
                         annotate("text", x=2003, y=0.70,colour="forestgreen", label="Age 5",size=9) +
                         annotate("text", x=2001, y=0.95,colour="purple", label="Age 7",size=9) +
                         geom_point(aes(shape=pyr,size=1.2)) + geom_point(aes(x=yr,y=obs),shape=8,size=4) +mytheme + guides(fill=FALSE,shape=FALSE,size=FALSE,colour=FALSE) 

mnpred[subn=="1"&age<8&age>3] %>%  ggplot(aes(x=yr,y=pred,shape=pyr,colour=as.factor(age))) + geom_point(aes(size=4)) + geom_point(aes(x=yr,y=obs),shape=3,size=4) +mytheme
mnpred 
str(repred)
repred
repred[subn<=15,subn:="No survey"]
/new

#-----------------------------------------------------------
 
(res[[1]]$yr)
dim(res[[8]]$wt_pre)
names(res[[1]] )
# Compile retrospectives
i=3
df <- data.frame()
yrs <- 1991:2014
for (i in 1:15)
{
  lastyr <- 2017-i-2
  yrsfit <- 1991:lastyr
  yrsprj <- (lastyr+1):(lastyr+3)
  nyrs   <- length(yrsfit) #+ length(yrsprj)
  t <- unlist(res[[i]]$data)
  t <- data.frame(t,yr=yrs,             run=i,src="obs")
  names(t) <-c(3:15,names(t[14:16]))
  r <- unlist(res[[i]]$W[c(yrsfit,yrsprj)-1969,])
  r
  r <- data.frame(r,yr=c(yrsfit,yrsprj),run=i,src=c(rep("est",length(yrsfit)),rep("projected",length(yrsprj)) ) )
  #r <- bind_rows(r,data.frame(t(replicate(3,res[[i]]$mnwt)),yr=yrsprj,run=rep(i,3),src=rep("mean",3)) )
  names(r) <-c(3:15,names(r[14:16]))
  names(r)
  df <- rbind(df,gather(t,age,wt,1:13), gather(r,age,wt,1:13))
}
odf <- df %>% filter(src=="obs",age>3,age<9) # ,run<14) 
tdf <- df %>% filter(src=="projected",age>3,age<9) # ,run<14) 
edf <- df %>% filter(run==1,src=="est",age>3,age<9) # ,run<14) 
tdf
ggplot(odf, aes(x=yr,y=wt)) + geom_point() + geom_line(data=edf,aes(x=yr,y=wt)) + geom_line(data=tdf,aes(x=yr,y=wt,colour=as.factor(run)))  + facet_grid(age ~ .,scale="free_y") + ylab("Weight (kg)") + xlab("Year") + scale_color_discrete(guide=FALSE) #+ geom_hline(data = df, yintercept=(mnwt))
length(res)

ggplot(tdf, aes(x=yr,y=wt,colour=as.factor(run))) + geom_line() + geom_point(data=odf,aes(x=yr,y=wt,colour=as.factor(run))) + facet_grid(. ~ age )


