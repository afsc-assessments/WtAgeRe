# wt.R
library(tidyverse)
library(data.table)
library(here)
source(here("R/readadmb.R"))
source(here("R/helper.R"))
# First 15 w/o survey
# Second 15 w/ survey
#system(paste0("cp arc/wt2_no_srv.dat wt_",mod_opt,".dat") )
#system(paste0("cp arc/wt2_no_srv.pin wt_",mod_opt,".pin") )
# system(paste0("cp arc/wt2_with_srv.pin wt_",mod_opt,".pin") )
#--------------------------------------------------------------------------------------
# get data (fishery)     
dfgoa<-read_dat(here("examples","goapollock","wtgoa.dat"))

df<-read_dat(here("examples","ebspollock","wt.dat"))
df_fsh <- data.frame(year=df$fshry_yrs,df$fishery)
df_srv <- data.frame(year=df$survey_yrs,df$survey)
df_bth <- rbind(data.frame(df_fsh,source="fishery"),data.frame(df_srv,source="survey"))
dd     <- read_rep("wt.rep") 
dd
df_tmp <- data.frame(year=dd$yr,dd$wt_pre,source="predicted") 
df_bth
names(df_tmp) <- c("year",3:15,"source")
names(df_bth) <- c("year",3:15,"source")
df_all <-rbind(df_bth,df_tmp) %>% mutate(source=fct_relevel(source, "fishery","survey","predicted"))
maxage=10 # For display only...
p1<-pivot_longer(df_all,cols = 2:14, names_to = "age", values_to = "wt") %>% group_by(age,source) %>% 
   mutate(age=as.numeric(age), mnwt=mean(wt)) %>% ungroup() %>% filter(age<=maxage) %>% mutate(anom=wt/mnwt-1,Anomaly=ifelse(abs(anom)>.5,NA,anom) ) %>%
   ggplot(aes(y=year,x=age,fill=Anomaly,label=round(wt,2))) + geom_tile() + 
   scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
   geom_text(size=3) + ylab("Year") + xlab("Age") + 
   scale_y_reverse() + theme_minimal(base_size=18) + facet_grid(.~source) ; p1
  ggsave("~/_mymods/ebswp/doc/figs/fsh_wtage_data_pred.pdf",plot=p1,width=12,height=15.0,units="in")
#--------------------------------------------------------------------------------------

  re_dat <- read_dat(here("examples", "ebspollock","wt.dat"))
names(re_dat)
re_dat
do_est = TRUE
ij=3
for (ij in 1:3)
{
  mod_opt <- ifelse(ij==1,"nosrv",
             ifelse(ij==2,"both",
             ifelse(ij==3,"yreff","coheff")))
  if (ij>1)
  {
    system(paste0("cp arc/wt2_with_srv.dat wt_",mod_opt,".dat") )
    system(paste0("cp wt_both.pin wt_",mod_opt,".pin") )
  }
  else
  {
    mod_opt="both"
    system(paste0("cp arc/wt2_no_srv.dat wt_",mod_opt,".dat") )
    system(paste0("cp arc/wt2_no_srv.pin wt_",mod_opt,".pin") )
  }
  d <- read_dat(paste0("wt_",mod_opt,".dat"))
  print(mod); print(mod_opt)
 df
 mod_opt="both"
 i=2021
 do_est=TRUE
  for (i in 2019:2010){
    df$cur_yr <- i
    print(i)
    write_dat(df,paste0("wt_",mod_opt,"_",i,".dat"))
    fname <- paste0("arc/wt_",mod_opt,"_",i)
    if (do_est)
    {
  	  system(paste0("./wt -ind wt_",mod_opt,"_",i,".dat -nox  "))
  	  system(paste0("cp wt.rep ", fname,".rep"))
      system(paste0("cp wt.std ", fname,".std"))
      system(paste0("cp wt.par ", fname,".par"))
    }
  	reres[[i]] <- read_admb(fname)
  }
}
# 1st 15 is w/o survey, second w/ and both, 3rd is w/ and yreff, 4th is w/ and coheff
system(paste0("./wt  "))
getwd()
length(reres)
df<-read_dat(here("examples","ebspollock","wt.dat"))
names(df)
df
df[8]


# Anomalies on the data
maxage=10 # For display only...
pivot_longer(df_bth,cols = 2:14, names_to = "age", values_to = "wt") %>% group_by(age,source) %>% 
   mutate(age=as.numeric(age), mnwt=mean(wt)) %>% ungroup() %>% filter(age<=maxage) %>% mutate(anom=wt/mnwt-1,Anomaly=ifelse(abs(anom)>.5,NA,anom) ) %>%
   ggplot(aes(y=year,x=age,fill=Anomaly,label=round(wt,2))) + geom_tile() + 
   scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
   geom_text(size=3) + ylab("Year") + xlab("Age") + 
   scale_y_reverse() + theme_minimal(base_size=18) + facet_grid(.~source)

df_fsh
maxage=13
pivot_longer(df_fsh,-c(year), names_to = "age", values_to = "wt") %>% group_by(age) %>% 
   mutate(age=as.numeric(age), mnwt=mean(wt)) %>% ungroup() %>% filter(age<=maxage) %>% mutate(anom=wt/mnwt-1,Anomaly=ifelse(abs(anom)>.5,NA,anom) ) %>%
   ggplot(aes(y=year,x=age,fill=Anomaly,label=round(wt,2))) + geom_tile() + 
   scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
   geom_text(size=2.3) + ylab("Year") + xlab("Age") +
   ylim(c(2020,1960))+  theme_minimal(base_size=18) + ggtitle("Saithe data")
function(file=c("wt"),dat=df,source=c("model")){
    df_tmp <- data.frame()
    #for (i in 1:length(file)){

#--------------------------
# Read in results/estimates
#--------------------------
  fn_get_pred
  cd <- getwd()
  setwd("examples/ebspollock/")
  df_pred <- fn_get_pred(file=c("examples/ebspollock/wt"),source=c("Pollock"))
  df_pred <- fn_get_pred(file=c("examples/ebspollock/wt"),dat=df_bth,source=c("Pollock"))
  tmp <- read_rep("examples/ebspollock/wt.rep")
  df_pred <- data.frame(year=tmp$yr,wt=tmp$wt_pre,source="new")
  names(df_pred) <- c("year",3:15,"source")
  fn_plot_anoms(df_pred)
  fn_plot_anoms
  df_pred

pivot_longer(df,-c(year,source), names_to = "age", values_to = "wt") %>% group_by(age) %>% 
   mutate(age=as.numeric(age), mnwt=mean(wt)) %>% ungroup() %>% filter(age<=maxage) %>% mutate(anom=wt/mnwt-1,Anomaly=ifelse(abs(anom)>.5,NA,anom) ) %>%
   ggplot(aes(y=year,x=age,fill=Anomaly,label=round(wt,2))) + geom_tile() + 
   scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
   geom_text(size=2.3) + ylab("Year") + xlab("Age") + 
   scale_y_reverse() + theme_minimal(base_size=18) + ggtitle("Saithe model")

#--------------------------
# Read in results/estimates
#--------------------------
  df_pred <- fn_get_pred(file=c("wt"),source=c("Pollock"))
  df_pred <- fn_get_pred(file=c("wt"),source=c("Saithe"))
  fn_plot_anoms(df_pred)
  df_pred


pivot_longer(df_pred,cols = 2:11, names_to = "age", values_to = "CV") %>% 
   mutate(age=as.numeric(age)) %>% filter(age<=maxage,type=="CV",Year>1970) %>% 
   ggplot(aes(y=Year,x=age,fill=CV,label=round(CV,2))) + geom_tile() + 
   scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
   geom_text(size=3) + ylab("Year") + xlab("Age") + ggtitle("Weight-at-age CV") +
   scale_y_reverse() + theme_minimal(base_size=16) #+ facet_grid(.~source)

   # show some examples
  tm <- paste0("arc/wt_both_",c(2010,2013,2017,2020))
  hm <- paste0("Model_",c(2010,2013,2017,2020))
  df_retro <- fn_get_pred(file=tm,source=hm)
  fn_plot_anoms(df_retro)

  tt<- read_admb("wt")  
  names(tt)
  tt$fit$names 
  wt_hist <- data.frame(matrix(tt$fit$est[grep("wt_hist",tt$fit$names )],ncol=13,byrow=TRUE),Year=tt$yr,"est")
  names(wt_hist) <- c(3:15,"Year","type")
  wt_std <- data.frame(matrix(tt$fit$std[grep("wt_hist",tt$fit$names )],ncol=13,byrow=TRUE),Year=tt$yr,"std")
  names(wt_std) <- c(3:15,"Year","type")
  wt_cv <- data.frame(wt_std[,1:13]/ wt_hist[,1:13],Year=tt$yr,"CV")
  names(wt_cv) <- c(3:15,"Year","type")
  wt_res<- rbind(wt_hist,wt_std,wt_cv)
  idx_coh <- grep("coh_eff",tt$fit$names )
  idx_coh[-1]
  names(tt)
  df_coh <- data.frame(Year=tt$cohort,est=tt$fit$est[idx_coh[-1]], se=tt$sigma_coh*tt$fit$std[idx_coh[-1]]) %>% mutate(lb=est-2*se, ub=est+2*se )
  ggplot(df_coh, aes(x=Year, y=est, ymax=ub, ymin=lb )) + geom_ribbon(fill="salmon",col="red") + geom_hline(yintercept=0,col="grey") + 
  geom_line(col="blue") + geom_point() + ggtitle("Cohort efffect, EBS pollock") + theme_minimal()
  ggplot(df_coh, aes(x=Year, y=est, ymax=ub, ymin=lb )) + geom_ribbon(fill="salmon",alpha=0.5,col="red") + 
  geom_hline(yintercept=0,col="grey") + 
  geom_line(col="blue") + geom_point() + ggtitle("Cohort efffect, saithe") + theme_minimal(base_size = 18)

  names(wt_hist)
  wt_hist %>% select(age3="3",Year) %>% transmute(wt=age3, coh=Year-3, anom=wt/mean(wt)-1) %>% 
      ggplot(aes(x=coh,y=anom)) + geom_line(size=1) + theme_minimal()

  df_rep <- read_rep(paste0("wt.rep"))
  df_pred<- data.frame(year=df_rep$yr,df_rep$wt_pre)
  names(df_pred) <- c("year",3:15)

  df_nosrv<-read_rep("wt_nosrv.rep")
  df_nosrv<- data.frame(year=df_nosrv$yr,df_nosrv$wt_pre)
  names(df_nosrv) <- c("year",3:15)

maxage=10
pivot_longer(df_pred,cols = 2:14, names_to = "age", values_to = "wt") %>% group_by(age) %>% 
   mutate(age=as.numeric(age), mnwt=mean(wt)) %>% ungroup() %>% filter(age<=maxage) %>% mutate(anom=wt/mnwt-1,Anomaly=ifelse(abs(anom)>.5,NA,anom) ) %>%
   ggplot(aes(y=year,x=age,fill=Anomaly,label=round(wt,2))) + geom_tile() + 
   scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
   geom_text(size=3) + ylab("Year") + xlab("Age") + 
   scale_y_reverse() + theme_minimal(base_size=18) #+ facet_grid(.~source)

df_coh <- data.frame(cohort = df_rep$cohort, Effect=df_rep$coh_eff[1:(length(df_rep$cohort))])
ggplot(df_coh,aes(x=cohort,y=Effect)) + geom_line(size=2) + theme_bw()
df_all <- rbind(data.frame(df_fsh,source="fishery"),data.frame(df_srv,source="index"),data.frame(df_pred,source="model"))
df_all <- rbind(data.frame(df_fsh,source="fishery data"),data.frame(df_srv,source="survey data"),data.frame(df_pred,source="model with survey"),data.frame(df_nosrv,source="model no survey"))
names(df_all) <- c("year",3:15,"source")
tail(df_all)


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


