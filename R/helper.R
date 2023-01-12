# Functions
   file=c("wt");source=c("model");i=1
fn_get_pred <- function(file=c("wt"),dat=df,source=c("model")){
    df_tmp <- data.frame()
    for (i in 1:length(file)){
      dd     <- read_rep(paste0(file[i],".rep"))
      df_tmp <- rbind(df_tmp, data.frame(year=dd$yr,dd$wt_pre,source=source[i]))
    }
    names(df_tmp) <- c("year",unlist(dat[7]):unlist(dat[8]),"source")
    return(df_tmp)
  }
  
fn_plot_anoms <- function(dfin, maxage=10,firstyr=1982,minage=3){
  p1<-pivot_longer(dfin,cols=2:11, names_to = "age", values_to = "wt") %>% group_by(age,source) %>% 
  mutate(age=as.numeric(age), mnwt=mean(wt)) %>% ungroup() %>% filter(year>=firstyr,age>=minage,age<=maxage) %>% 
   mutate(anom=wt/mnwt-1,Anomaly=ifelse(abs(anom)>.5,NA,anom) ) %>%
   ggplot(aes(y=year,x=age,fill=Anomaly,label=round(wt,2))) + geom_tile() + 
     scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
     geom_text(size=3) + ylab("Year") + xlab("Age") + 
     scale_y_reverse() + theme_minimal(base_size=18) 
     if (length(unique(dfin$source))>1) p1 <- p1 + facet_grid(.~source)
      return(p1)
 }
 #head(df_res)
#fn_plot_anoms(df_pred)
#dfin <- df_res
#names(dfin)[3:dim(dfin)[2]-1] 
#names(dfin)[2]