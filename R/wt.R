# Functions
   file=c("wt");source=c("model");i=1
fn_get_pred <- function(file=c("wt"),dat=df,source=c("model")){
    df_tmp <- data.frame()
    for (i in 1:length(file)){
      dd     <- read_rep(paste0(file[i],".rep"))
      df_tmp <- rbind(df_tmp, data.frame(year=dd$yr,dd$wt_pre,source=source[i]))
    }
    names(df_tmp) <- c("year",unlist(dat[8]):unlist(dat[9]),"source")
    return(df_tmp)
  }


fn_plot_anoms <- function(df, maxage=10,firstyr=1982,minage=3){
  pivot_longer(df,cols = 2:52, names_to = "age", values_to = "wt") %>% group_by(age,source) %>% 
   mutate(age=as.numeric(age), mnwt=mean(wt)) %>% ungroup() %>% filter(year>=firstyr,age>=minage,age<=maxage) %>% mutate(anom=wt/mnwt-1,Anomaly=ifelse(abs(anom)>.5,NA,anom) ) %>%
   ggplot(aes(y=year,x=age,fill=Anomaly,label=round(wt,2))) + geom_tile() + 
   scale_fill_gradient2(low = scales::muted("blue"), high = scales::muted("red"), na.value = "white") +
   geom_text(size=3) + ylab("Year") + xlab("Age") + 
   scale_y_reverse() + theme_minimal(base_size=18) + facet_grid(.~source)
 }