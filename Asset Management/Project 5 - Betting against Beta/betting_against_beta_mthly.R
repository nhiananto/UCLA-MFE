library(RPostgres)
library(tidyverse)
library(zoo)
require(data.table)
library(lubridate)
library(kableExtra)
library(dplyr)
library(foreign)
library(crunch)
library(moments)

setwd("~/UCLA MFE/Spring 2020/Quantitative Asset Management/Final Project")


wrds <- dbConnect(Postgres(),
                  host='wrds-pgdata.wharton.upenn.edu',
                  port=9737,
                  dbname='wrds',
                  sslmode='require')

crsp_sql <-dbSendQuery(wrds,"select a.permno,
                  c.shrcd,
                  c.exchcd,
                  a.date,
                  a.prc,
                  a.ret,
                  d.vwretd,
                  e.t90ret,
                  e.t30ret
           from crsp.msf a
           inner join crsp.msenames c
           on
            a.permno = c.permno
            and a.date between c.namedt and coalesce(c.nameendt, CURRENT_DATE)
            left join crsp.msi d
            on
            a.date=d.date
            left join crsp.mcti e
            on
              date_trunc('month',a.date)=date_trunc('month', e.caldt)
          where
            a.date >= '1925-12-31' and a.date<='2019-12-31'
           ")


crsp <- dbFetch(crsp_sql, n = -1)

dbClearResult(crsp_sql)

crsp <- data.table(crsp)



###########################################

#check uniqueness
crsp[, .N, by = .(permno, date)][N > 1]

###
crsp$date = as.Date(crsp$date)

#make sure all dates are the same (set as end of month)
crsp[,date := ceiling_date(date, unit = "month") - days(1)]

#### create monthly variable
crsp[, yrmo := make_date(year(date), month(date), 1)]

#lag 1, 5 year flag
crsp[, `:=`(lag_year_1 = date %m-% months(12),
            lag_year_5 = date %m-% months(60))]

#clean returns
crsp[ret %in% c(-66, -77, -88, -99), ret := NA]

#create excess return variable
crsp[, excess_ret := ret - t30ret]

#sharecode
crsp <- crsp[shrcd %in% c(10,11)]

#create log ret variable
crsp[, `:=`(log_ret = log(1 + ret),
            log_excess_ret = log(1 + excess_ret))]


#####################################################
#volatility asset
setorder(crsp, permno, date)

#date has to be >= lagged 1 year (of current record)
#use excess return as volatility
crsp[, `:=`(vol_asset = rollapplyr(.SD, 12,
                                   function(z) sd(fifelse(z[,2] >= z[12,3], 1, as.numeric(NA)) * as.numeric(z[,1])  , na.rm = T),
                                   by.column = F,
                                   fill = NA)
), by = .(permno), .SDcols = c("log_excess_ret", "date", "lag_year_1")]

#create valid vol estimate flag
crsp[ , `:=`(is_valid_vol_est = rollapplyr(.SD, 12,
                                           function(z) sum( !is.na(z[,1]) * as.numeric(z[,2] >= z[12,3])  ) >= 12,
                                           by.column = F,
                                           fill = NA) #ret not missing & within 1 year lag & 12 observations within the year
), by = .(permno), .SDcols = c("log_excess_ret", "date", "lag_year_1") ]

#set as False flag
crsp[is.na(is_valid_vol_est), is_valid_vol_est := F]


#####################################################
#volatility market, create list of dates and market vwretd (assume days are consecutive)
vol_m <- unique(crsp[!is.na(vwretd), .(date, vwretd, t30ret)])

setorder(vol_m, date)
#check any missing vol mkt
vol_m[is.na(vwretd)]


#check dups
vol_m[, .N, by = date][N > 1]

setorder(vol_m, date)

#excess mkt
vol_m[, excess_mkt := vwretd - t30ret]

#log mkt ret
vol_m[, `:=`(log_mkt_ret = log(1+vwretd),
             log_excess_mkt_ret = log(1+ excess_mkt))]


#calculate vol mkt log ret
#use excess mkt ret
vol_m[, vol_mkt := frollapply(log_excess_mkt_ret, 12, function(x) sd(x, na.rm = T))]


######### merge back to original dataset to get market vol
crsp <- merge(crsp, vol_m[, .(date, vol_mkt, log_mkt_ret, log_excess_mkt_ret)], by = c("date"), all.x = T)


#####################################################
setorder(crsp, permno, date)

###########
#correlation
crsp$log_excess_ret = as.numeric(crsp$log_excess_ret)
crsp$log_excess_mkt_ret = as.numeric(crsp$log_excess_mkt_ret)



#use excess return for correlation
#crsp[, correlation := NULL]
crsp[, correlation := as.numeric(NA)]
crsp[, correlation := rollapplyr(.SD, 60, #5 years
                              function(z) cor( ( fifelse(z[,3] >= z[60 ,4], as.numeric(1), as.numeric(NA), as.numeric(NA))  * as.numeric(z[,1]) ),
                                               as.numeric(z[,2]),
                                               use = "na.or.complete"), #only include records in the past 5 years
                                 by.column = F,
                                 fill = NA)
     ,by = .(permno), .SDcols = c("log_excess_ret", "log_excess_mkt_ret", "date", "lag_year_5")]

#create valid correlation flag
crsp[ , `:=`(is_valid_cor_est = rollapplyr(.SD, 60,
                                           function(z) sum( !is.na(z[,1]) * as.numeric(z[,2] >= z[60,3])  ) >= 36,
                                           by.column = F,
                                           fill = NA) #ret not missing & within 5 year lag & more than 36 observations (for past 60 months (5 years))
), by = .(permno), .SDcols = c("log_excess_ret", "date", "lag_year_5")]


#set valid flag to F for na cor
crsp[is.na(is_valid_cor_est), is_valid_cor_est := F]


###save dataset
write.csv.gz(crsp, "crsp_mthly_excess.csv.gz")


# library(crunch)
# crsp <- fread("~/UCLA MFE/Spring 2020/Quantitative Asset Management/Final Project/crsp_mthly_excess.csv.gz")

##### beta
crsp[, beta := fifelse(is_valid_cor_est & is_valid_vol_est, correlation * vol_asset / vol_mkt, as.numeric(NA))]

##### shrinked beta
crsp[, shrinked_beta := 0.6 * beta + 0.4]


###### get lagged beta & check valid last month flag
crsp[, yrmo := make_date(year(date), month(date), 1)]
crsp[, lag_shrinked_beta := fifelse(shift(yrmo, 1) == (yrmo %m-% months(1)), shift(shrinked_beta), as.numeric(NA) ), by = permno]

#delete with no lagged beta
crsp <- crsp[!is.na(lag_shrinked_beta)]

############################################################################
##### NYSE breakpoints by month using last observation beta
crsp[, beta_rank := frank(lag_shrinked_beta), by = date]

###deciles
#NYSE breakpoints using size (deciles)
crsp[, beta_decile := cut(beta_rank, breaks= quantile(fifelse(exchcd == 1, beta_rank, as.numeric(NA)),
                                                             probs=seq(0,1,0.1),na.rm=TRUE) , labels=F, right = F), by = date]

# crsp[beta_decile > 10, beta_decile:=NA]


#carry
setorder(crsp, date, beta_rank)
crsp[!is.na(beta_rank), beta_decile := na.locf(beta_decile, na.rm = F), by = .(date)]
#also backward
crsp[!is.na(beta_rank), beta_decile := na.locf(beta_decile, na.rm = F, fromLast = T), by = .(date)]


####################
#monthly portfolio

#create equally weighted portfolio
crsp_monthly = crsp[!is.na(beta_decile) & !is.na(ret), .(eq_ret_beta = mean(ret),
                                                         rf = last(t30ret)) ,by = .(beta_decile, date)]
crsp_monthly[, excess_ret := eq_ret_beta - rf]

crsp_monthly[ , `:=`(Year = year(date),
                     Month = month(date))]

setorder(crsp_monthly, date, beta_decile)
#####
summary_ret <- function(return_tbl, ret_name, start_date = "1925-01-01", end_date = "2012-03-31"){
  returns <- data.table::copy(return_tbl)
  
  returns <- returns[data.table::between(make_date(Year, Month, 1), start_date, end_date), ]
  
  result <-returns[, .(mean_ret = mean(get(ret_name)) * 100,
                       ann_ret = mean(get(ret_name)) * 12,
                       vol = sd(get(ret_name)) * sqrt(12),
                       skewness = skewness(get(ret_name)) ), by = beta_decile][, sr := ann_ret/vol]
  
  result <- t(result)
  result <- result[c(2,4,6,5), ]
  colnames(result) <- c(paste0("Decile ", 1:10))
  rownames(result) <- c("Mean Excess Return", "Standard Deviation", "Sharpe Ratio", "Skewness(m)") 
  return(result)
}




######################################################
#BAB factor
######################################################
#create high low beta weights
crsp[!is.na(beta_rank), `:=`(high_beta_weight = fifelse(beta_rank - mean(beta_rank) > 0 , beta_rank - mean(beta_rank), 0)/
                               sum( abs(beta_rank - mean(beta_rank)) )* 2,
                             low_beta_weight = fifelse(beta_rank - mean(beta_rank) < 0 , abs(beta_rank - mean(beta_rank)), 0)/
                               sum( abs(beta_rank - mean(beta_rank)) ) * 2
), by = date]




port_bab <- crsp[, .(low_beta = sum(low_beta_weight * lag_shrinked_beta, na.rm = T),
         low_beta_excess_ret = sum(low_beta_weight * excess_ret, na.rm = T),
         high_beta = sum(high_beta_weight * lag_shrinked_beta, na.rm = T),
         high_beta_excess_ret = sum(high_beta_weight * excess_ret, na.rm = T)), by = .(date)]

port_bab[, bab_excess_ret := 1/low_beta * low_beta_excess_ret - 1/high_beta * high_beta_excess_ret]
port_bab[, `:=`(Year = year(date),
            Month = month(date))]

bab_res <- function(start_date = "1925-01-01", end_date = "2012-12-31"){
  res <- port_bab[date >= start_date & date <= end_date, .(mean_ret = mean(bab_excess_ret) * 100,
                                                          ann_ret = mean(bab_excess_ret) * 12,
                                                          vol = sd(bab_excess_ret) * sqrt(12),
                                                          skewness = skewness(bab_excess_ret) )][, sr := ann_ret/vol]
  return(t(as.matrix(res))[c(1,3,5,4)])
}






#results table 3 1926-2012
tbl_3_res <- cbind(summary_ret(crsp_monthly, "excess_ret"), bab_res())
colnames(tbl_3_res) <- c(head(colnames(tbl_3_res), -1), "BAB Factor")

#results 2007-2020
output_2 <- cbind(summary_ret(crsp_monthly, "excess_ret", "2007-01-01", "2019-12-31"), bab_res("2007-01-01", "2019-12-31"))
colnames(output_2) <- c(head(colnames(output_2), -1), "BAB Factor")

#results 2015-2020
output_3 <- cbind(summary_ret(crsp_monthly, "excess_ret", "2015-01-01", "2019-12-31"), bab_res("2015-01-01", "2019-12-31"))
colnames(output_3) <- c(head(colnames(output_3), -1), "BAB Factor")


########################################################


########
#Fama French 3 factors
ff_3 = as.data.table(read.csv("F-F_Research_Data_Factors.CSV", stringsAsFactors = F, skip = 3, nrows = 1122))
colnames(ff_3) <- c("date", "mkt.rf", "smb","hml", "rf")
ff_3 <- ff_3[,  c( list("date" = date), lapply(.SD,  function(x) as.numeric(x)/100) ) , .SDcols = -c("date")]
ff_3[, date := as.character(date)][, date := make_date(as.numeric(substring(date, 1,4)), as.numeric(substring(date, 5,6)), day = 1)]
ff_3[, `:=`(Year = year(date),
          Month = month(date))]


#####Fama French 5 factors
ff_5 <- as.data.table(read.csv("F-F_Research_Data_5_Factors_2x3.csv", skip = 3, nrows = 682))
colnames(ff_5) <- c("date", tolower(colnames(ff_5)[-1]))
ff_5 <- ff_5[,  c( list("date" = date), lapply(.SD,  function(x) as.numeric(x)/100) ) , .SDcols = -c("date")]
ff_5[, date := as.character(date)][, date := make_date(as.numeric(substring(date, 1,4)), as.numeric(substring(date, 5,6)), day = 1)]
ff_5[, `:=`(Year = year(date),
          Month = month(date))]

######merge with original dataset 1963+
port_ff_1963 <- merge(crsp_monthly[date >= "1963-07-01" & date <= "2012-03-31", .(date, eq_ret_beta, beta_decile, Year, Month)],
      ff_5[, -c("date")], by = c("Year","Month"), all.x = T)
port_ff_1963[, excess_ret := eq_ret_beta - rf]

######merge with original dataset 1963+
port_bab_ff_1963 <- merge(port_bab[date >= "1963-07-01" & date <= "2012-03-31", .(date, bab_excess_ret, Year, Month)],
                          ff_5[, -c("date")], by = c("Year","Month"), all.x = T)
setnames(port_bab_ff_1963, "bab_excess_ret", "excess_ret")


#########################################################
#2007
######merge with original dataset 2007+
port_ff_2007 <- merge(crsp_monthly[date >= "2007-01-01" & date <= "2015-12-31", .(date, eq_ret_beta, beta_decile, Year, Month)],
                      ff_5[, -c("date")], by = c("Year","Month"), all.x = T)
port_ff_2007[, excess_ret := eq_ret_beta - rf]

######merge with original dataset 2007+
port_bab_ff_2007 <- merge(port_bab[date >= "2007-01-01" & date <= "2015-12-31", .(date, bab_excess_ret, Year, Month)],
                          ff_5[, -c("date")], by = c("Year","Month"), all.x = T)
setnames(port_bab_ff_2007, "bab_excess_ret", "excess_ret")


#########################################################
#2015
######merge with original dataset 2015+
port_ff_2015 <- merge(crsp_monthly[date >= "2015-01-01" & date <= "2019-12-31", .(date, eq_ret_beta, beta_decile, Year, Month)],
                      ff_5[, -c("date")], by = c("Year","Month"), all.x = T)
port_ff_2015[, excess_ret := eq_ret_beta - rf]

######merge with original dataset 2015+
port_bab_ff_2015 <- merge(port_bab[date >= "2015-01-01" & date <= "2019-12-31", .(date, bab_excess_ret, Year, Month)],
                          ff_5[, -c("date")], by = c("Year","Month"), all.x = T)
setnames(port_bab_ff_2015, "bab_excess_ret", "excess_ret")




reg<- function(data){
  
  model3<-lm(excess_ret~mkt.rf+smb+hml,data=data)
  model4<-lm(excess_ret~mkt.rf+smb+hml+rmw,data=data)
  model5<-lm(excess_ret~mkt.rf+smb+hml+rmw+cma,data=data)
  
  result3<-c(summary(model3)$coefficients[1]*100,coef(summary(model3))[1,"t value"])
  result4<-c(summary(model4)$coefficients[1]*100,coef(summary(model4))[1,"t value"])
  result5<-c(summary(model5)$coefficients[1]*100,coef(summary(model5))[1,"t value"])
  
  result<-c(result3,
            result4,
            result5,
            sd(data$excess_ret)/sd(data$mkt)*cor(data$mkt,data$excess_ret))
  return(result)
}



alpha_dt<- function(data){
  result <- matrix(data=NA,nrow=7,ncol=10)
  for (i in c(1:10)){
    result[,i]<-reg(data[beta_decile==i,])
  }
  result <- data.table(result)
  names(result)<- c(as.character(paste0("P",seq(1,10,1))))
  result<- data.table(" "=c("Three-factor alpha","t-value","Four-factor alpha","t-value","Five-factor alpha","t-value","Beta(realized)"),result)
  return(result)
  
}


###1963
output_alpha_1963 <- data.table(alpha_dt((port_ff_1963)),"BAB"=reg(port_bab_ff_1963))


###2007
output_alpha_2007 <- data.table(alpha_dt((port_ff_2007)),"BAB"=reg(port_bab_ff_2007))

###2015
output_alpha_2015 <- data.table(alpha_dt((port_ff_2015)),"BAB"=reg(port_bab_ff_2015))



################################################
#save results output
################################################

write.csv(tbl_3_res,"output_1926.csv")
write.csv(output_2,"output_2007.csv")
write.csv(output_3,"output_2015.csv")

######################
write.csv(output_alpha_1963, "output_alpha_1963.csv")
write.csv(output_alpha_2007, "output_alpha_2007.csv")
write.csv(output_alpha_2015, "output_alpha_2015.csv")



############### cum ret deciles
forplot=data.table::copy(crsp_monthly[Year >=1926 & Year <=2012])
forplot[,Dates:=as.Date(paste0(Year,"-",Month,"-","01"),"%Y-%m-%d")]
forplot=forplot[order(beta_decile,Dates)]
forplot[,cumret:=log(cumprod(1+excess_ret)),by=.(beta_decile)]

ggplot() +
  geom_line(aes(x = Dates, y = cumret, col = factor(beta_decile) ), data = forplot) +
  xlab("Time Period") + ylab("Cumulative Return") + ggtitle("Cumulative Return since 1926")+theme_bw()+
  theme(plot.title = element_text(size = 8, face = "bold",hjust=0.5),
        legend.position = "top",
        legend.justification = c("center", "top"),
        legend.box.just = "right",
        legend.margin = margin(1, 1, 1, 1))



################# cum ret bab

plot_bab<- port_bab[date >= "1926-01-01" & date <= "2012-12-31"]
plot_bab[,cum_ret:=log(cumprod(1+bab_excess_ret))]

ggplot(plot_bab,aes(x=date, y=cum_ret))+
  geom_line(col="navyblue")+
  geom_line(aes(y = mean(cum_ret)), color = "red", linetype = "dotted") +
  geom_line(aes(y = quantile(cum_ret, 0.75)), color = "black", linetype = "dashed") +
  geom_line(aes(y = quantile(cum_ret, 0.25)), color = "black", linetype = "dashed")+
  ggtitle(" Monthly Return of BAB")+
  xlab(" Investment Horizen")+
  ylab(" return")+
  ggtitle("Monthly Cumulative Returns for BAB factor")+
  theme(plot.title = element_text(size = 8, face = "bold",hjust=0.5))


ggplot(plot_bab,aes(x=date, y=bab_excess_ret))+
  geom_line(col="navyblue")+
  ggtitle(" Monthly Return of BAB")+
  xlab(" Investment Horizen")+
  ylab(" return")+
  ggtitle("Monthly Returns for BAB factor")+
  theme(plot.title = element_text(size = 8, face = "bold",hjust=0.5))
