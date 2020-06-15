library(tidyverse)
library(lubridate)
train_test_split <- ymd('2019-01-01')
kaub_data <- read_csv('data/kaub.csv')
mean_temp <- read_csv('data/mean_temp.csv')
precipitation <- read_csv('data/percipitation.csv')
mean_temp <- pivot_wider(mean_temp,names_from = id,values_from = value)
precipitation <- pivot_wider(precipitation,names_from = id,values_from = value)
colnames(mean_temp) <- paste(colnames(mean_temp),'temp',sep='_')
colnames(mean_temp)[1] <- 'date'
colnames(precipitation) <- paste(colnames(precipitation),'regen',sep='_')
colnames(precipitation)[1] <- 'date'
colnames(tsz) <- c('date','kaub')
covariates <- tsz %>% right_join(mean_temp,by='date')

covariates <- covariates %>% right_join(precipitation,by='date')
lags <- 0:12
names(lags) <- paste('lag',0:12,sep='_')
df<- map(covariates,function(col){
  map_df(lags,function(lg,x) {
    lag(x,lg)
  },col)
})
name <- names(df)  
for (i in 1:28) {
  colnames(df[[i]]) <- paste(name[i],colnames(df[[i]]),sep='_')
}
tbl <- df %>% flatten %>% as_tibble()
tbl <- tbl %>% select(-2:-13)
tbl<- tbl %>% mutate(date=ymd(date_lag_0))
tbl <- tbl %>% mutate(day=day(date),year=year(date),month=month(date))
tbl <- tbl %>% select(-date_lag_0)
