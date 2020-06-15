library(tidyverse)
library(lubridate)
library(parsnip)
train_test_split <- ymd('2019-01-01')
kaub_data <- read_csv('data/kaub.csv')
kaub <- kaub_data %>% mutate(date=ymd(date))

kaub <- kaub %>% mutate(day=day(date),
                        month=month(date),
                        year = year(date))

lags <- 1:12
names(lags) <- 1:12 %>% map(~{paste0('lag_',.)})
df <- lags %>% map_df(~{lag(kaub$level,.)})

kaub <- cbind(kaub,df)
kaub <- kaub[13:nrow(kaub),]
kaub_train <- kaub %>% filter(date<train_test_split)
kaub_test <-  kaub %>% filter(date>train_test_split)
kaub_train <- kaub_train %>% select(-date)
kaub_test <- kaub_test %>% select(-date)

model <- boost_tree(
  mode='regression',
  mtry= 20,
  trees=500,
  min_n=3,
  tree_depth = 8,
  learn_rate = 0.01,
  loss_reduction = 0.01) %>% set_engine('xgboost')

fit <- model %>% fit.model_spec(level ~ . , data=kaub_train)
pred <- predict(fit,kaub_test[,-1])$.pred
