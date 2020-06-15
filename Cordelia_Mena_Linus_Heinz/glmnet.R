library(tidyverse)
library(tidymodels)

## Data preparation
kaub <- read_csv("data/kaub.csv")
covariates <- read_csv("data/covariates.csv")

joined <- right_join(kaub,covariates) 
joined <- joined[12:nrow(joined),]## Train test split
traintest<- initial_time_split(joined, prop = 3/4)

ts_train <- training(traintest)
ts_test  <- testing(traintest)



### Preprocessing 
control <- recipe(level ~.,data = traintest) %>%
  step_rm(date) %>% 
  step_center(all_numeric())%>%
  step_knnimpute(all_numeric()) %>% 
  step_scale(all_numeric()) %>% 
  prep(training = ts_train, retain = TRUE)

train_data <- juice(control)
test_data  <- bake(control, ts_test)


## model building
glm_model <- linear_reg(mixture = 0.1) %>%  
  set_engine("glmnet",fit_intercept = F, class='gaussian')
## workflow


glm_model.fit <- parsnip::fit(glm_model,level ~. ,data = train_data)
pred <- multi_predict(glm_model.fit,test_data) %>% unnest() %>% 
  cbind(select(test_data,level))

# calculate mase
mase <- pred %>%
  group_by(penalty) %>%
  mase(.pred, level) %>% arrange(.estimate)

# select best performing model and size
results <- broom::tidy(glm_model.fit) 
results %>% group_by(lambda) %>% count() %>% arrange(n)
r <- results %>% filter(lambda<0.98 & lambda>0.96)
factors <- r$term[str_detect(r$term,'regen')] %>% str_remove_all('`')
covariates_glm_selected <- covariates[c(factors,'date')]
write_csv(covariates_glm_selected,'data/selectedcov.csv')
