library('bsts')
library(tidyverse)
model <- list()
kaub <- read_csv("data/kaub.csv")
covariates <- read_csv("data/covariates.csv")
glmnet_covariates <- read_csv("data/selectedcov.csv")


joined <- right_join(kaub,covariates) %>% select(-date)
joined <- joined[12:nrow(joined),]

joined[is.na(joined)] <- 0

joined_glm <- right_join(kaub,glmnet_covariates) %>% select(-date)
joined_glm <- joined_glm[12:nrow(joined_glm),]
joined_glm[is.na(joined_glm)] <- 0
model <- AddLocalLevel(model,joined$level)
model <- AddAr(model,joined$level)
model <- AddSeasonal(model,joined$level,51)

#fit <- bsts(joined$level,
#               state.specification = model,
#               niter = 1000)

#plot(fit)
#plot(fit, "components")

model_dynamic <- AddDynamicRegression(model,level~.,data = joined_glm)

fit_dynamic <- bsts(joined_glm$level,state.specification = model_dynamic,niter = 1000)

