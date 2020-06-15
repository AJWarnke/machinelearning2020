install.packages("tidyverse")
install.packages("tseries")
install.packages("forecast")
install.packages("pbapply")
install.packages("ggplot2")
install.packages("caret")
install.packages("glmnet")
install.packages("doParallel")
install.packages("rdwd")
install.packages("MCS")

library(tidyverse)
library(tseries)
library(lubridate)
library(forecast)
library(pbapply)
library(ggplot2)
library(caret)
library(glmnet)
library(rdwd)
library(MCS)

library(doParallel)


#import the water level data
water_levels = read_tsv("rawdata.txt")

plot(water_levels)
#stationary? ...not sure

#returns the first differences of a series. (first observation is dropped because there's no 0th observation to compute the difference with)
firstDifferences = function(series) {
  return(series[2:nrow(series), ] - series[1:(nrow(series)-1), ])
}


first_diffs = firstDifferences(water_levels[, 2])

#dates of first differences begin with the second date because the first observation got dropped (see above)
dates = water_levels[2:nrow(water_levels), 1]

plot(cbind(dates, first_diffs))
#looks definitely stationary now

#augmented Dickey-Fuller test to verify
adf.test(as.matrix(first_diffs))
#->null of non-stationarity refuted at >99% confidence level


#transforms a time series into a series of "observations":
#outcome variable is the first variable of the original series
#predictors are [num_lags] lagged values of each original variable, starting with the [horizon]'th lag
#...setting horizon to 10 equals 10-step-ahead forecasts
addPredictors = function(series, num_lags, horizon=10) {
  
  #add the lags (drops first [num_lags + horizon] rows ...those would have missing lags)
  lagged = embed(as.matrix(series), num_lags+horizon+1)
  
  #copy and repeat names from original series
  colnames(lagged) = rep(names(series), num_lags+horizon+1)
  #turns repeated names (name, name, name) into (name, name.1, name.2) etc.
  #->also denotes lag order
  colnames(lagged) = make.unique(colnames(lagged))
  
  #1 corresponds to t-0, 10 c. to t-9; t-10 is first included predictor for horizon=10
  return(as_tibble(lagged)[, -(2:(horizon*ncol(series)))])
}

#predictors for forecasting y_t are y_{t-10} to y_{t-365}
fd_365lags = addPredictors(first_diffs, 365)
dates_365lags = addPredictors(water_levels[-1, 1], 365)


#given a time series (as a tibble),
#returns a list of tibbles, each of which is a slice of the original series.
##slice time series into what will later be the "outer" windows
#note that horizon should be 10 for ARIMA models, but 1 for Lasso models, as those are transformed beforehand with addPredictors()

timeSlices = function(series, slice_length, last_of_first_slice=slice_length, horizon=1) {
  
  #...passing last_of_first_slice instead of just first_t is necessary because we'll later have different window lengths
  #-> unbiased comparison only if tested on exactly the same observations (shorter windows shouldn't be tested on earlier obs)
  #...conveniently, the last_of_first_slice that'll be passed to all calls later is simply the longest of all slice_length
  #(the longest window starts with obs 1; shorter windows start later so that the actuals for OOS testing are always the same)
  first_t = 1 + last_of_first_slice - slice_length
  last_start = nrow(series)-(slice_length-1)
  slices = lapply(first_t:(last_start), function(x) series[x:(x+(slice_length-1)), ])
  
  #for Lasso models, 10-step-ahead is handled in addPredictors 
  #so horizon should be set to 1; then just drops last slice
  
  #drop all slices without a corresponding actual to compare forecast against
  #1 step: drop length = length:length = length+1-1:length
  #2 step: drop length-1:length = length+1-2:length
  #etc.
  first_dropped = (length(slices)+1-horizon)
  slices = slices[-(first_dropped:length(slices))]
  
  return(slices)
}


#the forecasted observation corresponding to a slice is a single row that contains as its first value the actual the forecast is compared against (t),
#followed by the predictors (t-10, t-11, ...)


forecastedObservations = function(series, slice_length, last_of_first_slice, horizon=1) {
  
  #call with same arguments as timeSlices (!)
  #returns a vector of actuals
  #i. e., with length l and horizon h, 
  #for every slice (t, ..., T) -> y_(T+h)
  ##values of observations to pseudo-OOS-test on. one for every slice (...same order) <-> one forecast from every outer window
  
  first_actual = last_of_first_slice + horizon
  #assumes the first column is y_t
  actuals = series[first_actual:nrow(series), ]
  return(actuals)
}

#given a single slice (=outer window), choose a lambda using "inner window" validation
lambdaFromSlice = function(slice, val_windowlength, use_lasso=TRUE, horiz=1, rolling_validation=TRUE) {
  
  #rows: last 10 rows contain information later than t-10
  #->drop last 10 rows
  first_dropped = nrow(slice)-9
  truncated_slice = slice[-(first_dropped:nrow(slice)), ]
  
  #columns: predictors are the lags (...only information from t-10 or earlier because 1st to 9th lag have been dropped before in addPredictors())
  #easiest to just check if there's a "." in the name -> just lagged variables
  x = truncated_slice[, (stringi::stri_count(names(truncated_slice), fixed = ".") == 1)]
  #actuals that validation forecasts are compared against - simply the first column
  #coercion to double is necessary because otherwise train() assumes a classification problem (...why?!)
  y = as.double(as.matrix(truncated_slice[, 1]))
  
  validation_control <- trainControl(method = "timeslice", initialWindow = val_windowlength, horizon = 1,
                                     fixedWindow = rolling_validation, skip = 0, allowParallel = TRUE)
  
  #lambda grid - adapted from Introduction to Statistical Learning
  grid =c(10^seq(10, 1, length=20), 10^seq(1,-2, length=20))
  
  #fit on "inner" window/slice of [val_windowlength] obs, test (forecast and compare against actual)
  #repeat/walk window forward one step at a time, until last obs was tested on
  #do this for every lambda in the grid
  validation_lasso_fits <- train(x, y, method = "glmnet",
                                 trControl=validation_control, metric="RMSE",
                                 tuneGrid = expand.grid(alpha = ifelse(use_lasso, 1, 0), lambda = grid))
  
  #lambda that yielded the lowest mean error across all validation forecasts?
  return(validation_lasso_fits$bestTune$lambda)
  
}


#given a single slice (=outer window), fit on whole "outer" window, and return 1-step-ahead forecast
forecastFromSlice = function(slice, forecasted_obs, lambda_opt, use_lasso=TRUE, horiz=1) {
  #as above:
  #rows: last 10 rows contain information later than t-10
  #->drop last 10 rows
  first_dropped = nrow(slice)-9
  truncated_slice = slice[-(first_dropped:nrow(slice)), ]
  
  #columns: predictors are the lags (...only information from t-10 or earlier because 1st to 9th lag have been dropped before in addPredictors())
  #easiest to just check if there's a "." in the name -> just lagged variables
  x = as.matrix(truncated_slice[, (stringi::stri_count(names(truncated_slice), fixed = ".") == 1)])
  #actuals that validation forecasts are compared against - simply the first column
  #coercion to double is necessary because otherwise train() assumes a classification problem (...why?!)
  #print("berfore y")
  y = as.double(as.matrix(truncated_slice[, 1]))
  #print("after y")
  #for "1-step-ahead" forecast (still forecasts t+1 with (t-9, t-10, ...)),
  #predictors come form the "forecasted observation"
  x_new = as.matrix(forecasted_obs[, -1])
  #print(x_new)
  
  #print(x)
  #print(y)
  #fit on outer window, i. e., y_t = f_{t-1}(y_t-1, y_t-2, ...), where f minimizes RSS+lambda*absolute sum of coefs
  OOS_test_fit = glmnet(x, y, alpha = ifelse(use_lasso, 1, 0), lambda = lambda_opt)
  #print(coef(OOS_test_fit))
  #with the coefficients from that fit, forecast 1 step ahead, i. e., y_t+1 = f_{t-1}(y_t, y_t-1, ...)
  return(predict(OOS_test_fit, newx=x_new))
  
}


#return a vector of forecasts, one for every window   
autoArimaForecasts = function(window_length, last_of_first_window, horiz=10) {
  years = timeSlices(water_levels[, 2], window_length, last_of_first_window, horiz)
  #actuals to compare against
  actuals = forecastedObservations(water_levels[, 2], window_length, last_of_first_window, horiz)
  #dates of actuals
  dates = forecastedObservations(water_levels[, 1], window_length, last_of_first_window, horiz)
  #fit an ARIMA model - with order chosen by auto.arima - on every slice
  arima_fits = pbsapply(years, function(x) auto.arima(x))
  #generate 10-step-ahead forecasts from every fit (...recursive forecasts)
  arima_forecasts = pbsapply(arima_fits, function(x) forecast(x, h=10)$mean)
  #prediction error = y_hat - y
  PEs = arima_forecasts[10, ] - actuals
  results = tibble(
    date = dates$date, 
    y = actuals$level, 
    y_hat = arima_forecasts[10, ],
    PE = PEs$level
  )
  return(results)
}


lassoForecasts = function(winlength, val_winlength, last_of_first_win)
{
  #unfortunately hardcoded
  slices_365 = timeSlices(fd_365lags, winlength, last_of_first_win)
  date_slices_365 = timeSlices(dates_365lags, winlength, last_of_first_win)
  fobs_365 = forecastedObservations(fd_365lags, winlength, last_of_first_win)
  fc_dates_365 = forecastedObservations(dates_365lags, winlength, last_of_first_win)
  
  lambdas_365 = pblapply(slices_365, function(slice) lambdaFromSlice(slice, val_winlength))
  
  forecasts_365 = pblapply(1:length(slices_365),
                           function(i) forecastFromSlice(slices_365[[i]], fobs_365[i, ],lambda_opt = lambdas_365[[i]]))
  
  forecasts_365 = t(as_tibble(forecasts_365, .name_repair = "minimal"))
  
  #level forecast = first difference forecast + actual level one timestep earlier 
  #... y_hat_t =  y_t-1 + (y_t - y_t-1)_hat
  first_date = match(fc_dates_365[1, 1], as.character(water_levels[["date"]])) - 1
  
  level_forecasts_365 = water_levels[first_date:(nrow(water_levels)-1), 2] + forecasts_365
  
  PEs = level_forecasts_365 - water_levels[(first_date+1):nrow(water_levels), 2]
  
  frame = data.frame("date"=fc_dates_365[, 1], "actual"=water_levels[(first_date+1):nrow(water_levels), 2],
                     "fcast"=level_forecasts_365, "PE"=PEs,
                     "winlength"=winlength, "val_winlength"=val_winlength)
  names(frame) = c("date", "actual", "fcast", "PE", "winlength", "val_winlength")
  
  return(frame)
}

#parallel processing
my_cluster <- makePSOCKcluster(3)
registerDoParallel(my_cluster)

set.seed(1)

#this way the first ARIMA forecast is for the same time as the first Lasso forecast
several_arima_forecasts = pblapply(183*(1:4), function(winlength) autoArimaForecasts(winlength, 732+367))


several_lasso_forecasts_val20 = pblapply(183*(1:4), function(winlength) lassoForecasts(winlength, winlength-20, 732))

several_lasso_forecasts_val40 = pblapply(183*(1:4), function(winlength) lassoForecasts(winlength, winlength-40, 732))

several_lasso_forecasts_val80 = pblapply(183*(1:4), function(winlength) lassoForecasts(winlength, winlength-80, 732))

several_lasso_forecasts_val160 = pblapply(183*(1:4), function(winlength) lassoForecasts(winlength, winlength-160, 732))

print(length(several_lasso_forecasts_val160[[1]]$PE))
print(length(several_arima_forecasts[[1]]$PE))


allPEs = cbind(pbsapply(several_arima_forecasts, function(result) (result$PE)),
               pbsapply(several_lasso_forecasts_val20, function(result) (result$PE)),
               pbsapply(several_lasso_forecasts_val40, function(result) (result$PE)),
               pbsapply(several_lasso_forecasts_val80, function(result) (result$PE)),
               pbsapply(several_lasso_forecasts_val160, function(result) (result$PE)))


#model confidence set procedure at the 90% confidence level (alpha = 0.1)
MCSprocedure(allPEs^2, 0.1, B = 5000, cl = my_cluster)
#yields all Lasso models as the superior set of models
#so I decided to simply average all Lasso models for the final forecasts

#"real" OOS forecasts:
#I copied the data and added 10 rows at the end with value 0 (see rawdata_for_predictions.txt)


#this is necessary because of the hardcoded reference in lassoForecasts
#...adding an argument in lassoForecasts() would have been better but I ran out of time...
water_levels = read_tsv("rawdata_for_predictions.txt")


first_diffs = firstDifferences(water_levels[, 2])

dates = water_levels[2:nrow(water_levels), 1]

fd_365lags = addPredictors(first_diffs, 365)
dates_365lags = addPredictors(water_levels[-1, 1], 365)


real_lasso_forecasts = cbind(pbsapply(183*(1:4), function(winlength) lassoForecasts(winlength, winlength-20, 1235)$fcast),
                             pbsapply(183*(1:4), function(winlength) lassoForecasts(winlength, winlength-40, 1235)$fcast),
                             pbsapply(183*(1:4), function(winlength) lassoForecasts(winlength, winlength-80, 1235)$fcast),
                             pbsapply(183*(1:4), function(winlength) lassoForecasts(winlength, winlength-160, 1235)$fcast))

stopCluster(cl=my_cluster)

#average at each time = average across rows
ensemble_forecasts = rowMeans(real_lasso_forecasts)

#forecasts are previous actual plus forecasted first difference (growth). 
#so for 9 of the 10 "real" OOS forecasts, where I set the actual as 0, the forecast returned will be just the forecasted first difference
#-> need to sum manually starting with the first "real" OOS forecast
final_level_forecasts = sapply((10:19), function(i) sum(ensemble_forecasts[10:i]))

final_forecasts_with_dates = cbind(dates[(nrow(dates)-9):nrow(dates), ], final_level_forecasts)

write.csv(final_forecasts_with_dates,'final_forecasts.csv')