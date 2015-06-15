## read in and plot data

library(lubridate)
library(ggplot2)
library(dplyr)

sj <- read.csv('data/San_Juan_Training_Data.csv')

## add lag columns
sj <- sj %>% mutate(total_cases_lag3 = lag(total_cases, 3),
                    total_cases_lag1 = lag(total_cases),
                    week_start_date = ymd(week_start_date))

## simple time series plot
ggplot(sj, aes(x=week_start_date, y=total_cases)) + geom_line()
ggplot(sj, aes(x=week_start_date)) + geom_linerange(aes(ymin=0, ymax=total_cases))

## phase plots
(p <- ggplot(sj, aes(x=total_cases_lag1, y=total_cases, color=season)) + 
    geom_point() + geom_path() +
    scale_x_log10() + scale_y_log10() +
    geom_abline(b=1, a=0, color="grey"))

p + facet_wrap(~ season)

## smooth data
sm <- loess(total_cases ~ as.numeric(week_start_date), data=sj, span=1/80)
sj$smooth_cases <- sm$fitted
sj$smooth_cases_lag1 <- lag(sj$smooth_cases)

ggplot(sj, aes(x=week_start_date)) + 
    geom_linerange(aes(ymin=0, ymax=total_cases)) +
    geom_line(aes(y=smooth_cases), color="red", size=1)

## phase plot with smooth data
(p1 <- ggplot(sj, aes(x=smooth_cases_lag1, y=smooth_cases, color=season)) + 
    geom_point() + geom_path() +
    scale_x_log10() + scale_y_log10() +
    geom_abline(b=1, a=0, color="grey"))

p1 + facet_wrap(~ season)
