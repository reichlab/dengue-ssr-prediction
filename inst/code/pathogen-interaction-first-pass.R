## making graphs of preliminary data
## R21 grant
## February 2014
## Nicholas Reich

library(scales)
library(reshape2)

setwd("~/Dropbox/work/proposals/TSIR_R21/code")
flu <- read.csv("../data/WHO_NREVSS.csv", 
                colClasses=c("factor", "factor", rep("numeric", 11)),
                na.strings="X")
rsv <- read.csv("../data/RSV_CDC_national_1989_2010.csv")


## add date to flu data and plot
flu$date <- as.Date(paste0(flu$YEAR, "-", flu$WEEK, "-1"), "%Y-%U-%u")
flu$A <- rowSums(as.matrix(flu[,c(7:11,13)]))
flu1 <- subset(flu, select=c("A", "B", "date"))
flu1$type <- "flu"
ggplot(flu1) + geom_linerange(aes(date, ymin=0, ymax=A), alpha=.5) + theme_bw()
#qplot(date, A, data=flu1, geom="line")

## plot RSV data
rsv$date2 <- as.Date(as.character(rsv$RepWeekDate), format="%m/%d/%y")
ggplot(rsv) + geom_linerange(aes(date2, ymin=0, ymax=RSVpos), alpha=.5) + theme_bw()


## plot together
ggplot(rsv) + 
        geom_linerange(aes(date2, ymin=0, ymax=RSVpos), alpha=.5) + 
        geom_linerange(data=flu1, aes(date, ymin=0, ymax=A), color="red", alpha=.5) + 
        geom_linerange(data=flu1, aes(date, ymin=0, ymax=B), color="blue", alpha=.5) + 
        theme_bw() + guides(color=FALSE)


rsv$color <- "black"
flu2 <- melt(flu1, id.vars=c("date"))
flu2$value <- as.numeric(flu2$value)

## plot together
ggplot(rsv) + 
        geom_line(aes(date2, RSVpos, color=color), alpha=.8, lwd=1) + 
        geom_line(data=flu2, aes(date, value, size=1, color=variable), alpha=.5, lwd=1) + 
        theme_bw() + theme(legend.position=c(.1,.7), legend.title=element_blank()) +  
        #guides(size=FALSE) + 
        #xlim(as.Date("1998-05-01"), as.Date("2014-01-01")) +
        ylab("case counts") + xlab("") +
        scale_x_date(breaks=seq.Date(as.Date("1998-01-01"), as.Date("2011-01-01"), "year"),
                     limits=c(as.Date("1998-01-01"), as.Date("2010-01-01")),
                     labels = date_format("%Y")) + 
        scale_color_manual(values=c("black", "red", "blue", "white"), 
                          labels=c("Influenza A", "Influenza B", "RSV", ""))


## reading in and plotting dengue data
load("/Users/nick/Dropbox/work/research/dengueCrossProtection/data/bkk.dengue.all.cases.rda")
den <- melt(bkk.dengue.all.cases, id.vars="date", measure.vars=c("den.cases.str1","den.cases.str2","den.cases.str3","den.cases.str4"))

ggplot(den) + 
        geom_line(aes(date, value, color=variable), alpha=.8) + 
        theme_bw() + theme(legend.position=c(.1,.6)) + 
        ylab("case counts") + xlab("") +
        scale_color_brewer(palette="Dark2", name="Dengue",
                           labels=c("Serotype 1", "Serotype 2", "Serotype 3", "Serotype 4")) +
        scale_x_date(limits=c(as.Date("1974-08-01"), as.Date("2009-06-01")))
        

## smooth dengue data
for(i in 1:4){
       
}


## manifold plot of flu data
flu_wide <- spread(flu2[complete.cases(flu2),], key=variable, value=value)

## smooth of flu A
smooth_A <- loess(A ~ as.numeric(date), span=1/80, data=flu_wide)
flu_wide$sm_A <- smooth_A$fitted
flu_wide[flu_wide$sm_A<0, "sm_A"] <- 0

## smooth of flu B
smooth_B <- loess(B ~ as.numeric(date), span=1/80, data=flu_wide)
flu_wide$sm_B <- smooth_B$fitted
flu_wide[flu_wide$sm_B<0, "sm_B"] <- 0

## plots of A and B
qplot(x=date, y=A, data=flu_wide, geom="line")
qplot(x=date, y=B, data=flu_wide, geom="line")
qplot(x=date, y=sm_A, data=flu_wide, geom="line")
qplot(x=date, y=sm_B, data=flu_wide, geom="line")

## manifold plot
ggplot(flu_wide) + geom_path(aes(x=sm_A, y=sm_B)) + scale_x_log10() + scale_y_log10() +
    geom_abline(a=0, b=1, linetype=2, color="gray")

## lags for A
flu_wide$sm_A_lagx <- lag(flu_wide$sm_A, 52)
ggplot(flu_wide) + geom_path(aes(x=sm_A, y=sm_A_lagx)) + scale_x_log10() + scale_y_log10() +
    geom_abline(a=0, b=1, linetype=2, color="gray")

flu_wide$sm_B_lagx <- lag(flu_wide$sm_B, 52)
ggplot(flu_wide) + geom_path(aes(x=sm_B, y=sm_B_lagx)) + scale_x_log10() + scale_y_log10() +
    geom_abline(a=0, b=1, linetype=2, color="gray")


