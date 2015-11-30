library(ssr)
library(ggplot2)
library(grid)


last_train_ind <- which(San_Juan_test$season == "2008/2009" & San_Juan_test$season_week == 52)

train_label_x <- as.numeric(San_Juan_test$week_start_date[1]) + 0.8 * (as.numeric(San_Juan_test$week_start_date[last_train_ind]) - as.numeric(San_Juan_test$week_start_date[1]))

p <- ggplot() +
    geom_line(aes(x=week_start_date, y = total_cases), data=San_Juan_test) +
    geom_vline(xintercept = as.numeric(San_Juan_test$week_start_date[last_train_ind]), colour = "#E69F00") +
#    geom_text(data = NULL, label = "Train", x = San_Juan_test$week_start_date[floor(0.8 * last_train_ind)], y = 400) +
    geom_text(data = NULL, label = "Train", x = 0.8, y = 0.5) +
    xlab("Time") +
    ylab("Total Cases") +
    ggtitle("Dengue Fever Incidence in San Juan, PR") +
    theme_bw(base_size=22)


pdf("/media/evan/data/Reich/dengue-ssr-prediction/inst/intermediate-results/ssr-poster/plots/initial-data-plot.pdf", width=10, height=10)
print(p)
grid.text("Train Data", x = unit(0.56, "npc"), y = unit(0.91, "npc"), gp = gpar(fontsize = 20))
grid.text("Test Data", x = unit(0.895, "npc"), y = unit(0.91, "npc"), gp = gpar(fontsize = 20))
dev.off()

