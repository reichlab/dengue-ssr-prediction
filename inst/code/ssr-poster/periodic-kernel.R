library(ggplot2)
library(reshape2)

plot_df <- data.frame(t=seq_len(5 * 52))

kernel_center <- 130
rho <- pi / 52

h <- 0.1
plot_df$kernel_h0.1 <- exp( -0.5 * (sin(rho * (kernel_center - plot_df$t)) / h)^2)

h <- 1
plot_df$kernel_h1 <- exp( -0.5 * (sin(rho * (kernel_center - plot_df$t)) / h)^2)

h <- 10
plot_df$kernel_h10 <- exp( -0.5 * (sin(rho * (kernel_center - plot_df$t)) / h)^2)

plot_df <- melt(plot_df, id.vars = "t")
plot_df$variable <- as.character(plot_df$variable)
plot_df$bandwidth <- "0.1"
plot_df$bandwidth[plot_df$variable == "kernel_h1"] <- "1"
plot_df$bandwidth[plot_df$variable == "kernel_h10"] <- "10"

p <- ggplot(plot_df) +
    geom_line(aes(x = t, y = value, linetype = bandwidth, colour = bandwidth)) +
    geom_vline(xintercept = kernel_center) +
    scale_colour_manual("Bandwidth",
        breaks = c("0.1", "1", "10"),
        labels = c("0.1", "1", "10"),
        values = c("#E69F00", "#56B4E9", "#009E73")
    ) +
    scale_linetype("Bandwidth") +
    ylab("Kernel Function Value") +
    xlab("Week") +
    ggtitle("The Periodic Kernel") +
    theme_bw(base_size = 22)

pdf("/media/evan/data/Reich/dengue-ssr-prediction/inst/intermediate-results/ssr-poster/plots/periodic-kernel.pdf", width = 10, height = 4)
print(p)
dev.off()


