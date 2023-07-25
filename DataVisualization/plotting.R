library(ggplot2)

# correlation line plot
data_cors <- read.csv("sim_2023Jul18_101500/1C1_svm_random_cors_snp_yield.csv")
data_cors$X <- factor(data_cors$X, levels = data_cors$X)

plot_cors <- ggplot(data_cors, aes(x = X, y = X1, group = 1)) +
            geom_line() +
            geom_point()
ggsave("CorrelationLinePlot.pdf", plot = plot_cors)
