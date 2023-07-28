library(ggplot2)
library(reshape2)

nCycle <- 3
data_cors <- data.frame()

# correlation line plot
for (cycle in 1:nCycle) {
    if (cycle == 1) {
        data_cors <- read.csv(paste("sim_2023Jul18_101500/1C", cycle, "_svm_random_cors_snp_yield.csv", sep = ""))
        colnames(data_cors)[cycle + 1] <- paste("1C", cycle, sep = "")
    }else {
        temp_data <- read.csv(paste("sim_2023Jul18_101500/1C", cycle, "_svm_random_cors_snp_yield.csv", sep = ""))
        data_cors <- cbind(data_cors, temp_data$X1)
        colnames(data_cors)[cycle + 1] <- paste("1C", cycle, sep = "")
    }
}

data_cors$X <- factor(data_cors$X, levels = data_cors$X)
data_cors <- melt(data_cors, id.vars = "X", variable.name = "trials")

plot_cors <- ggplot(data_cors, aes(X, value, group = trials)) +
            geom_line(aes(colour = trials)) +
            geom_point(aes(colour = trials)) +
            labs(title = "Model Performance", x = "Generation", y = "Performance")
ggsave("CorrelationLinePlot.pdf", plot = plot_cors)
