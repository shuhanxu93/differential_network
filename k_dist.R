lcc <- function(my_graph) {
        group <- which.max(components(my_graph)$csize)
        ids <- which(components(my_graph)$membership == group)
        my_sub <- induced_subgraph(my_graph, ids)
        my_sub
}

library("igraph")
library("poweRlaw")
library("ggplot2")

g_epi <- read.graph("sample1.2_epithelialexp_005s.net", format = "pajek")
g_cancer <- read.graph("sample1.2_cancer_005s.net", format = "pajek")
g_both <- read.graph("sample1.2_epicanexp_005s.net", format = "pajek")

g_epi_lcc <- lcc(g_epi)
g_cancer_lcc <- lcc(g_cancer)
g_both_lcc <- lcc(g_both)

G.degrees <- degree(g_epi)

# Let's count the frequencies of each degree
G.degree.histogram <- as.data.frame(table(G.degrees))

# Need to convert the first column to numbers, otherwise
# the log-log thing will not work (that's fair...)
G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])

fit=lm(log10(G.degree.histogram[,2]) ~ log10(G.degree.histogram[,1]))

# Now, plot it!
ggplot(G.degree.histogram, aes(x = G.degrees, y = Freq)) +
        geom_point() +
        geom_abline(intercept = coef(fit)[1], slope = coef(fit)[2], col='red') +
        scale_x_continuous("Degree\n(nodes with this amount of connections)",
                           breaks = c(1, 3, 10, 30, 100, 300),
                           trans = "log10") +
        scale_y_continuous("Frequency\n(how many of them)",
                           breaks = c(1, 3, 10, 30, 100, 300, 1000),
                           trans = "log10") +
        ggtitle("Degree Distribution (log-log), Healthy Region FDR threshold -") +
        theme_bw()

