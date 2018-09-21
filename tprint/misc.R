lcc <- function(pajek.net) {
        my_graph <- read_graph(pajek.net, format = "pajek")
        group <- which.max(components(my_graph)$csize)
        ids <- which(components(my_graph)$membership == group)
        my_sub <- induced_subgraph(my_graph, ids)
        c(gorder(my_sub), gsize(my_sub))
}

library("igraph", lib.loc="/Library/Frameworks/R.framework/Versions/3.5/Resources/library")
lcc("sample1.2_epicanexp_2_005s.net")