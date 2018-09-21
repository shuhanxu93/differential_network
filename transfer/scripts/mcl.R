mcl_clu <- function(pajek.net) {
        
        prefix <- unlist(strsplit(pajek.net, "\\.(?=[^.]+$)", perl=TRUE))[1]
        clu.txt <- paste(prefix, '_mcl_clu.txt', sep='')
        
        my_graph <- read_graph(pajek.net, format = "pajek")
        adj <- get.adjacency(my_graph)
        mcl_comm <- mcl(adj, addLoops = TRUE)
        write(as.character(mcl_comm$Cluster), clu.txt)
}

library("igraph")
library("MCL")

mcl_clu("sample1.2_epithelialexp_005.net")
mcl_clu("sample1.2_cancer_005.net")
mcl_clu("sample1.2_epicanexp_005.net")

mcl_clu("sample1.2_epithelialexp_010.net")
mcl_clu("sample1.2_cancer_010.net")
mcl_clu("sample1.2_epicanexp_010.net")

mcl_clu("sample1.2_epithelialexp_015.net")
mcl_clu("sample1.2_cancer_015.net")
mcl_clu("sample1.2_epicanexp_015.net")

mcl_clu("sample1.2_epithelialexp_020.net")
mcl_clu("sample1.2_cancer_020.net")
mcl_clu("sample1.2_epicanexp_020.net")