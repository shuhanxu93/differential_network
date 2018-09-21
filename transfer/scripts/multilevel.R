multilevel_clu <- function(pajek.net) {
        
        prefix <- unlist(strsplit(pajek.net, "\\.(?=[^.]+$)", perl=TRUE))[1]
        clu.txt <- paste(prefix, '_multilevel_clu.txt', sep='')
        
        my_graph <- read_graph(pajek.net, format = "pajek")
        multilevel_comm <- cluster_louvain(my_graph)
        write(as.character(membership(multilevel_comm)), clu.txt)
        

}

library("igraph")

multilevel_clu("sample1.2_epithelialexp_005.net")
multilevel_clu("sample1.2_cancer_005.net")
multilevel_clu("sample1.2_epicanexp_005.net")

multilevel_clu("sample1.2_epithelialexp_010.net")
multilevel_clu("sample1.2_cancer_010.net")
multilevel_clu("sample1.2_epicanexp_010.net")

multilevel_clu("sample1.2_epithelialexp_015.net")
multilevel_clu("sample1.2_cancer_015.net")
multilevel_clu("sample1.2_epicanexp_015.net")

multilevel_clu("sample1.2_epithelialexp_020.net")
multilevel_clu("sample1.2_cancer_020.net")
multilevel_clu("sample1.2_epicanexp_020.net")