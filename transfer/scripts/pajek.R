write_pajek_pthres <- function(df, outfile, minsize = 2, threshold = 0.05) {
        
        # only consider links with at least minsize samplesize
        df_filtered <- df[df$samplesize >= minsize, ]
        
        # adjust pvalue
        pval_adj <- p.adjust(df_filtered$pvalue, method = "BH")
        select <- pval_adj < threshold
        
        from <- df_filtered$gene_1[select]
        to <- df_filtered$gene_2[select]
        
        node_name <- union(from, to)
        nnode <- length(node_name)
        node_list <- data.frame(node_index=1:nnode, node_name=node_name)
        
        nedge <- sum(select)
        from <- match(from, node_name)
        to <- match(to, node_name)
        edge_list <- cbind(from, to, 1)
        
        write(paste("*Vertices", nnode, sep=" " ), outfile)
        write.table(node_list, outfile, row.names = F, col.names = F, append = T)
        write(paste("*Edges", nedge, sep=" " ), outfile, append = T)
        write.table(edge_list, outfile, row.names = F, col.names = F, append = T)
        
}

write_clu <- function(infomap.clu) {
        
        prefix <- unlist(strsplit(infomap.clu, "\\.(?=[^.]+$)", perl=TRUE))[1]
        clu.txt <- paste(prefix, '_infomap_clu.txt', sep='')
        
        df <- read.table(infomap.clu, sep = " ", skip = 2,
                         stringsAsFactors = FALSE)
        order <- order(df[ , 1])
        clu_mem <- df[order, 2]
        
        write(as.character(clu_mem), clu.txt)

}

extract_module <- function(pajek.net, clu.txt, cluster_no) {
        
        con <- file(pajek.net, "r")
        first_line <- readLines(con, n=1)
        close(con)
        nnode <- as.numeric(unlist(strsplit(first_line, " "))[2])
        
        df_1 <- read.table(pajek.net, sep = " ", nrows = nnode, skip = 1,
                           stringsAsFactors = FALSE)
        gene <- df_1[, 2]
        
        df_2 <- read.table(clu.txt,
                           stringsAsFactors = FALSE)
        clu_mem <- df_2[ , 1]
        
        module_gene <- gene[clu_mem == cluster_no]
        smodule_gene <- gsub(" ENSG.*", "", module_gene)
        paste(smodule_gene, collapse = ',')
}

#extract_module("input/sample1.2_cancer.net", "output/sample1.2_cancer_clu.txt", 1)

extract_original <- function(pajek.net, clu.txt, cluster_no) {
        
        con <- file(pajek.net, "r")
        first_line <- readLines(con, n=1)
        close(con)
        nnode <- as.numeric(unlist(strsplit(first_line, " "))[2])
        
        df_1 <- read.table(pajek.net, sep = " ", nrows = nnode, skip = 1,
                           stringsAsFactors = FALSE)
        gene <- df_1[, 2]
        
        df_2 <- read.table(clu.txt,
                           stringsAsFactors = FALSE)
        clu_mem <- df_2[ , 1]
        
        gene[clu_mem == cluster_no]
}

#extract_original("input/sample1.2_cancer.net", "output/sample1.2_cancer_clu.txt", 1)

#./Infomap input/file output/ -N 100 --clu


write_pajek_pthres(df_epithelial, "input/sample1.2_epithelial_2_005.net", minsize = 2, threshold = 0.05)
#./Infomap input/sample1.2_epithelial_2_005.net output/ -N 100 --clu
write_clu("output/sample1.2_epithelial_2_005.clu")
clu_mem <- read.table("output/sample1.2_epithelial_2_005_clu.txt",
                      stringsAsFactors = FALSE)[ , 1]
table(clu_mem)

write_pajek_pthres(df_epithelialexp, "input/sample1.2_epithelialexp_2_005.net", minsize = 2, threshold = 0.05)
#./Infomap input/sample1.2_epithelialexp_2_005.net output/ -N 100 --clu
write_clu("output/sample1.2_epithelialexp_2_005.clu")
clu_mem <- read.table("output/sample1.2_epithelialexp_2_005_clu.txt",
                      stringsAsFactors = FALSE)[ , 1]
table(clu_mem)

write_pajek_pthres(df_cancer, "input/sample1.2_cancer_2_005.net", minsize = 2, threshold = 0.05)
#./Infomap input/sample1.2_cancer_2_005.net output/ -N 100 --clu
write_clu("output/sample1.2_cancer_2_005.clu")
clu_mem <- read.table("output/sample1.2_cancer_2_005_clu.txt",
                      stringsAsFactors = FALSE)[ , 1]
table(clu_mem)

write_pajek_pthres(df_epican, "input/sample1.2_epican_2_005.net", minsize = 2, threshold = 0.05)
#./Infomap input/sample1.2_epican_2_005.net output/ -N 100 --clu
write_clu("output/sample1.2_epican_2_005.clu")
clu_mem <- read.table("output/sample1.2_epican_2_005_clu.txt",
                      stringsAsFactors = FALSE)[ , 1]
table(clu_mem)

write_pajek_pthres(df_epicanexp, "input/sample1.2_epicanexp_2_005.net", minsize = 2, threshold = 0.05)
#./Infomap input/sample1.2_epicanexp_2_005.net output/ -N 100 --clu
write_clu("output/sample1.2_epicanexp_2_005.clu")
clu_mem <- read.table("output/sample1.2_epicanexp_2_005_clu.txt",
                      stringsAsFactors = FALSE)[ , 1]
table(clu_mem)


write_pajek_pthres(df_epithelialexp, "input/sample1.2_epithelialexp_005.net", threshold = 0.05)
write_pajek_pthres(df_cancer, "input/sample1.2_cancer_005.net", threshold = 0.05)
write_pajek_pthres(df_epicanexp, "input/sample1.2_epicanexp_005.net", threshold = 0.05)

write_pajek_pthres(df_epithelialexp, "input/sample1.2_epithelialexp_010.net", threshold = 0.10)
write_pajek_pthres(df_cancer, "input/sample1.2_cancer_010.net", threshold = 0.10)
write_pajek_pthres(df_epicanexp, "input/sample1.2_epicanexp_010.net", threshold = 0.10)

write_pajek_pthres(df_epithelialexp, "input/sample1.2_epithelialexp_015.net", threshold = 0.15)
write_pajek_pthres(df_cancer, "input/sample1.2_cancer_015.net", threshold = 0.15)
write_pajek_pthres(df_epicanexp, "input/sample1.2_epicanexp_015.net", threshold = 0.15)

write_pajek_pthres(df_epithelialexp, "input/sample1.2_epithelialexp_020.net", threshold = 0.20)
write_pajek_pthres(df_cancer, "input/sample1.2_cancer_020.net", threshold = 0.20)
write_pajek_pthres(df_epicanexp, "input/sample1.2_epicanexp_020.net", threshold = 0.20)

#./Infomap input/sample1.2_epithelialexp_005.net output/ -N 100 --clu
#./Infomap input/sample1.2_cancer_005.net output/ -N 100 --clu
#./Infomap input/sample1.2_epicanexp_005.net output/ -N 100 --clu

#./Infomap input/sample1.2_epithelialexp_010.net output/ -N 100 --clu
#./Infomap input/sample1.2_cancer_010.net output/ -N 100 --clu
#./Infomap input/sample1.2_epicanexp_010.net output/ -N 100 --clu

#./Infomap input/sample1.2_epithelialexp_015.net output/ -N 100 --clu
#./Infomap input/sample1.2_cancer_015.net output/ -N 100 --clu
#./Infomap input/sample1.2_epicanexp_015.net output/ -N 100 --clu

#./Infomap input/sample1.2_epithelialexp_020.net output/ -N 100 --clu
#./Infomap input/sample1.2_cancer_020.net output/ -N 100 --clu
#./Infomap input/sample1.2_epicanexp_020.net output/ -N 100 --clu

write_clu("output/sample1.2_epithelialexp_005.clu")
write_clu("output/sample1.2_cancer_005.clu")
write_clu("output/sample1.2_epicanexp_005.clu")

write_clu("output/sample1.2_epithelialexp_010.clu")
write_clu("output/sample1.2_cancer_010.clu")
write_clu("output/sample1.2_epicanexp_010.clu")

write_clu("output/sample1.2_epithelialexp_015.clu")
write_clu("output/sample1.2_cancer_015.clu")
write_clu("output/sample1.2_epicanexp_015.clu")

write_clu("output/sample1.2_epithelialexp_020.clu")
write_clu("output/sample1.2_cancer_020.clu")
write_clu("output/sample1.2_epicanexp_020.clu")
