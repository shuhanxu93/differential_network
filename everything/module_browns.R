row_logfold <- function(x, size1, size2) {
        
        control <- x[1:size1]
        perturb <- x[(size1+1):(size1+size2)]
        
        log10(median(perturb)/median(control))
        
}

row_wilcox <- function(x, size1, size2) {
        
        control <- x[1:size1]
        perturb <- x[(size1+1):(size1+size2)]

        wilcox.test(control, perturb)$p.value

}

stat_test <- function(st_data, control_coordinate, perturb_coordinate) {
        df_st <- read.csv(st_data, sep = "\t")
        mat_st <- as.matrix(df_st)
        
        all_coordinate <- names(df_st)
        
        select_control <- match(control_coordinate, all_coordinate)
        select_control <- select_control[!is.na(select_control)]
        control_size <- length(select_control)
        
        select_perturb <- match(perturb_coordinate, all_coordinate)
        select_perturb <- select_perturb[!is.na(select_perturb)]
        perturb_size <- length(select_perturb)
        
        mat_control <- mat_st[, select_control]
        mat_perturb <- mat_st[, select_perturb]
        mat_both <- cbind(mat_control, mat_perturb)
        mat_both <- round(mat_both, 1) - 1
        log_fold <- apply(mat_both, 1, row_logfold, size1 = control_size, size2 = perturb_size)
        p_values <- apply(mat_both, 1, row_wilcox, size1 = control_size, size2 = perturb_size)
        list(m = mat_both, l = log_fold, p = p_values)
}

lnames = load(file='coordinate.RData')
lnames

sample1.2_test <- stat_test("sample1.2_deconvolution.csv", sample1.2_epithelial_expanded, sample1.2_cancer)
sample1.2_pvalues <- sample1.2_test$p
sample1.2_pvalues_adj <- p.adjust(sample1.2_pvalues, method = "BH")


module_pvalue <- function(pajek.net, clu.txt, data_matrix, p_values) {
        
        #open library for Empirical Browns Method
        library(EmpiricalBrownsMethod)
        
        #extract information from pajek.net and clu.txt files
        con <- file(pajek.net, "r")
        first_line <- readLines(con, n=1)
        close(con)
        nnode <- as.numeric(unlist(strsplit(first_line, " "))[2])
        
        df_1 <- read.table(pajek.net, sep = " ", nrows = nnode, skip = 1)
        gene <- df_1[, 2]
        
        df_2 <- read.table(clu.txt)
        clu_mem <- df_2[ , 1]
        
        all_gene <- rownames(data_matrix)
        
        #create 3 vectors
        cluster <- 1:max(clu_mem)
        cluster_size <- numeric()
        browns <- numeric()
        
        #for each cluster_index in cluster
        for (cluster_index in cluster) {
                boo <- clu_mem == cluster_index
                cluster_size <- c(cluster_size, sum(boo))
                cluster_gene <- gene[boo]
                index <- match(cluster_gene, all_gene)
                browns <- c(browns, empiricalBrownsMethod(data_matrix[index, ], p_values[index]))

        }
        
        summary <- data.frame(cluster = cluster, cluster_size = cluster_size, browns = browns)
        summary
}

epithelial_modules <- module_pvalue("input/sample1.2_epithelialexp_2_005.net", "output/sample1.2_epithelialexp_2_005_clu.txt",
                                    sample1.2_test$m, sample1.2_pvalues_adj)
cancer_modules <- module_pvalue("input/sample1.2_cancer_2_005.net", "output/sample1.2_cancer_2_005_clu.txt",
                                sample1.2_test$m, sample1.2_pvalues_adj)
epican_modules <- module_pvalue("input/sample1.2_epicanexp_2_005.net", "output/sample1.2_epicanexp_2_005_clu.txt",
                                sample1.2_test$m, sample1.2_pvalues_adj)

gene <- rownames(sample1.2_test$m)
log_fold <- sample1.2_test$l
wilcox_padjust <- sample1.2_pvalues_adj

con <- file("input/sample1.2_epithelialexp_2_005.net", "r")
first_line <- readLines(con, n=1)
close(con)
nnode <- as.numeric(unlist(strsplit(first_line, " "))[2])
df_1 <- read.table("input/sample1.2_epithelialexp_2_005.net",
                   sep = " ", nrows = nnode, skip = 1)
cluster_gene <- df_1[, 2]
df_2 <- read.table("output/sample1.2_epithelialexp_2_005_clu.txt")
clu_mem <- df_2[ , 1]
epithelial_cluster <- clu_mem[match(gene, cluster_gene)]
epithelial_browns <- epithelial_modules[epithelial_cluster, 3]

con <- file("input/sample1.2_cancer_2_005.net", "r")
first_line <- readLines(con, n=1)
close(con)
nnode <- as.numeric(unlist(strsplit(first_line, " "))[2])
df_1 <- read.table("input/sample1.2_cancer_2_005.net",
                   sep = " ", nrows = nnode, skip = 1)
cluster_gene <- df_1[, 2]
df_2 <- read.table("output/sample1.2_cancer_2_005_clu.txt")
clu_mem <- df_2[ , 1]
cancer_cluster <- clu_mem[match(gene, cluster_gene)]
cancer_browns <- cancer_modules[cancer_cluster, 3]

con <- file("input/sample1.2_epicanexp_2_005.net", "r")
first_line <- readLines(con, n=1)
close(con)
nnode <- as.numeric(unlist(strsplit(first_line, " "))[2])
df_1 <- read.table("input/sample1.2_epicanexp_2_005.net",
                   sep = " ", nrows = nnode, skip = 1)
cluster_gene <- df_1[, 2]
df_2 <- read.table("output/sample1.2_epicanexp_2_005_clu.txt")
clu_mem <- df_2[ , 1]
both_cluster <- clu_mem[match(gene, cluster_gene)]
both_browns <- epican_modules[both_cluster, 3]

gene_summary <- data.frame(gene=gene, log_fold=log_fold, wilcox_padjust=wilcox_padjust,
                           epithelial_cluster=epithelial_cluster, epithelial_browns=epithelial_browns,
                           cancer_cluster=cancer_cluster, cancer_browns=cancer_browns,
                           both_cluster=both_cluster, both_browns=both_browns)

write.table(gene_summary, "sample1.2_genesummary.csv", sep = '\t', row.names = F, col.names = T)
gene_summary <- read.table("sample1.2_genesummary.csv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)

