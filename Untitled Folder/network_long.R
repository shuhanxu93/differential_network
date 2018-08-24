cor1 <- function(x, data_matrix) {
        spearman.test(data_matrix[x[1], ], data_matrix[x[2], ])$p.value
}

buildnetwork <- function(filename) {
        
        sample <- unlist(strsplit(filename, '_'))[1]
        
        df <- read.csv(filename, sep = '\t')
        gene <- rownames(df)
        gene_size <- length(gene)
        indx <- apply(df, 2, function(x) any(is.na(x) | is.infinite(x)))
        df <- df[, !indx]
        mat <- t(as.matrix(df))
        mat[mat < 1] <- NA
        
        boo <- lower.tri(meshgrid(1:gene_size)[[1]])
        gene_1 <- meshgrid(1:gene_size)[[1]][boo]
        gene_2 <- meshgrid(1:gene_size)[[2]][boo]
        
        results <- rcorr(mat, type = 'spearman')
        spearman <- results$r[boo]
        samplesize <- results$n[boo]
        pvalue <- results$P[boo]
        
        gene_combine <- cbind(gene_1, gene_2)
        boo <- samplesize > 1
        gene_boo <- gene_combine[boo, ]
        mat <- t(mat)
        detach("package:Hmisc", unload=TRUE)
        library("pspearman", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
        pvalue[boo] <- apply(gene_boo, 1, cor1, mat)
        
        pairwise <- data.frame(gene_1=gene[gene_1], gene_2=gene[gene_2], spearman=spearman,
                               samplesize=samplesize, pvalue=pvalue)
        
        write.table(pairwise, file = paste(sample, '_pairwise_exact22.csv', sep=''), 
                    sep = '\t', row.names = FALSE, col.names = TRUE)
}

buildsubnetwork <- function(filename, coordinate, subarea_name) {
        
        sample <- unlist(strsplit(filename, '_'))[1]
        sample <- paste(sample, subarea_name, sep='_')
        
        
        df <- read.csv(filename, sep = '\t')
        gene <- rownames(df)
        gene_size <- length(gene)
        indx <- apply(df, 2, function(x) any(is.na(x) | is.infinite(x)))
        df <- df[, !indx]
        all_coordinate <- names(df)
        select <- match(coordinate, all_coordinate)
        select <- select[!is.na(select)]
        print(length(select))
        mat <- t(as.matrix(df))
        mat <- mat[select, ]
        mat[mat < 1] <- NA
        
        boo <- lower.tri(meshgrid(1:gene_size)[[1]])
        gene_1 <- meshgrid(1:gene_size)[[1]][boo]
        gene_2 <- meshgrid(1:gene_size)[[2]][boo]
        
        results <- rcorr(mat, type = 'spearman')
        spearman <- results$r[boo]
        samplesize <- results$n[boo]
        pvalue <- results$P[boo]
        
        gene_combine <- cbind(gene_1, gene_2)
        boo <- samplesize > 1
        gene_boo <- gene_combine[boo, ]
        mat <- t(mat)
        detach("package:Hmisc", unload=TRUE)
        library("pspearman", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
        pvalue[boo] <- apply(gene_boo, 1, cor1, mat)
        
        pairwise <- data.frame(gene_1=gene[gene_1], gene_2=gene[gene_2], spearman=spearman,
                               samplesize=samplesize, pvalue=pvalue)

        write.table(pairwise, file = paste(sample, '_pairwise_exact22.csv', sep=''), 
                    sep = '\t', row.names = FALSE, col.names = TRUE)
}

setwd('/scratch/diffnet/shuhan/data')
library("pracma", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
library("Hmisc", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.4")
lnames = load(file='coordinate.RData')
lnames

#buildnetwork('sample1.2_deconvolution.csv')
#buildsubnetwork('sample1.2_deconvolution.csv', sample1.2_epithelial, "epithelial")
#buildsubnetwork('sample1.2_deconvolution.csv', sample1.2_epithelial_expanded, "epithelialexp")
#buildsubnetwork('sample1.2_deconvolution.csv', sample1.2_cancer, "cancer")
#combined = c(sample1.2_epithelial, sample1.2_cancer)
#buildsubnetwork('sample1.2_deconvolution.csv', combined, "epithelial_cancer")
combined = c(sample1.2_epithelial_expanded, sample1.2_cancer)
buildsubnetwork('sample1.2_deconvolution.csv', combined, "epithelialexp_cancer")