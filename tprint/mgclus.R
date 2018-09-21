mgclusin2pajek <- function(mgclusin) {
        
        prefix <- unlist(strsplit(mgclusin, "\\.(?=[^.]+$)", perl=TRUE))[1]
        outfile <- paste(prefix, '.net', sep='')
        
        df <- read.table(mgclusin, sep = "\t",
                         stringsAsFactors = FALSE)
        
        from <- df[, 1]
        to <- df[, 2]
        
        node_name <- union(from, to)
        nnode <- length(node_name)
        node_list <- data.frame(node_index=1:nnode, node_name=node_name)
        
        nedge <- dim(df)[1]
        from <- match(from, node_name)
        to <- match(to, node_name)
        if (dim(df)[2] > 2) {
                edge_list <- cbind(from, to, df[, 3])
        } else {
                edge_list <- cbind(from, to, 1)
        }

        write(paste("*Vertices", nnode, sep=" " ), outfile)
        write.table(node_list, outfile, row.names = F, col.names = F, append = T)
        write(paste("*Edges", nedge, sep=" " ), outfile, append = T)
        write.table(edge_list, outfile, row.names = F, col.names = F, append = T)
        
}

pajek2mgclusin <- function(pajek) {
        
        prefix <- unlist(strsplit(pajek, "\\.(?=[^.]+$)", perl=TRUE))[1]
        outfile <- paste(prefix, '.tsv', sep='')
        
        con <- file(pajek, "r")
        first_line <- readLines(con, n=1)
        close(con)
        nnode <- as.numeric(unlist(strsplit(first_line, " "))[2])
        
        node_list <- read.table(pajek, sep = " ", nrows = nnode, skip = 1,
                           stringsAsFactors = FALSE)
        node_name <- node_list[, 2]
        
        edge_list <- read.table(pajek, sep = " ", skip = nnode + 2,
                                stringsAsFactors = FALSE)
        from <- node_name[edge_list[, 1]]
        to <- node_name[edge_list[, 2]]
        if (dim(edge_list)[2] > 2) {
                weight <- edge_list[, 3]
        } else {
                weight <- 1
        }
        
        df <- data.frame(from, to, weight)
        write.table(df, outfile, quote = F, sep = "\t", row.names = F, col.names = F)
        
}

mgclusout2cluster <- function(pajek, mgclusout) {
        
        prefix <- unlist(strsplit(pajek, "\\.(?=[^.]+$)", perl=TRUE))[1]
        outfile <- paste(prefix, '_mgclus_clu.txt', sep='')
        
        con <- file(pajek, "r")
        first_line <- readLines(con, n=1)
        close(con)
        nnode <- as.numeric(unlist(strsplit(first_line, " "))[2])
        
        node_list <- read.table(pajek, sep = " ", nrows = nnode, skip = 1,
                                stringsAsFactors = FALSE)
        node_name <- node_list[, 2]
        
        df <- read.table(mgclusout, sep = "\t", fill = TRUE)
        cluster <- sapply(node_name, function(x, dataframe) which(dataframe==x,
                                                                  arr.ind=T)[1], df)
        
        write(as.character(cluster), outfile)
        
}








