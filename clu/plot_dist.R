clu_mem <- read.table("clu/sample1.2_epicanexp_020s_multilevel_clu.txt",
                      stringsAsFactors = FALSE)[ , 1]
sizes <- as.vector(table(clu_mem))
table(sizes)
x <- 1:max(sizes) + 0.5
hist(sizes, breaks = x)

