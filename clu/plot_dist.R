clu_mem <- read.table("sample1.2_epicanexp_020s_multilevel_clu.txt",
                      stringsAsFactors = FALSE)[ , 1]
sizes <- as.vector(table(clu_mem))
table(sizes)
x <- 1:max(sizes) + 0.5
hist(sizes, breaks = x)



i <- 1

l <- list();

clu_mem <- read.table("sample1.2_cancer_005s_multilevel_clu.txt",
                      stringsAsFactors = FALSE)[ , 1]
sizes <- as.vector(table(clu_mem))

l[["cancer"]] <- sizes


clu_mem <- read.table("sample1.2_epicanexp_005s_multilevel_clu.txt",
                      stringsAsFactors = FALSE)[ , 1]
sizes <- as.vector(table(clu_mem))

l[["epican"]] <- sizes

boxplot(l, las = 2)











l <- list()

l <-

        boxplot(data, las = 2, col = c(“red”,“sienna”,“palevioletred1”,“royalblue2”,“red”,“sienna”,“palevioletred1”,
                                        “royalblue2”,“red”,“sienna”,“palevioletred1”,“royalblue2”),
                at = c(1,2,3,4, 6,7,8,9, 11,12,13,14), par(mar = c(12, 5, 4, 2) + 0.1),
                names = c(“Station 1”,“Station 2”,“Station 3”,“Station 4”,“Station 1”,“Station 2”,“Station 3”,“Station 4”,“Station 1”,“Station 2”,“Station 3”,“Station 4”))


