library(png)
library(ggplot2)

# load ST data
raw = read.csv("/scratch/friedrich_stefanie/data/2_DIRINTH/18_STdata/all_deconvoluted/sample1.2_deconvolution.csv", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
head(raw)

# pre-processing
tmp = as.data.frame(noquote(unlist(lapply(strsplit(rownames(raw), split = " "),'[',1))))
data = cbind.data.frame (tmp, raw)
colnames(data) = c("gene", colnames(raw) )
head(data)

# search 
search = "YWHAE"
gene = t(data[grep(search, data[,1]), ][,-1]); #gene
data_search = c()

# pre-process earch result
list = strsplit(as.character(rownames(gene)),'x') ; xaxis = sapply(list, function (x) x[1]); xaxis = as.integer(substr(xaxis, 2, nchar(xaxis)) )
list = strsplit(as.character(rownames(gene)),'x') ; yaxis = sapply(list, function (x) x[2]); yaxis = as.integer(yaxis)
#data_search= cbind.data.frame(gene, xaxis, yaxis)
gene[is.na(gene)] <- 0;
head(gene)

#data = merge(x = (norm), y = (gene), by = 0 , all.x = TRUE)

#par(mfrow = c(3,4))
theme = theme(panel.grid.minor.x = element_line(colour = "grey"), panel.grid.minor.y = element_line(colour = "grey"),
              panel.grid.major.x = element_line(colour = "grey"), panel.grid.major.y = element_line(colour = "grey"),
              axis.title.x=element_blank() , axis.title.y=element_blank(), 
              axis.text.x=element_text(colour="grey20",size=8) , axis.text.y=element_text(colour="grey20",size=8), 
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              legend.position="none" ,aspect.ratio=1) 

p  =  ggplot(data = gene) + theme + ggtitle(search) +
      scale_x_continuous(limits = c(1, 33), breaks = seq(0, 33, 10) , position = "top") + 
      scale_y_reverse() +  
      geom_point(aes(x=gene$xaxis, y=gene$yaxis), size = 5, color = "grey89") ; p
p1 =  p + geom_point(aes(x=gene$xaxis, y=gene$yaxis, size = 5, color = gene[,1])) +
      scale_color_gradient(limit= c(0,max(gene[,1])), low="grey89", high="darkred") ; p1
      #scale_color_gradient(limit= c(2,max(gene[,1])), low="yellow", high="darkred")  
      #scale_size_continuous(range=c(0, 35)) #+ scale_size_area() 
      #geom_point(aes(x = gene$norm1, y = gene$norm2, size = ifelse(gene[,4] == "yes", 16, 0)), color="black") 

file = paste0("filename_gene.", search, ".pdf")
pdf(file)
dev.off()


grid.arrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, ncol=3)    


