extract <- function(filename) {
        
        subarea <- read.csv(filename, header=FALSE)
        subarea <- as.vector(t(subarea))
        subarea <- subarea[subarea != '']
        subarea <- subarea[order(subarea)]
        subarea <- paste('X', subarea, sep='')
        subarea
}

sample1.2_epithelial <- extract('1.2_epithelial.csv')
sample1.2_epithelial_expanded <- extract('1.2_epithelial_expanded.csv')
sample1.2_cancer <- extract('1.2_cancer.csv')
sample2.4_epithelial <- extract('2.4_epithelial.csv')
sample2.4_cancer <- extract('2.4_cancer.csv')
sample3.3_cancercore <- extract('3.3_cancercore.csv')
sample3.3_cancerperiphery <- extract('3.3_cancerperiphery.csv')

save(sample1.2_epithelial, sample1.2_epithelial_expanded, sample1.2_cancer, sample2.4_epithelial,
     sample2.4_cancer, sample3.3_cancercore, sample3.3_cancerperiphery, file='coordinate.RData')