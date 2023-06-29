# your input here
# cases should be a dataframe object with at least following columns:
# GeneID  vclass
# where:
#   GeneID should be ensemble ID
#   vclass should contain three categories: 
#   "syn", "LGD", "Dmis". Other annotation won't be recognized.
# samplenumber shoud be an integer referring to the trio number
cases = read.csv('data/Jin_2017.csv') # your input
samplenumber = 2645 # your input
# load reference files
reference = read.table("data/mutrate.3mer.txt", sep = "\t", header = TRUE, na.strings = ".")
blacklist <- read.table("data/GENCODEV19_blacklist.txt", sep = "\t", header = F)
reference = reference[!reference$GeneID %in% blacklist$V1,]

geneset = unique(reference$GeneID)
# load mouse heart dev files
mouse_dev = read.table('data/mouse_br_rnaseq3.rda.txt', header = TRUE, fill = TRUE)
mouse_dev = mouse_dev[!is.na(mouse_dev$e14.5_mean)&!is.na(mouse_dev$human.Ensembl.Gene.ID),]

rankPercentileData <- data.frame(gene=mouse_dev$human.Ensembl.Gene.ID, value=mouse_dev$e14.5_mean)
reference$mouse_dev <- rankPercentileData$value[match(reference$GeneID, rankPercentileData$gene)]
reference$mouse_dev_rank <- dplyr::percent_rank(reference$mouse_dev)
# if is na, set rank as 0.5
reference$mouse_dev_rank[is.na(reference$mouse_dev_rank)] = 0.5
# model
modelfile = dir('model/', '.R$')
for (ii in modelfile) {
  source(paste0('model/', ii))
}
result = gene_set_exttada(geneset, cases, samplenumber, reference, rankPercentileData)
saveRDS(result, file = "result.RDS")

