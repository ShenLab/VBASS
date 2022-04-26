DNVs <- read.csv('data/SPARK.WES1.csv', row.names = 1)
trios <- read.csv("data/SPARK.WES1.trios.csv", row.names = 1)

samplenumber <- length(unique(trios$IID))

reference <- read.delim('data/mutrate.3mer.txt', na.strings = c("NA", "."))
blacklist <- read.delim("data/GENCODEV19_blacklist.txt", header = F)
reference = reference[!reference$GeneID %in% blacklist$V1,]
geneset <- as.character(reference$GeneID)

# model
modelfile = dir('extTADA/model/', '.R$')
for (ii in modelfile) {
  source(paste0('extTADA/model/', ii))
}
result = gene_set_exttada(geneset, cases, samplenumber, reference)
saveRDS(result, file = "data/SPARK_extTADA.WES1.RDS")