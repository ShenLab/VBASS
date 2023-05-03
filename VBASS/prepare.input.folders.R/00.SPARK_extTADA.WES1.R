DNVs <- read.csv('data/SPARK.WES1.csv', row.names = 1)
trios <- read.csv("data/SPARK.WES1.trios.csv", row.names = 1)

samplenumber <- length(unique(trios$IID))

reference <- read.delim('data/mutrate.3mer.txt', na.strings = c("NA", "."))
blacklist <- read.delim("data/GENCODEV19_blacklist.txt", header = F)
reference = reference[!reference$GeneID %in% blacklist$V1,]
geneset <- as.character(reference$GeneID)
args <- commandArgs(trailingOnly = T)
seed <- as.numeric(args[1])
set.seed(seed)
# model
modelfile = dir('extTADA/model/', '.R$')
for (ii in modelfile) {
  source(paste0('extTADA/model/', ii))
}
result = gene_set_exttada(geneset, DNVs, samplenumber, reference)
saveRDS(result, file = paste0("data/SPARK_extTADA.WES1.seed.", seed, ".RDS"))

for (i in 0:9) {
  tmp <- (readRDS(paste0('data/SPARK_extTADA.WES1.seed.', i, '.RDS')))
  # that will be the optimized KL parameter
  print(log(1/tmp$pars0[1,1]-1)/2)
}
