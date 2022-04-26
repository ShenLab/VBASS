DNVs <- read.csv('data/SPARK.WES1.csv', row.names = 1)
trios <- read.csv("data/SPARK.WES1.trios.csv", row.names = 1)
seed <- 0
set.seed(as.numeric(seed))
discovery.trios <- as.character(sample(trios$IID, floor(dim(trios)[1]*2/3)))
validation.trios <- as.character(trios$IID[!trios$IID%in%discovery.trios])

discovery.DNVs <- DNVs[DNVs$IID%in%discovery.trios, ]
validation.DNVs <- DNVs[DNVs$IID%in%validation.trios, ]
cases <- validation.DNVs[!is.na(validation.DNVs$vclass),]

write.csv(discovery.DNVs, "data/DNVs.discovery.csv")
write.csv(validation.DNVs, "data/DNVs.validation.csv")

samplenumber <- length(unique(validation.trios))
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
saveRDS(result, file = paste0("data/SPARK_extTADA.WES1.validation.RDS"))
