cases <- read.csv('data/Jin_2017.csv')
samplenumber <- 2645
reference <- read.table("data/mutrate.3mer.txt", sep = "\t", header = TRUE, na.strings = ".")
blacklist <- read.table("data/GENCODEV19_blacklist.txt", sep = "\t", header = F)
reference = reference[!reference$GeneID %in% blacklist$V1,]
# load gnomAD score
expression <- read.table("data/mouse_br_rnaseq3.rda.txt", sep = "\t", header = TRUE, na.strings = ".")
# model
modelfile = dir('scalar_VBASS/', '.R$')
for (ii in modelfile) {
  sourcedir = paste0('scalar_VBASS/', ii)
  print(sourcedir)
  source(sourcedir)
}

geneset <- as.character(reference$GeneID)

covariatesData <- data.frame(gene=as.character(reference$GeneID),
                             expression=as.numeric(expression$`e14.5_rank`[match(reference$HGNC, expression$human.External.Gene.Name)]))
result = gene_set_scalarVBASS(geneset, cases, samplenumber, reference, covariatesData)
saveRDS(result, file = "result/CHD_VBASS_expression.RDS")

modelfile = dir('extTADA/', '.R$')
modelfile = modelfile[modelfile != "gene_set_exttada_single.R"]
modelfile = modelfile[modelfile != "sim_gene_set_extTADA.R"]

for (ii in modelfile) {
  sourcedir = paste0('extTADA/', ii)
  print(sourcedir)
  source(paste0('extTADA/', ii))
}
result = gene_set_exttada(geneset, cases, samplenumber, reference)
saveRDS(result, file = "result/CHD_extTADA.RDS")