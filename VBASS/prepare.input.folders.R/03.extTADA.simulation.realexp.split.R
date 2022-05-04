args = commandArgs(trailingOnly=TRUE)
input.batch <- args[1]

train.input <- read.csv(paste0('VBASS.simulation.with.RDS/input.x.',
                               input.batch, '.csv'),
                        row.names = 1)
train.input.var <- read.csv(paste0('VBASS.simulation.with.RDS/input.x_var.',
                                   input.batch, '.csv'),
                            row.names = 1)

input.dnv.table <- cbind(train.input[,c('dn.cls1', 'dn.cls2')],
                         train.input.var[,c('mut.cls1', 'mut.cls2')]/16616)

# reference <- input.dnv.table[,c('mut.cls1', 'mut.cls2')]

reference <- read.delim('data/mutrate.3mer.txt', na.strings = c("NA", "."))
blacklist <- read.delim("data/GENCODEV19_blacklist.txt", header = F)
reference = reference[!reference$GeneID %in% blacklist$V1,]

modelfile = dir('extTADA/model/', '.R$')
for (ii in modelfile) {
  source(paste0('extTADA/model/', ii))
}
samplenumber <- 16616
geneset <- as.character(row.names(input.dnv.table))
result = gene_set_dnvtable_exttada(geneset, input.dnv.table, samplenumber, reference)
saveRDS(result, file = paste0("VBASS.simulation.with.RDS/extTADA.simulation.realexp.",
                              input.batch, ".RDS"))
