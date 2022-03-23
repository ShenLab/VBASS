modelfile = dir('extTADA/', '.R$')
modelfile = modelfile[modelfile != "gene_set_exttada.R"]
modelfile = modelfile[modelfile != "gene_set_exttada_single.R"]
modelfile = modelfile[modelfile != "simulation_xTADA.2.R"]
for (ii in modelfile) {
  sourcedir = paste0('extTADA/', ii)
  source(sourcedir)
}

args = commandArgs(trailingOnly=TRUE)
C <- as.numeric(args[1])
C <- 0.28242408
samplesize <- as.numeric(args[2])
seed <- as.numeric(args[3])

reference <- read.table("~/Data/mutRate/EA.TEF.mutrate.3mer.txt", sep = "\t", header = TRUE, na.strings = ".")

trios <- read.delim("~/Data/PCGC.anno.hg19/PCGC.trios.txt", na.strings = c("NA", "."))
trios <- trios[grep("JinSC2017", trios$Publication),]

result = simulate_xTADA(reference, samplesize, effect_size_mean = c(19.94636306, 11.78832503),
                        effect_size_var = c(0.83709223, 0.89577342),
                        pi = 0.03728875,
                        A = 104.14597200, B = 0.73660024, C = C, seed=seed, repeat.size=1)

dnv_table <- result[[1]]
mutRate <- result[[2]]
covariates_array <- result[[3]]
risk_geneID <- result[[4]]
gene_set <- as.character(rownames(dnv_table))

simulation_result = sim_gene_set_exttada(dnv_table, mutRate, gene_set, samplesize)
saveRDS(simulation_result, file = paste0("RDS.files/extTADA_simulation_c.", C ,".size.",
                                         samplesize, ".seed.", seed, ".RDS"))

extTADAsim_post <- simulation_result$dataFDR

# Calculate Mutation Sum and Mutation Rank
mutRate$mutSum = mutRate$mut.cls1 + mutRate$mut.cls2
mutRate$mutRank = dplyr::percent_rank(mutRate$mutSum)
mutRate$Gene = rownames(mutRate)
rownames(mutRate) = NULL

extTADAsim_post$mutSum = extTADAsim_post$mut_LGD + extTADAsim_post$mut_Dmis

extTADAsim_post$true <- extTADAsim_post$Gene %in% risk_geneID
saveRDS(extTADAsim_post, file = paste0("RDS.files/extTADA_simulation_c.", C ,".size.",
                                       samplesize, ".seed.", seed, ".posterior.RDS"))


