library("rstan")
options(mc.cores = parallel::detectCores())
########################################
#######De novo only
DNscalarVBASS <- readChar(paste0("scalar_VBASS/", "/model.stan"),
                          file.info(paste0("scalar_VBASS/", "/model.stan"))$size)

