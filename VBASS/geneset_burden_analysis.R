get.mytada.data.mutrate <- function(geneset, cases, reference, Dmis_def) {
  # cls0: syn, cls1: LGD, cls2: Dmis(should give Dmis_def), cls3:mis
  # automatically check geneset is ID or Symbol
  if (sum(geneset %in% reference$GeneID) >= sum(geneset %in% reference$HGNC)) {
    index.to.use = "GeneID"
  } else {
    index.to.use = "HGNC"
  }
  mytada.data = data.frame(gene.id=geneset,
                           mut.cls0=numeric(length(geneset)),
                           mut.cls1=numeric(length(geneset)),
                           mut.cls2=numeric(length(geneset)),
                           mut.cls3=numeric(length(geneset)),
                           dn.cls0=numeric(length(geneset)),
                           dn.cls1=numeric(length(geneset)),
                           dn.cls2=numeric(length(geneset)),
                           dn.cls3=numeric(length(geneset)),
                           entrez.id=reference$EntrezID[match(geneset, reference[,index.to.use])])
  all.indexes <- seq(1, dim(reference)[1])
  for (i in 1:length(geneset)) {
    # check whether we found the gene
    index = all.indexes[reference[,index.to.use]==as.character(geneset[i])&!is.na(reference[,index.to.use])]
    # get mut rate
    if (length(index)!=0) {
      if (length(index)>1) {
        message(paste0("Warning: ", geneset[i], " has multiple rows in reference"))
      }
      mytada.data$mut.cls0[i] = sum(reference$Mu_Silent[index])
      mytada.data$mut.cls1[i] = sum(reference$Mu_LoF[index])
      if (!paste0("Mu_Dmis_", Dmis_def) %in% colnames(reference)) {
        message(paste0("Dmis def: ", Dmis_def, ", mutrate not found in reference!"))
        stop()
      } else {
        mytada.data$mut.cls2[i] = sum(reference[index, paste0("Mu_Dmis_", Dmis_def)])
      }
      mytada.data$mut.cls3[i] = sum(reference$Mu_Missense[index])
    } else {
      message(paste0("Gene: ", geneset[i], ", not found in reference!"))
      mytada.data$mut.cls0[i] = .Machine$double.eps
      mytada.data$mut.cls1[i] = .Machine$double.eps
      mytada.data$mut.cls2[i] = .Machine$double.eps
      mytada.data$mut.cls3[i] = .Machine$double.eps
    }
    # get case, note here Dmis just mean mis
    if (index.to.use == "HGNC") {
      gene_cases = cases[cases$Symbol==as.character(geneset[i]),]
      gene_cases = gene_cases[!is.na(gene_cases$Symbol),]
    } else {
      gene_cases = cases[grep(as.character(geneset[i]), cases$GeneID),]
      gene_cases = gene_cases[!is.na(gene_cases$GeneID),]
    }
    mytada.data$dn.cls0[i] = sum(!is.na(gene_cases$vclass) & gene_cases$vclass=="syn")
    mytada.data$dn.cls1[i] = sum(!is.na(gene_cases$vclass) & gene_cases$vclass=="LGD")
    mytada.data$dn.cls2[i] = sum(!is.na(gene_cases$vclass) & gene_cases$vclass=="Dmis")
    mytada.data$dn.cls3[i] = sum(!is.na(gene_cases$vclass) & gene_cases$vclass%in%c("Dmis", "mis"))
  }
  mytada.data
}

geneset_burden_analysis_mutrate <- function(geneset, cases, samplenumber_case, reference, Dmis_def) {
  # form a proper input for TADA
  mytada.data = get.mytada.data.mutrate(geneset, cases, reference, Dmis_def)
  message("sanity check")
  syn_oe = sum(mytada.data$dn.cls0) / (sum(mytada.data$mut.cls0) * samplenumber_case * 2)
  message(paste0("syn_oe = ", syn_oe))
  # mytada.data[,2:4]=mytada.data[,2:4]*syn_oe
  
  syn_oe = sum(mytada.data$dn.cls0) / (sum(mytada.data$mut.cls0) * samplenumber_case * 2)
  LGD_oe = sum(mytada.data$dn.cls1) / (sum(mytada.data$mut.cls1) * samplenumber_case * 2)
  Dmis_oe = sum(mytada.data$dn.cls2) / (sum(mytada.data$mut.cls2) * samplenumber_case * 2)
  mis_oe = sum(mytada.data$dn.cls3) / (sum(mytada.data$mut.cls3) * samplenumber_case * 2)
  
  pois_syn = 1-ppois(sum(mytada.data$dn.cls0)-1, (sum(mytada.data$mut.cls0) * samplenumber_case * 2))
  pois_LGD = 1-ppois(sum(mytada.data$dn.cls1)-1, (sum(mytada.data$mut.cls1) * samplenumber_case * 2))
  pois_Dmis = 1-ppois(sum(mytada.data$dn.cls2)-1, (sum(mytada.data$mut.cls2) * samplenumber_case * 2))
  pois_mis = 1-ppois(sum(mytada.data$dn.cls3)-1, (sum(mytada.data$mut.cls3) * samplenumber_case * 2))
  
  result = data.frame(obs_syn = sum(mytada.data$dn.cls0),
                      exp_syn = (sum(mytada.data$mut.cls0) * samplenumber_case * 2),
                      burden_syn = syn_oe,
                      p_syn = pois_syn,
                      obs_LGD = sum(mytada.data$dn.cls1),
                      exp_LGD = (sum(mytada.data$mut.cls1) * samplenumber_case * 2),
                      burden_LGD = LGD_oe,
                      p_LGD = pois_LGD,
                      obs_Dmis = sum(mytada.data$dn.cls2),
                      exp_Dmis = (sum(mytada.data$mut.cls2) * samplenumber_case * 2),
                      burden_Dmis = Dmis_oe,
                      p_Dmis = pois_Dmis,
                      obs_mis = sum(mytada.data$dn.cls3),
                      exp_mis = (sum(mytada.data$mut.cls3) * samplenumber_case * 2),
                      burden_mis = mis_oe,
                      p_mis = pois_mis
                      )
}

geneset_burden_analysis_dnvtable_mutrate <- function(geneset, dnv_table, samplenumber_case) {
  mytada.data <- dnv_table[match(geneset, rownames(dnv_table)), ]
  message("sanity check")
  if (is.null(mytada.data$dn.cls0)) {
    mytada.data$dn.cls0 <- 0
  }
  if (is.null(mytada.data$mut.cls0)) {
    mytada.data$mut.cls0 <- .Machine$double.eps
  }
  syn_oe = sum(mytada.data$dn.cls0) / (sum(mytada.data$mut.cls0) * samplenumber_case * 2)
  message(paste0("syn_oe = ", syn_oe))
  # mytada.data[,2:4]=mytada.data[,2:4]*syn_oe
  
  syn_oe = sum(mytada.data$dn.cls0) / (sum(mytada.data$mut.cls0) * samplenumber_case * 2)
  LGD_oe = sum(mytada.data$dn.cls1) / (sum(mytada.data$mut.cls1) * samplenumber_case * 2)
  Dmis_oe = sum(mytada.data$dn.cls2) / (sum(mytada.data$mut.cls2) * samplenumber_case * 2)
  if (!is.null(mytada.data$dn.cls3)) {
    mis_oe = sum(mytada.data$dn.cls3) / (sum(mytada.data$mut.cls3) * samplenumber_case * 2)
    pois_mis = 1-ppois(sum(mytada.data$dn.cls3)-1, (sum(mytada.data$mut.cls3) * samplenumber_case * 2))
  } else {
    mytada.data$dn.cls3 <- 0
    mytada.data$mut.cls3 <- .Machine$double.eps
    mis_oe <- 0
    pois_mis <- 1
  }
  pois_syn = 1-ppois(sum(mytada.data$dn.cls0)-1, (sum(mytada.data$mut.cls0) * samplenumber_case * 2))
  pois_LGD = 1-ppois(sum(mytada.data$dn.cls1)-1, (sum(mytada.data$mut.cls1) * samplenumber_case * 2))
  pois_Dmis = 1-ppois(sum(mytada.data$dn.cls2)-1, (sum(mytada.data$mut.cls2) * samplenumber_case * 2))
  
  result = data.frame(obs_syn = sum(mytada.data$dn.cls0),
                      exp_syn = (sum(mytada.data$mut.cls0) * samplenumber_case * 2),
                      burden_syn = syn_oe,
                      p_syn = pois_syn,
                      obs_LGD = sum(mytada.data$dn.cls1),
                      exp_LGD = (sum(mytada.data$mut.cls1) * samplenumber_case * 2),
                      burden_LGD = LGD_oe,
                      p_LGD = pois_LGD,
                      obs_Dmis = sum(mytada.data$dn.cls2),
                      exp_Dmis = (sum(mytada.data$mut.cls2) * samplenumber_case * 2),
                      burden_Dmis = Dmis_oe,
                      p_Dmis = pois_Dmis,
                      obs_mis = sum(mytada.data$dn.cls3),
                      exp_mis = (sum(mytada.data$mut.cls3) * samplenumber_case * 2),
                      burden_mis = mis_oe,
                      p_mis = pois_mis
  )
}


geneset_burden_analysis_casecontrol <- function(geneset, cases, controls, samplenumber_case, samplenumber_control) {
  # form a proper input for TADA
  mytada.data = data.frame(gene.id=geneset,
                           mut.cls0=numeric(length(geneset)),
                           mut.cls1=numeric(length(geneset)),
                           mut.cls2=numeric(length(geneset)),
                           mut.cls3=numeric(length(geneset)),
                           dn.cls0=numeric(length(geneset)),
                           dn.cls1=numeric(length(geneset)),
                           dn.cls2=numeric(length(geneset)),
                           dn.cls3=numeric(length(geneset)))
  for (i in 1:length(geneset)) {
    index = match(as.character(geneset[i]), reference$HGNC)
    # check whether we found the gene
    gene_controls = controls[controls$Symbol==as.character(geneset[i]),]
    mytada.data$mut.cls0[i] = sum(gene_controls$vclass=="syn")
    mytada.data$mut.cls1[i] = sum(gene_controls$vclass=="LGD")
    mytada.data$mut.cls2[i] = sum(gene_controls$vclass=="Dmis")
    mytada.data$mut.cls3[i] = sum(gene_controls$vclass=="Dmis") + sum(gene_controls$vclass=="mis")
    
    gene_cases = cases[cases$Symbol==as.character(geneset[i]),]
    mytada.data$dn.cls0[i] = sum(gene_cases$vclass=="syn")
    mytada.data$dn.cls1[i] = sum(gene_cases$vclass=="LGD")
    mytada.data$dn.cls2[i] = sum(gene_cases$vclass=="Dmis")
    mytada.data$dn.cls3[i] = sum(gene_cases$vclass=="Dmis") + sum(gene_cases$vclass=="mis")
  }
  message("sanity check before correction")
  syn_oe = sum(mytada.data$dn.cls0) / (sum(mytada.data$mut.cls0) * samplenumber_case / samplenumber_control)
  message(paste0("syn_oe = ", syn_oe))
  message("adjusting mutation rate according to syn_oe")
  # mytada.data[,2:4]=mytada.data[,2:4]*syn_oe
  
  syn_oe = sum(mytada.data$dn.cls0) / (sum(mytada.data$mut.cls0) * samplenumber_case / samplenumber_control)
  LGD_oe = sum(mytada.data$dn.cls1) / (sum(mytada.data$mut.cls1) * samplenumber_case / samplenumber_control)
  Dmis_oe = sum(mytada.data$dn.cls2) / (sum(mytada.data$mut.cls2) * samplenumber_case / samplenumber_control)
  mis_oe = sum(mytada.data$dn.cls3) / (sum(mytada.data$mut.cls3) * samplenumber_case / samplenumber_control)
  
  pois_syn = ppois(sum(mytada.data$dn.cls0), (sum(mytada.data$mut.cls0) * samplenumber_case / samplenumber_control), lower.tail = F)
  pois_LGD = ppois(sum(mytada.data$dn.cls1), (sum(mytada.data$mut.cls1) * samplenumber_case / samplenumber_control), lower.tail = F)
  pois_Dmis = ppois(sum(mytada.data$dn.cls2), (sum(mytada.data$mut.cls2) * samplenumber_case / samplenumber_control), lower.tail = F)
  pois_mis = ppois(sum(mytada.data$dn.cls3), (sum(mytada.data$mut.cls3) * samplenumber_case / samplenumber_control), lower.tail = F)
  
  result = data.frame(obs_syn = sum(mytada.data$dn.cls0),
                      exp_syn = (sum(mytada.data$mut.cls0) * samplenumber_case / samplenumber_control),
                      burden_syn = syn_oe,
                      p_syn = pois_syn,
                      obs_LGD = sum(mytada.data$dn.cls1),
                      exp_LGD = (sum(mytada.data$mut.cls1) * samplenumber_case / samplenumber_control),
                      burden_LGD = LGD_oe,
                      p_LGD = pois_LGD,
                      obs_Dmis = sum(mytada.data$dn.cls2),
                      exp_Dmis = (sum(mytada.data$mut.cls2) * samplenumber_case / samplenumber_control),
                      burden_Dmis = Dmis_oe,
                      p_Dmis = pois_Dmis,
                      obs_mis = sum(mytada.data$dn.cls3),
                      exp_mis = (sum(mytada.data$mut.cls3) * samplenumber_case / samplenumber_control),
                      burden_mis = mis_oe,
                      p_mis = pois_mis
  )
}