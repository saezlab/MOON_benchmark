library(readr)
library(cosmosR)

zscore_final_clean <- as.data.frame(
  read_csv("data/cytosig/zscore_raw.csv"))

zscore_meta_clean <- as.data.frame(
  read_csv("data/cytosig/zscore_meta_clean.csv"))

zscore_final_clean <- zscore_final_clean[,c(1,which(names(zscore_final_clean) %in% zscore_meta_clean$id))]

zscore_final_clean <- zscore_final_clean[-which(rowSums(is.na(zscore_final_clean[,-1])) > 1000),]

table(zscore_meta_clean$treatment)

zscore_meta_clean$treatment <- gsub("^41BBL$","TNFSF9",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^Activin A$","INHBA",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^CXCL4$","PF4",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^CD40L$","CD40LG",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^GMCSF$","CSF2",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^IFN1$","IFNA1",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^IL12$","IL12A",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^IL1$","IL1A",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^IL23$","IL23A",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^IL36$","IL36A",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^MCSF$","CSF1",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^NO$","Metab__HMDB0003378_c",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^OPGL$","TNFSF11",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^PGE2$","PTGES2",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^TNFA$","TNF",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^TRAIL$","TNFSF10",zscore_meta_clean$treatment)
zscore_meta_clean$treatment <- gsub("^TWEAK$","TNFSF12",zscore_meta_clean$treatment)


treatments <- unique(zscore_meta_clean$treatment)

data("meta_network")
nodes <- unique(c(meta_network$source, meta_network$target))

treatments[treatments %in% nodes]
treatments[!(treatments %in% nodes)]

write_csv(zscore_final_clean, file = "data/cytosig/zscore_final_clean_filtered.csv")
write_csv(zscore_meta_clean, file = "data/cytosig/zscore_meta_clean_renamed.csv")
