library(readr)
library(cosmosR)
library(decoupleR)
library(reshape2)
library(purrr)
library(dplyr)
library(ggplot2)
library(yardstick)
library(ggrepel)

source("scripts/support_functions.R")

zscores <- as.data.frame(
  read_csv(file = "data/cytosig/zscore_final_clean_filtered.csv"))
meta_data <- as.data.frame(
  read_csv(file = "data/cytosig/zscore_meta_clean_renamed.csv"))

# collectri <- decoupleR::get_collectri()
# save(collectri, file = "support/collectri.RData")
load("support/collectri.RData")

# data("meta_network")
# save(meta_network, file = "support/neta_network.RData")
load("support/neta_network.RData")

##let's try with an updatedo mnipath metaPKN
# load("../meta_PKN_BIGG/results/meta_PKN.RData")
# meta_network <- meta_PKN

nodes <- unique(c(meta_network$source, meta_network$target))

collectri <- collectri[collectri$source %in% nodes,]

row.names(zscores) <- zscores$gene
zscores <- zscores[,-1]

collectri <- collectri[collectri$target %in% row.names(zscores),]

#run this chunk to re-estimate the TF activties
# TF_activities <- apply(zscores,2,function(x, collectri){
#   TF_act <- decoupleR::run_ulm(na.omit(x), network = collectri)
#   return(TF_act[,c(2,4)])
# },collectri = collectri)
# 
# for(i in 1:length(TF_activities))
# {
#   names(TF_activities[[i]])[2] <- names(TF_activities)[i]
# }
# 
# save(TF_activities, file = "results/citosig_TF_activties.RData")
load("results/citosig_TF_activties.RData")

TF_activities_df <- TF_activities %>% reduce(full_join, by = "source")
TF_activities_df <- as.data.frame(TF_activities_df)

n_steps <- 10

#run this chunl to re-estimate the moon scores
# moon_res_list <- list()
# for(exp_counter in 1:length(zscores[1,]))
# {
#   print(exp_counter)
#   RNA_input <- zscores[,exp_counter]
#   names(RNA_input) <- row.names(zscores)
#   RNA_input <- RNA_input[-which(is.na(RNA_input))]
#   
#   meta_network_filtered <- meta_network_cleanup(meta_network)
#   meta_network_filtered <- cosmosR:::filter_pkn_expressed_genes(names(RNA_input), meta_pkn = meta_network_filtered)
#   
#   TF_input <- TF_activities_df[,exp_counter+1]
#   names(TF_input) <- TF_activities_df[,1]
#   
#   if(sum(is.na(TF_input)) > 0)
#   {
#     TF_input <- TF_input[-which(is.na(TF_input))]
#   }
#   TF_input_filtered <- cosmosR:::filter_input_nodes_not_in_pkn(TF_input, meta_network_filtered)
#   meta_network_filtered <- cosmosR:::keep_observable_neighbours(meta_network_filtered, n_steps = 10, observed_nodes = names(TF_input_filtered))
#   
#   meta_network_compressed_list <- compress_same_children(meta_network_filtered, sig_input = c(0),metab_input = TF_input_filtered)
#   
#   meta_network_compressed <- meta_network_compressed_list$compressed_network
#   
#   node_signatures <- meta_network_compressed_list$node_signatures
#   
#   duplicated_parents <- meta_network_compressed_list$duplicated_signatures
#   
#   meta_network_compressed <- meta_network_cleanup(meta_network_compressed)
#   
#   # test <- decoupleRnival(downstream_input = decoupleRnival_input, meta_network = meta_network_compressed, n_layers = 3, statistic = "ulm")
#   
#   meta_network_compressed_to_run <- meta_network_compressed
#   
#   #We run moon in a loop until TF-target coherence convergences
#   before <- 1
#   after <- 0
#   i <- 1
#   while (before != after & i < 10) {
#     before <- length(meta_network_compressed_to_run[,1])
#     moon_res <- cosmosR::moon(downstream_input = TF_input_filtered, 
#                               meta_network = meta_network_compressed_to_run, 
#                               n_layers = n_steps, 
#                               statistic = "ulm") 
#     
#     meta_network_compressed_to_run <- filter_incohrent_TF_target(moon_res, collectri, meta_network_compressed_to_run, RNA_input)
#     
#     after <- length(meta_network_compressed_to_run[,1])
#     i <- i + 1
#   }
#   
#   if(i < 10)
#   {
#     print(paste("Converged after ",paste(i-1," iterations", sep = ""),sep = ""))
#   } else
#   {
#     print(paste("Interupted after ",paste(i," iterations. Convergence uncertain.", sep = ""),sep = ""))
#   }
#   
#   moon_res <- decompress_moon_result(moon_res, meta_network_compressed_list, meta_network_compressed_to_run)
#   names(moon_res)[2] <- names(TF_activities_df)[exp_counter+1]
#   
#   moon_res_list[[exp_counter]] <- moon_res
# }
# 
# save(moon_res_list, file = "results/citosig_moon_activties.RData")
load("results/citosig_moon_activties.RData")

moon_res_list_filtered <- lapply(moon_res_list, function(x){x <- x[which(x$level != 0),c(4,2)]
x <- unique(x)
return(x)})

moon_res_df <- moon_res_list_filtered %>% reduce(full_join, by = "source_original")
moon_res_df <- as.data.frame(moon_res_df)

moon_res_df_long <- melt(moon_res_df)
moon_res_df_long$target <- gsub("@.*","",moon_res_df_long$variable)

meta_data$targets <- gsub("[@&].*","",meta_data$id)
targets <- unique(meta_data[,c(12,7)])
maping_vector <- setNames(targets$treatment, nm = targets$targets)

moon_res_df_long$target <- sapply(moon_res_df_long$target, function(x, maping_vector){
  return(maping_vector[x])
},maping_vector = maping_vector)

moon_res_df_long <- moon_res_df_long[,c(2,1,3,4)]

names(moon_res_df_long) <- c("id", "source", "score", "target")
moon_res_df_long$id <- as.character(moon_res_df_long$id)

roc <- decoupleRBench::calc_curve(moon_res_df_long)

ggplot(roc, aes(x = 1-specificity,
                y = sensitivity)) +
  geom_line() +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  xlab("FPR (1-specificity)") +
  ylab("TPR (sensitivity)")

moon_res_df_long_yarstick <- moon_res_df_long
moon_res_df_long_yarstick$groundtruth <- as.numeric(moon_res_df_long_yarstick$source == moon_res_df_long_yarstick$target)
bad_experiments <- moon_res_df_long_yarstick[moon_res_df_long_yarstick$groundtruth & is.na(moon_res_df_long_yarstick$score),"id"]
bad_experiments <- unique(bad_experiments)
moon_res_df_long_yarstick <- moon_res_df_long_yarstick[!(moon_res_df_long_yarstick$id %in% bad_experiments),]

moon_res_df_long_yarstick$groundtruth <- factor(moon_res_df_long_yarstick$groundtruth, levels = c(1,0))
moon_res_df_long_yarstick <- moon_res_df_long_yarstick[complete.cases(moon_res_df_long_yarstick),]

ligands <- unique(moon_res_df_long_yarstick$target)

for(ligand in ligands)
{
  print(ligand)
  moon_res_df_long_yarstick_sub <- moon_res_df_long_yarstick[moon_res_df_long_yarstick$target == ligand,]
  if(sum(as.numeric(as.character(moon_res_df_long_yarstick_sub$groundtruth))) == 0)
  {
    moon_res_df_long_yarstick <- moon_res_df_long_yarstick[!(moon_res_df_long_yarstick$target == ligand),]
  }
}

ROC_per_ligand <- moon_res_df_long_yarstick %>% as_tibble() %>% group_by(target) %>% roc_auc(groundtruth, score)

experiment_count <- unique(moon_res_df_long_yarstick[,c(1,4)])
experiment_count <- merge(experiment_count,data.frame(table(experiment_count$target)), by.x = "target", by.y = "Var1")
experiment_count <- unique(experiment_count[,c(1,3)])

ROC_per_ligand <- merge(ROC_per_ligand, experiment_count, by = "target")

ggplot(ROC_per_ligand, aes(x = Freq, y = .estimate, label = target)) + 
  geom_label_repel() +
  geom_point() +
  geom_hline(yintercept = 0.5) +
  theme_minimal()

ROC_per_ligand$Freq_binned <- as.character(floor(ROC_per_ligand$Freq / 20) * 20)
ROC_per_ligand <- ROC_per_ligand[order(ROC_per_ligand$Freq, decreasing = T),]
ROC_per_ligand$target <- factor(ROC_per_ligand$target, level = ROC_per_ligand$target)

ggplot(ROC_per_ligand, aes(x = target, y = .estimate, fill = Freq_binned)) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.5) +
  theme_minimal() + 
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) 

ligands_of_interest <- as.character(ROC_per_ligand[ROC_per_ligand$Freq > 6 & ROC_per_ligand$.estimate < 0.55,"target"])

moon_res_df_long_yarstick_receptors <- moon_res_df_long_yarstick[moon_res_df_long_yarstick$target %in% ligands_of_interest,]
receptors_to_merge <- meta_network[meta_network$source %in%ligands_of_interest,]
names(receptors_to_merge) <- c("target","sign","receptor")
moon_res_df_long_yarstick_receptors <- merge(moon_res_df_long_yarstick_receptors, receptors_to_merge, by = "target")
moon_res_df_long_yarstick_receptors$groundtruth <- as.numeric(moon_res_df_long_yarstick_receptors$source == moon_res_df_long_yarstick_receptors$target | moon_res_df_long_yarstick_receptors$source == moon_res_df_long_yarstick_receptors$receptor)
moon_res_df_long_yarstick_receptors$groundtruth <- factor(moon_res_df_long_yarstick_receptors$groundtruth, levels = c(1,0))
moon_res_df_long_yarstick_receptors[moon_res_df_long_yarstick_receptors$source == moon_res_df_long_yarstick_receptors$receptor,"score"] <- moon_res_df_long_yarstick_receptors[moon_res_df_long_yarstick_receptors$source == moon_res_df_long_yarstick_receptors$receptor,"score"] * moon_res_df_long_yarstick_receptors[moon_res_df_long_yarstick_receptors$source == moon_res_df_long_yarstick_receptors$receptor,"sign"]

ROC_per_receptor <- moon_res_df_long_yarstick_receptors %>% as_tibble() %>% group_by(receptor) %>% roc_auc(groundtruth, score)

experiment_count <- unique(moon_res_df_long_yarstick_receptors[,c(2,7)])
experiment_count <- merge(experiment_count,data.frame(table(experiment_count$receptor)), by.x = "receptor", by.y = "Var1")
experiment_count <- unique(experiment_count[,c(1,3)])

ROC_per_receptor <- merge(ROC_per_receptor, experiment_count, by = "receptor")
ROC_per_receptor <- merge(ROC_per_receptor, receptors_to_merge, by = "receptor",)

# ggplot(ROC_per_receptor, aes(x = Freq, y = .estimate, label = target)) + 
#   geom_label_repel() +
#   geom_point() +
#   geom_hline(yintercept = 0.5) +
#   theme_minimal()

ROC_per_receptor$Freq_binned <- as.character(floor(ROC_per_receptor$Freq / 20) * 20)
ROC_per_receptor <- ROC_per_receptor[order(ROC_per_receptor$Freq, decreasing = T),]
ROC_per_receptor$receptor <- factor(ROC_per_receptor$receptor, level = ROC_per_receptor$receptor)

ggplot(ROC_per_receptor, aes(x = receptor, y = .estimate, fill = target, group = target)) + 
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.5) +
  theme_minimal() + 
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  ylim(c(0,1))


#I should make networks connecting e.g. FGF10 to downstream TF to show why the estimation is so poor.
#Will make such network for a few interesting cases.


moon_res_df_long_target_only <- moon_res_df_long[which(moon_res_df_long$source == moon_res_df_long$target),]

#we can check for each ligand how well they are scored (the higher the better)
ggplot(moon_res_df_long_target_only, aes(x = source, y = score, group = source)) + 
  geom_boxplot() + 
  geom_jitter() + 
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + geom_hline(yintercept = 1.7)
  
moon_res_df_long_target_only <- merge(moon_res_df_long_target_only, meta_data, by = "id")
write_csv(moon_res_df_long_target_only, file = "results/moon_res_df_long_target_only.csv")

#chatgpt cleaned up the time and added a mean score column
moon_res_df_long_target_only_GPTcleanedup <- as.data.frame(
  read_csv("results/moon_res_df_long_target_only_GPTcleanedup_withMeanScore.csv"))
moon_res_df_long_target_only_GPTcleanedup <- moon_res_df_long_target_only_GPTcleanedup[!is.na(moon_res_df_long_target_only_GPTcleanedup$score),]
moon_res_df_long_target_only_GPTcleanedup <- moon_res_df_long_target_only_GPTcleanedup[!is.na(moon_res_df_long_target_only_GPTcleanedup$id),]
moon_res_df_long_target_only_GPTcleanedup$loghours <- log(moon_res_df_long_target_only_GPTcleanedup$inferred_time_all_conditions)
moon_res_df_long_target_only_GPTcleanedup <- moon_res_df_long_target_only_GPTcleanedup[order(moon_res_df_long_target_only_GPTcleanedup$median_score_by_target, decreasing = T),]
moon_res_df_long_target_only_GPTcleanedup$target <- factor(moon_res_df_long_target_only_GPTcleanedup$target, levels = unique(moon_res_df_long_target_only_GPTcleanedup$target))
moon_res_df_long_target_only_GPTcleanedup$source <- factor(moon_res_df_long_target_only_GPTcleanedup$source, levels = unique(moon_res_df_long_target_only_GPTcleanedup$source))
moon_res_df_long_target_only_GPTcleanedup$quality <- ifelse(moon_res_df_long_target_only_GPTcleanedup$median_score_by_target > 1.7, "lightgreen",
                                                            ifelse(moon_res_df_long_target_only_GPTcleanedup$median_score_by_target > 0.3, "orange",
                                                            "red"))

quality_mapping <- unique(moon_res_df_long_target_only_GPTcleanedup[,c("source","quality")])
quality_mapping <- quality_mapping[complete.cases(quality_mapping),]
#we can also add some of the metadata as colors for example
ggplot(moon_res_df_long_target_only_GPTcleanedup, aes(x = source, y = score, group = source, color = loghours, fill = platformType)) + 
  geom_boxplot(fill = quality_mapping$quality, coef = 6) + 
  geom_jitter() + 
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + geom_hline(yintercept = 0)

moon_res_df_long_target_only_GPTcleanedup <- moon_res_df_long_target_only_GPTcleanedup %>%
  group_by(source) %>%
  mutate(mean_ligand_score = mean(score, na.rm = TRUE)) %>%
  ungroup()
moon_res_df_long_target_only_GPTcleanedup <- as.data.frame(moon_res_df_long_target_only_GPTcleanedup)

sum(unique(moon_res_df_long_target_only_GPTcleanedup$mean_ligand_score) > 0)
sum(unique(moon_res_df_long_target_only_GPTcleanedup$mean_ligand_score) < 0)

moon_res_df_long_target_only_GPTcleanedup <- moon_res_df_long_target_only_GPTcleanedup %>%
  group_by(source) %>%
  mutate(sd_ligand_score = sd(score, na.rm = TRUE)) %>%
  ungroup()
moon_res_df_long_target_only_GPTcleanedup <- as.data.frame(moon_res_df_long_target_only_GPTcleanedup)


#get the levels of the ligands
moon_level_list_filtered <- lapply(moon_res_list, function(x){
  sample <- names(x)[2]
  x <- x[which(x$level != 0),c(4,3)]
  x <- unique(x)
  names(x)[2] <- sample
return(x)})

moon_level_df <- moon_level_list_filtered %>% reduce(full_join, by = "source_original")
moon_level_df <- as.data.frame(moon_level_df)

moon_level_df_long <- melt(moon_level_df)
moon_level_df_long$target <- gsub("@.*","",moon_level_df_long$variable)

moon_level_df_long$target <- sapply(moon_level_df_long$target, function(x, maping_vector){
  return(maping_vector[x])
},maping_vector = maping_vector)

moon_level_df_long <- moon_level_df_long[,c(2,1,3,4)]

names(moon_level_df_long) <- c("id", "source", "level", "target")
moon_level_df_long$id <- as.character(moon_level_df_long$id)

moon_level_df_long_reduced <- moon_level_df_long[which(moon_level_df_long$source %in% moon_res_df_long_target_only_GPTcleanedup$source),]

moon_res_df_long_target_only_GPTcleanedup_level <- merge(moon_res_df_long_target_only_GPTcleanedup, moon_level_df_long_reduced, by = c("id","source","target"))
moon_res_df_long_target_only_GPTcleanedup_level$level <- as.character(moon_res_df_long_target_only_GPTcleanedup_level$level)

#for paper
ggplot(moon_res_df_long_target_only_GPTcleanedup_level, aes(x = source, y = score, group = source, color = level)) + 
  geom_boxplot(fill = quality_mapping$quality, coef = 6, color = "black", alpha = 1) + 
  geom_jitter() + 
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + geom_hline(yintercept = 0)


ggplot(moon_res_df_long_target_only_GPTcleanedup_level, aes(x = level, y = score)) + 
  geom_boxplot(coef = 7) + geom_jitter(alpha = 0.3) + theme_minimal()

for(level in unique(moon_res_df_long_target_only_GPTcleanedup_level$level))
{
  print(paste("level: ",level, sep = ""))
  print(mean(moon_res_df_long_target_only_GPTcleanedup_level[moon_res_df_long_target_only_GPTcleanedup_level$level == level, "score"]))
  print(median(moon_res_df_long_target_only_GPTcleanedup_level[moon_res_df_long_target_only_GPTcleanedup_level$level == level, "score"]))
  print(sd(moon_res_df_long_target_only_GPTcleanedup_level[moon_res_df_long_target_only_GPTcleanedup_level$level == level, "score"]))
  print(length(moon_res_df_long_target_only_GPTcleanedup_level[moon_res_df_long_target_only_GPTcleanedup_level$level == level, "score"]))
}

cor.test(moon_res_df_long_target_only_GPTcleanedup_level$score, as.numeric(moon_res_df_long_target_only_GPTcleanedup_level$level), method = "kendall")

moon_res_df_ranks <- moon_res_df
moon_res_df_ranks[,-1] <- as.data.frame(apply(moon_res_df_ranks[,-1], 2, function(x){rank(as.numeric(x), na.last = "keep") / length(na.omit(x))}))
##IFNA1
ligand_to_rank <- "WNT3A"
ligand_to_rank_samples <- moon_res_df_long_target_only_GPTcleanedup_level[moon_res_df_long_target_only_GPTcleanedup_level$source == ligand_to_rank,"id"]

mean(moon_res_df_long_target_only_GPTcleanedup_level[moon_res_df_long_target_only_GPTcleanedup_level$target == ligand_to_rank,"score"])
mean(as.numeric(moon_res_df[moon_res_df$source_original == ligand_to_rank,which(names(moon_res_df) %in% ligand_to_rank_samples)]), na.rm = T)
mean(as.numeric(moon_res_df[moon_res_df$source_original == ligand_to_rank,-c(1,which(names(moon_res_df) %in% ligand_to_rank_samples))]), na.rm = T)
plot(density(as.numeric(moon_res_df[moon_res_df$source_original == ligand_to_rank,-c(1,which(names(moon_res_df) %in% ligand_to_rank_samples))]), na.rm = T))
#get the relative rank
mean(as.numeric(moon_res_df_ranks[moon_res_df_ranks$source_original == ligand_to_rank,-c(1,which(names(moon_res_df_ranks) %in% ligand_to_rank_samples))]), na.rm = T)
mean(as.numeric(moon_res_df_ranks[moon_res_df_ranks$source_original == ligand_to_rank,which(names(moon_res_df_ranks) %in% ligand_to_rank_samples)]), na.rm = T)

##all ligands
ligand_scores_benchamrk <- list()
for(ligand_to_rank in unique(moon_res_df_long_target_only_GPTcleanedup_level$target))
{
  print(ligand_to_rank)
  ligand_to_rank_samples <- moon_res_df_long_target_only_GPTcleanedup_level[moon_res_df_long_target_only_GPTcleanedup_level$source == ligand_to_rank,"id"]
  mean_quantile_in_trueExp <- mean(as.numeric(moon_res_df_ranks[moon_res_df_ranks$source_original == ligand_to_rank,which(names(moon_res_df_ranks) %in% ligand_to_rank_samples)]), na.rm = T)
  mean_quantile_in_falseExp <- mean(as.numeric(moon_res_df_ranks[moon_res_df_ranks$source_original == ligand_to_rank,-c(1,which(names(moon_res_df_ranks) %in% ligand_to_rank_samples))]), na.rm = T)
  mean_score_in_trueExp <- mean(as.numeric(moon_res_df[moon_res_df$source_original == ligand_to_rank,which(names(moon_res_df) %in% ligand_to_rank_samples)]), na.rm = T)
  mean_score_in_falseExp <- mean(as.numeric(moon_res_df[moon_res_df$source_original == ligand_to_rank,-c(1,which(names(moon_res_df) %in% ligand_to_rank_samples))]), na.rm = T)
  ligand_scores_benchamrk[[ligand_to_rank]] <- c(ligand_to_rank, mean_quantile_in_trueExp, mean_quantile_in_falseExp, mean_score_in_trueExp, mean_score_in_falseExp)
}

ligand_scores_benchamrk_df <- as.data.frame(do.call(rbind, ligand_scores_benchamrk))
ligand_scores_benchamrk_df[,-1] <- as.data.frame(apply(ligand_scores_benchamrk_df[,-1], 2, as.numeric))
names(ligand_scores_benchamrk_df) <- c("ligand_to_rank", "mean_quantile_in_trueExp", "mean_quantile_in_falseExp", "mean_score_in_trueExp", "mean_score_in_falseExp")

ggplot(ligand_scores_benchamrk_df, aes(x = mean_quantile_in_falseExp, y = mean_quantile_in_trueExp)) + geom_point() + theme_minimal() + ylim(c(0,1)) + xlim(c(0,1)) + geom_abline(intercept = 0)
sum(ligand_scores_benchamrk_df$mean_quantile_in_trueExp > ligand_scores_benchamrk_df$mean_quantile_in_falseExp)

ligand_quantiles_long <- melt(ligand_scores_benchamrk_df[order(ligand_scores_benchamrk_df[,2],decreasing = F),c(1,2,3)])
ligand_quantiles_long$ligand_to_rank <- factor(ligand_quantiles_long$ligand_to_rank, levels = unique(ligand_quantiles_long$ligand_to_rank))

#good sup figure
ggplot(ligand_quantiles_long, aes(x = value, y = ligand_to_rank, color = variable)) + geom_point(size = 3) + theme_minimal()

ligand_score_long <- melt(ligand_scores_benchamrk_df[order(ligand_scores_benchamrk_df[,4],decreasing = F),c(1,4,5)])
ligand_score_long$ligand_to_rank <- factor(ligand_score_long$ligand_to_rank, levels = unique(ligand_score_long$ligand_to_rank))

ggplot(ligand_score_long, aes(x = value, y = ligand_to_rank, color = variable)) + geom_point(size = 3) + theme_minimal()

#just chekc if platform or other variable has higher scores
variable_regulon <- moon_res_df_long_target_only_GPTcleanedup[,c(1,9)]
variable_regulon$mor <- 1
names(variable_regulon) <- c("target","source","mor")
measurments_variable <- moon_res_df_long_target_only_GPTcleanedup[,c(3),drop = F]
row.names(measurments_variable) <- moon_res_df_long_target_only_GPTcleanedup$id
result_ulm <- run_ulm(as.matrix(measurments_variable), variable_regulon, minsize = 1)
# mean(moon_res_df_long_target_only_GPTcleanedup[moon_res_df_long_target_only_GPTcleanedup$platformType == "MicroArray",3])
# mean(moon_res_df_long_target_only_GPTcleanedup[moon_res_df_long_target_only_GPTcleanedup$platformType == "RNASeq",3])


cor.test(moon_res_df_long_target_only_GPTcleanedup$score, moon_res_df_long_target_only_GPTcleanedup$loghours, method = "kendall")
ggplot(moon_res_df_long_target_only_GPTcleanedup, aes(x = loghours, y = score)) + geom_jitter() +
  geom_smooth(method='lm', formula= y~x) + theme_minimal() + geom_hline(yintercept = 0)


##Let'S redo the ROC plot with a dataset ocnsistent with the previous one
moon_res_df_long_yarstick_filtered <- moon_res_df_long_yarstick[moon_res_df_long_yarstick$id %in% unique(moon_res_df_long_target_only_GPTcleanedup$id),]
ROC_per_ligand <- moon_res_df_long_yarstick_filtered %>% as_tibble() %>% group_by(target) %>% roc_auc(groundtruth, score)

experiment_count <- unique(moon_res_df_long_yarstick_filtered[,c(1,4)])
experiment_count <- merge(experiment_count,data.frame(table(experiment_count$target)), by.x = "target", by.y = "Var1")
experiment_count <- unique(experiment_count[,c(1,3)])

ROC_per_ligand <- merge(ROC_per_ligand, experiment_count, by = "target")

ggplot(ROC_per_ligand, aes(x = Freq, y = .estimate, label = target)) + 
  geom_label_repel() +
  geom_point() +
  geom_hline(yintercept = 0.5) +
  theme_minimal()

ROC_per_ligand$Freq_binned <- as.character(floor(ROC_per_ligand$Freq / 20) * 20)
ROC_per_ligand <- ROC_per_ligand[order(ROC_per_ligand$Freq, decreasing = T),]
ROC_per_ligand$target <- factor(ROC_per_ligand$target, level = ROC_per_ligand$target)

ggplot(ROC_per_ligand, aes(x = target, y = .estimate, fill = Freq_binned)) + 
  geom_bar(stat = "identity", position = "dodge") +
  # geom_hline(yintercept = 0.5) +
  theme_minimal() + 
  geom_text(aes(label=Freq), position=position_dodge(width=0.9), vjust=-0.25) +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + 
  geom_hline(yintercept = 0.55) +  geom_hline(yintercept = 0.45)

experiment_id <- "BMP4@Condition:embryonic stem cell@E-MEXP-1192.MicroArray.HG-U133_Plus_2"
expressed_genes <- zscores[,experiment_id,drop = F]
expressed_genes <- expressed_genes[complete.cases(expressed_genes),,drop = F]
expressed_genes <- setNames(expressed_genes[,1], nm = row.names(expressed_genes))

moon_res_list_named <- moon_res_list
names(moon_res_list_named) <- lapply(moon_res_list_named, function(x){return(names(x)[2])})

moon_res_experiment <- moon_res_list_named[[experiment_id]] 

ligand <- as.character(
  moon_res_df_long_target_only_GPTcleanedup[moon_res_df_long_target_only_GPTcleanedup$id == experiment_id,"source"])

n_steps <- moon_res_experiment[moon_res_experiment$source == ligand,"level"]

meta_network_filtered <- cosmosR:::filter_pkn_expressed_genes(expressed_genes_entrez = names(expressed_genes), meta_network)

#This function allow us to get the subnetowrk that yielded the ligand score
moon_scoring_network <- get_moon_scoring_network(ligand, meta_network_filtered, moon_res_experiment, T)

names(moon_scoring_network$SIF)[2] <- "sign"
names(moon_scoring_network$ATT)[2] <- "moon_score"

write_csv(moon_scoring_network$SIF, file = "results/networktest_SIF.csv")
write_csv(moon_scoring_network$ATT, file = "results/networktest_ATT.csv")


#let's plot only ligand with experiements under 24h
# ggplot(moon_res_df_long_target_only_GPTcleanedup[moon_res_df_long_target_only_GPTcleanedup$inferred_time_all_conditions <= 24,], aes(x = source, y = score, group = source)) + 
#   geom_boxplot() + 
#   geom_jitter() + 
#   theme_minimal() +
#   theme(axis.text.x=element_text(angle=45, hjust=1)) + geom_hline(yintercept = 1.7)
# 
# #let's try to order them by mean score
# moon_res_df_long_target_only_GPTcleanedup_less_than_24h <- moon_res_df_long_target_only_GPTcleanedup[moon_res_df_long_target_only_GPTcleanedup$inferred_time_all_conditions <= 24,]
# moon_res_df_long_target_only_GPTcleanedup_less_than_24h <- moon_res_df_long_target_only_GPTcleanedup_less_than_24h[!is.na(moon_res_df_long_target_only_GPTcleanedup_less_than_24h$loghours),]
# moon_res_df_long_target_only_GPTcleanedup_less_than_24h$source <- as.character(moon_res_df_long_target_only_GPTcleanedup_less_than_24h$source)
# moon_res_df_long_target_only_GPTcleanedup_less_than_24h <- moon_res_df_long_target_only_GPTcleanedup_less_than_24h[order(moon_res_df_long_target_only_GPTcleanedup_less_than_24h$mean_score_by_target, decreasing = T),]
# moon_res_df_long_target_only_GPTcleanedup_less_than_24h$source <- factor(moon_res_df_long_target_only_GPTcleanedup_less_than_24h$source, levels = unique(moon_res_df_long_target_only_GPTcleanedup_less_than_24h$source))
# 
# ggplot(moon_res_df_long_target_only_GPTcleanedup_less_than_24h, aes(x = source, y = score, group = source)) + 
#   geom_boxplot() + 
#   geom_jitter() + 
#   theme_minimal() +
#   theme(axis.text.x=element_text(angle=45, hjust=1)) + geom_hline(yintercept = 1.7)
# 
# experiements_id_less_than_24h <- unique(moon_res_df_long_target_only_GPTcleanedup[moon_res_df_long_target_only_GPTcleanedup$inferred_time_all_conditions <= 24,"id"])


ligand <- "IL6"

ligand_targets <- meta_network[meta_network$source == ligand,"target"]
ligand_targets <- ligand_targets[ligand_targets %in% moon_res_df$source_original]

moon_res_df_long_ligand_targets <- moon_res_df_long[which(moon_res_df_long$source %in% c(ligand,ligand_targets) & moon_res_df_long$target == ligand),]

ggplot(moon_res_df_long_ligand_targets, aes(x = source, y = score, group = source)) + 
  geom_boxplot() + 
  geom_jitter() + 
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + geom_hline(yintercept = 1.7)









