library(readr)
library(ggplot2)
library(ggrepel)

source("scripts/support_functions.R")

moon_res_df_long_target_only_GPTcleanedup <- as.data.frame(
  read_csv("results/moon_res_df_long_target_only_GPTcleanedup.csv"))
moon_res_df_long_target_only_GPTcleanedup <- moon_res_df_long_target_only_GPTcleanedup[!(moon_res_df_long_target_only_GPTcleanedup$source %in% c("ANGPT1","TSLP")),]
moon_res_df_long_target_only_GPTcleanedup <- moon_res_df_long_target_only_GPTcleanedup[!is.na(moon_res_df_long_target_only_GPTcleanedup$score),]

stat_list <- list()
mean_list <- list()
for(ligand in unique(moon_res_df_long_target_only_GPTcleanedup$source))
{
  ligand_df <- moon_res_df_long_target_only_GPTcleanedup[moon_res_df_long_target_only_GPTcleanedup$source == ligand,]
  if(length(na.omit(ligand_df$score)) > 1)
  {
    t_stat <- t.test(ligand_df$score)
    stat_list[[ligand]] <- t_stat$p.value
    mean_list[[ligand]] <- mean(ligand_df$score, na.rm = T)
  } else
  {
    stat_list[[ligand]] <- 1
    mean_list[[ligand]] <- na.omit(ligand_df$score)
  }
}
stats_df <- as.data.frame(t(data.frame(stat_list)))
stats_df$source <- row.names(stats_df)
names(stats_df)[1] <- "p_value"

mean_df <- as.data.frame(t(data.frame(mean_list)))
mean_df$source <- row.names(mean_df)
names(mean_df)[1] <- "mean_score"

moon_res_df_long_target_only_GPTcleanedup <- merge(moon_res_df_long_target_only_GPTcleanedup,stats_df, by = "source")
moon_res_df_long_target_only_GPTcleanedup <- merge(moon_res_df_long_target_only_GPTcleanedup,mean_df, by = "source")

moon_res_df_long_target_only_GPTcleanedup$significant <- ifelse(moon_res_df_long_target_only_GPTcleanedup$p_value <= 0.05,"<=0.05",">0.05")
moon_res_df_long_target_only_GPTcleanedup <- moon_res_df_long_target_only_GPTcleanedup[order(moon_res_df_long_target_only_GPTcleanedup$mean_score, decreasing = T),]

moon_res_df_long_target_only_GPTcleanedup$source<- factor(moon_res_df_long_target_only_GPTcleanedup$source, levels = unique(moon_res_df_long_target_only_GPTcleanedup$source))

ggplot(moon_res_df_long_target_only_GPTcleanedup, aes(x = source, y = score, group = source, fill = significant, color = significant)) + 
  geom_boxplot(coef = 6) + 
  geom_jitter(alpha = 0.5) + 
  theme_minimal() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) + geom_hline(yintercept = 0) + scale_fill_manual(values = c("skyblue","lightgrey")) + scale_color_manual(values = c("darkblue","grey")) 

mean_pvalue_df <- merge(mean_df, stats_df)

ggplot(mean_pvalue_df, aes(x = mean_score, y = -log10(p_value), label = source)) + 
  geom_point() + 
  geom_label_repel() +
  geom_hline(yintercept = 1.3) + 
  xlim(-max(abs(mean_pvalue_df$mean_score)), max(abs(mean_pvalue_df$mean_score))) +
  theme_minimal()
