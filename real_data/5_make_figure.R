rm(list = ls())
require(data.table)
require(tidyverse)
library(latex2exp)
library(cowplot)

tissue_vector <- fread("./tissue.txt", header = F)$V1
N_PP <- c()
N_classic <- c()
N_PopInf <- c()
N_eff <- c()
df_cb <- c()
for (tissue in tissue_vector){
  print(tissue)
  # tissue <- "Artery_Aorta"
  df <- fread(paste0("./SexBSE/", tissue, ".txt.gz"))
  
  corr <- fread(paste0("./corr/", tissue, ".txt.gz"))
  df$corr <- corr$cor
  df$se_PP_ratio <- df$se_vec_PP/df$se_vec_classic
  df$se_PopInf_ratio <- df$se_vec_PopInf/df$se_vec_classic
  df$se_eff_ratio <- df$se_vec_eff/df$se_vec_classic
  df$Name <- c(1:nrow(df))
  print(head(df))
  df_cb <- rbind(df_cb, df)
  print(head(df_cb))
  
  # Adjust p-values outside the loop
  q_vec_PP <- p.adjust(df$P_vec_PP, method = "BH")
  q_vec_PopInf <- p.adjust(df$P_vec_PopInf, method = "BH")
  q_vec_classic <- p.adjust(df$P_vec_classic, method = "BH")
  q_vec_eff <- p.adjust(df$P_vec_eff, method = "BH")
  N_PP <- c(N_PP, sum(q_vec_PP < 0.05))
  N_PopInf <- c(N_PopInf, sum(q_vec_PopInf < 0.05))
  N_classic <- c(N_classic, sum(q_vec_classic < 0.05))
  N_eff <- c(N_eff, sum(q_vec_eff < 0.05))
  
}

N_summ <- data.frame(N_PP, N_PopInf, N_classic, N_eff)
df_cb$q_vec_PP <- p.adjust(df_cb$P_vec_PP, method = "BH")
df_cb$q_vec_PopInf <- p.adjust(df_cb$P_vec_PopInf, method = "BH")
df_cb$q_vec_classic <- p.adjust(df_cb$P_vec_classic, method = "BH")
df_cb$q_vec_eff <- p.adjust(df_cb$P_vec_eff, method = "BH")

N_PP_unique <- length(unique(df_cb$Name[df_cb$q_vec_PP < 0.05]))
N_PopInf_unique <- length(unique(df_cb$Name[df_cb$q_vec_PopInf < 0.05]))
N_classic_unique <- length(unique(df_cb$Name[df_cb$q_vec_classic < 0.05]))
N_eff_unique <- length(unique(df_cb$Name[df_cb$q_vec_eff < 0.05]))
df_N_genes <- data.frame(N <- c(N_PP_unique, N_PopInf_unique, N_classic_unique, N_eff_unique),
                         Method = c("PP", "POP-Inf", "Classic", "EIF*-based"))


df_plot <- data.frame(Classic = df_cb$beta_vec_classic, PP = df_cb$beta_vec_PP)
beta_pp <- ggplot(df_plot, aes(x = Classic, y = PP)) + 
  geom_point(size = 0.5, color = "#7fc97f") +                           
  geom_abline(slope = 1, intercept = 0, color = "#763262", linetype = "dashed") +
  theme_minimal() +
  xlim(-1.6, 1.6) +
  ylim(-1.6, 1.6) +
  labs(x = TeX("Classic $\\hat{\\theta}$"), y = TeX("PP $\\hat{\\theta}$")) +
  theme(axis.text.x = element_text(family = "Helvetica", size = 12),
        axis.text.y = element_text(family = "Helvetica", size = 12),
        axis.title.x = element_text(vjust = 1),
        axis.title.y = element_text(vjust = 2),
        text = element_text(family = "Helvetica", size=12, color = "black"),
        legend.position = ""
  )


df_plot <- data.frame(Classic = df_cb$beta_vec_classic, eff = df_cb$beta_vec_eff)
beta_eff <- ggplot(df_plot, aes(x = Classic, y = eff)) + 
  geom_point(size = 0.5, color = "#984EA3") +                           
  geom_abline(slope = 1, intercept = 0, color = "#763262", linetype = "dashed") +
  theme_minimal() +
  xlim(-1.6, 1.6) +
  ylim(-1.6, 1.6) +
  labs(x = TeX("Classic $\\hat{\\theta}$"), y = TeX("EIF*-based $\\hat{\\theta}$")) +
  theme(axis.text.x = element_text(family = "Helvetica", size = 12),
        axis.text.y = element_text(family = "Helvetica", size = 12),
        axis.title.x = element_text(vjust = 1),
        axis.title.y = element_text(vjust = 2),
        text = element_text(family = "Helvetica", size=12, color = "black"),
        legend.position = ""
  )


df_plot <- data.frame(Classic = df_cb$beta_vec_classic, PopInf = df_cb$beta_vec_PopInf)
beta_popinf <- ggplot(df_plot, aes(x = Classic, y = PopInf)) + 
  geom_point(size = 0.5, color = "#EF3B2C") +   
  geom_abline(slope = 1, intercept = 0, color = "#763262", linetype = "dashed") +
  theme_minimal() +
  xlim(-1.6, 1.6) +
  ylim(-1.6, 1.6) +
  labs(x = TeX("Classic $\\hat{\\theta}$"), y = TeX("POP-Inf $\\hat{\\theta}$")) +
  theme(axis.text.x = element_text(family = "Helvetica", size = 12),
        axis.text.y = element_text(family = "Helvetica", size = 12),
        axis.title.x = element_text(vjust = 1),
        axis.title.y = element_text(vjust = 2),
        text = element_text(family = "Helvetica", size=12, color = "black"),
        legend.position = ""
  )

# Create scatter plot with fitted line
se_pp <- ggplot(df_cb, aes(x = abs(df_cb$corr), y = df_cb$se_PP_ratio)) + 
  geom_point(size = 0.5, color = "#7fc97f") +                           
  geom_hline(yintercept = 1, linetype = "dashed", color = "#763262") + # Add horizontal line at Y = 1
  theme_minimal() +
  xlim(0, 0.8) +
  ylim(0.6, 2.3) +
  labs(y = expression(paste(hat("SE"),"(PP) /  ", hat("SE"), "(Classic)"
                            , sep = "")
  ),
  x =  "|Imputation correlation|"
  ) +
  theme(axis.text.x = element_text(family = "Helvetica", size = 12),
        axis.text.y = element_text(family = "Helvetica", size = 12),
        axis.title.x = element_text(vjust = 1),
        axis.title.y = element_text(vjust = 2, size = 8),
        text = element_text(family = "Helvetica", size=12, color = "black"),
        legend.position = ""
  )

# Create scatter plot with fitted line
se_eff <- ggplot(df_cb, aes(x = abs(df_cb$corr), y = df_cb$se_eff_ratio)) + 
  geom_point(size = 0.5, color = "#984EA3") +                          
  geom_hline(yintercept = 1, linetype = "dashed", color = "#763262") + # Add horizontal line at Y = 1
  theme_minimal() +
  xlim(0, 0.8) +
  ylim(0.6, 2.3) +
  labs(y = expression(paste(hat("SE"),"(EIF*-based) /  ", hat("SE"), "(Classic)"
                            , sep = "")
  ),
  x =  "|Imputation correlation|"
  ) +
  theme(axis.text.x = element_text(family = "Helvetica", size = 12),
        axis.text.y = element_text(family = "Helvetica", size = 12),
        axis.title.x = element_text(vjust = 1),
        axis.title.y = element_text(vjust = 2, size = 7),
        text = element_text(family = "Helvetica", size=12, color = "black"),
        legend.position = ""
  )

# Create scatter plot with fitted line
se_popinf <- ggplot(df_cb, aes(x = abs(df_cb$corr), y = df_cb$se_PopInf_ratio)) + 
  geom_point(size = 0.5, color = "#EF3B2C") +                          
  geom_hline(yintercept = 1, linetype = "dashed", color = "#763262") + # Add horizontal line at Y = 1
  theme_minimal() +
  xlim(0, 0.8) +
  ylim(0.6, 2.3) +
  labs(y = expression(paste(hat("SE"),"(POP-Inf) /  ", hat("SE"), "(Classic)"
                            , sep = "")
  ),
  x =  "|Imputation correlation|"
  ) +
  theme(axis.text.x = element_text(family = "Helvetica", size = 12),
        axis.text.y = element_text(family = "Helvetica", size = 12),
        axis.title.x = element_text(vjust = 1),
        axis.title.y = element_text(vjust = 2, size = 8),
        text = element_text(family = "Helvetica", size=12, color = "black"),
        legend.position = ""
  )


colnames(df_N_genes) <- c("N", "Method")
df_N_genes$Method <- factor(df_N_genes$Method, levels = c("POP-Inf", "EIF*-based", "Classic", "PP"))

N_genes_p <- ggplot(df_N_genes, aes(x=as.factor(Method), y = N, fill=as.factor(Method))) +
  geom_bar(stat="identity") +
  theme_classic() +
  coord_flip() +
  scale_fill_manual(values = c("#EF3B2C", "#984EA3", "#fdb462", "#7fc97f")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Number of sex−biased genes",
       x =  "Methods") +
  theme(axis.text.x = element_text(family = "Helvetica", size = 12),
        axis.text.y = element_text(family = "Helvetica", size = 12),
        axis.title.x = element_text(vjust = 1),
        axis.title.y = element_text(vjust = 2),
        text = element_text(family = "Helvetica", size=12, color = "black"),
        legend.position = ""
  )


p_top <- plot_grid(beta_pp, beta_eff, beta_popinf, se_pp, se_eff, se_popinf, labels = c('a', 'b', 'c', 'd', 'e', 'f'), nrow = 2, ncol = 3, label_size = 14)

bottom_row <- plot_grid(N_genes_p, labels = c('g'), label_size = 14)
p_all <- plot_grid(p_top, bottom_row, labels = c('', ''), rel_heights = c(2, 1), label_size = 14, ncol = 1)

jpeg_out_dir <-  paste0("./cow.jpeg")
jpeg(jpeg_out_dir, width = 7, height = 7, units = 'in', res = 300, quality=100)
par(mar=c(5,5,4,4)+.1)
print(p_all)
dev.off()