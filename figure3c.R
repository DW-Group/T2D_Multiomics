library(tidyverse)
library(data.table)

setwd("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig3")

# Model 1
load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/prediction/species_for_RF_T2D_ROC_basic2.RData")
auc_No_pool <- data.frame(group=gsub("cardis", "pool", auc_No_cardis$group),
                          study="Averaged",
                          auc=rowMeans(data.frame(auc1=auc_No_cardis$auc,
                                                  auc2=auc_No_direct$auc,
                                                  auc3=auc_No_hpfs$auc,
                                                  auc4=auc_No_nhs2$auc,
                                                  auc5=auc_No_sol$auc)))
auc_Yes_pool <- data.frame(group=gsub("cardis", "pool", auc_Yes_cardis$group),
                           study="Averaged",
                           auc=rowMeans(data.frame(auc1=auc_Yes_cardis$auc,
                                                   auc2=auc_Yes_direct$auc,
                                                   auc3=auc_Yes_hpfs$auc,
                                                   auc4=auc_Yes_nhs2$auc,
                                                   auc5=auc_Yes_sol$auc)))

auc_No_cardis$study <- "MetaCardis"
auc_No_direct$study <- "DIRECT-PLUS"
auc_No_hpfs$study <- "HPFS"
auc_No_nhs2$study <- "NHSII"
auc_No_sol$study <- "HCHS/SOL"

auc_Yes_cardis$study <- "MetaCardis"
auc_Yes_direct$study <- "DIRECT-PLUS"
auc_Yes_hpfs$study <- "HPFS"
auc_Yes_nhs2$study <- "NHSII"
auc_Yes_sol$study <- "HCHS/SOL"

auc_combine_long1 <- bind_rows(list(auc_No_cardis[, c(1,5,2)],
                                    auc_No_direct[, c(1,5,2)],
                                    auc_No_hpfs[, c(1,5,2)],
                                    auc_No_nhs2[, c(1,5,2)],
                                    auc_No_sol[, c(1,5,2)],
                                    auc_No_pool,
                                    
                                    auc_Yes_cardis[, c(1,5,2)],
                                    auc_Yes_direct[, c(1,5,2)],
                                    auc_Yes_hpfs[, c(1,5,2)],
                                    auc_Yes_nhs2[, c(1,5,2)],
                                    auc_Yes_sol[, c(1,5,2)],
                                    auc_Yes_pool))
auc_combine_long1$metf <- c(rep("Metformin nonuser", 600), 
                            rep("Metformin user", 600))
auc_combine_long1$Model <- "Model1"
gdata::keep(auc_combine_long1, sure=T)

# Model 2
load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/prediction/species_for_RF_T2D_ROC_spe2.RData")
auc_No_pool <- data.frame(group=gsub("cardis", "pool", auc_No_cardis$group),
                          study="Averaged",
                          auc=rowMeans(data.frame(auc1=auc_No_cardis$auc,
                                                  auc2=auc_No_direct$auc,
                                                  auc3=auc_No_hpfs$auc,
                                                  auc4=auc_No_nhs2$auc,
                                                  auc5=auc_No_sol$auc)))
auc_Yes_pool <- data.frame(group=gsub("cardis", "pool", auc_Yes_cardis$group),
                           study="Averaged",
                           auc=rowMeans(data.frame(auc1=auc_Yes_cardis$auc,
                                                   auc2=auc_Yes_direct$auc,
                                                   auc3=auc_Yes_hpfs$auc,
                                                   auc4=auc_Yes_nhs2$auc,
                                                   auc5=auc_Yes_sol$auc)))

auc_No_cardis$study <- "MetaCardis"
auc_No_direct$study <- "DIRECT-PLUS"
auc_No_hpfs$study <- "HPFS"
auc_No_nhs2$study <- "NHSII"
auc_No_sol$study <- "HCHS/SOL"

auc_Yes_cardis$study <- "MetaCardis"
auc_Yes_direct$study <- "DIRECT-PLUS"
auc_Yes_hpfs$study <- "HPFS"
auc_Yes_nhs2$study <- "NHSII"
auc_Yes_sol$study <- "HCHS/SOL"

auc_combine_long2 <- bind_rows(list(auc_No_cardis[, c(1,5,2)],
                                    auc_No_direct[, c(1,5,2)],
                                    auc_No_hpfs[, c(1,5,2)],
                                    auc_No_nhs2[, c(1,5,2)],
                                    auc_No_sol[, c(1,5,2)],
                                    auc_No_pool,
                                    
                                    auc_Yes_cardis[, c(1,5,2)],
                                    auc_Yes_direct[, c(1,5,2)],
                                    auc_Yes_hpfs[, c(1,5,2)],
                                    auc_Yes_nhs2[, c(1,5,2)],
                                    auc_Yes_sol[, c(1,5,2)],
                                    auc_Yes_pool))
auc_combine_long2$metf <- c(rep("Metformin nonuser", 600), 
                            rep("Metformin user", 600))
auc_combine_long2$Model <- "Model2"
gdata::keep(auc_combine_long1, auc_combine_long2, sure=T)

# Model 3
load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/prediction/species_for_RF_T2D_ROC_met2.RData")
auc_No_pool <- data.frame(group=gsub("cardis", "pool", auc_No_cardis$group),
                          study="Averaged",
                          auc=rowMeans(data.frame(auc1=auc_No_cardis$auc,
                                                  auc2=auc_No_direct$auc,
                                                  auc3=auc_No_hpfs$auc,
                                                  auc4=auc_No_nhs2$auc,
                                                  auc5=auc_No_sol$auc)))
auc_Yes_pool <- data.frame(group=gsub("cardis", "pool", auc_Yes_cardis$group),
                           study="Averaged",
                           auc=rowMeans(data.frame(auc1=auc_Yes_cardis$auc,
                                                   auc2=auc_Yes_direct$auc,
                                                   auc3=auc_Yes_hpfs$auc,
                                                   auc4=auc_Yes_nhs2$auc,
                                                   auc5=auc_Yes_sol$auc)))

auc_No_cardis$study <- "MetaCardis"
auc_No_direct$study <- "DIRECT-PLUS"
auc_No_hpfs$study <- "HPFS"
auc_No_nhs2$study <- "NHSII"
auc_No_sol$study <- "HCHS/SOL"

auc_Yes_cardis$study <- "MetaCardis"
auc_Yes_direct$study <- "DIRECT-PLUS"
auc_Yes_hpfs$study <- "HPFS"
auc_Yes_nhs2$study <- "NHSII"
auc_Yes_sol$study <- "HCHS/SOL"

auc_combine_long3 <- bind_rows(list(auc_No_cardis[, c(1,5,2)],
                                    auc_No_direct[, c(1,5,2)],
                                    auc_No_hpfs[, c(1,5,2)],
                                    auc_No_nhs2[, c(1,5,2)],
                                    auc_No_sol[, c(1,5,2)],
                                    auc_No_pool,
                                    
                                    auc_Yes_cardis[, c(1,5,2)],
                                    auc_Yes_direct[, c(1,5,2)],
                                    auc_Yes_hpfs[, c(1,5,2)],
                                    auc_Yes_nhs2[, c(1,5,2)],
                                    auc_Yes_sol[, c(1,5,2)],
                                    auc_Yes_pool))
auc_combine_long3$metf <- c(rep("Metformin nonuser", 600), 
                            rep("Metformin user", 600))
auc_combine_long3$Model <- "Model3"
gdata::keep(auc_combine_long1, auc_combine_long2, auc_combine_long3, sure=T)

# Model 4
load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/prediction/species_for_RF_T2D_ROC_spemet2.RData")
auc_No_pool <- data.frame(group=gsub("cardis", "pool", auc_No_cardis$group),
                          study="Averaged",
                          auc=rowMeans(data.frame(auc1=auc_No_cardis$auc,
                                                  auc2=auc_No_direct$auc,
                                                  auc3=auc_No_hpfs$auc,
                                                  auc4=auc_No_nhs2$auc,
                                                  auc5=auc_No_sol$auc)))
auc_Yes_pool <- data.frame(group=gsub("cardis", "pool", auc_Yes_cardis$group),
                           study="Averaged",
                           auc=rowMeans(data.frame(auc1=auc_Yes_cardis$auc,
                                                   auc2=auc_Yes_direct$auc,
                                                   auc3=auc_Yes_hpfs$auc,
                                                   auc4=auc_Yes_nhs2$auc,
                                                   auc5=auc_Yes_sol$auc)))

auc_No_cardis$study <- "MetaCardis"
auc_No_direct$study <- "DIRECT-PLUS"
auc_No_hpfs$study <- "HPFS"
auc_No_nhs2$study <- "NHSII"
auc_No_sol$study <- "HCHS/SOL"

auc_Yes_cardis$study <- "MetaCardis"
auc_Yes_direct$study <- "DIRECT-PLUS"
auc_Yes_hpfs$study <- "HPFS"
auc_Yes_nhs2$study <- "NHSII"
auc_Yes_sol$study <- "HCHS/SOL"

auc_combine_long4 <- bind_rows(list(auc_No_cardis[, c(1,5,2)],
                                    auc_No_direct[, c(1,5,2)],
                                    auc_No_hpfs[, c(1,5,2)],
                                    auc_No_nhs2[, c(1,5,2)],
                                    auc_No_sol[, c(1,5,2)],
                                    auc_No_pool,
                                    
                                    auc_Yes_cardis[, c(1,5,2)],
                                    auc_Yes_direct[, c(1,5,2)],
                                    auc_Yes_hpfs[, c(1,5,2)],
                                    auc_Yes_nhs2[, c(1,5,2)],
                                    auc_Yes_sol[, c(1,5,2)],
                                    auc_Yes_pool))
auc_combine_long4$metf <- c(rep("Metformin nonuser", 600), 
                            rep("Metformin user", 600))
auc_combine_long4$Model <- "Model4"
gdata::keep(auc_combine_long1, auc_combine_long2, 
            auc_combine_long3, auc_combine_long4, sure=T)

##############
auc_combine_long <- bind_rows(list(auc_combine_long1, auc_combine_long2, 
                                   auc_combine_long3, auc_combine_long4))

View(auc_combine_long %>% group_by(metf, study, Model) %>% 
       summarise(mean=mean(auc)))

auc_combine_long$metf <- 
  factor(auc_combine_long$metf,
         levels = c("Metformin user","Metformin nonuser"))

pal_sel <- paletteer::paletteer_d("ggsci::category20b_d3")
study_cols <- as.vector(pal_sel[1:4])
names(study_cols) <- 
  c("Model1","Model2","Model3","Model4")

auc_combine_long$study <- 
  factor(auc_combine_long$study,
         levels = c("Averaged", "DIRECT-PLUS", "HCHS/SOL", "HPFS", 
                    "MetaCardis", "NHSII"))

metf_group_names <- c(
  `Metformin user` = "Controls vs. Meformin-treatment T2D",
  `Metformin nonuser` = "Controls vs. Meformin-naive T2D"
)

library(ggplot2)
ggplot(auc_combine_long,aes(x=study,y=auc))+
  geom_boxplot(aes(color=Model),outlier.shape = NA)+
  geom_point(aes(color=Model),alpha=0.2,size = 0.15,
             position = position_jitterdodge(jitter.width = 0.1,
                                             jitter.height = 0.01,
                                             dodge.width = 0.8,
                                             seed = 1)
  )+
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), limits = c(0.2,1.01))+
  scale_color_manual(values = scales::alpha(study_cols,0.8),
                     labels=c("age+sex+BMI","age+sex+BMI+species",
                              "age+sex+BMI+metabolites",
                              "age+sex+BMI+species+metabolites"),
                     name="Models")+
  facet_grid(~metf,scales = "free",space = "free",
             labeller = as_labeller(metf_group_names))+
  theme_classic()+
  labs(y="AUC from random forest with \nleave-one-dataset-out cross-validation")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size=11,angle = 45,hjust = 1,vjust = 1,color = "black"),
        axis.text.y = element_text(size = 11,color = "black"),
        legend.position = "none",
        strip.text = element_text(size = 12)
  )
ggsave("./figure3c.pdf", width = 12,height = 4)

p_lgd <- ggplot(auc_combine_long,aes(x=study,y=auc))+
  geom_boxplot(aes(color=Model),outlier.shape = NA)+
  geom_point(aes(color=Model),alpha=0.2,size = 0.15,
             position = position_jitterdodge(jitter.width = 0.1,
                                             jitter.height = 0.01,
                                             dodge.width = 0.8,
                                             seed = 1)
  )+
  scale_color_manual(values = alpha(study_cols,0.8),
                     labels=c("age+sex+BMI","age+sex+BMI+species",
                              "age+sex+BMI+metabolites",
                              "age+sex+BMI+species+metabolites"),
                     name="Model")+
  facet_grid(~metf,scales = "free",space = "free",
             labeller = as_labeller(metf_group_names))+
  guides(color=guide_legend(ncol=4,title.position="left"))+
  theme_classic()+
  labs(y="AUC from random forest with \nleave-one-dataset-out cross-validation")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size=11,angle = 45,hjust = 1,vjust = 1,color = "black"),
        axis.text.y = element_text(size = 11,color = "black"))
p_lgd
lgd <- ggpubr::get_legend(p_lgd)

library(grid)
pdf("./figure3c_lgd.pdf",
    height = 0.5,width = 8)
grid.draw(lgd)
dev.off()