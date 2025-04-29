library(tidyverse)
load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/1.preprocess/pool_species_20250115_pcoa.RData")

# Fig. PCoA------
species_zero_adj0 <- sweep(species_zero_adj0, 2, colSums(species_zero_adj0), FUN="/")
summary(colSums(species_zero_adj0)) # should be all 1

species_all_zero <- sweep(species_all_zero, 2, colSums(species_all_zero), FUN="/")
summary(colSums(species_all_zero)) # should be all 1

# PCoA data before batch correction ----
library(vegan)
metadata3317 <- metadata3317 %>% 
  mutate(study=case_when(study=="MBS" ~ "NHSII",
                         study=="MLVS" ~ "HPFS",
                         study=="SOL" ~ "HCHS/SOL",
                         study=="Segal1" ~ "Talmor-Barkan_2022",
                         TRUE ~ study
  ))

table(metadata3317$study)
identical(colnames(species_all_zero), metadata3317$mgx_id) # TRUE
distance.matrix <- vegdist(t(species_all_zero), method="bray")

## do the MDS math (this is basically eigen value decomposition)
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE, k=20)

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- mds.stuff$points
mds.data <- data.frame(id=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2])
identical(mds.data$id, metadata3317$mgx_id) # TRUE
mds.data <- mds.data %>% 
  inner_join(metadata3317,by=c("id"="mgx_id"))

# PCoA data after batch correction  ------
species_zero_adj0 <- species_zero_adj0[, metadata3317$mgx_id]
identical(colnames(species_zero_adj0), metadata3317$mgx_id) # TRUE
distance.matrix2 <- vegdist(t(species_zero_adj0), method="bray")

## do the MDS math (this is basically eigen value decomposition)
mds.stuff2 <- cmdscale(distance.matrix2, eig=TRUE, x.ret=TRUE, k=20)

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per2 <- round(mds.stuff2$eig/sum(mds.stuff2$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values2 <- mds.stuff2$points
mds.data2 <- data.frame(id=rownames(mds.values2),
                        X=mds.values2[,1],
                        Y=mds.values2[,2])
identical(mds.data2$id, metadata3317$mgx_id) # TRUE
mds.data2 <- mds.data2 %>% 
  inner_join(metadata3317,by=c("id"="mgx_id"))

set.seed(1234)
R2_before <- adonis2(distance.matrix ~ study, data=metadata3317) # R2=0.0796
R2_before$R2
set.seed(1234)
R2_after <- adonis2(distance.matrix2 ~ study, data=metadata3317) # R2=0.0381
R2_after$R2

save.image(file="/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/figS2.RData")

library(ggsci)
(col_pal <- pal_nejm(palette = "default",alpha = 0.7)(6))
scales::show_col(col_pal)

studycol <- c(
  `DIRECT-PLUS`= col_pal[6],
  `HCHS/SOL` = col_pal[5],
  NHSII = col_pal[4],
  HPFS = col_pal[3], 
  MetaCardis = col_pal[2],
  `Talmor-Barkan_2022` = col_pal[1]
)

# draw PCoA plot
mds.data$study <- factor(mds.data$study, levels=c("DIRECT-PLUS", "HCHS/SOL", "HPFS", 
                                                  "MetaCardis", "NHSII", "Talmor-Barkan_2022"))
mds.data2$study <- factor(mds.data2$study, levels=c("DIRECT-PLUS", "HCHS/SOL", "HPFS", 
                                                    "MetaCardis", "NHSII", "Talmor-Barkan_2022"))

p_study1 <- ggplot(data=mds.data, aes(x=X, y=Y)) +
  # geom_point(size=3.5, aes(color=study), shape=19, alpha=0.3) +
  geom_point(size=3.5, aes(fill=study), shape=21) +
  scale_fill_manual(values=studycol) +
  xlab(paste("PCo1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("PCo2 - ", mds.var.per[2], "%", sep="")) +
  labs(title="") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5, size = 14, colour = "black"),
        legend.position="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size = 11, colour = "black"),
        legend.spacing.x = unit(0.1, 'cm')) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, override.aes = list(size = 4))) +
  annotate(geom="text", x=0.4, y=0.36, col="black", size=4.5,
           label="Before correction, R2=8.0%")

p_study2 <- ggplot(data=mds.data2, aes(x=X, y=-Y)) +
  # geom_point(size=3.5, aes(color=study), shape=19, alpha=0.3) +
  geom_point(size=3.5, aes(fill=study), shape=21) +
  scale_fill_manual(values=studycol) +
  xlab(paste("PCo1 - ", mds.var.per2[1], "%", sep="")) +
  ylab(paste("PCo2 - ", mds.var.per2[2], "%", sep="")) +
  labs(title="") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(size = 14, colour = "black"),
        axis.text.y = element_text(size = 14, colour = "black"),
        axis.title.x = element_text(size = 14, colour = "black"),
        axis.title.y = element_text(size = 14, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        plot.title = element_text(hjust = 0.5, size = 14, colour = "black"),
        legend.position = "none") +
  annotate(geom="text", x=0.44, y=0.32, col="black", size=4.5,
           label="After correction, R2=3.8%")

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/figS2.pdf",
    width = 10, height = 5.6, onefile = F) # Open a new pdf file
ggpubr::ggarrange(p_study1, p_study2,
                  ncol=2, legend = "bottom", common.legend = T)
dev.off() # Close the file














# ## species and PCoA correlation
# dim(mds.data2) #5465/30
# 
# spe_df <- rotate_df(species_zero_adj0) %>% 
#   rownames_to_column(var = "id")
# dim(spe_df) #5465/196
# 
# mds_spe <- mds.data2 %>% 
#   inner_join(spe_df,by="id")
# dim(mds_spe) #5465/223
# 
# cor.test(mds_spe$X,mds_spe$s__Prevotella_copri, method="spearman") # 0.78
# cor.test(mds_spe$Y,mds_spe$s__Prevotella_copri, method="spearman") # 0.08
# 
# pco1_cor <- c()
# for (i in 29:223){
#   pco1_cor <- rbind(pco1_cor,c(colnames(mds_spe)[i],
#                                cor(mds_spe$X,mds_spe[,i],method = "spearman")))
# }
# pco1_cor <- as.data.frame(pco1_cor) %>% 
#   setnames(c("species","pco1_cor"))
# 
# pco2_cor <- c()
# for (i in 29:223){
#   pco2_cor <- rbind(pco2_cor,c(colnames(mds_spe)[i],
#                                cor(mds_spe$Y,mds_spe[,i],method = "spearman")))
# }
# pco2_cor <- as.data.frame(pco2_cor) %>% 
#   setnames(c("species","pco2_cor"))
# 
# spe_pcoa_cor <- pco1_cor %>% 
#   inner_join(pco2_cor,by="species")
# spe_pcoa_cor[,2:3] <- map_df(spe_pcoa_cor[,2:3],as.numeric)
# write.xlsx(spe_pcoa_cor, "./pool_species_20250115_spe_pcoa_cor.xlsx")
