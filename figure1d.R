library(tidyverse)
library(openxlsx)

setwd("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/")
load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/prediction/prediction_rf_lodo.RData")
colnames(dat_comb)[1:300] # 30:224, species

# PCoA data after batch correction ----
library(vegan)
table(dat_comb$study)
distance.matrix <- vegdist(dat_comb[, 30:224], method="bray")

## do the MDS math (this is basically eigen value decomposition)
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE, k=20)

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- mds.stuff$points
mds.data <- data.frame(id=dat_comb$id,
                       X=mds.values[,1],
                       Y=mds.values[,2])
mds.data <- merge(dat_comb, mds.data, by="id")

met888 <- read.xlsx("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/mets_1090.xlsx",
                    sheet="met888")

# correlation between PCs and mets
PCmet_cor <- data.frame(CHEM_ID=met888$CHEM_ID,
                        corX=NA,
                        corY=NA)
PCmet_cor <- merge(met888, PCmet_cor, by="CHEM_ID")
for (i in 1:nrow(PCmet_cor)){
  PCmet_cor[i, "corX"]=cor(mds.data$X, mds.data[, PCmet_cor$CHEM_ID[i]], method="spearman")
  PCmet_cor[i, "corY"]=cor(mds.data$Y, mds.data[, PCmet_cor$CHEM_ID[i]], method="spearman")
}
PCmet_cor$add <- abs(PCmet_cor$corX)+abs(PCmet_cor$corY)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/mbx_t2d/mets_t2d_meta_combine.RData")
View(PCmet_cor %>% filter(CHEM_ID %in% meta_mets_t2d_sig$CHEM_ID))

addmet <- PCmet_cor %>% 
  filter(mets_name %in% c("cinnamoylglycine",
                          "imidazole propionate",
                          "phenylacetylglutamine",
                          "propionylcarnitine",
                          "3-ureidopropionate",
                          "acetylcarnitine",
                          "glycine",
                          "isoleucine",
                          "trigonelline",
                          "urate"))

library(ggplot2)
library(scales) 
library(ggsci)
library(ggrepel)
library(ggnewscale)

(col_pal2=pal_nejm("default",alpha=1)(8))
show_col(col_pal2)

col_super_path <- 
  c(
    "Lipids" = col_pal2[1],
    "Amino acids" = col_pal2[2],
    "Xenobiotics" = col_pal2[3],
    "Peptides" = col_pal2[4],
    "Nucleotides" = col_pal2[5],
    "Cofactors and vitamins" = col_pal2[6],
    "Carbohydrates and energy metabolism" = col_pal2[8],
    "Other metabolites"="grey70")

mds.data$BFratio <- mds.data$Bacteroidetes/mds.data$Firmicutes
mds.data$BFratio <- ifelse(mds.data$BFratio>10, 10, mds.data$BFratio)

pdf(file = "/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/figure1d.pdf",
    width = 8.4, height = 5)
ggplot(data=mds.data, aes(x=-X, y=Y)) +
  geom_point(size=3, aes(color=BFratio), shape=19, alpha=0.5) +
  xlab(paste("PCo1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("PCo2 - ", mds.var.per[2], "%", sep="")) +
  scale_x_continuous(limits = c(-0.35, 0.7), breaks = seq(-0.2, 0.6,by=0.2),
                     labels = label_number(accuracy = 0.1)) +
  scale_y_continuous(limits = c(-0.55, 0.4), breaks = seq(-0.4, 0.4,by=0.2),
                     labels = label_number(accuracy = 0.1)) +
  # scale_fill_gradient(low = "#56B1F7", high = "#132B43", limits=c(0,10), breaks=c(0,5,10)) +
  scale_color_viridis_c(direction=-1, limits=c(0,10), breaks=c(0,5,10)) +
  labs(color="Bacteroidetes/Firmicutes ratio") +
  new_scale_color() +
  geom_segment(aes(x = 0, y = 0, xend = corX*1.25, yend = corY*1.25,
                   color = SUPER_PATHWAY2),
               data = addmet,
               linewidth = 0.7,
               arrow = arrow(length = unit(0.3, "cm"))) +
  scale_color_manual(values = col_super_path) +
  geom_label_repel(aes(x = corX*1.25, y = corY*1.25, label = mets_name,
                 color = SUPER_PATHWAY2),
             data = addmet,
             size = 4,
             # vjust = 1,  # Adjust the vertical position of the label
             # hjust = 0.5,   # Adjust the horizontal position of the label
             show.legend = FALSE
  ) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, linewidth=1),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16, colour = "black"),
    legend.text = element_text(size = 16, colour = "black")
  )
dev.off()

save.image(file="/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/figure1d.RData")
