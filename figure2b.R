library(tidyverse)
library(openxlsx)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/prediction/prediction_rf_lodo.RData")
gdata::keep(met806, sure=T)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/figure1bc.RData")
mets01_t$present <- rowSums(mets01_t[, 2:7])

met93 <- intersect((mets01_t %>% filter(present==6) %>% pull(CHEM_ID)), met806)
gdata::keep(met93, sure=T)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/prediction/prediction_rf_lodo.RData")
colnames(dat_comb)[1:300] # 30:224, species

# PCoA data after batch correction ----
library(vegan)
met_all <- dat_comb[, met93]
met_all[is.na(met_all)] <- 0

# permanova
set.seed(1234)
pp <- adonis2(met_all ~ dat_comb$status_new, method = "euclidean")
# 1.11%

# PCoA plot
distance.matrix <- vegdist(met_all, method = "euclidean")

## do the MDS math (this is basically eigen value decomposition)
mds.stuff <- cmdscale(distance.matrix, eig=TRUE, x.ret=TRUE, k=20)

## calculate the percentage of variation that each MDS axis accounts for...
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)

## now make a fancy looking plot that shows the MDS axes and the variation:
mds.values <- mds.stuff$points
mds.data <- data.frame(id=dat_comb$id,
                       X=mds.values[,1],
                       Y=mds.values[,2])
mds.data <- merge(dat_comb, mds.data, by="id") %>% 
  mutate(status_new=case_when(status_new=="Con" ~ 'Control',
                              status_new=="Pre" ~ 'Prediabetes',
                              status_new=="T2D" ~ 'T2D'))

library(ggplot2)
library(scales) 
library(ggsci)
statuscol <- c("#0a74b2","#e18726","#b93f2b")
names(statuscol) <- c("Control", "Prediabetes", "T2D")

save.image(file="/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig2/figure2b.Rdata")

pdf(file = "/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig2/figure2b.pdf",
    width = 7, height = 5)
ggplot(data=mds.data, aes(x=X, y=Y)) +
  geom_point(size=4, aes(fill=status_new), alpha=0.7, shape=21) +
  xlab(paste("PCo1 - ", mds.var.per[1], "%", sep="")) +
  ylab(paste("PCo2 - ", mds.var.per[2], "%", sep="")) +
  scale_x_continuous(limits = c(-16.5, 12.5), breaks = seq(-15, 10, by=5)) +
  scale_y_continuous(limits = c(-12, 9.6), breaks = seq(-12, 8, by=4),
                     labels = label_number(accuracy = 0.1)) +
  labs(title="PCoA plot by diabetes status") +
  scale_fill_manual(values = statuscol,
                    name = "") +
  stat_ellipse(aes(color=status_new), linewidth=1.5)+
  scale_color_manual(values = statuscol) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 18, colour = "black"),
        axis.text.y = element_text(size = 18, colour = "black"),
        axis.title.x = element_text(size = 18, colour = "black"),
        axis.title.y = element_text(size = 18, colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        plot.title = element_text(size = 18, hjust = 0.5)) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 16, colour = "black"),
    legend.text = element_text(size = 16, colour = "black")
  )
dev.off()
