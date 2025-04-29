library(openxlsx)
library(tidyverse)

load(file="/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/interaction/dietPC_energy/meta_mbx_spe_int_t2d3_log.RData")
meta_mgx_mbx_int <- merge(met888[, c(2,6,7,8)], meta_mgx_mbx_int, 
                          by="CHEM_ID")

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/mgx_t2d/species_t2d_combine.RData")
spe_info_all <- spe_info_all %>% arrange(desc(ave_relab))

meta_mgx_mbx_int2 <- meta_mgx_mbx_int %>% 
  filter(q_fixed<0.25) # %>% 
  # filter(n>3) %>%
  # filter(species %in% spe_info_all$species[1:100])

meta_mgx_mbx_int3_0 <- meta_mgx_mbx_int %>% 
  group_by(SUPER_PATHWAY2, SUB_PATHWAY2) %>% summarise(n0=n())
meta_mgx_mbx_int3 <- meta_mgx_mbx_int2 %>% 
  group_by(SUPER_PATHWAY2, SUB_PATHWAY2) %>% summarise(n1=n())
identical(meta_mgx_mbx_int3_0[, 1:2], meta_mgx_mbx_int3[, 1:2])
meta_mgx_mbx_int3$por <- meta_mgx_mbx_int3$n1/meta_mgx_mbx_int3_0$n0*100

meta_mgx_mbx_int4 <- meta_mgx_mbx_int3 %>% group_by(SUPER_PATHWAY2) %>% summarise(n2=max(por)) %>% arrange(desc(n2))

meta_mgx_mbx_int3$SUPER_PATHWAY2 <- factor(meta_mgx_mbx_int3$SUPER_PATHWAY2, 
                                           levels=meta_mgx_mbx_int4$SUPER_PATHWAY2)
meta_mgx_mbx_int3 <- meta_mgx_mbx_int3 %>% arrange(SUPER_PATHWAY2, desc(por))
# meta_mgx_mbx_int3$SUB_PATHWAY2 <- factor(meta_mgx_mbx_int3$SUB_PATHWAY2,
#                                          levels=meta_mgx_mbx_int3$SUB_PATHWAY2[c(1:6,8:13,7,14,16:20,15,21:29)])
meta_mgx_mbx_int3$SUB_PATHWAY2 <- factor(meta_mgx_mbx_int3$SUB_PATHWAY2,
                                         levels=meta_mgx_mbx_int3$SUB_PATHWAY2[c(1:2,4,3,5:9,11:17,10,
                                                                                 18,19,21:25,20,26:29)])


library(pals)
library(ggsci)
library(scales)

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
    "Carbohydrates and energy" = col_pal2[8],
    "Other metabolites"="grey70")

library(ggplot2)
bar <- ggplot(aes(x=SUB_PATHWAY2, y=por, fill=SUPER_PATHWAY2), data=meta_mgx_mbx_int3) +
  geom_bar(stat = "identity", width=0.7 ) +
  scale_fill_manual(name="Metabolite super pathway",
                    values = col_super_path)+
  scale_y_continuous(limits = c(0, 22), breaks = c(0,5,10,15,20)) +
  labs(y="Proportion of interactions (%)")+
  theme_classic() +
  theme(
    panel.border = element_blank(), 
    axis.line.y = element_line(size = 0.8,color = "black"),
    axis.line.x = element_line(size = 0.8,color = "black"),
    axis.ticks = element_line(colour = "black", size=0.7),
    axis.ticks.length = unit(.15, "cm"),
    axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0,
                               size = 12,color = "black"),
    axis.text.y = element_text(size = 12,color = "black"),
    axis.title.y = element_text(size = 13,color = "black"),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0.5, 3, 0.5, 0.5), "cm")
  )

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig5/figure5d.pdf",
    width = 11, height = 5.5, onefile = F) # Open a new pdf file
print(bar)
dev.off() # Close the file

# legend
bar_legend <- ggplot(aes(x=SUB_PATHWAY2, y=por, fill=SUPER_PATHWAY2), data=meta_mgx_mbx_int3) +
  geom_bar(stat = "identity", width=0.7 ) +
  scale_fill_manual(name="Super pathway",
                    values = col_super_path)+
  labs(y="Number of interactions")+
  theme_classic() +
  theme(
    panel.border = element_blank(), 
    axis.line.y = element_line(size = 0.8,color = "black"),
    axis.line.x = element_line(size = 0.8,color = "black"),
    axis.ticks = element_line(colour = "black", size=0.7),
    axis.ticks.length = unit(.15, "cm"),
    axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0,
                               size = 12,color = "black"),
    axis.text.y = element_text(size = 12,color = "black"),
    axis.title.y = element_text(size = 13,color = "black"),
    axis.title.x = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 12,color = "black"),
    legend.title = element_text(size = 12,color = "black"),
    legend.key.height= unit(0.5, 'cm'),
    legend.key.width= unit(0.5, 'cm')) +
  guides(fill = guide_legend(title.position = "top"))+
  theme(legend.spacing.y = unit(0.3, "cm"))

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig5/figure5d_legend.pdf",
    width =8, height = 5, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(bar_legend))
dev.off() # Close the file
