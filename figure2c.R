library(openxlsx)
library(tidyverse)
library(ggplot2)
library(scales)

# Sex hormone
# Caffeine metabolism
# Phenol metabolism
# Secondary bile acids
# AAA metabolism
# Sphingomyelin
# Carnitines
# Food/nutrient biomarker

# dodge barplot ------
## by total ------
load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/variance/mbx_R2_bootstrapping_meta.RData")
meta_r2_spe_food_tmp <- meta_r2_t2
gdata::keep(meta_r2_spe_food_tmp, sure=T)
View(table(meta_r2_spe_food_tmp$SUB_PATHWAY2))

meta_r2_top <- meta_r2_spe_food_tmp %>% 
  filter(SUB_PATHWAY2 %in% c("Caffeine metabolism",
                             "Aromatic amino acids",
                             "Bile acids",
                             "Phenol metabolism",
                             "Branched-chain amino acids",
                             # "Carnitines and acyl carnitines",
                             "Sphingomyelins") | mets_name %in% c("hydroxy-CMPF",
                                                                 "beta-cryptoxanthin")) %>% 
  mutate(group=case_when(SUB_PATHWAY2=="Bile acids" ~ "Secondary bile acids",
                         SUB_PATHWAY2 %in% c("Cofactors and vitamins", "Fatty acids") ~ "Food/nutrient biomarker",
                         TRUE ~ SUB_PATHWAY2),
         mets_name=case_when(mets_name=="5-acetylamino-6-amino-3-methyluracil" ~ "AAMU",
                             TRUE ~ mets_name))

meta_r2_top <- meta_r2_top %>% mutate(group=factor(group,
                                                   levels=c("Phenol metabolism",
                                                            "Secondary bile acids",
                                                            "Aromatic amino acids",
                                                            "Branched-chain amino acids",
                                                            # "Carnitines and acyl carnitines",
                                                            "Sphingomyelins",
                                                            "Caffeine metabolism",
                                                            "Food/nutrient biomarker")),
                                      var_type=factor(var_type,
                                                      levels=c("Species_food","Species","Food")),
                                      R2=R2*100,
                                      LCL=LCL*100,
                                      UCL=UCL*100)
meta_r2_top_tmp <- meta_r2_top %>% group_by(group, mets_name) %>% summarise(max=max(R2)) %>% 
  arrange(group, desc(max))
meta_r2_top$mets_name <- 
  factor(meta_r2_top$mets_name, levels = (meta_r2_top_tmp$mets_name)) # rev

barplot <- ggplot(meta_r2_top,aes(x=mets_name,
                                  y=R2,
                                  fill=var_type))+
  geom_bar(stat = "identity",position = position_dodge(0.8))+
  geom_errorbar(aes(ymin=LCL,ymax=UCL),
                width=0.5,
                linewidth=0.7,
                position = position_dodge(0.8),
                color="grey40")+
  scale_fill_manual(name="Variable type",
                    values = c("#b93f2b","#0a74b2","#e18726"),
                    labels = c("Species + diet","Species","Diet"))+
  scale_y_continuous(breaks=c(0,10,20,30,40), limits = c(-0.5,40))+
  labs(y="Variance explained (%)")+
  theme_classic()+
  theme(
    panel.border = element_blank(), 
    axis.line.y = element_line(size = 1.5,color = "black"),
    axis.line.x = element_line(size = 0.8,color = "black"),
    axis.title.y = element_text(size = 15,color = "black"),
    axis.text.y = element_text(size = 14,color = "black"),
    axis.ticks.y = element_line(colour = "black", size=0.7),
    axis.ticks.length.y = unit(.15, "cm"),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = c(0.9,0.8),
    legend.text = element_text(size = 12,color = "black"),
    legend.title = element_text(size = 12,color = "black"),
    plot.margin = unit(c(0.5, 2, 0, 0.5), "cm")
  )

xx <- ggplot(meta_r2_top,aes(x=mets_name, y=1)) +
  geom_tile(aes(fill=group)) +
  scale_fill_manual(name="Metabolite group",
                    values = c("#749B58","#CC9900",
                               "#AE1F63","#377EB8",
                               "#FFD300","#924822", "#004CFF"),
                    labels = c("Phenol metabolism","Secondary bile acids",
                               "Aromatic amino acids", "Branched-chain amino acids",
                               "Sphingomyelin","Caffeine metabolism", "Food/nutrient biomarker")) +
  xlab("") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.ticks.x = element_line(colour = "black", size=0.7),
        axis.ticks.length.x = unit(.15, "cm"),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0,
                                   size = 14,color = "black"),
        axis.text.y = element_blank(),
        axis.title =element_blank(),
        legend.position = "none",
        plot.margin = unit(c(0, 2, 0.5, 0.5), "cm"))

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig2/figure2a.pdf",
    width = 14, height = 6, onefile = F)
egg::ggarrange(barplot, xx, 
               nrow=2, ncol=1,
               heights = c(1, 0.07))
dev.off()

# legend
xx_legend <- ggplot(meta_r2_top,aes(x=mets_name, y=1)) +
  geom_tile(aes(fill=group)) +
  scale_fill_manual(name="Metabolite group",
                    values = c("#749B58","#CC9900",
                               "#AE1F63","#377EB8",
                               "#FFD300","#924822", "#004CFF"),
                    labels = c("Phenol metabolism","Secondary bile acids",
                               "Aromatic amino acids", "Branched-chain amino acids",
                               "Sphingomyelin","Caffeine metabolism", "Food/nutrient biomarker")) +
  xlab("") +
  theme_bw() +
  theme(axis.text = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 12,color = "black"),
        legend.title = element_text(size = 12,color = "black"),
        legend.key.height= unit(0.4, 'cm'),
        legend.key.width= unit(0.4, 'cm')
  ) +
  guides(fill = guide_legend(title.position = "top", byrow = TRUE))+
  theme(legend.spacing.y = unit(0.1, "cm"))

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig2/figure2a_lgd.pdf",
    width =9, height = 1, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(xx_legend))
dev.off() # Close the file
