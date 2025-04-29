library(openxlsx)
library(stringr)
library(tidyverse)
library(data.table)
library(pheatmap)
library(sjmisc)
library(RColorBrewer)
library(grid)
library(ggpubr)
library(viridis)
library(openxlsx)

# correlation -----
load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/halla/mbx_mgx_for_halla_meta.RData")
met888 <- read.xlsx("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/mets_1090.xlsx",
                    sheet="met888")
setequal(meta_corr$pair, pair_final$pair)
meta_corr$pair_n <- pair_final$n
meta_corr <- merge(meta_corr, met888[, c("CHEM_ID", "mets_name", 
                                         "SUPER_PATHWAY2", "SUB_PATHWAY2")], by="CHEM_ID")
gdata::keep(meta_corr, sure=T)

# Mets with T2D Q<0.25 in either model -----
load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/mbx_t2d/mets_t2d_meta_combine.RData")

# Species with T2D -----
load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/mgx_t2d/species_t2d_combine.RData")
rownames(meta_spe_t2d) <- meta_spe_t2d$feature
table(meta_spe_t2d$type)
# Other T2D_depelted_nonsig    T2D_depleted_sig 
#    38                  62                   2 
# T2D_enriched_nonsig    T2D_enriched_sig 
#                  67                   7 

select <- dplyr::select

length(unique(meta_corr$species)) # 175
length(unique(meta_corr$CHEM_ID)) # 888

# only for T2D-related metabolites
meta_corr2 <- meta_corr %>% 
  filter(CHEM_ID %in% meta_mets_t2d_sig$CHEM_ID) %>% 
  mutate(rho_q=p.adjust(rho_p,method = "BH",n=length(rho_p)),
         pair=paste0(species, "_", CHEM_ID))

sum(meta_corr2$rho_q<0.25) #1849

meta_corr3 <- meta_corr2 %>% filter(rho_q<0.25)
data1 <- data.frame(table(meta_corr3$mets_name))
data1 <- merge(data1, meta_corr3[!duplicated(meta_corr3$mets_name), 
                                 c("mets_name", "SUPER_PATHWAY2", "SUB_PATHWAY2")],
               by.x="Var1", by.y="mets_name") %>% arrange(desc(Freq))

data2 <- data.frame(table(meta_corr3$species))
gdata::keep(data1, data2, sure=T)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/1.preprocess/pool_species_20250115_filtered_mmuphin.RData")
data2 <- merge(data2, taxa_all[,c(2,6,7)],
               by.x="Var1", by.y="species", all.x = T) %>% arrange(desc(Freq))
gdata::keep(data1, data2, sure=T)

table(data1[1:30, ]$SUB_PATHWAY2)
data1 <- data1 %>% mutate(
  Var1=case_when(Var1=="5alpha-androstan-3beta,17alpha-diol disulfate" ~ "5a-androstan-3b,17a-diol disulfate",
                 Var1=="5alpha-pregnan-3beta,20alpha-diol monosulfate (2)" ~ "5a-pregnan-3b,20a-diol monosulfate",
                 TRUE ~ Var1))
data1 <- data1[1:30, ]
data1$Var1 <- factor(data1$Var1, levels=data1$Var1)

col_pal3=ggsci::pal_nejm("default",alpha=1)(10)
scales::show_col(col_pal3)
col_super_path <- 
  c(
    "Lipids" = col_pal3[1],
    "Amino acids" = col_pal3[3],
    "Xenobiotics" = col_pal3[2],
    "Cofactors and vitamins" = col_pal3[6])

# (col_pal2 <- pals::glasbey(15))
# col_sub_path <- 
#   c("Phospholipids" = col_pal2[1],
#     "Chemicals and drugs" = col_pal2[2],
#     "Cofactors and vitamins" = col_pal2[3],
#     "Bile acids" = col_pal2[4],
#     "Other amino acids" = col_pal2[5],
#     "Other xenobiotics" = col_pal2[6],
#     "Fatty acids" = col_pal2[7],
#     "Sex hormones" = col_pal2[8],
#     "Phenol metabolism"=col_pal2[15])

library(ggplot2)
bar1 <- ggplot(aes(x=Var1, y=Freq, fill=SUPER_PATHWAY2), data=data1[1:30, ]) +
  geom_bar(stat = "identity", width=0.7 ) +
  scale_fill_manual(name="Super pathway",
                    values = col_super_path)+
  scale_y_continuous(limits = c(0, 100), breaks = c(0,20,40,60,80,100)) +
  labs(y="Number of significant correlations")+
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
    plot.margin = unit(c(0.5, 5, 0.5, 0.5), "cm")
  )

data2 <- data2[1:30, ] %>% mutate(
  phylum=gsub("p__", "", phylum),
  Var1=gsub("_", " ", gsub("s__", "", Var1))
)
data2$Var1 <- factor(data2$Var1, levels=data2$Var1)

table(data2[1:30, ]$phylum)
col_pal <- ggsci::pal_nejm(palette = "default",alpha = 1)(6)
scales::show_col(col_pal)
col_phylum <- 
  c(
    "Actinobacteria" = col_pal[3],
    "Bacteroidetes" = col_pal[1],
    "Euryarchaeota" = col_pal[5],
    "Firmicutes" = col_pal[2],
    "Proteobacteria" = col_pal[4],
    "Verrucomicrobia" = col_pal[6])

bar2 <- ggplot(aes(x=Var1, y=Freq, fill=phylum), data=data2) +
  geom_bar(stat = "identity", width=0.7) +
  scale_fill_manual(name="Phylum",
                    values = col_phylum)+
  scale_y_continuous(limits = c(0, 65), breaks = c(0,15,30,45,60)) +
  labs(y="Number of significant correlations")+
  theme_classic() +
  theme(
    panel.border = element_blank(), 
    axis.line.y = element_line(size = 0.8,color = "black"),
    axis.line.x = element_line(size = 0.8,color = "black"),
    axis.ticks = element_line(colour = "black", size=0.7),
    axis.ticks.length = unit(.15, "cm"),
    axis.text.x = element_text(angle = -45, vjust = 1, hjust = 0,
                               face = "italic",
                               size = 12,color = "black"),
    axis.text.y = element_text(size = 12,color = "black"),
    axis.title.y = element_text(size = 13,color = "black"),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0.5, 5, 0.5, 0.5), "cm")
  )

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/figS3.pdf",
    width = 11, height = 11, onefile = F) # Open a new pdf file
egg::ggarrange(bar1, bar2, nrow = 2)
dev.off() # Close the file

# legend
bar1_legend <- ggplot(aes(x=Var1, y=Freq, fill=SUPER_PATHWAY2), data=data1[1:30, ]) +
  geom_bar(stat = "identity", width=0.7 ) +
  scale_fill_manual(name="Super pathway",
                    values = col_super_path)+
  scale_y_continuous(limits = c(0, 100), breaks = c(0,20,40,60,80,100)) +
  labs(y="Number of significant correlations")+
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
    legend.key.height= unit(0.4, 'cm'),
    legend.key.width= unit(0.4, 'cm')) +
  guides(fill = guide_legend(title.position = "top",
                             nrow=2, byrow=TRUE))

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/figS3a_lgd.pdf",
    width = 5, height = 3, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(bar1_legend))
dev.off() # Close the file

bar2_legend <- ggplot(aes(x=Var1, y=Freq, fill=phylum), data=data2) +
  geom_bar(stat = "identity", width=0.7) +
  scale_fill_manual(name="Phylum",
                    values = col_phylum)+
  scale_y_continuous(limits = c(0, 65), breaks = c(0,15,30,45,60)) +
  labs(y="Number of significant correlations")+
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
    legend.key.height= unit(0.4, 'cm'),
    legend.key.width= unit(0.4, 'cm')) +
  guides(fill = guide_legend(title.position = "top",
                             nrow=2, byrow=TRUE))

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/figS3b_lgd.pdf",
    width =5, height = 3, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(bar2_legend))
dev.off() # Close the file
