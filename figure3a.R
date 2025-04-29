library(tidyverse)
library(openxlsx)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/mbx_t2d/mets_t2d_maaslin.RData")
direct <- direct_t2d3_mod4$results %>% filter(value=="status_con3") %>% mutate(study="direct")
mbs <- mbs_t2d3_mod4$results %>% filter(value=="status_con3") %>% mutate(study="mbs")
mlvs <- mlvs_t2d3_mod4$results %>% filter(value=="status_con3") %>% mutate(study="mlvs")
pedersen <- pedersen_t2d3_mod3$results %>% filter(value=="status_con3") %>% mutate(study="pedersen")
segal1 <- segal1_t2d3_mod4$results %>% filter(value=="status_con3") %>% mutate(study="segal1")
sol <- sol_t2d3_mod4$results %>% filter(value=="status_con3") %>% mutate(study="sol")
gdata::keep(direct, mbs, mlvs, pedersen, segal1, sol, sure=T)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/mbx_t2d/mets_t2d_meta_combine.RData")
met888 <- read.xlsx("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/mets_1090.xlsx",
                    sheet="met888")

meta_mets_t2d_sig <- merge(met888[, c("CHEM_ID", "mets_name", "SUPER_PATHWAY2", "SUB_PATHWAY2")],
                           meta_mets_t2d_sig[, -2], by="CHEM_ID", all.y = T)
meta_mets_t2d_sig2 <- meta_mets_t2d_sig %>% filter(n>3 & q<0.15) %>% arrange(desc(abs(beta)))
group1 <- meta_mets_t2d_sig2 %>% group_by(SUPER_PATHWAY2, SUB_PATHWAY2) %>% 
  summarise(n1=n()) %>% 
  filter(n1>5 & !(SUB_PATHWAY2 %in% c("Other xenobiotics", "Other amino acids")))
group2 <- group1 %>% group_by(SUPER_PATHWAY2) %>% 
  summarise(n2=sum(n1)) %>% arrange(desc(n2))
group1 <- group1 %>% 
  mutate(SUPER_PATHWAY2=factor(SUPER_PATHWAY2, levels=group2$SUPER_PATHWAY2)) %>%   
  arrange(SUPER_PATHWAY2, desc(n1))

meta_mets_t2d_sig3 <- meta_mets_t2d_sig2 %>% filter(SUB_PATHWAY2 %in% group1$SUB_PATHWAY2) %>% 
  mutate(SUPER_PATHWAY2=factor(SUPER_PATHWAY2, levels=group2$SUPER_PATHWAY2),
         SUB_PATHWAY2=factor(SUB_PATHWAY2, levels=group1$SUB_PATHWAY2),
         feature=CHEM_ID,
         study="pooled",
         ind=ifelse(beta>0, 1, 0)) %>% 
  arrange(SUPER_PATHWAY2, SUB_PATHWAY2, desc(ind), mets_name) %>% 
  select(feature, mets_name, SUPER_PATHWAY2, SUB_PATHWAY2, beta)
direct <- direct %>% filter(feature %in% meta_mets_t2d_sig3$feature) %>% 
  mutate(direct=coef) %>% select(feature, direct)
mbs <- mbs %>% filter(feature %in% meta_mets_t2d_sig3$feature) %>% 
  mutate(mbs=coef) %>% select(feature, mbs)
mlvs <- mlvs %>% filter(feature %in% meta_mets_t2d_sig3$feature) %>% 
  mutate(mlvs=coef) %>% select(feature, mlvs)
pedersen <- pedersen %>% filter(feature %in% meta_mets_t2d_sig3$feature) %>% 
  mutate(pedersen=coef) %>% select(feature, pedersen)
segal1 <- segal1 %>% filter(feature %in% meta_mets_t2d_sig3$feature) %>% 
  mutate(segal1=coef) %>% select(feature, segal1)
sol <- sol %>% filter(feature %in% meta_mets_t2d_sig3$feature) %>% 
  mutate(sol=coef) %>% select(feature, sol)

ppdata <- Reduce(function(x, y) merge(x, y, by="feature", all = T), 
                 list(meta_mets_t2d_sig3, direct, mbs, mlvs, pedersen, segal1, sol)) %>% 
  mutate(mets_name=factor(mets_name, levels = meta_mets_t2d_sig3$mets_name)) %>% 
  arrange(mets_name)

library(ComplexHeatmap)
library(circlize)
library(gridBase)
library(RColorBrewer)
library(scales)
library(pals)
library(ggsci)

dat.q.r <- 
  select(ppdata,
         c("mets_name",
           "beta",rev(c("direct","sol","mlvs","pedersen","mbs","segal1"))
         )) %>% 
  column_to_rownames(var = "mets_name")

# manually cut the super large/small values (maily in MBS)
sum(dat.q.r>0.4, na.rm = T) # 10
sum(dat.q.r<(-0.4), na.rm = T) # 3
dat.q.r[dat.q.r>0.4] <- 0.4
dat.q.r[dat.q.r< -0.4] <- -0.4

dat.q.r.mat <- as.matrix(dat.q.r)
colnames(dat.q.r.mat) 
colnames(dat.q.r.mat) <-
  c("Meta-analysis",rev(c("DIRECT-PLUS","HCHS/SOL","HPFS","Metacardis","NHSII","Talmor-Barkan_2022")))

split_superpath <- as.character(ppdata$SUPER_PATHWAY2)
(col_pal3=pal_nejm("default",alpha=1)(8))
show_col(col_pal3)
col_super_path <- 
  c(
    "Lipids" = col_pal3[1],
    "Amino acids" = col_pal3[2],
    "Xenobiotics" = col_pal3[3],
    "Peptides" = col_pal3[4],
    "Nucleotides" = col_pal3[5],
    "Carbohydrates and energy" = col_pal3[8])

split_subpath <- as.character(ppdata$SUB_PATHWAY2)
(col_pal2 <- pals::glasbey(15))
(col_pal2 <- alpha(col_pal2,0.7))
scales::show_col(col_pal2)
col_sub_path <- 
  c("Phospholipids" = col_pal2[1],
    "Sphingomyelins" = col_pal2[2],
    "Carnitines and acyl carnitines" = col_pal2[3],
    "Bile acids" = col_pal2[4],
    "Lysophospholipids" = col_pal2[5],
    "Cholines and acyl cholines" = col_pal2[6],
    "Fatty acids" = col_pal2[7],
    "Ceramides" = col_pal2[8],
    "Aromatic amino acids" = col_pal2[9],
    "Branched-chain amino acids" = col_pal2[10],
    "Arginine and proline"=col_pal2[11],
    "Carbohydrates and energy"=col_pal3[8],
    "Peptides"=col_pal2[13],
    "Nucleotides"=col_pal2[14],
    "Phenol metabolism"=col_pal2[15])

my_pal <- rev(brewer.pal(n=11,name="RdBu"))
breaklist <- seq(-0.4,0.4,by=0.001)
col_pal <- colorRampPalette(my_pal)(length(breaklist))
col_fun1 = colorRamp2(breaklist,col_pal)

gaps <- c(rep(1, 14),30)

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig3/figure3a.pdf",width = 10, height = 10)
circos.clear()
circos.par(gap.after = gaps,start.degree=270)
circos.heatmap.initialize(dat.q.r.mat, split = split_subpath,cluster = FALSE)
circos.heatmap(dat.q.r.mat,
               col = col_fun1,
               split = split,
               # show.sector.labels = TRUE,
               cluster = FALSE,
               # dend.side = "inside",
               rownames.side = "outside",
               rownames.cex = 0.7,
               track.height = 0.15,
               bg.border="grey70",
               bg.lwd = 1,
               bg.lty = 1)
circos.track(track.index = get.current.track.index(),
             panel.fun = function(x, y) {
               if(CELL_META$sector.numeric.index == 15) { # position to put colnames (the last sector)
                 cn = rev(colnames(dat.q.r.mat)[1:7])
                 n = length(cn)
                 circos.text(
                   rep(CELL_META$cell.xlim[2], n) + convert_x(2, "mm"),
                   1:n-0.5, cn,cex = 0.6, adj = c(1, 0.5),
                   facing = "outside") #1:n - 0.5
               }}, bg.border = NA)
set_track_gap(mm_h(0.5))
circos.heatmap(split_subpath,
               col = col_sub_path,
               cell.lty = 0,
               track.height=0.02)
set_track_gap(mm_h(1))
circos.heatmap(split_superpath,
               col = col_super_path,
               track.height=0.02)
dev.off()

# # Legends
# pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig3/figure3a_lgd_beta.pdf",
#     width = 1,height = 1.5)
# lgd_beta <- Legend(title = "Beta", col_fun = col_fun1)
# draw(lgd_beta)
# dev.off()
# 
# pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig3/figure3a_lgd_sub.pdf",
#     width = 10, height = 1.5)
# lgd_subpathway <- Legend(title = "Sub pathway", at = names(col_sub_path),
#                          legend_gp = gpar(fill = col_sub_path),
#                          ncol = 4)
# draw(lgd_subpathway)
# dev.off()
# 
# pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig3/figure3a_lgd_super.pdf",
#     width = 5, height = 2)
# lgd_superpathway <- Legend(title = "Super pathway", at = names(col_super_path),
#                            legend_gp = gpar(fill = col_super_path),
#                            ncol = 1)
# draw(lgd_superpathway)
# dev.off()
