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

select <- dplyr::select

length(unique(meta_corr$species)) # 175
length(unique(meta_corr$CHEM_ID)) # 888

# only for T2D-related metabolites
meta_corr2 <- meta_corr %>% 
  filter(CHEM_ID %in% meta_mets_t2d_sig$CHEM_ID) %>% 
  mutate(rho_q=p.adjust(rho_p,method = "BH",n=length(rho_p)),
         pair=paste0(species, "_", CHEM_ID))

sum(meta_corr2$rho_q<0.25) #1849
sum(meta_corr2$rho_q<0.1) #962

length(unique(meta_corr2[meta_corr2$rho_q<0.25, ]$species)) # 158
length(unique(meta_corr2[meta_corr2$rho_q<0.25, ]$CHEM_ID)) # 291

# check direction for spe-met-t2d
spe_met_t2d <- merge(meta_corr2, meta_spe_t2d[, c(1,3,8)], 
                     by.x="species", by.y="feature", all.x = T)
spe_met_t2d <- merge(spe_met_t2d,
                     meta_mets_t2d_sig[, c("CHEM_ID", "beta", "q", "group")],
                     by.x = "CHEM_ID",
                     by.y = "CHEM_ID")
spe_met_t2d <- spe_met_t2d %>% 
  mutate(check=case_when(rho>0 & coef>0 & beta>0 ~ "pos-pos-pos",
                         rho>0 & coef<0 & beta<0 ~ "neg-neg-pos",
                         rho<0 & coef>0 & beta<0 ~ "pos-neg-neg",
                         rho<0 & coef<0 & beta>0 ~ "neg-pos-neg",
                         TRUE ~ "other"))
table(spe_met_t2d$check)
# neg-neg-pos neg-pos-neg       other pos-neg-neg pos-pos-pos 
# 5525        5598       22811        6637        7016 

spe_info_all <- spe_info_all %>% arrange(desc(ave_relab))
spe1 <- meta_spe_t2d[meta_spe_t2d$type %in% 
                       c("T2D_enriched", "T2D_depleted"), ]$feature
spe2 <- spe_info_all[spe_info_all$ave_prev>30,]$species

### Heatmap filter:
# (met-T2D with q<0.1)
# (expected direction (spe_met_t2d$check))
# present >=4 studies
# strong correlation (abs>0.1 & q<0.1) or very significant (q<0.05)
# selected pathway
meta_corr_sig <- spe_met_t2d %>%
  # filter(rho_q<0.1) %>%
  # filter(check != "other") %>% 
  filter(pair_n>3) %>%
  filter((abs(rho)>=0.1 & rho_q<0.1) |
           (abs(rho)>=0.05 & rho_q<0.05)) %>% 
  filter(SUB_PATHWAY2 %in% c("Aromatic amino acids",
                             "Bile acids",
                             "Branched-chain amino acids",
                             # "Caffeine metabolism",
                             # "Carnitines and acyl carnitines",
                             # "Ceramides",
                             # "Cholines and acyl cholines",
                             # "Cofactors and vitamins",
                             "Phenol metabolism"
                             # "Sphingomyelins",
                             # "Sterols and steroids",
                             # "Sulfur-containing amino acids"
                             ))
table(meta_corr_sig[!duplicated(meta_corr_sig$CHEM_ID),]$SUB_PATHWAY2)
check_spe <- meta_corr_sig %>% group_by(species) %>% summarise(n=n())
meta_corr_sig <- meta_corr_sig %>% 
  filter(species %in% c(intersect(check_spe[check_spe$n>3, ]$species, spe2),
                        spe1))
# View(meta_corr2 %>% filter(species %in% c("s__Eubacterium_eligens",
#                                          "s__Eubacterium_siraeum",
#                                          "s__Eubacterium_sp_CAG_180",
#                                          "s__Faecalibacterium_prausnitzii",
#                                          "s__Coprococcus_eutactus") &
#                             mets_name %in% c("3-phenylpropionate",
#                                              "hippurate",
#                                              "4-methylcatechol sulfate",
#                                              "4-ethylphenylsulfate") &
#                              rho_q<0.01))

summary(meta_corr_sig$r)
dim(meta_corr_sig) # 208/20

table(abs(meta_corr_sig$r)>0.3) # 11
length(unique(meta_corr_sig$CHEM_ID)) # 23 mets
length(unique(meta_corr_sig$species)) # 34 spe
# 4 spe associated with T2D
length(intersect(unique(meta_corr_sig$species), 
                 meta_spe_t2d[meta_spe_t2d$type %in% 
                                c("T2D_enriched", "T2D_depleted"), ]$feature)) 

met_corr_final <- meta_corr2 %>% filter((species %in% meta_corr_sig$species) &
                                          (CHEM_ID %in% meta_corr_sig$CHEM_ID)) %>% 
  mutate(sig_sign=case_when(
    rho_q < 0.25 & rho_q >=0.1 ~ "+",
    rho_q < 0.1 & rho_q >=0.05 ~ "#",
    rho_q < 0.05 ~ "*"
  ))

save(meta_corr2, met_corr_final, 
     file="/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig4/figure4a.RData")

dat_htmp <- met_corr_final[,c("species", "rho", "mets_name")]
dat_htmp_wide <- dat_htmp %>% 
  pivot_wider(values_from = rho, names_from = species) %>% 
  column_to_rownames(var = "mets_name")
dim(dat_htmp_wide) # 23/34
dat_htmp_wide[is.na(dat_htmp_wide)] <- 0
dat_htmp_wide[dat_htmp_wide > 0.3] <- 0.3
dat_htmp_wide[dat_htmp_wide < -0.3] <- -0.3
red_blue <- rev(brewer.pal(n=11,name="RdBu"))
pur_green <- rev(brewer.pal(n=11,name="PRGn"))
breaklist <- seq(-0.3,0.3,by=0.001)

## significant sign
sig.mat.t <- met_corr_final[,c("species", "sig_sign", "mets_name")] %>% 
  pivot_wider(values_from = sig_sign, names_from = species) %>% 
  column_to_rownames(var = "mets_name") %>% 
  as.matrix()
sig.mat.t[is.na(sig.mat.t)] <- " "
identical(colnames(dat_htmp_wide),colnames(sig.mat.t)) # TRUE
identical(rownames(dat_htmp_wide),rownames(sig.mat.t)) # TRUE

## mets group
table(rownames(dat_htmp_wide) %in% rownames(meta_mets_t2d_sig))
meta_corr3 <- meta_corr_sig %>% filter(!duplicated(CHEM_ID))
rownames(meta_corr3) <- meta_corr3$mets_name
mets_group <- meta_corr3[rownames(dat_htmp_wide),] %>% 
  select("group","SUB_PATHWAY2") %>% 
  setnames(c("Metabolite status","Sub_path"))

mets_group <- mets_group %>%
  mutate(`Sub pathway`=`Sub_path`) %>% select(-"Sub_path")
identical(rownames(mets_group),rownames(dat_htmp_wide))
table(mets_group$`Metabolite status`)
table(mets_group$`Sub pathway`)

## species group
# species pcoa correlation
# spe_pcoa_cor <- read.xlsx("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/1.preprocess/spe_pcoa_cor.xlsx")
# meta_spe_t2d2 <- meta_spe_t2d %>% 
#   left_join(spe_pcoa_cor,c("feature"="species"))
# summary(meta_spe_t2d2$pco1_cor)
# summary(meta_spe_t2d2$pco2_cor)
# 
# meta_spe_t2d2 <- meta_spe_t2d2 %>%
#   mutate(
#     pco1_cor_grp=case_when(
#       pco1_cor>0.5 ~ "> 0.50",
#       pco1_cor<=0.5&pco1_cor>0.25 ~ "0.25 - 0.50",
#       pco1_cor<=0.25&pco1_cor>0 ~ "0 - 0.25",
#       pco1_cor<= 0&pco1_cor> -0.25 ~ "-0.25 - 0",
#       pco1_cor<= -0.25&pco1_cor>= -0.5 ~ "-0.5 - -0.25",
#       pco1_cor< -0.5 ~ "< -0.50"),
#     pco2_cor_grp=case_when(
#       pco2_cor>0.5 ~ "> 0.50",
#       pco2_cor<=0.5&pco2_cor>0.25 ~ "0.25 - 0.50",
#       pco2_cor<=0.25&pco2_cor>0 ~ "0 - 0.25",
#       pco2_cor<= 0&pco2_cor> -0.25 ~ "-0.25 - 0",
#       pco2_cor<= -0.25&pco2_cor>= -0.5 ~ "-0.5 - -0.25",
#       pco2_cor< -0.5 ~ "< -0.50"))
# table(meta_spe_t2d2$pco1_cor_grp)
# table(meta_spe_t2d2$pco2_cor_grp)
# rownames(meta_spe_t2d2) <- meta_spe_t2d2$feature
# 
# table(colnames(dat_htmp_wide) %in% meta_spe_t2d2$feature)
# species_group <- meta_spe_t2d2[colnames(dat_htmp_wide),] %>% 
#   select("type","pco2_cor_grp","pco1_cor_grp") %>% 
#   setnames(c("Species status","PCo2 correlation","PCo1 correlation"))

table(colnames(dat_htmp_wide) %in% meta_spe_t2d$feature)
species_group <- meta_spe_t2d[colnames(dat_htmp_wide),] %>%
  select("type") %>%
  setnames(c("Species status"))
rownames(species_group) <- gsub("s__","",rownames(species_group))
rownames(species_group) <- gsub("_"," ",rownames(species_group))
colnames(dat_htmp_wide) <- gsub("s__","",colnames(dat_htmp_wide))
colnames(dat_htmp_wide) <- gsub("_"," ",colnames(dat_htmp_wide))
identical(rownames(species_group), colnames(dat_htmp_wide)) # TRUE
table(species_group$`Species status`)

library(ggsci)
(col_pal2 <- pals::glasbey(15))
(col_pal2 <- alpha(col_pal2,0.7))
scales::show_col(col_pal2)
my_colour = list(
  `Species status` = 
    c(`T2D_depleted`="#0a74b2",
      `T2D_enriched` = "#b93f2b",
      `Nonsignificant` = "grey"
    ),
  `Metabolite status` = 
    c(`T2D_depleted`="#0a74b2",
      `T2D_enriched` = "#b93f2b"
    ),
  `PCo1 correlation` = 
    c("> 0.50" = "#833e13",
      "0.25 - 0.50" = "#b6561b",
      "0 - 0.25" = "#fae3c8",
      "-0.25 - 0" = "#e8fdff",
      "-0.5 - -0.25" = "#65b9c4",
      "< -0.50" = "#3b818d"),
  `PCo2 correlation` = 
    c("> 0.50" = "#833e13",
      "0.25 - 0.50" = "#b6561b",
      "0 - 0.25" = "#fae3c8",
      "-0.25 - 0" = "#e8fdff",
      "-0.5 - -0.25" = "#65b9c4",
      "< -0.50" = "#3b818d"),
  `Sub pathway` = 
    c("Aromatic amino acids" = col_pal2[9],
      "Branched-chain amino acids" = col_pal2[10],
      "Bile acids" = col_pal2[4],
      "Phenol metabolism"=col_pal2[15]))

## cluster method
cluster_method <- "ward.D2"

pdf(paste0("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig4/figure4a_",
             cluster_method,".pdf"),width = 16,height = 10)
  heatmap1 <- pheatmap(dat_htmp_wide,
                       clustering_method = cluster_method,
                       cluster_rows = T,
                       cluster_cols = T,
                       # scale = "column",
                       # scale = "row",
                       display_numbers = sig.mat.t,
                       # fontsize_col = 6,
                       # fontsize_row = 15,
                       fontsize_number = 10,
                       number_color = "grey33",
                       treeheight_col = 30,
                       treeheight_row = 30,
                       cellheight = 12,
                       cellwidth = 12,
                       color = colorRampPalette(pur_green)(length(breaklist)), #(100)
                       breaks = breaklist,
                       annotation_col = species_group,
                       annotation_row = mets_group,
                       annotation_names_col = T,
                       annotation_names_row = T,
                       annotation_colors = my_colour,
                       # labels_row = as.character(rownames(dat.q.r.t)),
                       # show_colnames = FALSE,
                       labels_col = as.character(colnames(dat_htmp_wide)),
                       labels_row = as.character(rownames(dat_htmp_wide)),
                       angle_col = "315"
                       # gaps_row = c(5,8),
                       # gaps_col = c(5)
  )
dev.off()
