# adding MBS into the T2D3 meta-analysis
library(tidyverse)
library(openxlsx)

pp <- c(
  "s__Lactobacillus_rogosae_X342", # glycocholate
  "s__Lactobacillus_rogosae_X1628", # glycochenodeoxycholate
  "s__Lactobacillus_rogosae_X1629"
)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/interaction/dietPC_energy/mbx_spe_int_t2d_log_sub_bymet.RData")
gdata::keep(mbs_mgx_mbx_int_0, mbs_mgx_mbx_int_1, pp, sure=T)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/interaction/dietPC_energy/mbx_spe_int_t2d3_log_sub_bymet.RData")
gdata::keep(direct_mgx_mbx_int_0, direct_mgx_mbx_int_1, 
            mbs_mgx_mbx_int_0, mbs_mgx_mbx_int_1, 
            mlvs_mgx_mbx_int_0, mlvs_mgx_mbx_int_1,
            pedersen_mgx_mbx_int_0, pedersen_mgx_mbx_int_1, 
            segal1_mgx_mbx_int_0, segal1_mgx_mbx_int_1, 
            sol_mgx_mbx_int_0, sol_mgx_mbx_int_1, 
            mets_info_all, pp, sure=T)

direct_mgx_mbx_int_0[,3:5] <- map_df(direct_mgx_mbx_int_0[,3:5],as.numeric)
direct_mgx_mbx_int_0$study <- "DIRECT-PLUS"
mbs_mgx_mbx_int_0[,3:5] <- map_df(mbs_mgx_mbx_int_0[,3:5],as.numeric)
mbs_mgx_mbx_int_0$study <- "NHSII"
mlvs_mgx_mbx_int_0[,3:5] <- map_df(mlvs_mgx_mbx_int_0[,3:5],as.numeric)
mlvs_mgx_mbx_int_0$study <- "HPFS"
pedersen_mgx_mbx_int_0[,3:5] <- map_df(pedersen_mgx_mbx_int_0[,3:5],as.numeric)
pedersen_mgx_mbx_int_0$study <- "MetaCardis"
segal1_mgx_mbx_int_0[,3:5] <- map_df(segal1_mgx_mbx_int_0[,3:5],as.numeric)
segal1_mgx_mbx_int_0$study <- "Talmor-Barkan_2022"
sol_mgx_mbx_int_0[,3:5] <- map_df(sol_mgx_mbx_int_0[,3:5],as.numeric)
sol_mgx_mbx_int_0$study <- "HCHS/SOL"

mgx_mbx_int_0 <- 
  bind_rows(mlvs_mgx_mbx_int_0,pedersen_mgx_mbx_int_0,mbs_mgx_mbx_int_0,
            segal1_mgx_mbx_int_0,sol_mgx_mbx_int_0,direct_mgx_mbx_int_0) 
mgx_mbx_int_0$pair <- paste0(mgx_mbx_int_0$species, "_", mgx_mbx_int_0$CHEM_ID)
pair_n_0 <- mgx_mbx_int_0 %>% group_by(pair, species, CHEM_ID) %>% summarize(n=n())
table(pair_n_0$n)
# 1     2     3     4     5     6 
# 43655 35369 35125 35114 11827  7296

direct_mgx_mbx_int_1[,3:5] <- map_df(direct_mgx_mbx_int_1[,3:5],as.numeric)
direct_mgx_mbx_int_1$study <- "DIRECT-PLUS"
mbs_mgx_mbx_int_1[,3:5] <- map_df(mbs_mgx_mbx_int_1[,3:5],as.numeric)
mbs_mgx_mbx_int_1$study <- "NHSII"
mlvs_mgx_mbx_int_1[,3:5] <- map_df(mlvs_mgx_mbx_int_1[,3:5],as.numeric)
mlvs_mgx_mbx_int_1$study <- "HPFS"
pedersen_mgx_mbx_int_1[,3:5] <- map_df(pedersen_mgx_mbx_int_1[,3:5],as.numeric)
pedersen_mgx_mbx_int_1$study <- "MetaCardis"
segal1_mgx_mbx_int_1[,3:5] <- map_df(segal1_mgx_mbx_int_1[,3:5],as.numeric)
segal1_mgx_mbx_int_1$study <- "Talmor-Barkan_2022"
sol_mgx_mbx_int_1[,3:5] <- map_df(sol_mgx_mbx_int_1[,3:5],as.numeric)
sol_mgx_mbx_int_1$study <- "HCHS/SOL"

mgx_mbx_int_1 <- 
  bind_rows(mlvs_mgx_mbx_int_1,pedersen_mgx_mbx_int_1,mbs_mgx_mbx_int_1,
            segal1_mgx_mbx_int_1,sol_mgx_mbx_int_1,direct_mgx_mbx_int_1) 
mgx_mbx_int_1$pair <- paste0(mgx_mbx_int_1$species, "_", mgx_mbx_int_1$CHEM_ID)
pair_n_1 <- mgx_mbx_int_1 %>% group_by(pair, species, CHEM_ID) %>% summarize(n=n())
table(pair_n_1$n)
# 1     2     3     4     5     6 
# 43655 35369 35125 35114 11827  7296

pair_final <- pair_n_1[pair_n_1$n>1, ]

gdata::keep(mgx_mbx_int_0, mgx_mbx_int_1, pp, sure=T)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/interaction/dietPC_energy/meta_mbx_spe_int_t2d3_log_sub_bymet.RData")
mgx_mbx_int_00 <- mgx_mbx_int_0 %>% filter(pair %in% pp) %>% 
  mutate(group="Low") %>% dplyr::select(pair, study, group, beta, se)
mgx_mbx_int_11 <- mgx_mbx_int_1 %>% filter(pair %in% pp) %>% 
  mutate(group="High") %>% dplyr::select(pair, study, group, beta, se)
meta_mgx_mbx_int_00 <- meta_mgx_mbx_int_0 %>% filter(pair %in% pp) %>% 
  mutate(study="Pooled", group="Low", beta=beta_fixed, se=se_fixed) %>% 
  dplyr::select(pair, study, group, beta, se)
meta_mgx_mbx_int_11 <- meta_mgx_mbx_int_1 %>% filter(pair %in% pp) %>% 
  mutate(study="Pooled", group="High", beta=beta_fixed, se=se_fixed) %>% 
  dplyr::select(pair, study, group, beta, se)

data_sub <- bind_rows(list(mgx_mbx_int_00, meta_mgx_mbx_int_00,
                           mgx_mbx_int_11, meta_mgx_mbx_int_11))

data_sub <- merge(meta_mgx_mbx_int_0[, c("pair", "species", "CHEM_ID")],
                  data_sub, by="pair")
met888 <- read.xlsx("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/mets_1090.xlsx",
                    sheet="met888")
data_sub <- merge(met888[, c("CHEM_ID", "mets_name")],
                  data_sub, by="CHEM_ID")
data_sub <- data_sub %>% 
  mutate(group=factor(group, levels=c("Low", "High")),
         species=gsub("_", " ",
                      gsub("s__", "", species)),
         study=factor(study, levels=c("DIRECT-PLUS", "HCHS/SOL", "HPFS", "MetaCardis", "NHSII", 
                                      "Talmor-Barkan_2022", "Pooled")),
         lci=ifelse(study=="Pooled", beta-1.96*se, NA),
         uci=ifelse(study=="Pooled", beta+1.96*se, NA),
         pair=case_when(pair=="s__Lactobacillus_rogosae_X342" ~ "Lactobacillus rogosae -\nglycocholate\n(FDR<0.001)",
                        pair=="s__Lactobacillus_rogosae_X1628" ~ "Lactobacillus rogosae -\nglycochenodeoxycholate\n(FDR<0.001)",
                        pair=="s__Lactobacillus_rogosae_X1629" ~ "Lactobacillus rogosae -\ntaurochenodeoxycholate\n(FDR<0.001)",
                        pair=="s__Lactobacillus_rogosae_X1648" ~ "Lactobacillus rogosae -\ntaurocholate\n(FDR<0.001)")
  ) %>% arrange(pair, group, study)
table(data_sub$species, data_sub$mets_name)

library(ggplot2)
stdy_shape <- c("DIRECT-PLUS"=0,
                "NHSII"=1,
                "HPFS"=2,
                "MetaCardis"=3,
                "Talmor-Barkan_2022"=4,
                "HCHS/SOL"=5,
                "Pooled"=16)

abd_col <- c("High"="#20854E",
             "Low"="#924822")

abd_col <- c("High"="#AD494AFF",
             "Low"="#5254A3FF")

ppa <- data_sub %>% 
  ggplot(aes(y=group, x=beta, xmin=lci, xmax=uci, 
             color=group, shape=study)) +
  geom_vline(xintercept = 0, linetype="dashed", linewidth = 0.7) +
  geom_point(size=3.5, aes(group=group), position=position_dodge(width=0.9)) +
  geom_errorbar(width=0.35, linewidth=1.2, aes(group=group), position=position_dodge(width=0.9)) +
  scale_color_manual(values = abd_col, name="Metabolite level",
                     guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(name="Study", values=stdy_shape) +
  # scale_x_continuous(limits = c(-1,1.1)) +
  labs(title="",
       x="Effect size")+
  facet_wrap(~ pair, dir = "v", ncol=1, 
             strip.position = "top") +
  theme_classic() +
  theme(
    axis.ticks.y = element_blank(),
    axis.line = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(size=12, hjust = 0),
    plot.title = element_text(size=12, hjust = 0.5),
    axis.text.x = element_text(size = 13,color = "black"),
    axis.text.y = element_blank(),
    # axis.text.y = element_text(size = 12,color = "black"),
    axis.title.x = element_text(size = 12,color = "black"),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(colour = "black", size=0.7),
    axis.ticks.length.x = unit(.15, "cm"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing.y = unit(0.4,"line"),
    legend.position = "none",
    legend.text = element_text(size = 12,color = "black"),
    legend.title = element_text(size = 12,color = "black"),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"))

ppa_lgd <- data_sub %>% 
  ggplot(aes(y=group, x=beta, xmin=lci, xmax=uci, 
             color=group, shape=study)) +
  geom_vline(xintercept = 0, linetype="dashed", linewidth = 0.7) +
  geom_point(size=3.5, aes(group=group), position=position_dodge(width=0.9)) +
  geom_errorbar(width=0.35, linewidth=1.2, aes(group=group), position=position_dodge(width=0.9)) +
  scale_color_manual(values = abd_col, name="Metabolite level",
                     guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(name="Study", values=stdy_shape) +
  # scale_x_continuous(limits = c(-1,1.1)) +
  labs(title="",
       x="Effect size")+
  facet_wrap(~ pair, dir = "v", ncol=1, 
             strip.position = "top") +
  theme(legend.position = "right",
        legend.text = element_text(size = 12,color = "black"),
        legend.title = element_text(size = 12,color = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white")) +
  guides(color = guide_legend(title.position = "top"),
         shape = guide_legend(title.position = "top"))

gdata::keep(ppa, ppa_lgd, sure=T)

#######################################################
pp <- c(
  "s__Clostridium_bolteae_X100003434", 
  "s__Coprococcus_comes_X100003434", 
  "s__Dorea_longicatena_X100003434",
  "s__Faecalibacterium_prausnitzii_X100003434"
)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/interaction/dietPC_energy/mbx_spe_int_t2d_log_sub_bymet.RData")
gdata::keep(mbs_mgx_mbx_int_0, mbs_mgx_mbx_int_1, pp, 
            ppa, ppa_lgd, sure=T)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/interaction/dietPC_energy/mbx_spe_int_t2d3_log_sub_bymet.RData")
gdata::keep(direct_mgx_mbx_int_0, direct_mgx_mbx_int_1, 
            mbs_mgx_mbx_int_0, mbs_mgx_mbx_int_1, 
            mlvs_mgx_mbx_int_0, mlvs_mgx_mbx_int_1,
            pedersen_mgx_mbx_int_0, pedersen_mgx_mbx_int_1, 
            segal1_mgx_mbx_int_0, segal1_mgx_mbx_int_1, 
            sol_mgx_mbx_int_0, sol_mgx_mbx_int_1, 
            mets_info_all, pp, 
            ppa, ppa_lgd, sure=T)

direct_mgx_mbx_int_0[,3:5] <- map_df(direct_mgx_mbx_int_0[,3:5],as.numeric)
direct_mgx_mbx_int_0$study <- "DIRECT-PLUS"
mbs_mgx_mbx_int_0[,3:5] <- map_df(mbs_mgx_mbx_int_0[,3:5],as.numeric)
mbs_mgx_mbx_int_0$study <- "NHSII"
mlvs_mgx_mbx_int_0[,3:5] <- map_df(mlvs_mgx_mbx_int_0[,3:5],as.numeric)
mlvs_mgx_mbx_int_0$study <- "HPFS"
pedersen_mgx_mbx_int_0[,3:5] <- map_df(pedersen_mgx_mbx_int_0[,3:5],as.numeric)
pedersen_mgx_mbx_int_0$study <- "MetaCardis"
segal1_mgx_mbx_int_0[,3:5] <- map_df(segal1_mgx_mbx_int_0[,3:5],as.numeric)
segal1_mgx_mbx_int_0$study <- "Talmor-Barkan_2022"
sol_mgx_mbx_int_0[,3:5] <- map_df(sol_mgx_mbx_int_0[,3:5],as.numeric)
sol_mgx_mbx_int_0$study <- "HCHS/SOL"

mgx_mbx_int_0 <- 
  bind_rows(mlvs_mgx_mbx_int_0,pedersen_mgx_mbx_int_0,mbs_mgx_mbx_int_0,
            segal1_mgx_mbx_int_0,sol_mgx_mbx_int_0,direct_mgx_mbx_int_0) 
mgx_mbx_int_0$pair <- paste0(mgx_mbx_int_0$species, "_", mgx_mbx_int_0$CHEM_ID)
pair_n_0 <- mgx_mbx_int_0 %>% group_by(pair, species, CHEM_ID) %>% summarize(n=n())
table(pair_n_0$n)
# 1     2     3     4     5     6 
# 43655 35369 35125 35114 11827  7296

direct_mgx_mbx_int_1[,3:5] <- map_df(direct_mgx_mbx_int_1[,3:5],as.numeric)
direct_mgx_mbx_int_1$study <- "DIRECT-PLUS"
mbs_mgx_mbx_int_1[,3:5] <- map_df(mbs_mgx_mbx_int_1[,3:5],as.numeric)
mbs_mgx_mbx_int_1$study <- "NHSII"
mlvs_mgx_mbx_int_1[,3:5] <- map_df(mlvs_mgx_mbx_int_1[,3:5],as.numeric)
mlvs_mgx_mbx_int_1$study <- "HPFS"
pedersen_mgx_mbx_int_1[,3:5] <- map_df(pedersen_mgx_mbx_int_1[,3:5],as.numeric)
pedersen_mgx_mbx_int_1$study <- "MetaCardis"
segal1_mgx_mbx_int_1[,3:5] <- map_df(segal1_mgx_mbx_int_1[,3:5],as.numeric)
segal1_mgx_mbx_int_1$study <- "Talmor-Barkan_2022"
sol_mgx_mbx_int_1[,3:5] <- map_df(sol_mgx_mbx_int_1[,3:5],as.numeric)
sol_mgx_mbx_int_1$study <- "HCHS/SOL"

mgx_mbx_int_1 <- 
  bind_rows(mlvs_mgx_mbx_int_1,pedersen_mgx_mbx_int_1,mbs_mgx_mbx_int_1,
            segal1_mgx_mbx_int_1,sol_mgx_mbx_int_1,direct_mgx_mbx_int_1) 
mgx_mbx_int_1$pair <- paste0(mgx_mbx_int_1$species, "_", mgx_mbx_int_1$CHEM_ID)
pair_n_1 <- mgx_mbx_int_1 %>% group_by(pair, species, CHEM_ID) %>% summarize(n=n())
table(pair_n_1$n)
# 1     2     3     4     5     6 
# 43655 35369 35125 35114 11827  7296

pair_final <- pair_n_1[pair_n_1$n>1, ]

gdata::keep(mgx_mbx_int_0, mgx_mbx_int_1, pp, 
            ppa, ppa_lgd, sure=T)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/interaction/dietPC_energy/meta_mbx_spe_int_t2d3_log_sub_bymet.RData")
mgx_mbx_int_00 <- mgx_mbx_int_0 %>% filter(pair %in% pp) %>% 
  mutate(group="Low") %>% dplyr::select(pair, study, group, beta, se)
mgx_mbx_int_11 <- mgx_mbx_int_1 %>% filter(pair %in% pp) %>% 
  mutate(group="High") %>% dplyr::select(pair, study, group, beta, se)
meta_mgx_mbx_int_00 <- meta_mgx_mbx_int_0 %>% filter(pair %in% pp) %>% 
  mutate(study="Pooled", group="Low", beta=beta_fixed, se=se_fixed) %>% 
  dplyr::select(pair, study, group, beta, se)
meta_mgx_mbx_int_11 <- meta_mgx_mbx_int_1 %>% filter(pair %in% pp) %>% 
  mutate(study="Pooled", group="High", beta=beta_fixed, se=se_fixed) %>% 
  dplyr::select(pair, study, group, beta, se)

data_sub <- bind_rows(list(mgx_mbx_int_00, meta_mgx_mbx_int_00,
                           mgx_mbx_int_11, meta_mgx_mbx_int_11))

data_sub <- merge(meta_mgx_mbx_int_0[, c("pair", "species", "CHEM_ID")],
                  data_sub, by="pair")
met888 <- read.xlsx("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/mets_1090.xlsx",
                    sheet="met888")
data_sub <- merge(met888[, c("CHEM_ID", "mets_name")],
                  data_sub, by="CHEM_ID")
data_sub <- data_sub %>% 
  mutate(group=factor(group, levels=c("Low", "High")),
         species=gsub("_", " ",
                      gsub("s__", "", species)),
         study=factor(study, levels=c("DIRECT-PLUS", "HCHS/SOL", "HPFS", "MetaCardis", "NHSII", 
                                      "Talmor-Barkan_2022", "Pooled")),
         lci=ifelse(study=="Pooled", beta-1.96*se, NA),
         uci=ifelse(study=="Pooled", beta+1.96*se, NA),
         pair=case_when(pair=="s__Clostridium_bolteae_X100003434" ~ "Clostridium bolteae -\nimidazole propionate\n(FDR=0.24)",
                        pair=="s__Coprococcus_comes_X100003434" ~ "Coprococcus comes -\nimidazole propionate\n(FDR=0.07)",
                        pair=="s__Dorea_longicatena_X100003434" ~ "Dorea longicatena -\nimidazole propionate\n(FDR=0.03)",
                        pair=="s__Faecalibacterium_prausnitzii_X100003434" ~ "Faecalibacterium prausnitzii -\nimidazole propionate\n(FDR=0.16)")
  ) %>% arrange(pair, group, study)
table(data_sub$species, data_sub$mets_name)

library(ggplot2)
stdy_shape <- c("DIRECT-PLUS"=0,
                "NHSII"=1,
                "HPFS"=2,
                "MetaCardis"=3,
                "Talmor-Barkan_2022"=4,
                "HCHS/SOL"=5,
                "Pooled"=16)

abd_col <- c("High"="#AD494AFF",
             "Low"="#5254A3FF")

data_sub1 <- data_sub[data_sub$species=="Clostridium bolteae", ]
data_sub11 <- data_sub1
data_sub111 <- data_sub1
data_sub11$pair <- "Clostridium bolteae -\nimidazole propionate\n(FDR=0.24) 2"
data_sub111$pair <- "Clostridium bolteae -\nimidazole propionate\n(FDR=0.24) 3"

data_sub2 <- data_sub[data_sub$species!="Clostridium bolteae", ]

data_sub3 <- bind_rows(list(data_sub1, data_sub11, data_sub111))

ppb <- data_sub2 %>% 
  ggplot(aes(y=group, x=beta, xmin=lci, xmax=uci, 
             color=group, shape=study)) +
  geom_vline(xintercept = 0, linetype="dashed", linewidth = 0.7) +
  geom_point(size=3.5, aes(group=group), position=position_dodge(width=0.9)) +
  geom_errorbar(width=0.35, linewidth=1.2, aes(group=group), position=position_dodge(width=0.9)) +
  scale_color_manual(values = abd_col, name="Metabolite level",
                     guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(name="Study", values=stdy_shape) +
  # scale_x_continuous(limits = c(-1,1.1)) +
  labs(title="",
       x="Effect size")+
  facet_wrap(~ pair, dir = "v", ncol=1, 
             strip.position = "top") +
  theme_classic() +
  theme(
    axis.ticks.y = element_blank(),
    axis.line = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(size=12, hjust = 0),
    plot.title = element_text(size=12, hjust = 0.5),
    axis.text.x = element_text(size = 13,color = "black"),
    axis.text.y = element_blank(),
    # axis.text.y = element_text(size = 12,color = "black"),
    axis.title.x = element_text(size = 12,color = "black"),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(colour = "black", size=0.7),
    axis.ticks.length.x = unit(.15, "cm"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing.y = unit(0.4,"line"),
    legend.position = "none",
    legend.text = element_text(size = 12,color = "black"),
    legend.title = element_text(size = 12,color = "black"),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"))

ppc <- data_sub3 %>% 
  ggplot(aes(y=group, x=beta, xmin=lci, xmax=uci, 
             color=group, shape=study)) +
  geom_vline(xintercept = 0, linetype="dashed", linewidth = 0.7) +
  geom_point(size=3.5, aes(group=group), position=position_dodge(width=0.9)) +
  geom_errorbar(width=0.35, linewidth=1.2, aes(group=group), position=position_dodge(width=0.9)) +
  scale_color_manual(values = abd_col, name="Metabolite level",
                     guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(name="Study", values=stdy_shape) +
  # scale_x_continuous(limits = c(-1,1.1)) +
  labs(title="",
       x="Effect size")+
  facet_wrap(~ pair, dir = "v", ncol=1, 
             strip.position = "top") +
  theme_classic() +
  theme(
    axis.ticks.y = element_blank(),
    axis.line = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(size=12, hjust = 0),
    plot.title = element_text(size=12, hjust = 0.5),
    axis.text.x = element_text(size = 13,color = "black"),
    axis.text.y = element_blank(),
    # axis.text.y = element_text(size = 12,color = "black"),
    axis.title.x = element_text(size = 12,color = "black"),
    axis.title.y = element_blank(),
    axis.ticks.x = element_line(colour = "black", size=0.7),
    axis.ticks.length.x = unit(.15, "cm"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    panel.spacing.y = unit(0.4,"line"),
    legend.position = "none",
    legend.text = element_text(size = 12,color = "black"),
    legend.title = element_text(size = 12,color = "black"),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"))

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/figS5.pdf",
    width = 6, height = 6, onefile = F) # Open a new pdf file
egg::ggarrange(ppa, ppc, ncol = 2)
dev.off() # Close the file

# pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/figS5_lgd.pdf",
#     width = 8, height = 5, onefile = F) # Open a new pdf file
# grid::grid.draw(ggpubr::get_legend(ppa_lgd))
# dev.off() # Close the file
