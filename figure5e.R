# adding MBS into the T2D3 meta-analysis
library(tidyverse)
library(openxlsx)

pp <- c(
  "s__Gemmiger_formicilis_X100001423", # 4-hydroxyhippurate
  "s__Coprococcus_comes_X100001423", # 4-hydroxyhippurate
  "s__Coprococcus_eutactus_X100001423", # 4-hydroxyhippurate
  
  "s__Roseburia_hominis_X100000463", # indolelactate
  "s__Faecalibacterium_prausnitzii_X100000463", # indolelactate
  "s__Coprococcus_comes_X100000463", # indolelactate
  
  "s__Roseburia_hominis_X266", # cholesterol
  "s__Ruminococcus_bicirculans_X266", # cholesterol
  "s__Faecalibacterium_prausnitzii_X266" # cholesterol
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
         pair=case_when(pair=="s__Coprococcus_comes_X100000463" ~ "Coprococcus comes - indolelactate",
                        pair=="s__Faecalibacterium_prausnitzii_X100000463" ~ "Faecalibacterium prausnitzii - indolelactate",
                        pair=="s__Roseburia_hominis_X100000463" ~ "Roseburia hominis - indolelactate",
                        pair=="s__Roseburia_hominis_X266" ~ "Roseburia hominis - cholesterol",
                        pair=="s__Ruminococcus_bicirculans_X266" ~ "Ruminococcus bicirculans - cholesterol",
                        pair=="s__Faecalibacterium_prausnitzii_X266" ~ "Faecalibacterium prausnitzii - cholesterol",
                        pair=="s__Gemmiger_formicilis_X100001423" ~ "Gemmiger formicilis - 4-hydroxyhippurate",
                        pair=="s__Coprococcus_eutactus_X100001423" ~ "Coprococcus eutactus - 4-hydroxyhippurate",
                        pair=="s__Coprococcus_comes_X100001423" ~ "Coprococcus comes - 4-hydroxyhippurate")
  ) %>% arrange(pair, group, study)
table(data_sub$species, data_sub$mets_name)

data_sub1 <- data_sub %>% 
  filter(mets_name %in% c("4-hydroxyhippurate")) %>% 
  mutate(species=factor(species, 
                     levels=c("Coprococcus comes", "Coprococcus eutactus", "Gemmiger formicilis")))
data_sub2 <- data_sub %>% 
  filter(mets_name %in% c("indolelactate")) %>% 
  mutate(species=factor(species, 
                        levels=c("Coprococcus comes", "Faecalibacterium prausnitzii", "Roseburia hominis")))
data_sub3 <- data_sub %>% 
  filter(mets_name %in% c("cholesterol")) %>% 
  mutate(species=factor(species, 
                        levels=c("Ruminococcus bicirculans", "Faecalibacterium prausnitzii", "Roseburia hominis")))

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

p1 <- data_sub1 %>% 
  ggplot(aes(y=group, x=beta, xmin=lci, xmax=uci, 
             color=group, shape=study)) +
  geom_vline(xintercept = 0, linetype="dashed", linewidth = 0.7) +
  geom_point(size=3.5, aes(group=group), position=position_dodge(width=0.9)) +
  geom_errorbar(width=0.35, linewidth=1.2, aes(group=group), position=position_dodge(width=0.9)) +
  scale_color_manual(values = abd_col, name="Metabolite level",
                     guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(name="Study", values=stdy_shape) +
  # scale_x_continuous(limits = c(-1,1.1)) +
  labs(title="4-hydroxyhippurate",
       x="Effect size")+
  facet_wrap(~ species, dir = "v", ncol=1, 
             strip.position = "top") +
  theme_classic() +
  theme(
        axis.ticks.y = element_blank(),
        axis.line = element_blank(), 
        strip.background = element_blank(),
        strip.text = element_text(size=12, hjust = 0, face = "italic"),
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

p2 <- data_sub2 %>% 
  ggplot(aes(y=group, x=beta, xmin=lci, xmax=uci, 
             color=group, shape=study)) +
  geom_vline(xintercept = 0, linetype="dashed", linewidth = 0.7) +
  geom_point(size=3.5, aes(group=group), position=position_dodge(width=0.9)) +
  geom_errorbar(width=0.35, linewidth=1.2, aes(group=group), position=position_dodge(width=0.9)) +
  scale_color_manual(values = abd_col, name="Metabolite level",
                     guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(name="Study", values=stdy_shape) +
  # scale_x_continuous(limits = c(-1,1.1)) +
  labs(title="indolelactate",
       x="Effect size")+
  facet_wrap(~ species, dir = "v", ncol=1, 
             strip.position = "top") +
  theme_classic() +
  theme(
    axis.ticks.y = element_blank(),
    axis.line = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(size=12, hjust = 0, face = "italic"),
    plot.title = element_text(size=12, hjust = 0.5),
    axis.text.x = element_text(size = 13,color = "black"),
    axis.text.y = element_blank(),
    # axis.text.y = element_text(size = 13,color = "black"),
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

p3 <- data_sub3 %>% 
  ggplot(aes(y=group, x=beta, xmin=lci, xmax=uci, 
             color=group, shape=study)) +
  geom_vline(xintercept = 0, linetype="dashed", linewidth = 0.7) +
  geom_point(size=3.5, aes(group=group), position=position_dodge(width=0.9)) +
  geom_errorbar(width=0.35, linewidth=1.2, aes(group=group), position=position_dodge(width=0.9)) +
  scale_color_manual(values = abd_col, name="Metabolite level",
                     guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(name="Study", values=stdy_shape) +
  # scale_x_continuous(limits = c(-1,1.1)) +
  labs(title="cholesterol",
       x="Effect size")+
  facet_wrap(~ species, dir = "v", ncol=1, 
             strip.position = "top") +
  theme_classic() +
  theme(
    axis.ticks.y = element_blank(),
    axis.line = element_blank(), 
    strip.background = element_blank(),
    strip.text = element_text(size=12, hjust = 0, face = "italic"),
    plot.title = element_text(size=12, hjust = 0.5),
    axis.text.x = element_text(size = 13,color = "black"),
    axis.text.y = element_blank(),
    # axis.text.y = element_text(size = 13,color = "black"),
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

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig5/figure5ee.pdf",
    width = 9, height = 4.5, onefile = F) # Open a new pdf file
egg::ggarrange(p1, p2, p3, nrow = 1)
dev.off() # Close the file

library(scales)
data_sub3$group <- factor(data_sub3$group, 
                          levels = c("High", "Low"))
p_legend <- data_sub3 %>% 
  ggplot(aes(y=group, x=beta, xmin=lci, xmax=uci, 
             color=group, shape=study)) +
  geom_vline(xintercept = 0, linetype="dashed", linewidth = 0.7) +
  geom_point(size=3, aes(group=group), position=position_dodge(width=0.9)) +
  geom_errorbar(width=0.3, linewidth=1, aes(group=group), position=position_dodge(width=0.9)) +
  scale_color_manual(values = abd_col, name="Metabolite level",
                     guide = guide_legend(reverse = TRUE)) +
  scale_shape_manual(name="Study", values=stdy_shape) +
  # scale_x_continuous(limits = c(-1,1.1)) +
  labs(title="cholesterol",
       x="Effect size")+
  facet_wrap(~ species, dir = "v", ncol=1, 
             strip.position = "top") +
  theme(legend.position = "right",
        legend.text = element_text(size = 12,color = "black"),
        legend.title = element_text(size = 12,color = "black"),
        legend.key = element_rect(colour = "transparent", fill = "white")) +
  guides(color = guide_legend(title.position = "top"),
         shape = guide_legend(title.position = "top"))

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig5/figure5e_lgd.pdf",
    width =8, height = 5, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(p_legend))
dev.off() # Close the file
