library(openxlsx)
library(stringr)
library(tidyverse)
library(data.table)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/mbx_t2d/mets_t2d_maaslin.RData")
gdata::keep(dat_direct_t, dat_mbs_t, dat_mlvs_t, dat_pedersen_t, dat_segal1_t, dat_sol_t, 
            mets_info_all, sure=T)

### figure 1b
metadata_all <- bind_rows(list(dat_direct_t[, c("study","age","sex","bmi","status_new")], 
                               dat_mbs_t[, c("study","age","sex","bmi","status_new")], 
                               dat_mlvs_t[, c("study","age","sex","bmi","status_new")], 
                               dat_pedersen_t[, c("study","age","sex","bmi","status_new")], 
                               dat_segal1_t[, c("study","age","sex","bmi","status_new")], 
                               dat_sol_t[, c("study","age","sex","bmi","status_new")])) %>% 
  mutate(study=case_when(study=="SOL" ~ "HCHS/SOL",
                         study=="MBS" ~ "NHSII",
                         study=="MLVS" ~ "HPFS",
                         study=="Segal1" ~ "Talmor-Barkan_2022",
                         TRUE ~ study))

library(ggplot2)
library(ggridges)
library(ggsci)

### age & bmi ------
metadata_all_long <- metadata_all %>%
  pivot_longer(cols = c(age, bmi), names_to = "age_bmi", values_to = "value") %>% 
  mutate(study=factor(study, 
                      levels=rev(c("DIRECT-PLUS","HCHS/SOL","HPFS","MetaCardis","NHSII","Talmor-Barkan_2022"))))

p_age_bmi <- ggplot(metadata_all_long, aes(x=value, y=study, fill=study)) +
  geom_density_ridges(alpha=0.7, scale=1.5, color="gray50",size=0.3) +
  scale_fill_nejm() +
  # scale_y_discrete(expand = c(0.02,0.02))+
  labs(x="Value") +
  facet_wrap(~ age_bmi, scales = "free_x",
             labeller = as_labeller(c(age = "Age (y)", bmi = "BMI (kg/m2)"))) +
  theme_classic() +
  theme(
    panel.background = element_blank(),
    axis.text = element_text(size = 11, colour = "black"),
    axis.title = element_blank(),
    # axis.title.x = element_text(size = 11,colour = "black"),
    strip.text = element_text(size=12,colour = "black"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text.x = element_blank()) +
  theme(axis.line.x.bottom=element_line(linewidth=0.5), 
        axis.line.y.left=element_line(linewidth=0.5))

### Sex ------
study_n <- metadata_all %>% 
  group_by(study) %>% 
  summarize(study_n=n())

sex_n <- metadata_all %>% 
  group_by(study,sex) %>% 
  summarize(n=n(),.groups = "keep")
sex_n <- sex_n %>% left_join(study_n,by="study")
sex_n$percent <- sex_n$n*100/sex_n$study_n 
sex_n$study <- factor(sex_n$study,
                      levels=rev(c("DIRECT-PLUS","HCHS/SOL","HPFS","MetaCardis","NHSII","Talmor-Barkan_2022")))

p_sex <- ggplot(sex_n,aes(x=percent,y=study,fill=sex))+
  geom_bar(stat = "identity",position = "stack",alpha = 0.7)+
  scale_fill_manual(breaks = c("Female","Male"),
                    values = c(pal_nejm("default")(6)[1],
                               pal_nejm("default")(6)[6]),
                    name = "Sex"
                    # values = c("#CE7777","#2B3A55")
  ) +
  labs(x="Sex (%)")+
  theme_classic()+
  theme(
    panel.background = element_blank(),
    axis.text = element_text(size = 11,colour = "black"),
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()) +
  theme(axis.line.x.bottom=element_line(linewidth=0.5), 
        axis.line.y.left=element_line(linewidth=0.5))

p_sex_lgd <- p_sex + theme(legend.position = "right")
p_sex_legend <- ggpubr::get_legend(p_sex_lgd)

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/figure1b_sex_legend.pdf",
    width = 1,height = 1)
grid::grid.draw(p_sex_legend)
dev.off()

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/figure1b.pdf", 
    width = 6,height = 4,onefile = F)
egg::ggarrange(p_age_bmi,p_sex, nrow = 1,ncol=2, widths = c(2,1))
dev.off()

### figure 1c
mbs_met01 <- tibble(
  CHEM_ID=colnames(dat_mbs_t)[51:249],
  MBS=1
)
mlvs_met01 <- tibble(
  CHEM_ID=colnames(dat_mlvs_t)[54:257],
  MLVS=1
)
pedersen_met01 <- tibble(
  CHEM_ID=colnames(dat_pedersen_t)[26:888],
  pedersen=1
)
sol_met01 <- tibble(
  CHEM_ID=colnames(dat_sol_t)[53:824],
  sol=1
)
segal1_met01 <- tibble(
  CHEM_ID=colnames(dat_segal1_t)[54:780],
  segal1=1
)
direct_met01 <- tibble(
  CHEM_ID=colnames(dat_direct_t)[46:916],
  direct=1
)

mets01_list <- list(segal1_met01,mbs_met01,pedersen_met01,mlvs_met01,sol_met01,direct_met01)

mets01_t <- mets01_list %>% 
  reduce(full_join,by="CHEM_ID") %>% 
  replace(is.na(.),0)

colnames(mets01_t)[2:7] <- 
  c("Talmor-Barkan_2022","NHSII","MetaCardis","HPFS","HCHS/SOL","DIRECT-PLUS")
mets01_t$n1 <- rowSums(mets01_t[, c("Talmor-Barkan_2022","NHSII","MetaCardis","HPFS","HCHS/SOL","DIRECT-PLUS")])
mets01_t$n2 <- rowSums(mets01_t[, c("Talmor-Barkan_2022","MetaCardis","HCHS/SOL","DIRECT-PLUS")])
View(mets01_t[mets01_t$n1==2 & mets01_t$n2==0, ])
met1090 <- read.xlsx("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/mets_1090.xlsx",
                     sheet="met1090")
View(met1090[met1090$CHEM_ID %in% mets01_t[mets01_t$n1==2 & mets01_t$n2==0, ]$CHEM_ID, ])

# met888 <- read.xlsx("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/mbx_t2d/mets_info_all888.xlsx",
#                     sheet="Sheet5")
# met_extra <- mets_info_all %>% filter(CHEM_ID %in% 
#                                         setdiff(mets01_t$CHEM_ID, met888$CHEM_ID))
# write.xlsx(list("met888"=met888,
#                 "met_extra"=met_extra),
#            file="/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/mets_1090.xlsx")

mets01_t <- mets01_t %>% 
  left_join(met1090[, c("CHEM_ID", "SUPER_PATHWAY2", "SUB_PATHWAY2")], by=c("CHEM_ID")) 

mets01_t$SUPER_PATHWAY2 <-
  factor(mets01_t$SUPER_PATHWAY2,
         levels = c("Lipids", "Amino acids", "Xenobiotics", "Peptides", "Nucleotides",
                    "Carbohydrates and energy metabolism", 
                    "Cofactors and vitamins", "Other metabolites"))

### upset plot -------
library(ComplexUpset)
library(RColorBrewer)
library(ggsci)
library(scales)

mets_pre <- c("Talmor-Barkan_2022","NHSII","MetaCardis","HPFS","HCHS/SOL","DIRECT-PLUS")

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

pdf(file = "/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/figure1c.pdf",
    width = 10, height = 5)
upset(mets01_t, mets_pre, 
      name = "Metabolites",
      set_sizes = (
        upset_set_size(position = "right")
      ),
      base_annotations = list(
        'Intersection size'=intersection_size(
          counts = FALSE,
          mapping=aes(fill=SUPER_PATHWAY2)) + 
          theme(legend.position = "none",
                axis.title = element_text(color="black")) +
          scale_fill_manual(values = col_super_path,
                            # guide="none",
                            name="Super Pathway")),
      min_size=5,
      width_ratio = 0.25
)
dev.off()

pdf(file = "/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/figure1c_legend.pdf",
    width = 10, height = 5)
upset(mets01_t,mets_pre,name = "Metabolites",
      set_sizes = (
        upset_set_size(position = "right")
      ),
      base_annotations = list(
        'Intersection size'=intersection_size(
          counts = FALSE,
          mapping=aes(fill=SUPER_PATHWAY2))+ 
          # theme(legend.position = "none")+
          scale_fill_manual(values = col_super_path,
                            # guide="none",
                            name="Super Pathway")),
      min_size=5,
      width_ratio = 0.2
)
dev.off()

save.image(file="/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/figure1bc.RData")

pdf(file = "/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig1/figure1c_all.pdf",
    width = 15, height = 5)
upset(mets01_t, mets_pre, 
      name = "Metabolites",
      set_sizes = (
        upset_set_size(position = "right")
      ),
      base_annotations = list(
        'Intersection size'=intersection_size(
          counts = T,
          mapping=aes(fill=SUPER_PATHWAY2)) + 
          theme(legend.position = "none",
                axis.title = element_text(color="black")) +
          scale_fill_manual(values = col_super_path,
                            # guide="none",
                            name="Super Pathway")),
      min_size=0,
      width_ratio = 0.25
)
dev.off()
