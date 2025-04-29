library(openxlsx)
library(tidyverse)
library(data.table)
library(ggplot2)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/prediction/prediction_rf_lodo.RData")
gdata::keep(dat_comb, sure=T)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/anpan/3.anpan_pglmm_c_mets_results.RData")
# re_pglmm_c, all PGLMM results
# re_pglmm_c2, add met name
# re_pglmm_c3, model == "base_fit" & elpd_diff != 0, the tree model is better than the base model
# re_pglmm_c4, elpd_diff < (-2)

length(unique(re_pglmm_c4$mets)) # 68
length(unique(re_pglmm_c4$species)) # 47

# Species with T2D -----
load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/mgx_t2d/species_t2d_combine.RData")
meta_spe_t2d$feature <- gsub("s__", "", meta_spe_t2d$feature)
rownames(meta_spe_t2d) <- meta_spe_t2d$feature
table(meta_spe_t2d$type)
# Other T2D_depleted_nonsig    T2D_depleted_sig T2D_enriched_nonsig    T2D_enriched_sig 
#    38                  62                   2                  67                  7 

# Species with T2D, by met -----
load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/interaction/dietPC_energy/check_sub_bymet.RData")
sub_t2d3_MBS <- sub_t2d3_MBS %>% mutate(check1=abs(beta_low - beta_high)/abs(beta_low),
                                        check2=abs(beta_low - beta_high)/abs(beta_high))

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/interaction/dietPC_energy/spe_met_pairs_all.RData")

re_pglmm_c4$pair <- paste0(re_pglmm_c4$species, "_", re_pglmm_c4$mets)
for(i in 1:nrow(re_pglmm_c4)) {
  temp <- dat_comb[, c("study", "status_new", re_pglmm_c4$species[i], re_pglmm_c4$mets[i])]
  temp <- temp[complete.cases(temp), ]
  re_pglmm_c4[i, "n_all"] <- nrow(temp)
  re_pglmm_c4[i, "n_con"] <- sum(temp$status_new=="Con")
  re_pglmm_c4[i, "n_pre"] <- sum(temp$status_new=="Pre")
  re_pglmm_c4[i, "n_t2d"] <- sum(temp$status_new=="T2D")
}

re_pglmm_c5 <- re_pglmm_c4 %>% 
  filter(pair %in% sub_t2d3_MBS[sub_t2d3_MBS$check1>1 | sub_t2d3_MBS$check2>1, ]$pair) %>% 
  filter(species %in% spe_info_all[spe_info_all$ave_prev>75, ]$species) %>% 
  filter(!(SUB_PATHWAY2 %in% c("Sphingomyelins", "Ceramides", "Cofactors and vitamins",
                               "Carnitines and acyl carnitines"))) %>% 
  filter(pair %in% mgx_mbx_int_pair[mgx_mbx_int_pair$n>4, ]$pair) %>% 
  filter(elpd_diff < (-4))

re_pglmm_c5$feature <- gsub("s__", "", re_pglmm_c5$species)
re_pglmm_c5 <- merge(re_pglmm_c5, meta_spe_t2d[, c("feature", "coef", "qval.fdr")],
                     by = "feature")
re_pglmm_c5 <- merge(re_pglmm_c5, sub_t2d3_MBS[, -2],
                     by = "pair")
re_pglmm_c5 <- re_pglmm_c5 %>% mutate(
  label=paste0(gsub("_", " ", feature),
               " - ",
               mets_name),
  delta=abs(elpd_diff),
  lci=delta - se_diff,
  uci=delta + se_diff) %>% arrange(delta)
re_pglmm_c5$label <- factor(re_pglmm_c5$label, levels=re_pglmm_c5$label)

### panel 1
p5data1 <- re_pglmm_c5[, c("label", "coef", "beta_low", "beta_high")] %>% 
  tidyr::pivot_longer(!label, names_to = "type", values_to = "coef") %>% 
  mutate(type=case_when(type=="coef" ~ "Overall",
                        type=="beta_low" ~ "Low metabolite",
                        type=="beta_high" ~ "High metabolite"),
         type=factor(type, levels=c("Overall", "Low metabolite", "High metabolite")))
summary(p5data1$coef)
sum(p5data1$coef>0.4) # 2
sum(p5data1$coef<(-0.6)) # 1
p5data1$coef[p5data1$coef>0.4] <- 0.4
p5data1$coef[p5data1$coef<(-0.6)] <- (-0.6)

p5data11 <- re_pglmm_c5[, c("label", "qval.fdr", "q_low", "q_high")] %>% 
  tidyr::pivot_longer(!label, names_to = "type", values_to = "fdr") %>% 
  mutate(type=case_when(type=="qval.fdr" ~ "Overall",
                        type=="q_low" ~ "Low metabolite",
                        type=="q_high" ~ "High metabolite"),
         type=factor(type, levels=c("Overall", "Low metabolite", "High metabolite")),
         star=ifelse(fdr<0.25, "*", ""))
identical(p5data1[, c("label", "type")], p5data11[, c("label", "type")]) # TRUE
p5data1$star <- p5data11$star

p6a1 <-ggplot(p5data1, aes(x=type, y=label)) +
  geom_tile(aes(fill=coef)) +
  scale_fill_gradient2(low="#16518E",
                       mid="white",
                       high="#950F26",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.6,0,0.4),
                       limits=c(-0.6,0.4)) +
  geom_text(aes(label=star), color="black", size=7, nudge_y = -0.2) +
  theme_classic() +
  theme(
    axis.line = element_blank(), 
    plot.title = element_text(size=12, hjust = 0.5),
    axis.text.y = element_text(size = 14,color = "black"),
    axis.text.x = element_text(size = 14,color = "black", angle = 45,
                               vjust = 1, hjust = 1),
    axis.title = element_blank(),
    axis.ticks = element_line(colour = "black", size=0.7),
    axis.ticks.length = unit(.15, "cm"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none",
    legend.text = element_text(size = 12,color = "black"),
    legend.title = element_text(size = 12,color = "black"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))

p6a2 <- ggplot(aes(x=delta, y=label, xmin=lci, xmax=uci), data=re_pglmm_c5) +
  geom_point(size=2) +
  geom_errorbar(width=0.3, linewidth=0.7) +
  xlab(expression(Delta*"ELPD")) +
  theme_classic() +
  theme(
    axis.line = element_blank(), 
    plot.title = element_text(size=12, hjust = 0.5),
    axis.text.x = element_text(size = 14,color = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 14,color = "black"),
    axis.title.y = element_blank(),
    axis.ticks = element_line(colour = "black", size=0.7),
    axis.ticks.length = unit(.15, "cm"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none",
    legend.text = element_text(size = 12,color = "black"),
    legend.title = element_text(size = 12,color = "black"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))

library(tidyr)
data6a3 <- re_pglmm_c5 %>% 
  pivot_longer(
    cols = c("n_con", "n_pre", "n_t2d"), 
    names_to = "status",
    values_to = "value"
  ) %>% select(label, status, value) %>% 
  mutate(status=case_when(status=="n_con" ~ "Control",
                          status=="n_pre" ~ "Prediabetes",
                          status=="n_t2d" ~ "T2D"))
statuscol <- c("#0a74b2","#e18726","#b93f2b")
names(statuscol) <- c("Control", "Prediabetes", "T2D")
data6a3$status <- factor(data6a3$status, levels=rev(c("Control", "Prediabetes", "T2D")))

p6a3 <- ggplot(aes(fill=status, y=label, x=value), data=data6a3) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=statuscol) +
  xlab("Sample size") +
  scale_x_continuous(limits=c(0, 2400),
                     breaks=c(0,500,1000,1500,2000)) +
  theme_classic() +
  theme(
    axis.line = element_blank(), 
    plot.title = element_text(size=12, hjust = 0.5),
    axis.text.x = element_text(size = 14,color = "black"),
    axis.text.y = element_blank(),
    axis.title.x = element_text(size = 14,color = "black"),
    axis.title.y = element_blank(),
    axis.ticks = element_line(colour = "black", size=0.7),
    axis.ticks.length = unit(.15, "cm"),
    panel.border = element_rect(color = "black", fill = NA, size = 1),
    legend.position = "none",
    legend.text = element_text(size = 12,color = "black"),
    legend.title = element_text(size = 12,color = "black"),
    plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig6/figure6a.pdf",
    width = 12.5, height = 4.5, onefile = F) # Open a new pdf file
egg::ggarrange(p6a1, p6a2, p6a3, 
               nrow = 1, widths = c(0.4, 1.3, 1.3))
dev.off() # Close the file

p6a1_lgd <- ggplot(p5data1, aes(x=type, y=label)) +
  geom_tile(aes(fill=coef)) +
  scale_fill_gradient2(low="#16518E",
                       mid="white",
                       high="#950F26",
                       midpoint = 0,
                       guide = "colourbar",
                       breaks=c(-0.6,0,0.4),
                       limits=c(-0.6,0.4)) +
  labs(fill="Species-T2D association") +
  theme(
    legend.text = element_text(size = 12,color = "black"),
    legend.title = element_text(size = 12,color = "black"),
    legend.position = "bottom") +
  guides(fill = guide_colourbar(title.position="top",
                                title.hjust = 0.5,
                                barheight = 1.5,
                                barwidth = 10))
pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig6/figure6a_lgd1.pdf",
    width = 3, height = 3, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(p6a1_lgd))
dev.off() # Close the file

data6a3$status <- factor(data6a3$status, levels=c("Control", "Prediabetes", "T2D"))
p6a3_lgd <- ggplot(aes(fill=status, y=label, x=value), data=data6a3) +
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(values=statuscol) +
  labs(fill="T2D status") +
  theme_classic() +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12,color = "black"),
    legend.title = element_text(size = 12,color = "black")) +
  guides(fill = guide_legend(title.position="top"))
pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig6/figure6a_lgd2.pdf",
    width = 5, height = 3, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(p6a3_lgd))
dev.off() # Close the file
