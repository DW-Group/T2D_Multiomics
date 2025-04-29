library(tidyverse)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/mgx_t2d/species_t2d_trend_data.RData")

# Fig. species distribution--------
table(spe_info_all$present)
# 1  2  3  4  5  6
# 20 34 23 12 30 76

nrow(spe_info_all) #195
length(spe_info_all$present[spe_info_all$present==1])/nrow(spe_info_all)
# 20; 10.26%
length(spe_info_all$present[spe_info_all$present>1&spe_info_all$present<6])/nrow(spe_info_all)
# 99; 50.77%
length(spe_info_all$present[spe_info_all$present==6])/nrow(spe_info_all)
# 76; 38.97%

spe_all <- merge(spe_info_all, spe_relab, by="species")
universal <- spe_all %>% filter(present==6) %>% arrange(desc(relab))
universal <- universal[1:25, ]
universal$group <- "Universal"
overlap <- spe_all %>% filter(present!=6 & present!=1) %>% arrange(desc(relab))
overlap <- overlap[1:25, ]
overlap$group <- "Overlapping"
solo <- spe_all %>% filter(present==1) %>% arrange(desc(relab))
solo <- solo[1:25,]
solo$group <- "Singular"

spe_plot <- rbind(universal,overlap,solo)
spe_plot2 <- spe_plot %>% dplyr::select(-"relab") %>%
  pivot_longer(cols = c("mbs_relab","mlvs_relab","pedersen_relab","sol_relab",
                        "segal1_relab","direct_relab"), 
               names_to = "study",values_to = "relab") %>%
  filter(relab>0)
spe_plot2$group <- factor(spe_plot2$group,
                          levels=c("Universal", "Overlapping", "Singular"))
spe_plot2$species <- gsub("s__", "", spe_plot2$species)
spe_plot2$species <- gsub("_", " ", spe_plot2$species)

spe_order <- spe_plot2 %>% group_by(group, species) %>%
  summarise(mean=mean(relab)) %>%
  ungroup(species) %>% arrange(group, desc(mean))
spe_plot2$species <- factor(spe_plot2$species, levels=spe_order$species)

spe_plot2 <- spe_plot2 %>%
  mutate(study=case_when(study=="mbs_relab" ~ "NHSII",
                         study=="mlvs_relab" ~ "HPFS",
                         study=="direct_relab" ~ "DIRECT-PLUS",
                         study=="sol_relab" ~ "HCHS/SOL",
                         study=="pedersen_relab" ~ "MetaCardis",
                         study=="segal1_relab" ~ "Talmor-Barkan_2022"
  ))
spe_plot2$study <-
  factor(spe_plot2$study,
         levels=c("DIRECT-PLUS", "HCHS/SOL", "HPFS", 
                  "MetaCardis", "NHSII", "Talmor-Barkan_2022"
         ))

plain <- function(x,...) {
  format(x, ..., scientific = FALSE, drop0trailing = TRUE)
}

library(ggsci)
(col_pal <- pal_nejm(palette = "default")(6))
scales::show_col(col_pal)

studycol <- c(
  `DIRECT-PLUS`= col_pal[6],
  `HCHS/SOL` = col_pal[5],
  NHSII = col_pal[4],
  HPFS = col_pal[3], 
  MetaCardis = col_pal[2],
  `Talmor-Barkan_2022` = col_pal[1]
)

spe_dist <- ggplot(spe_plot2, aes(x=species, y=relab)) +
  geom_boxplot() +
  geom_point(aes(fill=study), shape=21, size=3) +
  scale_fill_manual(values=studycol) +
  ylab("Relative abundance (%)") +
  xlab("") +
  scale_y_log10(breaks=c(10, 1, 0.1, 0.01), labels=plain) +
  facet_wrap( ~ group, scales = "free_x") +
  theme(strip.text.x = element_text(hjust = 0.5, color="black", size=14)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,
                                   size = 11, face = "italic", color="black"),
        axis.text.y = element_text(size = 14, color="black"),
        axis.title.y = element_text(size = 14, color="black")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black")) +
  theme(legend.position = c(.85, .8),
        legend.title = element_blank(),
        legend.text = element_text(size = 10, color="black"),
        legend.background=element_blank(),
        legend.key.height= unit(.3, 'cm'),
        legend.key.width= unit(.3, 'cm'),
        legend.key=element_rect(fill="white")) +
  guides(fill = guide_legend(nrow=5, byrow=TRUE))
ggsave(spe_dist, filename="/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/figS1.pdf",
       width =14, height = 6.4)
