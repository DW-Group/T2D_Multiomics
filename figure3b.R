library(tidyverse)
library(ggplot2)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/prediction/prediction_rf_lodo.RData")
met_data <- dat_comb[, c("study", "status_new", met1090[!is.na(met1090$order), ]$CHEM_ID)]
met_data2 <- met_data %>% group_by(status_new, study) %>% summarise(across(everything(),~mean(.,na.rm=TRUE)))
met_data2 <- met_data2 %>% mutate(study=case_when(study=="Segal1" ~ "Talmor-Barkan_2022",
                                                  study=="SOL" ~ "HCHS/SOL",
                                                  study=="MBS" ~ "NHSII",
                                                  study=="MLVS" ~ "HPFS",
                                                  TRUE ~ study),
                                  status_new=case_when(status_new=="Con" ~ "Control",
                                                       status_new=="Pre" ~ "Prediabetes",
                                                       status_new=="T2D" ~ "T2D"),
                                  study=factor(study, levels=c("DIRECT-PLUS",
                                                               "HCHS/SOL","HPFS", 
                                                               "MetaCardis","NHSII",
                                                               "Talmor-Barkan_2022")))

library(ggsci)
(col_pal <- pal_nejm(palette = "default",alpha = 0.7)(6))
scales::show_col(col_pal)

study_col <- c(
  `DIRECT-PLUS`= col_pal[6],
  `HCHS/SOL` = col_pal[5],
  NHSII = col_pal[4],
  HPFS = col_pal[3], 
  MetaCardis = col_pal[2],
  `Talmor-Barkan_2022` = col_pal[1]
)

# leucine, X397
p1 <- ggplot(met_data2[!is.na(met_data2$X397), ], aes(x=status_new, y=X397))+
  geom_boxplot(outlier.shape = NA,width=0.5)+
  geom_jitter(aes(color=study),size=3,width = 0.05,height = 0) +
  scale_color_manual(values = study_col)+
  scale_y_continuous(limits = c(-0.45, 0.8), breaks = c(-0.4, 0, 0.4, 0.8)) +
  labs(title = "Leucine",
       y="INT-transformed level")+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(size = 12, color = "black", angle = 30, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size = 12),
        legend.title = element_blank(),
        axis.line.x.bottom=element_line(linewidth=0.5), 
        axis.line.y.left=element_line(linewidth=1)) +
  annotate(geom="text", x=1.5, y=(0.8)-(0.8-(-0.45))/18, size=3.7, col="black",
           label="FDR<0.001")

# methionine sulfoxide, X100000039
p2 <- ggplot(met_data2[!is.na(met_data2$X100000039), ], aes(x=status_new, y=X100000039))+
  geom_boxplot(outlier.shape = NA,width=0.5)+
  geom_jitter(aes(color=study),size=3,width = 0.05,height = 0) +
  scale_color_manual(values = study_col)+
  scale_y_continuous(limits = c(-0.35, 0.8), breaks = c(-0.3, 0, 0.3, 0.6)) +
  labs(title = "Methionine sulfoxide",
       y="INT-transformed level")+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(size = 12, color = "black", angle = 30, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size = 12),
        legend.title = element_blank(),
        axis.line.x.bottom=element_line(linewidth=0.5), 
        axis.line.y.left=element_line(linewidth=1)) +
  annotate(geom="text", x=1.5, y=(0.8)-(0.8-(-0.35))/18, size=3.7, col="black",
           label="FDR=0.02")

# cinnamoylglycine, X100002253
p3 <- ggplot(met_data2[!is.na(met_data2$X100002253), ], aes(x=status_new, y=X100002253))+
  geom_boxplot(outlier.shape = NA,width=0.5)+
  geom_jitter(aes(color=study),size=3,width = 0.05,height = 0) +
  scale_color_manual(values = study_col)+
  scale_y_continuous(limits = c(-0.85, 0.4), breaks = c(-0.8, -0.4, 0, 0.4)) +
  labs(title = "Cinnamoylglycine",
       y="INT-transformed level")+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.text.x = element_blank(),
        # axis.text.x = element_text(size = 12, color = "black", angle = 30, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size = 12),
        legend.title = element_blank(),
        axis.line.x.bottom=element_line(linewidth=0.5), 
        axis.line.y.left=element_line(linewidth=1)) +
  annotate(geom="text", x=1.5, y=(-0.85)+(0.4-(-0.85))/18, size=3.7, col="black",
           label="FDR<0.001")

# linoleoylcholine, X100015760
p4 <- ggplot(met_data2[!is.na(met_data2$X100015760), ], aes(x=status_new, y=X100015760))+
  geom_boxplot(outlier.shape = NA,width=0.5)+
  geom_jitter(aes(color=study),size=3,width = 0.05,height = 0) +
  scale_color_manual(values = study_col)+
  scale_y_continuous(limits = c(-0.85, 0.6), breaks = c(-0.8, -0.4, 0, 0.4)) +
  labs(title = "Linoleoylcholine",
       y="INT-transformed level")+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.y = element_text(size = 12, color = "black"),
        # axis.text.x = element_blank(),
        axis.text.x = element_text(size = 12, color = "black", angle = 30, hjust = 1, vjust = 1),
        axis.title.y = element_text(size = 12, color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5,size = 12),
        legend.title = element_blank(),
        axis.line.x.bottom=element_line(linewidth=0.5), 
        axis.line.y.left=element_line(linewidth=1)) +
  annotate(geom="text", x=1.5, y=(-0.85)+(0.6-(-0.85))/18, size=3.7, col="black",
           label="FDR<0.001")

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig3/figure3b.pdf",
    width = 2.4, height = 9, onefile = F)
egg::ggarrange(p1,p2,p3,p4,nrow = 4)
dev.off()

pp_lgd <- ggplot(met_data2[!is.na(met_data2$X100002253), ], aes(x=status_new, y=X100002253))+
  geom_boxplot(outlier.shape = NA,width=0.5)+
  geom_jitter(aes(color=study),size=3,width = 0.05,height = 0) +
  scale_color_manual(values = study_col,
                     name="Study")+
  theme_classic()+
  theme(axis.text = element_blank(),
        legend.position = "top",
        legend.text = element_text(size = 12,color = "black"),
        legend.title = element_text(size = 12,color = "black"),
        legend.key.height= unit(0.4, 'cm'),
        legend.key.width= unit(0.4, 'cm')) +
  guides(color = guide_legend(title.position = "top", 
                              nrow=3, byrow = TRUE))+
  theme(legend.spacing.y = unit(0.1, "cm"))

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig3/figure3b_lgd.pdf",
    width =6, height = 2, onefile = F) # Open a new pdf file
grid::grid.draw(ggpubr::get_legend(pp_lgd))
dev.off() # Close the file
