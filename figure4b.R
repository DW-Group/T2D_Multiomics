library(tidyverse)
library(ggplot2)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/interaction/dietPC_energy/mbx_spe_int_t2d_log.RData")
# Log Transformation
LOG <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(log10(y)) # log2: maaslin
}
dat_spe_mbs_log <- dat_spe_mbs_pool[,c(1, 26:ncol(dat_spe_mbs_pool))] 
dat_spe_mbs_log[, 2:ncol(dat_spe_mbs_log)] <-
  map_df(dat_spe_mbs_log[, 2:ncol(dat_spe_mbs_log)], LOG)
dat_mbs_t <- dat_mets_mbs_int %>% 
  inner_join(dat_spe_mbs_log, by="id") %>% 
  inner_join(diet_mbs_pc10, by="id")
gdata::keep(dat_mbs_t, sure=T)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/interaction/dietPC_energy/mbx_spe_int_t2d3_log.RData")
# Log Transformation
LOG <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(log10(y)) # log2: maaslin
}
dat_spe_direct_log <- dat_spe_direct[,c(1,27:ncol(dat_spe_direct))] 
dat_spe_direct_log[,2:ncol(dat_spe_direct_log)] <-
  map_df(dat_spe_direct_log[,2:ncol(dat_spe_direct_log)], LOG)
dat_direct_t <- dat_mets_direct_int %>% 
  inner_join(dat_spe_direct_log, by="id") %>% 
  inner_join(diet_direct_pc10, by="id")

dat_spe_mlvs_log <- dat_spe_mlvs_pool[,c(1, 26:ncol(dat_spe_mlvs_pool))] 
dat_spe_mlvs_log[, 2:ncol(dat_spe_mlvs_log)] <-
  map_df(dat_spe_mlvs_log[, 2:ncol(dat_spe_mlvs_log)], LOG)
dat_mlvs_t <- dat_mets_mlvs_int %>% 
  inner_join(dat_spe_mlvs_log, by="id") %>% 
  inner_join(diet_mlvs_pc10, by="id")

dat_spe_pedersen_log <- dat_spe_pedersen[,c(1,26:ncol(dat_spe_pedersen))] 
dat_spe_pedersen_log[, 2:ncol(dat_spe_pedersen_log)] <-
  map_df(dat_spe_pedersen_log[, 2:ncol(dat_spe_pedersen_log)], LOG)
dat_pedersen_t <- dat_mets_pedersen_int %>% inner_join(dat_spe_pedersen_log, by="id")

dat_spe_segal1_log <- dat_spe_segal1[,c(1,26:ncol(dat_spe_segal1))] 
dat_spe_segal1_log[, 2:ncol(dat_spe_segal1_log)] <-
  map_df(dat_spe_segal1_log[, 2:ncol(dat_spe_segal1_log)], LOG)
dat_segal1_t <- dat_mets_segal1_int %>% 
  inner_join(dat_spe_segal1_log, by="id") %>% 
  inner_join(diet_segal1_pc10, by="id")

dat_spe_sol_log <- dat_spe_sol[,c(1,27:ncol(dat_spe_sol))] 
dat_spe_sol_log[, 2:ncol(dat_spe_sol_log)] <-
  map_df(dat_spe_sol_log[, 2:ncol(dat_spe_sol_log)], LOG)
dat_sol_t <- dat_mets_sol_int %>% 
  inner_join(dat_spe_sol_log, by="id") %>% 
  inner_join(diet_sol_pc10, by="id")

data_all <- bind_rows(list(dat_direct_t, dat_mbs_t, dat_mlvs_t, 
                           dat_pedersen_t, dat_segal1_t, dat_sol_t))
gdata::keep(data_all, sure=T)

# 3-methyl-2-oxobutyrate, s__Bacteroides_vulgatus
# 4-methyl-2-oxopentanoate, s__Bacteroides_vulgatus
# leucine, s__Bacteroides_vulgatus

met <- c("X100000936", "X100000551", "X100001485")
met0 <- c("3-methyl-2-oxobutyrate", "4-methyl-2-oxopentanoate", "leucine")
spe <- "s__Bacteroides_vulgatus"
spe0 <- "B. vulgatus"
plot(data_all$s__Bacteroides_vulgatus, data_all$X100000936)
plot(data_all$s__Bacteroides_vulgatus, data_all$X100000551)
plot(data_all$s__Bacteroides_vulgatus, data_all$X100001485)

pp <- list()
i=1
tempdata <- data.frame(x=data_all[, spe],
                       y=data_all[, met[i]]) %>% filter(!is.na(x) & !is.na(y))

pp[[i]] <- ggplot(tempdata, aes(x=x, y=y)) +
  geom_point(color="#5254A333") +
  geom_smooth(method = "lm", se = T, color="#AD494AFF") +
  scale_x_continuous(limits = c(-4,2), breaks = c(-4,-2,0,2)) +
  labs(y=met0[i], title=spe0, x="Log10(species)") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 13, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        # axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        # axis.title.x = element_text(size = 12, colour = "black"),
        legend.position = "none",
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(linewidth = 0.7),
        # axis.ticks = element_blank(),
        plot.title = element_text(size = 12,hjust = 0.5, 
                                  face = "italic", colour = "black"),
        axis.line.x.bottom=element_line(linewidth=0.5), 
        axis.line.y.left=element_line(linewidth=1))

i=2
tempdata <- data.frame(x=data_all[, spe],
                       y=data_all[, met[i]]) %>% filter(!is.na(x) & !is.na(y))

pp[[i]] <- ggplot(tempdata, aes(x=x, y=y)) +
  geom_point(color="#5254A333") +
  geom_smooth(method = "lm", se = T, color="#AD494AFF") +
  scale_x_continuous(limits = c(-4,2), breaks = c(-4,-2,0,2)) +
  labs(y=met0[i], title=spe0, x="Log10(species)") +
  theme_classic() +
  theme(axis.text.y = element_text(size = 13, colour = "black"),
        axis.title.y = element_text(size = 12, colour = "black"),
        # axis.text.x = element_blank(),
        # axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12, colour = "black"),
        axis.title.x = element_text(size = 12, colour = "black"),
        legend.position = "none",
        axis.ticks.length = unit(0.15, "cm"),
        axis.ticks = element_line(linewidth = 0.7),
        # axis.ticks = element_blank(),
        plot.title = element_blank(),
        axis.line.x.bottom=element_line(linewidth=0.5), 
        axis.line.y.left=element_line(linewidth=1))

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig4/figure4b.pdf",
    width=2.7, height=6.5, onefile = F)
egg::ggarrange(pp[[1]], pp[[2]], nrow=2)
dev.off()

# pp <- list()
# for (i in 1:3) {
#   tempdata <- data.frame(x=data_all[, spe],
#                          y=data_all[, met[i]]) %>% filter(!is.na(x) & !is.na(y))
#   
#   pp[[i]] <- ggplot(tempdata, aes(x=x, y=y)) +
#     geom_point(color="#5254A333") +
#     geom_smooth(method = "lm", se = T, color="#AD494AFF") +
#     labs(y=met0[i], title=spe0, x="Log10(species)") +
#     theme_classic() +
#     theme(axis.text.x = element_text(size = 12, colour = "black"),
#           axis.text.y = element_text(size = 12, colour = "black"),
#           plot.title = element_text(size = 12, hjust = 0.5, 
#                                     face = "italic", colour = "black"),
#           # axis.title.x = element_blank(),
#           axis.title = element_text(size = 12, colour = "black"),
#           axis.ticks.length = unit(0.15, "cm"),
#           axis.line.x.bottom=element_line(linewidth=0.5), 
#           axis.line.y.left=element_line(linewidth=1))
#   
# }
# 
# pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig4/figure4b.pdf",
#     width=2.5, height=8, onefile = F)
# egg::ggarrange(pp[[1]], pp[[2]], pp[[3]], nrow=3)
# dev.off()