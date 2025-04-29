library(tidyverse)
library(ggplot2)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig6/figure6b1.RData")
table(metafile0$group)
ggplot(metafile, aes(x=group, y=X100001423, color=group)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) +
  geom_jitter(alpha=0.2, size=0.3, width = 0.1)

data1 <- metafile0 %>% filter(group!=9) %>% 
  mutate(clade=case_when(group==0 ~ "Subclade 1",
                         group==2 ~ "Subclade 2",
                         group==1 ~ "Subclade 3"))

# Subclade 2 vs 1, p=0.0424
summary(lm(X100001423 ~ as.factor(group), data=data1[data1$group != 1, ]))  

# Subclade 3 vs 1, p=0.0465
summary(lm(X100001423 ~ as.factor(group), data=data1[data1$group != 2, ]))    
  
clade_col <- c("Subclade 1"="#924822",
               "Subclade 2"="#20854E",
               "Subclade 3"="#20854E")

pp1 <- ggplot(data1, aes(x=clade, y=X100001423, color=clade)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) +
  geom_jitter(alpha=0.2, size=0.3, width = 0.1) +
  scale_y_continuous(limits = c(-3.3, 3.3), breaks = c(-3,-2,-1,0,1,2,3)) +
  scale_color_manual(values = clade_col) +
  labs(x="",
       y="4-hydroxyhippurate",
       title="C. comes") +
  theme_classic()+
  theme(axis.text.y = element_text(size = 12,color = "black"),
        axis.text.x = element_text(size = 12,color = "black", 
                                   angle=30, hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size = 12,color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 12,color = "black", face = "italic", hjust = 0.5),
        axis.line = element_blank(),
        panel.border = element_rect(fill=NA, size=1.5),
        legend.position = "none")

gdata::keep(pp1, clade_col, sure=T)  

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig6/figure6b2.RData")
table(metafile0$group)
ggplot(metafile, aes(x=group, y=X266, color=group)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) +
  geom_jitter(alpha=0.2, size=0.3, width = 0.1)

data2 <- metafile0 %>% filter(group!=9) %>% 
  mutate(clade=case_when(group==0 ~ "Subclade 2",
                         group==1 ~ "Subclade 1",
                         group==2 ~ "Subclade 3"))

# Subclade 3 vs 2, p=0.0065
summary(lm(X266 ~ as.factor(group), data=data2[data2$group != 1, ]))  

# Subclade 1 vs 2, p=0.0018
summary(lm(X266 ~ as.factor(group), data=data2[data2$group != 2, ]))    

clade_col <- c("Subclade 2"="#924822",
               "Subclade 1"="#20854E",
               "Subclade 3"="#20854E")

pp2 <- ggplot(data2, aes(x=clade, y=X266, color=clade)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) +
  geom_jitter(alpha=0.2, size=0.3, width = 0.1) +
  scale_y_continuous(limits = c(-3.3, 3.3), breaks = c(-3,-2,-1,0,1,2,3)) +
  scale_color_manual(values = clade_col) +
  labs(x="",
       y="cholesterol",
       title="R. hominis") +
  theme_classic()+
  theme(axis.text.y = element_text(size = 12,color = "black"),
        axis.text.x = element_text(size = 12,color = "black", 
                                   angle=30, hjust=0.5, vjust=0.5),
        axis.title.y = element_text(size = 12,color = "black"),
        axis.title.x = element_blank(),
        plot.title = element_text(size = 12,color = "black", face = "italic", hjust = 0.5),
        axis.line = element_blank(),
        panel.border = element_rect(fill=NA, size=1.5),
        legend.position = "none")

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/figS6.pdf",
    width =6, height = 4.5, onefile = F) # Open a new pdf file
egg::ggarrange(pp1, pp2, ncol = 2)
dev.off() # Close the file
