# s__Eubacterium_hallii_X100001423 (4-hydroxyhippurate)
# low: 0.2549295, 0.07491236
# high: -0.1901971, 0.2225029
# FDR-inter: 0.204302

library(tidyverse)
library(readr)
library(openxlsx)
library(dplyr)
library(vegan)
library(data.table)
library(ape)
library(treeio)
library(ggtree)
library(ggtreeExtra)
library(phyloseq)
library(ggnewscale)
library(tidytree)
library(fastglm)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/1.preprocess/pool_species_20250115_filtered_mmuphin.RData")
gdata::keep(species_all_pool, metadata_spe_all_pool, sure=T)

load("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/anpan/1.metadata_final.RData")
setnames(metadata_final,"mgx_id","sample_id")
identical(colnames(species_all_pool), metadata_spe_all_pool$mgx_id) # TRUE
length(intersect(colnames(species_all_pool), metadata_final$sample_id)) # 2479
spe_final <- species_all_pool[, metadata_final$sample_id]
identical(colnames(spe_final), metadata_final$sample_id) # TRUE
metadata_final <- cbind(metadata_final, t(spe_final))
gdata::keep(metadata_final, sure=T)

select <- dplyr::select

pca = function(mat) {
  centered_mat = scale(mat, scale = FALSE)
  svd_res = svd(centered_mat)
  eigs = svd_res$d^2
  tenth_max_k = sum(eigs > (eigs[1]*.01)) + 1
  message(paste0("k = ", tenth_max_k))
  d = diag(svd_res$d)[,1:tenth_max_k]
  svd_res$u %*% d
}

metName <- "X100001423"
speName <- "s__Eubacterium_hallii"

phylo_effect <- read.xlsx(paste0("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/anpan/3.pglmm/pglmm_phylo_effect_post/phylo_effect_",
                                 metName, "_", speName, ".xlsx"))
phylo_effect <- select(phylo_effect,all_of(c("sample_id","mean")))

## filtered gene table -----
dat_filtered <- fread(paste0("/n/netscratch/dwang_lab/Lab/T2D_multiomics/anpan/genemodel/",
                             speName, "_", metName, "/filter_stats/filtered_", 
                             gsub("s__", "", speName),
                             ".tsv.gz"))

## Metadata ------
metafile <- metadata_final %>% 
  filter(sample_id %in% dat_filtered$sample_id) #%>% 
dim(metafile) 

metafile <- metafile %>% 
  mutate(metf=case_when(metformin==1 ~ "Yes",
                        metformin==0 | (study=="SOL" & is.na(metformin)) ~ "No",
                        TRUE ~ "Missing"))

colSums(is.na(metafile[,c("age", "sex", "bmi", "study","metf","status_new",metName)]))

## gene pre/abs table -----
genefile <- dat_filtered %>% 
  select(-all_of(c("age","sex","bmi","study","status_new",metName))) %>% 
  column_to_rownames(var = "sample_id")
genefile[genefile=="TRUE"] <- 1
genefile[genefile=="FALSE"] <- 0
summary(colSums(genefile))

## Build tree  -------
pca_df <- pca(genefile)
dim(pca_df) 
rownames(pca_df) <- rownames(genefile)

distfile <- vegdist(pca_df, method="euclidean")
# distfile <- vegdist(pca_df, method="jaccard")

treefile <- nj(distfile) |>
  ape::ladderize()
summary(treefile$edge.length)
treefile$edge.length[which(treefile$edge.length<0)]
table(treefile$tip.label %in% metafile$sample_id)

tree_df <- as_tibble(treefile)
tree_df <- tree_df %>% 
  left_join(metafile, by=c("label"="sample_id"))

treefile2 <- as.treedata(tree_df)

# # view internal node
# pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/anpan/5.callout/supp/check.pdf",
#     width = 50, height = 50, onefile = F) # Open a new pdf file
# ggtree(treefile2) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
#   layout_dendrogram() # + geom_tiplab()
# dev.off() # Close the file

# check clade
clade2 <- tree_subset(treefile2, node=1836, levels_back=100)
ggtree(clade2, aes(color=group)) + # geom_tiplab() + 
  xlim(0, 9) + scale_color_manual(values=c("black", "red")) +
  layout_dendrogram()

# group clade based on node number
zz11=groupClade(as_tibble(treefile), 1833) # 1
zz22=groupClade(as_tibble(treefile), 1836) # 1
identical(zz11[, 1:2], zz22[, 1:2])
zz11$group1 <- zz11$group
zz11$group2 <- zz22$group
zz11 <- zz11 %>% mutate(group=case_when(group1==1 ~ 1,
                                        group2==1 ~ 0,
                                        TRUE ~ 9),
                        clade=case_when(group1==1 ~ "Clade 1",
                                        group2==1 ~ "Clade 2",
                                        TRUE ~ "Others"))

zz1 <- zz11[!is.na(zz11$label), ]
metafile <- merge(metafile, zz1[, c("label", "group")],
                  by.x = "sample_id", by.y = "label")

metafile$group <- as.factor(metafile$group)
ggplot(metafile, aes(x=group, y=X100001423, color=group)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) +
  geom_jitter(alpha=0.2, size=0.3, width = 0.1)
ggplot(metafile, aes(x=group, y=s__Eubacterium_hallii, color=group)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) +
  geom_jitter(alpha=0.2, size=0.3, width = 0.1)

table(metafile$study, metafile$group)
# 0   1   2   9
# 0   1   9
# DIRECT-PLUS   6 222  33
# MBS           1  97  59
# MetaCardis  213  43 227
# Segal1        0 207  26
# SOL           1   4 316

# analysis for genes
genefile0 <- genefile
metafile0 <- metafile

metafile <- metafile %>% filter(group %in% c("0", "1")) %>% mutate(group=as.character(group))
genefile <- genefile[metafile$sample_id, ]
identical(rownames(genefile), metafile$sample_id) # TRUE
genefile <- genefile[, colSums(genefile==1)>5 & colSums(genefile==1)<(nrow(metafile)-5)]

library(foreach)
library(doParallel)

# Register a parallel backend
num_cores <- 66 # detectCores() - 1  # Use one less core than available to avoid overloading
cl <- makeCluster(num_cores)
registerDoParallel(cl)

res <- foreach(i = 1:ncol(genefile), .combine = rbind, .packages = 'fastglm') %dopar% {
  y = genefile[, i]
  
  covar <- c("group", "age", "sex", "bmi", "study", "status_new")
  glm_formula = as.formula(paste0(" ~ ", paste(covar, collapse = " + ")))
  x = model.matrix(glm_formula, data = metafile)
  mod <- fastglm(x = x, y = y, family = binomial(), method = 1)
  
  c(colnames(genefile)[i],
    colMeans(genefile == 1)[i],
    summary(mod)[['coefficients']][2, ])
}

# Stop the cluster after computation
stopCluster(cl)

res <- as.data.frame(res)
res[,2:6] <- map_df(res[,2:6],as.numeric)
res$q <- p.adjust(res$`Pr(>|z|)`, method="fdr", n=nrow(res))
colnames(res)[1:2] <- c("gene", "presence")

load("/n/holylfs05/LABS/dwang_lab/Users/zmei/database/uniref90_go_ko.RData")

res$uniref <- sapply(res$gene, function(x) sub(":.*", "", x))
res <- merge(go2uniref90, res, by.x="V1", by.y="uniref", all.y=T)
View(res[res$V1 %in% c("UniRef90_A0A174HHZ4",
                       "UniRef90_A0A174HHV7",
                       "UniRef90_C0EXY0"), ])

rm(go2uniref90, ko2uniref90)
save.image(file="/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/figS7.RData")

library(stringr)
res[res$V1=="UniRef90_C0EXY0", ]$TERM <- "oxidase"
res2 <- res %>% 
  filter(!str_detect(TERM, regex("DNA", ignore_case = F))) %>% 
  filter(!str_detect(TERM, regex("ATP", ignore_case = F))) %>% 
  filter(!str_detect(TERM, regex("ADP", ignore_case = F))) %>% 
  filter(!str_detect(TERM, regex("RNA", ignore_case = F))) %>% 
  filter(!str_detect(TERM, regex("nucleotide", ignore_case = F))) %>% 
  filter(!str_detect(TERM, regex("nucleoside", ignore_case = F))) %>% 
  filter(!str_detect(TERM, regex("TDP", ignore_case = T))) %>% 
  filter(!str_detect(TERM, regex("membrane", ignore_case = F))) %>% 
  filter(!str_detect(TERM, regex("ribosome", ignore_case = F))) %>% 
  filter(!str_detect(TERM, regex("cytoplasm", ignore_case = F))) %>% 
  filter(!str_detect(TERM, regex("NAD", ignore_case = T))) %>% 
  filter(!str_detect(gene, regex("NO_NAME", ignore_case = F))) %>% 
  filter(str_detect(TERM, regex("reductase", ignore_case = F)) | 
           str_detect(TERM, regex("oxidase", ignore_case = F)) |
           # str_detect(TERM, regex("catabolic", ignore_case = F)) |
           # str_detect(TERM, regex("anabolic", ignore_case = F)) |
           str_detect(TERM, regex("catalytic", ignore_case = F)))
res2 <- res2[!duplicated(res2$V1), ] %>% 
  mutate(fdr=p.adjust(`Pr(>|z|)`, method="fdr", n=length(`Pr(>|z|)`)))
sum(res2$fdr<0.1 & res2$Estimate>1)

genefile0$sample_id <- rownames(genefile0)
genes <- genefile0[, c("sample_id",
                       "UniRef90_A0A174HHZ4: 4-hydroxybenzoyl-CoA reductase subunit alpha",
                       "UniRef90_A0A174HHV7: 4-hydroxybenzoyl-CoA reductase subunit gamma",
                       "UniRef90_C0EXY0: Polyphenol oxidase"
)]
colnames(genes)[2:4] <- c("A0A174HHZ4: 4-hydro...unit a", 
                          "A0A174HHV7: 4-hydro...unit g", 
                          "C0EXY0: Polyphenol oxidase")
genes[,2:4] <- map_df(genes[,2:4],as.character)

# Convert from wide to long format
genes <- genes %>%
  pivot_longer(
    cols = c("A0A174HHZ4: 4-hydro...unit a", 
             "A0A174HHV7: 4-hydro...unit g", 
             "C0EXY0: Polyphenol oxidase"),  # Columns to pivot/melt
    names_to = "gene",                    # Name for the new key column
    values_to = "presence"              # Name for the new value column
  ) %>% mutate(gene=factor(gene, levels=c("A0A174HHZ4: 4-hydro...unit a", 
                                          "A0A174HHV7: 4-hydro...unit g", 
                                          "C0EXY0: Polyphenol oxidase")))

tree_df <- as_tibble(treefile)
identical(tree_df[, 1:2], zz11[, 1:2])
tree_df$group <- zz11$clade
tree_df <- tree_df %>% 
  left_join(metafile0[, c("sample_id", "study", 
                          "status_new", metName)], by=c("label"="sample_id")) %>% 
  mutate(study=case_when(study=="Segal1" ~ "Talmor-Barkan_2022",
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
treefile2 <- as.treedata(tree_df)

library(ggsci)
(col_pal <- pal_nejm(palette = "default")(6))
scales::show_col(col_pal)

study_col <- c(
  `DIRECT-PLUS`= col_pal[6],
  `HCHS/SOL` = col_pal[5],
  NHSII = col_pal[4],
  HPFS = col_pal[3], 
  MetaCardis = col_pal[2],
  `Talmor-Barkan_2022` = col_pal[1]
)
gene_col <- c(`1` = "#7fff00", `0` = "#104e8b")
clade_col <- c("Clade 2"="#924822",
               "Clade 1"="#20854E")
status_col <- c("Con"="#0a74b2",
                "Pre"="#e18726",
                "T2D"="#b93f2b")

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/figS7_3.pdf",
    width = 7, height = 10, onefile = F) # Open a new pdf file
ggtree(treefile2, size = 0.5,
       layout = "fan", open.angle = 10) + 
  xlim(0, NA) +
  geom_hilight(mapping=aes(subset = node %in% c(1836, 1833),
                           fill=group),
               # type = "gradient", gradient.direction = 'rt',
               alpha = 0.3) +
  scale_fill_manual(values=clade_col, name="Calde") +
  new_scale_fill() +
  geom_tippoint(shape=16, size=2.5, alpha=1, 
                mapping = aes_string(color=metName), 
                show.legend = T) +
  scale_color_gradient2(low = "blue", 
                        high = "red", 
                        mid = "white", 
                        midpoint = 0,
                        name="Metabolite level") +
  geom_fruit(
    data = phylo_effect,
    geom = geom_tile,
    aes(x="Phylo effect", y=sample_id, fill=mean),
    width = 5,
    axis.params=list(
      axis       = "",
      text.size  = 4,
      text.angle = 90,
      hjust      = 1,
      vjust      = 0.5),
    offset = 0.03, show.legend = T, inherit.aes = F) +
  scale_fill_continuous(name="Phylo effect") +
  new_scale_fill() +
  geom_fruit(
    geom = geom_tile,
    aes(x="Study", y=label, fill=study),
    width = 5,
    axis.params=list(
      axis       = "",
      text.size  = 4,
      text.angle = 90,
      hjust      = 1,
      vjust      = 0.5),
    offset = 0.04, show.legend = F, inherit.aes = F) +
  scale_fill_manual(values = study_col, 
                    name="Study") +
  # new_scale_fill() +
  # geom_fruit(
  #   geom = geom_tile,
  #   aes(x=0.00, y=label, fill=status_new),
  #   width = 5,
  #   offset = 0.08, show.legend = T, inherit.aes = F) +
  # scale_fill_manual(values = status_col, name="") +
  new_scale_fill() +
  geom_fruit(
    data = genes, geom = geom_tile,
    aes(x=gene, y=sample_id, fill=presence),
    width = 6.3,
    axis.params=list(
      axis       = "",
      text.size  = 4,
      text.angle = 90,
      hjust      = 1,
      vjust      = 0.5),
    offset = 0.045, show.legend = T, inherit.aes = F) +
  scale_fill_manual(values = gene_col, 
                    name="Gene presence") +
  layout_rectangular() +
  theme(
    # legend.position="none",
    plot.margin = unit(c(0.1, 0.1, 5, 0.1), "cm"))
dev.off()

######################################################################
data1 <- metafile %>% 
  mutate(clade=case_when(group==0 ~ "Subclade 1",
                         group==1 ~ "Subclade 2"))

# Subclade 2 vs 1, p<0.001
summary(lm(X100001423 ~ as.factor(group), data=data1))  

clade_col <- c("Subclade 1"="#924822",
               "Subclade 2"="#20854E")

pp1 <- ggplot(data1, aes(x=clade, y=X100001423, color=clade)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) +
  geom_jitter(alpha=0.2, size=0.3, width = 0.1) +
  scale_y_continuous(limits = c(-3.3, 3.5), breaks = c(-3,-2,-1,0,1,2,3)) +
  scale_color_manual(values = clade_col) +
  labs(x="",
       y="4-hydroxyhippurate",
       title="E. hallii") +
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

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/figS7_2.pdf",
    width =3, height = 4.5, onefile = F) # Open a new pdf file
pp1
dev.off() # Close the file
