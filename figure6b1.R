# s__Coprococcus_comes_X100001423 (4-hydroxyhippurate)
# low: 0.07183961, 0.6673714
# high: -0.5249144, 0.006942289
# FDR-inter: 0.08627229

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
speName <- "s__Coprococcus_comes"

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

# view internal node
ggtree(treefile2) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) +
  layout_dendrogram() # + geom_tiplab()

# check clade
clade2 <- tree_subset(treefile2, node=1825, levels_back=100)
ggtree(clade2, aes(color=group)) + # geom_tiplab() + 
  xlim(0, 9) + scale_color_manual(values=c("black", "red")) +
  layout_dendrogram()

# group clade based on node number
zz11=groupClade(as_tibble(treefile), 1825) # 1
zz22=groupClade(as_tibble(treefile), 2373) # 1
zz33=groupClade(as_tibble(treefile), 1671) # 1
identical(zz11[, 1:2], zz22[, 1:2])
identical(zz11[, 1:2], zz33[, 1:2])
zz11$group1 <- zz11$group
zz11$group2 <- zz22$group
zz11$group3 <- zz33$group
zz11 <- zz11 %>% mutate(group=case_when(group1==1 ~ 2,
                                        group2==1 ~ 1,
                                        group3==1 ~ 0,
                                        TRUE ~ 9),
                        clade=case_when(group1==1 ~ "Clade 1",
                                        group2==1 ~ "Clade 2",
                                        group3==1 ~ "Clade 3",
                                        TRUE ~ "Others"))

zz1 <- zz11[!is.na(zz11$label), ]
metafile <- merge(metafile, zz1[, c("label", "group")],
                  by.x = "sample_id", by.y = "label")

metafile$group <- as.factor(metafile$group)
ggplot(metafile, aes(x=group, y=X100001423, color=group)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) +
  geom_jitter(alpha=0.2, size=0.3, width = 0.1)
ggplot(metafile, aes(x=group, y=s__Coprococcus_comes, color=group)) + 
  geom_boxplot(outlier.shape = NA, width = 0.3) +
  geom_jitter(alpha=0.2, size=0.3, width = 0.1)

table(metafile$study, metafile$group)
# 0   1   2   9
# DIRECT-PLUS   3 205  27   6
# MBS           2  76   8  48
# MetaCardis  260  16  57 216
# Segal1        6   4 198   2
# SOL          25   0   9 281

# analysis for genes
genefile0 <- genefile
metafile0 <- metafile

metafile1 <- metafile %>% filter(group %in% c("0", "1")) %>% mutate(group=as.character(group))
genefile1 <- genefile[metafile1$sample_id, ]
identical(rownames(genefile1), metafile1$sample_id) # TRUE
genefile1 <- genefile1[, colSums(genefile1==1)>5 & colSums(genefile1==1)<(nrow(metafile1)-5)]

metafile2 <- metafile %>% filter(group %in% c("0", "2")) %>% mutate(group=as.character(group))
genefile2 <- genefile[metafile2$sample_id, ]
identical(rownames(genefile2), metafile2$sample_id) # TRUE
genefile2 <- genefile2[, colSums(genefile2==1)>5 & colSums(genefile2==1)<(nrow(metafile)-5)]

summary(lm(X100001423 ~ group, data=metafile1)) # 0.05
summary(lm(X100001423 ~ group, data=metafile2)) # 0.04

library(foreach)
library(doParallel)

# Register a parallel backend
num_cores <- 66 # detectCores() - 1  # Use one less core than available to avoid overloading
cl <- makeCluster(num_cores)
registerDoParallel(cl)

res1 <- foreach(i = 1:ncol(genefile1), .combine = rbind, .packages = 'fastglm') %dopar% {
  y = genefile1[, i]
  
  covar <- c("group", "age", "sex", "bmi", "study", "status_new")
  glm_formula = as.formula(paste0(" ~ ", paste(covar, collapse = " + ")))
  x = model.matrix(glm_formula, data = metafile1)
  mod <- fastglm(x = x, y = y, family = binomial(), method = 1)
  
  c(colnames(genefile1)[i],
    colMeans(genefile1 == 1)[i],
    summary(mod)[['coefficients']][2, ])
}

res2 <- foreach(i = 1:ncol(genefile2), .combine = rbind, .packages = 'fastglm') %dopar% {
  y = genefile2[, i]
  
  covar <- c("group", "age", "sex", "bmi", "study", "status_new")
  glm_formula = as.formula(paste0(" ~ ", paste(covar, collapse = " + ")))
  x = model.matrix(glm_formula, data = metafile2)
  mod <- fastglm(x = x, y = y, family = binomial(), method = 1)
  
  c(colnames(genefile2)[i],
    colMeans(genefile2 == 1)[i],
    summary(mod)[['coefficients']][2, ])
}

# Stop the cluster after computation
stopCluster(cl)

res1 <- as.data.frame(res1)
res1[,2:6] <- map_df(res1[,2:6],as.numeric)
res1$q <- p.adjust(res1$`Pr(>|z|)`, method="fdr", n=nrow(res1))
colnames(res1)[1:2] <- c("gene", "presence")

res2 <- as.data.frame(res2)
res2[,2:6] <- map_df(res2[,2:6],as.numeric)
res2$q <- p.adjust(res2$`Pr(>|z|)`, method="fdr", n=nrow(res2))
colnames(res2)[1:2] <- c("gene", "presence")

load("/n/holylfs05/LABS/dwang_lab/Users/zmei/database/uniref90_go_ko.RData")

res1$uniref <- sapply(res1$gene, function(x) sub(":.*", "", x))
res1 <- merge(go2uniref90, res1, by.x="V1", by.y="uniref", all.y=T)
View(res1[res1$V1 %in% c("UniRef90_C0BCN8",
                         "UniRef90_A0A173ZLJ0",
                         "UniRef90_A0A173Z1V4",
                         "UniRef90_A0A173YEX7",
                         "UniRef90_A0A173TJB9"), ])

res2$uniref <- sapply(res2$gene, function(x) sub(":.*", "", x))
res2 <- merge(go2uniref90, res2, by.x="V1", by.y="uniref", all.y=T)
View(res2[res2$V1 %in% c("UniRef90_C0BCN8",
                         "UniRef90_A0A173ZLJ0",
                         "UniRef90_A0A173Z1V4",
                         "UniRef90_A0A173YEX7",
                         "UniRef90_A0A173TJB9"), ])
rm(go2uniref90, ko2uniref90)
save.image(file="/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig6/figure6b1.RData")

library(stringr)
res1[res1$V1=="UniRef90_C0BCN8", ]$TERM <- "oxidase"
res2[res2$V1=="UniRef90_C0BCN8", ]$TERM <- "oxidase"
res11 <- res1 %>% 
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
res11 <- res11[!duplicated(res11$V1), ] %>% 
  mutate(fdr=p.adjust(`Pr(>|z|)`, method="fdr", n=length(`Pr(>|z|)`)))
sum(res11$fdr<0.1 & res11$Estimate>1)

res22 <- res2 %>% 
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

res22 <- res22[!duplicated(res22$V1), ] %>% 
  mutate(fdr=p.adjust(`Pr(>|z|)`, method="fdr", n=length(`Pr(>|z|)`)))
sum(res22$fdr<0.1 & res22$Estimate>1)
length(union(res11$V1, res22$V1))

# table S9
colnames(res11)[6:11] <- paste0("res11_", colnames(res11)[6:11])
colnames(res22)[6:11] <- paste0("res22_", colnames(res22)[6:11])
tableS9 <- merge(res11[, c(4,6,7,9,11)],
                 res22[, c(4,6,7,9,11)],
                 by="gene", all = T)
write.xlsx(tableS9, file="/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/tableS9.xlsx")

View(res11[res11$gene %in% c("UniRef90_C0BCN8: Polyphenol oxidase",
                             "UniRef90_A0A173ZLJ0: Benzoyl-CoA reductase subunit B",
                             "UniRef90_A0A173Z1V4: 4-hydroxybenzoyl-CoA reductase subunit alpha",
                             "UniRef90_A0A173YEX7: 4-hydroxybenzoyl-CoA reductase subunit gamma",
                             "UniRef90_A0A173TJB9: Anaerobic benzoate catabolism transcriptional regulator"),])

View(res22[res22$gene %in% c("UniRef90_C0BCN8: Polyphenol oxidase",
                           "UniRef90_A0A173ZLJ0: Benzoyl-CoA reductase subunit B",
                           "UniRef90_A0A173Z1V4: 4-hydroxybenzoyl-CoA reductase subunit alpha",
                           "UniRef90_A0A173YEX7: 4-hydroxybenzoyl-CoA reductase subunit gamma",
                           "UniRef90_A0A173TJB9: Anaerobic benzoate catabolism transcriptional regulator"),])

genefile0$sample_id <- rownames(genefile0)
genes <- genefile0[, c("sample_id",
                       "UniRef90_C0BCN8: Polyphenol oxidase",
                       "UniRef90_A0A173ZLJ0: Benzoyl-CoA reductase subunit B",
                       "UniRef90_A0A173Z1V4: 4-hydroxybenzoyl-CoA reductase subunit alpha",
                       "UniRef90_A0A173YEX7: 4-hydroxybenzoyl-CoA reductase subunit gamma"
)]
colnames(genes)[2:5] <- c("C0BCN8: Polyphenol oxidase", 
                          "A0A173ZLJ0: Benzoyl...unit B", 
                          "A0A173Z1V4: 4-hydro...unit a", 
                          "A0A173YEX7: 4-hydro...unit g")
genes[,2:5] <- map_df(genes[,2:5],as.character)

# Convert from wide to long format
genes <- genes %>%
  pivot_longer(
    cols = c("C0BCN8: Polyphenol oxidase", 
             "A0A173ZLJ0: Benzoyl...unit B", 
             "A0A173Z1V4: 4-hydro...unit a", 
             "A0A173YEX7: 4-hydro...unit g"),  # Columns to pivot/melt
    names_to = "gene",                    # Name for the new key column
    values_to = "presence"              # Name for the new value column
  ) %>% mutate(gene=factor(gene, levels=c("C0BCN8: Polyphenol oxidase", 
                                          "A0A173ZLJ0: Benzoyl...unit B", 
                                          "A0A173Z1V4: 4-hydro...unit a", 
                                          "A0A173YEX7: 4-hydro...unit g")))

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
clade_col <- c("Clade 3"="#924822",
               "Clade 2"="#20854E",
               "Clade 1"="#20854E")
status_col <- c("Con"="#0a74b2",
                "Pre"="#e18726",
                "T2D"="#b93f2b")

pdf("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig6/figure6b1.pdf",
    width = 7, height = 10, onefile = F) # Open a new pdf file
ggtree(treefile2, size = 0.3,
       layout = "fan", open.angle = 10) + 
  xlim(-5, NA) +
  geom_hilight(mapping=aes(subset = node %in% c(1671, 1825, 2373),
                           fill=group),
               # type = "gradient", gradient.direction = 'rt',
               alpha = 0.3) +
  scale_fill_manual(values=clade_col, name="Calde") +
  new_scale_fill() +
  geom_tippoint(shape=16, size=1.5, alpha=1, 
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
    width = 3,
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
    width = 3,
    axis.params=list(
      axis       = "",
      text.size  = 4,
      text.angle = 90,
      hjust      = 1,
      vjust      = 0.5),
    offset = 0.025, show.legend = F, inherit.aes = F) +
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
    width = 4,
    axis.params=list(
      axis       = "",
      text.size  = 4,
      text.angle = 90,
      hjust      = 1,
      vjust      = 0.5),
    offset = 0.03, show.legend = T, inherit.aes = F) +
  scale_fill_manual(values = gene_col, 
                    name="Gene presence") +
  layout_rectangular() +
  theme(
    # legend.position="none",
    plot.margin = unit(c(0.1, 0.1, 5, 0.1), "cm"))
dev.off() 
