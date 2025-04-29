library(tidyverse)
library(data.table)

load(file="/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/2.analysis/permanova/mbx_permanova_pc.RData")

# Combine results --------
### diet and species -----
re_perm_diet_spe <- tibble(
  study=c("MBS","MLVS","HCHS/SOL","Talmor_2022","DIRECT-PLUS"),
  R2_diet_spe=c(re_perm_diet_spe_mbs$R2[1],
                re_perm_diet_spe_mlvs$R2[1],
                re_perm_diet_spe_sol$R2[1],
                re_perm_diet_spe_segal1$R2[1],
                re_perm_diet_spe_direct$R2[1]
  ),
  P_diet_spe=c(re_perm_diet_spe_mbs$`Pr(>F)`[1],
               re_perm_diet_spe_mlvs$`Pr(>F)`[1],
               re_perm_diet_spe_sol$`Pr(>F)`[1],
               re_perm_diet_spe_segal1$`Pr(>F)`[1],
               re_perm_diet_spe_direct$`Pr(>F)`[1]
  )
)

### species -----
re_perm_spe <- tibble(
  study=c("MBS","MLVS","MetaCardis","HCHS/SOL","Talmor_2022","DIRECT-PLUS"),
  R2_spe=c(re_perm_spe_mbs$R2[1],
           re_perm_spe_mlvs$R2[1],
           re_perm_spe_ped$R2[1],
           re_perm_spe_sol$R2[1],
           re_perm_spe_segal1$R2[1],
           re_perm_spe_direct$R2[1]
  ),
  P_spe=c(re_perm_spe_mbs$`Pr(>F)`[1],
          re_perm_spe_mlvs$`Pr(>F)`[1],
          re_perm_spe_ped$`Pr(>F)`[1],
          re_perm_spe_sol$`Pr(>F)`[1],
          re_perm_spe_segal1$`Pr(>F)`[1],
          re_perm_spe_direct$`Pr(>F)`[1]
  )
)

### food -----
re_perm_food <- tibble(
  study=c("MBS","MLVS","HCHS/SOL","Talmor_2022","DIRECT-PLUS"),
  R2_food=c(re_perm_food_mbs$R2[1],
            re_perm_food_mlvs$R2[1],
            re_perm_food_sol$R2[1],
            re_perm_food_segal1$R2[1],
            re_perm_food_direct$R2[1]
  ),
  P_food=c(re_perm_food_mbs$`Pr(>F)`[1],
           re_perm_food_mlvs$`Pr(>F)`[1],
           re_perm_food_sol$`Pr(>F)`[1],
           re_perm_food_segal1$`Pr(>F)`[1],
           re_perm_food_direct$`Pr(>F)`[1])
)

### age -----
re_perm_age <- tibble(
  study=c("MBS","MLVS","MetaCardis","HCHS/SOL","Talmor_2022","DIRECT-PLUS"),
  R2_age=c(re_perm_age_mbs$R2[1],
           re_perm_age_mlvs$R2[1],
           re_perm_age_ped$R2[1],
           re_perm_age_sol$R2[1],
           re_perm_age_segal1$R2[1],
           re_perm_age_direct$R2[1]
  ),
  P_age=c(re_perm_age_mbs$`Pr(>F)`[1],
          re_perm_age_mlvs$`Pr(>F)`[1],
          re_perm_age_ped$`Pr(>F)`[1],
          re_perm_age_sol$`Pr(>F)`[1],
          re_perm_age_segal1$`Pr(>F)`[1],
          re_perm_age_direct$`Pr(>F)`[1])
)

### sex -----
re_perm_sex <- tibble(
  study=c(#"MBS","MLVS",
    "MetaCardis","HCHS/SOL","Talmor_2022","DIRECT-PLUS"),
  R2_sex=c(#re_perm_sex_mbs$R2[1],
    #re_perm_sex_mlvs$R2[1],
    re_perm_sex_ped$R2[1],
    re_perm_sex_sol$R2[1],
    re_perm_sex_segal1$R2[1],
    re_perm_sex_direct$R2[1]
  ),
  P_sex=c(#re_perm_sex_mbs$`Pr(>F)`[1],
    #re_perm_sex_mlvs$`Pr(>F)`[1],
    re_perm_sex_ped$`Pr(>F)`[1],
    re_perm_sex_sol$`Pr(>F)`[1],
    re_perm_sex_segal1$`Pr(>F)`[1],
    re_perm_sex_direct$`Pr(>F)`[1])
)

### bmi -----
re_perm_bmi <- tibble(
  study=c("MBS","MLVS","MetaCardis","HCHS/SOL","Talmor_2022","DIRECT-PLUS"),
  R2_bmi=c(re_perm_bmi_mbs$R2[1],
           re_perm_bmi_mlvs$R2[1],
           re_perm_bmi_ped$R2[1],
           re_perm_bmi_sol$R2[1],
           re_perm_bmi_segal1$R2[1],
           re_perm_bmi_direct$R2[1]
  ),
  P_bmi=c(re_perm_bmi_mbs$`Pr(>F)`[1],
          re_perm_bmi_mlvs$`Pr(>F)`[1],
          re_perm_bmi_ped$`Pr(>F)`[1],
          re_perm_bmi_sol$`Pr(>F)`[1],
          re_perm_bmi_segal1$`Pr(>F)`[1],
          re_perm_bmi_direct$`Pr(>F)`[1])
)

### metformin -----
re_perm_metf <- tibble(
  study=c("MBS","MLVS","MetaCardis","HCHS/SOL","DIRECT-PLUS"), #,"Segal1"
  R2_metf=c(re_perm_metf_mbs$R2[1],
            re_perm_metf_mlvs$R2[1],
            re_perm_metf_ped$R2[1],
            re_perm_metf_sol$R2[1],
            re_perm_metf_direct$R2[1]
            # re_perm_metf_segal1$R2[1]
  ),
  P_metf=c(re_perm_metf_mbs$`Pr(>F)`[1],
           re_perm_metf_mlvs$`Pr(>F)`[1],
           re_perm_metf_ped$`Pr(>F)`[1],
           re_perm_metf_sol$`Pr(>F)`[1],
           re_perm_metf_direct$`Pr(>F)`[1]#,
           # re_perm_metf_segal1$`Pr(>F)`[1]
  ))

### T2D -----
re_perm_t2d <- tibble(
  study=c("MBS","MLVS","MetaCardis","HCHS/SOL","Talmor_2022","DIRECT-PLUS"),
  R2_t2d=c(re_perm_t2d_mbs$R2[1],
           re_perm_t2d_mlvs$R2[1],
           re_perm_t2d_ped$R2[1],
           re_perm_t2d_sol$R2[1],
           re_perm_t2d_segal1$R2[1],
           re_perm_t2d_direct$R2[1]
  ),
  P_t2d=c(re_perm_t2d_mbs$`Pr(>F)`[1],
          re_perm_t2d_mlvs$`Pr(>F)`[1],
          re_perm_t2d_ped$`Pr(>F)`[1],
          re_perm_t2d_sol$`Pr(>F)`[1],
          re_perm_t2d_segal1$`Pr(>F)`[1],
          re_perm_t2d_direct$`Pr(>F)`[1])
)

### clinic -----
re_perm_clinic <- tibble(
  study=c("MBS","MLVS","MetaCardis","HCHS/SOL","Talmor_2022","DIRECT-PLUS"
  ),
  R2_clinic=c(re_perm_clinic_mbs$R2[1],
              re_perm_clinic_mlvs$R2[1],
              re_perm_clinic_ped$R2[1],
              re_perm_clinic_sol$R2[1],
              re_perm_clinic_segal1$R2[1],
              re_perm_clinic_direct$R2[1]
              # re_perm_clinic_segal2$R2[1]
  ),
  P_clinic=c(re_perm_clinic_mbs$`Pr(>F)`[1],
             re_perm_clinic_mlvs$`Pr(>F)`[1],
             re_perm_clinic_ped$`Pr(>F)`[1],
             re_perm_clinic_sol$`Pr(>F)`[1],
             re_perm_clinic_segal1$`Pr(>F)`[1],
             re_perm_clinic_direct$`Pr(>F)`[1]
             # re_perm_clinic_segal2$`Pr(>F)`[1]
  )
)

### all -----
re_perm_all <- tibble(
  study=c("MBS","MLVS","MetaCardis","HCHS/SOL","Talmor_2022", "DIRECT-PLUS"),
  R2_all=c(re_perm_all_mbs$R2[1],
           re_perm_all_mlvs$R2[1],
           re_perm_all_ped$R2[1],
           re_perm_all_sol$R2[1],
           re_perm_all_segal1$R2[1],
           re_perm_all_direct$R2[1]
  ),
  P_all=c(re_perm_all_mbs$`Pr(>F)`[1],
          re_perm_all_mlvs$`Pr(>F)`[1],
          re_perm_all_ped$`Pr(>F)`[1],
          re_perm_all_sol$`Pr(>F)`[1],
          re_perm_all_segal1$`Pr(>F)`[1],
          re_perm_all_direct$`Pr(>F)`[1]
  )
)

### meta-analysis -----
re_perm_list <- 
  list(re_perm_all,re_perm_diet_spe,re_perm_spe,re_perm_food,
       re_perm_age,re_perm_sex,re_perm_bmi,
       re_perm_metf,re_perm_t2d)
re_perm_t <- re_perm_list %>% 
  reduce(full_join,by="study") %>% 
  mutate(SE_all=R2_all/(abs(qnorm(P_all/2))),
         SE_diet_spe=R2_diet_spe/(abs(qnorm(P_diet_spe/2))),
         SE_spe=R2_spe/(abs(qnorm(P_spe/2))),
         SE_food=R2_food/(abs(qnorm(P_food/2))),
         SE_age=R2_age/(abs(qnorm(P_age/2))),
         SE_sex=R2_sex/(abs(qnorm(P_sex/2))),
         SE_bmi=R2_bmi/(abs(qnorm(P_bmi/2))),
         SE_metf=R2_metf/(abs(qnorm(P_metf/2))),
         SE_t2d=R2_t2d/(abs(qnorm(P_t2d/2)))
  )

library(meta)
meta_r2_all <- metagen(TE=re_perm_t$R2_all,seTE = re_perm_t$SE_all,
                       fixed = TRUE)
meta_r2_diet_spe <- metagen(TE=re_perm_t$R2_diet_spe,
                            seTE = re_perm_t$SE_diet_spe,
                            fixed = TRUE)
meta_r2_spe <- metagen(TE=re_perm_t$R2_spe,seTE = re_perm_t$SE_spe,
                       fixed = TRUE)
meta_r2_food <- metagen(TE=re_perm_t$R2_food,seTE = re_perm_t$SE_food,
                        fixed = TRUE)
meta_r2_age <- metagen(TE=re_perm_t$R2_age,seTE = re_perm_t$SE_age,
                       fixed = TRUE)
meta_r2_sex <- metagen(TE=re_perm_t$R2_sex,seTE = re_perm_t$SE_sex,
                       fixed = TRUE)
meta_r2_bmi <- metagen(TE=re_perm_t$R2_bmi,seTE = re_perm_t$SE_bmi,
                       fixed = TRUE)
meta_r2_metf <- metagen(TE=re_perm_t$R2_metf,seTE = re_perm_t$SE_metf,
                        fixed = TRUE)
meta_r2_t2d <- metagen(TE=re_perm_t$R2_t2d,seTE = re_perm_t$SE_t2d,
                       fixed = TRUE)
meta_r2_t <- tibble(
  vars=c("All","Diet+Species","Species","Diet","Age","Sex","BMI","Metformin","T2D"), 
  R2=c(meta_r2_all$TE.fixed,meta_r2_diet_spe$TE.fixed,
       meta_r2_spe$TE.fixed,meta_r2_food$TE.fixed,
       meta_r2_age$TE.fixed,meta_r2_sex$TE.fixed,
       meta_r2_bmi$TE.fixed,meta_r2_metf$TE.fixed,
       meta_r2_t2d$TE.fixed
  ),
  P=c(meta_r2_all$pval.fixed,meta_r2_diet_spe$pval.fixed,
      meta_r2_spe$pval.fixed,meta_r2_food$pval.fixed,
      meta_r2_age$pval.fixed,meta_r2_sex$pval.fixed,
      meta_r2_bmi$pval.fixed,meta_r2_metf$pval.fixed,
      meta_r2_t2d$pval.fixed
  ))

## Heatmap -----
re_perm_all2 <- re_perm_all %>% 
  setnames(c("Study","R2","P")) %>% 
  mutate(vartype="All")
re_perm_diet_spe2 <- re_perm_diet_spe %>% 
  setnames(c("Study","R2","P"))%>% 
  mutate(vartype="Diet+Species")
re_perm_spe2 <- re_perm_spe %>% 
  setnames(c("Study","R2","P"))%>% 
  mutate(vartype="Species")
re_perm_food2 <- re_perm_food %>% 
  setnames(c("Study","R2","P"))%>% 
  mutate(vartype="Diet")
re_perm_clinic2 <- re_perm_clinic %>% 
  setnames(c("Study","R2","P"))%>% 
  mutate(vartype="Clinic")
re_perm_age2 <- re_perm_age %>% 
  setnames(c("Study","R2","P"))%>% 
  mutate(vartype="Age")
re_perm_sex2 <- re_perm_sex %>% 
  setnames(c("Study","R2","P"))%>% 
  mutate(vartype="Sex")
re_perm_bmi2 <- re_perm_bmi %>% 
  setnames(c("Study","R2","P"))%>% 
  mutate(vartype="BMI")
re_perm_metf2 <- re_perm_metf %>% 
  setnames(c("Study","R2","P"))%>% 
  mutate(vartype="Metformin")
re_perm_t2d2 <- re_perm_t2d %>% 
  setnames(c("Study","R2","P"))%>% 
  mutate(vartype="T2D")

re_perm_long <- 
  bind_rows(re_perm_all2,re_perm_diet_spe2,re_perm_spe2,re_perm_food2,
            re_perm_age2,re_perm_sex2,re_perm_bmi2,re_perm_metf2,re_perm_t2d2)
re_perm_long$R2 <- re_perm_long$R2 * 100

meta_r2_t2 <- meta_r2_t %>% 
  mutate(Study="Meta-analysis") %>% 
  setnames("vars","vartype")
meta_r2_t2$R2 <- meta_r2_t2$R2*100

(re_perm_long <- bind_rows(re_perm_long,meta_r2_t2))

re_perm_long <- re_perm_long %>% mutate(
  Study=case_when(Study=="MBS" ~ "NHSII",
                  Study=="MLVS" ~ "HPFS",
                  Study=="Talmor_2022" ~ "Talmor-Barkan_2022",
                  TRUE ~ Study),
  Study=factor(Study, levels=rev(c("Meta-analysis","DIRECT-PLUS","HCHS/SOL","HPFS",
                                   "MetaCardis","NHSII","Talmor-Barkan_2022"))),
  vartype=factor(vartype, levels=c("All","Diet+Species","Species","Diet","T2D",
                                   "Metformin","BMI","Sex","Age")),
  sig_sign=case_when(
    P<=0.001 ~ "**",
    P<0.05 ~ "*")
)

library(ggsci)
ggplot(re_perm_long,aes(y=Study,x=vartype,fill=R2))+
  geom_tile(color="grey60")+
  geom_text(aes(label=format(round(R2,2),nsmall=2)),size=4)+
  geom_text(aes(label=sig_sign),size=4,vjust=-0.3)+
  # scale_fill_continuous(na.value = "grey40")+
  scale_fill_material(palette = "deep-orange",alpha = 1,reverse = F,
                      na.value="grey")+
  theme_minimal()+
  # labs(y="Variable type")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        # axis.title.y = element_text(size = 16,color = "black"),
        axis.text.y = element_text(size = 14,color = "black"),
        axis.text.x = element_text(size = 14,color = "black",angle = 45,
                                   hjust = 1,vjust = 1),
        legend.position = "right",
        panel.grid = element_blank())
ggsave("/n/holylfs05/LABS/dwang_lab/Lab/T2D_multiomics/3.figures/fig2/figure2a.pdf",
       width = 8,height = 4)
