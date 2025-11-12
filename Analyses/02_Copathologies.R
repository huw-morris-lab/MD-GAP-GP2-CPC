# Analysis of associations between clinical characteristics and co-pathologies
# Last updated in November 2025

# Import genetic PCs generated via PLINK2 for covariates
DIR=path/to/directory/

ALL_PC <- read.table("$DIR/ALL_PC.txt")%>%
  dplyr::rename(gp2_id = V1) %>%
  dplyr::select(gp2_id, PC1, PC2, PC3)


df <- df %>%
  mutate(copath1 = case_when(
    (is.na(path_autopsy_dx2) | path_autopsy_dx2 == "") &    #path_autopsy_dx is the path diagnosis
      (is.na(path_autopsy_dx3) | path_autopsy_dx3 == "") &
      (is.na(path_autopsy_dx4) | path_autopsy_dx4 == "") ~ 0,
    
    path_autopsy_dx1 != "" & !is.na(path_autopsy_dx1) &
      path_autopsy_dx2 != "" & !is.na(path_autopsy_dx2) &
      (is.na(path_autopsy_dx3) | path_autopsy_dx3 == "") &
      (is.na(path_autopsy_dx4) | path_autopsy_dx4 == "") ~ 1,
    
    is.na(path_autopsy_dx4) | path_autopsy_dx4 == "" ~ 2,
    
    TRUE ~ NA_real_
  ))

df <- left_join(df,ALL_PC, by = "gp2_id")

# Age at onset
AAO <- lm(AGE_ONSET ~ copath1 + sex + PC1 + PC2 + PC3, data = df)
summary(AAO)
# Disease duration
DX <- lm(DX1 ~ copath1 + sex + + AGE_ONSET + PC1 + PC2 + PC3, data = df)
summary(DX)
