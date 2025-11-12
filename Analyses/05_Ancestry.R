# Ancestry analysis examining associations between genetically defined ancestry and pathology
# Last updated in November 2025

# Read ancestry data generated via Genotools
# For demographics
ANCESTRY <- read.table("$DIR/all_ids_with_ancestry.txt", header = T) %>%
  dplyr::rename(gp2_id = FID) %>%
  dplyr::select(-IID, -BrainBank)

# df is dataframe including all clinical and path data
df <- df %>%
  left_join(ANCESTRY, by = "gp2_id")

# Define ancestry groups
ancestries <- list(
  AFR = c("AAC", "AFR"),
  AJ = "AJ",
  AMR = "AMR",
  CAH = "CAH",
  CAS = "CAS",
  EAS = "EAS",
  MDE = "MDE",
  SAS = "SAS",
  EUR = "EUR"
)

# Function to summarize per ancestry
summarise_ancestry <- function(df, ancestry_label, ancestry_codes) {
  df_sub <- df %>%
    filter(Ancestry %in% ancestry_codes) %>%
    filter(!is.na(age_at_death), !is.na(PATHOLOGY), !is.na(sex))
  
  list(
    Group = ancestry_label,
    PathologyCount = plyr::count(df_sub$PATHOLOGY),
    SexCount = plyr::count(df_sub$sex),
    MeanAge = mean(df_sub$age_at_death),
    SdAge = sd(df_sub$age_at_death)
  )
}

# Apply to all ancestries
ancestry_results <- lapply(names(ancestries), function(label) {
  summarise_ancestry(df, label, ancestries[[label]])
})

names(ancestry_results) <- names(ancestries)  # Label each result

# Comparing only ancestries with more than 10 cases
AJ_SAS_EUR <- df %>%
  filter(Ancestry == "EUR" |Ancestry == "AJ" | Ancestry == "SAS")

# Kruskal-Wallis
kruskal.test(age_at_death ~ Ancestry, data = AJ_SAS_EUR)

# Pairwise Wilcoxon
pairwise.wilcox.test(AJ_SAS_EUR$age_at_death, AJ_SAS_EUR$Ancestry,
                     p.adjust.method = "bonferroni")

# Comparing pathology type by ancestry
df <- df %>%
  filter(GENE_TYPE == "NO MUTATION",
         PATHOLOGY %in% c("LBD", "PSP")) %>%
  count(Ancestry, PATHOLOGY) %>%
  pivot_wider(names_from = PATHOLOGY1, values_from = n, values_fill = 0)

# Run chi-squared test
chi_res <- chisq.test(df)

prop_table <- prop.table(df, margin = 1)
print(prop_table)

pairwise_AJ_SAS <- df[c("AJ", "SAS"), ]
fisher.test(pairwise_AJ_SAS)

pairwise_EUR_SAS <- df[c("EUR", "SAS"), ]
fisher.test(pairwise_EUR_SAS)

# Create barplot 

png("$DIR/Ancestry.png", width = 8, height = 6, units = "in", res = 600, bg = "white")

bp <- barplot(
  t(prop_table),
  beside = TRUE,
  col = c("#CCE5F6", "#FF7F50"),
  ylim = c(0, 1),
  xlim = c(0, 10),  # adjust manually if needed
  main = "Proportions of LBD vs PSP",
  ylab = "Proportion"
)
dev.off()















