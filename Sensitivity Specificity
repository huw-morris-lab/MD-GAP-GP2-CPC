# Sensitivity & specificity analysis

# Loading library
library(tidyverse)

# Function to calculate diagnostic metrics from a confusion dataframe
calc_metrics <- function(df) {
  # Automatically identify the "positive" label
  positive_label <- df$Clinical[1] 
  patho_label <- df$Pathological[1]  
  
  mat <- xtabs(Count ~ Pathological + Clinical, data = df)
  
  # Extract counts using exact string matches from table
  TP <- mat[patho_label, positive_label]
  FP <- mat[setdiff(rownames(mat), patho_label), positive_label]
  FN <- mat[patho_label, setdiff(colnames(mat), positive_label)]
  TN <- mat[setdiff(rownames(mat), patho_label), setdiff(colnames(mat), positive_label)]
  
  sensitivity <- TP / (TP + FN)
  specificity <- TN / (TN + FP)
  PPV <- TP / (TP + FP)
  NPV <- TN / (TN + FN)
  
  data.frame(
    Label = positive_label,
    Sensitivity = round(sensitivity, 3),
    Specificity = round(specificity, 3),
    PPV = round(PPV, 3),
    NPV = round(NPV, 3)
  )
}

confusion_list <- list(
  PD = data.frame(
    Clinical = c("PD", "PD", "Not PD", "Not PD"),
    Pathological = c("LBD", "Not LBD", "LBD", "Not LBD"),
    Count = c(1654, 143, 143, 1606)
  ),
  DLB = data.frame(
    Clinical = c("DLB", "DLB", "Not DLB", "Not DLB"),
    Pathological = c("DLB", "Not DLB", "DLB", "Not DLB"),
    Count = c(215, 18, 1835, 1316)
  ),
  PSP = data.frame(
    Clinical = c("PSP", "PSP", "Not PSP", "Not PSP"),
    Pathological = c("PSP", "Not PSP", "PSP", "Not PSP"),
    Count = c(430, 61, 102, 2871)
  ),
  MSA = data.frame(
    Clinical = c("MSA", "MSA", "Not MSA", "Not MSA"),
    Pathological = c("MSA", "Not MSA", "MSA", "Not MSA"),
    Count = c(183, 61, 48, 3172)
  ),
  CBD = data.frame(
    Clinical = c("CBD", "CBD", "Not CBD", "Not CBD"),
    Pathological = c("CBD", "Not CBD", "CBD", "Not CBD"),
    Count = c(18, 58, 7, 3378)
  ),
  CTL = data.frame(
    Clinical = c("CTL", "CTL", "Not CTL", "Not CTL"),
    Pathological = c("CTL", "Not CTL", "CTL", "Not CTL"),
    Count = c(683, 62, 14, 2706)
  )
)

# Run and combine results
results <- do.call(rbind, lapply(confusion_list, calc_metrics))
print(results)


# Comparing PD vs PDD/DLB PPV
LBD <- df %>%
  filter(CLIN_DX == "PD" | CLIN_DX == "DLB" | CLIN_DX == "PDD") %>%
  mutate(PATH_LBD = case_when(
    PATHOLOGY %in% c("LBD", "DLB", "ILBD", "PD", "PDD") ~ 1,
    TRUE ~ 0))

LBD$GROUP <- case_when(
  LBD$CLIN_DX == "PD" ~ "PD",
  LBD$CLIN_DX %in% c("PDD", "DLB") ~ "PDD_DLB"
)



ppv_summary <- LBD %>%
  group_by(GROUP) %>%
  summarise(
    TP = sum(PATH_LBD == 1),     # Pathology confirms LBD
    TOTAL = n(),                      # All clinically diagnosed
    PPV = TP / TOTAL
  )

# Construct 2x2 table
table_ppv <- table(LBD$GROUP, LBD$PATH_LBD)
# Add labels
rownames(table_ppv) <- c("PD", "PDD_DLB")

print(table_ppv)
# Run Fisher’s exact test
fisher.test(table_ppv)
chisq_test(table_ppv)
