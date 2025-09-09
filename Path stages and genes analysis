########### Stacked bar chart ###########
### LBD SUBTYPE ###
library(dplyr)
library(stringr)
library(ggplot2)
library(RColorBrewer)
library(scales)

df <- df %>%
  filter(
    !is.na(LBD_Subtype),
    !LBD_Subtype %in% c("No information available", "Not applicable", "Unclassifiable")
  ) %>%
  count(GENE_TYPE, LBD_Subtype) %>%
  dplyr::rename(Count = n) %>%
  group_by(GENE_TYPE) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

df$LBD_Subtype <- factor(df$LBD_Subtype, levels = c("None", "Amygdala", "Brainstem", "Limbic", "Neocortical"))

# Define custom colors for McKeith stages 
custom_colors <- c(
  "None" = "#87CEEB",  # sky blue
  "Amygdala" = "#4682B4",  # steel blue
  "Brainstem" = "#FFE4E1",  # misty rose
  "Limbic" = "#FFC6B3",  
  "Neocortical" = "#FF7F50"  # coral   
)

# Plot
p <- ggplot(df, aes(x = GENE_TYPE, y = Proportion, fill = LBD_Subtype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Percent_Label), 
            position = position_stack(vjust = 0.5), 
            size = 3, color = "black") +
  scale_fill_manual(values = custom_colors, name = "LBD Subtype") +
  labs(
    x = "Mutation Group",
    y = "Proportion of Cases",
    title = "Distribution of LBD subtype by Mutation Group Across all Pathological Diagnoses"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )


### BRAAK NFT ###
# Group into Braak NFT categories including stage 0
df <- df %>%
  filter(!is.na(BRAAK_NFT)) %>%
  count(GENE_TYPE, BRAAK_NFT, name = "Count") %>%
  group_by(GENE_TYPE) %>%
  mutate(Proportion = Count / sum(Count)) %>%
  ungroup()

# Set stage order: from no pathology to severe
df$BRAAK_NFT <- factor(df$BRAAK_NFT,
                                      levels = c("0", "I–II", "III–IV", "V–VI"))

# Add percent labels
df <- df %>%
  mutate(Percent_Label = percent(Proportion, accuracy = 1))

# Define custom colors
braak_colors <- c(
  "0" = "#87CEEB",      # pale blue (no pathology)
  "I–II" = "#FFE4E1",   
  "III–IV" = "#FFC6B3",  
  "V–VI" =  "#FF7F50"    
)

# Plot
p_braak <- ggplot(df, aes(x = GENE_TYPE, y = Proportion, fill = BRAAK_NFT)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Percent_Label),
            position = position_stack(vjust = 0.5),
            size = 3, color = "black") +
  scale_fill_manual(values = braak_colors, name = "Braak NFT Stage") +
  labs(
    x = "Mutation Group",
    y = "Proportion of Cases",
    title = "Distribution of Braak NFT Stages by Mutation in Individuals with Lewy Body Diseases"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )

p_braak

########### Regressions ###########
library(MASS) 
# Set reference 
df$Mutation_status <- relevel(as.factor(df$Mutation_status), ref = "NO MUTATION") 

# Fit ordinal logistic regression model using proportional odds logistic regression
model <- polr(as.factor(Mckeith_coded) ~ Mutation_status + sex + age_at_death + DX1, 
              data = df, method = "logistic")
# The outcome Mckeith_coded is treated as an ordinal factor (e.g., stage 0 to 4).

# Extract the model summary 
model_summary <- summary(model)

########### SURVIVAL ###########
library(survival)
library(survminer)
library(broom)

# Make GENE_TYPE a factor with correct levels matching the palette
df$GENE_TYPE <- factor(df$GENE_TYPE, levels = c("GBA1 PD RISK", "GBA1 GD CAUSING","LRRK2", "NO MUTATION"))

# Create event variable
df$event <- 1 


# Create survival object
surv_obj <- Surv(time = df$DX1, event = df$event)
surv_obj
# Fit Kaplan-Meier model
km_fit <- survfit(surv_obj ~ GENE_TYPE + AGE_ONSET, data = df)

# Convert to dataframe using broom::tidy
km_df <- tidy(km_fit)

# Custom colors (named to match factor levels exactly)
custom_colors <- c("GBA1 PD RISK" = "#FFA07A", "GBA1 GD CAUSING" = "#E5673C", "LRRK2" = "#FFE4E1","NO MUTATION" = "#4682B4")
df$GENE_TYPE <- factor(df$GENE_TYPE, levels = names(custom_colors))
surv_obj <- Surv(time = df$DX1, event = df$event)
km_fit <- survfit(surv_obj ~ GENE_TYPE, data = df)
km_df <- tidy(km_fit)

# Clean strata and match colors
km_df$GENE_TYPE <- gsub("GENE_TYPE=", "", km_df$strata)
km_df$GENE_TYPE <- trimws(km_df$GENE_TYPE)
km_df$GENE_TYPE <- factor(km_df$GENE_TYPE, levels = names(custom_colors))

km_plot <- ggplot(km_df, aes(x = time, y = estimate, color = GENE_TYPE)) +
  geom_step(size = 1) +
  scale_color_manual(values = custom_colors) +
  xlim(0, 40) +
  labs(
    x = "Years since onset",
    y = "Survival Probability",
    color = "Gene Mutation",
    title = "Kaplan-Meier Survival by Mutation Status"
  ) +
  theme_minimal()

km_plot


# Set "NO MUTATION" as reference category
df$GENE_TYPE <- relevel(as.factor(df$GENE_TYPE), ref = "NO MUTATION")

# Kaplan-Meier log-rank test
survdiff(Surv(DX1, event) ~ GENE_TYPE, data = df)

# Cox proportional hazards model (adjusted for AGE_ONSET)
cox_fit <- coxph(Surv(DX1, event) ~ GENE_TYPE + AGE_ONSET, data = df)

# Summary output with HRs, 95% CI, p-values
summary(cox_fit)
