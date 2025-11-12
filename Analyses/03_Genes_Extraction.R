# Script to extract genes of interest
# Last updated in August 2025

# PART 1
---
# STEP 1: use plink2 to convert raw genotypes from each ancestry and brain bank into txt files after doing QC (--geno and --mind)
#.txt files used: "SNPlist.txt"
# The /txt file is a SNP list (chr start end) set up by a combination of Clinvar path/ likely path variants identified from WGS in the MDGAP cohort, filtered by a clinical or path diagnosis of movement disorders #and other variants of interest in PD genes)
# format: #CHROM START END 
# the NBA files are saved under each ancestry, for each brainbank
# output file: {ancestry}_{brainbank}_qc.txt (files are saved in brainbank folders under RDS/Lesley/CPC; the original files were derived from RDS/GP2_NBA/)
# last run 12 May by Lesley

# Set working directory
base_dir="path/to/folder"
cd "$base_dir" || exit

# Define ancestries and brainbanks
ancestries=("AAC" "AFR" "AJ" "AMR" "CAH" "CAS" "EAS" "EUR" "FIN" "MDE" "SAS")
brainbanks=("CAMBRIDGE" "EBB" "KINGS" "NWCSTL" "QSBB" "VBB" "DEMENTIA" "IBB" "MBB" "OXF" "SW" "BBDP" "SYD")

# Loop through each brainbank and ancestry
for brainbank in "${brainbanks[@]}"; do
    for ancestry in "${ancestries[@]}"; do
        echo "Processing $ancestry in $brainbank..."

        ancestry_dir="base_dir${brainbank}/raw_genotypes/${ancestry}"
        
        # Detect the base name of the pfile (e.g., MDGAP-QSBB_EUR_release10)
        pfile_base=$(ls "${ancestry_dir}"/*.pgen 2>/dev/null | head -n 1 | sed 's/\.pgen$//')

        if [ -z "${pfile_base}" ]; then
            echo "No .pgen file found for ${ancestry} in ${brainbank}, skipping..."
            continue
        fi

        # Create output directory for this brainbank/ancestry
        output_dir="${base_dir}/${brainbank}/${ancestry}"
        mkdir -p "$output_dir"

        # Define output base
        outbase="${output_dir}/${ancestry}_${brainbank}_qc"

        # Run PLINK2 steps
        plink2 --pfile "${pfile_base}" --geno 0.05 --make-pgen --out "${outbase}_geno"
        plink2 --pfile "${outbase}_geno" --mind 0.05 --make-pgen --out "${outbase}_mind"
        plink2 --pfile "${outbase}_mind" --extract range SNPlistnew_9.txt --recode vcf id-paste=iid --out "${outbase}"

        # Convert VCF to TXT
        sed '/^##/ d' < "${outbase}.vcf" > "${outbase}.txt"

        echo "Finished $ancestry in $brainbank, results in $output_dir"
    done
done

echo "All ancestries and brainbanks processed."

---
# STEP 2: Step 2: process the output in R; recode to genotypes to 0,1 and 2; transpose the data so that mutation probes are the column headings and IDs are the rows
#.txt files used: "mutation_names_14May_new.txt" in the format: NBA_name (Illumina probe name from the NeuroBooster20042459_A1.csv document (RDS/NeuroBooster) and SNP_id (easy to understand gene_proteinchange.suffix name for each probe) 
#output files:
#"summary_all_results_snps_transpose_all_cohorts_14May.csv"
#"summary_positive_results_snps_transpose_all_cohorts_14May.csv"
#"summary_results_NBA_name_snps_all_cohorts_14May.csv"
#Last run: 14 May 14:20

# --- Load Libraries ---
library(data.table)
library(tidyverse)
library(readxl)
library(dplyr)

# --- Load Mutation Names ---
mutation_names <- fread("$base_dir/mutation_names_14May_new.txt")

# --- Set Root Directory ---
brainbanks <- c("CAMBRIDGE", "EBB", "KINGS", "NWCSTL", "QSBB", "VBB", "DEMENTIA", "IBB", "MBB", "OXF", "SW")
ancestries <- c("AAC", "AFR", "AJ", "AMR", "CAH", "CAS", "EAS", "EUR", "FIN", "MDE", "SAS")

# --- Prepare Output Containers ---
all_data_t <- list()
positive_data_t <- list()
variant_summary <- list()

# --- Process Each File ---
for (brainbank in brainbanks) {
  for (ancestry in ancestries) {
    file_path <- file.path(base_dir, brainbank, ancestry, paste0(ancestry, "_", brainbank, "_qc.txt"))
    
    if (!file.exists(file_path)) next
    
    cat("Processing", ancestry, brainbank, "\n")
    
    # --- Read PLINK txt ---
    plink_vcf <- fread(file_path)
    names(plink_vcf)[3] <- "NBA_name"
    names(plink_vcf)[1] <- "CHROM"
    
    # Backup original NBA_name before merge
    plink_vcf$original_ID <- plink_vcf$NBA_name
    
    # --- Merge with mutation names ---
    plink_data <- merge(plink_vcf, mutation_names, by = "NBA_name", all.x = TRUE)
    
    # Replace NBA_name with SNP_id if available, otherwise keep original
    plink_data$NBA_name <- ifelse(is.na(plink_data$SNP_id), plink_data$original_ID, plink_data$SNP_id)
    
    # Drop unnecessary columns
    plink_data$SNP_id <- NULL
    plink_data$original_ID <- NULL
    
    # Drop PLINK metadata columns (columns 2 to 9)
    plink_data <- subset(plink_data, select = -c(2:9))
    
    # --- Convert genotype strings to numeric ---
    plink_data[plink_data == "0/0"] <- 0
    plink_data[plink_data == "0/1"] <- 1
    plink_data[plink_data == "1/1"] <- 2
    plink_data[plink_data == "./."] <- NA
    
    # --- Store variant summary (non-zero genotypes) ---
    any_vars <- filter_at(plink_data, vars(-NBA_name), any_vars(. >= 1))
    any_vars$brainbank <- brainbank
    any_vars$ancestry <- ancestry
    variant_summary[[paste0(brainbank, "_", ancestry)]] <- any_vars

    # --- Transpose data ---
    sample_ids <- names(plink_data)[-1]
    transposed_matrix <- t(as.matrix(plink_data[, -1, with = FALSE]))
    transposed <- as.data.table(transposed_matrix)
    setnames(transposed, plink_data$NBA_name)
    transposed[, sample_id := sample_ids]
    
    # Add metadata
    transposed$brainbank <- brainbank
    transposed$ancestry <- ancestry
    
    # Reorder columns: sample_id, brainbank, ancestry first, then mutations
    setcolorder(
      transposed,
      c("sample_id", "brainbank", "ancestry", setdiff(names(transposed), c("sample_id", "brainbank", "ancestry")))
    )
    
    # --- Filter positive genotype samples ---
    any_vars_t <- filter_at(transposed, vars(-sample_id, -brainbank, -ancestry), any_vars(. >= 1))
    any_vars_t$brainbank <- brainbank
    any_vars_t$ancestry <- ancestry

    # --- Store ---
    all_data_t[[paste0(brainbank, "_", ancestry)]] <- transposed
    positive_data_t[[paste0(brainbank, "_", ancestry)]] <- any_vars_t
  }
}

# --- Combine All Results ---
summary_all <- bind_rows(all_data_t)
summary_positive <- bind_rows(positive_data_t)
summary_variant <- bind_rows(variant_summary)

# --- Save Outputs with 14May Naming ---
fwrite(summary_all, file = "summary_all_results_snps_transpose_all_cohorts_14May.csv")
fwrite(summary_positive, file = "summary_positive_results_snps_transpose_all_cohorts_14May.csv")
fwrite(summary_variant, file = "summary_results_NBA_name_snps_all_cohorts_14May.csv")

---

# PART 2
---
# PART 2 looks for the same variants in WGS that was already extracted from NBA
#Then we compare the WGS and NBA output 
---
# Extract the same variant list (SNPlistnew_9.txt) from WGS data 
# output files: "UCL_WGS_Feb2025_SNPs.txt"

#Using the terminal code the following:
#path to directory where joint call CVF is saved
DIR=/path/to/data
#Check if there are .psam, .pvar and .pgen files in the UCL_VCF folder ready for use before running this part
nohup plink2 --vcf $DIR/UCL_WGS_Feb2025.vcf.gz --vcf-half-call missing --make-pgen --out UCL_WGS_Feb2025 &

#Last run: 14 May 
plink2 \
  --pfile $DIR/UCL_WGS_Feb2025 \
  --extract range $DIR/SNPlistnew_9.txt \
  --recode vcf id-paste=iid \
  --out UCL_WGS_Feb2025_SNPs


sed '/^##/ d' < UCL_WGS_Feb2025_SNPs.vcf > UCL_WGS_Feb2025_SNPs.txt 

---
#Convert the .txt file which is derived from WGS data into a readable .csv output 
#.txt files used:"mutation_names_build38new_14May.txt" in the format: SNP_id (same as before, an easy to understand gene_proteinchange.suffix name) and chr_pos (chr_pos_ref_alt) 
# output files:"summary_all_results_WGS_snps_transpose_14May.csv" and "summary_all_positive_results_WGS_snps_transpose_14May.csv"


#---Load libraries---####
library(data.table)
library(tidyverse)
library(readxl)
library(dplyr)

#Read in list of mutation names (must include NBA_id and variant name that i want it to be called)
mutation_names_build38new <- fread("$DIR/mutation_names_build38new_14May.txt")

#---Read in plink raw file---####
plink_vcf <- fread("UCL_WGS_Feb2025_SNPs.txt") #this will be the txt output of the plink commands 


#Rename NeuroChip SNP ID with mutation names (from VCF file - "--recode vcf")
names(plink_vcf)[3]<- "chr_pos" #same as the first column NBA_id in the mutation names file
names(plink_vcf)[1]<- "CHROM" #To remove '#' symbol
plink_data <- merge(plink_vcf, mutation_names_build38new, by = "chr_pos") 
plink_data$chr_pos <- plink_data$SNP_id
plink_data$SNP_id <- NULL
plink_data <- subset(plink_data, select = -c(2:9)) #check that no extra columns are deleted (keep SNPid and sample ids)
plink_data[plink_data == "0/0"] <- 0
plink_data[plink_data == "0/1"] <- 1
plink_data[plink_data == "1/1"] <- 2
plink_data[plink_data == "./."] <- NA

#This generates a list of all variants with a genotype != "0/0 in all individuals (columns)
#any_vars <- filter_at(plink_data, vars(-chr:pos), any_vars(.>=1)) 
#write_excel_csv(any_vars, "summary_results_MDGAP_PFP_SNPs.csv")

#This generates a list of all individuals (rows) with a variant != "0/0" and is more useful
plink_data_t <- setnames(plink_data[, data.table(t(.SD), keep.rownames=TRUE), .SDcols=-"chr_pos"], 
                         plink_data[, c('sample_id', chr_pos)])[] #change to the name of the NBA id column 

any_vars_t <- filter_at(plink_data_t, vars(-sample_id), any_vars(.>=1))
#ensure unique names by .x in txt file
write_excel_csv(any_vars_t, "summary_positive_results_WGS_snps_transpose_14May.csv")
write_excel_csv(plink_data_t, "summary_all_results_WGS_snps_transpose_14May.csv")

