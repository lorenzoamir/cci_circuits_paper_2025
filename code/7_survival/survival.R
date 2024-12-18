library(survival)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(tibble)

MIN_PATIENTS <- 10

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
tissuefile <- args[1]
output_dir <- args[2]

# Read the survival data
cat('Reading survival data\n')
df <- read_csv(tissuefile, col_types = cols())

# Check and fix column names if needed
if (names(df)[1] == "") {
  names(df)[1] <- "PatientID"  # Rename if the first column has no name
}

# Set the first column as row names
df <- df %>% column_to_rownames(var = names(df)[1])

# Print number of patients (rows) in the dataset
cat('Number of patients:', nrow(df), '\n')

# Read motifs file and filter for specific types
motifs <- read_csv('/home/lnemati/pathway_crosstalk/results/crosstalk/all_ccc_complex_pairs/adj/motifs/tumor/motifs.csv', col_types = cols())
# DEBUG
#motifs <- motifs %>% filter(Type %in% c("3_clique") %>% pull(Interaction))
#motifs <- motifs %>% filter(Type %in% c("3_clique", "4_clique", "4_no_crosstalk")) %>% pull(Interaction)
motifs <- motifs %>% pull(Interaction)

# Subset the motifs for testing purposes
cat('Number of cliques (subset):', length(motifs), '\n')

# Extract ccc interactions from motifs by splitting and flattening the list
ccc <- unique(unlist(str_split(motifs, '&')))  # Split interactions by '&' and flatten the list
cat('Number of ccc interactions:', length(ccc), '\n')

survival_analysis <- function(interaction, df) {
  cat(interaction, '\n')

  genes <- unique(unlist(str_split(interaction, '[+&_]')))
  
  # Check if multiple conditions are present
  multiple <- length(unique(df$condition)) > 1

  # Select relevant columns (OS.time, OS, condition, and the genes of interest)
  cols <- c("OS.time", "OS", "condition", genes)
  if (!multiple) {
    cols <- setdiff(cols, "condition")
  }
  df <- df %>% select(all_of(cols))

  # Convert OS time to years
  df <- df %>% mutate(OS.time = OS.time / 365)

  # Create empty dataframes for high and low expression groups
  high_expression_group <- data.frame()
  low_expression_group <- data.frame()

  # If multiple tumors are present, split by condition
  print('Splitting')
  if (multiple) {
    for (condition in unique(df$condition)) {
      tissue_df <- df %>% filter(condition == condition)
      medians <- sapply(genes, function(gene) median(tissue_df[[gene]], na.rm = TRUE))
      
      # Create binary flags for above/below median for each gene
      above_median_flags <- sapply(genes, function(gene) tissue_df[[gene]] > medians[gene])
      below_median_flags <- sapply(genes, function(gene) tissue_df[[gene]] <= medians[gene])
      
      # Identify patients that are above the median for all genes
      above_all_genes <- rowSums(above_median_flags) == length(genes)
      below_all_genes <- rowSums(below_median_flags) == length(genes)
      
      # Split patients into high and low expression groups based on all genes
      high_expression_group <- bind_rows(high_expression_group, tissue_df[above_all_genes, ])
      low_expression_group <- bind_rows(low_expression_group, tissue_df[below_all_genes, ])
    }
  } else {
    # If only one tissue, split into high and low expression groups
    medians <- sapply(genes, function(gene) median(df[[gene]], na.rm = TRUE))
    
    # Create binary flags for above/below median for all genes
    above_median_flags <- sapply(genes, function(gene) df[[gene]] > medians[gene])
    below_median_flags <- sapply(genes, function(gene) df[[gene]] <= medians[gene])
    
    # Identify patients that are above the median for all genes
    above_all_genes <- rowSums(above_median_flags) == length(genes)
    below_all_genes <- rowSums(below_median_flags) == length(genes)
    
    # Split patients into high and low expression groups based on all genes
    high_expression_group <- df[above_all_genes, ]
    low_expression_group <- df[below_all_genes, ]
  }
                                 
  # Add a 'group' column to indicate high (1) or low (0) expression groups
  high_expression_group <- high_expression_group %>% mutate(group = 1)
  low_expression_group <- low_expression_group %>% mutate(group = 0)

  # Combine the high and low expression groups into one dataframe
  df <- bind_rows(high_expression_group, low_expression_group)
  
  # Remove any rows that are not part of the high or low expression groups
  print('Subsetting')
  df <- df %>% filter(group %in% c(0, 1))

  # Count the number of patients in each group
  n_above_all <- nrow(high_expression_group)
  n_below_all <- nrow(low_expression_group)

  # If there are too few patients in either group, list of NA values
  if (n_above_all < MIN_PATIENTS || n_below_all < MIN_PATIENTS) {
    print('Too few patients')
    return(list(
      hr = NA,
      n_patients_low = n_below_all,
      n_patients_high = n_above_all,
      logrank_pval = NA,
      concordance_index = NA,
      ci_low = NA,
      ci_high = NA,
      se = NA
    ))
  }

  # Make sure the group and condition columns are factors
  print('Converting')
  df$group <- as.factor(df$group)
  if (multiple) {
    df$condition <- as.factor(df$condition)
  }

  # Remove all columns that are not needed and run the Cox model
  if (multiple) {
    df <- df %>% select(OS.time, OS, group, condition)
  } else {
    df <- df %>% select(OS.time, OS, group)
  }

  result <- tryCatch({
    if (multiple) {
      model <- coxph(formula = Surv(OS.time, OS) ~ group + strata(condition), data = df)
    } else {
      model <- coxph(formula = Surv(OS.time, OS) ~ group, data = df)
    }

    # Extract the metrics
    hr <- exp(coef(model)["group1"])
    logrank_pval <- summary(model)$sctest["pvalue"]
    concordance <- summary(model)$concordance[1]
    print('Getting confidence intervals')
    ci_low <- exp(confint(model)["group1", 1])
    ci_high <- exp(confint(model)["group1", 2])
    se <- summary(model)$coefficients["group1", "se(coef)"]
    print(summary(model)$coefficients)
    # Propagate error:
    # hr = exp(coef) => err_hr = exp(coef) * se(coef) = hr * se(coef)
    se <- hr * se

    print('Model summary')
    print(summary(model))

    return(list(
      hr = hr,
      n_patients_low = n_below_all,
      n_patients_high = n_above_all,
      logrank_pval = logrank_pval,
      concordance_index = concordance,
      ci_low = ci_low,
      ci_high = ci_high,
      se = se
    ))
  }, warning = function(w) {
    print("Warning occurred during Cox model fitting. Returning NAs.")
    return(list(
      hr = NA,
      n_patients_low = n_below_all,
      n_patients_high = n_above_all,
      logrank_pval = NA,
      concordance_index = NA,
      ci_low = NA,
      ci_high = NA,
      se = NA
    ))
  })

  return(result)
}

# Test all ccc interactions
ccc_results <- data.frame()
for (interaction in ccc) {
  result <- survival_analysis(interaction, df)
  
  # Convert the result into a data frame, and add the interaction and type
  result_df <- data.frame(
    interaction = interaction,
    hr = result$hr,
    n_patients_low = result$n_patients_low,
    n_patients_high = result$n_patients_high,
    logrank_pval = result$logrank_pval,
    concordance_index = result$concordance_index,
    ci_low = result$ci_low,
    ci_high = result$ci_high,
    se=result$se,
    #zph = result$zph,
    type = "ccc"
  )
  
  # Bind the new row to the results dataframe
  ccc_results <- bind_rows(ccc_results, result_df)
}

# Test all crosstalk interactions
crosstalk_results <- data.frame()
for (interaction in motifs) {
  result <- survival_analysis(interaction, df)
  
  # Convert the result into a data frame, and add the interaction and type
  result_df <- data.frame(
    interaction = interaction,
    hr = result$hr,
    n_patients_low = result$n_patients_low,
    n_patients_high = result$n_patients_high,
    logrank_pval = result$logrank_pval,
    concordance_index = result$concordance_index,
    ci_low = result$ci_low,
    ci_high = result$ci_high,
    se=result$se,
    #zph = result$zph,
    type = "crosstalk"
  )
  
  # Bind the new row to the results dataframe
  crosstalk_results <- bind_rows(crosstalk_results, result_df)
}

# Combine the results
all_results <- bind_rows(ccc_results, crosstalk_results) %>%
  arrange(logrank_pval)  # Sorting by logrank pval

# Save results
tissue_name <- tools::file_path_sans_ext(basename(tissuefile))
output_file <- file.path(output_dir, paste0(tissue_name, '.csv'))

cat('Saving results to:', output_file, '\n')
write_csv(all_results, output_file)

cat('Done: survival.R')
