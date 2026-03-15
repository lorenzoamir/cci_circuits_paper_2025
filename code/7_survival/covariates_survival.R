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

MIN_PATIENTS <- 10

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)

covariates_dir = '/home/lnemati/pathway_crosstalk/results/survival/covariates/'

# Find matching covariates file
all_covariates_files <- list.files(covariates_dir, full.names = TRUE)
tissue_name <- tools::file_path_sans_ext(basename(tissuefile))
tissue_name_pattern <- gsub(" ", "_", tissue_name, fixed = TRUE)

matching_file <- all_covariates_files[
    grepl(tissue_name_pattern, all_covariates_files, ignore.case = TRUE)
]

if (length(matching_file) == 0) {
    cat("No matching covariates file found for tissue:", tissue_name, "\n") 
    cat("Done: covariates_survival.R - no covariates file")
    stop("Done: covariates_survival.R - no covariates file")
} else if (length(matching_file) > 1) {
    cat("Multiple matching covariates files found for tissue:", tissue_name, "\n")
    cat("Done: covariates_survival.R - multiple covariates files")
    stop("Done: covariates_survival.R - multiple covariates files")
}

cat('Reading covariates from:', matching_file, '\n')
covariates_df <- read_csv(matching_file, col_types = cols())
circuits <- covariates_df %>% pull(interaction)             
                     
# Read the survival data
cat('Reading survival data\n')
df <- read_csv(tissuefile, col_types = cols())

# Check and fix column names if needed
if (names(df)[1] == "") {
  names(df)[1] <- "sample"  # Rename if the first column has no name
}

# Set the first column as row names
df <- df %>% column_to_rownames(var = names(df)[1])

# The ids are like this TCGA-3M-AB46-01, the last number after last - represents samples, remove it to get patient id
df$patient <- sapply(strsplit(rownames(df), '-'), function(x) paste(x[1:4], collapse = '-'))

# Print number of patients (rows) in the dataset
cat('Number of patients:', nrow(df), '\n')
                         
# Find matching covariates file
all_covariates_files <- list.files(covariates_dir, full.names = TRUE)
tissue_name <- tools::file_path_sans_ext(basename(tissuefile))
tissue_name_pattern <- gsub(" ", "_", tissue_name, fixed = TRUE)

matching_file <- all_covariates_files[
    grepl(tissue_name_pattern, all_covariates_files, ignore.case = TRUE)
]

if (length(matching_file) == 0) {
    cat("No matching covariates file found for tissue:", tissue_name, "\n") 
    cat("Done: covariates_survival.R - no covariates file")
    stop("Done: covariates_survival.R - no covariates file")
} else if (length(matching_file) > 1) {
    cat("Multiple matching covariates files found for tissue:", tissue_name, "\n")
    cat("Done: covariates_survival.R - multiple covariates files")
    stop("Done: covariates_survival.R - multiple covariates files")
}

cat('Reading covariates from:', matching_file, '\n')
covariates_df <- read_csv(matching_file, col_types = cols())
circuits <- covariates_df %>% pull(interaction)             
                     
# Extract ccc interactions from circuits by splitting and flattening the list
ccc <- unique(unlist(str_split(circuits, '&')))  # Split interactions by '&' and flatten the list
cat('Number of ccc interactions:', length(ccc), '\n')

multiple <- length(unique(df$condition)) > 1

if (multiple) {
  df <- df %>%
    group_by(condition) %>%
    mutate(
        binary_age = as.integer(age_at_diagnosis > median(age_at_diagnosis, na.rm = TRUE)),
        gender = as.integer(gender == "Female"),
        Lymphocytes = as.integer(Lymphocytes > median(Lymphocytes, na.rm = TRUE)),
        Neutrophils = as.integer(Neutrophils > median(Neutrophils, na.rm = TRUE)),
        Eosinophils = as.integer(Eosinophils > median(Eosinophils, na.rm = TRUE)),
        Mast_Cells = as.integer(`Mast_Cells` > median(`Mast_Cells`, na.rm = TRUE)),
        Dendritic_Cells = as.integer(`Dendritic_Cells` > median(`Dendritic_Cells`, na.rm = TRUE)),
        Macrophages = as.integer(Macrophages > median(Macrophages, na.rm = TRUE)),
        tumor_purity = as.integer(tumor_purity > median(tumor_purity, na.rm = TRUE))
    ) %>%
    ungroup()
} else {
  df <- df %>%
    mutate(
        binary_age = as.integer(age_at_diagnosis > median(age_at_diagnosis, na.rm = TRUE)),
        gender = as.integer(gender == "Female"),
        Lymphocytes = as.integer(Lymphocytes > median(Lymphocytes, na.rm = TRUE)),
        Neutrophils = as.integer(Neutrophils > median(Neutrophils, na.rm = TRUE)),
        Eosinophils = as.integer(Eosinophils > median(Eosinophils, na.rm = TRUE)),
        Mast_Cells = as.integer(`Mast_Cells` > median(`Mast_Cells`, na.rm = TRUE)),
        Dendritic_Cells = as.integer(`Dendritic_Cells` > median(`Dendritic_Cells`, na.rm = TRUE)),
        Macrophages = as.integer(Macrophages > median(Macrophages, na.rm = TRUE)),
        tumor_purity = as.integer(tumor_purity > median(tumor_purity, na.rm = TRUE))
    )
}

survival_analysis <- function(interaction, df, covariates = NULL) {
  genes <- unique(unlist(str_split(interaction, '[+&_]')))

  cols <- c("patient", "OS.time", "OS", "condition", genes, covariates)
  if (!multiple) {
    cols <- setdiff(cols, "condition")
  }
  df <- df %>% select(all_of(cols))

  df <- df %>% mutate(OS.time = OS.time / 365)

  high_expression_group <- data.frame()
  low_expression_group <- data.frame()

  if (multiple) {
    for (condition in unique(df$condition)) {
      tissue_df <- df %>% filter(.data$condition == condition)
      medians <- sapply(genes, function(gene) median(tissue_df[[gene]], na.rm = TRUE))

      above_median_flags <- sapply(genes, function(gene) tissue_df[[gene]] > medians[gene])
      below_median_flags <- sapply(genes, function(gene) tissue_df[[gene]] <= medians[gene])

      above_all_genes <- rowSums(above_median_flags) == length(genes)
      below_all_genes <- rowSums(below_median_flags) == length(genes)

      high_expression_group <- bind_rows(high_expression_group, tissue_df[above_all_genes, ])
      low_expression_group <- bind_rows(low_expression_group, tissue_df[below_all_genes, ])
    }
  } else {
    medians <- sapply(genes, function(gene) median(df[[gene]], na.rm = TRUE))
    above_median_flags <- sapply(genes, function(gene) df[[gene]] > medians[gene])
    below_median_flags <- sapply(genes, function(gene) df[[gene]] <= medians[gene])

    above_all_genes <- rowSums(above_median_flags) == length(genes)
    below_all_genes <- rowSums(below_median_flags) == length(genes)

    high_expression_group <- df[above_all_genes, ]
    low_expression_group <- df[below_all_genes, ]
  }

  patients_in_both_groups <- intersect(high_expression_group$patient, low_expression_group$patient)
  high_expression_group <- high_expression_group %>% filter(!patient %in% patients_in_both_groups)
  low_expression_group <- low_expression_group %>% filter(!patient %in% patients_in_both_groups)

  high_expression_group <- high_expression_group %>% mutate(group = 1)
  low_expression_group <- low_expression_group %>% mutate(group = 0)

  df <- bind_rows(high_expression_group, low_expression_group)
  df <- df %>% select(-patient)

  high_expression_ids <- paste(rownames(high_expression_group), collapse = ";")
  low_expression_ids <- paste(rownames(low_expression_group), collapse = ";")

  df <- df %>% filter(group %in% c(0, 1))

  n_above_all <- nrow(high_expression_group)
  n_below_all <- nrow(low_expression_group)

  if (n_above_all < MIN_PATIENTS || n_below_all < MIN_PATIENTS) {
    return(list(
      model_full = NA,
      model_reduced = NA,
      n_patients_low = n_below_all,
      n_patients_high = n_above_all
    ))
  }

  df$group <- as.factor(df$group)
  # Convert covariates to factors 
  if (!is.null(covariates)) df[covariates] <- lapply(df[covariates], as.factor)
  if (multiple) df$condition <- as.factor(df$condition)

  if (!is.null(covariates) && length(covariates) > 0) {
    valid_covariates <- c()
    for (cov in covariates) {
      if (!cov %in% names(df) || all(is.na(df[[cov]]))) next
      counts <- table(df[[cov]])
      if (all(counts >= MIN_PATIENTS)) valid_covariates <- c(valid_covariates, cov)
    }
    covariates <- if (length(valid_covariates) > 0) valid_covariates else NULL
  }
 
  # Block any covariate that has less than 2 categories
  if (!is.null(covariates)) {
   covariates <- covariates[sapply(covariates, function(cov) length(unique(df[[cov]])) >= 2)]
  }

  # If some covariate has less than MIN_PATIENTS in any category, exclude it from the model
  covariates_to_exclude <- c()

  for (cov in covariates) {
    counts <- table(df[[cov]])
    if (any(counts < MIN_PATIENTS)) {
      covariates_to_exclude <- c(covariates_to_exclude, cov)
    }
  }
  covariates <- setdiff(covariates, covariates_to_exclude)

  # Prepare formulas
  cov_string <- if (!is.null(covariates)) paste(covariates, collapse = " + ") else NULL

  if (multiple) {
    full_rhs <- paste(c("group", cov_string, "strata(condition)"), collapse = " + ")
    reduced_rhs <- paste(c(cov_string, "strata(condition)"), collapse = " + ")
  } else {
    full_rhs <- paste(c("group", cov_string), collapse = " + ")
    reduced_rhs <- cov_string
  }

  formula_full <- as.formula(paste("Surv(OS.time, OS) ~", full_rhs))
  formula_reduced <- if (!is.null(reduced_rhs)) as.formula(paste("Surv(OS.time, OS) ~", reduced_rhs)) else NULL

  result <- tryCatch({
    model_full <- coxph(formula_full, data = df)
    model_reduced <- if (!is.null(formula_reduced)) coxph(formula_reduced, data = df) else NULL

    list(
      model_full = model_full,
      model_reduced = model_reduced,
      n_patients_low = n_below_all,
      n_patients_high = n_above_all,
      covariates_used = covariates
    )
  }, warning = function(w) {
    message("Warning for motif ", interaction, ": ", w$message)
    # Try to fit the model excluding covariates one by one, start from last covariate in the list
    # Exclude last, if it fails, exclude last two etc.
    result <- NULL
    if (!is.null(covariates) && length(covariates) > 0) {
      for (i in seq_along(covariates)) {
        covariates_subset <- covariates[1:(length(covariates) - i)]
        cov_string_subset <- if (length(covariates_subset) > 0) paste(covariates_subset, collapse = " + ") else NULL
        full_rhs_subset <- paste(c("group", cov_string_subset, if (multiple) "strata(condition)"), collapse = " + ")
        reduced_rhs_subset <- paste(c(cov_string_subset, if (multiple) "strata(condition)"), collapse = " + ")
        formula_full_subset <- as.formula(paste("Surv(OS.time, OS) ~", full_rhs_subset))
        formula_reduced_subset <- if (!is.null(reduced_rhs_subset)) as.formula(paste("Surv(OS.time, OS) ~", reduced_rhs_subset)) else NULL
        result <- tryCatch({
          model_full_subset <- coxph(formula_full_subset, data = df)
          model_reduced_subset <- if (!is.null(formula_reduced_subset)) coxph(formula_reduced_subset, data = df) else NULL

          list(
            model_full = model_full_subset,
            model_reduced = model_reduced_subset,
            n_patients_low = n_below_all,
            n_patients_high = n_above_all,
            covariates_used = covariates_subset
          )
        }, warning = function(w) {
          message("Warning for motif ", interaction, " with covariates ", paste(covariates_subset, collapse = ", "), ": ", w$message)
          return(NULL)
        }, error = function(e) {
          message("Error for motif ", interaction, " with covariates ", paste(covariates_subset, collapse = ", "), ": ", e$message)
          return(NULL)
        })
        if (!is.null(result)) break
      }
    }
    return(result)
  }, error = function(e) {
    message("Error for motif ", interaction, ": ", e$message)
    return(NULL)
  })    
}
    
original_cols <- c("interaction", "tissue", "hr", "concordance_index", "logrank_pval", "motif")

cat('Subsetting columns in the original dataframe\n')
covariates_df <- covariates_df %>% select(all_of(original_cols))

cat('Initializing covariate results dataframe\n')
cov_results <- data.frame()

pb <- txtProgressBar(min = 0, max = length(circuits), style = 3)

# set covariates, 
covariates <- c("binary_tumor_stage", "binary_age", "gender", "Lymphocytes", "Neutrophils", "Eosinophils", "Mast_Cells", "Dendritic_Cells", "Macrophages", "tumor_purity")

# all covariates are binary, for each covariate get the median survival time of each category,
# sort covariates by the absolute difference of median survival time between categories
covariate_medians <- data.frame(covariate = covariates, median_diff = NA)
for (cov in covariates) {
  if (!cov %in% names(df) || all(is.na(df[[cov]]))) {
    covariate_medians <- covariate_medians %>% filter(covariate != cov)
    next
  }
  median_0 <- median(df$OS.time[df[[cov]] == 0], na.rm = TRUE)
  median_1 <- median(df$OS.time[df[[cov]] == 1], na.rm = TRUE)
  covariate_medians <- covariate_medians %>% mutate(median_diff = ifelse(covariate == cov, abs(median_0 - median_1), median_diff))
}
covariate_medians <- covariate_medians %>% arrange(desc(median_diff))
covariates <- covariate_medians$covariate   
cat('Covariates:', paste(covariates, collapse = ", "), '\n')

for (i in seq_along(circuits)) {
  interaction <- circuits[i]
  interaction_row <- data.frame(interaction = interaction)
    
  # Run survival for full and reduced model
  model_list <- tryCatch({
      survival_analysis(
          interaction,
          df,
          covariates = covariates
          )
  }, error = function(e) {
        message("Error in survival analysis for interaction ", interaction, ": ", e$message)
        return(NULL)
  })    
  
  # If models failed, fill with NA
  if (is.null(model_list) || 
      is.null(model_list$model_full) || 
      is.null(model_list$model_reduced)) {

    interaction_row <- interaction_row %>%
      mutate(
        lrt_pval = NA,
        aic_lower = NA,
        wald_pval = NA
      )

    cov_results <- bind_rows(cov_results, interaction_row)
    setTxtProgressBar(pb, i)
    next
  }

  # --- LRT ---
  lrt <- anova(model_list$model_reduced, model_list$model_full, test = "Chisq")
  pval_lrt <- lrt[2, 'Pr(>|Chi|)']

  # --- AIC ---
  aic_full <- AIC(model_list$model_full)
  aic_reduced <- AIC(model_list$model_reduced)
  aic_lower <- aic_full < aic_reduced   # TRUE / FALSE

  # --- Wald ---
  coef_summary <- summary(model_list$model_full)$coefficients
  if ("group1" %in% rownames(coef_summary)) {
    pval_wald <- coef_summary["group1", "Pr(>|z|)"]
  } else {
    pval_wald <- NA
  }

  # Add results to row
  interaction_row <- interaction_row %>%
    mutate(
      lrt_pval = pval_lrt,
      aic_lower = aic_lower,
      wald_pval = pval_wald,
      covariates_used = paste(model_list$covariates_used, collapse = ";"),
    )
  
  # Concatenate interactions
  cov_results <- bind_rows(cov_results, interaction_row)

  setTxtProgressBar(pb, i)
}
close(pb)

# Join results with original dataframe
covariates_df <- left_join(covariates_df, cov_results, by = "interaction")

# Save results to the same file
cat('Saving results to:', matching_file, '\n')
write_csv(covariates_df, matching_file)

cat('Done: covariates_survival.R')
