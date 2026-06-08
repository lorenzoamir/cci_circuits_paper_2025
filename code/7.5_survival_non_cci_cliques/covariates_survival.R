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

covariates_dir = '/home/lnemati/pathway_crosstalk/results/survival/non_cci_cliques/covariates/'

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

# Extract ccc interactions from circuits by splitting and flattening the list
ccc <- unique(unlist(str_split(circuits, '&')))  # Split interactions by '&' and flatten the list
cat('Number of ccc interactions:', length(ccc), '\n')

multiple <- length(unique(df$condition)) > 1

# Binarize gender only
if (multiple) {
  df <- df %>%
    group_by(condition) %>%
    mutate(gender = as.integer(gender == "Female")) %>%
    ungroup()
} else {
  df <- df %>%
    mutate(gender = as.integer(gender == "Female"))
}

survival_analysis <- function(interaction, df, covariates = NULL) {
  genes <- unique(unlist(str_split(interaction, '[+&_]')))

  cols <- c("patient", "OS.time", "OS", "condition", genes, covariates)
  if (!multiple) {
    cols <- setdiff(cols, "condition")
  }
  cols <- cols[cols %in% names(df)]
  df <- df %>% select(all_of(cols))

  df <- df %>% mutate(OS.time = OS.time / 365)

  if (multiple) df$condition <- as.factor(df$condition)

  # Filter valid covariates (must exist and not all NA)
  if (!is.null(covariates) && length(covariates) > 0) {
    valid_covariates <- covariates[
      covariates %in% names(df) &
      sapply(covariates, function(cov) !all(is.na(df[[cov]])))
    ]
    covariates <- if (length(valid_covariates) > 0) valid_covariates else NULL
  }

  gene_string <- paste(paste0("`", genes, "`"), collapse = " + ")
  cov_string <- if (!is.null(covariates) && length(covariates) > 0) paste(covariates, collapse = " + ") else NULL

  if (multiple) {
    full_rhs <- paste(c(gene_string, cov_string, "strata(condition)"), collapse = " + ")
    reduced_rhs <- paste(c(cov_string, "strata(condition)"), collapse = " + ")
  } else {
    full_rhs <- paste(c(gene_string, cov_string), collapse = " + ")
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
      covariates_used = covariates
    )
  }, warning = function(w) {
    message("Warning for motif ", interaction, ": ", w$message)
    result <- NULL
    if (!is.null(covariates) && length(covariates) > 0) {
      for (i in seq_along(covariates)) {
        covariates_subset <- covariates[1:(length(covariates) - i)]
        cov_string_subset <- if (length(covariates_subset) > 0) paste(covariates_subset, collapse = " + ") else NULL
        full_rhs_subset <- paste(c(gene_string, cov_string_subset, if (multiple) "strata(condition)"), collapse = " + ")
        reduced_rhs_subset <- paste(c(cov_string_subset, if (multiple) "strata(condition)"), collapse = " + ")
        formula_full_subset <- as.formula(paste("Surv(OS.time, OS) ~", full_rhs_subset))
        formula_reduced_subset <- if (!is.null(reduced_rhs_subset)) as.formula(paste("Surv(OS.time, OS) ~", reduced_rhs_subset)) else NULL
        result <- tryCatch({
          model_full_subset <- coxph(formula_full_subset, data = df)
          model_reduced_subset <- if (!is.null(formula_reduced_subset)) coxph(formula_reduced_subset, data = df) else NULL

          list(
            model_full = model_full_subset,
            model_reduced = model_reduced_subset,
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

# Set covariates; binary_tumor_stage and gender are binary, rest are numerical
covariates <- c("binary_tumor_stage", "age_at_diagnosis", "gender", "Lymphocytes", "Neutrophils", "Eosinophils", "Mast_Cells", "Dendritic_Cells", "Macrophages", "tumor_purity")

# Sort covariates by absolute correlation with OS.time
covariate_effects <- data.frame(covariate = covariates, effect_size = NA_real_)
for (cov in covariates) {
  if (!cov %in% names(df) || all(is.na(df[[cov]]))) {
    covariate_effects <- covariate_effects %>% filter(covariate != cov)
    next
  }
  effect_size <- abs(cor(df[[cov]], df$OS.time, use = "complete.obs"))
  covariate_effects <- covariate_effects %>% mutate(effect_size = ifelse(covariate == cov, effect_size, effect_size))
}
covariate_effects <- covariate_effects %>% arrange(desc(effect_size))
covariates <- covariate_effects$covariate
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
        aic_lower = NA
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

  # Add results to row
  interaction_row <- interaction_row %>%
    mutate(
      lrt_pval = pval_lrt,
      aic_lower = aic_lower,
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
