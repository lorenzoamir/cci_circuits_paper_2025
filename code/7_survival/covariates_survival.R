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

covariates_dir = '/home/lnemati/pathway_crosstalk/results/survival/covariates/'

# List
covariates_cols = c("gender", "age_at_diagnosis", "binary_tumor_stage")

# Read the survival data
cat('Reading survival data\n')
df <- read_csv(tissuefile, col_types = cols())

cat('Number of rows in the dataset:', nrow(df), '\n')

# Fix column name if needed
if (names(df)[1] == "") {
    names(df)[1] <- "sample"
}

# Set first column as row names
df <- df %>% column_to_rownames(var = names(df)[1])

# Create patient id if missing
if (!"patient" %in% colnames(df)) {
    df$patient <- sapply(strsplit(rownames(df), '-'),
                         function(x) paste(x[1:4], collapse = '-'))
}

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
motifs <- covariates_df %>% pull(interaction)

# Convert gender to binary
binary_gender <- function(gender) {
    if (gender == "Male") return(0)
    if (gender == "Female") return(1)
    return(NA)
}
df$gender <- sapply(df$gender, binary_gender)

cat('Number of motifs:', length(motifs), '\n')

ccc <- unique(unlist(str_split(motifs, '&')))
cat('Number of ccc interactions:', length(ccc), '\n')


survival_analysis <- function(interaction, df, cov) {

    cat(interaction, '-', cov, '\n')

    genes <- unique(unlist(str_split(interaction, '[+&_]')))
    multiple <- length(unique(df$condition)) > 1

    # If covariate is age, dichotomize it at the median, do it condition-wise if multiple conditions
    # Don't overwrite the original covariate column, create a new binary column for the model
    if (cov == "age_at_diagnosis") {
        if (multiple) {
            df <- df %>%
                group_by(condition) %>%
                mutate(age_binary = age_at_diagnosis > median(age_at_diagnosis, na.rm = TRUE)) %>%
                ungroup()
        } else {
            df <- df %>%
                mutate(age_binary = age_at_diagnosis > median(age_at_diagnosis, na.rm = TRUE))
        }
        cov <- "age_binary"
    }

    cols <- c("patient","OS.time","OS","condition",genes,cov)
    if (!multiple) cols <- setdiff(cols,"condition")

    df <- df %>% select(all_of(cols))
    df <- df %>% mutate(OS.time = OS.time / 365)

    if (!cov %in% colnames(df)) {
        cat("Skipping", cov, "- column missing\n")
        return(NULL)
    }

    if (all(is.na(df[[cov]]))) {
        cat("Skipping", cov, "- all NA\n")
        return(NULL)
    }

    high_expression_group <- data.frame()
    low_expression_group  <- data.frame()

    if (multiple) {

        for (condition in unique(df$condition)) {

            tissue_df <- df %>% filter(.data$condition == condition)

            medians <- sapply(genes, function(g)
                median(tissue_df[[g]], na.rm = TRUE)
            )

            above_flags <- sapply(genes, function(g)
                tissue_df[[g]] > medians[g]
            )

            below_flags <- sapply(genes, function(g)
                tissue_df[[g]] <= medians[g]
            )

            above_all <- rowSums(above_flags) == length(genes)
            below_all <- rowSums(below_flags) == length(genes)

            high_expression_group <- bind_rows(
                high_expression_group,
                tissue_df[above_all, ]
            )

            low_expression_group <- bind_rows(
                low_expression_group,
                tissue_df[below_all, ]
            )
        }

    } else {

        medians <- sapply(genes, function(g)
            median(df[[g]], na.rm = TRUE)
        )

        above_flags <- sapply(genes, function(g)
            df[[g]] > medians[g]
        )

        below_flags <- sapply(genes, function(g)
            df[[g]] <= medians[g]
        )

        above_all <- rowSums(above_flags) == length(genes)
        below_all <- rowSums(below_flags) == length(genes)

        high_expression_group <- df[above_all, ]
        low_expression_group  <- df[below_all, ]
    }

    # Remove overlap
    patients_in_both <- intersect(
        high_expression_group$patient,
        low_expression_group$patient
    )

    high_expression_group <- high_expression_group %>%
        filter(!patient %in% patients_in_both)

    low_expression_group <- low_expression_group %>%
        filter(!patient %in% patients_in_both)

    n_high <- nrow(high_expression_group)
    n_low  <- nrow(low_expression_group)

    if (n_high < MIN_PATIENTS || n_low < MIN_PATIENTS) {
        cat("Skipping — too few patients\n")
        return(NULL)
    }

    # Combine only selected patients
    df <- bind_rows(high_expression_group, low_expression_group)

    df$group <- as.factor(df[[cov]])

    if (nlevels(df$group) < 2) {
        cat("Skipping — only one level\n")
        return(NULL)
    }

    if (multiple) {
        df$condition <- as.factor(df$condition)
        model <- coxph(Surv(OS.time, OS) ~ group + strata(condition), data=df)
    } else {
        model <- coxph(Surv(OS.time, OS) ~ group, data=df)
    }

    coef_name <- rownames(summary(model)$coefficients)[1]

    hr  <- exp(coef(model)[coef_name])
    se  <- summary(model)$coefficients[coef_name,"se(coef)"]
    # Propagate error:
    # hr = exp(coef) => err_hr = exp(coef) * se(coef) = hr * se(coef)
    se <- hr * se

    result <- list(
        hr = hr,
        n_patients_low = n_low,
        n_patients_high = n_high,
        logrank_pval = summary(model)$sctest["pvalue"],
        concordance_index = summary(model)$concordance[1],
        se = se
    )

    names(result) <- paste0(names(result), "_", cov)

    return(result)
}


original_cols <- c("interaction", "tissue", "hr", "concordance_index", "logrank_pval", "motif")

cat('Subsetting columns in the original dataframe\n')
covariates_df <- covariates_df %>% select(all_of(original_cols))

cat('Initializing covariate results dataframe\n')
cov_results <- data.frame()

for (interaction in motifs) {

    cat("\nProcessing interaction:", interaction, "\n")

    interaction_row <- data.frame(interaction = interaction)

    for (cov in covariates_cols) {

        cat("   covariate:", cov, "\n")

        result <- survival_analysis(interaction, df, cov)
        
        # Handle case where result is NULL by filling with NA values
        if (is.null(result)) {
            # Print a warning message for debugging
            cat("   Warning: result is NULL for interaction:", interaction, "and covariate:", cov, "\n")
            result <- list(
                hr = NA,
                n_patients_low = NA,
                n_patients_high = NA,
                logrank_pval = NA,
                concordance_index = NA,
                se = NA
            )
            # If cov is age_at_diagnosis, change the suffix to age_binary in the column names
            if (cov == "age_at_diagnosis") {
                names(result) <- paste0(names(result), "_age_binary")
            } else {
            names(result) <- paste0(names(result), "_", cov)
            }
        } 

        result_df <- as.data.frame(result)

        interaction_row <- bind_cols(interaction_row, result_df)
    }

    cov_results <- bind_rows(cov_results, interaction_row)
}

cat('Joining covariate results with original motifs dataframe\n')
covariates_df <- left_join(covariates_df, cov_results, by = "interaction")

print(head(covariates_df))
print(colnames(covariates_df))

cat('Saving results to:', matching_file, '\n')
write_csv(covariates_df, matching_file)

cat('Done: covariates_survival.R')

