# ==============================================================================
# 1. Load required packages
# ==============================================================================
if (!require("readxl")) install.packages("readxl")
if (!require("dplyr")) install.packages("dplyr")
if (!require("coin")) install.packages("coin")
if (!require("ggplot2")) install.packages("ggplot2")

library(readxl)
library(dplyr)
library(coin)
library(ggplot2)

# ==============================================================================
# 2. Define maps and functions
# ==============================================================================
aa_mapping <- c(
  'ALA'='A', 'ARG'='R', 'ASN'='N', 'ASP'='D', 'CYS'='C',
  'GLU'='E', 'GLN'='Q', 'GLY'='G', 'HIS'='H', 'ILE'='I',
  'LEU'='L', 'LYS'='K', 'MET'='M', 'PHE'='F', 'PRO'='P',
  'SER'='S', 'THR'='T', 'TRP'='W', 'TYR'='Y', 'VAL'='V'
)

assign_region <- function(rel_pos) {
  if (rel_pos <= 0.33) return('N-Terminal')
  else if (rel_pos <= 0.66) return('Core/Middle')
  else return('C-Terminal')
}

# ==============================================================================
# 3. Read and process Excel files
# ==============================================================================
# Default to current directory or use command-line argument if provided
args <- commandArgs(trailingOnly = TRUE)
root_path <- if (length(args) > 0) args[1] else "."

cat("Searching for project folders in:", root_path, "
")
folders <- list.dirs(root_path, recursive = FALSE)
all_data_list <- list()

cat("Reading and extracting data from folders...
")

for (folder in folders) {
  file_path <- file.path(folder, "all_positions_summary.xlsx")
  if (file.exists(file_path)) {
    tryCatch({
      df <- read_excel(file_path)
      colnames(df) <- c('Mutation', 'col_1_complex', 'col_1_peptide',
                        'DeltaG_binding', 'DeltaDeltaG_binding',
                        'DeltaG_mutation', 'Position')
      
      df <- df %>% 
        filter(!grepl("WTref", Mutation, ignore.case = TRUE)) %>%
        mutate(
          DeltaDeltaG_binding = as.numeric(gsub("/", ".", as.character(DeltaDeltaG_binding))),
          Position = as.numeric(gsub("/", ".", as.character(Position)))
        ) %>%
        filter(!is.na(Position) & !is.na(DeltaDeltaG_binding))
      
      if (nrow(df) > 0) {
        df <- df %>%
          mutate(
            Introduced_AA_3L = toupper(sapply(strsplit(as.character(Mutation), "_"), `[`, 1)),
            Introduced_AA = aa_mapping[Introduced_AA_3L],
            Relative_Position = Position / max(Position, na.rm = TRUE),
            Structural_Region = sapply(Relative_Position, assign_region),
            
            # Categorize amino acids by chemical group
            AA_Group = case_when(
              Introduced_AA %in% c('W', 'H', 'Y', 'F')      ~ 'Aromatic',
              Introduced_AA %in% c('K', 'R')                ~ 'Positively Charged',
              Introduced_AA %in% c('D', 'E')                ~ 'Negatively Charged',
              Introduced_AA %in% c('S', 'T', 'N', 'Q', 'C') ~ 'Polar Uncharged',
              TRUE                                          ~ 'Aliphatic/Nonpolar'
            ),
            
            Binding_Improved = ifelse(DeltaDeltaG_binding < 0, "Beneficial", "Destructive"),
            Peptide = basename(folder)
          )
        
        all_data_list[[length(all_data_list) + 1]] <- df
      }
    }, error = function(e) {})
  }
}

if (length(all_data_list) == 0) {
  stop("No valid 'all_positions_summary.xlsx' files found in subdirectories.")
}

# Merge datasets
full_df <- bind_rows(all_data_list)
cat("Processing complete. Total input mutations:", nrow(full_df), "

")

# ==============================================================================
# 4. Statistical Analysis 1: Fisher's Exact Test & Chi-Square Test
# ==============================================================================
cat("=======================================================\n")
cat("   1. Fisher's Exact Test Results for All Groups       \n")
cat("=======================================================\n")

table_all_groups <- with(full_df, table(AA_Group, Binding_Improved))
cat("\n[Contingency table of mutations by group]:\n")
print(table_all_groups)

set.seed(123)
fisher_overall <- fisher.test(table_all_groups, simulate.p.value = TRUE, B = 100000)

cat("\n[Simulated Fisher's Exact Test Output (B = 100,000)]:\n")
cat("Overall P-value across groups:", fisher_overall$p.value, "\n")

chisq_res <- chisq.test(table_all_groups)
cat("\n[Chi-Square Test Output]:\n")
cat("Chi-Square Statistic (X-squared):", chisq_res$statistic, "\n")
cat("Chi-Square P-value:", chisq_res$p.value, "\n")

if (fisher_overall$p.value < 0.05) {
  cat("\nResult: Significant difference in binding success rate among AA groups.\n")
}

# ==============================================================================
# 5. Statistical Analysis 2: Permutation Test
# ==============================================================================
cat("\n=======================================================\n")
cat("   2. Permutation Test Results for All Groups          \n")
cat("=======================================================\n")

set.seed(123)
perm_overall <- oneway_test(DeltaDeltaG_binding ~ factor(AA_Group), 
                            data = full_df, 
                            distribution = approximate(nresample = 10000))

cat("\n[Permutation Test Output (coin package)]:\n")
print(perm_overall)

p_val_perm <- pvalue(perm_overall)
cat("\nOverall Permutation P-value (10,000 resamples):", p_val_perm, "\n")

if (p_val_perm < 0.05) {
  cat("Result: Distribution of DeltaDeltaG values differs significantly between groups.\n")
}

# ==============================================================================
# 6. Pairwise Post-hoc Test (Bonferroni adjustment)
# ==============================================================================
cat("\n=======================================================\n")
cat("   3. Pairwise Group Comparisons (Pairwise Chi-Square) \n")
cat("=======================================================\n")

pairwise_res <- pairwise.prop.test(table_all_groups, p.adjust.method = "bonferroni")
print(pairwise_res)

# ==============================================================================
# 7. Generate plots
# ==============================================================================
x_labels <- c(
  "Aromatic"           = "Aromatic",
  "Positively Charged" = "Positively\nCharged",
  "Negatively Charged" = "Negatively\nCharged",
  "Aliphatic/Nonpolar" = "Aliphatic/\nNonpolar",
  "Polar Uncharged"    = "Polar\nUncharged"
)

df_summary <- full_df %>%
  group_by(AA_Group) %>%
  summarise(
    Total = n(),
    Beneficial = sum(Binding_Improved == "Beneficial"),
    Success_Rate = (Beneficial / Total) * 100
  )

p1 <- ggplot(df_summary, aes(x = reorder(AA_Group, -Success_Rate), y = Success_Rate, fill = AA_Group)) +
  geom_bar(stat = "identity", width = 0.55, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.1f%%", Success_Rate)), vjust = -0.5, size = 3.5) +
  labs(title = "Beneficial Mutation Success Rate by Amino Acid Group",
       x = "Amino Acid Group",
       y = "Success Rate (%)") +
  theme_minimal(base_size = 12) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

plot_path <- file.path(root_path, "beneficial_mutation_success_rate.png")
ggsave(plot_path, plot = p1, width = 8, height = 6, dpi = 300)
cat("\nPlot successfully saved to:", normalizePath(plot_path), "\n")
