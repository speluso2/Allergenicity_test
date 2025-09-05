install.packages(
  c("ggplot2", "dplyr"),
  repos = "https://cloud.r-project.org"
)

# Example: install Bioconductor package if needed:
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")
