library(dplyr)

# Constants
models <- c("ann", "rrblup")
trainwith_options <- c("F2", "F5", "F2_and_F5")
traingen_parents_combos <- expand.grid(
  traingen = c("F2", "F5"),
  parents = c("F2", "F5"),
  stringsAsFactors = FALSE
)
crosses <- c("C1", "C2", "C3")

excluded_gens <- c("ParentPoolC1", "ParentPoolC2", "ParentPoolC3")

for (model in models) {
  for (trainwith in trainwith_options) {
    for (i in 1:nrow(traingen_parents_combos)) {
      traingen <- traingen_parents_combos$traingen[i]
      parents <- traingen_parents_combos$parents[i]

      for (cross in crosses) {
        # Construct file name
        file_name <- sprintf(
          "%s_%s_random_trainAt%s_trainWith%s_%sParents_pheno_snp_yield.rds",
          cross, model, traingen, trainwith, parents
        )

        if (!file.exists(file_name)) {
          message("File not found: ", file_name)
          next
        }

        # Load data
        df <- readRDS(file_name)

        # Subset & convert to numeric (just once here)
        df = as.data.frame(df[1])
        df <- df[, c("gen", "pheno", "ebv")]
        df$pheno <- as.numeric(as.character(df$pheno))
        df$ebv <- as.numeric(as.character(df$ebv))

        # Exclude parent pool gens
        df <- df[!df$gen %in% excluded_gens, ]

        # Compute correlations by gen
        correlation_data <- data.frame()
        for (g in unique(df$gen)) {
          sub <- df[df$gen == g, ]
          corr <- cor(sub$pheno, sub$ebv)
          correlation_data <- rbind(correlation_data, data.frame(
            gen = g,
            correlation = corr
          ))
        }

        # Construct output filename
        output_filename <- sprintf(
          "%s_%s_trainAt%s_trainWith%s_%sPhenoEBV_correlations.csv",
          cross, model, traingen, trainwith, parents
        )

        # Save correlation data
        write.csv(correlation_data, output_filename, row.names = FALSE)
        message("Saved correlations to: ", output_filename)
      }
    }
  }
}
