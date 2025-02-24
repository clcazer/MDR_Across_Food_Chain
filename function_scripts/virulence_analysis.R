library(stringr)
library(rlist)
library(ggplot2)
library(tidyr)
library(matrixStats)
library(readxl)
library(cowplot)
library(arules)
library(dplyr)
library(reshape2)
library(gridExtra)
library(grid)
library(lubridate)
library(dplyr)
library(tidyr)
library(purrr)
library(multcomp)
library(MCMCglmm)
library(coda)  # For convergence diagnostics










get_virulence_gene_desciptives <- function(df, use_bayesian = FALSE) {
    print("column names")
    print(colnames(df))
    # Get unique data sources
    data_sources <- unique(df$meta)
    print("Data sources found:")
    print(data_sources)
    
    # Create list to store results
    results <- list()
    
    for (source in data_sources) {
        print(paste("Processing source:", source))
        
        # Subset data for this source
        source_df <- df[df$meta == source, ]
        
        # Get total number of unique isolates (using GCA column)
        n_isolates <- length(unique(source_df$GCA))
        print(paste("Number of isolates:", n_isolates))
        
        # Get unique genes for this source - using Element.symbol
        unique_genes <- unique(source_df$Element.symbol)
        n_unique_genes <- length(unique_genes)
        print(paste("Number of unique genes:", n_unique_genes))
        print("Unique genes found:")
        print(unique_genes)
        
        # Calculate gene counts and prevalence
        gene_prevalence <- data.frame(
            Gene = unique_genes
        )
        
        # Calculate isolate-level counts for each gene
        gene_prevalence$Isolate_Count <- sapply(gene_prevalence$Gene, function(g) {
            length(unique(source_df$GCA[source_df$Element.symbol == g]))
        })
        
        # Calculate prevalence percentage
        gene_prevalence$Prevalence <- round(gene_prevalence$Isolate_Count / n_isolates * 100, 2)
        
        # Sort by prevalence (descending)
        gene_prevalence <- gene_prevalence[order(-gene_prevalence$Prevalence), ]
        
        # Add summary rows
        gene_prevalence <- rbind(
            gene_prevalence,
            data.frame(
                Gene = "Total_Unique_Genes",
                Isolate_Count = n_unique_genes,
                Prevalence = NA
            ),
            data.frame(
                Gene = "Total_Isolates",
                Isolate_Count = n_isolates,
                Prevalence = NA
            )
        )
        
        # Save to CSV
        write.csv(
            gene_prevalence,
            file = paste0("Virulence_Plasmid_data/virulence_descriptives/virulence_prevalence_", source, ".csv"),
            row.names = FALSE
        )
        
        # Store results
        results[[source]] <- list(
            n_unique_genes = n_unique_genes,
            n_isolates = n_isolates,
            prevalence_df = gene_prevalence
        )
    }
    
    # Create dataframe to store genes per isolate statistics
    genes_per_isolate_stats <- data.frame(
        Data_Source = character(),
        Mean_Genes = numeric(),
        SD_Genes = numeric(),
        Variance = numeric(),
        stringsAsFactors = FALSE
    )
    
    # Calculate genes per isolate for each source
    isolate_gene_counts <- list()  # Store counts for statistical testing
    
    for (source in data_sources) {
        source_df <- df[df$meta == source, ]
        
        # Count genes per isolate
        genes_per_isolate <- table(source_df$GCA)
        
        # Store counts for statistical testing
        isolate_gene_counts[[source]] <- as.numeric(genes_per_isolate)
        
        # Calculate mean and SD
        genes_per_isolate_stats <- rbind(
            genes_per_isolate_stats,
            data.frame(
                Data_Source = source,
                Mean_Genes = mean(genes_per_isolate),
                SD_Genes = sd(genes_per_isolate),
                Variance = var(genes_per_isolate)
            )
        )
    }
    
    # Create data frame for regression
    all_counts <- data.frame(
        Count = unlist(isolate_gene_counts),
        Source = factor(rep(names(isolate_gene_counts), sapply(isolate_gene_counts, length)))  # Convert to factor
    )
    
    # Check for overdispersion
    mean_count <- mean(all_counts$Count)
    var_count <- var(all_counts$Count)
    dispersion_ratio <- var_count / mean_count
    
    if (use_bayesian) {
        require(MCMCglmm)
        require(coda)  # For convergence diagnostics
        require(ggplot2)
        
        # Prepare data
        all_counts <- data.frame(
            Count = unlist(isolate_gene_counts),
            Source = factor(rep(names(isolate_gene_counts), 
                              sapply(isolate_gene_counts, length)))
        )
        
        # Check for overdispersion
        dispersion_ratio <- var(all_counts$Count) / mean(all_counts$Count)
        
        # Set weakly informative priors - now including prior for random effects
        prior <- list(
            R = list(V = 1, nu = 0.002),  # Residual variance prior
            G = list(G1 = list(V = 1, nu = 0.002))  # Random effects prior
        )
        
        # Fit Bayesian model
        model <- if(dispersion_ratio > 1.5) {
            MCMCglmm(
                Count ~ Source,
                random = ~units,  # adds observation-level random effects for overdispersion
                family = "poisson",
                data = all_counts,
                prior = prior,
                nitt = 50000,
                burnin = 10000,
                thin = 5
            )
        } else {
            MCMCglmm(
                Count ~ Source,
                family = "poisson",
                data = all_counts,
                prior = prior,
                nitt = 50000,
                burnin = 10000,
                thin = 5
            )
        }
        
        test_type <- if(dispersion_ratio > 1.5) "Bayesian Poisson with OLRE" else "Bayesian Poisson"
        
        # Extract posterior summaries
        posterior_summary <- summary(model)
        
        # Create statistical results dataframe for Bayesian analysis
        # Get the first non-intercept term's pMCMC value (representing the first source comparison)
        first_source_pMCMC <- posterior_summary$solutions[2, "pMCMC"]  # Skip intercept, take first source
        
        # Convergence diagnostics
        mcmc_chain <- as.mcmc(model$Sol)
        geweke_results <- geweke.diag(mcmc_chain)
        effective_sizes <- effectiveSize(mcmc_chain)
        
        # Save convergence diagnostics
        convergence_results <- data.frame(
            Parameter = colnames(model$Sol),
            Geweke_Z = geweke_results$z,
            Effective_Size = effective_sizes,
            Total_Samples = nrow(model$Sol)
        )
        
        write.csv(convergence_results,
                 "Virulence_Plasmid_data/virulence_descriptives/virulence_genes_bayesian_convergence.csv",
                 row.names = FALSE)
        
        # Create statistical results dataframe with added ESS
        statistical_results <- data.frame(
            Test_Type = test_type,
            P_Value = first_source_pMCMC,
            Effect_Size = NA,
            Effect_Size_Type = "NA - Using Posterior Probabilities",
            Dispersion_Ratio = dispersion_ratio,
            Model_Used = "Bayesian Poisson",
            Min_Effective_Size = min(effective_sizes),
            Mean_Effective_Size = mean(effective_sizes)
        )

        write.csv(statistical_results,
            "Virulence_Plasmid_data/virulence_descriptives/virulence_genes_bayesian_main_test.csv",
            row.names = FALSE)
        
        # Source comparisons with visualizations
        source_comparisons <- NULL
        if (length(levels(all_counts$Source)) > 1) {
            pairs <- combn(levels(all_counts$Source), 2)
            prob_diff <- matrix(NA, ncol(pairs), 3)
            
            # Create directory for plots if it doesn't exist
            dir.create("Virulence_Plasmid_data/virulence_descriptives/posterior_plots", 
                      showWarnings = FALSE, recursive = TRUE)
            
            for (i in 1:ncol(pairs)) {
                # Extract relevant posterior samples
                source1_coef <- if(pairs[1,i] == levels(all_counts$Source)[1]) {
                    "(Intercept)"
                } else {
                    paste0("Source", pairs[1,i])
                }
                source2_coef <- if(pairs[2,i] == levels(all_counts$Source)[1]) {
                    "(Intercept)"
                } else {
                    paste0("Source", pairs[2,i])
                }
                
                # Calculate differences
                if (source1_coef == "(Intercept)") {
                    diff_samples <- -model$Sol[,source2_coef]
                } else if (source2_coef == "(Intercept)") {
                    diff_samples <- model$Sol[,source1_coef]
                } else {
                    diff_samples <- model$Sol[,source1_coef] - model$Sol[,source2_coef]
                }
                
                # Create density plot
                png(sprintf("Virulence_Plasmid_data/virulence_descriptives/posterior_plots/comparison_%s_vs_%s.png", 
                    pairs[1,i], pairs[2,i]), width = 800, height = 600)
                d <- density(diff_samples)
                # Adjust x-limits to include 0
                xlim <- range(c(0, d$x))
                plot(d,
                     main=sprintf("Posterior Distribution: %s vs %s", pairs[1,i], pairs[2,i]),
                     xlab="Difference in gene counts",
                     ylab="Density",
                     ylim=c(0, max(d$y)),
                     xlim=xlim)  # Set x-axis to include 0
                abline(v=0, lty=2, col="red")
                abline(v=quantile(diff_samples, 0.025), lty=3, col="blue")
                abline(v=quantile(diff_samples, 0.975), lty=3, col="blue")
                legend("topright", 
                      legend=c("Zero difference", "95% credible interval"),
                      lty=c(2,3), 
                      col=c("red", "blue"))
                dev.off()
                
                # Create trace plot (without density)
                png(sprintf("Virulence_Plasmid_data/virulence_descriptives/posterior_plots/trace_%s_vs_%s.png", 
                    pairs[1,i], pairs[2,i]), width = 800, height = 400)
                # Reset plotting parameters
                par(mfrow=c(1,1), mar=c(5,4,4,2)+0.1)
                plot(seq_along(diff_samples), diff_samples, type="l",
                     main=sprintf("Trace Plot: %s vs %s", pairs[1,i], pairs[2,i]),
                     xlab="Iteration",
                     ylab="Difference in gene counts")
                dev.off()
                
                prob_diff[i,] <- c(
                    mean(diff_samples > 0),
                    quantile(diff_samples, 0.025),
                    quantile(diff_samples, 0.975)
                )
            }
            
            source_comparisons <- data.frame(
                Comparison = paste(pairs[1,], "vs", pairs[2,]),
                Prob_Positive_Diff = prob_diff[,1],
                CI_Lower = prob_diff[,2],
                CI_Upper = prob_diff[,3],
                ESS = sapply(1:ncol(pairs), function(i) {
                    effectiveSize(as.mcmc(diff_samples))
                })
            )
        }
        
        # Save Bayesian results
        write.csv(
            data.frame(
                Parameter = rownames(posterior_summary$solutions),
                Mean = posterior_summary$solutions[,"post.mean"],
                Lower_CI = posterior_summary$solutions[,"l-95% CI"],
                Upper_CI = posterior_summary$solutions[,"u-95% CI"],
                pMCMC = posterior_summary$solutions[,"pMCMC"]
            ),
            "Virulence_Plasmid_data/virulence_descriptives/virulence_genes_bayesian_estimates.csv",
            row.names = FALSE
        )
        
        if (!is.null(source_comparisons)) {
            write.csv(
                source_comparisons,
                "Virulence_Plasmid_data/virulence_descriptives/virulence_genes_bayesian_comparisons.csv",
                row.names = FALSE
            )
        }
        
        # Add convergence results to the returned list
        results$convergence_diagnostics <- convergence_results
        results$trace_plots <- "Saved in Virulence_Plasmid_data/virulence_descriptives/posterior_plots/"
        
        # Add new results to returned list
        results$genes_per_isolate_stats <- genes_per_isolate_stats
        results$statistical_results <- statistical_results
        
    } else {
        # Determine model type based on dispersion
        if(dispersion_ratio > 1.5) {
            require(MASS)  # for negative binomial
            model <- glm.nb(Count ~ Source, data = all_counts)
            test_type <- "Negative Binomial Regression"
        } else {
            model <- glm(Count ~ Source, family = poisson, data = all_counts)
            test_type <- "Poisson Regression"
        }
        
        test_result <- anova(model, test = "Chisq")
        
        # Calculate effect size (pseudo-R²)
        null_model <- update(model, . ~ 1)
        effect_size <- 1 - exp((logLik(model) - logLik(null_model))/nrow(all_counts))
        
        # Perform pairwise comparisons using LRT
        source_levels <- levels(all_counts$Source)
        pairwise_comparisons <- data.frame(
            Comparison = character(),
            Chi_Square = numeric(),
            P_Value = numeric(),
            stringsAsFactors = FALSE
        )
        
        if (length(source_levels) > 1) {
            pairs <- combn(source_levels, 2)
            for (i in 1:ncol(pairs)) {
                tryCatch({
                    # Subset data for this pair
                    pair_data <- all_counts[all_counts$Source %in% pairs[,i],]
                    pair_data$Source <- droplevels(pair_data$Source)
                    
                    # Fit models
                    if(dispersion_ratio > 1.5) {
                        pair_model <- glm.nb(Count ~ Source, data = pair_data)
                        pair_null <- glm.nb(Count ~ 1, data = pair_data)
                    } else {
                        pair_model <- glm(Count ~ Source, family = poisson, data = pair_data)
                        pair_null <- glm(Count ~ 1, family = poisson, data = pair_data)
                    }
                    
                    # Likelihood ratio test
                    lrt <- anova(pair_model, pair_null)
                    
                    # Extract test statistics
                    chi_square <- abs(2 * (logLik(pair_model) - logLik(pair_null)))
                    p_value <- pchisq(chi_square, df = 1, lower.tail = FALSE)
                    
                    new_row <- data.frame(
                        Comparison = paste(pairs[1,i], "vs", pairs[2,i]),
                        Chi_Square = chi_square,
                        P_Value = p_value,
                        stringsAsFactors = FALSE
                    )
                    pairwise_comparisons <- rbind(pairwise_comparisons, new_row)
                    
                }, error = function(e) {
                    warning(paste("Error in pairwise comparison for pair", i, ":", e$message))
                })
            }
            
            # Only adjust p-values if we have any comparisons
            if (nrow(pairwise_comparisons) > 0) {
                pairwise_comparisons$Adjusted_P_Value <- p.adjust(pairwise_comparisons$P_Value, 
                                                                method = "bonferroni")
            }
        }
        
        # Create statistical results dataframe for frequentist analysis
        statistical_results <- data.frame(
            Test_Type = test_type,
            P_Value = test_result$`Pr(>Chi)`[2],  # P-value for Source term
            Effect_Size = effect_size,
            Effect_Size_Type = "McFadden_Pseudo_R2",
            Dispersion_Ratio = dispersion_ratio,
            Model_Used = test_type
        )
        
        # Save frequentist results with distinct file names
        write.csv(statistical_results,
                 "Virulence_Plasmid_data/virulence_descriptives/virulence_genes_frequentist_main_test.csv",
                 row.names = FALSE)
        
        if (nrow(pairwise_comparisons) > 0) {
            write.csv(pairwise_comparisons,
                     "Virulence_Plasmid_data/virulence_descriptives/virulence_genes_frequentist_pairwise.csv",
                     row.names = FALSE)
        }
        
        # Add pairwise comparisons to results
        results$pairwise_comparisons <- pairwise_comparisons
    }
    
    # Add dispersion information to stats
    genes_per_isolate_stats$Dispersion_Ratio <- dispersion_ratio
    
    # Save results
    write.csv(genes_per_isolate_stats, 
              "Virulence_Plasmid_data/virulence_descriptives/virulence_genes_per_isolate_stats.csv", 
              row.names = FALSE)
    
    return(results)
}

virulence_df <- read.csv("Virulence_Plasmid_data/amrfinder_results-VIRULENCE-PLUS.point-muts-included.csv")
results <- get_virulence_gene_desciptives(virulence_df, use_bayesian = TRUE)
# print(results)
print("DONE")