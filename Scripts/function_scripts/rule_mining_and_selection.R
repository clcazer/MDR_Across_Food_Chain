

#load in necessary libraries
library(arules)
library(stringr)
library(rlist)
library(ggplot2)
library(tidyr)
library(matrixStats)
library(readxl)
library(cowplot)
library(arulesViz)
library(data.table)


library(here)

source(here("Scripts", "function_scripts", "data_wrangling.R"))



#function to mine rules and return quality measures
mine_associations <- function(df, year, title, all_measures = TRUE, use_measures = NULL, write_rules = FALSE, resistance_indicator, target) {
  #generate rules with the apriori algorithm
  rules <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, year = year, target = target)


  if (all_measures == TRUE) {
  measures <- interestMeasure(rules)
  } else {measures <- interestMeasure(x = rules, measure = use_measures)}

  measures <- do.call(data.frame, lapply(measures, 
                    function(value) { replace(value, is.infinite(value), NA)}))


  measures <- measures[, which(colMeans(is.na(measures)) < 0.1)] #keep only columns whare NA is less that %10
  measures <- measures[is.finite(rowSums(measures)), ] #drop remaining NA row-wise

    if (write_rules == TRUE) {
    write(rules,
        file = str_glue('ruleData/{title}_R_association_rules.csv'),
        sep = ",",
        quote = TRUE,
        row.names = FALSE)
        write.csv(measures, str_glue('ruleData/{title}_Rassociation_measures.csv'), row.names = FALSE)

    }

  return(measures)
}




#function to select quality measures and return the selected measures
select_quality_measures <- function(plothist = FALSE, resistance_indicator, df, target, data_source) {
    #create empty lists to hold the top loadings later
    pc1_top_loads <- list()
    pc2_top_loads <- list()
    pc3_top_loads <- list()
    pc4_top_loads <- list()

    #apply needed functions to the files (i.e., rule mining, pca, and selecting top loadings)
    for (year in c(min(df[, "Year"]):max(df[, "Year"]))) {
        # Check if there's any data for this year
        year_data <- df[df$Year == year, ]
        if (nrow(year_data) == 0) {
            warning(paste("Skipping year", year, "- no data available"))
            next
        }

        # Try to mine rules and handle potential errors
        result <- tryCatch({
            #mine rules and return quality measures as data variable
            rule_measures <- mine_associations(df = year_data, year = year, all_measures = TRUE, 
                                            resistance_indicator = resistance_indicator, target = target)

            # Remove constant columns
            non_constant_cols <- apply(rule_measures, 2, function(x) length(unique(x)) > 1)
            rule_measures <- rule_measures[, non_constant_cols, drop = FALSE]

            # Check if there are any non-constant columns left
            if (ncol(rule_measures) == 0) {
                return(NULL)
            }

            pca <- prcomp(rule_measures, scale = TRUE, center = TRUE, retx = TRUE, rank. = min(4, ncol(rule_measures)))

            #get the top loadings for each of the first 4 PCs
            load <- abs(pca$rotation)
            prop_load <- as.data.frame(apply(load, 2, function(x) x / sum(x)))

            pc1_loads <- prop_load |> sort_by(~ list(-PC1))
            pc2_loads <- prop_load |> sort_by(~ list(-PC2))
            pc3_loads <- prop_load |> sort_by(~ list(-PC3))
            pc4_loads <- prop_load |> sort_by(~ list(-PC4))

            list(
                pc1 = rownames(pc1_loads[1:5, ]),
                pc2 = rownames(pc2_loads[1:5, ]),
                pc3 = rownames(pc3_loads[1:5, ]),
                pc4 = rownames(pc4_loads[1:5, ])
            )
        }, error = function(e) {
            warning(paste("Error processing year", year, ":", e$message))
            return(NULL)
        })

        # If we got results, append them to our lists
        if (!is.null(result)) {
            pc1_top_loads <- list.append(pc1_top_loads, result$pc1)
            pc2_top_loads <- list.append(pc2_top_loads, result$pc2)
            pc3_top_loads <- list.append(pc3_top_loads, result$pc3)
            pc4_top_loads <- list.append(pc4_top_loads, result$pc4)
        }
    }

    # Check if we have any results
    if (length(pc1_top_loads) == 0) {
        stop("No valid data found for any year")
    }

    pca_df <- data.frame(PC1 = unlist(pc1_top_loads),
                        PC2 = unlist(pc2_top_loads),
                        PC3 = unlist(pc3_top_loads),
                        PC4 = unlist(pc4_top_loads))

    final_measures <- c(names(which.max(table(pca_df$PC1))),
                       names(which.max(table(pca_df$PC2))),
                       names(which.max(table(pca_df$PC3))),
                       names(which.max(table(pca_df$PC4))))

    if (plothist == TRUE) {
        pc1_plot <- ggplot(data = data.frame(pca_df$PC1), aes(x= pca_df$PC1)) + 
        geom_bar() + 
        labs(x = "quality measures", title = "PC1") +
        theme(axis.text.x = element_text(
            angle = 22.5,
            hjust = 1,
            size = 10
        ))
        ggsave(filename = str_glue("dataset_specific_outputs/{data_source}/figures/measure_pca_toploadings/{resistance_indicator}/{data_source}_{resistance_indicator}_pc1.png"), plot = pc1_plot)

        pc2_plot <- ggplot(data = data.frame(pca_df$PC2), aes(x= pca_df$PC2)) + 
        geom_bar() + 
        labs(x = "quality measures", title = "PC2") +
        theme( axis.text.x = element_text(
            angle = 22.5,
            hjust = 1,
            size = 10
          ))
        ggsave(filename = str_glue("dataset_specific_outputs/{data_source}/figures/measure_pca_toploadings/{resistance_indicator}/{data_source}_{resistance_indicator}_pc2.png"), plot = pc2_plot)
        

        pc3_plot <- ggplot(data = data.frame(pca_df$PC3), aes(x= pca_df$PC3)) + 
        geom_bar() + 
        labs(x = "quality measures", title = "PC3") +
        theme( axis.text.x = element_text(
            angle = 22.5,
            hjust = 1,
            size = 10
          ))

        ggsave(filename = str_glue("dataset_specific_outputs/{data_source}/figures/measure_pca_toploadings/{resistance_indicator}/{data_source}_{resistance_indicator}_pc3.png"), plot = pc3_plot)

        pc4_plot <- ggplot(data = data.frame(pca_df$PC4), aes(x= pca_df$PC4)) + 
        geom_bar() + 
        labs(x = "quality measures", title = "PC4") +
        theme( axis.text.x = element_text(
            angle = 22.5,
            hjust = 1,
            size = 10
          ))
        ggsave(filename = str_glue("dataset_specific_outputs/{data_source}/figures/measure_pca_toploadings/{resistance_indicator}/{data_source}_{resistance_indicator}_pc4.png"), plot = pc4_plot)
    }
    
    return(final_measures)
}



#function to mine rules and produce
#histrograms based on quality measure distributions
#for the purpose of selecting cutoff points
select_measure_cutoffs <- function(selected_measures, df, resistance_indicator, target, data_source) {
  #loop through each file (i.e, year)
  #get get a set of rules for each file
  #plot histogram for each of the files' 4 quality measures

  for (year in c(min(df[, "Year"]):max(df[, "Year"]))) {
    # Check if there's any data for this year
    year_data <- df[df$Year == year, ]
    if (nrow(year_data) == 0) {
      warning(paste("Skipping year", year, "- no data available"))
      next
    }

    print(paste("Processing year:", year))

    # Try to mine rules and handle potential errors
    result <- tryCatch({
      rule_measures <- mine_associations(df = year_data, year = year, 
                                       all_measures = FALSE, 
                                       use_measures = selected_measures,
                                       resistance_indicator = resistance_indicator, 
                                       target = target)

      if (is.null(rule_measures) || nrow(rule_measures) == 0) {
        warning(paste("No rules generated for year", year))
        return(NULL)
      }

      sample_size <- length(unique(year_data$ID))
      
      rule_measures <- as.data.frame(rule_measures)
      
      plot <- ggplot(gather(rule_measures), aes(value)) + 
        geom_histogram(bins = 50) + 
        facet_wrap(~key, scales = 'free_x') +
        labs(y = "number of rules", 
             title = str_glue('year = {year}; total rules = {nrow(rule_measures)}; sample size = {sample_size}'))
      
      ggsave(plot,
             filename = str_glue("dataset_specific_outputs/{data_source}/figures/r_measure_dists/{resistance_indicator}/{data_source}_{year}_{resistance_indicator}_measure_distributions.png"))
      
      TRUE  # Return TRUE to indicate successful processing
    }, error = function(e) {
      warning(paste("Error processing year", year, ":", e$message))
      return(NULL)
    })
    
    if (is.null(result)) {
      warning(paste("Failed to process year", year))
    }
  }
}



#create a function to make rarefaction curves
#with sample size on the x axis and number of rules of the y axis
#will need to shuffle all the data and take a subsample then
# run the mine_associations function on the subsample
#store the subsample size and the number of rules (repeat for each subsample size)
#then just graph
rarefaction_by_sample_size <- function(df, resistance_indicator, target, data_source) {
  print("started")
  #create empty list to hold matrices of data
  lom <- list()
  
  for (year in c(min(df[, "Year"]):max(df[, "Year"]))) {
    #read in the data
    year_df <- df[df$Year == year, ]
    
    # Skip years with no data
    if (nrow(year_df) == 0) {
      print(paste("Skipping year", year, "- no data available"))
      next
    }
    
    # Skip years with insufficient data
    if (nrow(year_df) < 20) {
      print(paste("Skipping year", year, "- insufficient data (less than 20 samples)"))
      next
    }

    #create matrix to hold data
    data_matrix <- matrix(0, 10, 10) #50 subsamples repeated and averaged 100 times

    #get evenly spaced subsample sizes
    subsample_size <- floor(seq(from = 20, to = nrow(year_df), length.out = 10))





    #also define some variables to index and update this matrix inside the inner forloop
    data_col <- 1
    while (data_col <= ncol(data_matrix)) {
    data_row <- 1
  print(year)
    #loop through all the subample sizes
    for (size in subsample_size){
    
    #get a random subsample for each subsample size in the sequence
    sub_sample <- year_df[sample(nrow(year_df), size), ]
    #do all the stuff to run apriori
    rules <- get_rules_or_itemsets(df = sub_sample, resistance_indicator = resistance_indicator,
                                   year = year, target = target)
  print(size)
  print(length(rules))
    data_matrix[data_row, data_col] <- length(rules)


    

    data_row <- data_row + 1
    }
    data_col <- data_col + 1
    }

    data_mean <- rowMeans2(data_matrix)
    data_stnd_err <- (rowSds(data_matrix))/(sqrt(ncol(data_matrix)))
    data_matrix <- cbind(data_matrix, data_mean)
    data_matrix <- cbind(data_matrix, data_stnd_err)
    data_matrix <- cbind(data_matrix, subsample_size)
    year <- c(replicate(10, as.numeric(year)))
    data_matrix <- cbind(data_matrix, year)
    data_matrix <- subset(data_matrix, select = c(data_mean, data_stnd_err, subsample_size, year))
    

    lom <- list.append(lom, data_matrix)
  }

  lom_df <- do.call(rbind, lapply(lom, data.frame))
  #lom_df <- lom_df[-which(lom_df$year == 2006),]
  #lom_df <- lom_df[-which(lom_df$year == 2013),]

  #mean_df <- gather(lom_df, key = "year", value = "data_mean")
  #(lom_df)
  lom_df$year <- as.factor(lom_df$year)
  plot <- ggplot(lom_df, aes(x = subsample_size, y = data_mean, 
  ymin = (data_mean - data_stnd_err), ymax = (data_mean + data_stnd_err), 
  colour = year, fill=year, linetype=year)) +
    geom_line(linewidth = 2) +
    geom_point(size = 7) +
    theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = 'white', colour = 'black')) +
    labs(y = str_glue("Number of {target}"), x = "Sample Size", title = str_glue("Rarefaction Curve ({resistance_indicator})")) +
    geom_ribbon(alpha=0.2)
    ggsave(str_glue("dataset_specific_outputs/{data_source}/figures/rarefaction_curves/{resistance_indicator}/{data_source}_{resistance_indicator}_rarefaction_by_sample_size_{target}.png"))


}







#function to make a rarefaction curve with number of rules on the x axis and
#number of unique class based rules on the y axis
rarefaction_rule_v_class <- function(df, resistance_indicator, target, data_source) {
  #read in all the files
  if (resistance_indicator == "phenotype") {
    class_mappings <- read_excel("EcoliBreakPoints.xlsx")
  }
  else if (resistance_indicator == "genotype") {
    class_mappings <- read.csv("gene_class_mappings.csv")
  } 
  else {stop("ERROR: resistance_indicator must be either phenotype or genotype")}

  #create empty list to hold matrices of data
  lom <- list()

  for (year in c(min(df[, "Year"]):max(df[, "Year"]))) {
    # Check if there's any data for this year
    year_data <- df[df$Year == year, ]
    if (nrow(year_data) == 0) {
      warning(paste("Skipping year", year, "- no data available"))
      next
    }

    # Try to get rules and handle potential errors
    result <- tryCatch({
      #do all the stuff to run apriori
      rules <- get_rules_or_itemsets(df = year_data, resistance_indicator = resistance_indicator,
                                    year = year, target = target)
      
      if (length(rules) == 0) {
        warning(paste("No rules generated for year", year))
        return(NULL)
      }
      
      pattern <- str_extract_all(labels(rules), "\\S+")

      if (resistance_indicator == "phenotype") {
        for (i in seq(class_mappings$abbreviation)) {
          pattern <- gsub(class_mappings$abbreviation[i], class_mappings$class[i], pattern, fixed = TRUE)
        }
      } else if (resistance_indicator == "genotype") {
        for (i in seq(class_mappings$Gene_family)) {
          pattern <- gsub(class_mappings$Gene_family[i], class_mappings$corrected_class[i], pattern, fixed = TRUE)
        }
      }

      #create matrix to hold data
      data_matrix <- matrix(0, 10, 10) #10 subsamples repeated and averaged 10 times

      #get 10 different evenly spaced subsample sizes
      subsample_size <- floor(seq(from = 5, to = length(pattern), length.out = 10))

      #also define some variables to index and update this matrix inside the inner forloop
      data_col <- 1
      while (data_col <= ncol(data_matrix)) {
        data_row <- 1
        #loop through all the subsample sizes
        for (size in subsample_size) {
          sub_sample <- pattern[sample(length(pattern), size)]
          data_matrix[data_row, data_col] <- length(unique(sub_sample))
          data_row <- data_row + 1
        }
        data_col <- data_col + 1
      }

      data_mean <- rowMeans2(data_matrix)
      data_stnd_err <- (rowSds(data_matrix))/(sqrt(ncol(data_matrix)))
      data_matrix <- cbind(data_matrix, data_mean)
      data_matrix <- cbind(data_matrix, data_stnd_err)
      data_matrix <- cbind(data_matrix, subsample_size)
      year_col <- c(replicate(10, as.numeric(year)))
      data_matrix <- cbind(data_matrix, year_col)
      
      data_matrix <- subset(data_matrix, select = c(data_mean, data_stnd_err, subsample_size, year_col))
      data_matrix
      
    }, error = function(e) {
      warning(paste("Error processing year", year, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(result)) {
      lom <- list.append(lom, result)
    }
  }

  # Check if we have any results
  if (length(lom) == 0) {
    stop("No valid data found for any year")
  }

  lom_df <- do.call(rbind, lapply(lom, data.frame))
  names(lom_df)[names(lom_df) == "year_col"] <- "year"  # Fix column name
  
  lom_df$year <- as.factor(lom_df$year)
  plot <- ggplot(lom_df, aes(x = subsample_size, y = data_mean, 
                            ymin = (data_mean - data_stnd_err), ymax = (data_mean + data_stnd_err), 
                            colour = year, fill=year, linetype=year)) +
    geom_line(linewidth = 2) +
    geom_point(size = 7) +
    theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = 'white', colour = 'black')) +
    labs(y = str_glue("Drug Class {target}"), 
         x = str_glue("Antimicrobial {target}"), 
         title = str_glue("Rarefaction Curve ({resistance_indicator})")) +
    geom_ribbon(alpha=0.2)
  
  ggsave(filename = str_glue("dataset_specific_outputs/{data_source}/figures/rarefaction_curves/{resistance_indicator}/{data_source}_{resistance_indicator}_rarefaction_by_class_{target}.png"))
}






unique_classes_represented <- function(df, resistance_indicator, target, data_source) {
  if (resistance_indicator == "phenotype") {
    class_mappings <- read_excel("EcoliBreakPoints.xlsx")
  }
  else if (resistance_indicator == "genotype") {
    class_mappings <- read.csv("gene_class_mappings.csv")
  } 
  else {stop("ERROR: resistance_indicator must be either phenotype or genotype")}

  #create empty list to hold matrices of data
  lom <- list()

  for (year in c(min(df[, "Year"]):max(df[, "Year"]))) {
    # Check if there's any data for this year
    year_data <- df[df$Year == year, ]
    if (nrow(year_data) == 0) {
      warning(paste("Skipping year", year, "- no data available"))
      next
    }

    # Try to process the year and handle potential errors
    result <- tryCatch({
      #create matrix to hold data
      data_matrix <- matrix(0, 10, 10) #10 subsamples repeated and averaged 10 times

      #do all the stuff to run apriori
      rules <- get_rules_or_itemsets(df = year_data, resistance_indicator = resistance_indicator,
                                    year = year, target = target)
      
      if (length(rules) == 0) {
        warning(paste("No rules generated for year", year))
        return(NULL)
      }
      
      pattern <- str_extract(labels(rules), "\\S+")

      if (resistance_indicator == "phenotype") {
        for (i in seq(class_mappings$abbreviation)) {
          pattern <- gsub(class_mappings$abbreviation[i], class_mappings$class[i], pattern, fixed = TRUE)
        }
      }
      else if (resistance_indicator == "genotype") {
        for (i in seq(class_mappings$Gene_family)) {
          pattern <- gsub(class_mappings$Gene_family[i], class_mappings$corrected_class[i], pattern, fixed = TRUE)
        }
      }

      #get 10 different evenly spaced subsample sizes
      subsample_size <- floor(seq(from = 5, to = length(pattern), length.out = 10))

      #also define some variables to index and update this matrix inside the inner forloop
      data_col <- 1
      while (data_col <= ncol(data_matrix)) {
        data_row <- 1
        #loop through all the subsample sizes
        for (size in subsample_size) {
          sub_sample <- pattern[sample(length(pattern), size)]
          
          sub_sample <- unique(unlist(strsplit(unique(sub_sample), ",")))
          sub_sample <- unique(unlist(strsplit(unique(sub_sample), "\\{")))
          sub_sample <- unique(unlist(strsplit(unique(sub_sample), "\\}")))
          
          data_matrix[data_row, data_col] <- length(unique(sub_sample))
          data_row <- data_row + 1
        }
        data_col <- data_col + 1
      }

      data_mean <- rowMeans2(data_matrix)
      data_stnd_err <- (rowSds(data_matrix))/(sqrt(ncol(data_matrix)))
      data_matrix <- cbind(data_matrix, data_mean)
      data_matrix <- cbind(data_matrix, data_stnd_err)
      data_matrix <- cbind(data_matrix, subsample_size)
      year_col <- c(replicate(10, as.numeric(year)))
      data_matrix <- cbind(data_matrix, year_col)
      
      data_matrix <- subset(data_matrix, select = c(data_mean, data_stnd_err, subsample_size, year_col))
      data_matrix
      
    }, error = function(e) {
      warning(paste("Error processing year", year, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(result)) {
      lom <- list.append(lom, result)
    }
  }

  # Check if we have any results
  if (length(lom) == 0) {
    stop("No valid data found for any year")
  }

  lom_df <- do.call(rbind, lapply(lom, data.frame))
  names(lom_df)[names(lom_df) == "year_col"] <- "year"  # Fix column name
  
  lom_df$year <- as.factor(lom_df$year)
  plot <- ggplot(lom_df, aes(x = subsample_size, y = data_mean, 
                            ymin = (data_mean - data_stnd_err), ymax = (data_mean + data_stnd_err), 
                            colour = year, fill=year, linetype=year)) +
    geom_line(linewidth = 2) +
    geom_point(size = 7) +
    theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5),
          panel.background = element_rect(fill = 'white', colour = 'black')) +
    labs(y = str_glue("Unique Classes Represented"), 
         x = str_glue("Antimicrobial {target}"), 
         title = str_glue("")) +
    geom_ribbon(alpha=0.2) +
    geom_vline(xintercept = 10) +
    geom_vline(xintercept = 1000)
  
  ggsave(filename = str_glue("dataset_specific_outputs/{data_source}/figures/rarefaction_curves/{resistance_indicator}/{data_source}_{resistance_indicator}_unique_classes_{target}.png"))
}





#plot measure cutoff on the x axis and number of rules or number unique
#class based rules on the y axis
rules_v_cutoffs <- function(df, resistance_indicator, target, measures, low, high, data_source) {
  #read in class mappings
  if (resistance_indicator == "phenotype") {
    class_mappings <- read_excel("EcoliBreakPoints.xlsx")
  }
  else if (resistance_indicator == "genotype") {
    class_mappings <- read.csv("gene_class_mappings.csv")
  } 
  else {stop("ERROR: resistance_indicator must be either phenotype or genotype")}

  for (measure in measures) {
    #create empty list to hold matrices of data
    lom <- list()

    for (year in c(min(df[, "Year"]):max(df[, "Year"]))) {
      # Check if there's any data for this year
      year_data <- df[df$Year == year, ]
      if (nrow(year_data) == 0) {
        warning(paste("Skipping year", year, "- no data available"))
        next
      }

      # Try to process the year and handle potential errors
      result <- tryCatch({
        #create matrix to hold data
        data_matrix <- matrix(0, 10, 2) #10 subsamples repeated and averaged 10 times

        #do all the stuff to run apriori
        rules <- get_rules_or_itemsets(df = year_data, resistance_indicator = resistance_indicator,
                                     year = year, target = target)
        
        if (length(rules) == 0) {
          warning(paste("No rules generated for year", year))
          return(NULL)
        }

        measures_df <- interestMeasure(rules, measure = measure)
        cut_off_value <- seq(from = low, to = high, length.out = 10)
        
        data_row <- 1
        for (cut_off in cut_off_value) {
          new_rules <- rules
          delete_index <- c(which((measures_df < cut_off)))

          if (length(delete_index != 0)) {
            new_rules <- rules[-delete_index]
          }

          pattern <- str_extract_all(labels(new_rules), "\\S+")
          
          if (resistance_indicator == "phenotype") {
            for (i in seq(class_mappings$abbreviation)) {
              pattern <- gsub(class_mappings$abbreviation[i], class_mappings$class[i], pattern, fixed = TRUE)
            }
          } else if (resistance_indicator == "genotype") {
            for (i in seq(class_mappings$Gene_family)) {
              pattern <- gsub(class_mappings$Gene_family[i], class_mappings$corrected_class[i], pattern, fixed = TRUE)
            }
          }

          data_matrix[data_row, 1] <- length(unique(pattern))
          data_matrix[data_row, 2] <- length(new_rules)
          data_row <- data_row + 1
        }

        colnames(data_matrix) <- c("class", "drug")
        data_matrix <- cbind(data_matrix, cut_off_value)
        year_col <- c(replicate(10, as.numeric(year)))
        data_matrix <- cbind(data_matrix, year_col)
        
        data_matrix <- subset(data_matrix, select = c(class, drug, cut_off_value, year_col))
        data_matrix
        
      }, error = function(e) {
        warning(paste("Error processing year", year, ":", e$message))
        return(NULL)
      })
      
      if (!is.null(result)) {
        lom <- list.append(lom, result)
      }
    }

    # Check if we have any results
    if (length(lom) == 0) {
      warning(paste("No valid data found for measure:", measure))
      next
    }

    lom_df <- do.call(rbind, lapply(lom, data.frame))
    names(lom_df)[names(lom_df) == "year_col"] <- "year"  # Fix column name
    
    lom_df$year <- as.factor(lom_df$year)

    plot1 <- ggplot(lom_df, aes(x = cut_off_value, y = class, 
                               colour = year, fill=year)) +
      geom_line(linewidth = 2) +
      geom_point(size = 3) +
      theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black')) +
      labs(y = str_glue("Class {target}"), x = NULL, title = str_glue("{measure} ({resistance_indicator})"))

    plot2 <- ggplot(lom_df, aes(x = cut_off_value, y = drug, 
                               colour = year, fill=year)) +
      geom_line(linewidth = 2) +
      geom_point(size = 3) +
      geom_hline(yintercept = 1000, linetype = "dotted", color = "red", linewidth = 1) +
      geom_vline(xintercept = case_when(
        measure == "cosine" ~ 0.5,
        measure == "jaccard" ~ 0.0,
        measure == "kulczynski" ~ 0.5,
        measure == "support" ~ 0.0,
        TRUE ~ 0.0
      ), color = "red", linewidth = 1) +
      theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5),
            panel.background = element_rect(fill = 'white', colour = 'black')) +
      labs(y = str_glue("number of {target}"), x = "Cut-off Value", title = str_glue(""))

    legend <- get_legend(plot2)
    pgrid <- plot_grid(plot1 + theme(legend.position="none"), 
                      plot2 + theme(legend.position="none"), 
                      ncol = 1)
    p <- plot_grid(pgrid, legend, ncol = 2, rel_widths = c(1, .12))

    ggsave(plot = plot2, 
           filename = str_glue("dataset_specific_outputs/{data_source}/figures/cut_off_selection/{resistance_indicator}/{data_source}_{resistance_indicator}_{target}_v_cutoffs_{measure}.png"), 
           height = 4, width = 6)
  }
}





#function to trim the rules based on the chosen cut-off values for the selected measures for a single year
implement_cuts <- function(df, resistance_indicator, target, cut_off, year, agg = FALSE, measures_used) {

    #do all the stuff to run apriori
    rules <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator,
                                   year = year, target = target, agg = agg)



    measures_df <- interestMeasure(rules, measure = c(measures_used))
    
    rules@quality[, measures_used[1]] <- measures_df[, measures_used[1]]
    rules@quality[, measures_used[2]] <- measures_df[, measures_used[2]]
    rules@quality[, measures_used[3]] <- measures_df[, measures_used[3]]
    rules@quality[, measures_used[4]] <- measures_df[, measures_used[4]]



   rules <- arules::subset(rules, subset = cosine > cut_off[1] & jaccard > cut_off[2] & kulczynski > cut_off[3] & support > cut_off[4])


    #write(rules,
        #file = str_glue("test.csv"),
       # sep = ",",
       # quote = TRUE,
       # row.names = FALSE)



    return(rules)

}


save_item_or_ruleset_aggregated <- function(genotype_df, phenotype_df, target, data_source) {

df_list <- list()
print(c(min(genotype_df[, "Year"]), max(genotype_df[, "Year"])))
for (year in c(min(genotype_df[, "Year"]):(max(genotype_df[, "Year"])))) {

    phenotype_rules_or_itemsets <- get_rules_or_itemsets(df = phenotype_df, resistance_indicator = "phenotype",
                                    year = year, target = target, agg = TRUE)



    genotype_rules_or_itemsets <- get_rules_or_itemsets(df = genotype_df, resistance_indicator = "genotype",
                                    year = year, target = target, agg = FALSE)
    if (target == "rules") {
    measures_df <- interestMeasure(phenotype_rules_or_itemsets, measure = c("cosine", "jaccard", "kulczynski"))
    
    phenotype_rules_or_itemsets@quality[, "cosine"] <- measures_df[, "cosine"]
    phenotype_rules_or_itemsets@quality[, "jaccard"] <- measures_df[, "jaccard"]
    phenotype_rules_or_itemsets@quality[, "kulczynski"] <- measures_df[, "kulczynski"]

    measures_df <- interestMeasure(genotype_rules_or_itemsets, measure = c("cosine", "jaccard", "kulczynski"))
    
    genotype_rules_or_itemsets@quality[, "cosine"] <- measures_df[, "cosine"]
    genotype_rules_or_itemsets@quality[, "jaccard"] <- measures_df[, "jaccard"]
    genotype_rules_or_itemsets@quality[, "kulczynski"] <- measures_df[, "kulczynski"]
    }
    phenotype_rule_item_df <- DATAFRAME(phenotype_rules_or_itemsets)
    phenotype_rule_item_df$Year <- rep(c(year), times = nrow(phenotype_rule_item_df))
    phenotype_rule_item_df$Resistance_Indicator <- rep(c("Phenotype"), times = nrow(phenotype_rule_item_df))

    genotype_rule_item_df <- DATAFRAME(genotype_rules_or_itemsets)
    genotype_rule_item_df$Year <- rep(c(year), times = nrow(genotype_rule_item_df))
    genotype_rule_item_df$Resistance_Indicator <- rep(c("Genotype"), times = nrow(genotype_rule_item_df))

    df_list <- list.append(df_list, phenotype_rule_item_df)
    df_list <- list.append(df_list, genotype_rule_item_df)
}

full_df <- bind_rows(df_list)

write.csv(full_df, str_glue("dataset_specific_outputs/{data_source}/ruleData/{data_source}_genotypeAndPhenotype_{target}.csv"), row.names = FALSE)

}








save_top_gene_rules_aggregated_v2 <- function(Retail_Meats_df, cecal_df, NAHLN_df, target) {
  # Pre-compile regex patterns and create lookup tables outside the loops
  gene_class_mappings <- read.csv("gene_class_mappings.csv")
  sanitized_genes <- make.names(gene_class_mappings$Genes)
  gene_name_map <- setNames(gene_class_mappings$Genes, sanitized_genes)
  split_pattern <- ",|\\{|\\}"
  
  # Process each dataset
  datasets <- list(
    "Retail_Meats" = Retail_Meats_df,
    "Cecal" = cecal_df,
    "NAHLN" = NAHLN_df
  )
  
  # Vectorized function to process gene names
  process_genes_vec <- function(genes_vec) {
    vapply(genes_vec, function(genes_str) {
      if (!is.character(genes_str)) return(NA_character_)
      
      genes <- strsplit(genes_str, split_pattern)[[1]]
      genes <- trimws(genes[genes != ""])
      fixed <- vapply(genes, function(g) {
        if (g %in% names(gene_name_map)) gene_name_map[g] else g
      }, character(1))
      paste0("{", paste(fixed, collapse = ","), "}")
    }, character(1))
  }
  
  # Function to standardize combinations
  standardize_combination <- function(lhs, rhs) {
    genes <- unique(c(
      strsplit(gsub("[{}]", "", lhs), ",")[[1]],
      strsplit(gsub("[{}]", "", rhs), ",")[[1]]
    ))
    paste0("{", paste(sort(trimws(genes)), collapse = ","), "}")
  }
  
  dataset_rules <- vector("list", length(datasets))
  names(dataset_rules) <- names(datasets)
  
  for (data_name in names(datasets)) {
    print(str_glue("Working on the {data_name} dataset"))
    current_df <- datasets[[data_name]]
    if (is.null(current_df) || nrow(current_df) == 0) next
    
    all_rules_for_dataset <- vector("list", length(unique(current_df$Year)))
    
    for (year in sort(unique(current_df$Year))) {
      print(str_glue("Working on the year {year}"))
      year_data <- current_df[current_df$Year == year, ]
      if (nrow(year_data) == 0) next
      
      tryCatch({
        genotype_rules <- get_rules_or_itemsets(
          df = year_data, 
          resistance_indicator = "genotype",
          year = year, 
          target = target, 
          agg = FALSE
        )
        
        if (length(genotype_rules) == 0) next
        
        # Calculate measures once
        measures_df <- interestMeasure(genotype_rules, measure = c("cosine", "support"))
        
        rule_df <- DATAFRAME(genotype_rules)
        
        # Convert factors to characters
        rule_df$LHS <- as.character(rule_df$LHS)
        rule_df$RHS <- as.character(rule_df$RHS)
        
        rule_df$cosine <- measures_df[, "cosine"]
        rule_df$support <- measures_df[, "support"]
        rule_df$Year <- year
        rule_df$Data_Source <- data_name
        
        # Filter empty rules
        valid_rules <- rule_df$LHS != "{}" & rule_df$RHS != "{}"
        rule_df <- rule_df[valid_rules, ]
        
        all_rules_for_dataset[[as.character(year)]] <- rule_df
        
      }, error = function(e) {
        warning(paste("Error processing", data_name, "year", year, ":", e$message))
      })
    }
    
    # Combine all years for this dataset
    combined_rules <- bind_rows(all_rules_for_dataset)
    
    if (nrow(combined_rules) > 0) {
      # Get top 10 rules by support and cosine across all years
      combined_rules <- combined_rules[order(-combined_rules$support, -combined_rules$cosine, -combined_rules$confidence), ]
      combined_rules <- head(combined_rules, 1000)
      
      # Fix gene names
      combined_rules$LHS <- process_genes_vec(combined_rules$LHS)
      combined_rules$RHS <- process_genes_vec(combined_rules$RHS)
      
      # Create standardized combinations and remove duplicates
      combined_rules$gene_combination <- mapply(standardize_combination, 
                                              combined_rules$LHS, 
                                              combined_rules$RHS)
      
      # Keep only the first occurrence of each gene combination
      combined_rules <- combined_rules[!duplicated(combined_rules$gene_combination), ]
      
      dataset_rules[[data_name]] <- combined_rules
    }
  }
  
  # Combine all datasets' top rules
  full_df <- bind_rows(dataset_rules)
  write.csv(full_df, "combined_csv_outputs/Unique_top_1000_rules_by_dataset.csv", row.names = FALSE)
  
  return(full_df)
}








save_top_gene_rules_aggregated_v3 <- function(Retail_Meats_df, cecal_df, NAHLN_df, target) {
  # Pre-compile regex patterns and create lookup tables outside the loops
  gene_class_mappings <- read.csv("gene_class_mappings.csv")
  sanitized_genes <- make.names(gene_class_mappings$Genes)
  gene_name_map <- setNames(gene_class_mappings$Genes, sanitized_genes)
  split_pattern <- ",|\\{|\\}"
  
  # Process each dataset
  datasets <- list(
    "Retail_Meats" = Retail_Meats_df
    # "Cecal" = cecal_df,
    # "NAHLN" = NAHLN_df
  )
  
  # Vectorized function to process gene names
  process_genes_vec <- function(genes_vec) {
    vapply(genes_vec, function(genes_str) {
      if (!is.character(genes_str)) return(NA_character_)
      
      genes <- strsplit(genes_str, split_pattern)[[1]]
      genes <- trimws(genes[genes != ""])
      fixed <- vapply(genes, function(g) {
        if (g %in% names(gene_name_map)) gene_name_map[g] else g
      }, character(1))
      paste0("{", paste(fixed, collapse = ","), "}")
    }, character(1))
  }
  
  # Function to standardize combinations
  standardize_combination <- function(lhs, rhs) {
    genes <- unique(c(
      strsplit(gsub("[{}]", "", lhs), ",")[[1]],
      strsplit(gsub("[{}]", "", rhs), ",")[[1]]
    ))
    paste0("{", paste(sort(trimws(genes)), collapse = ","), "}")
  }
  
  dataset_rules <- vector("list", length(datasets))
  names(dataset_rules) <- names(datasets)
  
  for (data_name in names(datasets)) {
    print(str_glue("Working on the {data_name} dataset"))
    current_df <- datasets[[data_name]]
    if (is.null(current_df) || nrow(current_df) == 0) next
    
    all_rules_for_dataset <- vector("list", length(unique(current_df$Year)))
    
    for (year in sort(unique(current_df$Year))) {
      print(str_glue("Working on the year {year}"))
      year_data <- current_df[current_df$Year == year, ]
      if (nrow(year_data) == 0) next
      
      tryCatch({
        genotype_rules <- get_rules_or_itemsets(
          df = year_data, 
          resistance_indicator = "genotype",
          year = year, 
          target = target, 
          agg = FALSE
        )
        
        if (length(genotype_rules) == 0) next
        
        # Calculate measures once
        measures_df <- interestMeasure(genotype_rules, measure = c("cosine", "support"))
        
        rule_df <- DATAFRAME(genotype_rules)
        
        # Convert factors to characters
        rule_df$LHS <- as.character(rule_df$LHS)
        rule_df$RHS <- as.character(rule_df$RHS)
        
        rule_df$cosine <- measures_df[, "cosine"]
        rule_df$support <- measures_df[, "support"]
        rule_df$Year <- year
        rule_df$Data_Source <- data_name
        
        # Filter empty rules
        valid_rules <- rule_df$LHS != "{}" & rule_df$RHS != "{}"
        rule_df <- rule_df[valid_rules, ]
        
        all_rules_for_dataset[[as.character(year)]] <- rule_df
        
      }, error = function(e) {
        warning(paste("Error processing", data_name, "year", year, ":", e$message))
      })
    }
    
    # Combine all years for this dataset
    combined_rules <- bind_rows(all_rules_for_dataset)
    
    if (nrow(combined_rules) > 0) {
      # Add class information to rules - vectorized version
      get_classes <- function(genes) {
        match_idx <- match(genes, gene_class_mappings$Genes)
        classes <- gene_class_mappings$corrected_class[match_idx]
        unique(classes[!is.na(classes)])
      }
      
      # Vectorized class assignment using data.table
      combined_rules_dt <- as.data.table(combined_rules)
      
      # Create gene_combination column first
      combined_rules_dt[, gene_combination := {
        Map(function(lhs, rhs) {
          genes <- unique(c(
            strsplit(gsub("[{}]", "", lhs), ",")[[1]],
            strsplit(gsub("[{}]", "", rhs), ",")[[1]]
          ))
          paste0("{", paste(sort(trimws(genes)), collapse = ","), "}")
        }, LHS, RHS)
      }]
      
      combined_rules_dt[, genes := {
        genes_list <- Map(function(lhs, rhs) {
          unique(unlist(strsplit(gsub("[{}]", "", c(lhs, rhs)), ",")))
        }, LHS, RHS)
        genes_list
      }]
      
      combined_rules_dt[, classes := {
        vapply(genes, function(x) {
          paste(sort(get_classes(trimws(x))), collapse = ",")
        }, character(1))
      }]
      
      # First, take top 200 rules by support/cosine using data.table
      setorder(combined_rules_dt, -support, -cosine)
      top_rules_dt <- head(combined_rules_dt, min(200, nrow(combined_rules_dt)))
      
      # Then select remaining rules using class diversity
      remaining_slots <- 1000 - nrow(top_rules_dt)
      if (remaining_slots > 0) {
        # Use data.table for remaining rules
        remaining_rules_dt <- combined_rules_dt[!gene_combination %in% top_rules_dt$gene_combination]
        
        # Pre-split classes for faster access
        all_classes <- unlist(strsplit(c(top_rules_dt$classes), ","))
        class_freq <- table(all_classes)
        
        selected_rules_dt <- data.table()
        
        # Process in batches for better performance
        batch_size <- 1000
        while(nrow(selected_rules_dt) < remaining_slots && nrow(remaining_rules_dt) > 0) {
          # Update class frequencies less frequently
          if(nrow(selected_rules_dt) %% batch_size == 0) {
            current_classes <- unlist(strsplit(c(top_rules_dt$classes, selected_rules_dt$classes), ","))
            class_freq <- table(current_classes)
          }
          
          # Calculate scores for current batch
          remaining_rules_dt[, score := {
            vapply(classes, function(class_str) {
              these_classes <- unlist(strsplit(class_str, ","))
              sum(1 / (class_freq[these_classes] + 1), na.rm = TRUE) * 
                .SD$support * .SD$cosine
            }, numeric(1))
          }]
          
          # Select top scoring rules from batch
          setorder(remaining_rules_dt, -score)
          n_select <- min(batch_size, remaining_slots - nrow(selected_rules_dt))
          selected_batch <- head(remaining_rules_dt, n_select)
          selected_rules_dt <- rbindlist(list(selected_rules_dt, selected_batch))
          remaining_rules_dt <- remaining_rules_dt[-(1:nrow(selected_batch))]
        }
        
        # Combine results
        combined_rules_dt <- rbindlist(list(top_rules_dt, selected_rules_dt))
      } else {
        combined_rules_dt <- top_rules_dt
      }
      
      # Convert back to data.frame and continue with existing logic
      combined_rules <- as.data.frame(combined_rules_dt)
      combined_rules$LHS <- process_genes_vec(combined_rules$LHS)
      combined_rules$RHS <- process_genes_vec(combined_rules$RHS)
      
      dataset_rules[[data_name]] <- combined_rules[!duplicated(combined_rules$gene_combination), ]
    }
  }
  
  # Combine all datasets' top rules
  full_df <- bind_rows(dataset_rules)
  write.csv(full_df, "combined_csv_outputs/ClassDiverse_Unique_top_1000_rules_by_dataset.csv", row.names = FALSE)
  
  return(full_df)
}


# Retail_Meats_df <- read.csv("Retail_Meats/Retail_Meats_wide_corrected_genotype.csv")
# cecal_df <- read.csv("cecal/cecal_wide_corrected_genotype.csv")
# NAHLN_df <- read.csv("NAHLN/NAHLN_wide_corrected_genotype.csv")


# save_top_gene_rules_aggregated_v3(Retail_Meats_df = Retail_Meats_df, cecal_df = cecal_df, NAHLN_df = NAHLN_df, target = "rules")


add_gene_family_combination <- function(rules_df) {
  # Read gene class mappings
  gene_class_mappings <- read.csv("gene_class_mappings.csv")
  
  # Create lookup table for faster mapping
  gene_to_class <- setNames(gene_class_mappings$corrected_class, gene_class_mappings$Genes)
  
  # Function to process a single gene combination
  process_combination <- function(gene_combo) {
    # Remove brackets and split by comma
    genes <- strsplit(gsub("[{}]", "", gene_combo), ",")[[1]]
    genes <- trimws(genes)  # Remove any whitespace
    
    # Look up classes for each gene
    classes <- gene_to_class[genes]
    
    # Remove any NA values and get unique classes
    classes <- unique(classes[!is.na(classes)])
    
    # Return formatted string
    if(length(classes) > 0) {
      return(paste0("{", paste(sort(classes), collapse = ","), "}"))
    } else {
      return("{}")
    }
  }
  
  # Add new column with class combinations
  rules_df$class_combination <- vapply(
    rules_df$gene_combination,
    process_combination,
    character(1)
  )
  
  return(rules_df)
}

# # Read your rules CSV
# rules_df <- read.csv("Unique_top_1000_rules_by_dataset.csv")

# # Add class combinations
# rules_df <- add_gene_family_combination(rules_df)

# # Write updated CSV
# write.csv(rules_df, "Unique_top_1000_rules_by_dataset_with_classes.csv", row.names = FALSE)
# print("DONE")



