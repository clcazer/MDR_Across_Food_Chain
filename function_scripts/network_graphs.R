library(igraph)
library(stringr)
library(here)
library(RColorBrewer)
library(pals)
library(data.table)
library(parallel)

source(here("function_scripts", "data_wrangling.R"))
source(here("function_scripts", "rule_mining_and_selection.R"))
source(here("function_scripts", "rule_analysis.R"))






graph_rules <- function(df, target, cut_off, resistance_indicator, measures_used, data_source, rules_selected, agg = FALSE, save_legend = TRUE){
    # Initial column name fixes
    colnames(df) <- gsub("BETA.LACTAM", "BETA-LACTAM", colnames(df))
    colnames(df) <- gsub("FOLATE.PATHWAY.INHIBITOR", "FOLATE-PATHWAY-INHIBITOR", colnames(df))

    # Get unique years from the data
    available_years <- sort(unique(df$Year))

    # Create a custom layout function
    custom_circular_layout <- function(graph, scale = 0.6) {
        node_colors <- V(graph)$color
        node_names <- V(graph)$name
        n <- vcount(graph)
        
        # Create a data frame with node names and colors
        nodes_df <- data.frame(name = node_names, color = node_colors, stringsAsFactors = FALSE)
        
        # Sort the data frame by color
        nodes_df_sorted <- nodes_df[order(nodes_df$color),]
        
        # Create circular coordinates 
        theta <- seq(0, 2*pi, length.out = n+1)[-1]
        layout <- cbind(cos(theta) * scale, sin(theta) * scale)
        
        # Create a named matrix for the layout, with rows ordered by color
        layout_sorted <- matrix(0, nrow = n, ncol = 2)
        rownames(layout_sorted) <- nodes_df_sorted$name
        layout_sorted[,1] <- layout[,1]
        layout_sorted[,2] <- layout[,2]
        
        return(list(layout = layout_sorted, order = nodes_df_sorted$name))
    }

    for (year in available_years){
        cat("Processing year:", year, "\n")

        if (rules_selected == "best"){
            rules <- implement_cuts(df = df, resistance_indicator = resistance_indicator, target = target, cut_off = cut_off, year = year, agg = FALSE, measures_used = measures_used)
        } else if (rules_selected == "all"){
            rules <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, year = year, target = target, agg = FALSE)
        }

        # Check if rules is empty or doesn't contain valid rules
        if (length(rules) == 0 || nrow(as(rules, "data.frame")) == 0) {
            cat("No valid rules found for year", year, "\n")
            next
        }

        # Convert rules to a graph, decomposing LHS
        rules_df <- DATAFRAME(rules, separate = TRUE)

        if (!all(c("RHS", "LHS") %in% colnames(rules_df))) {
            cat("Required columns 'RHS' and 'LHS' not found in rules for year", year, "\n")
            next
        }

        # Decompose rules and remove curly braces
        decomposed_rules <- data.frame()
        for (i in 1:nrow(rules_df)) {
            items <- c(unlist(strsplit(gsub("[{}]", "", as.character(rules_df$LHS[i])), ",")),
                       gsub("[{}]", "", as.character(rules_df$RHS[i])))
            items <- trimws(items)  # Trim whitespace from all items at once
            items <- items[items != ""]  # Remove any empty items
            
            if (length(items) < 2) {
                cat("Skipping rule", i, "due to insufficient items\n")
                next
            }
            
            for (j in 1:(length(items)-1)) {
                for (k in (j+1):length(items)) {
                    new_row <- data.frame(LHS = items[j],
                                          RHS = items[k])
                    decomposed_rules <- rbind(decomposed_rules, new_row)
                }
            }
        }

        # Check if decomposed_rules is empty
        if (nrow(decomposed_rules) == 0) {
            cat("No valid decomposed rules for year", year, "\n")
            next
        }

        # Count occurrences of each rule
        rule_counts <- table(paste(decomposed_rules$LHS, decomposed_rules$RHS, sep = "->"))
        decomposed_rules$count <- rule_counts[paste(decomposed_rules$LHS, decomposed_rules$RHS, sep = "->")]

        # Remove duplicate rules
        decomposed_rules <- unique(decomposed_rules)

        #load in the phenotype class mappings
        phenotype_class_mappings <- read_excel("EcoliBreakPoints.xlsx")
        genotype_class_mappings <- read.csv("gene_class_mappings.csv")
        genotype_class_mappings$Gene_family_original <- genotype_class_mappings$Gene_family
        genotype_class_mappings$Gene_family <- gsub(pattern = "[^A-Za-z0-9]", replacement = ".", genotype_class_mappings$Gene_family)


        # Before the main loop, create a fixed color palette for all possible classes
        all_possible_classes <- unique(c(phenotype_class_mappings$class, genotype_class_mappings$corrected_class))
        fixed_class_colors <- setNames(cols25(length(all_possible_classes)), all_possible_classes)
        #fixed_class_colors <- setNames(brewer.pal(length(all_possible_classes), "Dark2"), all_possible_classes)
        #fixed_class_colors <- setNames(rainbow(length(all_possible_classes)), all_possible_classes)

        #add a column that indicates the class of the item in the LHS and another column that indicates the class of the item in the RHS
        #only do this if agg is false

        if (!agg){
            
        if (resistance_indicator == "phenotype"){

        # Add LHS class
        decomposed_rules$LHS_class <- sapply(decomposed_rules$LHS, function(x) {
            match_row <- which(phenotype_class_mappings$abbreviation == x)
            if (length(match_row) > 0) {
                return(phenotype_class_mappings$class[match_row[1]])
            } else {
                return(NA)
            }
        })

        # Add RHS class
        decomposed_rules$RHS_class <- sapply(decomposed_rules$RHS, function(x) {
            match_row <- which(phenotype_class_mappings$abbreviation == x)
            if (length(match_row) > 0) {
                return(phenotype_class_mappings$class[match_row[1]])
            } else {
                return(NA)
            }
        })
        }

        #do the same things in the case of genotype
        if (resistance_indicator == "genotype"){

        # Add LHS class
        decomposed_rules$LHS_class <- sapply(decomposed_rules$LHS, function(x) {
            match_row <- which(genotype_class_mappings$Gene_family == x)
            if (length(match_row) > 0) {
                return(genotype_class_mappings$corrected_class[match_row[1]])
            } else {
                return(NA)
            }
        })

        # Add RHS class
        decomposed_rules$RHS_class <- sapply(decomposed_rules$RHS, function(x) {
            match_row <- which(genotype_class_mappings$Gene_family == x)   
            if (length(match_row) > 0) {
                return(genotype_class_mappings$corrected_class[match_row[1]])
            } else {
                return(NA)
            }
        })
        }










        # Before creating the graph, convert gene names back to original format if needed
        if (resistance_indicator == "genotype" && !agg) {
            # Create a lookup table for gene name conversion
            gene_lookup <- setNames(
                genotype_class_mappings$Gene_family_original,
                genotype_class_mappings$Gene_family
            )
            print("here's the gene lookup")
            print(gene_lookup)
            # Convert names in decomposed_rules
            decomposed_rules$LHS <- sapply(decomposed_rules$LHS, function(x) {
                if (x %in% names(gene_lookup)) return(gene_lookup[x])
                return(x)
            })
            decomposed_rules$RHS <- sapply(decomposed_rules$RHS, function(x) {
                if (x %in% names(gene_lookup)) return(gene_lookup[x])
                return(x)
            })
        }



        # Create a mapping of items to their classes
        item_classes <- unique(rbind(
            data.frame(item = decomposed_rules$LHS, class = decomposed_rules$LHS_class),
            data.frame(item = decomposed_rules$RHS, class = decomposed_rules$RHS_class)
        ))
        item_classes <- item_classes[!duplicated(item_classes$item), ]
        }

        rules_graph <- graph_from_data_frame(decomposed_rules, directed = FALSE)
        
        # Check if the graph is empty
        if (vcount(rules_graph) == 0) {
            cat("Empty graph for year", year, "\n")
            next
        }

        # Set edge attributes
       
        E(rules_graph)$width <- log(decomposed_rules$count) * 1.1 + 1  # Edge thickness proportional to rule count (normalized)    
        E(rules_graph)$color <- "lightblue"

        # Set vertex attributes
        if (agg) {
            # When agg is TRUE, use node labels (which are class names) for coloring
            V(rules_graph)$color <- ifelse(V(rules_graph)$name %in% names(fixed_class_colors), 
                               fixed_class_colors[V(rules_graph)$name], 
                               "gray")
        } else {
            # When agg is FALSE, use the previous method for coloring
            V(rules_graph)$color <- ifelse(!is.na(match(V(rules_graph)$name, item_classes$item)), 
                               fixed_class_colors[item_classes$class[match(V(rules_graph)$name, item_classes$item)]], 
                               "gray")
            
            # If genotype, check for sub_classes and modify labels
            if (resistance_indicator == "genotype") {
                V(rules_graph)$name <- sapply(V(rules_graph)$name, function(x) {
                    # Match node label to Gene_family
                    match_row <- which(genotype_class_mappings$Gene_family_original == x)
                    if (length(match_row) > 0) {
                        corrected_class <- genotype_class_mappings$corrected_class[match_row[1]]
                        sub_class <- genotype_class_mappings$sub_class[match_row[1]]
                        if (!is.na(sub_class) && !is.na(corrected_class) && sub_class != corrected_class) {
                            # Replace forward slashes with newlines in sub_class
                            sub_class <- gsub("/", "\n", sub_class)
                            return(paste0(x, "\n(", sub_class, ")"))
                        }
                    }
                    return(x)
                })
            }
        }

        # Node size proportional to degree (log normalized)
        node_degrees <- degree(rules_graph)
        V(rules_graph)$size <- (log(node_degrees) * 2 + 1) * 5

        # if (agg == FALSE){
        # # Before plotting, define legend_classes and legend_colors
        # legend_classes <- unique(c(decomposed_rules$LHS_class, decomposed_rules$RHS_class))
        # legend_classes <- legend_classes[!is.na(legend_classes)]
        # legend_colors <- fixed_class_colors[legend_classes]
        # }

        # Open a PNG device
        png(filename = str_glue("{data_source}/figures/network_graphs/{resistance_indicator}/{rules_selected}/agg{agg}/{data_source}_{resistance_indicator}_{target}_network_{year}_agg{agg}.png"), 
            width = 1800, height = 1800, res = 200)

        # Plot the graph
        par(mar=c(0,0,1,0))
        tryCatch({
            if (agg) {
                layout_info <- custom_circular_layout(rules_graph, scale = 0.85)
                layout <- layout_info$layout
            } else {
                layout_info <- custom_circular_layout(rules_graph, scale = 0.85)
                layout <- layout_info$layout
                new_order <- layout_info$order
                rules_graph <- permute(rules_graph, match(V(rules_graph)$name, new_order))
            }
            
            if (vcount(rules_graph) > 0) {
                plot(rules_graph, 
                     layout = layout,
                     rescale = FALSE,
                     vertex.label.cex = 0.8,
                     vertex.label.font = 4,
                     vertex.label.color = "black",
                     vertex.label.family = "sans",
                     edge.arrow.size = 0.5,
                     vertex.label.dist = 0,
                     edge.curved = 0.1,
                     main = str_to_title(str_glue("{data_source} ({year})")))
            } else {
                plot.new()
                text(0.5, 0.5, "No valid rules to plot", cex = 1.5)
            }

        }, error = function(e) {
            cat("Error plotting graph for year", year, ":", conditionMessage(e), "\n")
            plot.new()
            text(0.5, 0.5, paste("Error:", conditionMessage(e)), cex = 1.5)
        })

        # Close the PNG device
        dev.off()

        cat("Graph saved for year", year, "\n")

      
    }



    #SAVE THE LEGEND SEPARATELY
save_legend_image <- function(filename, 
                              all_possible_classes, 
                              fixed_class_colors, 
                              width = 650, 
                              height = 575, 
                              resolution = 150) {
  # Create a blank plot
  png(filename = filename, 
      width = width, 
      height = height, 
      res = resolution)



    # Legend
    legend_items <- c("Edge thickness: Normalized Rule count", 
                       "Node size: Normalized Degree")
    legend_colors <- c("lightblue", "black")  # lightblue for edges, black for nodes
    legend_pch <- c(NA, NA)  # Line for edge thickness, point for node size
    legend_lty <- c(NA, NA)   # Solid line for edge thickness, no line for node size

    legend_items <- c(legend_items, all_possible_classes)
    class_colors <- fixed_class_colors[all_possible_classes]
    legend_colors <- c(legend_colors, class_colors)
    legend_pch <- c(NA, NA, rep(19, length(all_possible_classes)))
    legend_lty <- c(NA, NA, rep(NA, length(all_possible_classes)))
  
  # Plot the legend
  plot.new()
  legend("center", 
         legend = legend_items, 
         col = legend_colors, 
        pch = legend_pch,  # Use the legend_pch values 
         pt.cex = 2, 
         cex = 0.9, 
         bty = "o", 
         bg = "white", 
         box.lwd = 2, 
         xpd = TRUE)
  

  # Draw a line in the legend
  #segments(x0 = 0.0, x1 = 0.1, y0 = 1.02, y1 = 1.1, lty = 1, col = "lightblue", lwd = 2)

  # Close the PNG device
  dev.off()
}
if (save_legend == TRUE) {
save_legend_image(filename = "legend.png", 
                  all_possible_classes = all_possible_classes, 
                  fixed_class_colors = fixed_class_colors)}

}








graph_rules_optimized <- function(df, target, cut_off, resistance_indicator, measures_used, data_source, rules_selected, agg = FALSE, save_legend = FALSE){
    #gsub BETA.LACTAM with BETA-LACTAM in df colnames
    colnames(df) <- gsub("BETA.LACTAM", "BETA-LACTAM", colnames(df))
    #gsub FOLATE.PATHWAY.INHIBITOR with FOLATE-PATHWAY-INHIBITOR in df colnames
    colnames(df) <- gsub("FOLATE.PATHWAY.INHIBITOR", "FOLATE-PATHWAY-INHIBITOR", colnames(df))

    # Create a custom layout function
    custom_circular_layout <- function(graph, scale = 0.6) {
        node_colors <- V(graph)$color
        node_names <- V(graph)$name
        n <- vcount(graph)
        
        # Create a data frame with node names and colors
        nodes_df <- data.frame(name = node_names, color = node_colors, stringsAsFactors = FALSE)
        
        # Sort the data frame by color
        nodes_df_sorted <- nodes_df[order(nodes_df$color),]
        
        # Create circular coordinates 
        theta <- seq(0, 2*pi, length.out = n+1)[-1]
        layout <- cbind(cos(theta) * scale, sin(theta) * scale)
        
        # Create a named matrix for the layout, with rows ordered by color
        layout_sorted <- matrix(0, nrow = n, ncol = 2)
        rownames(layout_sorted) <- nodes_df_sorted$name
        layout_sorted[,1] <- layout[,1]
        layout_sorted[,2] <- layout[,2]
        
        return(list(layout = layout_sorted, order = nodes_df_sorted$name))
    }

    # Load mappings outside the loop
    phenotype_class_mappings <- read_excel("EcoliBreakPoints.xlsx")
    genotype_class_mappings <- read.csv("wNAHLN_gene_class_mappings.csv")
    genotype_class_mappings$Gene_family_original <- genotype_class_mappings$Gene_family
    genotype_class_mappings$Gene_family <- gsub(pattern = "[^A-Za-z0-9]", replacement = ".", genotype_class_mappings$Gene_family)

    # Create fixed color palette
    all_possible_classes <- unique(c(phenotype_class_mappings$class, genotype_class_mappings$corrected_class))
    fixed_class_colors <- setNames(cols25(length(all_possible_classes)), all_possible_classes)

    # Process all years using Windows-compatible parallel processing
    if (.Platform$OS.type == "windows") {
        # Create cluster for Windows
        cl <- makeCluster(detectCores() - 1)
        
        # Export required packages and source files to cluster
        clusterEvalQ(cl, {
            library(igraph)
            library(stringr)
            library(data.table)
            library(here)
            library(readxl)
            library(pals)
            library(arules)
            
            # Source required function files
            source(here("function_scripts", "data_wrangling.R"))
            source(here("function_scripts", "rule_mining_and_selection.R"))
            source(here("function_scripts", "rule_analysis.R"))
            source(here("function_scripts", "network_graphs.R"))
        })
        
        # Export all necessary objects
        clusterExport(cl, c(
            "df", "target", "cut_off", "resistance_indicator", 
            "measures_used", "data_source", "rules_selected", 
            "agg", "custom_circular_layout", 
            "phenotype_class_mappings", "genotype_class_mappings",
            "all_possible_classes", "fixed_class_colors"
        ), envir = environment())
        
        # Use available_years instead of sequence
        invisible(parLapply(cl, available_years, function(year) {
            cat("Processing year:", year, "\n")
            print("about to mine the rules")
            # Get rules based on selection method
            rules <- if (rules_selected == "best") {
                implement_cuts(df = df, resistance_indicator = resistance_indicator, 
                             target = target, cut_off = cut_off, year = year, 
                             agg = FALSE, measures_used = measures_used)
            } else {
                get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, 
                                    year = year, target = target, agg = FALSE)
            }
            print("mined the rules!")

            # Early returns for invalid rules
            if (length(rules) == 0 || nrow(as(rules, "data.frame")) == 0) {
                cat("No valid rules found for year", year, "\n")
                return(NULL)
            }

            rules_df <- DATAFRAME(rules, separate = TRUE)
            if (!all(c("RHS", "LHS") %in% colnames(rules_df))) {
                cat("Required columns not found for year", year, "\n")
                return(NULL)
            }

            print("about to decompose rules!")

            # Optimize rule decomposition
            decomposed_rules <- do.call(rbind, lapply(1:nrow(rules_df), function(i) {
                items <- unique(trimws(unlist(strsplit(gsub("[{}]", "", as.character(c(rules_df$LHS[i], rules_df$RHS[i]))), ","))))
                if (length(items) < 2) return(NULL)
                
                # Create all unique pairs of items
                combinations <- t(combn(items, 2))
                data.frame(LHS = combinations[, 1], RHS = combinations[, 2], stringsAsFactors = FALSE)
            }))
        print("done decomposing the rules")

            if (is.null(decomposed_rules) || nrow(decomposed_rules) == 0) {
                cat("No valid decomposed rules for year", year, "\n")
                return(NULL)
            }

            # Convert decomposed_rules to data.table
            decomposed_rules <- as.data.table(decomposed_rules)

            # Count occurrences and remove duplicates
            decomposed_rules[, count := .N, by = .(LHS, RHS)]
            decomposed_rules <- unique(decomposed_rules)

            # Add class information
            if (!agg) {
                if (resistance_indicator == "phenotype") {
                    class_lookup <- setNames(phenotype_class_mappings$class, 
                                          phenotype_class_mappings$abbreviation)
                    decomposed_rules$LHS_class <- class_lookup[decomposed_rules$LHS]
                    decomposed_rules$RHS_class <- class_lookup[decomposed_rules$RHS]
                } else if (resistance_indicator == "genotype") {
                    class_lookup <- setNames(genotype_class_mappings$corrected_class, 
                                          genotype_class_mappings$Gene_family)
                    decomposed_rules$LHS_class <- class_lookup[decomposed_rules$LHS]
                    decomposed_rules$RHS_class <- class_lookup[decomposed_rules$RHS]
                }

                # Convert gene names back if needed
                if (resistance_indicator == "genotype") {
                    gene_lookup <- setNames(
                        genotype_class_mappings$Gene_family_original,
                        genotype_class_mappings$Gene_family
                    )
                    decomposed_rules$LHS <- gene_lookup[decomposed_rules$LHS]
                    decomposed_rules$RHS <- gene_lookup[decomposed_rules$RHS]
                }
            }

            # Create graph
            rules_graph <- graph_from_data_frame(decomposed_rules, directed = FALSE)
            
            if (vcount(rules_graph) == 0) {
                cat("Empty graph for year", year, "\n")
                return(NULL)
            }

            # Set edge attributes
            E(rules_graph)$width <- log(decomposed_rules$count) * 1.1 + 1
            E(rules_graph)$color <- "lightblue"

            # Set vertex attributes
            if (agg) {
                V(rules_graph)$color <- ifelse(V(rules_graph)$name %in% names(fixed_class_colors), 
                                             fixed_class_colors[V(rules_graph)$name], 
                                             "gray")
            } else {
                item_classes <- unique(rbind(
                    data.frame(item = decomposed_rules$LHS, class = decomposed_rules$LHS_class),
                    data.frame(item = decomposed_rules$RHS, class = decomposed_rules$RHS_class)
                ))
                item_classes <- item_classes[!duplicated(item_classes$item), ]
                
                V(rules_graph)$color <- ifelse(!is.na(match(V(rules_graph)$name, item_classes$item)), 
                                             fixed_class_colors[item_classes$class[match(V(rules_graph)$name, item_classes$item)]], 
                                             "gray")
                
                if (resistance_indicator == "genotype") {
                    V(rules_graph)$name <- sapply(V(rules_graph)$name, function(x) {
                        match_row <- which(genotype_class_mappings$Gene_family_original == x)
                        if (length(match_row) > 0) {
                            corrected_class <- genotype_class_mappings$corrected_class[match_row[1]]
                            sub_class <- genotype_class_mappings$sub_class[match_row[1]]
                            if (!is.na(sub_class) && !is.na(corrected_class) && sub_class != corrected_class) {
                                sub_class <- gsub("/", "\n", sub_class)
                                return(paste0(x, "\n(", sub_class, ")"))
                            }
                        }
                        return(x)
                    })
                }
            }

            # Set node sizes
            node_degrees <- degree(rules_graph)
            V(rules_graph)$size <- (log(node_degrees) * 2 + 1) * 5

            # Create and save plot
            png(filename = str_glue("{data_source}/figures/network_graphs/{resistance_indicator}/{rules_selected}/agg{agg}/{data_source}_{resistance_indicator}_{target}_network_{year}_agg{agg}.png"), 
                width = 1800, height = 1800, res = 200)

            par(mar=c(0,0,1,0))
            tryCatch({
                layout_info <- custom_circular_layout(rules_graph, scale = 0.85)
                layout <- layout_info$layout
                if (!agg) {
                    new_order <- layout_info$order
                    rules_graph <- permute(rules_graph, match(V(rules_graph)$name, new_order))
                }
                
                if (vcount(rules_graph) > 0) {
                    plot(rules_graph, 
                         layout = layout,
                         rescale = FALSE,
                         vertex.label.cex = 0.8,
                         vertex.label.font = 4,
                         vertex.label.color = "black",
                         vertex.label.family = "sans",
                         edge.arrow.size = 0.5,
                         vertex.label.dist = 0,
                         edge.curved = 0.1,
                         main = str_to_title(str_glue("{data_source} ({year})")))
                } else {
                    plot.new()
                    text(0.5, 0.5, "No valid rules to plot", cex = 1.5)
                }
            }, error = function(e) {
                cat("Error plotting graph for year", year, ":", conditionMessage(e), "\n")
                plot.new()
                text(0.5, 0.5, paste("Error:", conditionMessage(e)), cex = 1.5)
            })

            dev.off()
            cat("Graph saved for year", year, "\n")
        }))
        
        # Stop cluster
        stopCluster(cl)
    } else {
        # Use mclapply for Unix/Linux/Mac
        invisible(mclapply(available_years, function(year) {
            cat("Processing year:", year, "\n")
            print("about to mine the rules")
            # Get rules based on selection method
            rules <- if (rules_selected == "best") {
                implement_cuts(df = df, resistance_indicator = resistance_indicator, 
                             target = target, cut_off = cut_off, year = year, 
                             agg = FALSE, measures_used = measures_used)
            } else {
                get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator, 
                                    year = year, target = target, agg = FALSE)
            }
            print("mined the rules!")

        # Early returns for invalid rules
        if (length(rules) == 0 || nrow(as(rules, "data.frame")) == 0) {
            cat("No valid rules found for year", year, "\n")
            return(NULL)
        }

        rules_df <- DATAFRAME(rules, separate = TRUE)
        if (!all(c("RHS", "LHS") %in% colnames(rules_df))) {
            cat("Required columns not found for year", year, "\n")
            return(NULL)
        }

        print("about to decompose rules!")

        # Optimize rule decomposition
        decomposed_rules <- do.call(rbind, lapply(1:nrow(rules_df), function(i) {
            items <- unique(trimws(unlist(strsplit(gsub("[{}]", "", as.character(c(rules_df$LHS[i], rules_df$RHS[i]))), ","))))
            if (length(items) < 2) return(NULL)
            
            # Create all unique pairs of items
            combinations <- t(combn(items, 2))
            data.frame(LHS = combinations[, 1], RHS = combinations[, 2], stringsAsFactors = FALSE)
        }))
print("done decomposing the rules")

        if (is.null(decomposed_rules) || nrow(decomposed_rules) == 0) {
            cat("No valid decomposed rules for year", year, "\n")
            return(NULL)
        }

        # Convert decomposed_rules to data.table
        decomposed_rules <- as.data.table(decomposed_rules)

        # Count occurrences and remove duplicates
        decomposed_rules[, count := .N, by = .(LHS, RHS)]
        decomposed_rules <- unique(decomposed_rules)

        # Add class information
        if (!agg) {
            if (resistance_indicator == "phenotype") {
                class_lookup <- setNames(phenotype_class_mappings$class, 
                                      phenotype_class_mappings$abbreviation)
                decomposed_rules$LHS_class <- class_lookup[decomposed_rules$LHS]
                decomposed_rules$RHS_class <- class_lookup[decomposed_rules$RHS]
            } else if (resistance_indicator == "genotype") {
                class_lookup <- setNames(genotype_class_mappings$corrected_class, 
                                      genotype_class_mappings$Gene_family)
                decomposed_rules$LHS_class <- class_lookup[decomposed_rules$LHS]
                decomposed_rules$RHS_class <- class_lookup[decomposed_rules$RHS]
            }

            # Convert gene names back if needed
            if (resistance_indicator == "genotype") {
                gene_lookup <- setNames(
                    genotype_class_mappings$Gene_family_original,
                    genotype_class_mappings$Gene_family
                )
                decomposed_rules$LHS <- gene_lookup[decomposed_rules$LHS]
                decomposed_rules$RHS <- gene_lookup[decomposed_rules$RHS]
            }
        }

        # Create graph
        rules_graph <- graph_from_data_frame(decomposed_rules, directed = FALSE)
        
        if (vcount(rules_graph) == 0) {
            cat("Empty graph for year", year, "\n")
            return(NULL)
        }

        # Set edge attributes
        E(rules_graph)$width <- log(decomposed_rules$count) * 1.1 + 1
        E(rules_graph)$color <- "lightblue"

        # Set vertex attributes
        if (agg) {
            V(rules_graph)$color <- ifelse(V(rules_graph)$name %in% names(fixed_class_colors), 
                                         fixed_class_colors[V(rules_graph)$name], 
                                         "gray")
        } else {
            item_classes <- unique(rbind(
                data.frame(item = decomposed_rules$LHS, class = decomposed_rules$LHS_class),
                data.frame(item = decomposed_rules$RHS, class = decomposed_rules$RHS_class)
            ))
            item_classes <- item_classes[!duplicated(item_classes$item), ]
            
            V(rules_graph)$color <- ifelse(!is.na(match(V(rules_graph)$name, item_classes$item)), 
                                         fixed_class_colors[item_classes$class[match(V(rules_graph)$name, item_classes$item)]], 
                                         "gray")
            
            if (resistance_indicator == "genotype") {
                V(rules_graph)$name <- sapply(V(rules_graph)$name, function(x) {
                    match_row <- which(genotype_class_mappings$Gene_family_original == x)
                    if (length(match_row) > 0) {
                        corrected_class <- genotype_class_mappings$corrected_class[match_row[1]]
                        sub_class <- genotype_class_mappings$sub_class[match_row[1]]
                        if (!is.na(sub_class) && !is.na(corrected_class) && sub_class != corrected_class) {
                            sub_class <- gsub("/", "\n", sub_class)
                            return(paste0(x, "\n(", sub_class, ")"))
                        }
                    }
                    return(x)
                })
            }
        }

        # Set node sizes
        node_degrees <- degree(rules_graph)
        V(rules_graph)$size <- (log(node_degrees) * 2 + 1) * 5

        # Create and save plot
        png(filename = str_glue("{data_source}/figures/network_graphs/{resistance_indicator}/{rules_selected}/agg{agg}/{data_source}_{resistance_indicator}_{target}_network_{year}_agg{agg}.png"), 
            width = 1800, height = 1800, res = 200)

        par(mar=c(0,0,1,0))
        tryCatch({
            layout_info <- custom_circular_layout(rules_graph, scale = 0.85)
            layout <- layout_info$layout
            if (!agg) {
                new_order <- layout_info$order
                rules_graph <- permute(rules_graph, match(V(rules_graph)$name, new_order))
            }
            
            if (vcount(rules_graph) > 0) {
                plot(rules_graph, 
                     layout = layout,
                     rescale = FALSE,
                     vertex.label.cex = 0.8,
                     vertex.label.font = 4,
                     vertex.label.color = "black",
                     vertex.label.family = "sans",
                     edge.arrow.size = 0.5,
                     vertex.label.dist = 0,
                     edge.curved = 0.1,
                     main = str_to_title(str_glue("{data_source} ({year})")))
            } else {
                plot.new()
                text(0.5, 0.5, "No valid rules to plot", cex = 1.5)
            }
        }, error = function(e) {
            cat("Error plotting graph for year", year, ":", conditionMessage(e), "\n")
            plot.new()
            text(0.5, 0.5, paste("Error:", conditionMessage(e)), cex = 1.5)
        })

            dev.off()
            cat("Graph saved for year", year, "\n")
        }, mc.cores = detectCores() - 1))
    }

    # Save legend if requested
    if (save_legend) {
        save_legend_image <- function(filename, 
                                    all_possible_classes, 
                                    fixed_class_colors, 
                                    width = 650, 
                                    height = 575, 
                                    resolution = 150) {
            png(filename = filename, 
                width = width, 
                height = height, 
                res = resolution)

            legend_items <- c("Edge thickness: Normalized Rule count", 
                            "Node size: Normalized Degree")
            legend_colors <- c("lightblue", "black")
            legend_pch <- c(NA, NA)
            legend_lty <- c(NA, NA)

            legend_items <- c(legend_items, all_possible_classes)
            class_colors <- fixed_class_colors[all_possible_classes]
            legend_colors <- c(legend_colors, class_colors)
            legend_pch <- c(NA, NA, rep(19, length(all_possible_classes)))
            legend_lty <- c(NA, NA, rep(NA, length(all_possible_classes)))

            plot.new()
            legend("center", 
                   legend = legend_items, 
                   col = legend_colors, 
                   pch = legend_pch,
                   pt.cex = 2, 
                   cex = 0.9, 
                   bty = "o", 
                   bg = "white", 
                   box.lwd = 2, 
                   xpd = TRUE)

            dev.off()
        }
        
        save_legend_image(filename = "NAHLN_legend.png", 
                         all_possible_classes = all_possible_classes, 
                         fixed_class_colors = fixed_class_colors)
    }
}

# Rprof("profile.out")
# # Run your function
# Rprof(NULL)
# summaryRprof("profile.out")