library(igraph)
library(stringr)
library(here)

source(here("function_scripts", "data_wrangling.R"))
source(here("function_scripts", "rule_mining_and_selection.R"))
source(here("function_scripts", "rule_analysis.R"))






graph_rules <- function(df, target, cut_off, resistance_indicator, measures_used, data_source, rules_selected, agg = FALSE){
    # Create a custom layout function
    custom_circular_layout <- function(graph) {
        node_colors <- V(graph)$color
        node_names <- V(graph)$name
        n <- vcount(graph)
        
        # Create a data frame with node names and colors
        nodes_df <- data.frame(name = node_names, color = node_colors, stringsAsFactors = FALSE)
        
        # Sort the data frame by color
        nodes_df_sorted <- nodes_df[order(nodes_df$color),]
        
        # Create circular coordinates
        theta <- seq(0, 2*pi, length.out = n+1)[-1]
        layout <- cbind(cos(theta), sin(theta))
        
        # Create a named matrix for the layout, with rows ordered by color
        layout_sorted <- matrix(0, nrow = n, ncol = 2)
        rownames(layout_sorted) <- nodes_df_sorted$name
        layout_sorted[,1] <- layout[,1]
        layout_sorted[,2] <- layout[,2]
        
        return(list(layout = layout_sorted, order = nodes_df_sorted$name))
    }

    for (year in c(min(df$Year):max(df$Year))){
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
        genotype_class_mappings$Gene_family <- gsub(pattern = "[^A-Za-z0-9]", replacement = ".", genotype_class_mappings$Genes)


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







        # Create a mapping of items to their classes
        item_classes <- unique(rbind(
            data.frame(item = decomposed_rules$LHS, class = decomposed_rules$LHS_class),
            data.frame(item = decomposed_rules$RHS, class = decomposed_rules$RHS_class)
        ))
        item_classes <- item_classes[!duplicated(item_classes$item), ]
        }
        # Before the main loop, create a fixed color palette for all possible classes
        all_possible_classes <- unique(c(phenotype_class_mappings$class, genotype_class_mappings$corrected_class))
        fixed_class_colors <- setNames(rainbow(length(all_possible_classes)), all_possible_classes)

        rules_graph <- graph_from_data_frame(decomposed_rules, directed = FALSE)
        
        # Check if the graph is empty
        if (vcount(rules_graph) == 0) {
            cat("Empty graph for year", year, "\n")
            next
        }

        # Set edge attributes
        if (max(decomposed_rules$count) <= 5000) {
        E(rules_graph)$width <- ((decomposed_rules$count-1)/(5000-1) + 0.1)* 10  # Edge thickness proportional to rule count (normalized)    
        } else {E(rules_graph)$width <- 1}
        E(rules_graph)$color <- "black"

        # Set vertex attributes
        if (agg) {
            # When agg is TRUE, use node labels (which are class names) for coloring
            V(rules_graph)$color <- fixed_class_colors[V(rules_graph)$name]
        } else {
            # When agg is FALSE, use the previous method
            V(rules_graph)$color <- fixed_class_colors[item_classes$class[match(V(rules_graph)$name, item_classes$item)]]
        }
        # Node size proportional to degree (normalized, but avoids dividing by zero in the case of min=max)
        node_degrees <- degree(rules_graph)
        min_degree <- min(node_degrees)
        max_degree <- max(node_degrees)
        if (max_degree <= 50) {
           V(rules_graph)$size <- 30 * ( 1 +(node_degrees - 1) / (50 - 1))
        } else {
             V(rules_graph)$size <- 70  # Default size if all degrees are the same
        }

        if (agg == FALSE){
        # Before plotting, define legend_classes and legend_colors
        legend_classes <- unique(c(decomposed_rules$LHS_class, decomposed_rules$RHS_class))
        legend_classes <- legend_classes[!is.na(legend_classes)]
        legend_colors <- fixed_class_colors[legend_classes]
        }

        # Open a PNG device
        png(filename = str_glue("{data_source}/figures/network_graphs/{resistance_indicator}/{rules_selected}/agg{agg}/{data_source}_{resistance_indicator}_{target}_network_{year}_agg{agg}.png"), 
            width = 1500, height = 1000, res = 150)

        # Plot the graph
        par(mar=c(1,1,4,4))
        tryCatch({
            if (agg) {
                layout <- layout_in_circle(rules_graph)
            } else {
                layout_info <- custom_circular_layout(rules_graph)
                layout <- layout_info$layout
                new_order <- layout_info$order
                rules_graph <- permute(rules_graph, match(V(rules_graph)$name, new_order))
            }
            
            if (vcount(rules_graph) > 0) {
                plot(rules_graph, 
                     layout = layout,
                     vertex.label.cex = 0.8,
                     vertex.label.font = 2,
                     vertex.label.color = "black",
                     vertex.label.family = "sans",
                     edge.arrow.size = 0.5,
                     vertex.label.dist = 0,
                     edge.curved = 0.1,
                     main = str_glue("{data_source} Decomposed Association Rules Network - {year}"))

                # Legend
                legend_items <- c("Edge thickness: Normalized Rule count", 
                                  "Node size: Normalized Degree")
                legend_colors <- c("lightblue", "black")  # lightblue for edges, black for nodes
                legend_pch <- c(NA, 19)  # Line for edge thickness, point for node size
                legend_lty <- c(1, NA)   # Solid line for edge thickness, no line for node size

                if (!agg) {
                    # Only add class information to legend if there are classes
                    legend_classes <- unique(c(decomposed_rules$LHS_class, decomposed_rules$RHS_class))
                    legend_classes <- legend_classes[!is.na(legend_classes)]
                    if (length(legend_classes) > 0) {
                        legend_items <- c(legend_items, legend_classes)
                        class_colors <- fixed_class_colors[legend_classes]
                        legend_colors <- c(legend_colors, class_colors)
                        legend_pch <- c(legend_pch, rep(19, length(legend_classes)))
                        legend_lty <- c(legend_lty, rep(NA, length(legend_classes)))
                    }
                }

                legend("bottomright", 
                       legend = legend_items,
                       col = legend_colors,
                       pch = legend_pch,
                       lty = legend_lty,
                       pt.cex = 2,  # Increase point size in legend
                       cex = 0.9,
                       bty = "o",
                       bg = "white",
                       box.lwd = 2,
                       xpd = TRUE,
                       inset = c(-0.08, 0))
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
}


