

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


library(here)

source(here("function_scripts", "data_wrangling.R"))



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

    #create empty lists to hold the top lodaings later
    pc1_top_loads <- list()
    pc2_top_loads <- list()
    pc3_top_loads <- list()
    pc4_top_loads <- list()



    #apply needed functions to the files (i.e., rule mining, pca, and selecting top lodaings)
    for (year in c(min(df[, "Year"]):max(df[, "Year"]))) {

                #mine rules and return quality measures as data variable
                rule_measures <- mine_associations(df = df, year = year, all_measures = TRUE, 
                resistance_indicator = resistance_indicator, target = target)

                # Remove constant columns
                non_constant_cols <- apply(rule_measures, 2, function(x) length(unique(x)) > 1)
                rule_measures <- rule_measures[, non_constant_cols, drop = FALSE]

                # Check if there are any non-constant columns left
                if (ncol(rule_measures) > 0) {
                    pca <- prcomp(rule_measures, scale = TRUE, center = TRUE, retx = TRUE, rank. = min(4, ncol(rule_measures)))

                    #get the top loadings for each of the first 4 PCs
                    load <- abs(pca$rotation)
                    prop_load <- as.data.frame(apply(load, 2, function(x) x / sum(x)))

                    pc1_loads <- prop_load |> sort_by(~ list(-PC1))
                    pc2_loads <- prop_load |> sort_by(~ list(-PC2))
                    pc3_loads <- prop_load |> sort_by(~ list(-PC3))
                    pc4_loads <- prop_load |> sort_by(~ list(-PC4))

                    #append the year general (i.e, defined outside of lapply) top load list with the top loads for this particular year
                    pc1_top_loads <- list.append(pc1_top_loads, rownames(pc1_loads[1:5, ]))
                    pc2_top_loads <- list.append(pc2_top_loads, rownames(pc2_loads[1:5, ]))
                    pc3_top_loads <- list.append(pc3_top_loads, rownames(pc3_loads[1:5, ]))
                    pc4_top_loads <- list.append(pc4_top_loads, rownames(pc4_loads[1:5, ]))
                } else {
                    warning(paste("No non-constant columns for year", year))
                    next
                }
      }


    pca_df <- data.frame(PC1 = unlist(pc1_top_loads),
                    PC2 = unlist(pc2_top_loads),
                    PC3 = unlist(pc3_top_loads),
                    PC4 = unlist(pc4_top_loads))

    #pca_df[pca_df=="table.n11"] <- "table"
    final_measures <- c(names(which.max(table(pca_df$PC1))),
                      names(which.max(table(pca_df$PC2))),
                      names(which.max(table(pca_df$PC3))),
                      names(which.max(table(pca_df$PC4))))


    if (plothist == TRUE){

    pc1_plot <- ggplot(data = data.frame(pca_df$PC1), aes(x= pca_df$PC1)) + 
    geom_bar() + 
    labs(x = "quality measures", title = "PC1") +
    theme( axis.text.x = element_text(
        angle = 22.5,
        hjust = 1,
        size = 10
      ))
    ggsave(filename = str_glue("{data_source}/figures/measure_pca_toploadings/{resistance_indicator}/{data_source}_{resistance_indicator}_pc1.png"), plot = pc1_plot)

    pc2_plot <- ggplot(data = data.frame(pca_df$PC2), aes(x= pca_df$PC2)) + 
    geom_bar() + 
    labs(x = "quality measures", title = "PC2") +
    theme( axis.text.x = element_text(
        angle = 22.5,
        hjust = 1,
        size = 10
      ))
    ggsave(filename = str_glue("{data_source}/figures/measure_pca_toploadings/{resistance_indicator}/{data_source}_{resistance_indicator}_pc2.png"), plot = pc2_plot)
    

    pc3_plot <- ggplot(data = data.frame(pca_df$PC3), aes(x= pca_df$PC3)) + 
    geom_bar() + 
    labs(x = "quality measures", title = "PC3") +
    theme( axis.text.x = element_text(
        angle = 22.5,
        hjust = 1,
        size = 10
      ))

    ggsave(filename = str_glue("{data_source}/figures/measure_pca_toploadings/{resistance_indicator}/{data_source}_{resistance_indicator}_pc3.png"), plot = pc3_plot)

    pc4_plot <- ggplot(data = data.frame(pca_df$PC4), aes(x= pca_df$PC4)) + 
    geom_bar() + 
    labs(x = "quality measures", title = "PC4") +
    theme( axis.text.x = element_text(
        angle = 22.5,
        hjust = 1,
        size = 10
      ))
    ggsave(filename = str_glue("{data_source}/figures/measure_pca_toploadings/{resistance_indicator}/{data_source}_{resistance_indicator}_pc4.png"), plot = pc4_plot)
    }
    #final_measures <- append(final_measures, c("support", "lift"))
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
    print("###")

    rule_measures <- mine_associations(df = df, year, all_measures = FALSE, use_measures = selected_measures,
                resistance_indicator = resistance_indicator, target = target)

  #incidence_matrix <- incidence_matrix[incidence_matrix$Year == year, "ID"]
  sample_size <- length(c(unique(df[df$Year == year, "ID"]))) 




  rule_measures <- as.data.frame(rule_measures)
  #rule_measures <- subset(rule_measures, select = -c(table.n00, table.n10, table.n01)) # nolint
  
  
  
  ggplot(gather(rule_measures), aes(value)) + 
    geom_histogram(bins = 50) + 
    facet_wrap(~key, scales = 'free_x') +
    labs(y = "number of rules", title = str_glue('year = {year}; total rules = {nrow(rule_measures)}; sample size = {sample_size}'))
    #geom_density(color = 2)
    ggsave(str_glue("{data_source}/figures/r_measure_dists/{resistance_indicator}/{data_source}_{year}_{resistance_indicator}_measure_distributions.png"))

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
  #read in all the files

  #create empty list to hold matrices of data
  lom <- list()
  
  
  
  for (year in c(min(df[, "Year"]):max(df[, "Year"]))) {

    #create matrix to hold data
    data_matrix <- matrix(0, 10, 10) #50 subsamples repeated and averaged 100 times


    #read in the data
    year_df <- df[df$Year == year, ]
    #get 50 different evenly spaced subsample sizes

    subsample_size <- floor(seq(from = 20,to = nrow(year_df), length.out = 10))





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
    ggsave(str_glue("{data_source}/figures/rarefaction_curves/{resistance_indicator}/{data_source}_{resistance_indicator}_rarefaction_by_sample_size_{target}.png"))


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


    #year_df <- df[df$Year == year, ]





    #create matrix to hold data
    data_matrix <- matrix(0, 10, 10) #10 subsamples repeated and averaged 10 times

    #do all the stuff to run apriori
    rules <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator,
                                   year = year, target = target)
    
    
    
    pattern <- str_extract_all(labels(rules), "\\S+")



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
      else {stop("ERROR: resistance_indicator must be either phenotype or genotype")}


    #get 10 different evenly spaced subsample sizes
    subsample_size <- floor(seq(from = 5,to = length(pattern), length.out = 10))


    #also define some variables to index and update this matrix inside the inner forloop
    data_col <- 1
    while (data_col <= ncol(data_matrix)) {
    data_row <- 1

        #loop through all the subample sizes
    for (size in subsample_size){
    

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
    year <- c(replicate(10, as.numeric(year)))
    data_matrix <- cbind(data_matrix, year)

   
    data_matrix <- subset(data_matrix, select = c(data_mean, data_stnd_err, subsample_size, year))
   

    lom <- list.append(lom, data_matrix)


  }

  lom_df <- do.call(rbind, lapply(lom, data.frame))
  #lom_df <- lom_df[-which(lom_df$year == 2006),]
  #lom_df <- lom_df[-which(lom_df$year == 2013),]

  lom_df$year <- as.factor(lom_df$year)
  plot <- ggplot(lom_df, aes(x = subsample_size, y = data_mean, 
  ymin = (data_mean - data_stnd_err), ymax = (data_mean + data_stnd_err), 
  colour = year, fill=year, linetype=year)) +
    geom_line(linewidth = 2) +
    geom_point(size = 7) +
    theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = 'white', colour = 'black')) +
    labs(y = str_glue("Drug Class {target}"), x = str_glue("Antimicrobial {target}"), title = str_glue("Rarefaction Curve ({resistance_indicator})")) +
    geom_ribbon(alpha=0.2) 
    #geom_vline(xintercept = 1000)
    print("DONE!")
    ggsave(filename = str_glue("{data_source}/figures/rarefaction_curves/{resistance_indicator}/{data_source}_{resistance_indicator}_rarefaction_by_class_{target}.png"))
    print("DONE!")
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


    #year_df <- df[df$Year == year, ]





    #create matrix to hold data
    data_matrix <- matrix(0, 10, 10) #10 subsamples repeated and averaged 10 times

    #do all the stuff to run apriori
    rules <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator,
                                   year = year, target = target)
    
    
    
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
      else {stop("ERROR: resistance_indicator must be either phenotype or genotype")}


    #get 10 different evenly spaced subsample sizes
    subsample_size <- floor(seq(from = 5,to = length(pattern), length.out = 10))


    #also define some variables to index and update this matrix inside the inner forloop
    data_col <- 1
    while (data_col <= ncol(data_matrix)) {
    data_row <- 1

        #loop through all the subample sizes
    for (size in subsample_size){
    

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
    year <- c(replicate(10, as.numeric(year)))
    data_matrix <- cbind(data_matrix, year)

   
    data_matrix <- subset(data_matrix, select = c(data_mean, data_stnd_err, subsample_size, year))
   

    lom <- list.append(lom, data_matrix)


  }

  lom_df <- do.call(rbind, lapply(lom, data.frame))
  #lom_df <- lom_df[-which(lom_df$year == 2006),]
  #lom_df <- lom_df[-which(lom_df$year == 2013),]

  print(lom_df)
  lom_df$year <- as.factor(lom_df$year)
  plot <- ggplot(lom_df, aes(x = subsample_size, y = data_mean, 
  ymin = (data_mean - data_stnd_err), ymax = (data_mean + data_stnd_err), 
  colour = year, fill=year, linetype=year)) +
    geom_line(linewidth = 2) +
    geom_point(size = 7) +
    theme(text = element_text(size = 20), plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = 'white', colour = 'black')) +
    labs(y = str_glue("Unique Classes Represented"), x = str_glue("Antimicrobial {target}"), title = str_glue("Rarefaction Curve ({resistance_indicator})")) +
    geom_ribbon(alpha=0.2) +
    geom_vline(xintercept = 10) +
    geom_vline(xintercept = 1000)
    print("DONE!")
    ggsave(filename = str_glue("{data_source}/figures/rarefaction_curves/{resistance_indicator}/{data_source}_{resistance_indicator}_unique_classes_{target}.png"))
    print("DONE!")
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


    #create matrix to hold data
    data_matrix <- matrix(0, 10, 2) #50 subsamples repeated and averaged 100 times

    #do all the stuff to run apriori

    rules <- get_rules_or_itemsets(df = df, resistance_indicator = resistance_indicator,
                                   year = year, target = target)



    measures_df <- interestMeasure(rules, measure = measure)


 

    cut_off_value <- seq(from = low,to = high, length.out = 10)

    
    data_row <- 1

    for (cut_off in cut_off_value) {
      new_rules <- rules


    delete_index <- c(which((measures_df < cut_off)))

    if (length(delete_index != 0)) {new_rules <- rules[-delete_index]}
    else { new_rules <- rules}
    


    pattern <- str_extract_all(labels(new_rules), "\\S+")
    

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
      else {stop("ERROR: resistance_indicator must be either phenotype or genotype")}

   

    data_matrix[data_row, 1] <- length(unique(pattern))
    data_matrix[data_row, 2] <- length(new_rules)

    data_row <- data_row + 1

    }
    
    
    


     
    

 
    colnames(data_matrix) <- c("class", "drug")
    data_matrix <- cbind(data_matrix, cut_off_value)
    year <- c(replicate(10, as.numeric(year)))
    data_matrix <- cbind(data_matrix, year)

   
    data_matrix <- subset(data_matrix, select = c(class, drug, cut_off_value, year))
   

    lom <- list.append(lom, data_matrix)


  }

  lom_df <- do.call(rbind, lapply(lom, data.frame))
  #lom_df <- lom_df[-which(lom_df$year == 2006),]
  #lom_df <- lom_df[-which(lom_df$year == 2013),]


  lom_df$year <- as.factor(lom_df$year)


  plot1 <- ggplot(lom_df, aes(x = cut_off_value, y = class, 
  colour = year, fill=year)) +
    geom_line(linewidth = 2) +
    geom_point(size = 3) +
    theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = 'white', colour = 'black')) +
    labs(y = str_glue("Class {target}"), x = NULL, title = str_glue("{measure} ({resistance_indicator})"))
    #geom_hline(yintercept = 1000)
   

  plot2 <- ggplot(lom_df, aes(x = cut_off_value, y = drug, 
  colour = year, fill=year)) +
    geom_line(linewidth = 2) +
    geom_point(size = 3) +
    geom_vline(xintercept = case_when(
      measure == "cosine" ~ 0.5,
      measure == "jaccard" ~ 0.0,
      measure == "kulczynski" ~ 0.5,
      measure == "support" ~ 0.0,
      TRUE ~ 0.0
    ), color = "red", linewidth = 1) +
    theme(text = element_text(size = 15), plot.title = element_text(hjust = 0.5),
    panel.background = element_rect(fill = 'white', colour = 'black')) +
    labs(y = str_glue("number of {target}"), x = "Cut-off Value", title = str_glue("{measure} ({resistance_indicator})"))
    #geom_hline(yintercept = 1000)
    


  legend <- get_legend(plot2)

  # build grid without legends
  pgrid <- plot_grid(plot1 + theme(legend.position="none"), 
  plot2 + theme(legend.position="none"), 
  ncol = 1)
 
 # add legend
  p <- plot_grid(pgrid, legend, ncol = 2, rel_widths = c(1, .12))


#save plot with height of 4 inches and width of 6 inches
  ggsave(plot = plot2, filename = str_glue("{data_source}/figures/cut_off_selection/{resistance_indicator}/{data_source}_{resistance_indicator}_{target}_v_cutoffs_{measure}.png"), height = 4, width = 6)
}
}
#cecal_phenotype_df <- read.csv("cecal/cecal_wide_resStatus_phenotype.csv")
#rules_v_cutoffs(df = cecal_phenotype_df, resistance_indicator = "phenotype", target = "rules", measures = "jaccard", low = 0, high = 1, data_source = "cecal")





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

write.csv(full_df, str_glue("{data_source}/ruleData/{data_source}_genotypeAndPhenotype_{target}.csv"), row.names = FALSE)

}

