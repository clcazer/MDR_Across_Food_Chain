#' Generate BibTeX citations for R packages and R itself
#' 
#' @param packages Character vector of package names
#' @param output_dir Directory where BibTeX files should be saved (defaults to current directory)
#' @return Invisible list of file paths where citations were saved
#' @examples
#' get_package_citations(c("dplyr", "ggplot2"))
get_package_citations <- function(packages, output_dir = ".") {
  # Create output directory if it doesn't exist
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Initialize list to store file paths
  citation_files <- list()
  
  # First, get R citation
  tryCatch({
    r_cite <- citation()
    r_bibtex <- toBibtex(r_cite)
    r_filename <- file.path(output_dir, "R_citation.bib")
    writeLines(r_bibtex, r_filename)
    citation_files[["R"]] <- r_filename
    message(sprintf("Citation for R saved to %s", r_filename))
  }, error = function(e) {
    warning(sprintf("Could not generate citation for R: %s", e$message))
  })
  
  # Process each package
  for (pkg in packages) {
    # Check if it's an installed package
    if (!pkg %in% rownames(installed.packages())) {
      warning(sprintf("'%s' is not an installed package - skipping", pkg))
      next
    }
    
    tryCatch({
      # Get citation for package
      cite <- citation(pkg)
      
      # Convert to BibTeX format
      bibtex <- toBibtex(cite)
      
      # Create filename
      filename <- file.path(output_dir, paste0(pkg, "_citation.bib"))
      
      # Write to file
      writeLines(bibtex, filename)
      
      # Store file path
      citation_files[[pkg]] <- filename
      
      message(sprintf("Citation for package '%s' saved to %s", pkg, filename))
    }, error = function(e) {
      warning(sprintf("Could not generate citation for package '%s': %s", pkg, e$message))
    })
  }
  
  invisible(citation_files)
}



# Example usage
packages <- c("arules", "MCMCglmm", "coda")
get_package_citations(packages, output_dir = "Write_up/refs")