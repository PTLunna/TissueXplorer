library(shiny)
library(DBI)
library(RSQLite)
library(VennDiagram)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(dplyr)
library(shinycssloaders)
# For enrichment analysis
library(enrichR)
library(DT)
library(rsconnect)

#print(getwd())
#rsconnect::writeManifest()


db_path <- file.path(getwd(), "TissueXplorer.db")

# Função para mostrar nomes sem underscores
nice_name <- function(x) gsub("_", " ", x)

ui <- fluidPage(
  tags$div(
    style = "position: absolute; top: 10px; left: 9px;",
    tags$img(src = "logo1.png", height = "60px")
  ),
  div(
    style = "text-align: center;",
    titlePanel(
      div("TissueXplorer", style = "font-family: 'Raleway', sans-serif; font-weight: 700;")
    ),
    tags$h4("A user-friendly bioinformatic tool for transcriptomic data analysis",
            style = "margin-top: -3px; color: #666; font-family: 'Raleway', sans-serif;")
  ),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("n_species", "Number of Species to Compare", 2, min = 2, max = 3),
      uiOutput("species_select"),
      uiOutput("tissue_select"),
      checkboxInput("use_tpm", "Use TPM threshold", value = TRUE),
      conditionalPanel(
        condition = "input.use_tpm == true",
        numericInput("tpm_thresh", "TPM Threshold", value = 1, min = 0)
      ),
      textInput("gene_search", "Search for a Specific Gene Symbol (optional)", placeholder = "e.g., TNMD"),
      actionButton("goButton", "Compare")
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Home",
                 fluidRow(
                   column(12,
                          style = "padding-top: 30px;",
                          
                          div(style = "text-align: center;",
                              tags$img(src = "logo1.png", height = "100px"),
                              h2("Welcome to TissueXplorer!", style = "font-size: 28px; font-weight: bold; margin-bottom: 30px;")
                          ),
                          style = "font-size: 18px; line-height: 1.6;",
                          p("A web tool for transcriptomic data analysis that enables you to compare gene expression data between different species of mammals in the same tissue."),
                          p("This small tutorial will help you get started with the application:"),
                          tags$ol(
                            tags$li("Select the number of species you want to compare."),
                            tags$li("Choose the species and the tissue you are interested in."),
                            tags$li("(Optional) Set a TPM threshold to filter low-expression genes."),
                            tags$li("Click 'Compare' to generate a Venn diagram and gene intersection results."),
                            tags$li("Explore an overview of gene expression similarities across all species for the selected tissue in the Heatmap tab"),
                            tags$li("Use the Gene Search tab to look for one or more specific genes."),
                            tags$li("Run GO enrichment analysis in the Enrichment tab using any subset of genes."),
                            tags$li("Upload your own gene list or dataset under 'Data Integration', just make sure your uploaded CSV files follow the required format described in the upload section.")
                          ),
                          p("Enjoy TissueXplorer!"),
                          p("TissueXplorer uses public transcriptomic data from baseline experiments available in Expression atlas (https://www.ebi.ac.uk/gxa/home) and Bgee (https://www.bgee.org/). If you use tissueXplorer tool in your research, please cite  tissueXplorer and the  corresponding studies in your research for this consult our application note (Article in preparation). This tool was developed for use in evolutionary biology studies and is not intended for use in medical diagnosis.")
                   )
                 )
        ),
        tabPanel("Venn Diagram & Gene List",
                 div(style = "display: flex; justify-content: center; align-items: center;",
                     plotOutput("vennPlot", width = "600px", height = "500px")
                 ),
                 selectInput("intersection_choice", "Select Intersection to View Genes", choices = NULL),
                 dataTableOutput("single_intersection_table"),
                 downloadButton("download_gene_table", "Download Selected Genes"),
                 downloadButton("download_venn", "Download Venn Diagram")
        ),
        tabPanel("Heatmap",
                 withSpinner(plotOutput("heatmapPlot"), type = 6),
                 br(), br(),
                 # Add the plot showing the intersection percentages
                 div(style = "padding: 10px; background-color: #f0f0f0; border-radius: 8px;",
                     h4("Percentage of Shared Genes Between Species (for Selected Tissue)", style = "font-weight: bold;"),
                     withSpinner(plotOutput("intersectionPercentagesPlot"), type = 6),
                     br(),
                     p("The heatmap above shows the percentage of genes shared between each pair of species for the selected tissue and TPM threshold."),
                     p("To read the chart correctly, first locate Species 1 on the vertical axis (rows), then find Species 2 on the horizontal axis (columns). The corresponding value represents the percentage of genes from Species 1 that are shared with Species 2."),
                     p("For example, if the value is 77%, it means that 77% of the total genes expressed in Species 1 of the selected tissue, are also found in Species 2.")
                 ),
                 downloadButton("download_heatmap", "Download Heatmap"),
                 downloadButton("download_percentage_heatmap", "Download Percentage Heatmap")
        ),
        tabPanel("Gene Search",
                 verbatimTextOutput("geneSearchResult")
        ),
        tabPanel("Enrichment Analysis",
                 selectInput("enrichment_group", "Select Gene Set for Enrichment",
                             choices = NULL),
                 selectInput("enrichment_database", "Select Database for Enrichment",
                             choices = c("GO_Biological_Process_2023" = "GO_Biological_Process_2023",
                                         "GO_Molecular_Function_2023" = "GO_Molecular_Function_2023",
                                         "GO_Cellular_Component_2023" = "GO_Cellular_Component_2023",
                                         "KEGG_2021_Human" = "KEGG_2021_Human",
                                         "WikiPathway_2023_Human" = "WikiPathway_2023_Human",
                                         "Reactome_2022" = "Reactome_2022"),
                             selected = "GO_Biological_Process_2023"),
                 numericInput("max_terms", "Maximum number of terms to show", value = 15000, min = 5, max = 20000),
                 checkboxInput("include_nonsignificant", "Include non-significant terms (p > 0.05)", value = FALSE),
                 actionButton("run_enrichment", "Run Enrichment Analysis"),
                 br(), br(),
                 withSpinner(plotOutput("enrichment_barplot"), type = 6),
                 br(),
                 withSpinner(DT::dataTableOutput("enrichment_table"), type = 6),
                 downloadButton("download_enrichment_plot", "Download Enrichment Plot"),
                 downloadButton("download_enrichment_table", "Download Enrichment Table")
        ),
        tabPanel("Data Integration",
                 div(
                   style = "margin-bottom: 20px; padding: 20px; background-color: #f9f9f9; border: 1px solid #ddd; border-radius: 8px;",
                   
                   tags$p("You can upload two different types of files:"),
                   tags$ul(
                     tags$li("Full datasets"),
                     tags$li("Gene lists")
                   ),
                   
                   tags$h4("📄 Full dataset requirements:", style = "font-weight: bold; margin-top: 20px;"),
                   tags$ul(
                     tags$li("The first column must contain the Gene symbols and be named 'GeneSymbol'. The remaining columns should correspond to the names of the tissues and contain the TPM values for each gene."),
                     tags$li("The file must be in CSV format (.csv)."),
                     tags$li("Expression values must be numeric and in TPM (Transcripts per million)."),
                     tags$li("Replace the missing values with 0 and ensure all columns have headers.")
                   ),
                   
                   tags$h4("📄 Gene list requirements:", style = "font-weight: bold; margin-top: 20px;"),
                   tags$ul(
                     tags$li("A single column named 'GeneSymbol' is required. You do not need to specify TPMs."),
                     tags$li("You should specify the tissue in the dataset name."),
                     tags$li("The file must be in CSV format (.csv).")
                   )
                 ),
                 
                 textInput("custom_dataset_name", "Dataset Name", placeholder = "e.g., Mus musculus"),
                 fileInput("custom_file", "Upload CSV File", accept = ".csv"),
                 actionButton("add_dataset_btn", "Add Dataset"),
                 br(),
                 tableOutput("custom_data_preview"),
                 verbatimTextOutput("add_dataset_status")
                 
        )
      )
    )
  )
)


server <- function(input, output, session) {
  conn <- dbConnect(SQLite(), dbname = db_path)
  
  base_species_list <- dbListTables(conn)
  sample_table <- dbReadTable(conn, base_species_list[1])
  tissues <- colnames(sample_table)[!colnames(sample_table) %in% "GeneSymbol"]
  
  custom_datasets <- reactiveValues(data = list())
  intersection_data <- reactiveVal(list())
  
  species_list_all <- reactive({
    c(base_species_list, names(custom_datasets$data))
  })
  
  output$species_select <- renderUI({
    species_choices <- species_list_all()
    species_labels <- nice_name(species_choices)
    names(species_choices) <- species_labels
    lapply(1:input$n_species, function(i) {
      selectInput(paste0("species", i), paste("Select Species", i),
                  choices = species_choices,
                  selected = species_choices[i])
    })
  })
  
  output$tissue_select <- renderUI({
    req(input$n_species)
    
    selected_species <- sapply(1:input$n_species, function(i) input[[paste0("species", i)]])
    available_tissues <- NULL
    
    for (species in selected_species) {
      if (is.null(species)) next
      
      df <- if (species %in% names(custom_datasets$data)) {
        custom_datasets$data[[species]]
      } else {
        dbReadTable(conn, species)
      }
      
      species_tissues <- setdiff(colnames(df), "GeneSymbol")
      
      if (is.null(available_tissues)) {
        available_tissues <- species_tissues
      } else {
        available_tissues <- intersect(available_tissues, species_tissues)
      }
    }
    
    if (length(available_tissues) == 0) {
      return(tags$p("No common tissues found among selected species."))
    }
    
    tissue_labels <- nice_name(available_tissues)
    names(available_tissues) <- tissue_labels
    
    selectInput("tissue", "Select Tissue", choices = available_tissues)
  })
  
  get_gene_sets <- reactive({
    req(input$goButton)
    isolate({
      gene_sets <- list()
      for (i in 1:input$n_species) {
        species <- input[[paste0("species", i)]]
        if (!is.null(species)) {
          df <- if (species %in% names(custom_datasets$data)) {
            custom_datasets$data[[species]]
          } else {
            dbReadTable(conn, species)
          }
          if (input$use_tpm && input$tissue %in% colnames(df)) {
            df <- df %>% filter(.data[[input$tissue]] >= input$tpm_thresh)
          }
          gene_sets[[species]] <- df$GeneSymbol
        }
      }
      gene_sets
    })
  })
  
  # NEW: Create a reactive expression for heatmap gene sets that responds to TPM changes
  get_heatmap_gene_sets <- reactive({
    req(input$tissue)
    
    withProgress(message = "Processing data for heatmaps...", value = 0, {
      all_species <- species_list_all()
      tissue <- input$tissue
      gene_sets <- list()
      
      incProgress(0.1, detail = "Filtering genes by TPM threshold...")
      
      for (i in seq_along(all_species)) {
        species <- all_species[i]
        df <- if (species %in% names(custom_datasets$data)) {
          custom_datasets$data[[species]]
        } else {
          dbReadTable(conn, species)
        }
        
        if (tissue %in% colnames(df)) {
          df_filtered <- df %>% filter(.data[[tissue]] >= input$tpm_thresh)
          gene_sets[[species]] <- df_filtered$GeneSymbol
        }
        
        incProgress(0.7 / length(all_species), detail = paste("Processed:", nice_name(species)))
      }
      
      incProgress(0.9, detail = "Data processing complete")
      
      return(gene_sets)
    })
  })
  
  output$vennPlot <- renderPlot({
    gene_sets <- get_gene_sets()
    n <- length(gene_sets)
    if (n >= 2 && n <= 3) {
      category_names <- nice_name(names(gene_sets))
      cat_pos <- if (n == 2) c(-30, 30) else c(-30, 30, 180)
      cat_dist <- if (n == 2) 0.05 else 0.08
      venn.plot <- venn.diagram(
        gene_sets,
        category.names = category_names,
        filename = NULL,
        fill = c("lightblue", "lightgreen", "lightyellow")[1:n],
        alpha = 0.5,
        cex = 2.5,
        cat.cex = 2,
        cat.pos = cat_pos,
        cat.dist = cat_dist
      )
      grid.draw(venn.plot)
    }
  })
  
  observeEvent(input$goButton, {
    gene_sets <- get_gene_sets()
    species_names <- names(gene_sets)
    n <- length(gene_sets)
    map <- list()
    
    if (n == 2) {
      A <- gene_sets[[1]]; B <- gene_sets[[2]]
      nameA <- nice_name(species_names[1])
      nameB <- nice_name(species_names[2])
      map[[paste0(nameA)]] <- setdiff(A, B)
      map[[paste0(nameB)]] <- setdiff(B, A)
      map[[paste0(nameA, " & ", nameB)]] <- intersect(A, B)
      
    } else if (n == 3) {
      A <- gene_sets[[1]]; B <- gene_sets[[2]]; C <- gene_sets[[3]]
      nameA <- nice_name(species_names[1])
      nameB <- nice_name(species_names[2])
      nameC <- nice_name(species_names[3])
      map[[paste0(nameA)]] <- setdiff(A, union(B, C))
      map[[paste0(nameB)]] <- setdiff(B, union(A, C))
      map[[paste0(nameC)]] <- setdiff(C, union(A, B))
      map[[paste0(nameA, " & ", nameB)]] <- setdiff(intersect(A, B), C)
      map[[paste0(nameA, " & ", nameC)]] <- setdiff(intersect(A, C), B)
      map[[paste0(nameB, " & ", nameC)]] <- setdiff(intersect(B, C), A)
      map[[paste0(nameA, " & ", nameB, " & ", nameC)]] <- intersect(intersect(A, B), C)
    }
    
    intersection_data(map)
    updateSelectInput(session, "intersection_choice", choices = names(map))
  })
  
  output$single_intersection_table <- renderDataTable({
    req(input$intersection_choice)
    genes <- intersection_data()[[input$intersection_choice]]
    data.frame(GeneSymbol = genes)
  })
  
  output$download_gene_table <- downloadHandler(
    filename = function() {
      paste0("genes_", gsub(" ", "_", input$intersection_choice), "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      genes <- intersection_data()[[input$intersection_choice]]
      clean_gene_list <- as.character(na.omit(unique(unlist(genes))))
      write.csv(data.frame(GeneSymbol = clean_gene_list), file, row.names = FALSE, quote = TRUE)
    }
  )
  
  output$download_venn <- downloadHandler(
    filename = function() {
      paste0("venn_diagram_", Sys.Date(), ".png")
    },
    content = function(file) {
      gene_sets <- get_gene_sets()
      n <- length(gene_sets)
      if (n >= 2 && n <= 3) {
        category_names <- nice_name(names(gene_sets))
        cat_pos <- if (n == 2) c(-30, 30) else c(-30, 30, 180)
        cat_dist <- if (n == 2) 0.05 else 0.08
        venn.plot <- venn.diagram(
          gene_sets,
          category.names = category_names,
          filename = NULL,
          fill = c("lightblue", "lightgreen", "lightyellow")[1:n],
          alpha = 0.5,
          cex = 2.5,
          cat.cex = 2,
          cat.pos = cat_pos,
          cat.dist = cat_dist
        )
        png(file, width = 800, height = 600, bg = "transparent")
        grid.draw(venn.plot)
        dev.off()
      }
    }
  )
  
  output$geneSearchResult <- renderPrint({
    req(input$gene_search, input$goButton)
    isolate({
      result_text <- ""
      found_any <- FALSE
      
      # Step 1: Clean and split multiple gene symbols
      query_genes <- unlist(strsplit(input$gene_search, ","))
      query_genes <- trimws(query_genes)  # Remove leading/trailing spaces
      
      for (i in 1:input$n_species) {
        species <- input[[paste0("species", i)]]
        if (!is.null(species)) {
          df <- if (species %in% names(custom_datasets$data)) {
            custom_datasets$data[[species]]
          } else {
            dbReadTable(conn, species)
          }
          
          matches <- df %>%
            filter(GeneSymbol %in% query_genes)
          
          if (input$use_tpm && input$tissue %in% colnames(df)) {
            matches <- matches %>% filter(.data[[input$tissue]] >= input$tpm_thresh)
          }
          
          if (nrow(matches) > 0) {
            found_any <- TRUE
            result_text <- paste0(result_text, "\nSpecies: ", nice_name(species), "\n")
            for (j in 1:nrow(matches)) {
              result_text <- paste0(
                result_text,
                "  - Gene: ", matches$GeneSymbol[j],
                " | Tissue: ", nice_name(input$tissue),
                " | TPM: ", matches[[input$tissue]][j], "\n"
              )
            }
          }
        }
      }
      
      if (!found_any) {
        return("No matching gene(s) found in selected species/tissue with the given threshold.")
      } else {
        cat(result_text)
      }
    })
  })
  
  # FIXED: Make the first heatmap reactive to TPM threshold changes
  output$heatmapPlot <- renderPlot({
    gene_sets <- get_heatmap_gene_sets()  # Use the new reactive expression
    
    if (length(gene_sets) < 2) {
      return()
    }
    
    withProgress(message = "Generating shared genes heatmap...", value = 0, {
      incProgress(0.3, detail = "Building heatmap matrix...")
      
      species_names <- names(gene_sets)
      mat <- matrix(0, nrow = length(species_names), ncol = length(species_names))
      rownames(mat) <- colnames(mat) <- nice_name(species_names)
      
      for (i in seq_along(gene_sets)) {
        for (j in seq_along(gene_sets)) {
          mat[i, j] <- length(intersect(gene_sets[[i]], gene_sets[[j]]))
        }
      }
      
      incProgress(0.8, detail = "Rendering heatmap...")
      
      pheatmap::pheatmap(
        mat,
        display_numbers = matrix(as.character(round(mat, 0)), nrow = nrow(mat)),
        number_format = "%.0f",
        number_color = "black",
        fontsize_number = 14,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        main = paste("Shared Genes in", nice_name(input$tissue))
      )
      
      incProgress(1, detail = "Shared genes heatmap complete")
    })
  })
  
  # NEW: Make the intersection percentages plot independently reactive
  output$intersectionPercentagesPlot <- renderPlot({
    gene_sets <- get_heatmap_gene_sets()  # Use the same reactive expression
    
    if (length(gene_sets) < 2) {
      return()
    }
    
    withProgress(message = "Generating intersection percentages heatmap...", value = 0, {
      species_names <- names(gene_sets)
      
      incProgress(0.2, detail = "Calculating intersection percentages...")
      
      # Calculate the gene intersection percentages for each pair of species
      intersection_percentages <- data.frame(matrix(ncol = length(species_names), nrow = length(species_names)))
      colnames(intersection_percentages) <- rownames(intersection_percentages) <- nice_name(species_names)
      
      total_pairs <- length(species_names)^2
      current_pair <- 0
      
      for (i in seq_along(species_names)) {
        for (j in seq_along(species_names)) {
          current_pair <- current_pair + 1
          A <- gene_sets[[species_names[i]]]
          B <- gene_sets[[species_names[j]]]
          
          # Calculate the number of common and exclusive genes, after applying the TPM threshold
          common_genes <- length(intersect(A, B))
          exclusive_A <- length(setdiff(A, B))
          exclusive_B <- length(setdiff(B, A))
          
          # Calculate the percentage for species A
          if ((exclusive_A + common_genes) > 0) {
            percentage_A <- common_genes / (exclusive_A + common_genes) * 100
          } else {
            percentage_A <- 100
          }
          
          # Calculate the percentage for species B
          if ((exclusive_B + common_genes) > 0) {
            percentage_B <- common_genes / (exclusive_B + common_genes) * 100
          } else {
            percentage_B <- 100
          }
          
          # Assign percentages to the matrix
          intersection_percentages[i, j] <- round(percentage_B, 2)
          intersection_percentages[j, i] <- round(percentage_A, 2)
          
          # Update progress
          if (current_pair %% 3 == 0) {  # Update every 3 pairs to avoid too frequent updates
            incProgress(0.6 * (current_pair / total_pairs), 
                        detail = paste("Processed", current_pair, "of", total_pairs, "pairs"))
          }
        }
      }
      
      incProgress(0.9, detail = "Rendering percentage heatmap...")
      
      # Transpose the matrix to switch axis
      intersection_percentages_t <- t(intersection_percentages)
      
      pheatmap::pheatmap(
        intersection_percentages_t,
        display_numbers = TRUE,
        number_format = "%.2f%%",
        number_color = "black",
        fontsize_number = 12,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        main = paste("Gene Intersection Percentages in", nice_name(input$tissue)),
        legend = FALSE
      )
      
      incProgress(1, detail = "Intersection percentages heatmap complete")
    })
  })
  
  # Heatmap de genes compartilhados
  output$download_heatmap <- downloadHandler(
    filename = function() {
      paste0("shared_genes_heatmap_", Sys.Date(), ".png")
    },
    content = function(file) {
      gene_sets <- get_heatmap_gene_sets()
      species_names <- names(gene_sets)
      mat <- matrix(0, nrow = length(species_names), ncol = length(species_names))
      rownames(mat) <- colnames(mat) <- nice_name(species_names)
      
      for (i in seq_along(gene_sets)) {
        for (j in seq_along(gene_sets)) {
          mat[i, j] <- length(intersect(gene_sets[[i]], gene_sets[[j]]))
        }
      }
      
      png(file, width = 900, height = 700,  bg = "transparent")
      pheatmap::pheatmap(
        mat,
        display_numbers = matrix(as.character(round(mat, 0)), nrow = nrow(mat)),
        number_format = "%.0f",
        number_color = "black",
        fontsize_number = 14,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        main = paste("Shared Genes in", nice_name(input$tissue))
      )
      dev.off()
    }
  )
  
  # Heatmap de porcentagens
  output$download_percentage_heatmap <- downloadHandler(
    filename = function() {
      paste0("intersection_percentages_heatmap_", Sys.Date(), ".png")
    },
    content = function(file) {
      gene_sets <- get_heatmap_gene_sets()
      species_names <- names(gene_sets)
      intersection_percentages <- matrix(0, ncol = length(species_names), nrow = length(species_names))
      colnames(intersection_percentages) <- rownames(intersection_percentages) <- nice_name(species_names)
      
      for (i in seq_along(species_names)) {
        for (j in seq_along(species_names)) {
          A <- gene_sets[[species_names[i]]]
          B <- gene_sets[[species_names[j]]]
          common_genes <- length(intersect(A, B))
          exclusive_A <- length(setdiff(A, B))
          if ((exclusive_A + common_genes) > 0) {
            percentage_A <- common_genes / (exclusive_A + common_genes) * 100
          } else {
            percentage_A <- 100
          }
          intersection_percentages[i, j] <- round(percentage_A, 2)
        }
      }
      
      png(file, width = 900, height = 700,  bg = "transparent")
      pheatmap::pheatmap(
        t(intersection_percentages),
        display_numbers = TRUE,
        number_format = "%.2f%%",
        number_color = "black",
        fontsize_number = 12,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        main = paste("Gene Intersection Percentages in", nice_name(input$tissue)),
        legend = FALSE
      )
      dev.off()
    }
  )
  
  
  
  output$custom_data_preview <- renderTable({
    req(input$custom_file)
    head(read.csv(input$custom_file$datapath))
  })
  
  observeEvent(input$add_dataset_btn, {
    req(input$custom_file, input$custom_dataset_name)
    new_data <- read.csv(input$custom_file$datapath)
    
    if (!"GeneSymbol" %in% colnames(new_data)) {
      output$add_dataset_status <- renderText("❌ ERROR: Your CSV must include a 'GeneSymbol' column.")
      return()
    }
    
    # Clean the name for SQLite compatibility (no spaces, avoid special chars)
    dataset_name <- gsub("[^a-zA-Z0-9_]", "_", input$custom_dataset_name)
    
    if (ncol(new_data) == 1) {
      # Add placeholder tissue for compatibility
      new_data$PlaceholderTissue <- 1
      output$add_dataset_status <- renderText(
        paste("✔️ Gene list", dataset_name, "added with placeholder tissue.")
      )
    } else {
      output$add_dataset_status <- renderText(
        paste("✔️ Dataset", dataset_name, "added successfully!")
      )
    }
    
    # Save to SQLite database (overwrite if table already exists)
    tryCatch({
      dbWriteTable(conn, dataset_name, new_data, overwrite = TRUE)
    }, error = function(e) {
      output$add_dataset_status <- renderText(
        paste("❌ Failed to save dataset to database:", e$message)
      )
    })
  })
  
  
  # === ENRICHMENT ANALYSIS WITH enrichR ===
  
  # Update enrichment group choices when user clicks "Compare"
  observeEvent(input$goButton, {
    # Use the intersection_data that was already calculated for the Venn diagram
    map <- intersection_data()
    if (length(map) > 0) {
      choices <- names(map)
      updateSelectInput(session, "enrichment_group", choices = choices)
    }
  })
  
  # Perform enrichment analysis when user clicks "Run Enrichment Analysis"
  enrich_result <- eventReactive(input$run_enrichment, {
    req(input$enrichment_group, input$enrichment_database)
    
    withProgress(message = "Running enrichment analysis...", value = 0, {
      
      incProgress(0.2, detail = "Getting gene list...")
      
      # Get the intersection data that was calculated for the Venn diagram
      map <- intersection_data()
      selected_genes <- map[[input$enrichment_group]]
      
      if (is.null(selected_genes) || length(selected_genes) == 0) {
        return(list(error = "No genes found for selected group"))
      }
      
      incProgress(0.5, detail = paste("Running enrichment with", length(selected_genes), "genes..."))
      
      # Use enrichR for enrichment analysis
      tryCatch({
        # Get available databases to ensure our selection is valid
        dbs <- listEnrichrDbs()
        
        # Check if the selected database is available
        if (!input$enrichment_database %in% dbs$libraryName) {
          return(list(error = paste("Database", input$enrichment_database, "not available")))
        }
        
        # Run enrichment
        result <- enrichr(selected_genes, input$enrichment_database)
        
        incProgress(1, detail = "Analysis complete!")
        
        if (is.null(result) || length(result) == 0) {
          return(list(error = "No enrichment results returned"))
        }
        
        # Extract the results for the selected database
        enrichment_df <- result[[input$enrichment_database]]
        
        if (is.null(enrichment_df) || nrow(enrichment_df) == 0) {
          return(list(error = "No significant enrichment found"))
        }
        
        # Filter results based on significance setting
        if (input$include_nonsignificant) {
          # Include all results (both significant and non-significant)
          filtered_df <- enrichment_df
          significance_note <- "Showing all terms (including non-significant)"
        } else {
          # Only significant results (p-value < 0.05)
          filtered_df <- enrichment_df[enrichment_df$P.value < 0.05, ]
          significance_note <- "Showing only significant terms (p < 0.05)"
          
          if (nrow(filtered_df) == 0) {
            return(list(error = "No significant terms found (p < 0.05). Try checking 'Include non-significant terms' to see all results."))
          }
        }
        
        # Sort by p-value and take top results
        filtered_df <- filtered_df[order(filtered_df$P.value), ]
        filtered_df <- head(filtered_df, input$max_terms)
        
        return(list(
          data = filtered_df, 
          genes_analyzed = length(selected_genes),
          significance_note = significance_note
        ))
        
      }, error = function(e) {
        return(list(error = paste("Error in enrichment analysis:", e$message)))
      })
    })
  })
  
  # Render barplot of top enriched terms
  output$enrichment_barplot <- renderPlot({
    result <- enrich_result()
    
    if (!is.null(result$error)) {
      # Create an empty plot with error message
      plot(1, 1, type = "n", xlab = "", ylab = "", main = result$error, axes = FALSE)
      return()
    }
    
    if (is.null(result$data) || nrow(result$data) == 0) {
      plot(1, 1, type = "n", xlab = "", ylab = "", main = "No enrichment results to display", axes = FALSE)
      return()
    }
    
    df <- result$data
    
    # Prepare data for plotting
    df$neg_log_p <- -log10(df$P.value)
    df$Term_short <- ifelse(nchar(df$Term) > 50, 
                            paste0(substr(df$Term, 1, 47), "..."), 
                            df$Term)
    
    # Take top terms for plotting
    n_terms <- min(15, nrow(df))
    df_plot <- head(df, n_terms)
    df_plot$Term_short <- factor(df_plot$Term_short, levels = rev(df_plot$Term_short))
    
    # Create barplot
    ggplot(df_plot, aes(x = Term_short, y = neg_log_p)) +
      geom_col(fill = "steelblue", alpha = 0.7) +
      coord_flip() +
      labs(
        title = paste("Enrichment Analysis Results -", input$enrichment_database),
        subtitle = paste("Analyzed", result$genes_analyzed, "genes |", result$significance_note),
        x = "Terms",
        y = "-log10(P-value)"
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10)
      )
  })
  
  output$download_enrichment_plot <- downloadHandler(
    filename = function() {
      paste0("enrichment_plot_", Sys.Date(), ".png")
    },
    content = function(file) {
      result <- enrich_result()
      
      if (!is.null(result$error) || is.null(result$data) || nrow(result$data) == 0) {
        png(file, width = 800, height = 600)
        plot(1, 1, type = "n", xlab = "", ylab = "", main = "No enrichment results to export", axes = FALSE)
        dev.off()
        return()
      }
      
      df <- result$data
      df$neg_log_p <- -log10(df$P.value)
      df$Term_short <- ifelse(nchar(df$Term) > 50, 
                              paste0(substr(df$Term, 1, 47), "..."), 
                              df$Term)
      n_terms <- min(15, nrow(df))
      df_plot <- head(df, n_terms)
      df_plot$Term_short <- factor(df_plot$Term_short, levels = rev(df_plot$Term_short))
      
      p <- ggplot(df_plot, aes(x = Term_short, y = neg_log_p)) +
        geom_col(fill = "steelblue", alpha = 0.7) +
        coord_flip() +
        labs(
          title = paste("Enrichment Analysis Results -", input$enrichment_database),
          subtitle = paste("Analyzed", result$genes_analyzed, "genes |", result$significance_note),
          x = "Terms",
          y = "-log10(P-value)"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 12),
          axis.text.y = element_text(size = 10),
          axis.text.x = element_text(size = 10)
        )
      
      ggsave(file, plot = p, width = 9, height = 6)
    }
  )
  
  
  # Render enrichment results table
  output$enrichment_table <- DT::renderDataTable({
    result <- enrich_result()
    
    if (!is.null(result$error)) {
      return(data.frame(Message = result$error))
    }
    
    if (is.null(result$data) || nrow(result$data) == 0) {
      return(data.frame(Message = "No enrichment results found."))
    }
    
    df <- result$data
    
    # Select relevant columns for display
    display_cols <- c("Term", "P.value", "Adjusted.P.value", "Combined.Score", "Genes")
    df_display <- df[, intersect(display_cols, colnames(df)), drop = FALSE]
    
    # Round P.value and Adjusted.P.value to 6, Combined.Score to 2
    if ("P.value" %in% colnames(df_display)) {
      df_display$P.value <- round(df_display$P.value, 6)
    }
    if ("Adjusted.P.value" %in% colnames(df_display)) {
      df_display$Adjusted.P.value <- round(df_display$Adjusted.P.value, 6)
    }
    if ("Combined.Score" %in% colnames(df_display)) {
      df_display$Combined.Score <- round(df_display$Combined.Score, 2)
    }
    
    
    DT::datatable(
      df_display, 
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        autoWidth = TRUE,
        columnDefs = list(list(width = '300px', targets = c(0, 4))) # Make Term and Genes columns wider
      ), 
      rownames = FALSE,
      caption = paste("Enrichment results for", input$enrichment_database, "-", result$significance_note)
    )
  })
  
  output$download_enrichment_table <- downloadHandler(
    filename = function() {
      paste0("enrichment_table_", Sys.Date(), ".csv")
    },
    content = function(file) {
      result <- enrich_result()
      
      if (!is.null(result$error) || is.null(result$data) || nrow(result$data) == 0) {
        write.csv(data.frame(Message = result$error), file, row.names = FALSE)
        return()
      }
      
      df <- result$data
      display_cols <- c("Term", "P.value", "Adjusted.P.value", "Combined.Score", "Genes")
      df_filtered <- df[, intersect(display_cols, colnames(df)), drop = FALSE]
      
      # Optionally rename columns for CSV (remove dots)
      colnames(df_filtered) <- gsub("\\.", " ", colnames(df_filtered))
      
      write.csv(df_filtered, file, row.names = FALSE)
    }
  )
}
shinyApp(ui = ui, server = server)