library(shiny)
library(DBI)
library(RSQLite)
library(VennDiagram)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(dplyr)
library(shinycssloaders)
#For enrichment analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(DT)
library(rsconnect)

#print(getwd())

db_path <- file.path(getwd(), "TissueXplorer.db")

#rsconnect::writeManifest()

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")


#remove.packages("fgsea")


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
        numericInput("tpm_thresh", "TPM Threshold", value = 0.2, min = 0)
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
                            tags$li("Upload your own gene list or dataset under 'Add User Data', just make sure your uploaded CSV files follow the required format described in the upload section.")
                          ),
                          p("Enjoy TissueXplorer!")
                   )
                 )
        ),
        tabPanel("Venn Diagram & Gene List",
                 div(style = "display: flex; justify-content: center; align-items: center;",
                     plotOutput("vennPlot", width = "600px", height = "500px")
                 ),
                 selectInput("intersection_choice", "Select Intersection to View Genes", choices = NULL),
                 dataTableOutput("single_intersection_table"),
                 downloadButton("download_gene_table", "Download Selected Genes")
        ),
        tabPanel("Heatmap",
                 withSpinner(plotOutput("heatmapPlot"), type = 6)
        ),
        tabPanel("Gene Search",
                 verbatimTextOutput("geneSearchResult")
        ),
        tabPanel("Enrichment Analysis",
                 selectInput("enrichment_group", "Select Gene Set for Enrichment",
                             choices = NULL),
                 actionButton("run_enrichment", "Run GO Enrichment"),
                 plotOutput("go_barplot"),
                 DT::dataTableOutput("go_table")
                 
        ),
        tabPanel("Add User Data",
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
        cex = 2.3,
        cat.cex = 1.5,
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
      
      # FORCE: flatten to a 1D vector, no matter what
      clean_gene_list <- as.character(na.omit(unique(unlist(genes))))
      
      # Now save as a proper 1-column CSV
      write.csv(data.frame(GeneSymbol = clean_gene_list), file, row.names = FALSE, quote = TRUE)
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
  
  
  output$heatmapPlot <- renderPlot({
    req(input$tissue)
    
    withProgress(message = "Generating heatmap...", value = 0, {
      all_species <- species_list_all()
      tissue <- input$tissue
      gene_sets <- list()
      
      incProgress(0.1, detail = "Filtering genes by TPM...")
      
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
        
        incProgress(0.3 / length(all_species), detail = paste("Processed:", nice_name(species)))
      }
      
      if (length(gene_sets) < 2) {
        return()
      }
      
      incProgress(0.5, detail = "Building heatmap matrix...")
      
      species_names <- names(gene_sets)
      mat <- matrix(0, nrow = length(species_names), ncol = length(species_names))
      rownames(mat) <- colnames(mat) <- nice_name(species_names)
      
      for (i in seq_along(gene_sets)) {
        for (j in seq_along(gene_sets)) {
          mat[i, j] <- length(intersect(gene_sets[[i]], gene_sets[[j]]))
        }
      }
      
      incProgress(0.9, detail = "Rendering heatmap...")
      
      pheatmap::pheatmap(
        mat,
        display_numbers = matrix(as.character(round(mat, 0)), nrow = nrow(mat)),
        number_format = "%.0f",
        number_color = "black",
        fontsize_number = 14,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        main = paste("Shared Genes in", nice_name(tissue))
      )
      
      incProgress(1, detail = "Done")
    })
  })
  
  
  
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
    
    if (ncol(new_data) == 1) {
      # Gene list only — just add with placeholder tissue
      new_data$PlaceholderTissue <- 1  # Will ensure compatibility with filtering
      output$add_dataset_status <- renderText(
        paste("✔️ Gene list", input$custom_dataset_name, "added with placeholder tissue.")
      )
    } else {
      output$add_dataset_status <- renderText(
        paste("✔️ Dataset", input$custom_dataset_name, "added successfully!")
      )
    }
    
    custom_datasets$data[[input$custom_dataset_name]] <- new_data
  })
  
  
  session$onSessionEnded(function() {
    dbDisconnect(conn)
  })
  
  
  # === GO ENRICHMENT ANALYSIS ===
  
  # Update enrichment group choices when user clicks "Compare"
  observeEvent(input$goButton, {
    gene_sets <- get_gene_sets()
    species_names <- names(gene_sets)
    n <- length(gene_sets)
    choices <- list()
    
    if (n == 2) {
      nameA <- nice_name(species_names[1])
      nameB <- nice_name(species_names[2])
      choices <- c(
        species_names,
        paste0(species_names[1], "_AND_", species_names[2])
      )
      names(choices) <- c(
        nice_name(species_names),
        paste0(nameA, " & ", nameB)
      )
    } else if (n == 3) {
      nameA <- nice_name(species_names[1])
      nameB <- nice_name(species_names[2])
      nameC <- nice_name(species_names[3])
      
      choices <- c(
        species_names,
        paste0(species_names[1], "_AND_", species_names[2]),
        paste0(species_names[1], "_AND_", species_names[3]),
        paste0(species_names[2], "_AND_", species_names[3]),
        paste0(species_names[1], "_AND_", species_names[2], "_AND_", species_names[3])
      )
      
      names(choices) <- c(
        nice_name(species_names),
        paste0(nameA, " & ", nameB),
        paste0(nameA, " & ", nameC),
        paste0(nameB, " & ", nameC),
        paste0(nameA, " & ", nameB, " & ", nameC)
      )
    }
    
    updateSelectInput(session, "enrichment_group", choices = choices)
  })
  
  
  
  # Perform GO enrichment when user clicks "Run GO Enrichment"
  enrich_result <- eventReactive(input$run_enrichment, {
    req(input$enrichment_group, input$tissue)
    
    gene_sets <- get_gene_sets()
    group_key <- input$enrichment_group
    selected_genes <- NULL
    
    if (group_key %in% names(gene_sets)) {
      selected_genes <- gene_sets[[group_key]]
    } else if (grepl("_AND_", group_key)) {
      group_parts <- strsplit(group_key, "_AND_")[[1]]
      selected_genes <- Reduce(intersect, gene_sets[group_parts])
    }
    
    if (length(selected_genes) == 0) return(NULL)
    
    # Wrap everything in withProgress
    withProgress(message = "Running enrichment analysis...", value = 0, {
      incProgress(0.2, detail = "Converting gene symbols...")
      
      entrez_ids <- bitr(selected_genes,
                         fromType = "SYMBOL",
                         toType = "ENTREZID",
                         OrgDb = org.Hs.eg.db)
      
      if (is.null(entrez_ids) || nrow(entrez_ids) == 0) return(NULL)
      
      incProgress(0.5, detail = "Running enrichGO...")
      
      result <- enrichGO(
        gene = entrez_ids$ENTREZID,
        OrgDb = org.Hs.eg.db,
        ont = "BP",  # Biological Process
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
      
      incProgress(1, detail = "Done.")
      result
    })
  })
  
  
  
  # Render barplot of top GO terms
  output$go_barplot <- renderPlot({
    ego <- enrich_result()
    if (!is.null(ego) && nrow(ego) > 0) {
      barplot(ego, showCategory = 10, title = "Top GO Terms (Biological Process)")
    }
  })
  
  # Render GO table
  output$go_table <- DT::renderDataTable({
    ego <- enrich_result()
    if (is.null(ego) || nrow(ego) == 0) {
      return(data.frame(Message = "No enriched GO terms found."))
    }
    
    df <- as.data.frame(ego)[, c("ID", "Description", "GeneRatio", "p.adjust", "Count")]
    DT::datatable(df, options = list(pageLength = 10), rownames = FALSE)
  })
  
  
}

shinyApp(ui = ui, server = server)

