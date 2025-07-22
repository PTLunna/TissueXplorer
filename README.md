TissueXplorer is a user-friendly web application for comparative transcriptomic data analysis. It enables researchers to analyze gene expression across different mammalian species for a given tissue

Available online in https://tissuexplorer.shinyapps.io/TissueXplorer/
&
You can run TissueXplorer locally using RStudio by running the file "app.R" or "storedata.R" in this repository.

app.R          Standard version used for deployment (in https://tissuexplorer.shinyapps.io/TissueXplorer/). For security, it does not store custom datasets.

storedata.R    Stores user-uploaded datasets in your local SQLite database. After adding your own data, simply restart the app, and your custom datasets will appear automatically.
 

To run TissueXplorer locally, you need to install the following R packages:
install.packages(c(
  "shiny", "DBI", "RSQLite", "VennDiagram", "ggplot2", "reshape2",
  "pheatmap", "dplyr", "shinycssloaders", "DT", "rsconnect"
))
install.packages("enrichR")

TissueXplorer has the following features:
- Compare gene expression between 2 or 3 different species within the same tissue
- Venn diagrams for visualizing shared and unique genes
- Heatmaps showing shared gene counts and intersection percentages
- Gene search functionality
- Allows for a selection of a TPM threshold
- Enrichment analysis feature
- Custom dataset upload (TPM tables or gene lists)
- Downloadable results for all visualizations and data tables

