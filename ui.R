# metaAML ui.R
#
# Armon Azizi, Brooks Benard
# aazizi@stanford.edu
# bbenard@stanford.edu

# Define main shiny ui
shinyUI(fluidPage(
  navbarPage(inverse=TRUE, title=NULL,
             
             tabPanel("metaAML",
                      
                      titlePanel("Meta analysis of mutations in AML"),
                      
                      # Sidebar panel for user input
                      sidebarPanel(
                        
                        actionButton("analysis_start", "Run Analysis"),
                        
                        # dropdown selection with multiple choice
                        h4("Select Cohorts Of Interest:"),
                        selectizeInput("cohort_choices", 
                                       label=NULL, 
                                       choices = c("All", "De novo", "Secondary", "Relapse", "Therapy related", "Other"), 
                                       selected="All", 
                                       multiple=FALSE),
                        
                        # Slider bar for choosing frequency cutoff
                        h4("Choose Gene Frequency Inclusion Cutoff:"),
                        sliderInput("inclusion_cutoff", label=NULL,
                                    min = 1, max = 100, value = 10),
                        
                        # Select variant or SNV
                        h4("Visualize Variants or VAFs"),
                        selectizeInput("subset_to_variant_of_vaf", 
                                       label=NULL, 
                                       choices = c("Variants", "VAFs"), 
                                       selected="Variants", 
                                       multiple=FALSE),
                        
                        
                        # Gene subset selection
                        h4("Include or Exclude Specific Mutations on plots?"),
                        selectizeInput("subset_to_genes", 
                                       label=NULL, 
                                       choices = c("Default (all genes)", "Include (all rows)", "Include (select rows)", "Exclude"), 
                                       selected="Default (all)", 
                                       multiple=FALSE),
                        
                        # Gene subset selection
                        h4("Select Genes To Subset Mutation Landscape Plot:"),
                        selectizeInput("gene_subsets", 
                                       label=NULL, 
                                       choices = c("All Genes", "ASXL1","ATRX","BCOR", "BCORL1", "BRAF","CBL","CBLB","CDKN2A","CEBPA","CREBBP","CDCF"," CUX1","CSF3R","DNMT3A","EP300","ETV6","EZH2","FBXW7","FLT3","GATA2","GNAS","IDH1","IDH2","IKZF1","JAK2","KDM5A","KDM6A","KIT","KRAS","MLL","MLL2","MLL3","MLL5","MPL","MYC","NF1","NPM1","NRAS","PHF6","PRPF40B","PTEN","PTPN11","RAD21","RB1","RUNX1","SETBP1", "SF1","SF3A1","SF3B1", "SMC1A", "SMC3", "SRSF2","SH2B3","STAG1","STAG2","TET2","TP53","U2AF1","U2AF2","WT1","ZBTB33","ZRSR2"), 
                                       selected="All Genes", 
                                       options = list(maxItems = 10),
                                       multiple=TRUE),
                        
                        # Select variants for VAF correlation
                        h4("Genes for VAF correlation (x axis vs. y axis)"),
                        selectizeInput("gene_scatter", 
                                       label=NULL, 
                                       choices = c("All Genes", "ASXL1","ATRX","BCOR", "BCORL1", "BRAF","CBL","CBLB","CDKN2A","CEBPA","CREBBP","CDCF"," CUX1","CSF3R","DNMT3A","EP300","ETV6","EZH2","FBXW7","FLT3","GATA2","GNAS","IDH1","IDH2","IKZF1","JAK2","KDM5A","KDM6A","KIT","KRAS","MLL","MLL2","MLL3","MLL5","MPL","MYC","NF1","NPM1","NRAS","PHF6","PRPF40B","PTEN","PTPN11","RAD21","RB1","RUNX1","SETBP1", "SF1","SF3A1","SF3B1", "SMC1A", "SMC3", "SRSF2","SH2B3","STAG1","STAG2","TET2","TP53","U2AF1","U2AF2","WT1", "ZBTB33", "ZRSR2"), 
                                       selected= c("NRAS", "KRAS"), 
                                       options = list(maxItems = 2),
                                       multiple=TRUE),
                        
                        # Select variants for 3 sequential VAF ordering
                        h4("Genes for VAF subclonal ordering"),
                        selectizeInput("genes_vaf", 
                                       label=NULL, 
                                       choices = c("All Genes", "ASXL1","ATRX","BCOR", "BCORL1", "BRAF","CBL","CBLB","CDKN2A","CEBPA","CREBBP","CDCF"," CUX1","DNMT3A","EP300","ETV6","EZH2","FBXW7","FLT3","GATA2","GNAS","IDH1","IDH2","IKZF1","JAK2","KDM5A","KDM6A","KIT","KRAS","MLL","MLL2","MLL3","MLL5","MPL","MYC","NF1","NPM1","NRAS","PHF6","PRPF40B","PTEN","PTPN11","RAD21","RB1","RUNX1","SF1","SF3A1","SF3B1", "SMC1A", "SMC3", "SRSF2","SH2B3","STAG1","STAG2","TET2","TP53","U2AF1","U2AF2","WT1","ZBTB33","ZRSR2"), 
                                       selected= c("DNMT3A", "NPM1", "FLT3"), 
                                       options = list(maxItems = 3),
                                       multiple=TRUE),
                        
                        # Select genetic background for survival
                        h4("Select Genetic Background for Survival Analysis"),
                        selectizeInput("genetic_background", 
                                       label=NULL, 
                                       choices = c("None", "DNMT3A", "FLT3", "NPM1", "Complex Cytogenetics", "Male", "Female", "Favorable", "Intermediate", "Adverse"), 
                                       selected="None", 
                                       multiple=FALSE),
                        
                        # Gene selection for survival
                        h4("Select Genes For Survival Analysis:"),
                        selectizeInput("gene_survival", 
                                       label=NULL, 
                                       choices = c("All Genes", "ASXL1","ATRX","BCOR", "BCORL1", "BRAF","CBL","CBLB","CDKN2A","CEBPA","CREBBP","CDCF"," CUX1","DNMT3A","EP300","ETV6","EZH2","FBXW7","FLT3","GATA2","GNAS","IDH1","IDH2","IKZF1","JAK2","KDM5A","KDM6A","KIT","KRAS","MLL","MLL2","MLL3","MLL5","MPL","MYC","NF1","NPM1","NRAS","PHF6","PRPF40B","PTEN","PTPN11","RAD21","RB1","RUNX1","SF1","SF3A1","SF3B1", "SMC1A", "SMC3", "SRSF2","SH2B3","STAG1","STAG2","TET2","TP53","U2AF1","U2AF2","WT1","ZBTB33","ZRSR2"), 
                                       selected= c("BCOR", "RUNX1"),
                                       options = list(maxItems = 2),
                                       multiple=TRUE)
                        
                      ),
                      
                      # Main panel for displaying outputs
                      mainPanel(
                        
                        tabsetPanel(type = "tabs",
                                    tabPanel("Mutation Landscape",
                                             plotOutput("comut_plot",
                                                        height = "750px")
                                    ),
                                    tabPanel("Co-Occurance Plot",
                                             plotOutput("cooccurrance_plot",
                                                        height = "750px")
                                    ),
                                    tabPanel("VAF Scatterplot",
                                             plotOutput("scatter_plot",
                                                        height = "750px")
                                    ),
                                    tabPanel("Survival from VAF Scatterplot",
                                             plotOutput("surv_scatter_plot",
                                                        height = "750px")
                                    ),
                                    tabPanel("3 Gene Co-occuring VAF plot",
                                             plotOutput("three_gene_vaf_plot",
                                                        height = "750px")
                                    ),
                                    tabPanel("Survival Analysis",
                                             plotOutput("survival_plot",
                                                        height = "750px"),
                                             plotOutput("res",
                                                        height = "750px")
                                    )
                        )
                      )
             ),
             tabPanel("About",
                      wellPanel(
                        h2("About metaAML", align="center")
                      ),
                      p("MetaAML is a combination of the largest genomic studdies avaliable in AML. This potal was developed as a recourse to perform simple analysis on genetic interactions and clinical outcomes.")
                      
             ),
             tabPanel("Tutorial",
                      p("insert tutorial here")
             ),
             tabPanel("Data Access",
                      h2("Links to studies include:"),
                      wellPanel(
                        a(h5("Beat AML"), href = "https://www.nature.com/articles/s41586-018-0623-z", target="_blank"),
                        a(h5("AML Multistage"), href = "https://www.nature.com/articles/ng.3756", target="_blank"),
                        a(h5("TCGA AML"), href = "https://www.nejm.org/doi/full/10.1056/NEJMoa1301689", target="_blank")
                        
                      )
             )
  )
))

