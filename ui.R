##########################################
##   Shiny for -omics data presentation ##
##   Allison Barry                      ##
##   University of Oxford               ##
##   allimariebarry@gmail.com           ##
##   for non-commercial use only        ##
##########################################

#load("C:/Users/allim/AppData/Local/Packages/CanonicalGroupLimited.UbuntuonWindows_79rhkp1fndgsc/LocalState/rootfs/home/amb/Database/drg-directory/drg-directory.RData")

#rsconnect::deployApp()
#setRepositories(addURLs = c(BioC = "https://bioconductor.org/packages/3.15/bioc"))
library(cowplot)
library(shiny)
library(data.table)
library(shinythemes)
library(shinydashboard)
library(shinyFeedback)
library(plotly)
library(tidyverse)
library("DESeq2")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library(genefilter)
library("vsn")
library("BiocParallel")
library("AnnotationDbi")
library(ggrepel)
library(stringr)
library(ggplot2)
library(viridis)
library("pacman")
library("leaflet")
library(gridExtra)
library("cowplot")
library("plyr")
library(shinymanager)
library("fontawesome")
library(Seurat)
library("shinyWidgets")


# global variables 
subpopulations = c(
  "Nociceptors 3D" = "TDNV_3D.csv",
  "Nociceptors 4W" = "TDNV_4W.csv",
  "PEP 3D" = "CGRT_3D.csv",
  "PEP 4W" = "CGRT_4W.csv", 
  "NP 3D" = "MRTD_3D.csv",
  "NP 4W" = "MRTD_4W.csv",
  "C-LTMR 3D" = "CRTH_3D.csv",
  "C-LTMR 4W" = "CRTH_4W.csv",
  "Ad- AB-RA LTMRs 3D" = "TBAC_3D.csv",
  "Ad- AB-RA LTMRs 4W" = "TBAC_4W.csv"
)
### Modules
plotdot_ui <- function(id, dataset, combined=FALSE) {
  plotlyOutput(NS(id, "bulkseq_dots"))

}

plotscdot_ui <- function(id) {
  plotlyOutput(NS(id, "scrna_dots"))
}

plotfeature_ui <- function(id) {
  plotOutput(NS(id, "scrna_feature"))
}

plotumap_ui <- function(id) {
  plotlyOutput(NS(id, "scrna_umap"))
}

plothomescdot_ui <- function(id, dataset, combined=FALSE) {
  fluidRow(
    column(3,
           downloadButton(NS(id,"downloadscDot"), "Download")),
    column(10, offset = 1, 
           plotOutput(NS(id, "home_scrna_dots")))
  )
}

plotcombine_ui <- function(id) {
  fluidRow(
    column(12, 
           downloadButton(NS(id,"downloadDot"), "Download"),
           plotOutput(NS(id,"dot"))
           )
  )
}

deg_combine_ui <- function(id) {
  fluidRow(
    column(12, 
           downloadButton(NS(id,"downloadDeg"), "Download"),
           plotOutput(NS(id,"deg")))
  )
}

plotline_ui <- function(id, dataset) {
  box(title = "Injury", status = "primary", 
      plotOutput(NS(id, "bulkseq_lines"))
  )
}

plotsubtype_ui <- function(id, dataset=FALSE) {
  box(width = 12, 
      title = "Subtype Results", status = "primary", 
      plotOutput(NS(id, "bulkseq_lines_subtype"))
  )
}

deg_plot_ui <- function(id) {
  plotlyOutput(NS(id, "deg_plot"))
}

goi_table_ui <- function(id, dataset) {
  DT::dataTableOutput(NS(id,"goi_table"))
}

volcano_plot_ui <- function(id) {
  plotOutput(NS(id, "volcano")
  )
}

contrast_table_ui <- function(id) {
  DT::dataTableOutput(NS(id, "contrast_table"))
}


ui = function(req) {fluidPage(
  
  includeCSS("www/style.css"),
  
  shinyFeedback::useShinyFeedback(),
  
  dashboardPage(
    dashboardHeader(title="DRG Directory", titleWidth = 225
                    ),
    dashboardSidebar(width = 225,
                     sidebarMenu(
                       id = "tabs",
                       
                       HTML(paste0( # oxford logo + ndcn link    
                         "<br><br>",
                         "<a href='https://www.ndcn.ox.ac.uk/research/neural-injury-group' target='_blank'> <img style = 'display: block; margin-left: auto; margin-right: auto;' src='oxfordlogo2.png' width = '105'></a>",
                         "<br><br>")
                       ),
                       
                       selectizeInput(
                         inputId = "geneid", 
                         label = "Search Genes:", 
                         multiple = TRUE,
                         choices = NULL
                       ),
                       
                       fileInput("file", 
                                 label = "Upload File:",
                                 accept = c(
                                 'text/csv',
                                 'text/comma-separated-values',
                                 'text/tab-separated-values',
                                 'text/plain',
                                 '.csv', '.txt',
                                 '.tsv')), 
                       actionButton("reset", "Clear"),
                      
                       selectizeInput(
                         inputId = "sex", 
                         label = "Select Sex:", 
                         choices = c('Both', 'Separate'), 
                         selected = 'Both'),
                       
                       # sidebar menu for tabs (pages)    
                       menuItem("Home", tabName = "tabhome", icon = icon("home")),
                       menuItem("Datasets", tabName = "tabdatasets", icon = icon("database"), 
                                menuItem("Mouse DRG Subtype data", tabName = "tabsubtype"),
                                menuItem("Mouse DRG Bulk data", tabName = "tabmouse"), 
                                menuItem("Rat DRG data", tabName = "tabrat"), 
                                menuItem("Human iPSC data", tabName = "tabhuman"),
                                menuItem("Human Diabetes Skin data", tabName = "tabdb"),
                                menuItem("Human Skin CTS data", tabName = "tabcts"), 
                                menuItem("Spatial Seq data", tabName = "tabspat")
                                ),
                       
                       menuItem("User Guide + Data", tabName = "tabcode", icon = icon("folder-open")),
                       menuItem("Contact", tabName = "tabguide", icon = icon("info-circle")),
                       br(),
                       br()
                     ),
                     
                     #Footer (icons + social media links)
                     HTML(paste0(
                       "<table style='clear: both;padding: 0;text-align: center;vertical-align: middle;line-height: normal;
                margin: 30px;position: fixed;bottom: 0px;width: 165px'>", # start table
                       
                       # icons placed in 6 column table
                       "<tr >",
                       "<td style='padding: 10px;'></td>",
                       "<td style='padding: 5px;'><a href='https://twitter.com/aliibarry' target='_blank'><i class='fab fa-twitter fa-lg'></i></a></td>",
                       "<td style='padding: 5px;'><a href='https://github.com/aliibarry' target='_blank'><i class='fab fa-github fa-lg'></i></a></td>",
                       "<td style='padding: 5px;'><a href='https://orcid.org/0000-0002-6787-6889' target='_blank'><i class='fab fa-orcid fa-lg'></i></a></td>",
                       "<td style='padding: 5px;'><a href='https://www.ndcn.ox.ac.uk/team/allison-barry' target='_blank'><i class='fas fa-brain fa-lg'></i></a></td>",
                       "<td style='padding: 10px;'></td>",
                       "</tr>",
                       
                       # second table row, merged columns    
                       "<tr>",
                       "<script>","var today = new Date();","var yyyy = today.getFullYear();","</script>",
                       "<td  colspan='6', style = 'text-align: center;'><small>&copy; - <a href='https://github.com/aliibarry?tab=projects' target='_blank'>github.com/aliibarry</a> - <script>document.write(yyyy);</script></small></td>",
                       "</tr>",
                       
                       "</table>", #end table
                       "<br>")
                     )
    ),
    
    # Main panels for output, each TabItem is a different page, same sidebar menu
    dashboardBody(
      tabItems(
        tabItem(tabName="tabhome", 
                #h4("Home"),
                fluidRow(
                  tableOutput("contents"),
                  # text summary for page
                  box(width=12,
                      status = "primary", 
                      solidHeader = TRUE,
                      #height = 275,
                      title = "Overview",
                      p("This database provides an interface to explore murine -omics datasets.
                        Full descriptions and project-specific analyses are highlighed in the 'Datasets' tab.
                        All results are plotted as median transcripts per million (TPM).
                        Search result data for each dataset is available for download
                        on their respective pages (see: Datasets). All raw data is available
                        on GEO. Hypothesis testing with FDR < 0.05 and log2 fold
                        change (LFC) > 1.
                        ")
                  )
                ), #fluidrow
                br(),
                fluidRow(
                  box(status = "primary",
                      width = 12, 
                      title = "Dataset Overview",
                      column(width = 10, 
                             includeHTML("/Users/lynnezhao/Desktop/table.html")
                             ),
                      column(width = 2,
                             br(),br(),
                             box(status="primary", width = 12,
                               actionLink("link_to_subtype", "Browse Subtype Data")
                             )
                             ),
                      column(width = 2,
                             box( status = "primary", width = 12,
                                  actionLink("link_to_mouse", "Browse Mouse Bulk Data")
                             )),
                      column(width = 2,
                             box(status = "primary", width = 12,
                                  actionLink("link_to_rat", "Browse Rat DRG Data")
                             )
                             ),
                      column(width = 2,
                             box(status = "primary", width = 12,
                                 actionLink("link_to_human", "Browse iPSC Data")
                             )),
                      column(width = 2,
                             box(status = "primary", width = 12,
                                 actionLink("link_to_db", "Browse Diabetes Data")
                             )),
                      column(width = 2,
                             box(status = "primary", width = 12,
                                 actionLink("link_to_cts", "Browse CTS Data")
                             )),
                      column(width = 2,
                             box(status = "primary", width = 12,
                                 actionLink("link_to_zheng", "Browse Zheng Data")
                             )),
                      column(width = 2,
                             box(status = "primary", width = 12,
                                 actionLink("link_to_scrna", "Browse Spatial Seq Data")
                             )),
                      collapsible = TRUE),
                 
                ),
                fluidRow(
                  br(),
                  column(3,offset = 0, 
                         actionButton("load", "Plot Graphs"), br()
                         ), 
                  br()
                ),
                
                fluidRow(
                  br(),
                  column(width = 12,
                         p("All results are plotted as median vst-transformed count data. 
                                     Search result data is available for download in in the 'Tables' tab. 
                                     Interactive queries for hypothesis testing with FDR and log2 fold 
                                     change (LFC) details are linked below.")
                  ),
                  box(width = 12,
                      title = "Naive",
                      collapsible = TRUE,
                      status = "primary", 
                      plotcombine_ui("dot"), 
                      plothomescdot_ui("homespat")
                  ),
                  box(width = 12, 
                      title = "Differential Gene Analysis",
                      collapsible = TRUE,
                      status = "primary",
                      deg_combine_ui("deg_plot"))
                )
                
                
        ), # tabItem HOME
      
        tabItem(tabName="tabsubtype", 
                fluidRow(
                  box(title = "DRG subpopulation RNAseq",
                      width = 12,
                      status = "primary", 
                      solidHeader = TRUE,
                      includeMarkdown("datasetsummary.md"),
                      img(src = "schematic.png", height = 150, width = 400), 
                  )
                ),
                
                br(),
                actionLink("link_to_home", "Back to Home", icon = icon("home")),
                 
                br(),br(),
                actionButton("load2", "Plot Graphs"),
                br(),br(),
                fluidRow(
                  box(title = "Naive",
                      status = "primary", 
                         plotdot_ui("dot_subtype", FALSE)
                  ), 
                  plotline_ui("line", "subtype"),
                  plotsubtype_ui("lines_subtype", "subtype"), 
                  box(title = "Differential Gene Analysis",
                      status = "primary", 
                         deg_plot_ui("deg_plot_subtype")), 
                  box(title = "Ipsilateral vs Contralateral",
                      status = "primary", 
                         actionButton("volc", "Plot Volcano Graphs"),
                         selectInput("volca", "", 
                                     choices = c(
                                       "Nociceptors 3D",
                                       "Nociceptors 4W",
                                       "PEP 3D",
                                       "PEP 4W", 
                                       "NP 3D",
                                       "NP 4W",
                                       "C-LTMR 3D",
                                       "C-LTMR 4W",
                                       "Ad- AB-RA LTMRs 3D",
                                       "Ad- AB-RA LTMRs 4W"),
                                     selected = ""),
                         volcano_plot_ui("volcano")
                  )

                ), # line plots, line subtype, vocano plots 
                fluidRow(
                  box(
                    width = 12,
                    title = "Result Table", 
                    status = "primary", 
                    includeMarkdown("datatable_notes.Rmd"),
                    downloadButton("downloadData", "Download"),
                    br(),
                    goi_table_ui("goi_table", "subtype")
                  )
                ),
                fluidRow(
                  box(
                    width = 12,
                    title = "Differential Analysis Table", 
                    status = "primary", 
                    selectInput("contrast", "", 
                                choices = c(
                                  "Nociceptors 3D",
                                  "Nociceptors 4W",
                                  "PEP 3D",
                                  "PEP 4W", 
                                  "NP 3D",
                                  "NP 4W",
                                  "C-LTMR 3D",
                                  "C-LTMR 4W",
                                  "Ad- AB-RA LTMRs 3D",
                                  "Ad- AB-RA LTMRs 4W"),
                                selected = ""), 
                    contrast_table_ui("contrast_table")
                  )
                )
          ), 
        
        # for mouse data page 
        tabItem(tabName = "tabmouse", 
                actionButton("load3", "Plot Graphs"), # action button
                actionLink("link_to_home2", "Back to Home", icon = icon("home")),
                br(),br(),
                fluidRow(
                  box(title = "Naive", status = "primary", 
                      plotdot_ui("dot_mouse", FALSE)
                  ),
                  
                  plotline_ui("mouse_line", "mouse"),
                  plotsubtype_ui("mouse_lines_subtype", "mouse"), 
                  box (title = "Differential Gene Analysis", status = "primary", 
                        deg_plot_ui("deg_plot_mouse")), 
                  box(
                    title = "SHAM vs SNI", status = "primary", 
                    actionButton("volcmouse", "Plot Volcano Graphs"),
                    selectInput("volcamouse", "", 
                                choices = c("B10D2", "BALB"),
                                selected = ""),
                    volcano_plot_ui("volcano_mouse")
                  )
                ),
                fluidRow(
                  box(
                    width = 12,
                    title = "Result Table", 
                    status = "primary",
                    goi_table_ui("goi_table_mouse")
                  )
                ), 
                fluidRow(
                  box(width = 12, 
                      title = "Differential Analysis Table",
                      status = "primary",
                      selectInput("contrastm", "", 
                                 choices = c("B10D2", "BALB"),
                                 selected = ""), 
                      contrast_table_ui("mouse_contrast_table")
                      
                  )
                )
                
        ), 
        
        tabItem(tabName = "tabrat", 
                actionButton("load4", "Plot Graphs"), # action button 
                actionLink("link_to_home3", "Back to Home", icon = icon("home")),
                fluidRow(
                  box(title = "Naive", status = "primary", 
                      plotdot_ui("dot_rat", FALSE)
                  ), 
                  plotline_ui("rat_line", "rat")), 
                fluidRow(
                  box(width = 6, 
                      title = "Differential Gene Analysis",
                      status = "primary",
                      deg_plot_ui("deg_plot_rat")), 
                  box(width=6,
                    title = "SHAM vs SNT", status = "primary", 
                    actionButton("volcrat", "Plot Volcano Graphs"),
                    volcano_plot_ui("volcano_rat")
                  )
                ), 
                fluidRow(
                  box(
                    width = 12,
                    title = "Result Table", 
                    status = "primary",
                    goi_table_ui("goi_table_rat")
                  )
                ), 
                fluidRow(
                  box(width = 12, 
                      title = "Differential Analysis Table",
                      status = "primary",
                      DT::dataTableOutput("contrast_table_rat")
                  )
                )
        ),
        
        tabItem(tabName = "tabhuman", 
                fluidRow(
                  box(title = "iPSC RNAseq",
                      status = "primary", 
                      width=12,
                      solidHeader = TRUE,
                      column(width=9,
                             # includeMarkdown("ipscsummary.md")
                             ),
                      column(width=3,
                             img(src = "ipsc.jpeg", height = 320, width = 320))   
                  )
                ),
                br(),
                actionButton("load5", "Plot Graphs"), # action button 
                actionLink("link_to_home4", "Back to Home", icon = icon("home")),
                br(),br(),
                fluidRow(
                  # can plot by cell-line or differentiation state, choose it 
                  box(
                    status="primary",
                    title = "Naive",
                    width = 6,
                    tabsetPanel(
                      tabPanel(
                        title = "By Differentiation",
                        plotdot_ui("dot_human", FALSE)
                      ),
                      tabPanel(
                        title = "By Cell Line",
                        plotdot_ui("dot_human2", FALSE)
                      )  
                    )
                  ),
                  plotline_ui("human_line", "human")
                  ),
                fluidRow(
                  box(width = 6, 
                      title = "Differential Gene Analysis",
                      status = "primary",
                      deg_plot_ui("deg_plot_human")),
                  box(
                    title = "Healthy vs HSN1", status = "primary", 
                    actionButton("volch", "Plot Volcano Graphs"),
                    selectInput("volcahuman", "", 
                                choices = c("iPSCDN_young", "iPSCDN_old"),
                                selected = ""),
                    volcano_plot_ui("volcano_human")
                  )
                ), 
                fluidRow(
                  box(
                    width = 12,
                    title = "Result Table", 
                    status = "primary",
                    goi_table_ui("goi_table_human")
                  )
                ),
                fluidRow(
                  box(
                    width = 12,
                    title = "Differential Analysis Table", 
                    status = "primary", 
                    selectInput("contrasth", "", 
                                choices = c("iPSCDN_young", "iPSCDN_old"),
                                selected = ""), 
                    contrast_table_ui("contrast_table_human")
                  )
                )
        ),
        
        tabItem(tabName = "tabdb",
                fluidRow(
                  box(title = "Human Diabetes RNAseq",
                      status = "primary", 
                      width=12,
                      solidHeader = TRUE,
                      # includeMarkdown("ipscsummary.md")
                  )
                ),
                br(),
                actionButton("load6", "Plot Graphs"), # action button 
                actionLink("link_to_home5", "Back to Home", icon = icon("home")),
                br(),br(),
                fluidRow(
                  box(title = "Naive", status = "primary", 
                      plotdot_ui("dot_db", FALSE)
                  ), 
                  plotline_ui("db_line", "skin")), 
                fluidRow(
                  box(width = 6, 
                      title = "Differential Gene Analysis",
                      status = "primary",
                      deg_plot_ui("deg_plot_db")), 
                  box(width = 6,
                      title = "Painful vs Painless", status = "primary", 
                      actionButton("volcd", "Plot Volcano Graphs"),
                      selectInput("volcadb", "", 
                                  choices = c("Diabetes_skin", "Diabetes_skin_females", 
                                              "Diabetes_skin_males"),
                                  selected = ""),
                      volcano_plot_ui("volcano_db")
                  )
                ), 
                 
                fluidRow(
                  box(
                    width = 12,
                    title = "Differential Analysis Table", 
                    status = "primary", 
                    selectInput("contrastd", "", 
                                choices = c("Diabetes_skin", "Diabetes_skin_females", 
                                            "Diabetes_skin_males"),
                                selected = ""), 
                    contrast_table_ui("contrast_table_db")
                  )
                )
        ),
        
        tabItem(tabName = "tabcts", 
                fluidRow(
                  box(title = "Human Skin CTS RNAseq",
                      status = "primary", 
                      width=12,
                      solidHeader = TRUE,
                      # includeMarkdown("ipscsummary.md")
                  )
                ),
                br(),
                actionButton("load7", "Plot Graphs"), # action button 
                actionLink("link_to_home6", "Back to Home", icon = icon("home")),
                br(),br(),
                fluidRow(
                  box(title = "Naive", status = "primary", 
                      plotdot_ui("dot_cts", FALSE)
                  ), 
                  plotline_ui("cts_line", "skin")), 
                fluidRow(
                  box(width = 6, 
                      title = "Differential Gene Analysis",
                      status = "primary",
                      deg_plot_ui("deg_plot_cts")), 
                  box(width = 6,
                      title = "Painful vs Painless", status = "primary", 
                      actionButton("volcc", "Plot Volcano Graphs"),
                      volcano_plot_ui("volcano_cts")
                  )
                ), 
                fluidRow(
                  box(width = 12, 
                      title = "Differential Analysis Table",
                      status = "primary",
                      DT::dataTableOutput("contrast_table_cts")
                  )
                )
        ),
        
        tabItem(tabName = "tabspat", 
                actionButton("load8", "Plot Graphs"), 
                actionLink("link_to_home7", "Back to Home", icon = icon("home")),
                fluidRow(
                  box(width=6, 
                    title = "Naive", status = "primary", 
                      plotscdot_ui("dotspat")
                  ),
                  box(width=6,
                      title = "UMAP", status = "primary", 
                      plotumap_ui("umapspat")
                  )
                ), 
               
                fluidRow(
                  box(width=6,
                      height = 550,
                      title = "Feature", status = "primary", 
                      column(width=2,br(),
                             actionButton("load9", "Plot Graphs")),
                      column(width=10,
                             selectizeInput(
                               inputId = "scgeneid", 
                               label = "Search Genes:", 
                               multiple = FALSE,
                               choices = NULL
                             )),
                      plotfeature_ui("featurespat")
                  )
                )
                
        ),
        
        ## Supply simple links for each paper + supplementary repository
        tabItem(tabName="tabcode",
                h4("Data Access"),
                includeMarkdown('codedata.md'),
                br(),
                
                # HTML links to funders, embedded in a table for alignment
                HTML(paste0(
                  "<table style='margin-left:auto; margin-right:auto; bottom=0;'>",
                  
                  "<tr>",
                  "<td><a href='https://www.gtc.ox.ac.uk/' target='_blank'> <img style= 'display: center;' src='gtclogo.png', width = '100'></a></td>",
                  "<td><a href='https://wellcome.org/' target='_blank'> <img style = 'display: center;' src='wt-logo.svg'. width = '100'></a></td>",
                  "</tr>",
                  
                  "</table>"
                )
                )
                
        ), #tabitem code
        
        tabItem(tabName="tabguide",
                h4("Contact"),
                includeMarkdown('userguide.md'),
                
                # Lab location map, in a column solely for aesthetic
                column(10, offset = 1,
                       leafletOutput('myMap', width = "100%", height = 350)
                )
                
        ) # tabItem help guide
      ) # tabItems (all) 
    )  #mainpanel dashboardbody close
  ) # dashboad page close
) #fluid page
}#ui

# ui = secure_app(ui)
