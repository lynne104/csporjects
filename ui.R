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
# library(Seurat)
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
plotdot_ui <- function(id, dataset) {
  plotlyOutput(NS(id, "bulkseq_dots"))
}

plotline_ui <- function(id, dataset) {
  box(title = "Injury", status = "primary", 
         if (dataset == "subtype"){
           plotlyOutput(NS(id, "bulkseq_lines"))
         }, 
         if (dataset == "mouse") {
           plotlyOutput(NS(id, "mouse_lines"))
         }
  )
}

plotsubtype_ui <- function(id, dataset) {
  box(width = 12, 
      title = "Subtype Results", status = "primary", 
         if (dataset == "subtype"){
           plotlyOutput(NS(id, "bulkseq_lines_subtype"))
         }, 
         if (dataset == "mouse") {
           plotlyOutput(NS(id, "mouse_lines_subtype"))
         }
  )
}

deg_plot_ui <- function(id) {
  plotlyOutput(NS(id, "deg_plot"))
}

goi_table_ui <- function(id, dataset) {
  DT::dataTableOutput(NS(id,"goi_table"))
}

volcano_plot_ui <- function(id) {
  plotOutput(NS(id, "volcano"), 
  )
}

contrast_table_ui <- function(id) {
  DT::dataTableOutput(NS(id, "contrast_table"))
}

### IU 
shinyUI(fluidPage(
  
  #CSS style sheet
  includeCSS("www/style.css"),
  
  shinyFeedback::useShinyFeedback(),
  
  # UI built with shiny dashboard, requires boxes/columns in tabItems 
  dashboardPage(
    dashboardHeader(title="DRG Directory", titleWidth = 225),
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
                       selectizeInput(
                         inputId = "sex", 
                         label = "Select Sex:", 
                         choices = c('Both', 'Separate'), 
                         selected = 'Both'),
                       
                       # sidebar menu for tabs (pages)    
                       menuItem("Home", tabName = "tabhome", icon = icon("home")),
                       menuItem("Datasets", tabName = "tabdatasets", icon = icon("database"), 
                                menuItem("DRG subpopulation RNAseq", tabName = "tabsubtype", icon=icon("dna")),
                                menuItem("Mouse data", tabName = "tabmouse", icon=icon("frog")), 
                                menuItem("Rat data", tabName = "tabrat", icon=icon("frog"))
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
                  # text summary for page
                  box(width=12,
                      status = "primary", 
                      solidHeader = TRUE,
                      #height = 275,
                      title = "Overview",
                      p("This database provides an interface to explore murine -omics datasets. 
                                  Deep RNA-seq of mouse DRG subpopulations after spare nerve injury (SNI) was performed 
                                  to interrogate subtype-specific and shared injury signatures. 
                                  Specific transgenic details and methodologies will be available in the published reports, 
                                  and are currently available by request. 
                                  All work presented here is currently unpublished, and is the work of Ali Barry, 
                                  Giorgos Baskozos, and Dave Bennett at the University of Oxford (NDCN).
                                Full descriptions and relevant citations can be found in the 'Dataset' tab.  
                                ")
                  )
                ), #fluidrow
                br(),
                fluidRow(
                  box(status = "primary",
                    width = 4, actionLink("link_to_subtype", "Browse DRG subpopulation RNAseq Data", icon=icon("dna"))
                  ),
                  box(status = "primary",
                      width = 4, actionLink("link_to_mouse", "Browse Mouse Data")),
                  box(status = "primary",
                      width = 4, actionLink("link_to_rat", "Browse Rat Data"))
                ),
                
                fluidRow(
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
                  box(width = 10,
                      title = "Naive",
                      status = "primary",
                      plotdot_ui("dot")
                  ),
                  box(width = 10, 
                      title = "Differential Gene Analysis",
                      status = "primary",
                      deg_plot_ui("deg_plot"))
                )
                
                
        ), # tabItem HOME
      
        tabItem(tabName="tabsubtype", 
                fluidRow(
                  box(title = "DRG subpopulation RNAseq",
                      width = 12,
                      status = "primary", 
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
                         plotdot_ui("dot_subtype")
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
                                     choices = subpopulations,
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
                    goi_table_ui("goi_table", "subtype")
                  )
                ),
                fluidRow(
                  box(
                    width = 12,
                    title = "Differential Analysis Table", 
                    status = "primary", 
                    selectInput("contrast", "", 
                                choices = subpopulations,
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
                      plotdot_ui("dot_mouse")
                  ),
                  
                  plotline_ui("mouse_line", "mouse"),
                  plotsubtype_ui("mouse_lines_subtype", "mouse"), 
                  box (title = "Differential Gene Analysis", status = "primary", 
                        deg_plot_ui("deg_plot_mouse")), 
                  box(
                    title = "SHAM vs SNI", status = "primary", 
                    actionButton("volcmouse", "Plot Volcano Graphs"),
                    selectInput("volcamouse", "", 
                                choices = c("B10D2_Mouse_SNI_vs_SHAM.csv", "BALBc_Mouse_SNI_vs_SHAM.csv"),
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
                                 choices = c("B10D2_Mouse_SNI_vs_SHAM.csv", "BALBc_Mouse_SNI_vs_SHAM.csv"),
                                 selected = ""), 
                      contrast_table_ui("mouse_contrast_table")
                  )
                )
                
        ), 
        
        tabItem(tabName = "tabrat", 
                actionButton("load4", "Plot Graphs"), # action button 
                plotdot_ui("dot_rat")
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
)#ui