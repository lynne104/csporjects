##########################################
##   Shiny for -omics data presentation ##
##   Allison Barry                      ##
##   University of Oxford               ##
##   allimariebarry@gmail.com           ##
##   for non-commercial use only        ##
##########################################
library(profvis)
# load data  
data_dir = "/Users/lynnezhao/Desktop/drg/drg-directory.RData"
load(data_dir)

data_dir = "/Users/lynnezhao/Desktop/drg/datamining.RData"
load(data_dir)

# data_dir = "/Users/lynnezhao/Desktop/drg/ratsnt.RData"
# load(data_dir)

TABLE_PATH = "/Users/lynnezhao/Desktop/data/" 

# path for volcano plots dfs 
PATH = "/Users/lynnezhao/Desktop/drg/"

theme_line = theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(family = "Gadugi", size=12),
        axis.text.x = element_text(family = "Gadugi", size=10, angle = 45, hjust= 1),
        axis.text.y = element_text(family = "Gadugi", size=10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(), legend.justification = c(0,0.3)) 

th = theme_bw() + theme(panel.grid = element_blank(),
                               axis.title = element_blank(),
                               axis.text.x = element_text(size=8, angle = 45, hjust= 1),
                               axis.text.y = element_text(size=10),
                               axis.ticks.x = element_blank(),
                               axis.ticks.y = element_blank(), legend.justification = c(0,0.3))

th_vol = theme_bw() + theme(panel.grid = element_blank(),
                            axis.ticks.x = element_blank(),
                            axis.ticks.y = element_blank(), legend.justification = c(0,0.3))

population_labels = c("TBAC" = "A\u03b4-LTMR + A\u03b2-RA-LTMR",
                      "CRTH" = "C-LTMR", "MRTD" = "NP",
                      "CGRT" = "PEP", "TDNV" = "Nociceptors", 
                      "b10b2" = "b10b2", "balb" = "balb")
sexlabels = c("F" = "Female", "M" = "Male")
subpopulations = c(
  "Nociceptors 3D" = "TDNV_3D.csv", "Nociceptors 4W" = "TDNV_4W.csv",
  "PEP 3D" = "CGRT_3D.csv", "PEP 4W" = "CGRT_4W.csv", 
  "NP 3D" = "MRTD_3D.csv","NP 4W" = "MRTD_4W.csv",
  "C-LTMR 3D" = "CRTH_3D.csv","C-LTMR 4W" = "CRTH_4W.csv",
  "Ad- AB-RA LTMRs 3D" = "TBAC_3D.csv","Ad- AB-RA LTMRs 4W" = "TBAC_4W.csv"
)

subpopulation_labels = c(
  "TDNV_3D" = "Nociceptors 3D", "TDNV_4W" = "Nociceptors 4W",
  "CGRT_3D" = "PEP 3D","CGRT_4W" = "PEP 4W", 
  "MRTD_3D" = "NP 3D", "MRTD_4W" = "NP 4W",
  "CRTH_3D" = "C-LTMR 3D","CRTH_4W" = "C-LTMR 4W",
  "TBAC_3D" = "Ad- AB-RA LTMRs 3D", "TBAC_4W" = "Ad- AB-RA LTMRs 4W"
)

df <- function(data, type) {
  matfilt <- data()[,1:154]
  tcounts <- t(matfilt) %>%
    base::merge(bulkseq_colData, ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-nrow(matfilt)+1):(ncol(.)))
  tcounts$symbol <- gene_data[tcounts$gene,]$mgi_symbol # add gene symbols for goi
  return(tcounts)
}

prep <- function(df) {
  tcounts = df()
  tcounts_med <- tcounts %>% dplyr::group_by(Condition, symbol) %>% 
    dplyr::summarise(expression=median(expression))
  df_rat <- tcounts_med[!tcounts_med$Condition %in% "SNT", ]
  return(df_rat)
}
# data preprocessing; return median expression of a given gene list 
preprocess <- function(df, rownum, paintype, colData, sex, dataset, graphtype) {
  matfilt <- df()[,1:rownum]
  tcounts <- t(matfilt) %>%
    base::merge(colData, ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-nrow(matfilt)+1):(ncol(.)))
  tcounts$symbol <- gene_data[tcounts$gene,]$mgi_symbol # add gene symbols for goi
  if (graphtype == "dot") {
    if (sex == "Both") {
      tcounts_med <- tcounts %>% dplyr::group_by(Condition, Population, symbol) %>% 
        dplyr::summarise(expression=median(expression))
    }
    if (sex == "Separate") {
      tcounts_med <- tcounts %>% dplyr::group_by(Condition, Population, symbol, Sex) %>% 
        dplyr::summarise(expression=median(expression))
    }
    tcounts_med <- tcounts_med[!tcounts_med$Condition %in% paintype, ]
  }
  if (graphtype == "line") {
    if (dataset == "subtype") {
      if (sex == "Both") {
        tcounts_med <- tcounts %>% dplyr::group_by(Condition, symbol, Timepoint) %>% 
          dplyr::summarise(expression=median(expression))
      }
      if (sex == "Separate") {
        tcounts_med <- tcounts %>% dplyr::group_by(Condition, symbol, Sex, Timepoint) %>% 
          dplyr::summarise(expression=median(expression))
      }
    }
    if (dataset == "mouse") {
      if (sex == "Both") {
        tcounts_med <- tcounts %>% dplyr::group_by(Condition, symbol) %>% 
          dplyr::summarise(expression=median(expression))
      }
      if (sex == "Separate") {
        tcounts_med <- tcounts %>% dplyr::group_by(Condition, symbol, Sex) %>% 
          dplyr::summarise(expression=median(expression))
      }
    }
  }

  if (graphtype == "line_subtype") {
    if (dataset == "subtype") {
      if (sex == "Both") {
        tcounts_med <- tcounts %>% dplyr::group_by(Condition, symbol, Timepoint, Population) %>% 
          dplyr::summarise(expression=median(expression))
      }
      if (sex == "Separate") {
        tcounts_med <- tcounts %>% dplyr::group_by(Condition, symbol, Sex, Timepoint, Population) %>% 
          dplyr::summarise(expression=median(expression))
      }
    }
    if (dataset == "mouse") {
      if (sex == "Both") {
        tcounts_med <- tcounts %>% dplyr::group_by(Condition, symbol, Population) %>% 
          dplyr::summarise(expression=median(expression))
      }
      if (sex == "Separate") {
        tcounts_med <- tcounts %>% dplyr::group_by(Condition, symbol, Sex, Population) %>% 
          dplyr::summarise(expression=median(expression))
      }
      }
   
  }
  tcounts_med$Dataset = rep(dataset, nrow(tcounts_med))
  return(tcounts_med)
  
}

# plot line plots for subtype data and mouse data   
plotline_server <- function(id, df, sex, dataset) {
  moduleServer(id, function(input, output, session) {
    tcounts_med = df
    if (dataset == "subtype") {
      g = ggplot(data=tcounts_med, aes(x=Condition, y=expression, group=interaction(symbol, Timepoint))) + 
        scale_colour_viridis(discrete=TRUE, end = .80) + 
        geom_line(aes(color=symbol, linetype=Timepoint)) + geom_point(aes(col=symbol)) + theme_line + ylab("Expression") +
        guides(col=FALSE, linetype=guide_legend(label.hjust=0.5, ncol=8))
      if (sex == 'Both') {
        output$bulkseq_lines <- renderPlotly({g})
      }
      if (sex =='Separate'){
        output$bulkseq_lines <- renderPlotly({
          g + facet_wrap(~Sex, ncol=2, labeller = labeller(Sex = sexlabels)) + labs(col="") 
        })
      }
    }
    if (dataset == "mouse") {
      g = ggplot(data=tcounts_med, aes(x=Condition, y=expression, group=symbol)) + 
        scale_colour_viridis(discrete=TRUE, end = .80) + 
        geom_line(aes(col=symbol)) + geom_point(aes(col=symbol)) + theme_line + ylab("Expression") +
        guides(col=FALSE, linetype=guide_legend(label.hjust=0.5, ncol=8))
      if (sex == "Both") {
        output$mouse_lines <- renderPlotly({g})
      }
      else {
        output$mouse_lines <- renderPlotly({g + facet_wrap(~Sex, ncol=2, labeller = labeller(Sex = sexlabels)) + labs(col="")})
      }
    }
   
  })
}


# ploting dot plots for all data 
plotdot_server <- function(id, df, sex) {
  moduleServer(id, function(input, output, session) {
    tcounts_med = df
    g = ggplot(data = tcounts_med, aes(x= Population, y = symbol)) + 
      scale_colour_viridis_c(option = "magma", end = .90) + 
      geom_point(aes(col=expression, size=expression)) + th + 
      facet_wrap(~Dataset, scales = 'free_x') +
      scale_x_discrete(labels=population_labels)
    if (sex == 'Both') {
      output$bulkseq_dots <- renderPlotly({g})
    }
    if (sex == "Separate") {
      output$bulkseq_dots <- renderPlotly({ 
        g + facet_wrap(~Sex, ncol=2, labeller = labeller(Sex = sexlabels))
      })
    }
  })
}

# ploting subtype plots; 
plotsubtype_server <- function(id, df, sex, dataset) {
  moduleServer(id, function(input, output, session) {
    tcounts_med = df
    if (dataset == "subtype"){
      g =  ggplot(data=tcounts_med, aes(x=Condition, y=expression, group=interaction(symbol, Timepoint))) +  
        scale_colour_viridis(discrete=TRUE, end = .80) + 
        geom_line(aes(col=symbol, linetype=Timepoint)) + geom_point(aes(col=symbol)) + 
        # facet_grid(symbol~Population, labeller=labeller(Population=population_labels)) +
        facet_wrap(~Population, ncol=5, labeller=labeller(Population=population_labels))+ 
        th + ylab("Expression") + labs(col="")
      if (sex == "Both") {
        output$bulkseq_lines_subtype <- renderPlotly({g})
      }
      else {
        output$bulkseq_lines_subtype <- renderPlotly({g + 
            facet_grid(Sex~Population, labeller=labeller(Population=population_labels, Sex=sexlabels)) + labs(col="")
          })
      }
    }
    if (dataset == "mouse"){
      g =  ggplot(data=tcounts_med, aes(x=Condition, y=expression, group=symbol)) +  
        scale_colour_viridis(discrete=TRUE, end = .80) + 
        geom_line(aes(col=symbol)) + geom_point(aes(col=symbol)) + 
        facet_wrap(~Population)+ th + ylab("Expression") + labs(col="")
      if (sex == "Both") {
        output$mouse_lines_subtype <- renderPlotly({g})
      }
      else {
        output$mouse_lines_subtype <- renderPlotly({g + 
            facet_grid(Sex~Population, labeller=labeller(Sex=sexlabels)) + labs(col="")
        })
      }
    }
  })
}

# a server for rendering contrast tables 
contrast_table_server <- function(id, df) {
  moduleServer(id,function(input,output,session) {
    output$contrast_table <- DT::renderDataTable({
      req(df())
      DT::datatable(
        df(),
        width = 12,
        class = 'nowrap',
        options = list(scrollX = TRUE, pageLength = 10)
      )
    })
  })
}

# a server for rendering goi tables 
goi_table_server <- function(id, df, dataset) {
  moduleServer(id, function(input, output, session){
    output$goi_table <- DT::renderDataTable({
      req(df())
      datatable <- df()
      datatable$geneID <- rownames(datatable)
      DT::datatable(
        datatable,rownames=datatable$symbol,
        style="default",
        class = 'nowrap',
        options = list(scrollX = TRUE, pageLength = 5)
      )
    })
  })
}

# plot deg servers
deg_plot_server <- function(id, df) {
  moduleServer(id, function(input, output, session){
    output$deg_plot <- renderPlotly({
      ggplot(data = df(), aes(x=interaction(Population), y = symbol, 
                              text = paste('padj: ',padj))) + 
        # facet_grid(species~dataset, scale="free", space="free") +
        scale_colour_viridis_c(option = "viridis", end = .90) +
        geom_point(aes(col=padj, shape=sig, size=0.5)) + scale_shape_manual(values=c(1, 19)) +
        th + labs(shape = "", size = "") + facet_wrap(~Dataset, scales = 'free_x')  
    })
    
  })
}

volcano_plot_server <- function(id, geneids, df) {
  moduleServer(id, function(input, output, session){
    output$volcano <- renderPlot({
      il_genes <- df %>% filter(symbol %in% geneids) 
      ggplot(data = df, aes(x = log2FoldChange, y = log10fdr)) + 
        geom_point(colour = "grey", alpha = 0.5) +
        geom_point(data = il_genes, # New layer containing data subset il_genes       
                   size = 2,
                   shape = 21,
                   fill = "firebrick",
                   colour = "black") + 
        geom_hline(yintercept = -log10(0.05),
                   linetype = "dashed") + 
        geom_vline(xintercept = c(log2(0.5), log2(2)),
                   linetype = "dashed") +
        geom_label_repel(data = il_genes,   
                         aes(label = symbol),
                         force = 2,
                         nudge_y = 1) + 
        scale_colour_manual(values = cols) + 
        scale_x_continuous(breaks = c(seq(-10, 10, 2)),     
                           limits = c(-10, 10)) + th_vol + 
        labs(y="-log10(FDR)")
    })
  })
} 

################################################ SERVER ##############################################################################
shinyServer(function(input, output, session) {
  load(data_dir)
  
  # select genes 
  updateSelectizeInput(session, 
                       inputId = "geneid", 
                       label = "Search Genes:",
                       choices = bulkseq_mat[,155], 
                       server = TRUE,
                       selected = "Atf3"
                       #options = list(placeholder = 'select a gene', 'plugins' = list('remove_button'))
  ) 
  
  
  ### linking to other tabs
  observeEvent(input$link_to_subtype, {
    newvalue <- "tabsubtype"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$link_to_mouse, {
    newvalue <- "tabmouse"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$link_to_rat, {
    newvalue <- "tabrat"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$link_to_home, {
    newvalue <- "tabhome"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$link_to_home2, {
    newvalue <- "tabhome"
    updateTabItems(session, "tabs", newvalue)
  })
  

  # data preprocessing for deg plots 
  b10d2 = b10d2[c("log2FoldChange", "padj", "symbol")]
  b10d2$Population = rep("b10d2", nrow(b10d2))
  b10d2$Dataset = rep("mouse", nrow(b10d2))
  
  balb = balb[c("log2FoldChange", "padj", "symbol")]
  balb$Population = rep("balb", nrow(balb))
  balb$Dataset = rep("mouse", nrow(balb))
  
  mouse_deg_df = rbind(b10d2, balb)
  mouse_deg_df = mutate(mouse_deg_df, sig=ifelse(mouse_deg_df$padj<0.05, "SIG", "NS"))
  
  # for subtype dataset 
  PATH = "/Users/lynnezhao/Desktop/drg/"
  subtype_deg_df = data.frame()
  for (pop in subpopulations) {
    res <- fread(paste0(PATH, pop))
    colnames(res) <- c('symbol', 'log2FoldChange', 'padj', 'ID')
    rownames(res) <- res$ID
    res$Population <- rep(substring(pop, 1, 7), nrow(res)) # add the population label
    res$Dataset = rep("subtype", nrow(res))
    subtype_deg_df = bind_rows(subtype_deg_df, res)
  }
  subtype_deg_df = mutate(subtype_deg_df, sig=ifelse(subtype_deg_df$padj<0.05, "SIG", "NS"))
  subtype_deg_df = subtype_deg_df[c("log2FoldChange", "padj", "symbol", "Population","Dataset", "sig")]
  combined_deg_df = rbind(mouse_deg_df, subtype_deg_df)
  
  # plot homepage plots 
  observeEvent(input$load, {
    data1 <- reactive({
      TPM_mouse[TPM_mouse[,21] %in% input$geneid,]
    }) 
    data2 <- reactive({
      bulkseq_mat[bulkseq_mat[,155] %in% input$geneid,]
    })
    data3 <- reactive({
      counts_rat[counts_rat[,6] %in% input$geneid,]
    })

    r = rbind(preprocess(reactive({data1()}), 20, "SNI", TPM_mouse_colData, input$sex, "mouse", "dot"), 
              preprocess(reactive({data2()}), 154, "ipsi", bulkseq_colData, input$sex,"subtype", "dot"))
    plotdot_server("dot",r, input$sex) #dotplot
 
    # ploting differential gene analysis 
    deg <- reactive({
      df <- combined_deg_df %>% filter(symbol %in% input$geneid)
      return(df)
    }) %>% bindCache(input$geneid)
    
    deg_plot_server("deg_plot", reactive({deg()}))
  })
  

  # plots for subtype dataset page   
  observeEvent(input$load2, {
    data <- reactive({
      bulkseq_mat[bulkseq_mat[,155] %in% input$geneid,]
    }) %>% bindCache(input$geneid)
    plotdot_server("dot_subtype", preprocess(reactive({data()}), 154, "ipsi", bulkseq_colData, input$sex,"subtype", "dot"), input$sex)
    plotline_server("line", preprocess(reactive({data()}), 154, "ipsi", bulkseq_colData, input$sex,"subtype", "line"), input$sex, "subtype") # lineplot
    plotsubtype_server("lines_subtype", preprocess(reactive({data()}), 154, "ipsi", bulkseq_colData, input$sex,"subtype", "line_subtype"), input$sex, "subtype") 
    goi_table_server("goi_table", reactive({data()}), "subtype")
    
    # ploting differential gene analysis 
    deg <- reactive({
      df <- subtype_deg_df %>% filter(symbol %in% input$geneid)
      return(df)
    }) %>% bindCache(input$geneid)
    
    deg_plot_server("deg_plot_subtype", reactive({deg()}))
  })

  ############################# need change 
  observeEvent(input$load3, {
    data <- reactive({
      TPM_mouse[TPM_mouse[,21] %in% input$geneid,]
    }) %>% bindCache(input$geneid)
    plotdot_server("dot_mouse", preprocess(reactive({data()}), 20, "SNI", TPM_mouse_colData, input$sex, "mouse", "dot"), input$sex) 
    plotline_server("mouse_line", preprocess(reactive({data()}), 20, "SNI", TPM_mouse_colData, input$sex, "mouse", "line"), input$sex, "mouse") 
    plotsubtype_server("mouse_lines_subtype", preprocess(reactive({data()}), 20, "SNI", TPM_mouse_colData, input$sex, "mouse", "line_subtype"), input$sex, "mouse") 
    goi_table_server("goi_table_mouse", reactive({data()}), "mouse")
    
    # ploting differential gene analysis 
    deg <- reactive({
      df <- mouse_deg_df %>% filter(symbol %in% input$geneid)
      return(df)
    }) %>% bindCache(input$geneid)
    
    deg_plot_server("deg_plot_mouse", reactive({deg()}))
    
  })
  ##################################

  # display table in 'hypothesis testing' page 
  observeEvent(input$load4, {
    data <- reactive({
      counts_rat[counts_rat[,6] %in% input$geneid,]
    })
    

    plotdot_server("dot_rat", prep(reactive({data()})),input$sex)

  })
  usercontrast <- reactive({
    req(input$contrast)
    res <- fread(input$contrast)
    return(res)
  }) %>% bindCache(input$contrast)
  
  contrast_table_server("contrast_table",reactive({usercontrast()})) 
  
  # ploting volcano plots
  observeEvent(input$volc, {
    res = fread(paste0(PATH, input$volca))
    colnames(res) <- c('symbol', 'log2FoldChange', 'padj', 'ID')
    res = mutate(res, log10fdr=-log10(padj))
    volcano_plot_server("volcano", input$geneid, res)
  })
  
  ################################ volcano plots for mouse data
  observeEvent(input$volcmouse, {
    res = fread(paste0(TABLE_PATH, input$volcamouse))
    res = mutate(res, log10fdr=-log10(padj))
    volcano_plot_server("volcano_mouse", input$geneid, res)
  })
  ############
  
  usercontrast_mouse <- reactive({
    req(input$contrastm)
    res <- fread(paste0(TABLE_PATH, input$contrastm))
    return(res)
  }) 
  contrast_table_server("mouse_contrast_table",reactive({usercontrast_mouse()})) 
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("drgdirectory_search", ".csv", sep = "")
    },
    
    content = function(file) {
      datatable <- data()
      write.csv(datatable, file, row.names = FALSE)
    })
  
  PlotHeight = reactive(
    return(length(data()))
  )
  
  ### leaflet map for contact details
  output$myMap <- renderLeaflet({
    m <- leaflet() %>% addTiles()
    m <- m %>% setView( -1.238233, 51.756192, zoom = 13)
    m %>% addPopups(-1.2217, 51.76529, "Neural Injury Group")
  })
  
})

