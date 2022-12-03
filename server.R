##########################################
##   Shiny for -omics data presentation ##
##   Allison Barry                      ##
##   University of Oxford               ##
##   allimariebarry@gmail.com           ##
##   for non-commercial use only        ##
##########################################
library(profvis)
# load data  
load("/Users/lynnezhao/Desktop/drg/database.RData")
pbmc = readRDS("/Users/lynnezhao/Desktop/drg.combined.rds")
dataset = read.csv("/Users/lynnezhao/Desktop/dataset.csv",header = TRUE)
colnames(dataset) = c("Title", "Pain Model", "Citation", "Release Date")

source("/Users/lynnezhao/Desktop/labels.R") # include the theme and labels 

# data preprocessing; return median expression of a given gene list 
preprocess <- function(df, rownum, paintype, colData, sex, dataset, graphtype) {
  
  # first, make sure the input colData has 'Timepoint', 'Dataset', 'Species', 'symbol', 'Population', 'Sex'
  matfilt <- df()[,1:rownum]
  tcounts <- t(matfilt) %>%
    base::merge(colData, ., by="row.names") %>%
    gather(gene, expression, (ncol(.)-nrow(matfilt)+1):(ncol(.)))
  if ((dataset == "mouse" )|(dataset == "subtype")) {
    tcounts$symbol <- gene_data[tcounts$gene,]$mgi_symbol # add gene symbols for goi
  }
  if (dataset == "rat") {
    tcounts$symbol <- rat_gene_data[tcounts$gene,]$rgd_symbol 
  }
  if ((dataset == "iPSC_SN")|(dataset == "skin")) {
    tcounts$symbol <- human_gene_data[tcounts$gene,]$hgnc_symbol
  }
  if (graphtype == "dot") {
    if (sex == "Both") {
      tcounts_med <- tcounts %>% dplyr::group_by(Condition, Population, symbol, Dataset, Species) %>% 
        dplyr::summarise(expression=median(expression))
    }
    if (sex == "Separate") {
      tcounts_med <- tcounts %>% dplyr::group_by(Condition, Population, symbol, Sex, Dataset, Species) %>% 
        dplyr::summarise(expression=median(expression))
    }
    tcounts_med <- tcounts_med[!tcounts_med$Condition %in% paintype, ]
  }
  if (graphtype == "line") {
    if (sex == "Both") {
      tcounts_med <- tcounts %>% dplyr::group_by(Condition, symbol, Timepoint, Dataset, Species) %>% 
        dplyr::summarise(expression=median(expression))
    }
    if (sex == "Separate") {
      tcounts_med <- tcounts %>% dplyr::group_by(Condition, symbol, Sex, Timepoint, Dataset, Species) %>% 
        dplyr::summarise(expression=median(expression))
    }
  }
  if (graphtype == "line_subtype") {
    if (sex == "Both") {
      tcounts_med <- tcounts %>% dplyr::group_by(Condition, symbol, Timepoint, Population, Dataset, Species) %>% 
        dplyr::summarise(expression=median(expression))
    }
    if (sex == "Separate") {
      tcounts_med <- tcounts %>% dplyr::group_by(Condition, symbol, Sex, Timepoint, Population, Dataset, Species) %>% 
        dplyr::summarise(expression=median(expression))
    }
  }
  return(tcounts_med)
 
}
# plot the combined dot plot 
plotcombine_server <- function(id, df, sex, genes) {
  moduleServer(id, function(input, output, session) {
    df_rat = df[df$Species == "Rat (DRG)",]
    df_mouse = df[df$Species == "Mouse (DRG)",]
    df_human = df[df$Species == "human",]
    
    g1 = ggplot(df_rat, aes(x= Population, y = symbol)) + scale_colour_viridis_c(option = "magma", end = .90) +
      scale_fill_continuous(limits=c(-10,10)) + scale_size_continuous(limits=c(-10,10)) +
      geom_point(aes(col=log(expression), size=log(expression))) + thc + guides(col = FALSE, size = FALSE) +
      facet_grid(.~Dataset, scale = "free", space='free') + ggtitle("Rat (DRG)") + scale_x_discrete(labels=population_labels)
    g2 = ggplot(df_human, aes(x= Population, y = symbol)) + scale_colour_viridis_c(option = "magma", end = .90) +
      scale_fill_continuous(limits=c(-10,10)) + 
      geom_point(aes(col=log(expression), size=log(expression))) + thc + guides(col = FALSE, size = FALSE) +
      facet_grid(.~Dataset, scale = "free", space='free') + ggtitle("Human") + scale_x_discrete(labels=population_labels)
    g3 = ggplot(df_mouse, aes(x= Population, y = symbol)) + scale_colour_viridis_c(option = "magma", end = .90) +
      scale_fill_continuous(limits=c(-10,10)) +
      geom_point(aes(col=log(expression), size=log(expression))) + thc + guides(size = FALSE) +
      facet_grid(.~Dataset, scale = "free", space='free') + ggtitle("Mouse (DRG)") + 
      scale_x_discrete(labels=population_labels) + labs(col="log(TPM)", size = "")
    
    if (sex == "Both"){
       output$dot <- renderCachedPlot({
         plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))
      }, cacheKeyExpr = {list(df, genes)})
       
      
    }
    if (sex == "Separate") {
      g1 = g1 + facet_grid(.~Dataset+Sex, labeller = labeller(Sex = sexlabels), scale ="free",space='free')
      g2 = g2 + facet_grid(.~Dataset+Sex, labeller = labeller(Sex = sexlabels), scale ="free",space='free')
      g3 = g3 + facet_grid(.~Dataset+Sex, labeller = labeller(Sex = sexlabels), scale ="free",space='free')
      output$dot <- renderCatchedPlot({
        plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))
      }, cacheKeyExpr = {list(df,genes)})
      
    }
    output$downloadDot <- downloadHandler(
      filename = function() {
        paste("plot1", ".png", sep = "")
      },
      
      content = function(file) {
        p <- plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))
        ggsave(p, filename = file, width = 14, height = 4, dpi = 300, units = "in", device='png')
      }
      )
 
  })
}

# plot line plots for subtype data and mouse data   
plotline_server <- function(id, df, sex, dataset) {
  moduleServer(id, function(input, output, session) {
    tcounts_med = df
    g = ggplot(data=tcounts_med, aes(x=Condition, y=expression, group=interaction(symbol, Timepoint))) + 
      scale_colour_viridis(discrete=TRUE, end = .80) + 
      geom_line(aes(color=symbol, linetype=Timepoint)) + geom_point(aes(col=symbol)) + theme_line + ylab("Expression") +
      guides(linetype=guide_legend(label.hjust=0.5)) 
    if (sex == 'Both') {
      output$bulkseq_lines <- renderPlot({g})
    }
    if (sex =='Separate'){
      output$bulkseq_lines <- renderPlot({
        g + facet_wrap(~Sex, labeller = labeller(Sex = sexlabels), scales = "free_x") + labs(col="") 
      })
    }
    
  })
}

plotdot_server <- function(id, df, sex) {
  moduleServer(id, function(input, output, session) {
    tcounts_med = df
    g = ggplot(data = tcounts_med, aes(x= Population, y = symbol)) + 
      scale_colour_viridis_c(option = "magma", end = .90) +
      geom_point(aes(col=log(expression), size=log(expression))) + th + 
      facet_grid(.~Dataset, scales = "free", space="free") + 
      scale_x_discrete(labels=population_labels) + labs(col="log(TPM)", size = "")
    
    if (sex == 'Both') {
      output$bulkseq_dots <- renderPlotly({g })
    }
    if (sex == "Separate") {
      output$bulkseq_dots <- renderPlotly({ 
        g + facet_wrap(~Sex, labeller = labeller(Sex = sexlabels), scales ="free_x")
      })
    }
  })
}

plotscdot_server <- function(id, pbmc, genes) {
  moduleServer(id, function(input, output, session) {
    output$scrna_dots <- renderPlotly({
      hg = subset(human_gene_data, mgi_symbol %in% genes)$hgnc_symbol
      DotPlot(pbmc, features = hg) + 
        theme(axis.title = element_blank(),
              axis.text.x = element_text(size=8, angle = 45, hjust= 1),
              axis.text.y = element_text(size=8),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(), legend.justification = c(0,0.3),
              legend.title = element_text(size=10), legend.key.size = unit(0.4, "cm"))
    })
  })
}

plotumap_server <- function(id, pbmc) {
  moduleServer(id, function(input, output, session) {
    output$scrna_umap <- renderPlotly({
      DimPlot(pbmc) 
    })
  })
}

plotfeature_server <- function(id, pbmc, gene) {
  moduleServer(id, function(input, output, session) {
    output$scrna_feature <- renderCachedPlot({
      hg = subset(human_gene_data, mgi_symbol==gene)$hgnc_symbol
      FeaturePlot(pbmc, features = hg)
    }, cacheKeyExpr = {list(pbmc,gene)})
  })
}

plothomescdot_server <- function(id, pbmc, genes) {
  moduleServer(id, function(input, output, session) {
    hg = subset(human_gene_data, mgi_symbol %in% genes)$hgnc_symbol
    g = DotPlot(pbmc, features = hg) + thc + ggtitle("Human DRG subtypes", subtitle="Spatial-seq (Tavares-Ferrelra 2021)")
    output$home_scrna_dots <- renderCachedPlot({g}, cacheKeyExpr = {list(genes, pbmc)})
    
    output$downloadscDot <- downloadHandler(
      filename = function() {
        paste("scplot", ".png", sep = "")
      },
      
      content = function(file) {
        ggsave(g, filename = file, width = 10, height = 8, dpi = 300, units = "in", device='png')
      }
    )
    
  })
}

# ploting subtype plots; 
plotsubtype_server <- function(id, df, sex, dataset) {
  moduleServer(id, function(input, output, session) {
    tcounts_med = df
    pop_num = n_distinct(df$Population)
    g =  ggplot(data=tcounts_med, aes(x=Condition, y=expression, group=interaction(symbol, Timepoint))) +  
      scale_colour_viridis(discrete=TRUE, end = .80) + 
      geom_line(aes(col=symbol, linetype=Timepoint)) + geom_point(aes(col=symbol)) + 
      # facet_grid(symbol~Population, labeller=labeller(Population=population_labels)) +
      facet_wrap(~Population, ncol = pop_num, labeller=labeller(Population=population_labels))+ 
      theme_line + ylab("Expression") 
    if (sex == "Both") {
      output$bulkseq_lines_subtype <- renderPlot({g})
    }
    else {
      output$bulkseq_lines_subtype <- renderPlot({g + 
          facet_grid(Sex~Population, labeller=labeller(Population=population_labels, Sex=sexlabels)) + labs(col="")
      })
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

# plot deg plots 
deg_plot_server <- function(id, df) {
  moduleServer(id, function(input, output, session){
    output$deg_plot <- renderPlotly({
      ggplot(data = df(), aes(x=interaction(Population), y = symbol, 
                              text = paste('padj: ',padj))) + 
        scale_colour_viridis_c(option = "viridis", end = .90) +
        geom_point(aes(col=log2FoldChange, shape=sig, size=0.5)) + scale_shape_manual(values=c(1, 19)) +
        th + labs(shape = "", size = "") + facet_wrap(~Dataset, scales = "free_x") + scale_x_discrete(labels =subpopulation_labels)
    })
    
  })
}
# plot the combined deg plots
deg_combine_server <- function(id, datar, datam, datah) {
  moduleServer(id, function(input, output, session){ 
    df_rat = datar()
    df_mouse = datam()
    df_human = datah()

    
    g1 = ggplot(data = df_rat, aes(x=interaction(Population), y = symbol, text = paste('padj: ',padj))) +
      geom_point(aes(col=log2FoldChange, shape=sig, size=0.3)) + ggtitle("Rat (DRG)") +
      labs(shape = "", size = "") + facet_grid(.~Dataset, scale = "free", space='free')  +
      guides(col = FALSE, shape = FALSE, sig=FALSE, size = FALSE) +
      scale_colour_continuous(limits=c(0,1)) + scale_x_discrete(labels=subpopulation_labels) +
      scale_colour_viridis_c(option = "viridis", end = .90) +
      scale_shape_manual(values=c(1, 19)) + thc + labs(shape = "", size = "")
    
    
    g2 = ggplot(data = df_human, aes(x=interaction(Population), y = symbol, text = paste('padj: ',padj))) + ggtitle("Human") +
      geom_point(aes(col=log2FoldChange, shape=sig, size=0.3)) + facet_grid(.~Dataset, scale = "free", space='free') +
      guides(col = FALSE, shape = FALSE, sig=FALSE, size = FALSE) + scale_x_discrete(labels=subpopulation_labels) +
      scale_colour_continuous(limits=c(0,1)) +
      scale_colour_viridis_c(option = "viridis", end = .90) +
      scale_shape_manual(values=c(1, 19)) + thc + labs(shape = "", size = "")
    
    g3 = ggplot(data = df_mouse, aes(x=interaction(Population), y = symbol, text = paste('padj: ',padj))) +  
      ggtitle("Mouse (DRG)") + geom_point(aes(col=log2FoldChange, shape=sig, size=0.3)) + 
    facet_grid(.~Dataset, scale = "free", space='free') + scale_x_discrete(labels=subpopulation_labels) +
      scale_colour_continuous(limits=c(0,1)) +
      scale_colour_viridis_c(option = "viridis", end = .90) +
      scale_shape_manual(values=c(1, 19)) + thc + labs(shape = "", size = "") + guides(size=FALSE)
    
    output$deg <- renderCachedPlot({plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))},
                             cacheKeyExpr = {list(datar, datah, datam)})
    
    output$downloadDeg <- downloadHandler(
      filename = function() {
        paste("deg_plot", ".png", sep = "")
      },
      
      content = function(file) {
        p <- plot_grid(g1,g2, g3,ncol=3, rel_widths = c(1/11,4/11,6/11))
        ggsave(p, filename = file, width = 14, height = 4, dpi = 300, units = "in", device='png')
      })

  })
}
# plot volcano plots
volcano_plot_server <- function(id, geneids, df) {
  moduleServer(id, function(input, output, session){
    output$volcano <- renderCachedPlot({
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
                           limits = c(-10, 10)) + theme_line + 
        labs(y="-log10(FDR)")}, cacheKeyExpr = {list(df, geneids)})
  })
} 

credentials <- data.frame(
  user = c("shiny", "ndcn"),
  password = c("azerty", "ndcn-rnaseq"),
  stringsAsFactors = FALSE
)

################################################ SERVER ##############################################################################
shinyServer(function(input, output, session) {
  load(data_dir)
  
  # select genes 
  updateSelectizeInput(session, 
                       inputId = "geneid", 
                       label = "Search Genes:",
                       choices = TPM_subtype[,80], 
                       server = TRUE,
                       selected = "Atf3"
  ) 
  
  updateSelectizeInput(session, 
                       inputId = "scgeneid", 
                       label = "Search Genes:",
                       choices = TPM_subtype[,80], 
                       server = TRUE,
                       selected = "Atf3"
  )
  
  # get a list of genes, separated by comma 
  updateTextInput(session, 
                       inputId = "genes", 
                       label = "Search Genes:"
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
  
  observeEvent(input$link_to_scrna, {
    newvalue <- "tabspat"
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
  
  observeEvent(input$link_to_home3, {
    newvalue <- "tabhome"
    updateTabItems(session, "tabs", newvalue)
  })
  observeEvent(input$link_to_home4, {
    newvalue <- "tabhome"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$link_to_home5, {
    newvalue <- "tabhome"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$link_to_home6, {
    newvalue <- "tabhome"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$link_to_home7, {
    newvalue <- "tabhome"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$link_to_human, {
    newvalue <- "tabhuman"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$link_to_db, {
    newvalue <- "tabdb"
    updateTabItems(session, "tabs", newvalue)
  })
  
  observeEvent(input$link_to_cts, {
    newvalue <- "tabcts"
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
  
  # for rat dataset
  rat = ratdrg[c("log2FoldChange", "padj", "symbol")]
  rat$Population = rep("DRG", nrow(rat))
  rat$Dataset = rep("rat", nrow(rat))
  rat_deg_df = mutate(rat, sig=ifelse(rat$padj<0.05, "SIG", "NS"))
  
  # for human dataset
  young = young[c("log2FoldChange", "padj", "hsymbol")]
  colnames(young) = c("log2FoldChange", "padj", "symbol")
  young$Population = rep("young", nrow(young))
  young$Dataset = rep("iPSC_SN", nrow(young))
  
  old = old[c("log2FoldChange", "padj", "hsymbol")]
  colnames(old) = c("log2FoldChange", "padj", "symbol")
  old$Population = rep("old", nrow(old))
  old$Dataset = rep("iPSC_SN", nrow(old))
  human_deg_df = rbind(young, old)
  human_deg_df = mutate(human_deg_df, sig=ifelse(human_deg_df$padj<0.05, "SIG", "NS"))
  
  # for two human datasets 
  DB = DB[c("log2FoldChange", "padj", "symbol")]
  DB$Population = rep("Diabetes", nrow(DB))
  DB$Dataset = rep("pdn", nrow(DB))
  
  DB_female = DB_female[c("log2FoldChange", "padj", "symbol")]
  DB_female$Population = rep("Diabetes_male", nrow(DB_female))
  DB_female$Dataset = rep("pdn", nrow(DB_female))
  
  DB_male = DB_male[c("log2FoldChange", "padj", "symbol")]
  DB_male$Population = rep("Diabetes_female", nrow(DB_male))
  DB_male$Dataset = rep("pdn", nrow(DB_male))
  
  # for skin 
  HS = HS[c("log2FoldChange", "padj", "symbol")]
  HS$Population = rep("skin", nrow(HS))
  HS$Dataset = rep("cts", nrow(HS))
  
  db_deg_df = rbind(DB, DB_female, DB_male)
  db_deg_df = mutate(db_deg_df, sig=ifelse(db_deg_df$padj<0.05, "SIG", "NS"))
  
  cts_deg_df = HS
  cts_deg_df = mutate(cts_deg_df, sig=ifelse(cts_deg_df$padj<0.05, "SIG", "NS"))
  
  # for subtype dataset 
  subtype_deg_df = data.frame()
  subpopulations_df = list(TDNV_3D,TDNV_4W,CGRT_3D,CGRT_4W,MRTD_3D, MRTD_4W,CRTH_3D,CRTH_4W,TBAC_3D,TBAC_4W)
  labels = c(
    "Nociceptors 3D",
    "Nociceptors 4W",
    "PEP 3D",
    "PEP 4W", 
    "NP 3D",
    "NP 4W",
    "C-LTMR 3D",
    "C-LTMR 4W",
    "Ad- AB-RA LTMRs 3D",
    "Ad- AB-RA LTMRs 4W")
  for (i in c(1:length(subpopulations_df))) {
    res = data.frame(subpopulations_df[i])[c("log2FoldChange", "padj", "symbol")]
    res$Population <- rep(labels[i], nrow(res)) # add the population label
    res$Dataset = rep("subtype", nrow(res))
    subtype_deg_df = bind_rows(subtype_deg_df, res)
  }
  subtype_deg_df = mutate(subtype_deg_df, sig=ifelse(subtype_deg_df$padj<0.05, "SIG", "NS"))
  subtype_deg_df = subtype_deg_df[c("log2FoldChange", "padj", "symbol", "Population","Dataset", "sig")]
  
  mouse_all_deg_df = rbind(mouse_deg_df, subtype_deg_df)
  human_all_deg_df = rbind(human_deg_df, db_deg_df, cts_deg_df)
  
  combined_deg_df = rbind(mouse_all_deg_df, rat_deg_df, human_all_deg_df)

  # a reactive variable that records whether a file is uploaded. 
  rv <- reactiveValues(
    clear = FALSE,
    data = FALSE
  )
  
  observeEvent(input$file, {
    rv$clear <- FALSE
    rv$data = TRUE
  }, priority = 1000)
  
  observeEvent(input$reset, {
    rv$clear = TRUE
    rv$data = FALSE
  })
  
  # plot homepage plots 
  observeEvent(input$load, {
    
    file_input <- reactive({
        if (rv$clear == TRUE) {
          return(NULL)
        }
      if(rv$clear==FALSE && rv$data == TRUE) {
        goi = read.table(input$file$datapath)
        rownames(goi) <- goi[,1]
        goi <- goi[which(rownames(goi) %in% TPM_subtype[,80]==TRUE),]
        return(goi)}
    })
    
    if (is.null(file_input())) {
      genes = input$geneid
    }
    else {
      genes = file_input()
    }
    
    data1 <- reactive({
      TPM_mouse[TPM_mouse[,21] %in% genes,]
    }) 
    data2 <- reactive({
      TPM_subtype[TPM_subtype[,80] %in% genes,]
    })
    data3 <- reactive({
      TPM_rat[TPM_rat[,9] %in% genes,]
    })
    
    data4 <- reactive({
      TPM_ipsc[TPM_ipsc[,30] %in% genes,]
    })
    
    data5 <- reactive({
      TPM_HS_CTS[TPM_HS_CTS[,96] %in% genes,]
    })
    
    data6 <- reactive({
      TPM_HS_diabetes[TPM_HS_diabetes[,108] %in% genes,]
    })
    
    data7 <- reactive({
      TPM_zheng[TPM_zheng[,27] %in% genes,]
    })
    
    data8 <- reactive({
      hg = subset(human_gene_data, mgi_symbol %in% genes)$hgnc_symbol
      TPM_humandrg[TPM_humandrg[,52] %in% hg,]
    })
    
    df_mouse = preprocess(reactive({data1()}), 20, "SNI", TPM_mouse_colData, input$sex, "mouse", "dot")
    df_subtype = preprocess(reactive({data2()}), 79, "Ipsi", bulkseq_colData, input$sex,"subtype", "dot")
    df_rat = preprocess(reactive({data3()}), 8, "SNT", TPM_rat_colData, input$sex, "rat", "dot")
    df_ipsc = preprocess(reactive({data4()}), 28, "patient", ipsc_colData, input$sex, "iPSC_SN", "dot")
    df_skin_CTS = preprocess(reactive({data5()}), 94, "pre_Surgery", skin_colData, input$sex, "skin", "dot")
    df_diabetes = preprocess(reactive({data6()}), 106, "PDPN", db_colData, input$sex, "skin", "dot")
    df_zheng = preprocess(reactive({data7()}), 26, "", zheng_colData,input$sex, "mouse", "dot")
    df_humandrg = preprocess(reactive({data8()}), 51, "P", humandrg_colData,input$sex, "skin", "dot")
    
    r = rbind(df_mouse, df_subtype, df_rat, df_ipsc, df_skin_CTS, df_diabetes, df_zheng, df_humandrg)
    
    plotcombine_server("dot", r, input$sex, genes) #dotplot
    
    plothomescdot_server("homespat", pbmc, genes)

    degr <- reactive({
      rg = subset(rat_gene_data, mgi_symbol %in% genes)$rgd_symbol
      df <- rat_deg_df %>% filter(symbol %in% rg)
      return(df)
    }) %>% bindCache(genes)
    
    degm <- reactive({
      df <- mouse_all_deg_df %>% filter(symbol %in% genes)
      return(df)
    }) %>% bindCache(genes)
    
    degh <- reactive({
      hg = subset(human_gene_data, mgi_symbol %in% genes)$hgnc_symbol
      df <- human_all_deg_df %>% filter(symbol %in% hg)
      return(df)
    }) %>% bindCache(genes)
    
    deg_combine_server("deg_plot", reactive({degr()}), reactive({degm()}), reactive({degh()}))
    
  })
  
  # plots for subtype dataset page   
  observeEvent(input$load2, {
    file_input <- reactive({
      if (rv$clear == TRUE) {
        return(NULL)
      }
      if(rv$clear==FALSE && rv$data == TRUE) {
        goi = read.table(input$file$datapath)
        rownames(goi) <- goi[,1]
        goi <- goi[which(rownames(goi) %in% TPM_subtype[,80]==TRUE),]
        return(goi)}
    })
    
    if (is.null(file_input())) {
      genes = input$geneid
      
    }
    else {
      genes = file_input()
    }
    data <- reactive({
      TPM_subtype[TPM_subtype[,80] %in% genes,]
    }) %>% bindCache(genes)
    data2 <- reactive({
      bulkseq_mat[bulkseq_mat[,155] %in% genes,]
    }) %>% bindCache(genes)
    plotdot_server("dot_subtype", preprocess(reactive({data()}), 79, "ipsi", bulkseq_colData, input$sex,"subtype", "dot"), input$sex)
    plotline_server("line", preprocess(reactive({data2()}), 154, "ipsi", bulkseq_colData, input$sex,"subtype", "line"), input$sex, "subtype") # lineplot
    plotsubtype_server("lines_subtype", preprocess(reactive({data2()}), 154, "ipsi", bulkseq_colData, input$sex,"subtype", "line_subtype"), input$sex, "subtype") 
    goi_table_server("goi_table", reactive({data()}), "subtype")
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("drgdirectory_search", ".csv", sep = "")
      },
      
      content = function(file) {
        datatable <- data()
        write.csv(datatable, file, row.names = FALSE)
      })
    deg <- reactive({
      df <- subtype_deg_df %>% filter(symbol %in% genes)
      return(df)
    }) %>% bindCache(genes)
    
    deg_plot_server("deg_plot_subtype", reactive({deg()}))
  })

  observeEvent(input$load3, {
    file_input <- reactive({
      if (rv$clear == TRUE) {
        return(NULL)
      }
      if(rv$clear==FALSE && rv$data == TRUE) {
        goi = read.table(input$file$datapath)
        rownames(goi) <- goi[,1]
        goi <- goi[which(rownames(goi) %in% TPM_subtype[,80]==TRUE),]
        return(goi)}
    })
    
    if (is.null(file_input())) {
      genes = input$geneid
      
    }
    else {
      genes = file_input()
    }
    data <- reactive({
      TPM_mouse[TPM_mouse[,21] %in% genes,]
    }) %>% bindCache(genes)
    plotdot_server("dot_mouse", preprocess(reactive({data()}), 20, "SNI", TPM_mouse_colData, input$sex, "mouse", "dot"), input$sex) 
    plotline_server("mouse_line", preprocess(reactive({data()}), 20, "SNI", TPM_mouse_colData, input$sex, "mouse", "line"), input$sex, "mouse") 
    plotsubtype_server("mouse_lines_subtype", preprocess(reactive({data()}), 20, "SNI", TPM_mouse_colData, input$sex, "mouse", "line_subtype"), input$sex, "mouse") 
    goi_table_server("goi_table_mouse", reactive({data()}), "mouse")
    deg <- reactive({
      df <- mouse_deg_df %>% filter(symbol %in% genes)
      return(df)
    }) %>% bindCache(genes)
    
    deg_plot_server("deg_plot_mouse", reactive({deg()}))
    
  })
  # for rat page 
  observeEvent(input$load4, {
    file_input <- reactive({
      if (rv$clear == TRUE) {
        return(NULL)
      }
      if(rv$clear==FALSE && rv$data == TRUE) {
        goi = read.table(input$file$datapath)
        rownames(goi) <- goi[,1]
        goi <- goi[which(rownames(goi) %in% TPM_subtype[,80]==TRUE),]
        return(goi)}
    })
    
    if (is.null(file_input())) {
      genes = input$geneid
      
    }
    else {
      genes = file_input()
    }
    data <- reactive({
      TPM_rat[TPM_rat[,9] %in% genes,]
    })
    plotdot_server("dot_rat", preprocess(reactive({data()}), 8, "SNT", TPM_rat_colData, input$sex, "rat", "dot"),input$sex)
    plotline_server("rat_line", preprocess(reactive({data()}), 8, "SNT", TPM_rat_colData, input$sex, "rat", "line"), input$sex, "rat")
    
    deg <- reactive({
      rg = subset(rat_gene_data, mgi_symbol %in% genes)$rgd_symbol
      df <- rat_deg_df %>% filter(symbol %in% rg)
      return(df)
    }) %>% bindCache(genes)
    
    deg_plot_server("deg_plot_rat", reactive({deg()}))
    goi_table_server("goi_table_rat", reactive({data()}), "rat")
  })
  
  # for ipsc page 
  observeEvent(input$load5, {
    file_input <- reactive({
      if (rv$clear == TRUE) {
        return(NULL)
      }
      if(rv$clear==FALSE && rv$data == TRUE) {
        goi = read.table(input$file$datapath)
        rownames(goi) <- goi[,1]
        goi <- goi[which(rownames(goi) %in% TPM_subtype[,80]==TRUE),]
        return(goi)}
    })
    
    if (is.null(file_input())) {
      genes = input$geneid
      
    }
    else {
      genes = file_input()
    }
    data <- reactive({
      TPM_ipsc[TPM_ipsc[,30] %in% genes,]
    })
    
    plotdot_server("dot_human", preprocess(reactive({data()}), 28, "patient", ipsc_colData, input$sex, "iPSC_SN", "dot"),input$sex)
    new_colData = ipsc_colData
    new_colData$Population = ipsc_colData$cell_line
    plotdot_server("dot_human2", preprocess(reactive({data()}), 28, "patient", new_colData, input$sex, "iPSC_SN", "dot"),input$sex)
    plotline_server("human_line", preprocess(reactive({data()}), 28, "patient", ipsc_colData, input$sex, "iPSC_SN", "line"), input$sex, "iPSC_SN")
    deg <- reactive({
      hg = subset(human_gene_data, mgi_symbol %in% genes)$hgnc_symbol 
      df <- human_deg_df %>% filter(symbol %in% hg)
      return(df)
    }) %>% bindCache(genes)
    
    deg_plot_server("deg_plot_human", reactive({deg()}))
    goi_table_server("goi_table_human", reactive({data()}), "iPSC_SN")
    
  })
  
  # for diabetes skin plots 
  observeEvent(input$load6, {
    file_input <- reactive({
      if (rv$clear == TRUE) {
        return(NULL)
      }
      if(rv$clear==FALSE && rv$data == TRUE) {
        goi = read.table(input$file$datapath)
        rownames(goi) <- goi[,1]
        goi <- goi[which(rownames(goi) %in% TPM_subtype[,80]==TRUE),]
        return(goi)}
    })
    
    if (is.null(file_input())) {
      genes = input$geneid
      
    }
    else {
      genes = file_input()
    }
    data <- reactive({
      TPM_HS_diabetes[TPM_HS_diabetes[,108] %in% genes,]
    })
    
    plotdot_server("dot_db", preprocess(reactive({data()}), 106, "PDPN", db_colData, input$sex, "skin", "dot"),input$sex)
    plotline_server("db_line", preprocess(reactive({data()}), 106, "PDPN", db_colData, input$sex, "skin", "line"), input$sex, "skin")
    
    deg <- reactive({
      hg = subset(human_gene_data, mgi_symbol %in% genes)$hgnc_symbol
      df <- db_deg_df %>% filter(symbol %in% hg)
      return(df)
    }) %>% bindCache(genes)
    
    deg_plot_server("deg_plot_db", reactive({deg()}))
    
    
  })
  
  # for cts plots 
  observeEvent(input$load7, {
    file_input <- reactive({
      if (rv$clear == TRUE) {
        return(NULL)
      }
      if(rv$clear==FALSE && rv$data == TRUE) {
        goi = read.table(input$file$datapath)
        rownames(goi) <- goi[,1]
        goi <- goi[which(rownames(goi) %in% TPM_subtype[,80]==TRUE),]
        return(goi)}
    })
    
    if (is.null(file_input())) {
      genes = input$geneid
      
    }
    else {
      genes = file_input()
    }
    data <- reactive({
      TPM_HS_CTS[TPM_HS_CTS[,96] %in% genes,]
    })
    
    plotdot_server("dot_cts", preprocess(reactive({data()}), 94, "pre_Surgery", skin_colData, input$sex, "skin", "dot"),input$sex)
    plotline_server("cts_line", preprocess(reactive({data()}), 94, "pre_Surgery", skin_colData, input$sex, "skin", "line"), input$sex, "cts")
    
    deg <- reactive({
      hg = subset(human_gene_data, mgi_symbol %in% genes)$hgnc_symbol
      df <- cts_deg_df %>% filter(symbol %in% hg)
      return(df)
    }) %>% bindCache(genes)
    
    deg_plot_server("deg_plot_cts", reactive({deg()}))
    
    
  })
  
  # spatial seq plots 
  observeEvent(input$load8, {
    file_input <- reactive({
      if (rv$clear == TRUE) {
        return(NULL)
      }
      if(rv$clear==FALSE && rv$data == TRUE) {
        goi = read.table(input$file$datapath)
        rownames(goi) <- goi[,1]
        goi <- goi[which(rownames(goi) %in% TPM_subtype[,80]==TRUE),]
        return(goi)}
    })
    
    if (is.null(file_input())) {
      genes = input$geneid
      
    }
    else {
      genes = file_input()
    }
    plotscdot_server("dotspat", pbmc, genes)
    plotumap_server("umapspat",pbmc)
  })
  
  observeEvent(input$load9,{
    plotfeature_server("featurespat",pbmc,input$scgeneid)
  })
  
  # contrast table for subtype
  usercontrast <- reactive({
    req(input$contrast)
    if (input$contrast == "Nociceptors 3D"){res = TDNV_3D}
    if (input$contrast == "Nociceptors 4W"){res = TDNV_4W}
    if (input$contrast == "PEP 3D"){res = CGRT_3D}
    if (input$contrast == "PEP 4W"){res = CGRT_4W}
    if (input$contrast == "NP 3D"){res = MRTD_3D}
    if (input$contrast == "NP 4W"){res = MRTD_4W}
    if (input$contrast == "C-LTMR 3D"){res = CRTH_3D}
    if (input$contrast == "C-LTMR 4W"){res = CRTH_4W}
    if (input$contrast == "Ad- AB-RA LTMRs 3D"){res = TBAC_3D}
    if (input$contrast == "Ad- AB-RA LTMRs 4W"){res = TBAC_4W}
    return(res)
  }) %>% bindCache(input$contrast)
  contrast_table_server("contrast_table",reactive({usercontrast()})) 
  
  usercontrast_mouse <- reactive({
    req(input$contrastm)
    if (input$contrastm == "B10D2") {
      res = b10d2
    }
    if (input$contrastm == "BALB") {
      res = balb
    }
    return(res)
  })
  contrast_table_server("mouse_contrast_table",reactive({usercontrast_mouse()})) 
  
  # for rat
  output$contrast_table_rat <- DT::renderDataTable({
    DT::datatable(
      ratdrg,
      width = 12,
      class = 'nowrap',
      options = list(scrollX = TRUE, pageLength = 10)
    )
  })
  
  output$contrast_table_cts <- DT::renderDataTable({
    DT::datatable(
      HS,
      width = 12,
      class = 'nowrap',
      options = list(scrollX = TRUE, pageLength = 10)
    )
  })
  
  usercontrast_human <- reactive({
    req(input$contrasth)
    if (input$contrasth == "iPSCDN_young"){res = young}
    if (input$contrasth == "iPSCDN_old"){res = old}
    return(res)
  })
  contrast_table_server("contrast_table_human",reactive({usercontrast_human()})) 
  
  usercontrast_db <- reactive({
    req(input$contrastd)
    if (input$contrastd == "Diabetes_skin") {res = DB}
    if (input$contrastd == "Diabetes_skin_females") {res = DB_female}
    if (input$contrastd == "Diabetes_skin_male") {res = DB_male}
    return(res)
  })
  contrast_table_server("contrast_table_db",reactive({usercontrast_db()})) 
  
  
  # ploting volcano plots for subtype
  observeEvent(input$volc, {

    if (input$volca == "Nociceptors 3D"){res = TDNV_3D}
    if (input$volca == "Nociceptors 4W"){res = TDNV_4W}
    if (input$volca == "PEP 3D"){res = CGRT_3D}
    if (input$volca == "PEP 4W"){res = CGRT_4W}
    if (input$volca == "NP 3D"){res = MRTD_3D}
    if (input$volca == "NP 4W"){res = MRTD_4W}
    if (input$volca == "C-LTMR 3D"){res = CRTH_3D}
    if (input$volca == "C-LTMR 4W"){res = CRTH_4W}
    if (input$volca == "Ad- AB-RA LTMRs 3D"){res = TBAC_3D}
    if (input$volca == "Ad- AB-RA LTMRs 4W"){res = TBAC_4W}
    res = mutate(res, log10fdr=-log10(padj))
    volcano_plot_server("volcano", input$geneid, res)
  })
  
  observeEvent(input$volcmouse, {
    if (input$volcamouse == "B10D2") {
      res = b10d2
    }
    if (input$volcamouse == "BALB") {
      res = balb
    }
    res = mutate(res, log10fdr=-log10(padj))
    volcano_plot_server("volcano_mouse", input$geneid, res)
  })
  #for rat 
  observeEvent(input$volcrat, {
    res = ratdrg
    res = mutate(res, log10fdr=-log10(padj))
    rg = subset(rat_gene_data, rgd_symbol %in% input$geneid)$mgi_symbol 
    volcano_plot_server("volcano_rat", rg, res)
  })
  # for human 
  observeEvent(input$volch, {
    if (input$volcahuman == "iPSCDN_young"){res = young}
    if (input$volcahuman == "iPSCDN_old"){res = old}
    res = mutate(res, log10fdr=-log10(padj))
    hg = subset(human_gene_data, mgi_symbol %in% input$geneid)$hgnc_symbol
    volcano_plot_server("volcano_human", hg, res)
  })
  #for diabetes
  observeEvent(input$volcd, {
    if (input$volcadb == "Diabetes_skin") {res = DB}
    if (input$volcadb == "Diabetes_skin_females") {res = DB_female}
    if (input$volcadb == "Diabetes_skin_male") {res = DB_male}
    res = mutate(res, log10fdr=-log10(padj))
    hg = subset(human_gene_data, mgi_symbol %in% input$geneid)$hgnc_symbol
    volcano_plot_server("volcano_db", hg, res)
  })
  
  observeEvent(input$volcc, {
    res = HS
    res = mutate(res, log10fdr=-log10(padj))
    hg = subset(human_gene_data, mgi_symbol %in% input$geneid)$hgnc_symbol
    volcano_plot_server("volcano_cts", hg, res)
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

options(warn=-1) # remove warnings 

