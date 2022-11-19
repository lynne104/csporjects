##########################################
##   Shiny for -omics data presentation ##
##   Allison Barry                      ##
##   University of Oxford               ##
##   allimariebarry@gmail.com           ##
##   for non-commercial use only        ##
##########################################
library(profvis)
# load data  
load("/Users/lynnezhao/Desktop/drg/drg-directory.RData")
load("/Users/lynnezhao/Desktop/drg/datamining.RData")
load("/Users/lynnezhao/Desktop/drg/ratsnt.RData")
load("/Users/lynnezhao/Desktop/drg/ipsc.RData")
load("/Users/lynnezhao/Desktop/drg/diabetes.RData")

TABLE_PATH = "/Users/lynnezhao/Desktop/data/" 

# path for volcano plots dfs 
PATH = "/Users/lynnezhao/Desktop/drg/"

# themes 
theme_line = theme_bw() + 
  theme(panel.grid = element_blank(),
        axis.title = element_text(size=12),
        axis.text.x = element_text(size=10, angle = 45, hjust= 1),
        axis.text.y = element_text(size=10),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(), legend.justification = c(0,0.3)) 

th = theme_bw() + theme(panel.grid = element_blank(),
                               axis.title = element_blank(),
                               axis.text.x = element_text(size=8, angle = 45, hjust= 1),
                               axis.text.y = element_text(size=8),
                               axis.ticks.x = element_blank(),
                               axis.ticks.y = element_blank(), legend.justification = c(0,0.3))

thc = theme_bw() + theme(plot.title = element_text(hjust = 0.5, size = 12, vjust=0),
  panel.grid = element_blank(), plot.margin=unit(c(10,2,10,2), 'mm'),
                        axis.title = element_blank(),
                        axis.text.x = element_text(size=8, angle = 45, hjust= 1),
                        axis.text.y = element_text(size=8),
                        axis.ticks.x = element_blank(),
                        axis.ticks.y = element_blank(), legend.justification = c(0,0.3))
  

population_labels = c("TBAC" = "A\u03b4-LTMR + A\u03b2-RA-LTMR",
                      "CRTH" = "C-LTMR", "MRTD" = "NP",
                      "CGRT" = "PEP", "TDNV" = "Nociceptors", 
                      "B10.D2" ="b10d2", 
                      "BALB.c" = "balb")
sexlabels = c("F" = "Female", "M" = "Male", "mixed" = "Mixed")
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
  "TBAC_3D" = "Ad- AB-RA LTMRs 3D", "TBAC_4W" = "Ad- AB-RA LTMRs 4W", 
  "B10.D2" ="b10d2", 
  "BALB.c" = "balb"
)
# data preprocessing; return median expression of a given gene list 
preprocess <- function(df, rownum, paintype, colData, sex, dataset, graphtype) {
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

plotcombine_server <- function(id, df, sex) {
  moduleServer(id, function(input, output, session) {
    df_rat = df[df$Species == "Rat (DRG)",]
    df_mouse = df[df$Species == "Mouse (DRG)",]
    df_human = df[df$Species == "human",]
    
    g1 = ggplot(df_rat, aes(x= Population, y = symbol)) + scale_colour_viridis_c(option = "magma", end = .90) +
      scale_fill_continuous(limits=c(-10,10)) + scale_size_continuous(limits=c(-10,10)) +
      geom_point(aes(col=expression, size=expression)) + thc + guides(col = FALSE, size = FALSE) +
      facet_wrap(~Dataset, scales = "free_x") + ggtitle("Rat (DRG)")
    g2 = ggplot(df_human, aes(x= Population, y = symbol)) + scale_colour_viridis_c(option = "magma", end = .90) +
      scale_fill_continuous(limits=c(-10,10)) + 
      geom_point(aes(col=expression, size=expression)) + thc + guides(col = FALSE, size = FALSE) +
      facet_wrap(~Dataset, scales = "free_x") + ggtitle("Human")
    g3 = ggplot(df_mouse, aes(x= Population, y = symbol)) + scale_colour_viridis_c(option = "magma", end = .90) +
      scale_fill_continuous(limits=c(-10,10)) +
      geom_point(aes(col=expression, size=expression)) + thc +
      facet_wrap(~Dataset, scales = "free_x") + ggtitle("Mouse (DRG)") + 
      scale_x_discrete(labels=population_labels) + labs(col="log(TPM)")
    
    if (sex == "Both"){
      output$dot1 <- renderPlotly({g1})
      output$dot2 <- renderPlotly({g2})
      output$dot3 <- renderPlotly({g3})
      
    }
    if (sex == "Separate") {
      g1 = g1 + facet_wrap(~Sex, labeller = labeller(Sex = sexlabels), scales ="free_x")
      g2 = g2 + facet_wrap(~Sex, labeller = labeller(Sex = sexlabels), scales ="free_x")
      g3 = g3 + facet_wrap(~Sex, labeller = labeller(Sex = sexlabels), scales ="free_x")
      output$dot1 <- renderPlotly({g1})
      output$dot2 <- renderPlotly({g2})
      output$dot3 <- renderPlotly({g3})
    }
  })
}

# plot line plots for subtype data and mouse data   
plotline_server <- function(id, df, sex, dataset) {
  moduleServer(id, function(input, output, session) {
    tcounts_med = df
    g = ggplot(data=tcounts_med, aes(x=Condition, y=expression, group=interaction(symbol, Timepoint))) + 
      scale_colour_viridis(discrete=TRUE, end = .80) + 
      geom_line(aes(color=symbol, linetype=Timepoint)) + geom_point(aes(col=symbol)) + theme_line + ylab("Expression") +
      guides(col=FALSE, linetype=guide_legend(label.hjust=0.5, ncol=8)) 
    if (sex == 'Both') {
      output$bulkseq_lines <- renderPlotly({g})
    }
    if (sex =='Separate'){
      output$bulkseq_lines <- renderPlotly({
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
      geom_point(aes(col=expression, size=expression)) + th + 
      facet_grid(.~Dataset, scales = "free", space="free") + 
      scale_x_discrete(labels=population_labels)
    
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

# ploting subtype plots; 
plotsubtype_server <- function(id, df, sex, dataset) {
  moduleServer(id, function(input, output, session) {
    tcounts_med = df
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
        scale_colour_viridis_c(option = "viridis", end = .90) +
        geom_point(aes(col=log2FoldChange, shape=sig, size=0.5)) + scale_shape_manual(values=c(1, 19)) +
        th + labs(shape = "", size = "") + facet_wrap(~Dataset, scales = "free_x")
    }) 
  })
}

deg_combine_server <- function(id, datar, datam, datah) {
  moduleServer(id, function(input, output, session){ 
    df_rat = datar()
    df_mouse = datam()
    df_human = datah()
    
    g1 = ggplot(data = df_rat, aes(x=interaction(Population), y = symbol, 
                                 text = paste('padj: ',padj))) + scale_colour_continuous(limits=c(0,1)) +
      scale_colour_viridis_c(option = "viridis", end = .90) + guides(col = FALSE, shape = FALSE, sig=FALSE) +
      geom_point(aes(col=log2FoldChange, shape=sig, size=0.3)) + scale_shape_manual(values=c(1, 19)) + ggtitle("Rat (DRG)") +
      thc + labs(shape = "", size = "") + facet_wrap(~Dataset, scales = "free_x") 
    g2 = ggplot(data = df_human, aes(x=interaction(Population), y = symbol, 
                                   text = paste('padj: ',padj))) + ggtitle("Human") +
      scale_colour_viridis_c(option = "viridis", end = .90) + guides(col = FALSE, shape = FALSE, sig=FALSE) + 
      geom_point(aes(col=log2FoldChange, shape=sig, size=0.3)) + scale_shape_manual(values=c(1, 19)) +
      thc + labs(shape = "", size = "")+ facet_wrap(~Dataset, scales = "free_x") 
    g3 = ggplot(data = df_mouse, aes(x=interaction(Population), y = symbol, 
                                   text = paste('padj: ',padj))) + 
      scale_colour_viridis_c(option = "viridis", end = .90) + ggtitle("Mouse (DRG)") +
      geom_point(aes(col=log2FoldChange, shape=sig, size=0.3)) + scale_shape_manual(values=c(1, 19)) +
      thc + labs(shape = "", size = "")+ facet_wrap(~Dataset, scales = "free_x") + scale_x_discrete(labels=subpopulation_labels)
    
    output$deg1 <- renderPlotly({g1})
    output$deg2 <- renderPlotly({g2})
    output$deg3 <- renderPlotly({g3})
    
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
                           limits = c(-10, 10)) + theme_line + 
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
                       choices = logTPM_subtype[,80], 
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
  rat = rat[c("log2FoldChange", "padj", "symbol")]
  rat$Population = rep("rat", nrow(rat))
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
  
  mouse_all_deg_df = rbind(mouse_deg_df, subtype_deg_df)
  human_all_deg_df = rbind(human_deg_df, db_deg_df, cts_deg_df)
  
  combined_deg_df = rbind(mouse_all_deg_df, rat_deg_df, human_all_deg_df)

  # plot homepage plots 
  observeEvent(input$load, {
    data1 <- reactive({
      logTPM_mouse[logTPM_mouse[,21] %in% input$geneid,]
    }) 
    data2 <- reactive({
      logTPM_subtype[logTPM_subtype[,80] %in% input$geneid,]
    })
    data3 <- reactive({
      logTPM_rat[logTPM_rat[,9] %in% input$geneid,]
    })
    
    data4 <- reactive({
      logTPM_ipsc[logTPM_ipsc[,30] %in% input$geneid,]
    })
    
    data5 <- reactive({
      logTPM_HS_CTS[logTPM_HS_CTS[,96] %in% input$geneid,]
    })
    
    data6 <- reactive({
      logTPM_HS_diabetes[logTPM_HS_diabetes[,108] %in% input$geneid,]
    })
    
    df_mouse = preprocess(reactive({data1()}), 20, "SNI", TPM_mouse_colData, input$sex, "mouse", "dot")
    df_subtype = preprocess(reactive({data2()}), 79, "Ipsi", bulkseq_colData, input$sex,"subtype", "dot")
    df_rat = preprocess(reactive({data3()}), 8, "SNT", TPM_rat_colData, input$sex, "rat", "dot")
    df_ipsc = preprocess(reactive({data4()}), 28, "patient", ipsc_colData, input$sex, "iPSC_SN", "dot")
    df_skin_CTS = preprocess(reactive({data5()}), 94, "pre_Surgery", skin_colData, input$sex, "skin", "dot")
    df_diabetes = preprocess(reactive({data6()}), 106, "PDPN", db_colData, input$sex, "skin", "dot")
    
    r = rbind(df_mouse, df_subtype, df_rat, df_ipsc, df_skin_CTS, df_diabetes)
    
    plotcombine_server("dot", r, input$sex) #dotplot
 
    degr <- reactive({
      rg = subset(rat_gene_data, mgi_symbol %in% input$geneid)$rgd_symbol
      df <- rat_deg_df %>% filter(symbol %in% rg)
      return(df)
    }) %>% bindCache(input$geneid)
    
    degm <- reactive({
      df <- mouse_all_deg_df %>% filter(symbol %in% input$geneid)
      return(df)
    }) %>% bindCache(input$geneid)
    
    degh <- reactive({
      hg = subset(human_gene_data, mgi_symbol %in% input$geneid)$hgnc_symbol
      df <- human_all_deg_df %>% filter(symbol %in% hg)
      return(df)
    }) %>% bindCache(input$geneid)
    
    deg_combine_server("deg_plot", reactive({degr()}), reactive({degm()}), reactive({degh()}))
  })
  
  # plots for subtype dataset page   
  observeEvent(input$load2, {
    data <- reactive({
      logTPM_subtype[logTPM_subtype[,80] %in% input$geneid,]
    }) %>% bindCache(input$geneid)
    data2 <- reactive({
      bulkseq_mat[bulkseq_mat[,155] %in% input$geneid,]
    }) %>% bindCache(input$geneid)
    plotdot_server("dot_subtype", preprocess(reactive({data()}), 79, "ipsi", bulkseq_colData, input$sex,"subtype", "dot"), input$sex)
    plotline_server("line", preprocess(reactive({data2()}), 154, "ipsi", bulkseq_colData, input$sex,"subtype", "line"), input$sex, "subtype") # lineplot
    plotsubtype_server("lines_subtype", preprocess(reactive({data2()}), 154, "ipsi", bulkseq_colData, input$sex,"subtype", "line_subtype"), input$sex, "subtype") 
    goi_table_server("goi_table", reactive({data()}), "subtype")
    deg <- reactive({
      df <- subtype_deg_df %>% filter(symbol %in% input$geneid)
      return(df)
    }) %>% bindCache(input$geneid)
    
    deg_plot_server("deg_plot_subtype", reactive({deg()}))
  })

  observeEvent(input$load3, {
    data <- reactive({
      logTPM_mouse[logTPM_mouse[,21] %in% input$geneid,]
    }) %>% bindCache(input$geneid)
    plotdot_server("dot_mouse", preprocess(reactive({data()}), 20, "SNI", TPM_mouse_colData, input$sex, "mouse", "dot"), input$sex) 
    plotline_server("mouse_line", preprocess(reactive({data()}), 20, "SNI", TPM_mouse_colData, input$sex, "mouse", "line"), input$sex, "mouse") 
    plotsubtype_server("mouse_lines_subtype", preprocess(reactive({data()}), 20, "SNI", TPM_mouse_colData, input$sex, "mouse", "line_subtype"), input$sex, "mouse") 
    goi_table_server("goi_table_mouse", reactive({data()}), "mouse")
    deg <- reactive({
      df <- mouse_deg_df %>% filter(symbol %in% input$geneid)
      return(df)
    }) %>% bindCache(input$geneid)
    
    deg_plot_server("deg_plot_mouse", reactive({deg()}))
    
  })
  # for rat page 
  observeEvent(input$load4, {
    data <- reactive({
      logTPM_rat[logTPM_rat[,9] %in% input$geneid,]
    })
    plotdot_server("dot_rat", preprocess(reactive({data()}), 8, "SNT", TPM_rat_colData, input$sex, "rat", "dot"),input$sex)
    plotline_server("rat_line", preprocess(reactive({data()}), 8, "SNT", TPM_rat_colData, input$sex, "rat", "line"), input$sex, "rat")
    
    deg <- reactive({
      rg = subset(rat_gene_data, mgi_symbol %in% input$geneid)$rgd_symbol
      df <- rat_deg_df %>% filter(symbol %in% rg)
      return(df)
    }) %>% bindCache(input$geneid)
    
    deg_plot_server("deg_plot_rat", reactive({deg()}))
    goi_table_server("goi_table_rat", reactive({data()}), "rat")
  })
  
  # for ipsc page 
  observeEvent(input$load5, {
    data <- reactive({
      logTPM_ipsc[logTPM_ipsc[,30] %in% input$geneid,]
    })
    plotdot_server("dot_human", preprocess(reactive({data()}), 28, "patient", ipsc_colData, input$sex, "iPSC_SN", "dot"),input$sex)
    plotline_server("human_line", preprocess(reactive({data()}), 28, "patient", ipsc_colData, input$sex, "iPSC_SN", "line"), input$sex, "iPSC_SN")
    
    deg <- reactive({
      hg = subset(human_gene_data, mgi_symbol %in% input$geneid)$hgnc_symbol 
      df <- human_deg_df %>% filter(symbol %in% hg)
      return(df)
    }) %>% bindCache(input$geneid)
    
    deg_plot_server("deg_plot_human", reactive({deg()}))
    goi_table_server("goi_table_human", reactive({data()}), "iPSC_SN")
    
  })
  
  # for diabetes skin plots 
  observeEvent(input$load6, {
    data <- reactive({
      logTPM_HS_diabetes[logTPM_HS_diabetes[,108] %in% input$geneid,]
    })
    
    plotdot_server("dot_db", preprocess(reactive({data()}), 106, "PDPN", db_colData, input$sex, "skin", "dot"),input$sex)
    plotline_server("db_line", preprocess(reactive({data()}), 106, "PDPN", db_colData, input$sex, "skin", "line"), input$sex, "skin")
    
    deg <- reactive({
      hg = subset(human_gene_data, mgi_symbol %in% input$geneid)$hgnc_symbol
      df <- db_deg_df %>% filter(symbol %in% hg)
      return(df)
    }) %>% bindCache(input$geneid)
    
    deg_plot_server("deg_plot_db", reactive({deg()}))
    goi_table_server("goi_table_db", reactive({data()}), "skin")
    
  })
  
  # for cts plots 
  observeEvent(input$load7, {
    data <- reactive({
      logTPM_HS_CTS[logTPM_HS_CTS[,96] %in% input$geneid,]
    })
    
    plotdot_server("dot_cts", preprocess(reactive({data()}), 94, "pre_Surgery", skin_colData, input$sex, "skin", "dot"),input$sex)
    plotline_server("cts_line", preprocess(reactive({data()}), 94, "pre_Surgery", skin_colData, input$sex, "skin", "line"), input$sex, "cts")
    
    deg <- reactive({
      hg = subset(human_gene_data, mgi_symbol %in% input$geneid)$hgnc_symbol
      df <- cts_deg_df %>% filter(symbol %in% hg)
      return(df)
    }) %>% bindCache(input$geneid)
    
    deg_plot_server("deg_plot_cts", reactive({deg()}))
    goi_table_server("goi_table_cts", reactive({data()}), "skin")
    
  })
  
  # contrast table 
  usercontrast <- reactive({
    req(input$contrast)
    res <- fread(input$contrast)
    return(res)
  }) %>% bindCache(input$contrast)
  contrast_table_server("contrast_table",reactive({usercontrast()})) 
  
  usercontrast_mouse <- reactive({
    req(input$contrastm)
    res <- fread(paste0(TABLE_PATH, input$contrastm))
    return(res)
  })
  contrast_table_server("mouse_contrast_table",reactive({usercontrast_mouse()})) 
  
  res = fread(paste0(TABLE_PATH, "Rat_SNT_vs_SHAM.csv")) 
  # for rat
  output$contrast_table_rat <- DT::renderDataTable({
    DT::datatable(
      res,
      width = 12,
      class = 'nowrap',
      options = list(scrollX = TRUE, pageLength = 10)
    )
  })
  
  cts = fread(paste0(TABLE_PATH, "HS_skin_CTS_pre_post_surgery.csv"))
  output$contrast_table_cts <- DT::renderDataTable({
    DT::datatable(
      cts,
      width = 12,
      class = 'nowrap',
      options = list(scrollX = TRUE, pageLength = 10)
    )
  })
  
  usercontrast_human <- reactive({
    req(input$contrasth)
    res <- fread(paste0("/Users/lynnezhao/Desktop/data/", input$contrasth))
    return(res)
  })
  contrast_table_server("contrast_table_human",reactive({usercontrast_human()})) 
  
  usercontrast_db <- reactive({
    req(input$contrastd)
    res <- fread(paste0("/Users/lynnezhao/Desktop/data/", input$contrastd))
    return(res)
  })
  contrast_table_server("contrast_table_db",reactive({usercontrast_db()})) 
  
  
  # ploting volcano plots
  observeEvent(input$volc, {
    res = fread(paste0(PATH, input$volca))
    colnames(res) <- c('symbol', 'log2FoldChange', 'padj', 'ID')
    res = mutate(res, log10fdr=-log10(padj))
    volcano_plot_server("volcano", input$geneid, res)
  })
  
  observeEvent(input$volcmouse, {
    res = fread(paste0(TABLE_PATH, input$volcamouse))
    res = mutate(res, log10fdr=-log10(padj))
    volcano_plot_server("volcano_mouse", input$geneid, res)
  })
  
  observeEvent(input$volcrat, {
    res = fread(paste0(TABLE_PATH, "Rat_SNT_vs_SHAM.csv"))
    res = mutate(res, log10fdr=-log10(padj))
    rg = subset(rat_gene_data, rgd_symbol %in% input$geneid)$mgi_symbol 
    volcano_plot_server("volcano_rat", rg, res)
  })
  
  observeEvent(input$volch, {
    res = fread(paste0(TABLE_PATH, input$volcahuman))
    res = mutate(res, log10fdr=-log10(padj))
    hg = subset(human_gene_data, mgi_symbol %in% input$geneid)$hgnc_symbol
    volcano_plot_server("volcano_human", hg, res)
  })
  
  observeEvent(input$volcd, {
    res = fread(paste0(TABLE_PATH, input$volcadb))
    res = mutate(res, log10fdr=-log10(padj))
    hg = subset(human_gene_data, mgi_symbol %in% input$geneid)$hgnc_symbol
    volcano_plot_server("volcano_db", hg, res)
  })
  
  observeEvent(input$volcc, {
    res = fread(paste0(TABLE_PATH, "HS_skin_CTS_pre_post_surgery.csv"))
    res = mutate(res, log10fdr=-log10(padj))
    hg = subset(human_gene_data, mgi_symbol %in% input$geneid)$hgnc_symbol
    volcano_plot_server("volcano_cts", hg, res)
  })

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

