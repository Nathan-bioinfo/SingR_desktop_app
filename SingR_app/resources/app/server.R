source("func.R")
load_packages()

server = function(input, output, session){
  
  #####################################################################
  
  ######################   INTRODUCTION    ############################
  
  #####################################################################
  
  # intro text
  output$introduction = renderText("
                                   Introduction
                                   Here is an application that i've designed during my trainee period at Nantes's CRTI. 
                                   You'll find several single cell data visualizations.
                                   The program is designed for Seurat build datasets. 
                                   It's absolutely necessary to prepare youre data as explained here: 
                                   https://satijalab.org/seurat/articles/multimodal_vignette.html. 
                                   Help buttons are here to provide more informations.
                                   
                                   Download
                                   All plot and table are downloadable. Notice that the first two plots (interactive quality control)
                                   have several options on the top right corner (zoom, cut...) as the download option.
                                   All the other figures have a download button on the left bottom corner.
                                   
                                   Nichenet
                                   Nichenet uses different calculation and method to help you determine potential ligand and target cells.
                                   If you're not familiar with nichenet, please go to: 
                                   https://github.com/saeyslab/nichenetr/blob/master/vignettes/seurat_wrapper.md
                                   
                                   Disclaimer
                                   I build this app for the special use of Sophie Brouard's team 4 at CRTI.
                                   As R, it's of course open source, but all results and interpretations you'll obtain should be taken carefully.
                                   Staistics and plots are simplified to offer a light and speed app, 
                                   and you should be aware that published results are on you're own responsability.
                                   Please wait until white loading circle on the left bottom part dissapear after pressing buttons.
                                   
                                   Thank you and enjoy
                                   Nathan Caceres")

  #####################################################################

  ######################   UPLOADING       ############################

  #####################################################################

  # to get user path
  volumes = getVolumes()

  # intro tuto bubble
  observeEvent(input$intro_help, {
    shinyalert("Tutorial ", "Choose the example preloaded dataframe or upload your own:
               For personal dataframe, please press \"Choose a file\" button and select a Seurat file to upload.
               Then press load in both case.
               Maximum dataset size is 5 GB. Consider using DietSeurat to reduce youre dataset size.
               If text appears in the \"metadata to keep\" box, you can select those you want and press \"trim\".
               Now youre ready to continue", type = "info", size = "m")
  })

  observe({
    toggle("user_file_path", condition = "dataf" %in% input$file_type)
  })

  observe({
    shinyFileChoose(input, "user_file_path", roots = volumes, session = session)
    if(!is.null(input$user_file_path)){
      file_selected <<-parseFilePaths(volumes, input$user_file_path)
    }
  })

  observeEvent(
    c(input$file_type,input$user_file_path),{
      if("dataf" %in% input$file_type)
      {
        if(length(file_selected$datapath) != 0)
        {
          if (tolower(tools::file_ext(file_selected$datapath)) != "rds")
          {
            disable("preloaded_go")
            shinyalert("Data has to be Seurat object with .rds extension", type = "error")
          }
          else
          {
            enable("preloaded_go")
          }
        }
        else
        {
          disable("preloaded_go")
        }
      }
      else if ("ex" %in% input$file_type)
      {
        enable("preloaded_go")
      }
    })

  observe({
    print(paste("file path = ",input$file_path))
  })

  #get uploaded data if selected
  data = eventReactive(input$preloaded_go,{
    tryCatch(
      {
        rm(list=ls())
        gc()
        if(input$file_type == "dataf")
        {
          predf = readRDS(file_selected$datapath)
        }
        else if (input$file_type == "ex")
        {
          predf = readRDS("./data/data_example.rds")
        }
        reduc = names(predf@reductions)
        df = DietSeurat(
          predf,
          counts = TRUE,
          data = TRUE,
          #scale.data = TRUE,
          dimreducs = reduc
        )
        rm(predf)
        gc()
        # add mitochondrial percent metadata
        if(!exists("percent.mito")){
          mito.genes <- grep(pattern = "^MT-", x = rownames(df@assays[["RNA"]]), value = TRUE)
          percent.mito <- Matrix::colSums(df@assays[["RNA"]][mito.genes, ])/Matrix::colSums(df@assays[["RNA"]])
          df <- AddMetaData(object = df, metadata = percent.mito, col.name = "percent.mito")
          df$percent.mito <- percent.mito
        }
        # choose default assay
        Idents(df) = df$seurat_clusters
        Idents(df) = gsub(" ","_",Idents(df))
        DefaultAssay(df) <- "RNA"
        gc()
        return(df)
      },
      error = function(e)
      {
        shinyalert("Problem uploading data. Please use Seurat Object", type = "error")
        return(NULL)
      }
    )
  })

  observe({
    validate(need(!is.null(data()),"data null"))
    #if(input$file_type == "ex" && input$example_df == "tiny")
    if(input$file_type == "ex")
    {
      updateSelectInput(session, inputId = "trim", choices = colnames(data()@meta.data), selected = c("seurat_clusters","prediction","patient","condition"))
    }
    else
    {
      updateSelectInput(session, inputId = "trim", choices = colnames(data()@meta.data), selected = "seurat_clusters")
    }
  })

  meta_data = eventReactive(input$trim_go,{
    my_list = input$trim
    shinyalert("data set ready", type = "success")
    gc()
    return(my_list)
  })

  observe({
    print(paste("meta_data[1] =",meta_data()[1]))
  })

  # modified data with new idents
  data_2 = reactive({
    validate(need(!is.null(data()),"data null"))
    df = data()
    DefaultAssay(df) = "RNA"
    tryCatch(
      {
        if(input$tabs == "geneset")
        {
          Idents(df) = input$gs_idents
        }
        else if(input$tabs == "boxplot" && input$boxplot_group_by != "")
        {
          Idents(df) = input$boxplot_group_by
        }
        else if(input$tabs == "gene_expression" && input$expression_group_by_choice != "")
        {
          input_meta = input$expression_group_by_choice
          Idents(df) = input_meta
          expr = paste0("df$",input_meta," = gsub(\" \",\"_\",df$",input_meta,")")
          eval(parse(text = expr))
        }
        else if(input$tabs == "gene_expression_2" && input$expression_group_by_choice_2 != "")
        {
          Idents(df) = input$expression_group_by_choice_2
        }
        else if(input$tabs == "niche_calc" && input$niche_idents != "")
        {
          Idents(df) = input$niche_idents
        }
        else if(input$tabs == "differential_calculation" && input$diff_idents != "")
        {
          Idents(df) = input$diff_idents
        }
        else if(input$tabs == "ridge" && input$ridge_group_by != "")
        {
          Idents(df) = input$ridge_group_by
        }
        else if(input$tabs == "dot_plot" && input$dot_plot_group_by != "")
        {
          Idents(df) = input$dot_plot_group_by
        }
        Idents(df) = gsub(" ","_",Idents(df))
        gc()
        return(df)
      },
      error = function(e)
      {
        shinyalert("Error: problem with Idents. Please use Seurat Object", type = "error")
        return(NULL)
      }
    )
  })

  # display gene_list in inputs needed
  observeEvent(input$preloaded_go,{
    validate(need(!is.null(data_2()),"data 2 null"))
    updateSelectizeInput(session, inputId = "dot_plot_genes", choices = rownames(data_2()), server = TRUE, selected = FALSE)
    updateSelectizeInput(session, inputId = "density_gene_1", choices = rownames(data_2()), server = TRUE, selected = FALSE)
    updateSelectizeInput(session, inputId = "density_gene_2", choices = rownames(data_2()), server = TRUE, selected = FALSE)
    updateSelectizeInput(session, inputId = "density_gene_3", choices = rownames(data_2()), server = TRUE, selected = FALSE)
    updateSelectizeInput(session, inputId = "density_gene_4", choices = rownames(data_2()), server = TRUE, selected = FALSE)
    updateSelectizeInput(session, inputId = "ridge_gene", choices = rownames(data_2()), server = TRUE, selected = FALSE)
    updateSelectizeInput(session, inputId = "expression_gene_choice", choices = rownames(data_2()), server = TRUE, selected = FALSE)
    updateSelectizeInput(session, inputId = "expression_gene_choice_2", choices = rownames(data_2()), server = TRUE, selected = FALSE)
  })

  # upload confirmation
  observeEvent(input$user_file,{
    shinyalert("Download complete. Please select which metadata to keep and press \"trim\" button.", type = "success")
  })

  active_idents = reactive({
    validate(need(!is.null(data_2()),"data 2 null"))
    tryCatch(
      {
        if(all.is.numeric(levels(Idents(data_2()))))
        {
          gc()
          return(sort(as.numeric(levels(Idents(data_2())))))
        }
        else
        {
          gc()
          return(sort(levels(Idents(data_2()))))
        }
      },
      warning = function(w)
      {
        message(w)
        return(NULL)
      },
      error = function(e)
      {
        message(e)
        return(NULL)
      }
    )
  })

  observe({
    print(paste("active idents[1] = ",active_idents()[1]))
  })

  observe({
    toggle(id= "example_df", condition = {"ex" %in%  input$file_type})
  })


  #####################################################################

  ##############  QUALITY CONTROL : BOXPLOT ###########################

  #####################################################################

  # boxplot tutorial
  observeEvent(input$boxplot_help, {
    shinyalert("Tutorial ", "First choose display option : \"RNA count\" for number of RNA/cell, \"RNA Feature\" for total different RNA/cell, \"Percent mitochondrial\" for mitochondiral RNA percentage.
               Group data by choosing specific metadata in your dataset.
               Try different color option for customization.
               Other customization features available under the plot box.", type = "info", size = "m")
  })

  # get option choosed for boxplot display
  boxplot_option_input = eventReactive(input$boxplot_go,{
    switch(input$boxplot_option,
           "RNA count" = "nCount_RNA",
           "RNA feature"  = "nFeature_RNA",
           "% mitochondrial" = "percent.mito"
    )
  })

  # display otpions for group by
  observe({
    validate(need(!is.null(data()),"data null"))
    updateSelectInput(session, inputId = "boxplot_group_by", choice = meta_data(), selected = FALSE)
  })

  boxplot_n_color = eventReactive(
    c(input$boxplot_group_by,input$boxplot_idents,input$boxplot_idents_activation),{
    validate(need(!is.null(data_2()),"data 2 null"))
    n_color(input$boxplot_group_by,input$boxplot_idents,input$boxplot_idents_activation,data_2())
  })

  observe({
    updateSelectizeInput(session,inputId = "boxplot_colors", options = list(maxItems = boxplot_n_color()))
  })

  observe({
    toggle("boxplot_colors", condition = "custom" %in% input$boxplot_custom_colors)
  })

  boxplot_col_options = reactive({
    col_options(input$boxplot_custom_colors,input$boxplot_colors,boxplot_n_color())
  })

  observe({
    if(!input$boxplot_idents_activation)
    {
      updateSelectInput(session, inputId = "boxplot_idents", choices = active_idents(), selected = FALSE)
    }
    updateSelectInput(session,"boxplot_idents", choices = active_idents())
  })

  observe({
    toggle("boxplot_idents", condition = input$boxplot_idents_activation)
  })

  # drawing plot
  draw_boxplot = eventReactive(
    c(input$boxplot_go,input$boxplot_hide_dots,input$boxplot_add_box,input$boxplot_white),{
      p = VlnPlot(data_2(),
                  boxplot_option_input(),
                  group.by = input$boxplot_group_by,
                  idents = input$boxplot_idents,
                  cols = boxplot_col_options()) +
        ylab(input$boxplot_option) +
        ggtitle(paste(input$boxplot_option,"among",input$boxplot_group_by)) +
        xlab("group") +
        if(input$boxplot_add_box)
        {
          geom_boxplot()
        }
      else
      {
        NULL
      }
      if(input$boxplot_hide_dots)
      {
        p$layers[[2]] = NULL
      }
      return(p)
    })

  # plot output
  output$boxplot_ui = renderUI({
    plotOutput("boxplot", width = paste0(input$boxplot_width, "%"), height = input$boxplot_height)
  })

  output$boxplot = renderPlot({
    draw_boxplot()
  })

  # download plot
  output$down_boxplot = downloadHandler(
    filename = "boxplot.png",
    content = function(file){
      png(file,width = 1000,height = 1000)
      print(draw_boxplot())
      dev.off()
    }
  )

  #####################################################################

  ##############  QUALITY CONTROL : BARPLOT ###########################

  #####################################################################

  # barplot tutorial
  observeEvent(input$barplot_help, {
    shinyalert("Tutorial ", "Choose which metadata to group by (x axis).
               Choose subgroup in your metadata.
               Switch between percent of cells in each group or total cells.
               Other customization features available under the plot box.",
               type = "info", size = "m")
  })

  # display reactive values
  observe({
    validate(need(!is.null(data_2()),"data 2 null"))
    updateSelectInput(session, inputId = "barplot_meta", choice = meta_data(), selected = FALSE)
    updateSelectInput(session, inputId = "barplot_group_by", choice = meta_data(), selected = FALSE)
  })

  barplot_n_color = eventReactive(input$barplot_meta,{
    validate(need(!is.null(data_2()),"data 2 null"))
    n_color(input$barplot_meta,0,F,data_2())
  })

  observe({
    updateSelectizeInput(session,inputId = "barplot_colors", options = list(maxItems = barplot_n_color()))
  })

  observe({
    toggle("barplot_colors", condition = "custom" %in% input$barplot_custom_colors)
  })

  barplot_col_options = reactive({
    col_options(input$barplot_custom_colors,input$barplot_colors,barplot_n_color())
  })

  # drawing plot
  draw_barplot = eventReactive(input$barplot_go,{
    if("default" %in% input$barplot_custom_colors)
    {
      dittoBarPlot(data(), input$barplot_meta, group.by = input$barplot_group_by, scale = input$barplot_scale) + theme(axis.text.x = element_text(angle = 70, hjust=1))
    }
    else
    {
      dittoBarPlot(data(), input$barplot_meta, group.by = input$barplot_group_by, scale = input$barplot_scale, color.panel = barplot_col_options()) + theme(axis.text.x = element_text(angle = 70, hjust=1))
    }
  })

  # plot output
  output$barplot_ui = renderUI({
    plotOutput("barplot", width = paste0(input$barplot_width, "%"), height = input$barplot_height)
  })

  output$barplot = renderPlot({
    draw_barplot()
  })

  # download plot
  output$down_barplot = downloadHandler(
    filename = "barplot",
    content = function(file){
      png(file,width = 1000,height = 1000)
      print(draw_barplot())
      dev.off()
    }
  )

  #####################################################################

  #############  GENE EXPRESSION : VIOLIN STATS / T.TEST   ############

  #####################################################################

  # stat tutorial
  observeEvent(input$gene_expression_help, {
    shinyalert("Tutorial ", "Select a gene by writing it's name (auto-completion).
               Group data by choosing specific metadata in your dataset.
               Check \"isolate specific classes\" to remove some classes on the plot.
               Try different color option for customization.
               Check \"stats\" box to add student tests.
               Warning : if you don't see stats line on your plot, consider change y axis limits
               Other customization features available under the plot box.",
               type = "info", size = "m")
  })

  # display reactive values
  observe({
    validate(need(!is.null(data()),"data null"))
    updateSelectInput(session, inputId = "expression_group_by_choice", choice = meta_data(), selected = FALSE)
  })

  observe({
    updateSelectInput(session, inputId = "student_idents", choices = active_idents())
    if(!input$student_idents_activation)
    {
      updateSelectInput(session, inputId = "student_idents", choices = active_idents(), selected = FALSE)
    }
  })

  observe({
    toggle("student_stats", condition = input$student_stats_activation)
  })

  observe({
    toggle(id = "student_idents", condition = input$student_idents_activation)
  })

  # display either p-value or ***
  format = eventReactive(input$p_format,{
    if(input$p_format){
      return("p.format")
    }else{
      return("p.signif")
    }
  })

  observe({
    toggle(id = "expression_colors", condition = "custom" %in% input$expression_custom_color)
  })

  expression_n_color = eventReactive(
    c(input$expression_group_by_choice,input$student_idents,input$student_idents_activation),{
      validate(need(!is.null(data_2()),"data 2 null"))
      n_color(input$expression_group_by_choice,input$student_idents,input$student_idents_activation,data_2())
    })

  observe({
    updateSelectizeInput(session,inputId = "expression_colors", options = list(maxItems = expression_n_color()))
  })

  observe({
    toggle("expression_colors", condition = "custom" %in% input$expression_custom_colors)
  })

  expression_col_options = reactive({
    col_options(input$expression_custom_colors,input$expression_colors,expression_n_color())
  })

  student_stats_list = reactive({
    if(input$student_idents_activation && length(input$student_idents)>1)
    {
      return(get_pairs(input$student_idents))
    }
    else
    {
      return(get_pairs(active_idents()))
    }
  })

  observe({
    validate(need(!is.null(student_stats_list()),"student list null"))
    updateSelectizeInput(session,"student_stats", choices = student_stats_list())
  })

  student_x_lab = reactive({
    if(input$student_idents_activation)
    {
      return(input$student_idents)
    }
    else
    {
      return(active_idents())
    }
  })

  # make plot + stats
  draw_gene_expression_plot = eventReactive(
    c(input$gene_expression_go,input$gene_expression_y_lim,input$gene_expression_hide_dots,input$p_format,input$expression_add_box),{
      tryCatch(
        {
          p = VlnPlot(data_2(), features = input$expression_gene_choice,
                      pt.size = 0.1,
                      group.by = input$expression_group_by_choice,
                      y.max = input$gene_expression_y_lim,
                      idents = input$student_idents,
                      cols = expression_col_options()) +
            stat_compare_means(step.increase = 0.15,
                               method = "t.test",
                               label = format(),
                               comparisons = str_split(input$student_stats, " VS "))+
            scale_x_discrete(limits = factor(student_x_lab()))+
            if(input$expression_add_box)
            {
              geom_boxplot()
            }
          else
          {
            NULL
          }
          if(input$gene_expression_hide_dots)
          {
            p$layers[[2]] = NULL
          }
          return(p)
        },
        error = function(e)
        {
          shinyalert(message(e), type= "error")
          return(NULL)
        },
        warning = function(w)
        {
          shinyalert(message(w), type= "error")
          return(NULL)
        }

      )
    })

  # render plot
  output$gene_expression_ui = renderUI({
    plotOutput("gene_expression", width = paste0(input$gene_expression_width, "%"), height = input$gene_expression_height)
  })

  output$gene_expression = renderPlot({
    validate(
      need(input$gene_expression_go, "")
    )
    draw_gene_expression_plot()
  })

  # download plot
  output$down_gene_expression = downloadHandler(
    filename = "gene_expression_plot.png",
    content = function(file){
      png(file,width = 1000,height = 1000)
      print(draw_gene_expression_plot())
      dev.off()
    }
  )

  #####################################################################

  #############  GENE EXPRESSION : VIOLIN STATS / ANOVA   ############

  #####################################################################

  # Anova tutotial
  observeEvent(input$gene_expression_help_2, {
    shinyalert("Tutorial ", "Select a gene by writing it's name (auto-completion).
               Group data by choosing specific metadata in your dataset.
               Check \"split by\" box to subdvide youre data.
               Check \"isolate specific classes\" to remove some classes on the plot.
               Try different color option for customization.
               Other customization features available under the plot box.",
               type = "info", size = "m")
  })

  # display reactive value
  observe({
    validate(need(!is.null(data()),"data null"))
    updateSelectInput(session, inputId = "expression_group_by_choice_2",  choice = meta_data(), selected = FALSE)
  })

  observe({
    validate(need(!is.null(data()),"data null"))
    updateSelectInput(session, inputId = "expression_split_by_choice_2", choice = meta_data(), selected = FALSE)
  })

  observe({
    validate(need(!is.null(data()),"data null"))
    updateSelectInput(session,inputId = "anova_idents", choices = active_idents(), selected = FALSE)
  })

  # hide or show split by list
  observe({
    toggle(id = "expression_split_by_choice_2", condition = input$split_by_anova_activation)
    toggle(id = "anova_idents", condition = input$anova_idents_activation)
  })

  # return user's choice for split by otpion
  anova_split_by_option = eventReactive(
    c(input$expression_split_by_choice_2,input$split_by_anova_activation),{
      if(input$split_by_anova_activation == FALSE){
        return(NULL)
      }
      else{
        return(input$expression_split_by_choice_2)
      }
    })

  # return user's choice for idents otpion
  anova_idents_option = eventReactive(
    c(input$anova_idents,input$anova_idents_activation),{
      if(!input$anova_idents_activation){
        return(NULL)
      }else{
        return(input$anova_idents)
      }
    })

  observe({
    toggle(id = "expression_colors_2", condition = "custom" %in% input$expression_custom_color_2)
  })

  expression_n_color_2 = eventReactive(
    c(input$expression_group_by_choice_2,input$anova_idents,input$anova_idents_activation),{
      validate(need(!is.null(data_2()),"data 2 null"))
      n_color(input$expression_group_by_choice_2,input$anova_idents,input$anova_idents_activation,data_2())
    })

  observe({
    updateSelectizeInput(session,inputId = "expression_colors_2", options = list(maxItems = expression_n_color_2()))
  })

  observe({
    toggle("expression_colors_2", condition = "custom" %in% input$expression_custom_colors_2)
  })

  expression_col_options_2 = reactive({
    col_options(input$expression_custom_colors_2,input$expression_colors_2,expression_n_color_2())
  })

  # draw plot + stats
  draw_anova = eventReactive(
    c(input$gene_expression_go_2,input$gene_expression_y_lim_2,input$gene_expression_hide_dots_2,input$expression2_add_box),{
      p = VlnPlot(data_2(), features = input$expression_gene_choice_2,
                  pt.size = 0.1,
                  group.by = input$expression_group_by_choice_2,
                  split.by = anova_split_by_option(),
                  idents = anova_idents_option(),
                  y.max = input$gene_expression_y_lim_2,
                  cols = expression_col_options_2()
      ) + stat_compare_means(method = "anova", label.y = input$gene_expression_y_lim_2, label.x = 2) +
        if(input$expression2_add_box)
        {
          geom_boxplot()
        }
      else
      {
        NULL
      }
      if(input$gene_expression_hide_dots_2)
      {
        p$layers[[2]] = NULL
      }
      return(p)
    })

  # display plot
  output$gene_expression_ui_2 = renderUI({
    plotOutput("gene_expression_2", width = paste0(input$gene_expression_width_2, "%"), height = input$gene_expression_height_2)
  })

  output$gene_expression_2 = renderPlot({
    validate(
      need(input$gene_expression_go_2, "")
    )
    draw_anova()
  })

  # download plot
  output$down_gene_expression_2 = downloadHandler(
    filename = "anova.png",
    content = function(file){
      png(file,width = 1000,height = 1000)
      print(draw_anova())
      dev.off()
    }
  )

  #####################################################################

  ##############  GENE EXPRESSION : RIDGE  ############################

  #####################################################################

  # ridge tutorial
  observeEvent(input$ridge_help, {
    shinyalert("Tutorial ", "Select a gene by writing it's name (auto-completion).
               Group data by choosing specific metadata in your dataset.
               Check \"isolate specific classes\" to remove some classes on the plot.
               Try different color option for customization.
               Other customization features available under the plot box.",  type = "info", size = "m")
  })

  # display reactive value
  observe({
    validate(need(!is.null(data()),"data null"))
    updateSelectInput(session, inputId = "ridge_group_by", choice = meta_data(), selected = FALSE)
  })

  observe({
    validate(need(!is.null(data()),"data null"))
    updateSelectInput(session,inputId = "ridge_idents", choices = active_idents())
  })

  # hide or show options
  observe({
    toggle(id = "ridge_idents", condition = input$ridge_idents_activation)
  })

  # return user's choice for idents otpion
  ridge_idents_option = eventReactive(
    c(input$ridge_idents,input$ridge_idents_activation),{
      if(input$ridge_idents_activation == FALSE){
        return(NULL)
      }else{
        return(input$ridge_idents)
      }
    })

  observe({
    toggle(id = "ridge_colors", condition = "custom" %in% input$ridge_custom_color)
  })

  ridge_n_color = eventReactive(
    c(input$ridge_group_by,input$ridge_idents,input$ridge_idents_activation),{
      validate(need(!is.null(data_2()),"data 2 null"))
      n_color(input$ridge_group_by,input$ridge_idents,input$ridge_idents_activation,data_2())
  })

  observe({
    updateSelectizeInput(session,inputId = "ridge_colors", options = list(maxItems = ridge_n_color()))
  })

  observe({
    toggle("ridge_colors", condition = "custom" %in% input$ridge_custom_colors)
  })

  ridge_col_options = reactive({
    col_options(input$ridge_custom_colors,input$ridge_colors,ridge_n_color())
  })

  # draw plot
  draw_ridge_plot = eventReactive(input$ridge_go,{
    RidgePlot(data_2(), feature = input$ridge_gene, group.by = input$ridge_group_by, idents = ridge_idents_option(),cols = ridge_col_options())
  })

  # display plot
  output$ridge_plot_ui = renderUI({
    plotOutput("ridge_plot", width = paste0(input$ridge_width, "%"), height = input$ridge_height)
  })

  output$ridge_plot = renderPlot({
    draw_ridge_plot()
  })

  # download plot
  output$down_ridge = downloadHandler(
    filename = "ridge_plot.png",
    content = function(file){
      png(file,width = 1000,height = 1000)
      print(draw_ridge_plot())
      dev.off()
    }
  )

  #####################################################################

  ##############  GENE EXPRESSION : DOT PLOT  #########################

  #####################################################################

  # dotplot tutorial
  observeEvent(input$dotplot_help, {
    shinyalert("Tutorial ", "Select some genes by writing their names (auto-completion).
               Group data by choosing specific metadata in your dataset.
               Check \"isolate specific classes\" to remove some classes on the plot.
               Other customization features available under the plot box.",  type = "info", size = "m")
  })

  # display reactive value
  observe({
    validate(need(!is.null(data_2()),"data 2 null"))
    updateSelectInput(session, inputId = "dot_plot_group_by", choice = meta_data(), selected = FALSE)
    updateSelectInput(session, "dot_plot_idents", choices = active_idents())
  })

  observe({
    toggle(id = "dot_plot_idents", condition = input$dot_plot_idents_activation)
  })

  # get gene list from user choice
  dot_plot_gene_list = eventReactive(input$dot_plot_genes,{
    list = c()
    list <- unlist(strsplit(input$dot_plot_genes, ","))
    return(list)
  })

  # draw plot
  draw_dot_plot = eventReactive(input$dot_plot_go,{
    DotPlot(data_2(), features =  dot_plot_gene_list() ,group.by = input$dot_plot_group_by, idents = input$dot_plot_idents)
  })

  # display plot
  output$dot_plot_ui = renderUI({
    plotOutput("dot_plot", width = paste0(input$dot_plot_width, "%"), height = input$dot_plot_height)
  })

  output$dot_plot = renderPlot({
    draw_dot_plot()
  })

  # download plot
  output$down_dot_plot = downloadHandler(
    filename = "dot_plot.png",
    content = function(file){
      png(file,width = 1000,height = 1000)
      print(draw_dot_plot())
      dev.off()
    }
  )

  #####################################################################

  ##########################   GENESETS  ##############################

  #####################################################################

  # dotplot tutorial
  observeEvent(input$gs_help, {
    shinyalert("Tutorial ", "Choose the example geneset provided by Hallmark (general biological processes) or uplod youre own.
               Choose best options for your dataset if you upload one (right extension and separator), and verify on the example panel that everything is ok.
               Press \"create gene list\" button when everything is ready (could take a while).
               \"Idents\" selection allow to display which metadata you want to analyze.
               Then select two condition to compare (\"condition 1\" and \"condition 2\").
               \"sub category\" selection allow to keep only some value (eg: NK cells and B cells in the \"keep\" box).
               When youre ready press \"run\", the calculation could take a while, specially if you have provided a big geneset.
               Please wait until white loading circle on the left bottom part dissapear.
               Other customization features available under the plot box.",  type = "info", size = "m")
  })

  observe({
    shinyFileChoose(input, "gs_user_file_path", roots = volumes, session = session)
    if(!is.null(input$gs_user_file_path)){
      gs_file_selected <<-parseFilePaths(volumes, input$gs_user_file_path)
    }
  })

  observe({
    toggle("gs_user_file_path", condition = "upload" %in% input$gs_file_type)
    toggle("gs_file_extension", condition = "upload" %in% input$gs_file_type)
    toggle("gs_file_sep", condition = "upload" %in% input$gs_file_type)
  })

  observeEvent(
    c(input$gs_file_type,input$gs_user_file_path),{
      if("upload" %in% input$gs_file_type)
      {
        if(length(gs_file_selected$datapath) != 0)
        {
          ext = tolower(tools::file_ext(gs_file_selected$datapath))
          if (ext != "xls" & ext != "xlsx" & ext != "csv" & ext != "ods" & ext !="txt" & ext !="rds")
          {
            disable("gs_go")
            shinyalert("Unsupported file extension", type = "error")
          }
          else if(ext != input$gs_file_extension)
          {
            disable("gs_go")
            shinyalert("Conflict between extension choose and file uploaded", type = "warning")
          }
          else
          {
            enable("gs_go")
            shinyalert("File ready", type = "success")
          }
        }
        else
        {
          disable("gs_go")
        }
      }
      else if ("example" %in% input$gs_file_type)
      {
        enable("gs_go")
      }
    })

  gs_df = eventReactive(
    c(input$gs_file_type,input$gs_user_file_path,input$gs_file_sep),{
      if("upload" %in% input$gs_file_type)
      {
        if (length(gs_file_selected$datapath) != 0)
        {
          ext_2 = tolower(tools::file_ext(gs_file_selected$datapath))
          separator = input$gs_file_sep
          if("tab" %in% input$gs_file_sep)
          {
            separator = ""
          }
          if(ext_2 == "csv" | ext_2 == "txt")
          {
            df = read.csv(gs_file_selected$datapath, sep = separator, fileEncoding =  "UTF-8-BOM")
          }
          else if(ext_2 == "rds")
          {
            df = readRDS(gs_file_selected$datapath)
          }
          else if(ext_2 == "xls" | ext_2 == "xlsx")
          {
            df = read.table(gs_file_selected$datapath,sep = "\t", header=TRUE)
          }
          else
          {
            shinyalert("Cannot read file", type = "error")
          }
        }
      }
      else
      {
        df = read.csv("./data/geneset_example.csv", sep = ",", fileEncoding =  "UTF-8-BOM")
      }
      return(df)
    })

  observeEvent(
    c(input$gs_file_type,input$gs_user_file_path,input$gs_file_sep),{
      print(paste("lenght path = ",length(gs_file_selected$datapath)))
      toggle("gs_table", condition =  length(gs_file_selected$datapath) != 0)
      toggle("gs_table_box", condition =  length(gs_file_selected$datapath) != 0)
    })

  output$gs_table = renderTable({
    if(length(colnames(gs_df()))> 5)
    {
      head(gs_df()[,1:5])
    }
    else
    {
      head(gs_df())
    }
  })

  geneset_h = eventReactive(input$gs_go,{
    tryCatch(
      {
        if("upload" %in% input$gs_file_type)
        {
          genesets <<- gs_df()
          df = create_list(genesets)
          shinyalert("Gene list created", type = "success")
          gc()
          return(df)
        }
        else
        {
          df = readRDS("./data/geneset_h_example.RData")
          shinyalert("Gene list created", type = "success")
          return(df)
        }

      },
      error = function(e)
      {
        shinyalert("Something went wrong
                    Press help for datatype info", type = "error")
        shiny:::reactiveStop(e)
      }
    )
  })

  observe({
    validate(need(!is.null(data_2()),"data 2 null"))
    updateSelectInput(session,inputId = "gs_cat", choice = meta_data(), selected = FALSE)
    updateSelectInput(session,inputId = "gs_sub_cat", choice = meta_data(), selected = FALSE)
  })

  observe({
    updateSelectInput(session,inputId = "gs_cond_1", choices = active_idents())
    updateSelectInput(session,inputId = "gs_cond_2", choices = active_idents())
  })

  observe({
    validate(need(input$gs_sub_cat != "","gs sub cat null"))
    updateSelectInput(session,inputId = "gs_sub_cat_choices", choices = data_2()[[input$gs_sub_cat]])
  })

  observe({
    validate(need(!is.null(gs_df()),"gs df null"))
    updateSelectInput(session,inputId = "gs_plot_meta", choices = colnames(gs_df()))
  })

  observe({
    validate(need(!is.null(gs_stats_list()),"gs list null"))
    updateSelectizeInput(session,"gs_stats", choices = gs_stats_list())
  })

  observe({
    validate(need(!is.null(data()),"data null"))
    updateSelectizeInput(session, inputId = "gs_idents", choice = meta_data(), selected = FALSE)
  })

  gs_cats = reactive({
    if(!is.null(data_2()))
    {
      expr = sort(unique(eval(parse(text = paste0("data_2()$",input$gs_cat)))))
      return(expr)
    }
    else
    {
      return(NULL)
    }
  })

  meta_list = eventReactive(input$gs_go_2,{
    tryCatch(
      {
        cluster_ident = FALSE
        df = data_2()
        my_idents = levels(Idents(df))
        for(i in 1:length(my_idents))
        {
          if(grepl("^[[:digit:]]+", my_idents[i]))
          {
            levels(Idents(df)) = paste0("n_",levels(Idents(df)))
            cluster_ident = TRUE
            break
          }
        }
        cond_1 = input$gs_cond_1
        cond_2 = input$gs_cond_2
        cond_1 = gsub(" ","_",cond_1)
        cond_2 = gsub(" ","_",cond_2)
        assign(cond_1,subset(df, idents = eval(parse(text = paste0("\"",cond_1,"\"")))), envir = .GlobalEnv)
        assign(cond_2,subset(df, idents = eval(parse(text = paste0("\"",cond_2,"\"")))), envir = .GlobalEnv)
        sub_cat = input$gs_sub_cat
        expr = paste0("Idents(",cond_1,")=",cond_1,"$",sub_cat)
        eval(parse(text = expr))
        expr2 = paste0("Idents(",cond_2,")=",cond_2,"$",sub_cat)
        eval(parse(text = expr2))
        list_sub_cond = input$gs_sub_cat_choices
        meta_list = c()
        for (i in 1:length(list_sub_cond)) {
          name_cond_1 = paste0(list_sub_cond[i],"_",cond_1)
          name_cond_2 = paste0(list_sub_cond[i],"_",cond_2)
          name_cond_1 = gsub(" ","_",name_cond_1)
          name_cond_2 = gsub(" ","_",name_cond_2)
          assign(name_cond_1,subset(eval(parse(text = cond_1)), idents = list_sub_cond[i]), envir = .GlobalEnv)
          assign(name_cond_2,subset(eval(parse(text = cond_2)), idents = list_sub_cond[i]), envir = .GlobalEnv)
          meta_list = append(meta_list,name_cond_1)
          meta_list = append(meta_list,name_cond_2)
        }
        if(cluster_ident)
        {
          levels(Idents(df)) = gsub("n_","",levels(Idents(df)))
        }
        gc()
        return(meta_list)
      }
      ,
      error = function(e)
      {
        shinyalert(paste("Error: ",e), type = "error")
        return(NULL)
      }
    )
  })

  list_briscoe = eventReactive(input$gs_go_2,{
    validate(need(!is.null(meta_list()),"meta list null"))
    tryCatch(
      {
        print(paste("meta list = ",meta_list()))
        add_module_score(meta_list(),geneset_h())
        my_list = regroup(gs_df(),meta_list())
        gc()
        return(my_list)
      },
      error = function(e)
      {
        shinyalert("Error: genes not found in your dataframe. Consider using an ohter dataframe as reference", type = "error")
        return(NULL)
      }
    )
  })

  observe({
    print(paste("length gs_df = ",length(gs_df())))
    print(paste("length geneset = ",length(geneset_h())))
    print(paste("list briscoe = ",length(list_briscoe())))
  })

  observe({
    toggle("gs_file", condition = input$gs_file_type != "example")
    toggle("gs_conversion_from", condition = input$gs_conversion == TRUE)
    toggle("gs_conversion_to", condition = input$gs_conversion == TRUE)
    toggle("gs_go2", condition = input$gs_conversion == TRUE)
    toggle("geneset_ref_specie", condition = input$gs_file_type != "example")
  })

  gs_plot_boxes = eventReactive(input$gs_plot_add_boxes,{
    if(input$gs_plot_add_boxes)
    {
      return("boxplot")
    }
    else
    {
      return(NULL)
    }
  })

  gs_plot_legend = eventReactive(input$gs_plot_legend,{
    if(input$gs_plot_legend)
    {
      return("top")
    }
    else
    {
      return("none")
    }
  })

  gs_plot_angle = eventReactive(input$gs_plot_angle,{
    if(input$gs_plot_angle)
    {
      return(45)
    }
    else
    {
      return(0)
    }
  })

  gs_stats_list = reactive({
    validate(need(!is.null(meta_list()),"meta list null"))
    return(get_pairs(meta_list()))
  })

  align_step = reactive({
    if(input$gs_align_bar)
    {
      return(0)
    }
    else
    {
      return(0.15)
    }
  })

  geneset_n_color = reactive({
    return(length(meta_list()))
  })

  observe({
    updateSelectizeInput(session,inputId = "geneset_colors", options = list(maxItems = geneset_n_color()))
  })

  observe({
    toggle("geneset_colors", condition = "custom" %in% input$geneset_custom_colors)
  })

  geneset_col_options = reactive({
    col_options(input$geneset_custom_colors,input$geneset_colors,geneset_n_color())
  })

  # draw plot
  draw_gs_plot = eventReactive(
    c(input$gs_go_2,input$gs_go_3,input$gs_plot_angle,input$gs_plot_legend,input$gs_plot_add_boxes,input$gs_plot_meta,align_step(),input$gs_stats),{
      my_labels = meta_list()
      my_labels = gsub("_"," ",my_labels)
      p<-ggviolin(list_briscoe(),
                  x = paste0("V",length(list_briscoe())),
                  y = input$gs_plot_meta,
                  fill = paste0("V",length(list_briscoe())),
                  add = gs_plot_boxes(),
                  add.params = list(color = "black",
                                    fill = "white",
                                    width = 0.05),
                  ylab = input$gs_plot_meta,
                  xlab = "group",
                  palette = geneset_col_options())+
        theme(axis.text.x = element_text(angle = gs_plot_angle(), vjust = 0.5, hjust=0.5))+
        theme(legend.position = gs_plot_legend())+
        scale_x_discrete(labels= my_labels)+
        stat_compare_means(step.increase = align_step(),
                           symnum.args = list(cutpoints = c(0, 0.0001,  0.001, 0.01, 0.05 ,1),
                                              symbols = c("****", "***", "**", "*", "ns")),
                           method = "t.test",
                           comparisons = str_split(input$gs_stats, " VS "),
                           label.y = 0.5)
      return(p)
    })

  # display plot
  output$gs_plot_ui = renderUI({
    plotOutput("gs_plot", width = paste0(input$gs_width, "%"), height = input$gs_height)
  })

  output$gs_plot = renderPlot({
    draw_gs_plot()
  })

  #####################################################################

  ##########################    REDUCTIONS  ###########################

  #####################################################################

  # cluster tutorial
  observeEvent(input$cluster_help, {
    shinyalert("Tutorial ", "If you have reduction calculated in youre dataset, choose one in \"reduction\" box, else, check \"calculate\" to let the program add reductions for you.
               You have to choose the number of reductions for each calculation, values have been setted by R default choices.
               ", type = "info", size = "m")
  })

  data_3 = eventReactive(input$reduction_go,{
    df = data_2()
    Idents(df) = df$seurat_clusters
    DefaultAssay(df) <- "RNA"
    df <- NormalizeData(object = df, normalization.method = "LogNormalize",  scale.factor = 10000)
    DefaultAssay(df) <- "integrated"
    df = ScaleData(df)
    df <- RunPCA(df, features = VariableFeatures(object = df),npcs = as.numeric(input$pca_reductions))
    print(paste("lenght reduction = ",length(df@reductions)))
    df <- RunUMAP(df, dims = 1:as.numeric(input$umap_reductions))
    print(paste("lenght reduction = ",length(df@reductions)))
    df <- RunTSNE(df, dims = 1:as.numeric(input$tsne_reductions))
    shinyalert("reduction calculated", type = "success")
    return(df)
  })

  observe({
    validate(need(!is.null(data_3()),"data 3 null"))
    print(paste("data 3 reductions = ",names(data_3()@reductions)))
  })

  # display reactive value
  observe({
    validate(need(!is.null(data_2()),"data 2 null"))
    updateSelectInput(session, inputId = "cluster_group_by", choice = meta_data(), selected = FALSE)
  })

  # get option choosed for culster label
  cluster_label_input = eventReactive(input$cluster_go,{
    switch(input$cluster_labels,
           "yes" = "TRUE",
           "no"  = "FALSE"
    )
  })

  observe({
    validate(need(length(data_2()@reductions)>0,"no reductions"))
    updateSelectInput(session,"cluster_reduction", choices = names(data_2()@reductions))
    validate(need(length(data_3()@reductions)>0,"no reductions"))
    updateSelectInput(session,"cluster_reduction", choices = names(data_3()@reductions))
  })

  observe({
    toggle("pca_reductions", condition = input$reduction_activation)
    toggle("tsne_reductions", condition = input$reduction_activation)
    toggle("umap_reductions", condition = input$reduction_activation)
    toggle("reduction_go", condition = input$reduction_activation)
  })

  # plot creation function
  draw_clusters = eventReactive(input$cluster_go,{
    print("enter")
    if(length(data_2()@reductions)>0)
    {
      print("data 2")
      DimPlot(data_2(), reduction = input$cluster_reduction, group.by = input$cluster_group_by, label = cluster_label_input())
    }
    else if(length(data_3()@reductions)>0)
    {
      print("data_3")
      DimPlot(data_3(), reduction = input$cluster_reduction, group.by = input$cluster_group_by, label = cluster_label_input())
    }
    else
    {
      print("error, no reduction found")
    }
  })

  # display plot
  output$Clusters_plot_ui = renderUI({
    plotOutput("cluster_plot", width = paste0(input$cluster_width, "%"), height = input$cluster_height)
  })

  output$cluster_plot = renderPlot({
    draw_clusters()
  })

  # download plot
  output$down_clusters = downloadHandler(
    filename = "Visualization.png",
    content = function(file){
      png(file, width = 1000, height = 1000)
      print(draw_clusters())
      dev.off()
    }
  )

  #####################################################################

  ######################     DENSITY PLOT  ############################

  #####################################################################

  # density tutorial
  observeEvent(input$density_help, {
    shinyalert("Tutorial ", "Select a gene by writing it's name (auto-completion)
               Add other genes to perform coexpression analysis", type = "info", size = "m")
  })

  # display multiple gene box
  observe({
    toggle(id = "density_gene_2", condition = input$density_gene_2_activation)
    toggle(id = "density_gene_3", condition = input$density_gene_3_activation)
    toggle(id = "density_gene_4", condition = input$density_gene_4_activation)
    toggle(id = "density_gene_3_activation", condition = input$density_gene_2_activation)
    toggle(id = "density_gene_4_activation", condition = input$density_gene_3_activation)
  })

  # plot creation function
  draw_density = eventReactive(input$density_go,{
    if(!input$density_gene_2_activation){
      FeaturePlot(data(), features = input$density_gene_1, cols = c("green", "blue"))
    }
    else{
      if(!input$density_gene_3_activation){
        p1 = FeaturePlot(data(), features = input$density_gene_1, cols = c("green", "blue"))
        p2 = FeaturePlot(data(), features = input$density_gene_2, cols = c("green", "blue"))
        plot_grid(p1,p2)
      }
      else{
        if(!input$density_gene_4_activation){
          p1 = FeaturePlot(data(), features = input$density_gene_1, cols = c("green", "blue"))
          p2 = FeaturePlot(data(), features = input$density_gene_2, cols = c("green", "blue"))
          p3 = FeaturePlot(data(), features = input$density_gene_3, cols = c("green", "blue"))
          plot_grid(p1,p2,p3)
        }
        else{
          p1 = FeaturePlot(data(), features = input$density_gene_1, cols = c("green", "blue"))
          p2 = FeaturePlot(data(), features = input$density_gene_2, cols = c("green", "blue"))
          p3 = FeaturePlot(data(), features = input$density_gene_3, cols = c("green", "blue"))
          p4 = FeaturePlot(data(), features = input$density_gene_4, cols = c("green", "blue"))
          plot_grid(p1,p2,p3,p4)
        }
      }
    }
  })

  # display plot
  output$Density_plot_ui = renderUI({
    plotOutput("Density_plot", width = paste0(input$density_width, "%"), height = input$density_height)
  })

  output$Density_plot = renderPlot({
    draw_density()
  })

  # download plot
  output$down_density = downloadHandler(
    filename = "density_plot.png",
    content = function(file){
      png(file, width = 1000, height = 1000)
      print(draw_density())
      dev.off()
    }
  )

  #####################################################################

  ######################     DENSITY PLOT  2  #########################

  #####################################################################

  # density tutorial
  observeEvent(input$density_help_2, {
    shinyalert("Tutorial ", "Select two genes by writing their names (auto-completion)", type = "info", size = "m")
  })

  # display gene list
  observe({
    validate(need(!is.null(data_2()),"data 2 null"))
    updateSelectizeInput(session, inputId = "density_gene_1_2", choices = rownames(data_2()), server = TRUE,selected = FALSE)
    updateSelectizeInput(session, inputId = "density_gene_2_2", choices = rownames(data_2()), server = TRUE,selected = FALSE)
  })

  # plot creation function
  draw_density_2 = eventReactive(input$density_go_2,{
    FeaturePlot(data(), features = c(input$density_gene_1_2,input$density_gene_2_2), cols = c("green", "blue"), blend = TRUE)
  })

  # display plot
  output$Density_plot_ui_2 = renderUI({
    plotOutput("Density_plot_2", width = paste0(input$density_width_2, "%"), height = input$density_height_2)
  })

  output$Density_plot_2 = renderPlot({
    draw_density_2()
  })

  # download plot
  output$down_density_2 = downloadHandler(
    filename = "density_plot.png",
    content = function(file){
      png(file, width = 1000, height = 1000)
      print(draw_density_2())
      dev.off()
    }
  )

  #####################################################################

  ###############    DIFFERENTIAL CALCULATION  ########################

  #####################################################################

  # differential calulation tutorial
  observeEvent(input$diff_help, {
    shinyalert("Tutorial ", "Differential analysis helps to understand differential gene expression between two conditions.
               \"Idents box helps to determine which condition to bu used for calculation.
               \"Method box allow to switch between 3 methods, regarding of the type of analysis (1 condition against all the other, multiple condition against multiple conditon or each condition against each other conditon.
               Regarding which method you have choossed, you choose group 1 or group 1 and 2 for comparision.
               \"Metadata to look up\" allow to keep only one cellular type. Please choose the cellular prediction in this box to keep one cellular type in the following box \"keep\".
               The log fold change treshold determine the precision. The lower is the treshold, the more differential genes you will found.
               FInally, the p-value treshold help to discrimine gene under this limit.
               ", type = "info", size = "m")
  })

  diff_keep_options = reactive({
    validate(need(!is.null(data_2()),"data 2 null"))
    if(input$diff_meta == "")
    {
      return(NULL)
    }
    else
    {
      my_text = paste0("unique(data_2()$",input$diff_meta,")")
      return(eval(parse(text = my_text)))
    }
  })

  observe({
    updateSelectInput(session, "diff_idents", choices = meta_data())
  })

  # display reactive value
  observe({
    updateSelectInput(session,inputId = "diff_calcul_cluster_1", choices = active_idents())
    updateSelectInput(session,inputId = "diff_calcul_cluster_1_2",choices = active_idents())
    updateSelectInput(session,inputId = "diff_calcul_cluster_2", choices = active_idents())
  })

  observe({
    updateSelectInput(session, inputId = "diff_meta", choices = meta_data())
  })

  observe({
    updateSelectInput(session, inputId = "diff_keep", choices = diff_keep_options())
  })

  # display and hide reactive values
  observe({
    toggle(id = "diff_calcul_cluster_1", condition = {"1 VS all" %in%  input$diff_calcul_type})
    toggle(id = "diff_calcul_cluster_1_2", condition = {"multiple VS multiple" %in%  input$diff_calcul_type})
    toggle(id = "diff_calcul_cluster_2", condition = {"multiple VS multiple" %in%  input$diff_calcul_type})
    toggle(id = "diff_calcul_split_by", condition = input$diff_calcul_split_by_activation)
    toggle(id = "diff_calcul_split_by_choices", condition = input$diff_calcul_split_by_activation)
    toggle(id = "split", condition = input$diff_calcul_split_by_activation)
  })

  # calculation
  diff_genes = eventReactive(input$diff_calcul_go,{
    df = data_2()
    my_text = paste0("Idents(df) = df$",input$diff_meta)
    eval(parse(text = my_text))
    df = subset(df, idents = input$diff_keep)
    my_text2 = paste0("Idents(df) = df$",input$diff_idents)
    eval(parse(text = my_text2))
    if(input$diff_calcul_type == "each VS all" ){
      DE = FindAllMarkers(df, logfc.threshold = input$diff_calcul_fold_change)
      DE = DE[,-7]
      DE = DE[,-6]
    } else if (input$diff_calcul_type == "1 VS all"){
      DE = FindMarkers(df, ident.1 = input$diff_calcul_cluster_1, logfc.threshold = input$diff_calcul_fold_change,test.use = "MAST", min.pct  = 0.1,)
    } else {
      DE = FindMarkers(df, ident.1 = input$diff_calcul_cluster_1_2, ident.2 = input$diff_calcul_cluster_2, logfc.threshold = input$diff_calcul_fold_change,test.use = "MAST", min.pct  = 0.1, )
    }
    gc()
    return(DE)
  })

  # render differential table
  output$diff_calcul_table = renderTable({
    validate(need(!is.null(diff_genes()),"diff genes null"))
    head(diff_genes())
    }, rownames =  TRUE)

  # download table
  output$down_diff_table = downloadHandler(
    filename = "differential_table.csv",
    content = function(file){
      write.csv(diff_genes(),file)
    }
  )


  #####################################################################

  ######################    VOLCANO PLOT   ############################

  #####################################################################

  # volcano tutorial
  observeEvent(input$volcano_help, {
    shinyalert("Tutorial ", "Here, all the differential genes calculated befor are shown.
               You have to resize youre sample to obtain readable plot (next two panels.
               To resize, move the \"f-value cutoff\" slider. All genes inside those two vertical bars will be removed from the differential table.
               When you think you have the right amout of gene, press \"resize differential table\" button.
               ", type = "info", size = "m")
  })

  # volcano warning message
  output$warning_volcano = renderText("Please, go to Find differential gene tab first")

  # display and hide boxes
  observe({
    toggle(id = "warning_volcano_box", condition = !input$diff_calcul_go)
    toggle(id = "volcano_box", condition = input$diff_calcul_go)
  })

  # confirmation message for differential resize
  observeEvent(input$volcano_go, {
    shinyalert("Nice !", "You've set a new table with top differential genes, go check it on differential dot plot tab", type = "success")
  })

  # resize gene tab
  top_n_genes = eventReactive(
    c(input$volcano_go,input$diff_calcul_go),{
      if(input$volcano_go){
        DE = diff_genes()[which(abs(diff_genes()$avg_log2FC) > input$volcano_f_value),]
        return(DE)
      }else{
        DE = diff_genes()
        return(DE)
      }
    })

  # plot creation function
  draw_volcano_plot = eventReactive(
    c(input$diff_calcul_go,input$volcano_x_lim,input$volcano_f_value),{
      EnhancedVolcano(diff_genes(),
                      lab = rownames(diff_genes()),
                      x = 'avg_log2FC',
                      y = 'p_val_adj',
                      pCutoff = 0.05,
                      xlim =c(-input$volcano_x_lim,input$volcano_x_lim),
                      FCcutoff = input$volcano_f_value,
                      labSize = 4,
                      titleLabSize = 10,
                      subtitleLabSize = 0,
                      captionLabSize = 10,
                      axisLabSize =10
      )
    })

  # display plot
  output$volcano_plot_ui = renderUI({
    plotOutput("volcano_plot", width = paste0(input$volcano_width, "%"), height = input$volcano_height)
  })
  output$volcano_plot = renderPlot({
    draw_volcano_plot()
  })

  # download plot
  output$down_volcano = downloadHandler(
    filename = "volcano_plot.png",
    content = function(file){
      png(file, width = 1000, height = 1500)
      print(draw_volcano_plot())
      dev.off()
    }
  )

  #####################################################################

  ################      DIFF DOT PLOT + TABLE      ####################

  #####################################################################

  # differential dot plot table tutorial
  observeEvent(input$diff_dot_help, {
    shinyalert("Tutorial ", "Group data by choosing specific metadata in your dataset.
                Notice the new table with top differential genes you've selected on the previous plot.
                If dotplot is not readable, consider going on volcano panel and move f-value cutoff slider and press resize to get a smaller top differential table.
               ", type = "info", size = "m")
  })

  # warning message if user hasn't resize gene tab
  output$diff_dot_warning = renderText("Please go to volcano tab and press resize differential table")
  output$diff_dot_warning_findall = renderText("Sorry, this plot is not available for\"all cluster VS all cluster\" differential calculation option")

  # display and hide boxes
  observe({
    toggle(id = "diff_dot_box", condition = input$volcano_go & input$diff_calcul_type != "each cluster VS all clusters")
    toggle(id = "diff_dot_table_box", condition = input$volcano_go & input$diff_calcul_type != "each cluster VS all clusters")
    toggle(id = "diff_dot_inputs_box", condition = input$volcano_go & input$diff_calcul_type != "each cluster VS all clusters")
    toggle(id = "diff_dot_warning_box", condition = !input$volcano_go & input$diff_calcul_type != "each cluster VS all clusters")
    toggle(id = "diff_dot_warning_box_findall", condition = input$diff_calcul_type == "each cluster VS all clusters" )
  })

  # display reactive values
  observe({
    validate(need(!is.null(data_2()),"data 2 null"))
    updateSelectInput(session, inputId = "diff_dot_plot_group_by", choice = meta_data(), selected = FALSE)
  })

  # render table
  output$diff_dot_table = renderTable({
    validate(need(!is.null(top_n_genes()),"top genes null"))
    rownames(head(top_n_genes()))
    }, rownames =  TRUE)

  # download table
  output$down_diff_dot_table = downloadHandler(
    filename = "top_n_differential_table.csv",
    content = function(file){
      write.csv(top_n_genes(),file)
    }
  )

  # draw plot
  draw_diff_dot_plot = eventReactive(input$diff_dot_plot_go,{
    dittoDotPlot(data(), vars = rownames(top_n_genes()) ,group.by = input$diff_dot_plot_group_by )
  })

  # display plot
  output$diff_dot_plot_ui = renderUI({
    plotOutput("diff_dot_plot", width = paste0(input$diff_dot_plot_width, "%"), height = input$diff_dot_plot_height)
  })
  output$diff_dot_plot = renderPlot({
    draw_diff_dot_plot()
  })

  # download plot
  output$down_diff_dot_plot = downloadHandler(
    filename = "diff_dot_plot.png",
    content = function(file){
      png(file,width = 1000,height = 1000)
      print(draw_diff_dot_plot())
      dev.off()
    }
  )

  #####################################################################

  ######################      HEATMAP      ############################

  #####################################################################

  # heatmap tutorial
  observeEvent(input$heatmap_help,{
    shinyalert("Tutorial",
               "Press \"draw heatmap\" button to get the plot.
               If heatmap is not readable, consider going on volcano panel and move f-value cutoff slider and press resize to get a smaller top differential table.",
               type = "info", size = "m")
  })

  output$heatmap_warning = renderText("Please go to volcano tab and press resize differential table")

  # display and hide boxes
  observe({
    toggle(id = "heatmap_warning_box", condition = !input$volcano_go)
    toggle(id = "heatmap_box", condition = input$volcano_go)
  })

  # display plot

  draw_heatmap = eventReactive(input$heatmap_go,{
    df = data_2()
    DefaultAssay(df) = "integrated"
    DoHeatmap(df, features = rownames(top_n_genes())) + NoLegend()
  })

  output$Heatmap_plot_ui = renderUI({
    plotOutput("Heatmap_plot", width = paste0(input$heatmap_width, "%"), height = input$heatmap_height)
  })

  output$Heatmap_plot = renderPlot({
    draw_heatmap()
  })

  # download plot
  output$down_heatmap = downloadHandler(
    filename = "Heatmap.png",
    content = function(file){
      png(file, width = 1000, height = 1000)
      print(DoHeatmap(data(), features = rownames(top_n_genes())) + NoLegend())
      dev.off()
    }
  )


  #####################################################################

  #####################     CLUSTER PROFILER    #######################

  #####################################################################

  # CP tutorial
  observeEvent(input$CP_help, {
    shinyalert("Tutorial ", "This tab based on cluster profiler package allows to estimate biological pathways, connected to all differential genes you've calculated before.
               It's based on estimations and predictions, and results should be taken carefully.
               Select species you're studying.
               Select which genes to analyze, down regulated or upregulated", type = "info", size = "m")
  })

  # display table
  output$cluster_table_output = renderTable(top_n_genes(), rownames = TRUE)

  # get genes choosed
  CP_genes = eventReactive(input$CP_go,{
    DE = diff_genes()
    DE_up = as.data.frame(DE)
    DE_up = DE_up[0,]
    DE_down = as.data.frame(DE)
    DE_down = DE_down[0,]
    for (i in 1:length(DE$avg_log2FC)) {
      if(DE$avg_log2FC[i] > 0){
        DE_up = rbind(DE_up, DE[i,])
      }else{
        DE_down = rbind(DE_down, DE[i,])
      }
    }
    gc()
    if(input$CP_up_down == "all genes"){
      return(DE)
    }else if (input$CP_up_down == "up regulated"){
      return(DE_up)
    }else if(input$CP_up_down == "down regulated"){
      return(DE_down)
    }
  })

  # get species choosed
  CP_specie_choice = eventReactive(input$CP_species,{
    switch(input$CP_species,
           "human" = org.Hs.eg.db,
           "mouse"  = org.Mm.eg.db,
           "rat" = org.Rn.eg.db
    )
  })

  # Gene onthology calculation
  CP_enrich_go = eventReactive(input$CP_go,{
    df = enrichGO(gene = rownames(CP_genes()),
                  universe = rownames(data()),
                  OrgDb         = CP_specie_choice(),
                  keyType       = 'SYMBOL',
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.5,
                  qvalueCutoff  = 0.5,
                  minGSSize=5)
    gc()
    return(df)
  })

  # draw plot
  draw_CP_barplot = eventReactive(
    c(input$CP_go,input$CP_barplot_categories),{
      barplot(CP_enrich_go(), showCategory = input$CP_barplot_categories,font.size =12, label_format = 50)
    })

  # display plot
  output$CP_barplot_ui = renderUI({
    plotOutput("CP_barplot", width = paste0(input$CP_barplot_width, "%"), height = input$CP_barplot_height)
  })
  output$CP_barplot = renderPlot({
    draw_CP_barplot()
  })

  # download plot
  output$down_CP_barplot = downloadHandler(
    filename = "Cluster_profiler.png",
    content = function(file){
      png(file, width = 2000, height = 1000)
      print(draw_CP_barplot())
      dev.off()
    }
  )

  # download table
  output$down_diff_CP_table = downloadHandler(
    filename = "top_n_differential_table.csv",
    content = function(file){
      write.csv(top_n_genes(),file)

    }
  )

  #####################################################################

  ###################     CLUSTER PROFILER GSEA #######################

  #####################################################################

  # CP tutorial
  observeEvent(input$CP_help_2, {
    shinyalert("Tutorial ", "This tab is smilar to the previous one, expect that the differential genes are ordered.
                             That means the first pathways (on top of the plot) are based on the more differential genes, the more representativ ones.
                             Dotplot could be divided into two part: suppressed or activated pathway, but you could obtain only supressed or only repressed genes.
                             After organism is choosed, the sliders for \"GSS\" indicates minimun and maximum size of genesets per pathways. That means you constrains the analysis between those two limits.
                             \"P-value cutoff\" is the discrimination power of your tests.
               ", type = "info", size = "m")
  })

  # get genes choosed
  CP_genes_2 = eventReactive(input$CP_go_2,{
    DE = diff_genes()
    DE_list = DE$avg_log2FC
    names(DE_list) = rownames(DE)
    DE_list = sort(DE_list, decreasing = TRUE)
    return(DE_list)
  })

  # get species choosed
  CP_specie_choice_2 = eventReactive(input$CP_species_2,{
    switch(input$CP_species_2,
           "human" = org.Hs.eg.db,
           "mouse"  = org.Mm.eg.db,
           "rat" = org.Rn.eg.db
    )
  })

  # Gene onthology calculation
  CP_GSEA = eventReactive(input$CP_go_2,{
    tryCatch({
      df  = gseGO(geneList = CP_genes_2(),
                  ont = "BP",
                  keyType = "SYMBOL",
                  minGSSize = input$GSS_min,
                  maxGSSize = input$GSS_max,
                  pvalueCutoff = input$GSEA_p_val,
                  verbose = TRUE,
                  OrgDb = CP_specie_choice_2(),
                  pAdjustMethod = "bonferroni",
                  eps = 0)
      gc()
      return(df)
    }, error = function(e) {
      shinyalert("Warning", "some error advice in calc", type ="error")
      return(NULL)
    })
  })

  # draw plot
  draw_CP_dotplot = eventReactive(
    c(!is.null(CP_GSEA()),input$CP_dotplot_categories),{
      tryCatch({
        gene_count<- CP_GSEA()@result %>% group_by(ID) %>% summarise(count = sum(str_count(core_enrichment, "/")) + 1)
        dot_df<- left_join(CP_GSEA()@result, gene_count, by = "ID") %>% mutate(GeneRatio = count/setSize)
        dot_df = dot_df[1:input$CP_dotplot_categories,]
        dot_df$type = "upregulated"
        dot_df$type[dot_df$NES < 0] = "downregulated"
        p <- ggplot(dot_df, aes(x = GeneRatio, y = fct_reorder(Description, GeneRatio)), font.size =12,label_format = 5) +
          geom_point(aes(size = GeneRatio, color = p.adjust)) +
          scale_colour_gradient(limits=c(0, 0.10), low="red") +
          ylab(NULL) +
          ggtitle("GO pathway enrichment")+
          facet_grid(.~type)
        return(p)
      }, warning = function(w) {
        shinyalert("Warning", "some warning advice in plot", type ="warning")
        return(NULL)
      }, error = function(e) {
        if(is.null(CP_GSEA())){
          return(NULL)
        }else{
          shinyalert("Warning", "some error advice in plot", type ="error")
          return(NULL)
        }
      })

    })

  # display plot
  output$CP_dotplot_ui = renderUI({
    validate(
      need(!is.null(draw_CP_dotplot()), "error : plot problem")
    )
    if(is.null(draw_CP_dotplot()))
    {
      return(NULL)
    }
    else
    {
      plotOutput("CP_dotplot", width = paste0(input$CP_dotplot_width, "%"), height = input$CP_dotplot_height)
    }
  })

  output$CP_dotplot = renderPlot({
    validate(
      need(input$CP_go_2, "")
    )
    if(is.null(draw_CP_dotplot()))
    {
       return(NULL)
    }
    else
    {
      draw_CP_dotplot()
    }
  })

  # download plot
  output$down_CP_dotplot = downloadHandler(
    filename = "GSEA.png",
    content = function(file){
      png(file, width = 2000, height = 1000)
      print(draw_CP_dotplot())
      dev.off()
    }
  )

  #####################################################################

  ##################   DETAILED GSEA PLOT     #########################

  #####################################################################

  # CP tutorial
  observeEvent(input$CP_help_3, {
    shinyalert("Tutorial ", "If you have calculated some pathways in the previous panel, here you can plot them to see precisely how many genes are involved in the pathway and how much the pathway is suppressed or activated.
               Select a pathways in the box and press run.
               ", type = "info", size = "m")
  })

  observe({
    if(!is.null(CP_GSEA())){
      updateSelectInput(session,inputId = "CP_detailed_pathway", choices = CP_GSEA()@result$Description )
    }
  })

  observeEvent(input$tabs,{
    if(input$tabs == "cluster_profiler_3" & is.null(draw_CP_dotplot())) {
      shinyalert("Warning", "Go to GSEA and run calculation", type ="warning")
    }
  })

  observeEvent(input$tabs,{
    if(input$tabs == "cluster_profiler_3" & !input$CP_go_2) {
      shinyalert("Warning", "Go to GSEA and run calculation", type ="warning")
    }
  })

  n_pathway_item = eventReactive(input$CP_detailed_pathway,{
    n = which(CP_GSEA()@result$Description == input$CP_detailed_pathway)
    return(n)
  })

  # draw plot
  draw_CP_detailed = eventReactive(input$CP_go_3,{
    gseaplot(CP_GSEA(),n_pathway_item(),title = CP_GSEA()@result$Description[n_pathway_item()])
  })

  # display plot
  output$CP_detailed_ui = renderUI({
    plotOutput("CP_detailed", width = paste0(input$CP_detailed_width, "%"), height = input$CP_detailed_height)
  })
  output$CP_detailed = renderPlot({
    draw_CP_detailed()
  })

  # download plot
  output$down_CP_detailed = downloadHandler(
    filename = "GSEA_detailed.png",
    content = function(file){
      png(file, width = 2000, height = 1000)
      print(draw_CP_detailed())
      dev.off()
    }
  )

  # download table
  output$down_diff_CP_table_2 = downloadHandler(
    filename = "top_n_differential_table.csv",
    content = function(file){
      write.csv(top_n_genes(),file)
    }
  )


  #####################################################################

  ##################      NICHENET     ################################

  #####################################################################

  # nichenet tutorial
  observeEvent(input$niche_help, {
    shinyalert("Tutorial ", "Nichenet allows to estimate ligand/target relationship in your data set.
               It's based on estimation, calculation and known relations. That means results should be taken carefully.
               First, choose which group should be consider at potential senders and receiver (\"sender group(s)\" and \"receiver group(s)\")
               Then you have to indicate two conditions (\"receiver condition\" and \"reference condition\") , to perform a differential analysis between those.
               Differential genes on those conditions will be analyse to find potential ligands.
               The table given as an output show genes names, statistical values and the last column \"bona fide ligand\" indicate if the relation is known (\"True\") or predicted (\"False\")", type = "info", size = "m")
  })

  # display reactive values
  observe({
    validate(need(!is.null(data()),"data null"))
    updateSelectInput(session, inputId = "niche_condition", choice = meta_data(), selected = FALSE)
  })

  observe({
    validate(need(!is.null(data()),"data null"))
    updateSelectInput(session, inputId = "niche_idents", choice = meta_data(), selected = FALSE)
  })


  observe({
    validate(need(!is.null(data()),"data null"))
    updateSelectInput(session, inputId = "niche_sender", choices = active_idents())
    updateSelectInput(session, inputId = "niche_receiver", choices = active_idents())
  })

  # hide and display boxes
  observe({
    toggle(id = "niche_1", condition = input$niche_go_2)
    toggle(id = "niche_2", condition = input$niche_go_2)
    toggle(id = "niche_3", condition = input$niche_go_2)
    toggle(id = "niche_warning_1", condition = !input$niche_go_2)
    toggle(id = "niche_warning_2", condition = !input$niche_go_2)
    toggle(id = "niche_warning_3", condition = !input$niche_go_2)
  })

  # warning messages
  output$niche_message_1 = renderText("Please go to Nichenet calculation and press run")
  output$niche_message_2 = renderText("Please go to Nichenet calculation and press run")
  output$niche_message_3 = renderText("Please go to Nichenet calculation and press run")

  # upload reference files
  ligand_target_matrix = eventReactive(input$niche_go_2,{
    f = readRDS("./ressources/ligand_target_matrix.rds")
    return(f)
  })

  lr_network = eventReactive(input$niche_go_2,{
    f = readRDS("./ressources/lr_network.rds")
    return(f)
  })

  weighted_networks = eventReactive(input$niche_go_2,{
    f = readRDS("./ressources/weighted_networks.rds")
    return(f)
  })

  # display meta
  niche_choices = eventReactive(input$niche_condition,{
    validate(need(!is.null(data()),"data null"))
    expr = sort(unique(eval(parse(text = paste0("data()$",input$niche_condition)))))
  })

  observe({
    updateSelectInput(session, inputId = "niche_condition_oi", choices = niche_choices() )
    updateSelectInput(session, inputId = "niche_condition_ref", choices = niche_choices() )
  })

  # set receiver cluster
  niche_clusters_receiver = eventReactive(input$niche_receiver,{
    if(input$niche_idents == "seurat_clusters")
    {
      my_list = as.numeric(unlist(strsplit(input$niche_receiver,",")))
    }
    else
    {
      my_list = input$niche_receiver
    }
    return(my_list)
  })

  # set sender cluster
  niche_clusters_sender = eventReactive(input$niche_sender,{
    if(input$niche_idents == "seurat_clusters")
    {
      my_list = as.numeric(unlist(strsplit(input$niche_sender,",")))
    }
    else
    {
      my_list = input$niche_sender
    }
    return(my_list)
  })

  # nichenet calculation
  nichenet_output = eventReactive(input$niche_go_2,{
    nichenet_seuratobj_aggregate(
      seurat_obj = data_2(),
      receiver = niche_clusters_receiver(),
      sender =  niche_clusters_sender(),
      condition_colname = input$niche_condition, condition_oi = input$niche_condition_oi, condition_reference = input$niche_condition_ref,
      ligand_target_matrix = ligand_target_matrix(), lr_network = lr_network(), weighted_networks = weighted_networks(), organism = "human")
  })

  # display table
  output$niche_table = renderTable(nichenet_output()$ligand_activities, rownames = TRUE)

  # download table
  output$down_niche_table = downloadHandler(
    filename = "niche_table.csv",
    content = function(file){
      write.csv(nichenet_output()$ligand_activities,file)
    }
  )

  ##################      DOTPLOT      ################################

  # nichenet tutorial
  observeEvent(input$niche_dot_help, {
    shinyalert("Tutorial ", "This plot show expression levels of the ligands calculated among youre clusters", type = "info", size = "m")
  })

  # draw plot
  draw_niche_dot = eventReactive(input$niche_go_2,{
    DotPlot(data(), features = nichenet_output()$top_ligands %>% rev(), cols = "RdYlBu") + RotatedAxis()
  })

  # display plot
  output$niche_dot_ui = renderUI({
    plotOutput("niche_dot", width = paste0(input$niche_dot_width, "%"), height = input$niche_dot_height)
  })
  output$niche_dot = renderPlot({
    draw_niche_dot()
  })

  # download plot
  output$down_niche_dot = downloadHandler(
    filename = "niche_dotplot.png",
    content = function(file){
      png(file, width = 1000, height = 1000)
      print(draw_niche_dot())
      dev.off()
    }
  )

  ##################      HEATMAP      ################################

  # nichenet tutorial
  observeEvent(input$niche_map_help, {
    shinyalert("Tutorial ", "This plot show ligand activity on target genes", type = "info", size = "m")
  })

  # draw plot
  draw_niche_map = eventReactive(input$niche_go_2,{
    nichenet_output()$ligand_target_heatmap
  })

  # display plot
  output$niche_map_ui = renderUI({
    plotOutput("niche_map", width = paste0(input$niche_map_width, "%"), height = input$niche_map_height)
  })
  output$niche_map = renderPlot({
    draw_niche_map()
  })

  # download plot
  output$down_niche_map = downloadHandler(
    filename = "niche_heatmap.png",
    content = function(file){
      png(file, width = 1000, height = 1000)
      print(draw_niche_map())
      dev.off()
    }
  )

  ##################      TARGET       ################################

  # nichenet tutorial
  observeEvent(input$niche_target_help, {
    shinyalert("Tutorial ", "This plot show ligand activity on target genes, ligand expression among clusters and differential potential.", type = "info", size = "m")
  })

  # draw plot
  draw_niche_target = eventReactive(input$niche_go_2,{
    nichenet_output()$ligand_activity_target_heatmap
  })

  # display plot
  output$niche_target_ui = renderUI({
    plotOutput("niche_target", width = paste0(input$niche_target_width, "%"), height = input$niche_target_height)
  })
  output$niche_target = renderPlot({
    draw_niche_target()
  })

  # download plot
  output$down_niche_target = downloadHandler(
    filename = "niche_target_plot.png",
    content = function(file){
      png(file, width = 1000, height = 1000)
      print(draw_niche_target())
      dev.off()
    }
  )

  ##############################                    END                    ######################################################

  session$onSessionEnded(function() {
    print("Ending session")
    rm(list=ls())
    gc()
    print("Memory cleaned")
  })

}
