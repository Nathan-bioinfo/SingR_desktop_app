source("func.R")
load_packages()

ui = dashboardPage(
  # title
  dashboardHeader(title = "singR"),
  # lateral menu
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem("Welcome", tabName = "Intro"),
      menuItem("Upload data", tabName = "Uploading"),
      menuItem("Quality control", tabName = "quality",
               menuSubItem("Violin plot", tabName = "boxplot"),
               menuSubItem("Barplot", tabName = "barplot")
      ),
      menuItem("Gene expression", tabName = "expression",
               menuSubItem("Student test", tabName = "gene_expression"),
               menuSubItem("ANOVA", tabName = "gene_expression_2"),
               menuSubItem("Ridge plot", tabName = "ridge"),
               menuSubItem("Dot plot", tabName = "dot_plot")
      ),
      menuItem("Gene sets", tabName = "geneset"),
      menuItem("Reductions", tabName = "Clusters"),
      menuItem("Density per marker", tabName ="Density_gen",
               menuSubItem("1 to 4 markers", tabName = "Density"),
               menuSubItem("2 markers (color threshold)", tabName = "Density_2")
      ),
      menuItem("Differential genes", tabName = "Differential",
               menuItem("Find differential genes", tabName = "differential_calculation"),
               menuItem("Volcano plot", tabName = "Volcano"), 
               menuItem("Differential dotplot", tabName = "diff_dot"),
               menuItem("Heatmap", tabName = "Heatmap"),
               menuItem("Gene ontology", tabName = "GO",
                        menuSubItem("enrichment GO", tabName = "cluster_profiler"),
                        menuSubItem("GSEA global pathways", tabName = "cluster_profiler_2"),
                        menuSubItem("GSEA detailed pathway", tabName = "cluster_profiler_3")
                        )
      ),
      menuItem("Nichenet", tabName = "Nichenet",
               menuSubItem("Calculation", tabName = "niche_calc"),
               menuSubItem("Dot Plot", tabName = "niche_dot"),
               menuSubItem("Heatmap", tabName = "niche_map"),
               menuSubItem("Target Heatmap", tabName = "niche_target")
      ) 
    )
  ),
  # corpus
  dashboardBody(
    fluidRow(
      tags$head(
        tags$style(
          HTML(".shiny-output-error-validation
          {
            color: #ff0000;
            font-weight: bold;
          }"
          )
        )
      ),
      useShinyjs(),
      add_busy_spinner(spin = "fading-circle" , position = "bottom-left", color = "white", height = "100px", width = "100px"),
      tabItems(
         tabItem("Intro",
                shinydashboardPlus::box(
                  title = "Hello",
                  status = "primary",
                  verbatimTextOutput("introduction"),
                  img(src = "UN.png", height = 150, width = 330),
                  img(src = "CRTI.png", height = 150, width = 400),
                  img(src = "INSERM.png", height = 150, width = 330),
                  width = 9
                )
        ), 
        tabItem("Uploading",
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Uploading area",
                  width = 2,
                  actionButton("intro_help", "Help"),
                  radioButtons("file_type", "Upload dataframe or use example ?", c("example dataframe" = "ex","upload personnal dataframe" = "dataf"), selected = "ex"),
                  shinyFilesButton("user_file_path", "Choose a file" ,title = "Please select a file:", multiple = FALSE, buttonType = "default", class = NULL),
                  disabled(actionButton("preloaded_go", "load")),
                  selectInput("trim","metadata to keep", choices = NULL, multiple = T),
                  actionButton("trim_go", "trim")
                )
        ),
        tabItem("boxplot",
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Boxplot Inputs",
                  width = 2,
                  actionButton("boxplot_help", "help"),
                  selectInput(inputId = "boxplot_option", label = "display", choices = c("RNA count", "RNA feature", "% mitochondrial")),
                  selectInput(inputId = "boxplot_group_by", label = "group by" , choices = c()),
                  checkboxInput("boxplot_idents_activation", label = "isolate specific classes", value = FALSE),
                  selectInput(inputId = "boxplot_idents", label = "classe(s) to keep", choice = c(), multiple = TRUE, selectize = TRUE),
                  radioButtons("boxplot_custom_colors","colors",choices = c("default","custom","white","black")),
                  selectizeInput(inputId = "boxplot_colors", label = "colors", choices = colors(), multiple = TRUE),
                  actionButton(inputId = "boxplot_go", label = "run")
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Boxplot",
                  id = "boxplot_box",
                  uiOutput("boxplot_ui"),
                  fluidRow(
                    column(3,sliderInput(inputId = "boxplot_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                  ),
                  fluidRow(
                    column(3,sliderInput(inputId = "boxplot_height", label = "Change height (px)", min = 0, max = 800, value = 400))
                  ),
                  fluidRow(
                    column(1,checkboxInput(inputId = "boxplot_hide_dots", label = "Hide dots", value = FALSE)),
                    column(1,checkboxInput(inputId = "boxplot_add_box", label = "Show boxplot", value = FALSE)),
                    column(4,downloadButton("down_boxplot", label = "Download"))
                  ),
                  width = 10
                )
        ),
        tabItem("barplot",
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Barplot inputs",
                  width = 2,
                  actionButton("barplot_help", "help"),
                  selectInput(inputId = "barplot_group_by", label = "group by" , choices = c()),
                  selectInput(inputId = "barplot_meta", label = "metadata", choices = c()),
                  selectInput(inputId = "barplot_scale", label = "scale", choices = c("percent", "count")),
                  radioButtons("barplot_custom_colors", "colors", choices = c("default", "custom")),
                  selectizeInput(inputId = "barplot_colors", label = "colors", choices = colors(), multiple = TRUE),
                  actionButton(inputId = "barplot_go", label = "run")
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Barplot",
                  id = "barplot_box",
                  uiOutput("barplot_ui"),
                  fluidRow(
                    column(3,sliderInput(inputId = "barplot_width", label = "Change width (%)", min = 0, max = 100, value = 60))
                  ),
                  fluidRow(
                    column(3,sliderInput(inputId = "barplot_height", label = "Change height (px)", min = 0, max = 800, value = 400))
                  ),
                  downloadButton("down_barplot", label = "Download"),
                  width = 10
                )
        ),
        tabItem("gene_expression",
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Gene expression inputs",
                  actionButton("gene_expression_help", "help"),
                  selectizeInput(inputId = "expression_gene_choice",label = 'choose one gene (auto-completion)',choices = NULL,selected = NULL, multiple = TRUE, options = list(maxItems = 1)),
                  selectInput(inputId = "expression_group_by_choice", label = "group by", choices = c()),
                  checkboxInput("student_idents_activation", label = "isolate specific classes", value = FALSE),
                  selectInput(inputId = "student_idents", label = "classe(s) to keep", choice = c(), multiple = TRUE, selectize = TRUE),
                  radioButtons("expression_custom_colors","colors",choices = c("default","custom","white","black")),
                  selectizeInput(inputId = "expression_colors", label = "colors", choices = colors(), multiple = TRUE),
                  checkboxInput("student_stats_activation", "stats"),
                  selectizeInput("student_stats", "statistics", multiple = TRUE, choices = NULL),
                  actionButton(inputId = "gene_expression_go", label = "run"),
                  width = 2
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Gene expression plot",
                  uiOutput("gene_expression_ui"),
                  fluidRow(
                    column(3,sliderInput(inputId = "gene_expression_y_lim", label = "change y axis limits", min = 1, max = 20, value = 10))
                  ),
                  fluidRow(
                   column(3,sliderInput(inputId = "gene_expression_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                   column(3,sliderInput(inputId = "gene_expression_height", label = "Change height (px)", min = 0, max = 800, value = 400)),
                  ),
                  fluidRow(
                    column(1,checkboxInput(inputId = "gene_expression_hide_dots", label = "Hide dots", value = FALSE)),
                    column(1,checkboxInput(inputId = "expression_add_box", label = "Show boxplot", value = FALSE)),
                    column(1,checkboxInput(inputId = "p_format", label = "display p value", value = FALSE)),
                    column(3,downloadButton("down_gene_expression", label = "Download"))
                  ),
                  width = 10
                )
        ),
        tabItem("gene_expression_2",
                shinydashboardPlus::box(
                  status = "primary",
                  title= "Anova inputs",
                  actionButton("gene_expression_help_2", "help"),
                  selectizeInput(inputId = "expression_gene_choice_2",label = 'choose one gene (auto-completion)',choices = NULL,selected = NULL, multiple = TRUE, options = list(maxItems = 1)),
                  selectInput(inputId = "expression_group_by_choice_2", label = "group by", choices = c()),
                  checkboxInput(inputId = "split_by_anova_activation", label = "split by", value = FALSE),
                  selectInput(inputId = "expression_split_by_choice_2", label = "split by", choices = c()),
                  checkboxInput(inputId = "anova_idents_activation", label = "isolate specific classes", value = FALSE),
                  selectInput(inputId = "anova_idents", label = "classe(s) to keep", choice = c(),selected = NULL, multiple = TRUE, selectize = TRUE),
                  radioButtons("expression_custom_colors_2","colors",choices = c("default","custom","white","black")),
                  selectizeInput(inputId = "expression_colors_2", label = "colors", choices = colors(), multiple = TRUE),
                  actionButton(inputId = "gene_expression_go_2", label = "run"),
                  width = 2
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Anova plot",
                  uiOutput("gene_expression_ui_2"),
                  fluidRow(
                    column(3,sliderInput(inputId = "gene_expression_width_2", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "gene_expression_height_2", label = "Change height (px)", min = 0, max = 800, value = 400))
                  ),
                  fluidRow(
                     column(3,sliderInput(inputId = "gene_expression_y_lim_2", label = "change y axis limits", min = 1, max = 20, value = 7)),
                     column(1, style = "margin-top:50px;",checkboxInput(inputId = "gene_expression_hide_dots_2", label = "Hide dots", value = FALSE)),
                     column(1, offset = 0, style = "margin-top:50px;",checkboxInput(inputId = "expression2_add_box", label = "Boxplot", value = FALSE)),
                     column(1, style = "margin-top:50px;",downloadButton("down_gene_expression_2", label = "Download"))
                   ),
                  width = 10
                )
        ),
        tabItem("ridge",
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Ridge inputs",
                  actionButton("ridge_help","help"),
                  selectizeInput(inputId = "ridge_gene",label = 'choose one gene (auto-completion)',choices = NULL, multiple = FALSE, selected = NULL),
                  selectInput(inputId = "ridge_group_by", label = "group by", choices = c()),
                  checkboxInput(inputId = "ridge_idents_activation", label = "isolate specific classes", value = FALSE),
                  selectInput(inputId = "ridge_idents", label = "classe(s) to keep", choice = c(),selected = NULL, multiple = TRUE, selectize = TRUE),
                  radioButtons("ridge_custom_colors","colors",choices = c("default","custom","white","black")),
                  selectizeInput(inputId = "ridge_colors", label = "colors", choices = colors(), multiple = TRUE),
                  actionButton(inputId = "ridge_go", label = "run"),
                  width = 2
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Ridgeplot",
                  uiOutput("ridge_plot_ui"),
                  fluidRow(
                    column(3,sliderInput(inputId = "ridge_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "ridge_height", label = "Change height (px)", min = 0, max = 800, value = 400))
                  ),
                  downloadButton("down_ridge", label = "Download"),
                  width = 10
                )
        ),
        tabItem("dot_plot",
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Dotplot inputs",
                  actionButton("dotplot_help", "help"),
                  selectizeInput(inputId = 'dot_plot_genes',label = 'Genes (as many as you want)',choices = NULL,selected = NULL, multiple = TRUE,options = list(create = FALSE)),
                  selectInput(inputId = "dot_plot_group_by", label = "group by", choices = c()),
                  checkboxInput("dot_plot_idents_activation", label = "isolate specific classes", value = FALSE ),
                  selectInput(inputId = "dot_plot_idents", label = "idents", choices = c(), selectize = TRUE, multiple = TRUE),
                  actionButton(inputId = "dot_plot_go", label = "run"),
                  width = 2
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Dotplot",
                  uiOutput("dot_plot_ui"),
                  fluidRow(
                    column(3,sliderInput(inputId = "dot_plot_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "dot_plot_height", label = "Change height (px)", min = 0, max = 800, value = 400))
                  ),
                  downloadButton("down_dot_plot", label = "Download"),
                  width = 10
                )
        ),
        tabItem("geneset",
                column(width = 2,
                       shinydashboardPlus::box(
                         status = "primary",
                         title = "Geneset inputs",
                         actionButton("gs_help", "help"),
                         radioButtons("gs_file_type", "Geneset dataframe", c("example","upload")),
                         radioButtons("gs_file_extension", "File type", c("csv","xls","ods","txt","rds")),
                         shinyFilesButton("gs_user_file_path", "Choose a file" ,title = "Please select a file:", multiple = FALSE, buttonType = "default", class = NULL),
                         radioButtons("gs_file_sep", "Separator", c(",",";","tab")),
                         actionButton("gs_go", "create gene list"),
                         selectInput("gs_idents", "idents", choices = c()),
                         selectInput("gs_cond_1", "condition 1", choices = c()),
                         selectInput("gs_cond_2", "condition 2", choices = c()),
                         selectInput("gs_sub_cat", "sub category", choices =c()),
                         selectInput("gs_sub_cat_choices", "keep", choices = c(),multiple = TRUE),
                         actionButton("gs_go_2", "run"),
                         width = NULL
                       ),
                       shinydashboardPlus::box(
                         status = "primary",
                         title = "Geneset colors",
                         radioButtons("geneset_custom_colors","colors",choices = c("default","custom","white","black")),
                         selectizeInput(inputId = "geneset_colors", label = "colors", choices = colors(), multiple = TRUE),
                         actionButton("gs_go_3", "run"),
                         width = NULL
                       )
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  id = "gs_table_box",
                  title = "Table quick view (first 5 columns), choose best separator on left panel",
                  tableOutput("gs_table"),
                  width = 10
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "GS plot",
                  uiOutput("gs_plot_ui"),
                  fluidRow
                  (
                    column(3,selectInput("gs_plot_meta", "condition", choices = NULL,width = "400px")),
                    column(3,selectizeInput("gs_stats", "statistics", multiple = TRUE, choices = NULL, width = "400px"))
                  ),
                  fluidRow
                  (
                    column(3,sliderInput(inputId = "gs_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "gs_height", label = "Change height (px)", min = 0, max = 800, value = 400))
                  ),
                  fluidRow(
                    column(1,checkboxInput("gs_plot_legend", "add legend")),
                    column(1,checkboxInput("gs_plot_add_boxes", "add boxes")),
                    column(1,checkboxInput("gs_plot_angle", "labels angle")),
                    column(1,checkboxInput("gs_align_bar", "align stats bar")),
                    column(2,actionButton("gs_download", "download"))
                  ),
                  width = 10
                )
        ),
        tabItem("Clusters",
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Reductions inputs",
                  actionButton("cluster_help", "help"),
                  selectInput(inputId = "cluster_reduction" , label = "reduction", choices = NULL),
                  selectInput(inputId = "cluster_group_by", label = "group by", choices = c()),
                  selectInput(inputId = "cluster_labels", label = "add label", choices = c("yes","no")),
                  checkboxInput("reduction_activation", "calculate"),
                  textInput("pca_reductions","pca reductions", value =50),
                  textInput("umap_reductions","umap reductions", value = 20),
                  textInput("tsne_reductions","tsne reductions", value = 5),
                  actionButton("reduction_go", "calculate reduction"),
                  actionButton(inputId = "cluster_go", label = "plot"),
                  width = 2
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Clustering",
                  uiOutput("Clusters_plot_ui"), 
                  fluidRow(
                    column(3,sliderInput(inputId = "cluster_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "cluster_height", label = "Change height (px)", min = 0, max = 800, value = 400))
                  ),
                  downloadButton("down_clusters", label = "Download"),
                  width = 10
                )
        ),
        tabItem("Density",
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Density inputs",
                  actionButton("density_help", "help"),
                  selectizeInput(inputId = "density_gene_1",label = 'choose first gene (auto-completion)',choices = NULL,selected = NULL, multiple = TRUE, options = list(maxItems = 1)),
                  checkboxInput(inputId = "density_gene_2_activation", label = "add gene"),
                  selectizeInput(inputId = "density_gene_2",label = 'choose second gene',choices = NULL,selected = NULL, multiple = TRUE, options = list(maxItems = 1)),
                  checkboxInput(inputId = "density_gene_3_activation", label = "add gene"),
                  selectizeInput(inputId = "density_gene_3",label = 'choose third gene',choices = NULL,selected = NULL, multiple = TRUE, options = list(maxItems = 1)),
                  checkboxInput(inputId = "density_gene_4_activation", label = "add gene"),
                  selectizeInput(inputId = "density_gene_4",label = 'choose fourth gene',choices = NULL,selected = NULL, multiple = TRUE, options = list(maxItems = 1)),
                  actionButton(inputId = "density_go", label = "run"),
                  width = 2
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Density plot",
                  uiOutput("Density_plot_ui"), 
                  fluidRow(
                    column(3,sliderInput(inputId = "density_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "density_height", label = "Change height (px)", min = 0, max = 800, value = 400))
                  ),
                  downloadButton("down_density", label = "Download"),
                  width = 10
                )
        ),
        tabItem("Density_2",
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Density inputs",
                  actionButton("density_help_2", "help"),
                  selectizeInput(inputId = "density_gene_1_2",label = 'choose first gene (auto-completion)',choices = c(),selected = NULL, multiple = TRUE, options = list(maxItems = 1)),
                  selectizeInput(inputId = "density_gene_2_2",label = 'choose second gene',choices = c(),selected = NULL, multiple = TRUE, options = list(maxItems = 1)),
                  actionButton(inputId = "density_go_2", label = "run"),
                  width = 2
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Density plot",
                  uiOutput("Density_plot_ui_2"), 
                  fluidRow(
                    column(3,sliderInput(inputId = "density_width_2", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "density_height_2", label = "Change height (px)", min = 0, max = 800, value = 400))
                  ),
                  downloadButton("down_density_2", label = "Download"),
                  width = 10
                )
        ),
        tabItem("differential_calculation",
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Differential inputs",
                  actionButton("diff_help", "help"),
                  selectInput("diff_idents", "idents", choices = NULL),
                  selectInput(inputId = "diff_calcul_type", label = "choose a method", choices = c( "1 VS all","multiple VS multiple","each VS all")),
                  selectInput(inputId = "diff_calcul_cluster_1", label = "group 1",choices = c(),selected = NULL, multiple = FALSE, selectize = TRUE),
                  selectInput(inputId = "diff_calcul_cluster_1_2", label = "group 1",choices = c(),selected = NULL, multiple = TRUE, selectize = TRUE),
                  selectInput(inputId = "diff_calcul_cluster_2", label = "group 2",choices = c(),selected = NULL, multiple = TRUE, selectize = TRUE),
                  selectInput(inputId = "diff_meta", label = "metadata to look up", choices = c()),
                  selectInput(inputId = "diff_keep", label = "keep", choices = c()),
                  sliderInput(inputId = "diff_calcul_fold_change", label = "log fold change", value = 0.8, min = 0.01, max = 1),
                  sliderInput(inputId = "diff_calcul_p_value", label = "p-value cutoff", min = 0.005, max = 1, value = 0.05),
                  actionButton(inputId = "diff_calcul_go", label = "run"),
                  width = 2
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Head differential general table",
                  downloadButton("down_diff_table", label = "Download full table"),
                  tableOutput(outputId = "diff_calcul_table"),
                  width = 10
                )
        ),
        tabItem("Volcano",
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Warning",
                  id = "warning_volcano_box",
                  width = 4,
                  verbatimTextOutput("warning_volcano")
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  id = "volcano_box",
                  title = "volcano plot",
                  uiOutput("volcano_plot_ui"),
                  fluidRow(
                    column(3,sliderInput(inputId = "volcano_f_value", label = "f-value cutoff", min = 0.001, max = 5, value = 0.2)),
                    column(3,sliderInput(inputId = "volcano_x_lim", label = "Change x axis limits", min = 0.1, max = 50, value = 2))
                  ),
                  fluidRow(
                    column(3,sliderInput(inputId = "volcano_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "volcano_height", label = "Change height (px)", min = 0, max = 1000, value = 600))
                  ),
                  fluidRow(
                    column(6,
                     actionButton("volcano_help", "help"),
                     actionButton("volcano_go", label = "resize differential table"),
                     downloadButton("down_volcano", label = "Download plot")
                    )
                  ),
                  width = 12
                )
        ),
        tabItem("diff_dot",
                column(width = 2,
                       shinydashboardPlus::box(
                         title = "Differential dotplot inputs",
                         id = "diff_dot_inputs_box",
                         width = NULL,
                         status = "primary",
                         actionButton("diff_dot_help", "help"),
                         selectInput(inputId = "diff_dot_plot_group_by", label = "1: group by", choices = c()),
                         actionButton(inputId = "diff_dot_plot_go", label = "run")
                       ),
                       shinydashboardPlus::box(
                         title = "Head top differential genes",
                         id = "diff_dot_table_box",
                         width = NULL,
                         status = "primary",
                         tableOutput("diff_dot_table"),
                         downloadButton("down_diff_dot_table", label = "Downlad")
                       )
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Differential dotplot",
                  id = "diff_dot_box",
                  width = 10,
                  uiOutput("diff_dot_plot_ui"),
                  fluidRow(
                    column(3,sliderInput(inputId = "diff_dot_plot_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "diff_dot_plot_height", label = "Change height (px)", min = 0, max = 1000, value = 600))
                  ),
                  downloadButton("down_diff_dot_plot", label = "Download")
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Warning",
                  id = "diff_dot_warning_box",
                  width = 4,
                  verbatimTextOutput("diff_dot_warning")
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Warning",
                  id = "diff_dot_warning_box_findall",
                  width = 3,
                  verbatimTextOutput("diff_dot_warning_findall")
                )
        ),
        tabItem("Heatmap",
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Warning",
                  id = "heatmap_warning_box",
                  width = 4,
                  verbatimTextOutput("heatmap_warning")
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Heatmap",
                  id = "heatmap_box",
                  uiOutput("Heatmap_plot_ui"),
                  fluidRow(
                    column(3,sliderInput(inputId = "heatmap_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "heatmap_height", label = "Change height (px)", min = 0, max = 1000, value = 600))
                  ),
                  fluidRow(
                    column(1,actionButton("heatmap_go", "Draw Heatmap")),
                    column(1,actionButton("heatmap_help", "help")),
                    column(1,downloadButton("down_heatmap", label = "Download"))
                  ),
                  width = 12
                )
        ),
        tabItem("cluster_profiler",
                column(width = 2,
                       shinydashboardPlus::box(
                         title = "Gene onthology inputs",
                         id = "cluster_profiler_inputs_box",
                         actionButton("CP_help", "help"),
                         selectInput(inputId = "CP_species", label = "organism", choices = c("human","mouse","rat")),
                         selectInput(inputId = "CP_up_down", label = "genes", choices = c("all genes","up regulated", "down regulated")),
                         actionButton(inputId = "CP_go", label = "run"),
                         width = NULL,
                         status = "primary"
                       )
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Gene onthology barplot",
                  id = "cluster_profiler_box",
                  uiOutput("CP_barplot_ui"),
                  fluidRow(
                    column(3,sliderInput(inputId = "CP_barplot_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "CP_barplot_height", label = "Change height (px)", min = 0, max = 1000, value = 600))
                  ),
                  fluidRow(
                    column(3,sliderInput(inputId = "CP_barplot_categories", label = "number of category to display", value = 25, min = 1, max = 100)),
                    column(1, style ="margin-left:200px; margin-top:45px;",downloadButton("down_CP_barplot", label = "Download"))
                  ),
                  width = 10
                )
        ),
        tabItem("cluster_profiler_2",
                column(width = 2,
                       shinydashboardPlus::box(
                         title = "Gene onthology inputs",
                         id = "cluster_profiler_inputs_box_2",
                         actionButton("CP_help_2", "help"),
                         selectInput(inputId = "CP_species_2", label = "organism", choices = c("human","mouse","rat")),
                         sliderInput("GSS_min", label = "GSS min", max = 400, min = 1, value = 40),
                         sliderInput("GSS_max", label = "GSS max", max = 2000, min = 50, value = 800),
                         sliderInput("GSEA_p_val", label = "p-value cutoff", min = 0.001, max = 1, value = 0.05),
                         actionButton(inputId = "CP_go_2", label = "run"),
                         width = NULL,
                         status = "primary"
                       )
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Gene onthology_dotplot",
                  id = "cluster_profiler_box_2",
                  uiOutput("CP_dotplot_ui"),
                  fluidRow(
                    column(3,sliderInput(inputId = "CP_dotplot_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "CP_dotplot_height", label = "Change height (px)", min = 0, max = 1000, value = 600))
                  ),
                  fluidRow(
                    column(3,sliderInput(inputId = "CP_dotplot_categories", label = "number of category to display", value = 25, min = 1, max = 100)),
                    column(1, style ="margin-left:200px; margin-top:45px;",downloadButton("down_CP_dotplot", label = "Download"))
                  ),
                  width = 10
                )
        ),
        tabItem("cluster_profiler_3",
                column(width = 2,
                       shinydashboardPlus::box(
                         title = "Detailed pathway inputs",
                         id = "cluster_profiler_inputs_box_3",
                         actionButton("CP_help_3", "help"),
                         selectInput(inputId = "CP_detailed_pathway", label = "pathway", choices = NULL),
                         actionButton(inputId = "CP_go_3", label = "run"),
                         width = NULL,
                         verbatimTextOutput("testi"),
                         status = "primary"
                       )
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "GO_detailed",
                  id = "cluster_profiler_box_3",
                  uiOutput("CP_detailed_ui"),
                  fluidRow(
                    column(3,sliderInput(inputId = "CP_detailed_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "CP_detailed_height", label = "Change height (px)", min = 0, max = 800, value = 400))
                  ),
                  downloadButton("down_CP_detailed", label = "Download"),
                  width = 10
                )
        ),
        tabItem("niche_calc",
                shinydashboardPlus::box(
                  title = "Ligands calculation inputs",
                  status = "primary",
                  actionButton("niche_help", "help"),
                  selectInput("niche_idents", "idents", choices = c()),
                  selectInput(inputId = "niche_sender", label = "sender group(s)", choices = c(), selected = FALSE, multiple = TRUE, selectize = TRUE),
                  selectInput(inputId = "niche_receiver", label = "receiver group(s)", choices = c(), selected = FALSE, multiple = TRUE, selectize = TRUE),
                  selectInput(inputId = "niche_condition", label = "conditions to compare", choices = c()),
                  selectInput(inputId = "niche_condition_oi", label = "receiver condition", choices = c()),
                  selectInput(inputId = "niche_condition_ref", label = "reference condition", choices = c()),
                  actionButton("niche_go_2", label = "run"),
                  width = 2
                ),
                shinydashboardPlus::box(
                  status = "primary",
                  title = "Ligands table",
                  downloadButton("down_niche_table", label = "Download table"),
                  tableOutput(outputId ="niche_table"),
                  width = 10
                )
        ),
        tabItem("niche_dot",
                shinydashboardPlus::box(
                  id = "niche_1",
                  title = "Ligand dotplot",
                  status = "primary",
                  uiOutput("niche_dot_ui"),
                  fluidRow(
                    column(3,sliderInput(inputId = "niche_dot_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "niche_dot_height", label = "Change height (px)", min = 0, max = 800, value = 400))
                  ),
                  downloadButton("down_niche_dot", label = "Download"),
                  actionButton("niche_dot_help", label = "Help"),
                  width = 12 
                ),
                shinydashboardPlus::box(
                  id = "niche_warning_1",
                  title = "Warning",
                  status = "primary",
                  width = 4,
                  verbatimTextOutput("niche_message_1")
                )
        ),
        tabItem("niche_map",
                shinydashboardPlus::box(
                  status = "primary",
                  id = "niche_2",
                  title = "Ligand heatmap",
                  uiOutput("niche_map_ui"),
                  fluidRow(
                    column(3,sliderInput(inputId = "niche_map_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "niche_map_height", label = "Change height (px)", min = 0, max = 800, value = 400))
                  ),
                  downloadButton("down_niche_map", label = "Download"),
                  actionButton("niche_map_help", label = "Help"),
                  width = 12
                ),
                shinydashboardPlus::box(
                  id = "niche_warning_2",
                  title = "Warning",
                  status = "primary",
                  width = 4,
                  verbatimTextOutput("niche_message_2")
                )
        ),
        tabItem("niche_target",
                shinydashboardPlus::box(
                  status = "primary",
                  id = "niche_3",
                  title = "Ligand heatmap/scale",
                  uiOutput("niche_target_ui"),
                  fluidRow(
                    column(3,sliderInput(inputId = "niche_target_width", label = "Change width (%)", min = 0, max = 100, value = 60)),
                    column(3,sliderInput(inputId = "niche_target_height", label = "Change height (px)", min = 0, max = 800, value = 400))
                  ),
                  downloadButton("down_niche_target", label = "Download"),
                  actionButton("niche_target_help", label = "Help"),
                  width = 12
                ),
                shinydashboardPlus::box(
                  id = "niche_warning_3",
                  title = "Warning",
                  status = "primary",
                  width = 4,
                  verbatimTextOutput("niche_message_3")
                )
        )
      )
    )
  )
)