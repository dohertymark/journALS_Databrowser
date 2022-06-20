source("/Users/markdoherty/Documents/GitHub/journALS/App/functions.R")
library(shiny)
require(shinyWidgets)
require(shinyjs)

js <- "$(document).ready(function(){$('#tabset').tabCollapse();});"

ui <- 

useShinyjs()
# These tags ensure images scale correctly on different screens
tagList(tags$head(
  # tags$style(
  # type="text/css",
  # # This is really important to scale plots with screen size 
  # "#server.comparison_plots img {max-width: 80%; width: 80%; height: auto}
  # #server.age_plot img {max-width: 100%; width: 100%; height: auto}
  # #server.age_plot_ecdf img {max-width: 100%; width: 100%; height: auto}
  # #server.age_plot_kde img {max-width: 100%; width: 100%; height: auto}
  # #server.pedigree_plot img {max-width: 100%; width: 100%; height: auto}
  # #server.heterogeneity_plots img {max-width: 100%; width: 100%; height: auto}
  # #server.acmg_pathogenicity_plots img {max-width: 100%; width: 100%; height: auto;image-rendering: auto;}
  # #server.regions_analysis_plot img {max-width: 100%; width: 100%; height: auto}
  # #server.summary.proportion_explained.plot img {max-width: 100%; width: 100%; height: auto}"


  # #"#faf95 img {max-width: 80%; width: 80%; height: auto}
  # #age_plot_ecdf img {max-width: 100%; width: 100%; height: auto}
  # #age_plot_kde img {max-width: 100%; width: 100%; height: auto}
  # #pedigree_plot img {max-width: 100%; width: 100%; height: auto}
  # #heterogeneity_plots img {max-width: 100%; width: 100%; height: auto}
  # #regions_analysis_plot img {max-width: 100%; width: 100%; height: auto}"
  # ),
  tags$link(rel="shortcut icon", href="favicon.ico"),
  # This is required to get size of user screen to scale images
  tags$script('$(document).on("shiny:connected", function(e) {
    Shiny.onInputChange("innerWidth", window.innerWidth);
    });
    $(window).resize(function(e) {
      Shiny.onInputChange("innerWidth", window.innerWidth);
      });
      '),
  tags$style(HTML(
    "a.js-tabcollapse-panel-heading {
     display: block;
     text-align: center;
     }"
     )
  ),
  tags$script(src = "bootstrap-tabcollapse.js"),
  tags$script(HTML(js))),

navbarPage(h4("journALS"),
  collapsible = TRUE,
  theme = "modified.css",
  windowTitle = "journALS",


  ###################
  ### Home Page   ###
  ###################

  tabPanel(h4("Home"),

  br(),
  column(8,offset=2,
    h2(mainPanel(htmlOutput("home_header_text"),width=12,align="center")),
    br(),
    br(),
    h3(mainPanel(htmlOutput("home_body_text"),width=12,align="justify"))),
  column(12,
    br(),
    br(),    
    br(),
    br(),    
    br(),
    br(),    
    br(),
    br(),
    hr(),
    print("journALS is funded by Science Foundation Ireland and FutureNeuro"),
    br(),
    print("This work was supported by TCHPC (Research IT, Trinity College Dublin)"),
    br(),
    br(),
    br(),
    print("This website is currently under development and Beta testing"))
     ),


  
  #######################
  ### Variant Browser ###
  #######################

  tabPanel(h4("Variant Browser"),
    #id = "form",
    br(),
    select_variant_browser_input(),
    br(),
    br(),
    tabsetPanel(type = "tabs",

      tabPanel(h4("Pathogenicity"),
        fluidRow(
          column(6,offset=3,
            plotOutput("server.acmg_pathogenicity_plots")
            ))
        ),
      tabPanel(h4("General Information"),
        br(),
        variant_browser_general_information_input(),        
        DT::dataTableOutput("variant_table_single_variant"),
        br(),
        ),
      tabPanel(h4("Phenotype Information"),
        br(),
        variant_browser_phenotype_information_input(),
        DT::dataTableOutput("variant_table_single_variant_phenotype"),
        br(),
        br(),
        h6(mainPanel(htmlOutput("phenotype_table_text_line1"),width=12))
        ),
      tabPanel(h4("Geographic Heterogeneity"),
        br(),
        ui.heterogeneity_plots(),
        br(),
        br(),        
        br(),
        br(),
        h6(mainPanel(htmlOutput("heterogeneity_caveat1"),width=12)),
        h6(mainPanel(htmlOutput("heterogeneity_caveat2"),width=12))
        ),
      tabPanel(h4("Age of Onset"),
        br(),
        #variant_age_tab_input(),
        #ui.age_plot.output(),
        ui.age_plot(),
        h6(mainPanel(htmlOutput("age_text_line2"),width=12)),
        h6(mainPanel(htmlOutput("age_text_line1"),width=12)),
        h6(mainPanel(htmlOutput("age_text_line3"),width=12)),

        ),
      tabPanel(h4("Pedigrees"),
        br(),
        variant_pedigree_tab_input(),
        plotOutput("server.pedigree_plot",inline=T),
        # Adding breaks to give more room for larger pedigrees
        br(),
        br(),
        br(),
        br(),
        br(),
        br(),
        br(),
        br(),
        h6(mainPanel(htmlOutput("pedigree_text_line0"),width=12)),
        h6(mainPanel(htmlOutput("pedigree_text_line1"),width=12)),
        h6(mainPanel(htmlOutput("pedigree_text_line2"),width=12)),
        h6(mainPanel(htmlOutput("pedigree_text_line3"),width=12)),
        fluidRow(
          column(1,br(),downloadButton("download.pedigree_plot", label=h6("Save Figure"))),
          column(1,br(),downloadButton("downloadPedFile", label=h6("Save Ped File")))
          )
        )
      )
    ),

  #######################
  ### Region Browser  ###
  #######################
  tabPanel(h4("Region Browser"),
    chooseSliderSkin(c("Flat"),"#01184E"),
    br(),
    select_region_browser_input(),
    br(),
    br(),
    tabsetPanel(type = "tabs",
      tabPanel(h4("Individuals"),
        br(),
        DT::dataTableOutput("region_browser_individuals_table")),
      tabPanel(h4("Analysis"),
        br(),
        ui.regions_analysis_plot(),
        h6(mainPanel(htmlOutput("regions_analysis_text_line1"),width=12)),
        conditionalPanel(
          "input.regions_history_selection=='Overall'",
          h6(mainPanel(htmlOutput("regions_analysis_text_line2"),width=12)))
        ),
      tabPanel(h4("Population Studies"),
        br(),
        DT::dataTableOutput("region_browser_pop_studies_table"))

      )
    ),

  ####################
  ### Gene Browser ###
  ####################
  tabPanel(h4("Gene Browser"),
    chooseSliderSkin(c("Flat"),"#01184E"),
    fluidRow(column(8, offset=2,align="center",
      selectizeInput(
        "gene", 
        h4("Select Gene of Interest:",style="font-weight: bold;"),
        choices=NULL,
        options=list(
          maxItems=1,
          placeholder='Select Gene',
          create=TRUE,
          #class="dropdown-menu"),selected="TARDBP"
          class="dropdown-menu",selected="TARDBP")
        )
      )
    ),
    br(),
    tabsetPanel(type = "tabs",

      tabPanel(h4("Comparison Plot"),
        ui.comparison_plots(),
        #gene_comparison_output(),
        #br(),
        #gene_comparison_input()
        ),
      tabPanel(h4("Gene Plot"),
        ui.gene_plot(),
        ),
      tabPanel(h4("Variant Table"),
        br(),
        gene_table_input(),
        h5(DT::dataTableOutput("variant_table"))
        )
      )
    ),

  #######################
  ### Summary         ###
  #######################
  tabPanel(h4("Summary"),
    br(),
    column(10,offset=1,
      h2(mainPanel(htmlOutput("summary.header.text"),width=12,align="left")),
      br(),
      br(),
      br(),
      h4(mainPanel(htmlOutput("summary.overview.text1"),width=12,align="justify")),
       br(),
       br(),
       fluidRow(
        imageOutput("server.summary.proportion_explained.plot"),width=12,offset=1,align="left"),
        hr(),     
    selectizeInput("summary_gene_selection", "Select Gene with Pathogenic or Likely Pathogenic Variant:",
      choices=NULL,
      options=list(
        maxItems=1,
        placeholder='Select Gene:',
          create=TRUE)),#align="left"
   h4(mainPanel(htmlOutput("summary.gene.text1"),width=12,align="justify"))

      )
      ),
  #######################
  ### Downloads       ###
  #######################
  tabPanel(h4("Downloads"),
    verbatimTextOutput("That")
    ),
  #######################
  ### Annotate        ###
  #######################
  tabPanel(h4("Annotate"),
    br(),
    column(8,offset=2,
      h2(mainPanel(htmlOutput("upload_header_text"),width=12,align="center")),
      br(),
      br(),
      br(),
      br()),
    fluidRow(
      column(6,align="right",
        fileInput("upload_file", "Choose identifier.txt file",
          multiple = TRUE,
          accept = c("text",
           "text/plain",
           ".txt"))),
      column(6,align="left",downloadButton("downloadUploadResults", label=h6("Download journALS annotations"))
        )),
    br(),
    br(),
    column(8,offset=2,
      h4(mainPanel(htmlOutput("upload_body_text"),width=12,align="justify"))
      ),
    fluidRow(
      column(6,align="right",
        br(),
        br(),
        downloadButton("downloadExampleTxt", label=h6("Download example.txt"))
        ),
      column(6,align="left",
        br(),
        br(),
        downloadButton("downloadHeaderExplanation", label=h6("Download column explanations"))
        )
      )
    ),
  #######################
  ### About           ###
  #######################
  tabPanel(h4("About"),
    br(),
    column(8,offset=2,
      h2(mainPanel(htmlOutput("about.FAQ.header"),width=12,align="left")),
      br(),
      br(),
      br(),
      h3(mainPanel(htmlOutput("about.FAQ.citation.header"),width=12,align="justify")),
      br(),
      h4(mainPanel(htmlOutput("about.FAQ.citation.body"),width=12,align="justify")),
      h3(mainPanel(htmlOutput("about.FAQ.gnomAD_AFs.header"),width=12,align="justify")),
      br(),
      h4(mainPanel(htmlOutput("about.FAQ.gnomAD_AFs.body"),width=12,align="justify")),
      h3(mainPanel(htmlOutput("about.FAQ.resources.header"),width=12,align="justify")),
      br(),
      h4(mainPanel(htmlOutput("about.FAQ.resources.body"),width=12,align="justify"))


      )

    )
  )
)
