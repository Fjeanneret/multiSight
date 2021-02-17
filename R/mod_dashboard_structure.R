#' dashboard_structure UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd 
#'
#' @importFrom shiny NS tagList column sliderInput
#' textInput actionButton htmlOutput fluidRow icon h1
#' @import shinydashboard  
mod_dashboard_structure_ui <- function(id){
  ns <- NS(id)
  dashboardPage(
    dashboardHeader(title = "multiSight"),
    
    dashboardSidebar(    
      sidebarMenu(
        menuItem("Home", 
                 tabName = "home", 
                 icon = icon("home")),
        menuItem("Classification", 
                 tabName = "models", 
                 icon = icon("bullseye")),
        menuItem("Biological Insights",
                 tabName = "understand", 
                 icon = icon("database")),
        menuItem("Assumptions", 
                 tabName = "hypothesis", 
                 icon = icon("book"))
      )),
    
    dashboardBody(
      fluidRow(
        tabItems(
          # Home tab content
          tabItem(
            box(title = h1("Home"),
                width = 10, status = "warning",
                  "This Home tab is dedicated to input omic datasets and 
                    to save results"
            ),
            # h1("Home", align = "center"), 
            tabName = "home",
            # mod_home_ui("home_ui_1")
            mod_user_input_ui("user_input_ui_1")
          ),
          # Discrimination tab content
          tabItem(
            box(title = h1("Classification"),
                width = 10, status = "warning",
                "This Classification tab is dedicated to classification models
                performances and features selected."
            ),
            tabName = "models",
            mod_MLmodels_ui("MLmodels_ui_1")
          ),
          # Understand tab content
          tabItem(
            box(title = h1("Biological Insights"),
                width = 10, status = "warning",
                "This Biological Insights tab is dedicated to single-omic and 
                multi-omic functional enrichment with diablo or DESeq2 
                selected features.\n In case of 2 or more omic datasets 
                enriched, Stouffer's values and relative enrichment map for 
                each pathways or ontology are provided to see at a glance which
                ones are meaningful between several omics data.
                "
            ),
            tabName = "understand", 
            mod_understand_ui("understand_ui_1")
          ),
          # Hypothesis tab content
          tabItem(
            box(title = h1("Assumptions"),
                width = 10, status = "warning",
                "This Assumptions tab is dedicated to Newtwork inferences and 
                a integrated PubMed query module.",
                "Networks are build according to 3 methods: correlation, 
                partial correlation and mutual information (see ", 
                tags$a(href="https://doi.org/10.3389/fgene.2019.00535", 
                       "Inferring Interaction Networks From Multi-Omics Data", 
                       target="_blank"), ")."
                
            ),
                  
              tabName = "hypothesis", 
              box(title ="Diablo features", width = 10, 
                  status = "warning", collapsible = TRUE,
                  solidHeader = TRUE, 
                  actionButton(ns("diabloNetworkStart"),
                         "Build networks", icon("project-diagram")),
                  rep_br(times = 2),
                  mod_hypothesisGenerator_ui(ns("hypothesisGenerator_ui_1"))
            ),
            box(title ="Biosigner features", width = 10, 
                status = "warning", collapsible = TRUE,
                solidHeader = TRUE,
                actionButton(ns("biosignerNetworkStart"),
                         "Build networks", icon("project-diagram")),
                rep_br(times = 2),
                mod_hypothesisGenerator_ui(ns("hypothesisGenerator_ui_2"))
          )
        )
      )
    )
  )
)
}
    
#' dashboard_structure Server Function
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#' @param obj R6 object to wrap all data from different analysis.
#'
#' @importFrom shiny observeEvent req 
#'
#' @noRd 
mod_dashboard_structure_server <- function(input,
                                           output, 
                                           session)
{
  
    obj <- obj
    ns <- session$ns
    
    observeEvent(input$diabloNetworkStart,{
      req(length(obj$classification$diabloResult$biosignature) > 0)
      biosignature <- obj$classification$diabloResult$biosignature
      callModule(mod_hypothesisGenerator_server, 
                 "hypothesisGenerator_ui_1", 
                 obj, 
                 biosignature,
                 "diabloNetwork")
    })
    
    observeEvent(input$biosignerNetworkStart, {
      req(length(obj$classification$biosignerResult$biosignature) > 0)
      biosignature <- obj$classification$biosignerResult$biosignature
      callModule(mod_hypothesisGenerator_server, 
                 "hypothesisGenerator_ui_2", 
                 obj, 
                 biosignature,
                 "biosignerNetwork")
    })
  # callModule(mod_user_input_server, "user_input_ui_1")  
}
    
## To be copied in the UI
# mod_dashboard_structure_ui("dashboard_structure_ui_1")
    
## To be copied in the server
# callModule(mod_dashboard_structure_server, "dashboard_structure_ui_1")
 
