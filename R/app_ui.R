#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @importFrom shiny tagList
#' @noRd
app_ui <- function(request) {
    tagList(
        # List the first level UI elements here 
        mod_dashboard_structure_ui("dashboard_structure_ui_1"),
    )
}

