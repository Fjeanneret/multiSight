#' Select only data values from diablo features selected. 
#'
#' @param input,output,session Internal parameters for shiny
#' @param query character vector to query pubmed database
#' @param featMethod character string to specify 
#' feature selection method (e.g. "diablo")
#'
#' @return uiOutput, a html content to embed in client side.
#'
#' @importFrom shiny renderUI a
#' @importFrom easyPubMed get_pubmed_ids fetch_pubmed_data  articles_to_list 
pubmedInsight <- function(session, input, output, query, featMethod)
{ 
    if (query != "")
    {
        message("PubMed query: ", query)
        suppressWarnings(my_query <- 
                             get_pubmed_ids(pubmed_query_string = query))
        obj$networkInference[[featMethod]]$biblio$request <- 
            my_query$QueryTranslation
        if (my_query$Count != 0)
        {
              my_data <- fetch_pubmed_data(my_query, 
                                           encoding = "ASCII",  
                                           retmax = 10)
              my_abstracts_txt <- fetch_pubmed_data(my_query, 
                                                    format = "abstract")
              articleList <- articles_to_list(my_data) # list of articles
              articleNb <- length(articleList)
              if (articleNb > 5) articleNb <- 5 # print 5 articles
              output$pubmed <- renderUI({
                ## for each article creates a box with title and abstract
                  displayPubmed(articleList, articleNb, featMethod)
                
              })
        }else
        {
          output$pubmed <- renderUI({
            box(width = 10, 
                title = "0 article found, try with another query", 
                collapsed = TRUE
            )
          })
        }
    }
}

#' UI Pubmed articles 
#'
#' @description Selects and displays pubmed articles according to
#' nnumber provided.
#'
#' @param articleList articleList from articles_to_list() and fetch_pubmed_data
#' functions
#' @param articleNb number of articles to display
#' @param featMethod character string to specify 
#' feature selection method (e.g. "diablo")
#' @return uiOutput, a html content to embed in client side.
#'
#' @noRd
#'
#' @importFrom shiny renderUI a
#' @importFrom shinydashboard box
#' @importFrom easyPubMed custom_grep article_to_df 
displayPubmed <- function(articleList, articleNb, featMethod)
{
    lapply(seq(1, articleNb), function(i){
        article <- articleList[[i]]
        title <-  custom_grep(article, "ArticleTitle", "char" )
        ## A dataframe to retrieve abstract without html tags
        article_df <- article_to_df(article, 
                                    max_chars = -1, 
                                    getAuthors = FALSE)
        abstract <- article_df$abstract
        ## creating a universal url with doi 
        doi <- article_df$doi
        url <- paste0("https://doi.org/", doi)
        html_url <- a("Read online", href = url, target="_blank")
        title_url <- list(title = title, url = url)
        ## Add in obj object pubmed results
        obj$networkInference[[featMethod]]$biblio$title <-
            c(obj$networkInference[[featMethod]]$biblio$title, 
              title_url$title)
        obj$networkInference[[featMethod]]$biblio$url <-
            c(obj$networkInference[[featMethod]]$biblio$url,
              title_url$url)
        ## UI
        box(width = 10, 
            title = HTML(title),
            collapsible = TRUE, 
            collapsed = TRUE,
            HTML(abstract),
            html_url
        )
    })
}
