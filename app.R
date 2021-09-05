#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(shinycssloaders)
library(shinyjs)
library(plotly)
library(tidyverse)
library(ggsci)
library(cowplot)
library(patchwork)
library(ggrepel)
library(leaflet)
library(jsonlite)


# Constants -------------------------------------------------------------------------------------------------------

DEFAULT_SPINNER_COLOR = "#3C8DBC"

# Helper functions ------------------------------------------------------------------------------------------------

small_legend <- function(.plot, pointSize = 1, textSize = 8, spaceLegend = 0.8) {
    .plot <- .plot +
        guides(shape = guide_legend(override.aes = list(size = pointSize)),
               color = guide_legend(override.aes = list(size = pointSize))) +
        theme(legend.title = element_text(size = textSize), 
              legend.text  = element_text(size = textSize),
              legend.key.size = unit(spaceLegend, "lines"))
    return(.plot)
}


# Pages -------------------------------------------------------------------


# Header ------------------------------------------------------------------

Header <- dashboardHeader(
    title = span("CoVMutIT", class="logo", style="font-size: 20px;")
)


# Sidebar ----------------------------------------------------------------

Sidebar <- dashboardSidebar(
    sidebarMenu(
        id = "main_sidebar_tabs",
        menuItem("Home", tabName = "home", icon = icon("home")),
        menuItem("Details", tabName = "details", icon = icon("globe")),
        menuItem("Statistics", tabName = "statistics", icon = icon("chart-line")),
        menuItem("Predict", tabName = "predict", icon = icon("compass")),
        menuItem("Doucments", tabName = "documents", icon = icon("question-circle")),
        menuItem("About us", tabName = "about", icon = icon("address-card"))
    )
)


# Body --------------------------------------------------------------------


# HomeTab ---------------------------------------------------------------------------------------------------------


HomeTab <- tabItem(
    tabName = "home",
    fluidRow(
        div(
            class = "col-sm-12 col-lg-8",
            box(
                width=12,
                leafletOutput("main_map")
            )
        ),
        div(
            class = "col-sm-12 col-lg-4",
            box(
                width=12,
                selectInput(
                    "home__gene_select",
                    label = "Gene",
                    choices = c("S"),
                    selected = "S",
                    multiple = FALSE
                ),
                selectInput(
                    "home__mutation_select",
                    label = "Mutation",
                    choices = c("A23403G"),
                    selected = "A23403G",
                    multiple = FALSE
                ),
                actionButton("home__show_details", "Show details")
            )
        ),
        div(
            class = "col-sm-12",
            box(
                title = "line plot showing trend of mutant frequency changes",
                plotlyOutput("mutation_freq_monthly_trend")
                
            )
        )
    )
)


# DetailsTab ------------------------------------------------------------------------------------------------------


DetailsTab <- tabItem(
    tabName = "details",
    fluidRow(
        div(
            class = "col-sm-8 col-lg-9",
            box(
                width = 12,
                title = "Basic information",
                "hello world!"
            ),
            box(
                width = 12,
                closable = FALSE,
                collapsible = TRUE,
                "hello world!"
            ),
            box(
                width = 12,
                closable = FALSE,
                collapsible = TRUE,
                "hello world!"
            ),
            box(
                width = 12,
                closable = FALSE,
                collapsible = TRUE,
                "hello world!"
            )
        ),
        div(
            class = "col-sm-4 col-lg-3",
            box(
                width = 12,
                selectInput(
                    "details__gene_select",
                    label = "Select Gene",
                    choices = c("S"),
                    selected = "S",
                    multiple = FALSE
                ),
                selectInput(
                    "details__mutation_select",
                    label = "Select Gene",
                    choices = c("S"),
                    selected = "S",
                    multiple = FALSE
                ),
                selectInput(
                    "details__country_select",
                    label = "Select Gene",
                    choices = c("S"),
                    selected = "S",
                    multiple = FALSE
                ),
                selectInput(
                    "details__date_select",
                    label = "Select Gene",
                    choices = c("S"),
                    selected = "S",
                    multiple = FALSE
                )
            )
            
        )
    )
)


# StatisticsTab ---------------------------------------------------------------------------------------------------



StatisticsTab <- tabItem(
    tabName = "statistics",
    fluidRow(
        div(
            class = "col-sm-12",
            box(
                width=12,
                solidHeader = FALSE,
                background = NULL,
                status = "primary",
                fluidRow(
                    column(
                        width=3,
                        descriptionBlock(
                            header = "2,518,369",
                            text = "Total assemblies",
                            rightBorder = TRUE,
                            marginBottom = FALSE
                        )
                    ),
                    column(
                        width=3,
                        descriptionBlock(
                            header = "2,518,369",
                            text = "Total assemblies",
                            rightBorder = TRUE,
                            marginBottom = FALSE
                        )
                    ),
                    column(
                        width=3,
                        descriptionBlock(
                            header = "2,518,369",
                            text = "Total assemblies",
                            rightBorder = TRUE,
                            marginBottom = FALSE
                        )
                    ),
                    column(
                        width=3,
                        descriptionBlock(
                            header = "2,518,369",
                            text = "Total assemblies",
                            rightBorder = TRUE,
                            marginBottom = FALSE
                        )
                    )
                )
            )
        )
    ),
    fluidRow(
        
        ###
        div(
            class = "col-lg-4 col-md-6 col-sm-12", 
            box(
                width = 12,
                title = "Statistics",
                plotlyOutput("trend_d_plot") %>% withSpinner(color = DEFAULT_SPINNER_COLOR)
            )
        ),
        
        ###
        div(
            class = "col-lg-4 col-md-6 col-sm-12", 
            box(
                width = 12,
                title = "Statistics",
                plotlyOutput("trend_m_plot") %>% withSpinner(color = DEFAULT_SPINNER_COLOR)
            )
        ),
        
        ###
        div(
            class = "col-lg-4 col-md-6 col-sm-12", 
            box(
                width = 12,
                title = "Statistics",
                plotOutput("top10_countries_sequences_count_piechart") %>% withSpinner(color = DEFAULT_SPINNER_COLOR)
            )
        ),
        
        ###
        div(
            class = "col-lg-4 col-md-6 col-sm-12", 
            box(
                width = 12,
                title = "Statistics",
                plotlyOutput("mutations_accumulation_trends") %>% withSpinner(color = DEFAULT_SPINNER_COLOR)
            )
        )
    ) 
)


# PredictTab --------------------------------------------------------------

PredictTab <- tabItem(
    tabName = "predict",
    fluidRow(
        box(
            title = "Predict",
            width = 12,
            fluidRow(
                div(
                    class = "col-sm-12",
                    "[description holdplace]",
                ),
                div(
                    class = "col-sm-12 col-lg-6",
                    fileInput(
                        "predict_upload", 
                        label = "Upload File", 
                        multiple = FALSE, 
                        accept = c(".csv", ".tsv", ".xlsx")),
                ) 
            ),
            footer = actionButton("run_predict", "Predict")
        )
        
        
    )
    
)


# DocumentsTab ----------------------------------------------------------------------------------------------------



DocumentsTab <- tabItem(
    tabName = "documents",
    fluidRow(
       box(
           "hello world",
           width = 12
       ) 

    )
    
)


# AboutUsTab ------------------------------------------------------------------------------------------------------


AboutUsTab <- tabItem(
    tabName = "about",
    fluidRow(
        ### 
        div(
            class = "col-lg-3 col-md-6 col-sm-12",
            box(
                width = 12,
                boxProfile(
                    image = "user-profile-m.png",
                    title = "Huang Kai, Ph.D.",
                    subtitle = "Professor of Cardiology",
                    bordered = TRUE,
                    boxProfileItem(
                        title = "Email",
                        description = "huangkai1@hust.edu.cn"
                    )
                )
            )
        ),
        
        ### 
        div(
            class = "col-lg-3 col-md-6 col-sm-12",
            box(
                width = 12,
                boxProfile(
                    image = "user-profile-m.png",
                    title = "Wang Chaolong, Ph.D.",
                    subtitle = "Professor of Bioinformatics",
                    bordered = TRUE,
                    boxProfileItem(
                        title = "Email",
                        description = "chaolong@hust.edu.cn"
                    )
                )
            )
        ),
        
        ###
        div(
            class = "col-lg-3 col-md-6 col-sm-12",
            box(
                width = 12,
                boxProfile(
                    image = "user-profile-m.png",
                    title = "Liu Wei, Ph.D.",
                    subtitle = "Routine maintenance and functional update",
                    bordered = TRUE,
                    boxProfileItem(
                        title = "Email",
                        description = "liuweiathust@foxmail.com"
                    )
                )
            )
        ),
        
        ###
        div(
            class = "col-lg-3 col-md-6 col-sm-12",
            box(
                width = 12,
                boxProfile(
                    image = "user-profile-f.png",
                    title = "Li Yanze, Ph.D. Candidate",
                    subtitle = "Routine maintenance",
                    bordered = TRUE,
                    boxProfileItem(
                        title = "Email",
                        description = "d202081632@hust.edu.cn"
                    )
                )
            )
        )
        
    )
)

Body <- dashboardBody(
    tabItems(
        HomeTab,
        DetailsTab,
        PredictTab,
        DocumentsTab,
        StatisticsTab,
        AboutUsTab
    )
)



# Footer ------------------------------------------------------------------

Footer <- dashboardFooter(
    left = span("Version: ", textOutput("version", inline = TRUE), style = "font-size: 12px; color: #888888;"),
    right = span(
        "Clinic Center of Human Gene Research, ",
        a("Union Hospital", href="http://www.whuh.com/", target="_blank"), 
        ", ",
        a("Tongji Medical College", href="http://www.tjmu.edu.cn/", target="_blank"),
        ", ",
        a("HUST", href="https://www.hust.edu.cn/", target="_blank"),
        ", ",  
        "China", 
        style = "font-size: 12px"
    )
)


# Load Data ---------------------------------------------------------------

version <- readr::read_rds("data/version.rds")
trend_d <- readr::read_rds("data/GISAID_daily_sequences_count.rds")
trend_m <- readr::read_rds("data/GISAID_sequences_count_trends.rds")
top10_countries_sequences_count <- readr::read_rds("data/top10_countries_sequences_count.rds")
mutations_accumulation_trends <- readr::read_rds("data/mutations_accumutation_trends.rds")
balding_nichols_model_resutls <- readr::read_rds("data/BN_results_aggregated.rds")
mutations_monthly_count_all <- readr::read_rds("data/mutations_monthly_count_table_all.rds")
mutations_monthly_count_each <- readr::read_rds("data/mutations_monthly_count_table_each.rds")

# -------------------------------------------------------------------------

# Define UI for application that draws a histogram
ui <- dashboardPage(
    header = Header,
    sidebar = Sidebar,
    body = Body,
    footer = Footer,
    skin = "blue-light"
)



# Define server logic required to draw a histogram
server <- function(input, output, session) {
    observeEvent(input$home__show_details, {
        updateTabItems(session, "main_sidebar_tabs", "details")
    })
    
    output$version <- renderText(version$version)
    
    output$trend_d_plot <- renderPlotly({
        g <- trend_d %>% filter(group == "slide") %>% 
            ggplot(aes(x=date, y=value)) +
            geom_line(color=pal_aaas()(1), size=0.7) +
            scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m") +
            scale_y_continuous(labels = scales::label_number_si()) +
            ylab("7-Day Moving average") +
            theme(
                axis.title = element_text(size=9),
                axis.title.x = element_blank(),
                axis.text = element_text(size=8),
                axis.text.x = element_text(angle=45, hjust=0.8, vjust=1.0),
            )
        ggplotly(g)
    })
    
    output$trend_m_plot <- renderPlotly({
        g <- trend_m %>%  
            mutate(date=lubridate::ymd(date, truncated = TRUE)) %>% 
            ggplot(aes(x=date, y=count)) +
            geom_bar(stat="identity", fill=pal_aaas()(1), alpha=0.7) +
            scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m") +
            scale_y_continuous(labels = scales::label_number_si()) +
            ylab("Monthly average") +
            theme(
                axis.title = element_text(size=9),
                axis.title.x = element_blank(),
                axis.text = element_text(size=8),
                axis.text.x = element_text(angle=45, hjust=0.8, vjust=1.0),
            )
        ggplotly(g)
    })
    
    output$top10_countries_sequences_count_piechart <- renderPlot({
        g <- top10_countries_sequences_count %>% 
            ggplot(aes(x=2, y=count, fill=factor(label, levels=label))) +
            geom_bar(stat="identity") +
            coord_polar("y", start=0) +
            scale_fill_manual(
                values=c(pal_npg()(10), "#CCCCCC"), 
                name=paste("Countries (GISAID: ", format(version$last_update, "%b/%d/%Y"), ")")) +
            ylab("") +
            xlim(0.2, 2.5) +
            theme_void() +
            theme(
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                axis.line.x = element_blank(),
                axis.ticks.x = element_blank(),
                legend.text = element_text(family = "mono"),
            )
        return(small_legend(g, textSize = 9, pointSize = 8))
    })
    
    output$mutations_accumulation_trends <- renderPlotly({
        g <- mutations_accumulation_trends %>%  
            ggplot(aes(x=date, y=mean, group=group, color=group)) +
            geom_line() +
            geom_point(size=0.8)+
            geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.3,
                          position=position_dodge(0.05)) +
            ylab("#(mutantions)") +
            scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m") +
            scale_color_aaas() +
            theme(
                axis.title = element_text(size=9),
                axis.title.x = element_blank(),
                axis.text = element_text(size=8),
                axis.text.x = element_text(angle=45, hjust=0.8, vjust=1.0),
                legend.position = "bottom",
            )
        ggplotly(g)%>% layout(legend = list(orientation = "h"))
    })
    
    output$main_map <- renderLeaflet({
        leaflet() %>% addTiles()
    })
    
    output$mutation_freq_monthly_trend <- renderPlotly({
        g <- mutations_monthly_count_each %>% 
            filter(mutation == input$home__mutation_select) %>% 
            gather(date, count, colnames(mutations_monthly_count_each)[-1]) %>% 
            left_join(mutations_monthly_count_all) %>% 
            mutate(freq = count / total, date = lubridate::ymd(date, truncated = TRUE)) %>% 
            ggplot(aes(x=date, y=freq)) +
            geom_line(color=pal_aaas()(1), size=0.7) +
            scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m") +
            scale_y_continuous(labels = scales::percent) +
            ylab("Mutant frequency") +
            theme(
                axis.title = element_text(size=9),
                axis.title.x = element_blank(),
                axis.text = element_text(size=8),
                axis.text.x = element_text(angle=45, hjust=0.8, vjust=1.0),
            )
        ggplotly(g)
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
