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
library(shinyjs)


# Pages -------------------------------------------------------------------


# Header ------------------------------------------------------------------

Header <- dashboardHeader(
    title = span("CoVMutIT", class="logo", style="font-size: 20px;")
)


# Sidebar ----------------------------------------------------------------

Sidebar <- dashboardSidebar(
    sidebarMenu(
        menuItem("Home", tabName = "home", icon = icon("home")),
        menuItem("Details", tabName = "details", icon = icon("globe")),
        menuItem("Statistics", tabName = "statistics", icon = icon("chart-line")),
        menuItem("Doucments", tabName = "documents", icon = icon("question-circle")),
        menuItem("About us", tabName = "about", icon = icon("address-card"))
    )
)


# Body --------------------------------------------------------------------

HomeTab <- tabItem(
    tabName = "home",
    fluidRow(
        box(
            "home"
        )
    )
)

DetailsTab <- tabItem(
    tabName = "details",
    fluidRow(
        box(
            "details"
        )
    )
)

StatisticsTab <- tabItem(
    tabName = "statistics",
    fluidRow(
        box(
            "statistics"
        )
    ) 
)

DocumentsTab <- tabItem(
    tabName = "documents",
    fluidRow(
       box(
           "hello world",
           width = 12
       ) 

    )
    
)

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


# -------------------------------------------------------------------------

# Define UI for application that draws a histogram
ui <- dashboardPage(
    header = Header,
    sidebar = Sidebar,
    body = Body,
    footer = Footer,
    skin = "red-light"
)

version <- readr::read_rds("data/version.rds")

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    output$version <- renderText(version$version)
}

# Run the application 
shinyApp(ui = ui, server = server)
