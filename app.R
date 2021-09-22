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
# library(shinyWidgets)
library(shinyjs)
library(plotly)
library(tidyverse)
library(ggsci)
library(cowplot)
library(patchwork)
library(ggrepel)
library(leaflet)
library(jsonlite)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)



# Load modules ------------------------------------------------------------

source("./predict.R")

# Set ggplot2 default theme -----------------------------------------------

# theme_set(theme_bw())

# Constants -------------------------------------------------------------------------------------------------------

DEFAULT_SPINNER_COLOR = "#3C8DBC"
DEFAULT_LINE_COLOR = "#4393C3"

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


# Load Data ---------------------------------------------------------------

version <- readr::read_rds("data/version.rds")
trend_d <- readr::read_rds("data/GISAID_daily_sequences_count.rds")
trend_m <- readr::read_rds("data/GISAID_sequences_count_trends.rds")
top10_countries_sequences_count <- readr::read_rds("data/top10_countries_sequences_count.rds")
mutations_accumulation_trends <- readr::read_rds("data/mutations_accumutation_trends.rds")
balding_nichols_model_results <- readr::read_rds("data/BN_results_aggregated.rds")
mutations_monthly_count_all <- readr::read_rds("data/mutations_monthly_count_table_all.rds")
mutations_monthly_count_each <- readr::read_rds("data/mutations_monthly_count_table_each.rds")
mutation_position_table <- readr::read_rds("data/mutation_position_table.rds")
country_code_table <- readr::read_rds("data/country_code_table.rds")
mutation_with_annotation <- readr::read_rds("data/mutations_with_annotation.rds")
candidate_mutations <- readr::read_rds("data/candidate_mutations.rds")
candidate_mutations_count_table_country <-  readr::read_rds("data/candidate_mutations_count_table_country.rds")
candidate_mutations_total_table_country <-  readr::read_rds("data/candidate_mutations_total_table_country.rds")
candidate_mutations_count_table_global <-  readr::read_rds("data/candidate_mutations_count_table_global.rds")
candidate_mutations_total_table_global <-  readr::read_rds("data/candidate_mutations_total_table_global.rds")
predict_example_data <- readr::read_rds("data/predict_example_data.rds")
colnames(predict_example_data) <- c("mutation", "position", "freq_prev", "freq_next")


world_sf <- ne_countries(scale = "medium", returnclass = "sf")


# Genes -------------------------------------------------------------------

GENE_LIST <- c("ORF1ab", "S", "ORF3a", "E", "M", "ORF6", "ORF7a", "ORF7b", "ORF8", "N", "ORF10")

MUTATION_LIST <- candidate_mutations %>% distinct(mutation) %>% pull()


# initialize data ---------------------------------------------------------

DATE_LIST <- colnames(candidate_mutations_total_table_global)

ISO3C_LIST <- candidate_mutations_count_table_country %>% distinct(iso3c) %>% pull()

COUNTRY_LIST <- world_sf %>% 
    left_join(distinct(candidate_mutations_count_table_country, iso3c), by=c("iso_a3" = "iso3c")) %>% 
    select(geounit, iso_a3) %>% 
    st_drop_geometry() %>% 
    drop_na() %>% 
    arrange(geounit) %>% 
    deframe()


# Preprocess mutation table for HTML displaying ---------------------------

prepare_mutation_table <- function(mutation_table) {
    mutation_table %>% mutate(min_pvalue = ifelse(pvalue_min == 0, "0", scales::scientific(pvalue_min, digits = 2))) %>% 
        mutate(nt_position = scales::label_comma(accuracy = 1)(position)) %>% 
        mutate(aa_position = scales::label_comma(accuracy = 1)(aa_pos)) %>% 
        mutate(ref = ifelse(nchar(ref) > 6, "...", ref)) %>% 
        mutate(alt = ifelse(nchar(alt) > 6, "...", alt)) %>% 
        mutate(aa_change = ifelse(nchar(aa_change) > 10, "...", aa_change)) %>% 
        select(
            Mutation=mutation, 
            `p-value(min)`=min_pvalue, 
            NT_Position=nt_position, 
            Gene=gene, 
            REF=ref, 
            ALT=alt, 
            # AA_Position=aa_position, 
            AA_Change=aa_change 
            # Class=class
        )
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
            class = "col-sm-12 col-lg-2",
            fluidRow(
                box(
                    width=12,
                    selectInput(
                        "home__gene_select",
                        label = "Gene",
                        choices = GENE_LIST,
                        selected = "S",
                        multiple = FALSE
                    ),
                    selectInput(
                        "home__pvalue_select",
                        label = "p-value",
                        choices = c("<= 1e-5", "<= 1e-6", "<= 1e-7", "<= 1e-8", "<= 1e-9", "<= 1e-10"),
                        selected = "<= 1e-10",
                        multiple = FALSE
                    )
                )
            )
        ),
        
        div(
            class = "col-sm-12 col-lg-6",
            fluidRow(
                box(
                    width = 12,
                    DT::dataTableOutput("home__mutations_table")
                )
            )

        ),

        div(
            class = "col-sm-12 col-lg-4",
            fluidRow(
                box(
                    width = 12,
                    status = "danger",
                    title = "Selected Mutation",
                    descriptionBlock(
                        header = textOutput("home__mutation_selected"),
                        text = textOutput("home__mutation_selected_aachange")
                    )
                ),
                
                box(
                    width = 12,
                    closable = FALSE,
                    collapsible = TRUE,
                    plotOutput("home__global_mutation_count_plot", height = 250)
                ),
                
                box(
                    width = 12,
                    selectInput(
                        "home__country_select",
                        label = "Country",
                        choices = COUNTRY_LIST,
                        selected = "USA",
                        multiple = FALSE
                    ),
                    plotOutput("home__country_mutation_freq_plot", height = 300)
                )
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
            fluidRow(
                div(
                    class = "col-sm-12 col-lg-6",
                    box(
                        width = 12,
                        closable = FALSE,
                        collapsible = TRUE,
                        plotOutput("mutation_count_geographic_distribution", height = "250px") %>% withSpinner(color = DEFAULT_SPINNER_COLOR)
                    )
                ),
                
                div(
                    class = "col-sm-12 col-lg-6",
                    box(
                        width = 12,
                        closable = FALSE,
                        collapsible = TRUE,
                        plotOutput("mutation_freq_geographic_distribution", height = "250px") %>% withSpinner(color = DEFAULT_SPINNER_COLOR)
                    )
                ),
                
                div(
                    class = "col-sm-12 col-lg-6",
                    box(
                        width = 12,
                        closable = FALSE,
                        collapsible = TRUE,
                        plotOutput("mutation_count_geographic_distribution_within_date", height = "250px") %>% withSpinner(color = DEFAULT_SPINNER_COLOR)
                    )
                ),
                
                div(
                    class = "col-sm-12 col-lg-6",
                    box(
                        width = 12,
                        closable = FALSE,
                        collapsible = TRUE,
                        plotOutput("mutation_freq_geographic_distribution_within_date", height = "250px") %>% withSpinner(color = DEFAULT_SPINNER_COLOR)
                    )
                ),
                
                div(
                  class = "col-sm-12 col-lg-12",
                  box(
                      width = 12,
                      closable = FALSE,
                      collapsible = TRUE,
                      plotlyOutput("mutation_monthly_freq_line_plot", height = "300px") %>% withSpinner(color = DEFAULT_SPINNER_COLOR)
                  )
                ),
                
                div(
                  class = "col-sm-12 col-lg-12",
                  box(
                      width = 12,
                      closable = FALSE,
                      collapsible = TRUE,
                      plotOutput("balding_nichols_manhattan_plot", height = "250px") %>% withSpinner(color = DEFAULT_SPINNER_COLOR)
                  )
                ),
                
                div(
                  class = "col-sm-12 col-lg-12",
                  box(
                      width = 12,
                      closable = FALSE,
                      collapsible = TRUE,
                      plotOutput("mutation_freq_bimonthly_scatter", height = "250px") %>% withSpinner(color = DEFAULT_SPINNER_COLOR)
                  )
                )  
            ),
        ),
        div(
            class = "col-sm-4 col-lg-3",
            box(
                width = 12,
                selectInput(
                    "details__gene_select",
                    label = "Select Gene",
                    choices = GENE_LIST,
                    selected = "S",
                    multiple = FALSE
                ),
                selectInput(
                    "details__mutation_select",
                    label = "Select Mutation",
                    choices = MUTATION_LIST,
                    selected = "A23063T",
                    multiple = FALSE
                ),
                selectInput(
                    "details__country_select",
                    label = "Select Country",
                    choices = COUNTRY_LIST,
                    selected = "USA",
                    multiple = FALSE
                ),
                selectInput(
                    "details__date_select",
                    label = "Select Date",
                    choices = DATE_LIST[4:length(DATE_LIST)],
                    selected = DATE_LIST[length(DATE_LIST)],
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
                    class = "col-sm-12 col-lg-12",
                    fileInput(
                        "predict_upload", 
                        label = "Upload File", 
                        multiple = FALSE, 
                        accept = c(".csv", ".tsv", ".xlsx")),
                ) 
            ),
            footer = div(
                actionButton("predict__run_predict", "Predict"),
                actionButton("predict__load_example_data", "Show Example")
            )
        ),
    ),
    # box(
    #     title = "Uploaded File Path",
    #     width = 12,
    #     DT::dataTableOutput("upload_file_path")
    # ),
    hidden(
        fluidRow(
            id = "predict__exmpale_data_container", 
            box(
                width = 12,
                collapsible = TRUE,
                title = "Example data for prediction",
                tableOutput("predict__example_data_table"),
                footer = actionButton("predict__run_predict_with_example_data", "Predict (Example data)")
            )
        )
    ),
    # div(
    #     box(
    #         width = 12,
    #         title = "Debug",
    #         textOutput("predict__debug_output")
    #     )
    # ),
    hidden(
        fluidRow(
            id = "predict__progress_bar_container",
            box(
                width = 12,
                title = NULL,
                headerBorder = FALSE,
                shinyWidgets::progressBar(id = "predict__progress_bar", value = 0, total = 4000, title = "Predicting ...", display_pct = TRUE) 
            )
        )
    ),
    hidden(
        fluidRow(
          id = "predict__result_container",
          box(
              div(
                  align = "center",
                  plotOutput("predict__result_F_estimate_plot", width = 350, height = 300) %>% withSpinner(color = DEFAULT_SPINNER_COLOR)
              )
          ),
          box(
              plotOutput("predict__result_mutation_freq_scatter_plot", height = 300) %>% withSpinner(color = DEFAULT_SPINNER_COLOR)
          ),
          box(
            width = 12,
            plotOutput("predict__result_pvalue_manhattan_plot", height = 300) %>% withSpinner(color = DEFAULT_SPINNER_COLOR)
          ),
          box(
              width = 12,
              DT::dataTableOutput("predict__result_table")
          )
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
    useShinyjs(),
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

    # Footer ------------------------------------------------------------------

    output$version <- renderText(version$version)


    # Home --------------------------------------------------------------------

    home__CandidateMutations <- reactive({
        candidate_mutations %>% 
            filter(gene == input$home__gene_select) %>% 
            filter(pvalue_min <= switch(
                input$home__pvalue_select,
                "<= 1e-5"  = 1e-5, 
                "<= 1e-6"  = 1e-6, 
                "<= 1e-7"  = 1e-7, 
                "<= 1e-8"  = 1e-8, 
                "<= 1e-9"  = 1e-9, 
                "<= 1e-10" = 1e-10
            )) %>% 
            prepare_mutation_table()
    })
    
    home__MutationSelected <- reactive({
        ifelse(
            length(input$home__mutations_table_rows_selected),
            home__CandidateMutations()[input$home__mutations_table_rows_selected,] %>% pull(Mutation),
            "A23403G"
        )
    })
    
    home__MutationSelectedAAChange <- reactive({
        ifelse(
            length(input$home__mutations_table_rows_selected),
            home__CandidateMutations()[input$home__mutations_table_rows_selected,] %>% pull(AA_Change),
            "D614G"
        )
    })
    
    output$home__mutation_selected <- renderText(home__MutationSelected())
    
    output$home__mutation_selected_aachange <- renderText({
        sprintf("%s:%s", input$home__gene_select, home__MutationSelectedAAChange())
    })
    
    
    output$home__mutations_table <- DT::renderDataTable(
        home__CandidateMutations(), 
        selection = 'single',
        options = list(
            pageLength = 25,
            columnDefs = list(list(className = 'dt-center', targets = "_all"))
        )
    )

    
    output$home__global_mutation_count_plot <- renderPlot({
        # geographic heatmap
        mutationCountTable <- candidate_mutations_count_table_country %>% 
            filter(mutation == home__MutationSelected()) %>% 
            gather("date", "count", colnames(candidate_mutations_count_table_country)[-c(1, 2)]) %>% 
            group_by(mutation, iso3c) %>% 
            summarise(total=sum(count)) %>% 
            ungroup()
        mutationCountSF <- world_sf %>% left_join(mutationCountTable, by=c("iso_a3" = "iso3c"))
        ggplot(data = mutationCountSF) +
            geom_sf(aes(fill=total)) +
            scale_fill_distiller(palette = "RdBu", trans = "log10")
    })
    
    output$home__country_mutation_freq_plot <- renderPlot({
        # line plot
        mutationCountTable <- candidate_mutations_count_table_country %>%
            filter(mutation == home__MutationSelected()) %>% 
            filter(iso3c == input$home__country_select) %>% 
            gather("date", "count", colnames(candidate_mutations_count_table_country)[-c(1, 2)]) %>% 
            mutate(date = lubridate::ymd(date, truncated = TRUE))
        
        mutationTotalTable <- candidate_mutations_total_table_country %>% 
            filter(iso3c == input$home__country_select) %>% 
            gather("date", "total", colnames(candidate_mutations_count_table_country)[-c(1, 2)]) %>% 
            mutate(date = lubridate::ymd(date, truncated = TRUE))
        
        table <- mutationCountTable %>% 
            left_join(mutationTotalTable, by=c("iso3c", "date")) %>% 
            filter(date >= as.Date("2020-03-01")) %>% 
            mutate(freq=ifelse(total==0, 0, count / total)) 
        

        scale_factor <- 1 / max(table$total) 
        
        ggplot() +
            geom_line(aes(x = table$date, y=table$freq), color="#EE0000", alpha=0.85, size=1) +
            geom_bar(aes(x = table$date, y=table$total * scale_factor), stat="identity", fill="#3B4992", alpha=0.5) +
            scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m") +
            scale_y_continuous(
                name = "Frequency", labels = scales::percent, limits = c(0, 1),
                sec.axis = sec_axis(~ . / scale_factor, name = "#(Assemblies)", labels = scales::label_number_si())
            ) +
            theme(
                axis.title.x = element_blank(),
                axis.text.x = element_text(angle=45, hjust=0.8, vjust=1.0),
            )
    })
    
    # Details -----------------------------------------------------------------
    
    output$mutation_count_geographic_distribution <- renderPlot({
        mutationCountTable <- candidate_mutations_count_table_country %>% 
            filter(mutation == input$details__mutation_select) %>% 
            gather("date", "count", colnames(candidate_mutations_count_table_country)[-c(1, 2)]) %>% 
            group_by(mutation, iso3c) %>% 
            summarise(total=sum(count)) %>% 
            ungroup()
        
        mutationCountSF <- world_sf %>% left_join(mutationCountTable, by=c("iso_a3" = "iso3c"))
        
        ggplot(data = mutationCountSF) +
            geom_sf(aes(fill=total)) +
            scale_fill_distiller(name = "count", palette = "RdBu", trans = "log10")
    })
    
    output$mutation_freq_geographic_distribution <- renderPlot({
        mutationCountTable <- candidate_mutations_count_table_country %>% 
            filter(mutation == input$details__mutation_select) %>% 
            gather("date", "count", colnames(candidate_mutations_count_table_country)[-c(1, 2)]) %>% 
            group_by(mutation, iso3c) %>% 
            summarise(count=sum(count)) %>% 
            ungroup()
        
        mutationTotalTable <- candidate_mutations_total_table_country %>% 
            gather("date", "total", colnames(candidate_mutations_total_table_country)[-c(1)]) %>% 
            group_by(iso3c) %>% 
            summarise(total=sum(total)) %>% 
            ungroup()
        
        mutationFreqTable <- mutationCountTable %>% 
            left_join(mutationTotalTable, by="iso3c") %>% 
            mutate(freq = ifelse((total == 0), NA, count / total)) %>% 
            mutate(percent = freq *100)
        
        mutationCountSF <- world_sf %>% left_join(mutationFreqTable, by=c("iso_a3" = "iso3c"))
        
        ggplot(data = mutationCountSF) +
            geom_sf(aes(fill=percent)) +
            scale_fill_distiller(name = "percent", palette = "RdBu", limits=c(0, 100))
    })
    
    output$mutation_count_geographic_distribution_within_date <- renderPlot({
        mutationCountTable <- candidate_mutations_count_table_country %>% 
            filter(mutation == input$details__mutation_select) %>% 
            gather("date", "count", colnames(candidate_mutations_count_table_country)[-c(1, 2)]) %>% 
            filter(date == input$details__date_select) %>% 
            ungroup()
        
        mutationCountSF <- world_sf %>% left_join(mutationCountTable, by=c("iso_a3" = "iso3c"))
        
        ggplot(data = mutationCountSF) +
            geom_sf(aes(fill=count)) +
            scale_fill_distiller(name = "count", palette = "RdBu", trans = "log10")
    })
    
    output$mutation_freq_geographic_distribution_within_date <- renderPlot({
        mutationCountTable <- candidate_mutations_count_table_country %>% 
            filter(mutation == input$details__mutation_select) %>% 
            gather("date", "count", colnames(candidate_mutations_count_table_country)[-c(1, 2)]) %>% 
            filter(date == input$details__date_select) %>% 
            ungroup()
        
        mutationTotalTable <- candidate_mutations_total_table_country %>% 
            gather("date", "total", colnames(candidate_mutations_total_table_country)[-c(1)]) %>% 
            filter(date == input$details__date_select)
        
        mutationFreqTable <- mutationCountTable %>% 
            left_join(mutationTotalTable, by="iso3c") %>% 
            mutate(freq = ifelse((total == 0), NA, count / total)) %>% 
            mutate(percent = freq *100)
        
        mutationCountSF <- world_sf %>% left_join(mutationFreqTable, by=c("iso_a3" = "iso3c"))
        
        ggplot(data = mutationCountSF) +
            geom_sf(aes(fill=percent)) +
            scale_fill_distiller(name = "percent", palette = "RdBu", limits=c(0, 100))
    })
    
    output$mutation_freq_monthly_trend <- renderPlot({
        g <- mutations_monthly_count_each %>% 
            filter(mutation == input$home__mutation_select) %>% 
            gather(date, count, colnames(mutations_monthly_count_each)[-1]) %>% 
            left_join(mutations_monthly_count_all) %>% 
            mutate(freq = count / total, date = lubridate::ymd(date, truncated = TRUE)) %>% 
            ggplot(aes(x=date, y=freq)) +
            geom_line(color=DEFAULT_LINE_COLOR, size=0.7) +
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
    
    output$mutation_monthly_freq_line_plot <- renderPlotly({
        date_next <- lubridate::ymd(input$details__date_select, truncated = TRUE)
        date_prev <- date_next - months(1)
        mutationCountTable <- candidate_mutations_count_table_country %>%
            filter(mutation == input$details__mutation_select) %>% 
            filter(iso3c == input$details__country_select) %>% 
            gather("date", "count", colnames(candidate_mutations_count_table_country)[-c(1, 2)]) %>% 
            mutate(date = lubridate::ymd(date, truncated = TRUE))
        
        mutationTotalTable <- candidate_mutations_total_table_country %>% 
            filter(iso3c == input$details__country_select) %>% 
            gather("date", "total", colnames(candidate_mutations_count_table_country)[-c(1, 2)]) %>% 
            mutate(date = lubridate::ymd(date, truncated = TRUE))
        
        table_geounit <- mutationCountTable %>% 
            left_join(mutationTotalTable, by=c("iso3c", "date")) %>% 
            filter(date >= as.Date("2020-03-01")) %>% 
            mutate(freq=ifelse(total==0, 0, count / total)) 
        
        table_global <- mutations_monthly_count_each %>% 
            filter(mutation == input$details__mutation_select) %>% 
            gather(date, count, colnames(mutations_monthly_count_each)[-1]) %>% 
            left_join(mutations_monthly_count_all) %>% 
            mutate(freq = count / total, date = lubridate::ymd(date, truncated = TRUE)) 
        
        g <- select(table_geounit, date, selected_geounit=freq) %>% 
            left_join(select(table_global, date, global=freq)) %>% 
            gather("group", "freq", 2:3) %>% 
            ggplot(aes(x=date, y=freq, group=group, color=group)) +
            # geom_vline(xintercept=date_prev, color="#EE0000", size=1.0, alpha=0.2, linetype=2) +
            # geom_vline(xintercept=date_next, color="#EE0000", size=1.0, alpha=0.2, linetype=2) +
            # geom_rect(aes(xmin=date_prev, xmax=date_next, ymin=-Inf, ymax=Inf), fill="#EE0000", alpha=0.01) +
            geom_line(aes(text=sprintf("Date: %s\nFrequency: %.2f%%", format(date, "%Y-%m"), freq * 100))) +
            scale_color_aaas() +
            scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m") +
            scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
            ylab("Mutant frequency") +
            theme(
                axis.title = element_text(size=9),
                axis.title.x = element_blank(),
                axis.text = element_text(size=8),
                axis.text.x = element_text(angle=45, hjust=0.8, vjust=1.0),
                legend.position = "bottom",
            )
        
        ggplotly(g, tooltip = 'text') %>%  
            layout(legend = list(orientation = "h"), hovermode = 'x unified', height = 300)
    })
    
    output$balding_nichols_manhattan_plot <- renderPlot({
        max_value <- 20
        
        bn_manhattan_data <- balding_nichols_model_results %>% 
            filter(iso3c == input$details__country_select & month_next == input$details__date_select) %>% 
            left_join(mutation_position_table, by="mutation") %>% 
            mutate(log10pvalue=ifelse(pvalue == 0, max_value, -log10(pvalue)))
        
        bn_manhattan_highlight <- bn_manhattan_data %>% filter(mutation == input$details__mutation_select)
        
        ggplot() +
            geom_vline(xintercept=21563, color="#CCCCCC", size=1.0, alpha=0.8, linetype=2) +
            geom_vline(xintercept=25384, color="#CCCCCC", size=1.0, alpha=0.8, linetype=2) +
            geom_point(aes(x=position, y=log10pvalue), size=2, color="#999999", alpha=0.85, data=bn_manhattan_data) +
            geom_point(aes(x=position, y=log10pvalue), size=2, color="#EE0000", data=bn_manhattan_highlight) +
            scale_y_continuous(limits = c(0, max_value)) +
            geom_text_repel(data=bn_manhattan_highlight, mapping=aes(x=position, y=log10pvalue, label=mutation), color="#EE0000", box.padding = 0.5) +
            ylab(expression("-log"[10]~"(p-value)")) +
            xlab("SARS-CoV-2 genome position") +
            theme(
                legend.position = "none"
            )
        
    })
    
    output$mutation_freq_bimonthly_scatter <- renderPlot({
        bn_scatter_data <- balding_nichols_model_results %>% 
            filter(iso3c == input$details__country_select & month_next == input$details__date_select)
        bn_scatter_highlight <- bn_scatter_data %>% filter(mutation == input$details__mutation_select)
        
        ggplot() +
            geom_abline(slope=1, intercept = 0, color="#CCCCCC", linetype=2, size=1) +
            geom_point(aes(x=freq_prev, y=freq_next), data=bn_scatter_data, size=2, alpha=0.8, color="#999999") +
            geom_point(aes(x=freq_prev, y=freq_next), data=bn_scatter_highlight, size=2, alpha=0.8, color="#EE0000") +
            geom_text_repel(aes(x=freq_prev, y=freq_next, label=mutation), data=bn_scatter_highlight, color="#EE0000", box.padding = 0.5) +
            scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
            scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
            xlab("mutation frequency (prev-month)") +
            ylab("mutation frequency (next-month)") +
            coord_fixed() 
    })

    # Statistics --------------------------------------------------------------
    
    output$trend_d_plot <- renderPlotly({
        g <- trend_d %>% filter(group == "slide") %>% 
            ggplot(aes(x=date, y=value)) +
            geom_line(color=DEFAULT_LINE_COLOR, size=0.7) +
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
            geom_bar(stat="identity", fill=DEFAULT_LINE_COLOR, alpha=0.7) +
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
            geom_line(aes(text=sprintf("Date: %s\nCount: %.2f", format(date, "%Y-%m"), mean))) +
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
        ggplotly(g, tooltip = 'text')%>% layout(legend = list(orientation = "h"), hovermode = 'x unified')
    })
    

    # Predict -----------------------------------------------------------------

    observeEvent(input$predict__load_example_data, {
        if (input$predict__load_example_data %% 2 == 1) {
            shinyjs::show("predict__exmpale_data_container")
            shinyjs::hide("predict__result_container")
            updateActionButton(session, inputId = "predict__load_example_data", label = "Hide Example")
        } else {
            shinyjs::hide("predict__exmpale_data_container")
            updateActionButton(session, inputId = "predict__load_example_data", label = "Show Example")
        }

    })
    
    observeEvent(input$predict__run_predict, {
        req(input$predict_upload)
        
        shinyjs::show("predict__progress_bar_container")
        
        shinyWidgets::updateProgressBar(session = session, id = "predict__progress_bar", value = 0, total = 4000, title = "Predicting ...")
        
        uploaded_data <- readxl::read_xlsx(input$predict_upload$datapath)
        colnames(uploaded_data) <- c("mutation", "position", "freq_prev", "freq_next")
        
        uploaded_data %<>% filter_at(vars(c(freq_prev, freq_next)), any_vars(. != 0 ))
        
        predict_result <- covmutit_predict(uploaded_data, session = session)
        
        shinyjs::show("predict__result_container")
        output$predict__result_F_estimate_plot <- renderPlot(predict_result$F_estimate_plot)
        output$predict__result_mutation_freq_scatter_plot <- renderPlot(predict_result$scatter_plot)
        output$predict__result_pvalue_manhattan_plot <- renderPlot(predict_result$manhattan_plot)
        output$predict__result_table <- DT::renderDataTable(
            predict_result$table,
            server = FALSE,
            extensions = 'Buttons',
            selection = "none",
            options = list(
                pageLength = 25,
                columnDefs = list(list(className = 'dt-center', targets = "_all")),
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel')
            ),
            class = "display"
        )
        
        shinyjs::hide("predict__progress_bar_container")
    })
    
    # debugMessage <- reactive("hello world!")
    # 
    # output$predict__debug_output <- renderText({
    #     debugMessage()
    # })
    # 
    # output$upload_file_path <- DT::renderDataTable(input$predict_upload)
    
    predictExampleData <- reactive({
        predict_example_data[sample(nrow(predict_example_data), ceiling(nrow(predict_example_data) / 10)), ] %>% 
            filter_at(vars(c(freq_prev, freq_next)), any_vars(. != 0 )) %>% 
            mutate(
                freq_prev = ifelse(freq_prev == 0 , "0", formatC(freq_prev, format = "e", digits = 2)), 
                freq_next = ifelse(freq_next == 0 , "0", formatC(freq_next, format = "e", digits = 2))
            ) %>% 
            head(n=10) %>% 
            arrange(position)
    })

    output$predict__example_data_table <- renderTable(predictExampleData(), align = "c")
    
    observeEvent(input$predict__run_predict_with_example_data, {
        shinyjs::show("predict__progress_bar_container")
        
        shinyWidgets::updateProgressBar(session = session, id = "predict__progress_bar", value = 0, total = 4000, title = "Predicting ...")
        
        example_data <- predict_example_data[sample(nrow(predict_example_data), ceiling(nrow(predict_example_data) / 100)), ] %>%
            filter_at(vars(c(freq_prev, freq_next)), any_vars(. != 0 ))
        
        example_result <- covmutit_predict(example_data, session = session)
        
        shinyjs::show("predict__result_container")
        output$predict__result_F_estimate_plot <- renderPlot(example_result$F_estimate_plot)
        output$predict__result_mutation_freq_scatter_plot <- renderPlot(example_result$scatter_plot)
        output$predict__result_pvalue_manhattan_plot <- renderPlot(example_result$manhattan_plot)
        output$predict__result_table <- DT::renderDataTable(
            example_result$table,
            server = FALSE,
            extensions = 'Buttons',
            selection = "none",
            options = list(
                pageLength = 25,
                columnDefs = list(list(className = 'dt-center', targets = "_all")),
                dom = 'Bfrtip',
                buttons = c('copy', 'csv', 'excel')
            ),
            class = "display"
        )
        
        shinyjs::hide("predict__progress_bar_container")
    })
    
}

# Run the application
options(shiny.port = 80)
shinyApp(ui = ui, server = server)

