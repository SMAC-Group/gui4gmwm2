# source("R/plot_wv_and_datasheet.R")
# this is a specific branch of package gmwm!!!, namely the branch gmwm2 !
# remotes::install_github(repo = "https://github.com/SMAC-Group/gmwm", ref = "gmwm2")
# remotes::install_github(repo = "https://github.com/SMAC-Group/wv")



library(wv)
library(ggplot2)
library(gmwm)
library(scales)
library(reshape)
library(shiny)
library(shinydashboard)
library(bslib)


# load data wavelet
load("data/imudata.RData")

# source
source("my_plot_wvar.R")
source("my_plot_gmwm.R")
source("transform_parameters.R")



################
# CONSTANTS
# Purpose: Centralized defaults and dataset-specific parameters used by UI,
#          transformations, and datasheet computations.
################
const.RENDER_PLOT_WIDTH <- 1000
const.RENDER_PLOT_HEIGHT <- 600
const.RENDER_PLOT_RES <- 100 # default is 72

const.FIGURE_PLOT_HEIGHT <- "60vh"
const.FIGURE_PLOT_HEIGHT_REDUCED <- "45vh"
const.FIGURE_PLOT_HEIGHT_LOGO <- "100px"

const.nb_of_digits <- 7

# convert degrees-per-second to radians-per-second
const.degps_2_radps <- 1 / 360 * 2 * pi

# constant default frequency for custom data
const.DEFAULT_FREQ <- 1 # [Hz]
const.DATASET_FREQ <- c(
  MTiG = 100,
  navchip = 250,
  imar = 400,
  ln200 = 400
)

# constants for the custom datasheet
const.DEFAULT_WN <- 5.304377e-07
const.DEFAULT_QN <- 1.681227e-06
const.DEFAULT_SIGMA2_GM <- 7.348293e-09
const.DEFAULT_BETA_GM <- 3.526037e-01
const.DEFAULT_RW <- 1.268314e-12
const.DEFAULT_DR <- 3.913529e-09
const.DEFAULT_BI <- const.DEFAULT_WN
const.DEFAULT_BIF0 <- NA

# source: https://www.xsens.com/wp-content/uploads/2013/11/MTi-G_User_Manual_and_Technical_Documentation.pdf
# https://www.xsens.com/tags/accelerometers/
# https://www.xsens.com/tags/gyroscopes/
# the frequency here is 100, because the dataset was acquired at this rate
const.MTIG.GYRO_WN <- (0.05 * const.degps_2_radps * sqrt(100))^2
const.MTIG.GYRO_BI <- (20 / 3600 * const.degps_2_radps)^2
const.MTIG.GYRO_BIF0 <- 1 / 40 # Hz
const.MTIG.ACC_WN <- (0.002 * sqrt(100))^2
const.MTIG.ACC_BI <- (30 * 1e-6 * 10)^2
const.MTIG.ACC_BIF0 <- 1 / 0.5 # Hz

# source: http://cdn-docs.av-iq.com/dataSheet//NavChip_Product_Brief.pdf,
# the frequency here is 250, because the dataset was acquired at this rate
const.NAVCHIP.GYRO_WN <- (0.003 * const.degps_2_radps * sqrt(250))^2 # [(rad/s)^2]
const.NAVCHIP.GYRO_BI <- (10 / 3600 * const.degps_2_radps)^2
const.NAVCHIP.GYRO_BIF0 <- 1 / 2 # Hz
const.NAVCHIP.ACC_WN <- (50 * 1e-6 * 10 * sqrt(250))^2 # [(m/s^2)^2]
const.NAVCHIP.ACC_BI <- (0.05 * 1e-3 * 10)^2
const.NAVCHIP.ACC_BIF0 <- 1 / 5 # Hz

# source: http://www.northropgrumman.com/Capabilities/LN200FOG/Documents/ln200.pdf
# the frequency here is 400, because the dataset was acquired at this rate
const.LN200.GYRO_WN <- (0.05 / 60 * sqrt(400) * const.degps_2_radps)^2
const.LN200.GYRO_BI <- NA
const.LN200.GYRO_BIF0 <- NA # Hz
const.LN200.ACC_WN <- const.DEFAULT_WN
const.LN200.ACC_BI <- NA
const.LN200.ACC_BIF0 <- NA # Hz

# source: http://www.imar-navigation.de/downloads/IMU_FSAS.pdf
# the frequency here is 400, because the dataset was acquired at this rate
const.IMAR.GYRO_WN <- (0.15 / 60 * sqrt(400) * const.degps_2_radps)^2
const.IMAR.GYRO_BI <- (0.1 / 3600 * const.degps_2_radps)^2
const.IMAR.GYRO_BIF0 <- 1 / 100 # Hz
const.IMAR.ACC_WN <- const.DEFAULT_WN
const.IMAR.ACC_BI <- NA
const.IMAR.ACC_BIF0 <- NA # Hz

################
# FUNCTIONS for COSINE INTEGRAL calculations
# Purpose: Helpers for bias instability / drift computations using cosine
#          integral terms.
################
cos_function <- function(t) {
  cos(t) / t
}

Ci <- function(x) {
  -integrate(f = cos_function, lower = x, upper = 2e3, subdivisions = 10000)$value
}

VCi <- Vectorize(Ci, c("x"))

sigma2_T <- function(T, f0, B) {
  2 * B * B / pi * (log(2) - ((sin(pi * f0 * T))^3) / (2 * (pi * f0 * T)^2) * (sin(pi * f0 * T) + 4 * pi * f0 * T * cos(pi * f0 * T)) + VCi(2 * pi * f0 * T) - VCi(4 * pi * f0 * T))
}



smac_url <- a("https://github.com/SMAC-Group/gui4gmwm2", href = "https://github.com/SMAC-Group/gui4gmwm2")
smac_url_description <- "gui4gmwm2 GitHub repository:"

# Increase file upload limit from default 5MB to 100MB
options(shiny.maxRequestSize = 100 * 1024^2)

################
# UI
# Purpose: Defines all layout, inputs, and outputs for the app.
################
ui <- shinyUI(fluidPage(
  theme = bs_theme(
    version = 5,
    bootswatch = "flatly",
    base_font = font_google("IBM Plex Sans"),
    heading_font = font_google("IBM Plex Sans Condensed")
  ),
  shinyjs::useShinyjs(),
  tags$style(HTML(".js-irs-0 .irs-single, .js-irs-0 .irs-bar-edge, .js-irs-0 .irs-bar {background: #111111}")),
  tags$style(HTML(".js-irs-1 .irs-single, .js-irs-1 .irs-bar-edge, .js-irs-1 .irs-bar {background: #333333}")),
  # App-wide CSS customizations
  tags$style(type = "text/css", "
    body {
      background: linear-gradient(180deg, #f7f9fb 0%, #eef2f6 100%);
    }
    .card .table thead th {
      background-color: #e2e8f0;
      color: #0f172a;
    }
    .github-link {
      display: inline-block;
      padding: 4px 8px;
      line-height: 1;
      background: #f8fafc;
      border: 1px solid #e2e8f0;
      border-radius: 10px;
      text-decoration: none;
      color: #0f172a;
    }
    .github-row {
      display: flex;
      align-items: center;
      justify-content: flex-start;
      gap: 2px;
      width: fit-content;
    }
    .github-row .shiny-image-output {
      width: auto;
      flex: 0 0 auto;
    }
    #render_logo_unige img,
    #tabhelpplotlogo_epfl img {
      width: auto !important;
      max-width: 100%;
      height: auto;
      display: block;
    }
    .affiliation-logos {
      display: flex;
      align-items: center;
      justify-content: space-between;
      gap: 12px;
      overflow: hidden;
    }
    .affiliation-logos .shiny-image-output {
      flex: 0 0 auto;
      display: flex;
      align-items: center;
      width: auto !important;
      height: auto !important;
    }
    .affiliation-logos img {
      height: 70px;
      width: auto;
      display: block;
    }
    .github-link img {
      height: 22px;
      width: auto;
      display: block;
    }
    .tab-content {
      padding-top: 10px;
    }
    .well, .card, .panel {
      border-radius: 14px;
      box-shadow: 0 8px 24px rgba(15, 23, 42, 0.08);
    }
    #summ {
      background-color: #ffffff;
      color: #0f172a;
      width: 100%;
      font-size: 14px;
      border-radius: 12px;
      padding: 12px;
    }
    .form-control, .selectize-input, .btn {
      border-radius: 10px;
    }
    .nav-tabs .nav-link.active,
    .nav-pills .nav-link.active {
      background-color: #111111 !important;
      border-color: #111111 !important;
      color: #ffffff !important;
    }
    .nav-tabs .nav-link {
      color: #111111;
    }
    .btn-primary {
      background-color: #111111;
      border-color: #111111;
    }
    .card .table {
      margin-top: 10px;
      width: 100%;
    }
    .card .table th, .card .table td {
      vertical-align: middle;
    }
  "),
  title = "GMWM GUI",
  # Top-level navigation tabs
  tabsetPanel(
    id = "tabs",
    # Wavelet variance tab: input selection + plot
    tabPanel("Wavelet Variance",
             div(style = "width:75%; margin: 0 auto;",
                 plotOutput(outputId = "plot_wv",
                            height = const.FIGURE_PLOT_HEIGHT,
                            width = "100%")),
             
             radioButtons("data_input_choice", "Select data input:", choices = c("From library" = "library", "Custom" = "custom")),
             # Sampling frequency (Hz); synced with Summary tab
             numericInput(
               inputId = "freq_input",
               label = "Sampling frequency (Hz):",
               min = 1,
               value = const.DATASET_FREQ[[names(data)[1]]],
               step = 1
             ),
             # Library dataset controls
             conditionalPanel(
               condition = "input.data_input_choice == 'library'",
               
               column(
                 3,
                 selectInput("imu_obj", "Select IMU file:",
                             names(data),
                             selected = 1
                 ),
                 selectInput("selected_sensor", "Select sensor", choices = "Gyro.X", selected = "Gyro.X"),
               )

             ),
             # Custom dataset controls
             conditionalPanel(
               condition = "input.data_input_choice == 'custom'",
               
                 column(
                   3,
                  
                   
          
                   radioButtons("user_specified_separator", "Specify the field separator character",
                                choices = c(Comma = ",",
                                            Semicolon = ";",
                                            Tab = "\t"),
                                selected = ","),
                   checkboxInput(inputId = "user_specified_header", label="The file contains the names of the variables as its first line (TRUE if checked).", value = FALSE, width = NULL),

                   

                   fileInput("user_defined_txt_file", "Select input file (max 100MB):",
                             accept = c(
                               "text/txt",
                               "text/comma-separated-values,text/plain",
                               ".txt",
                               ".imu",
                               ".csv",
                               placeholder = "No file selected"
                             )
                   ),
                   tags$small("Preview (first 5 rows) based on selected separator and header."),
                   tableOutput(outputId = "file_preview"),
                   span(
                     "",
                     div(style = "display:inline-block;",
                         title = "Specify the column that contains the signal in the provided file. The frequency is fixed to 1.",
                         icon("info-circle"))),
                   numericInput(inputId = "user_defined_txt_file_column", label = "Select column number:",
                                min = 1, max = 10, value = 1)

                 )

               
             ),
             

        

             ),
    # Model fit visualization tab
    tabPanel("GMWM fit",
             div(style = "width:75%; margin: 0 auto;",
                 plotOutput(outputId = "plot_fit", height = const.FIGURE_PLOT_HEIGHT, width = "100%"))
    ),

    # Summary tab: parameter estimates + KF-transformed parameters
    tabPanel(
      "Summary",
      fluidRow(
        column(
          6,
          bslib::card(
            bslib::card_header(h4("Estimated parameters")),
            tags$p(textOutput("summary_obj")),
            tableOutput(outputId = "summ")
          )
        ),
        column(
          6,
          bslib::card(
            bslib::card_header(h4("Kalman Filter estimated parameters")),
            numericInput(
              inputId = "freq_input_summary",
              label = "Sampling frequency (Hz):",
              min = 1,
              value = const.DATASET_FREQ[[names(data)[1]]],
              step = 1
            ),
            tableOutput(outputId = "summ_kf")
          )
        )
      )
    ),
    # Help & About tab
    tabPanel(
      "Help",
      fluidRow(
        column(
          7,
          bslib::card(
            bslib::card_header(h4("Help & About")),
            tags$p("This app provides tools to estimate composite stochastic models for IMU datatsets using the GMWM estimator."),
            # GitHub repo link with logo
            tags$div(
              class = "github-row",
              imageOutput(outputId = "logo_github", height = "60px"),
              tags$a(
                class = "github-link",
                href = "https://github.com/SMAC-Group/gui4gmwm2",
                target = "_blank",
                "SMAC-Group/gui4gmwm2"
              )
            ),
            tags$hr(),
            tags$p(strong("How to use")),
            tags$p("1. Select a dataset (library or custom)."),
            tags$p("2. Choose a sensor and model components."),
            tags$p("3. Click “Fit Model” to estimate and review results.")
          )
        ),
        column(
          5,
          bslib::card(
            bslib::card_header(h4("Credits")),
            tags$p(strong("Application developed by:")),
            tags$div(
              tags$a(href="https://stephaneguerrier.com/", "Stéphane Guerrier",  target="_blank"),
              tags$br(),
              tags$a(href="https://lionelvoirol.com/", "Lionel Voirol",  target="_blank"),
              tags$br(),
              "Philipp Clausen",
              tags$br(),
              "Justin Lee",
              tags$br(),
              tags$a(href="https://robertomolinari.github.io/", "Roberto Molinari",  target="_blank"),
              tags$br(),
              tags$a(href="https://people.epfl.ch/jan.skaloud", "Jan Skaloud",  target="_blank")
            )
          )
        )
      ),
      br(),
      div(style = "width:90%; margin: 0 auto;",
          bslib::card(
            bslib::card_header(h4("Affiliations")),
            tags$div(
              class = "affiliation-logos",
              imageOutput(outputId = "render_logo_unige"),
              imageOutput(outputId = "tabhelpplotlogo_epfl")
            )
          )
      )
    )
  ),
  
  
  conditionalPanel(condition = "input.tabs == 'Wavelet Variance' || input.tabs == 'GMWM fit'",
                   
                   
                   
                   column(
                     4,
                     checkboxGroupInput("model", "Select Model",
                                        c(
                                          "Quantization Noise" = "QN",
                                          "White Noise" = "WN",
                                          "Random Walk" = "RW",
                                          "Drift" = "DR",
                                          "Gauss-Markov" = "GM"
                                        ),
                                        selected = "WN"
                     ),
                     conditionalPanel(
                       condition = "input.model.indexOf('GM')>-1",
                       sliderInput("gm_nb", "Number of Gauss-Markov Processes", 1, 5, 2)
                     ),
                     actionButton("fit3", label = "Fit Model")
                   )
                   ),
  

))



################
# SERVER
# Purpose: Reactive logic, model fitting, plotting, and table rendering.
################
server <- function(input, output, session) {
  # === HELP TAB OUTPUTS ===
  # Render the GitHub repo text link (legacy output; kept for compatibility)
  output$tabhelpurl <- renderUI({
    tagList(smac_url_description, smac_url)
  })



  # Help tab: EPFL logo (uses local file in ./logo)
  output$tabhelpplotlogo_epfl <- renderImage(
    {
      # Prefer the specific EPFL variant; fall back if missing
      filename <- file.path("./logo", paste("logo_epfl_6", ".png", sep = ""))
 
      filename <- normalizePath(filename)
      list(src = filename, height = "80px")
    },
    deleteFile = FALSE
  )
  # Help tab: UNIGE logo (local file)
  output$render_logo_unige <- renderImage(
    {
      filename <- normalizePath(file.path("./logo", paste("logo_unige", ".png", sep = "")))
      list(src = filename, height = const.FIGURE_PLOT_HEIGHT_LOGO)
    },
    deleteFile = FALSE
  )

  # Help tab: GitHub logo (local file)
  output$logo_github <- renderImage(
    {
      filename <- normalizePath(file.path("./logo", paste("logo_github", ".png", sep = "")))
      list(src = filename, height = "60px")
    },
    deleteFile = FALSE
  )


  # === INPUT DERIVATIONS ===
  # Update sensor choices when a dataset changes
  observe({
    updateSelectInput(session, "selected_sensor", choices = names(data[[input$imu_obj]])[-length(names(data[[input$imu_obj]]))])
  })

  # === WAVELET VARIANCE PLOT ===
  # Render the empirical wavelet variance for either library or custom data
  output$plot_wv <- renderPlot({
    
    
    if ("library" %in% input$data_input_choice){ 
      # Using library data: plot selected sensor from selected IMU dataset
      # par(mar = c(4, 5, 3, 2))
      my_plot_wvar(data[[input$imu_obj]][[input$selected_sensor]], legend_position = "bottomleft")
      # plot()
    } else{ 
      # Using custom data: parse file, build WV, then plot
      inFile <- input$user_defined_txt_file
      if (is.null(inFile)){
        # No file selected yet: show placeholder plot
        plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
        text(.5,.5 , "Please provide a data file", cex=2)
      }else{
        # Robust file import with error handling (e.g., wrong separator)
        my_data <- tryCatch(
          read.csv(inFile$datapath, header = input$user_specified_header, sep = input$user_specified_separator),
          error = function(e) {
            plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
            text(0.5, 0.5, paste("File read error:", e$message), cex = 1)
            return(NULL)
          }
        )
        if (is.null(my_data)) {
          return()
        }

        # Extract and validate selected column
        raw_col <- my_data[, input$user_defined_txt_file_column]
        x <- suppressWarnings(as.numeric(raw_col))
        if (!any(is.finite(x))) {
          plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
          text(0.5, 0.5, "Selected column is not numeric.\nCheck separator, header, and column index.", cex = 1)
          return()
        }
        x <- x[is.finite(x)]
        if (length(x) < 2) {
          plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
          text(0.5, 0.5, "Not enough numeric data to compute wavelet variance.", cex = 1)
          return()
        }

        # Compute wavelet variance from selected column
        wv_obj <- tryCatch(
          wv::wvar(x),
          error = function(e) {
            plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
            text(0.5, 0.5, paste("Wavelet variance error:", e$message), cex = 1)
            return(NULL)
          }
        )
        if (is.null(wv_obj)) {
          return()
        }
        tryCatch(
          my_plot_wvar(wv_obj, legend_position = "bottomleft"),
          error = function(e) {
            plot(c(0, 1), c(0, 1), ann = FALSE, bty = "n", type = "n", xaxt = "n", yaxt = "n")
            text(0.5, 0.5, paste("Plot error:", e$message), cex = 1)
            return(NULL)
          }
        )
        }
      

    }

    
    
    

  })

  # === CUSTOM FILE PREVIEW ===
  # Show a small preview table to help users validate separator/header choices
  output$file_preview <- renderTable({
    if (!"custom" %in% input$data_input_choice) {
      return(NULL)
    }
    inFile <- input$user_defined_txt_file
    if (is.null(inFile)) {
      return(NULL)
    }
    df <- tryCatch(
      read.csv(inFile$datapath, header = input$user_specified_header, sep = input$user_specified_separator),
      error = function(e) {
        return(NULL)
      }
    )
    if (is.null(df)) {
      return(NULL)
    }
    head(df, 5)
  }, striped = FALSE, bordered = TRUE, spacing = "s")




  # === MODEL FITTING ===
  # Fit model (triggered by "Fit Model" action button)
  fit <- eventReactive(input$fit3, {
    withProgress(message = "Fitting desired model...", value = 0, {

    
      # Build model incrementally from selected components
      first <- TRUE
      counter_model_size <- 0

      if ("GM" %in% input$model) {
        for (i in 1:input$gm_nb) {
          counter_model_size <- counter_model_size + 1
          if (first == TRUE) {
            model <- GM()
            first <- FALSE
          } else {
            model <- model + GM()
          }
        }
      }

      if ("WN" %in% input$model) {
        counter_model_size <- counter_model_size + 1
        if (first == TRUE) {
          model <- WN()
          first <- FALSE
        } else {
          model <- model + WN()
        }
      }

      if ("QN" %in% input$model) {
        counter_model_size <- counter_model_size + 1
        if (first == TRUE) {
          model <- QN()
          first <- FALSE
        } else {
          model <- model + QN()
        }
      }


      if ("RW" %in% input$model) {
        counter_model_size <- counter_model_size + 1
        if (first == TRUE) {
          model <- RW()
          first <- FALSE
        } else {
          model <- model + RW()
        }
      }

      if ("DR" %in% input$model) {
        counter_model_size <- counter_model_size + 1
        if (first == TRUE) {
          model <- DR()
          first <- FALSE
        } else {
          model <- model + DR()
        }
      }

      # Safety fallback: if no model selected, use 3x GM
      if (is.null(model)) {
        model <- 3 * GM()
      }

      # Compute gmwm fit using either raw data (library) or WV object (custom)
      if ("library" %in% input$data_input_choice){ 
        gmwm::gmwm(model, data[[input$imu_obj]][[input$selected_sensor]], robust = F)
      }else{
        inFile <- input$user_defined_txt_file
        my_data <- tryCatch(
          read.csv(inFile$datapath, header = input$user_specified_header, sep = input$user_specified_separator),
          error = function(e) {
            showNotification(paste("File read error:", e$message), type = "error", duration = 6)
            return(NULL)
          }
        )
        if (is.null(my_data)) {
          return(NULL)
        }
        raw_col <- my_data[, input$user_defined_txt_file_column]
        x <- suppressWarnings(as.numeric(raw_col))
        if (!any(is.finite(x))) {
          showNotification("Selected column is not numeric. Check separator, header, and column index.", type = "error", duration = 6)
          return(NULL)
        }
        x <- x[is.finite(x)]
        if (length(x) < 2) {
          showNotification("Not enough numeric data to compute wavelet variance.", type = "error", duration = 6)
          return(NULL)
        }
        wv_obj <- tryCatch(
          wv::wvar(x),
          error = function(e) {
            showNotification(paste("Wavelet variance error:", e$message), type = "error", duration = 6)
            return(NULL)
          }
        )
        if (is.null(wv_obj)) {
          return(NULL)
        }
        gmwm::gmwm(model, wv_obj , robust = F)
      }
     

    })


  })
  
  # Navigate to "GMWM fit" tab after fitting
  observeEvent(input$fit3, {
    updateNavbarPage(session, "tabs", selected = "GMWM fit")
  })

  # Plot estimated fit
  output$plot_fit <- renderPlot({
    # plot(fit())
    my_plot_gmwm(fit())
  })
  
  
  # === FREQUENCY SYNC ===
  # Single source of truth for frequency across UI tabs
  freq_selected <- reactiveVal(const.DATASET_FREQ[[names(data)[1]]])

  # When user edits frequency on the Wavelet tab, sync to Summary
  observeEvent(input$freq_input, {
    if (!is.null(input$freq_input) && input$freq_input != freq_selected()) {
      freq_selected(input$freq_input)
      updateNumericInput(session, "freq_input_summary", value = input$freq_input)
    }
  }, ignoreInit = TRUE)

  # When user edits frequency on the Summary tab, sync to Wavelet
  observeEvent(input$freq_input_summary, {
    if (!is.null(input$freq_input_summary) && input$freq_input_summary != freq_selected()) {
      freq_selected(input$freq_input_summary)
      updateNumericInput(session, "freq_input", value = input$freq_input_summary)
    }
  }, ignoreInit = TRUE)

  # When a library dataset is selected, update frequency automatically
  observeEvent(input$imu_obj, {
    if (is.null(input$imu_obj)) {
      return()
    }
    if (!"library" %in% input$data_input_choice) {
      return()
    }
    freq <- const.DATASET_FREQ[[input$imu_obj]]
    if (is.null(freq)) {
      return()
    }
    if (freq != freq_selected()) {
      freq_selected(freq)
      updateNumericInput(session, "freq_input", value = freq)
      updateNumericInput(session, "freq_input_summary", value = freq)
    }
  }, ignoreInit = FALSE)

  # === SUMMARY OUTPUTS ===
  # Objective function summary line
  output$summary_obj <- renderText({
    gmwm_fit <- fit()
    paste("Objective Function:", formatC(gmwm_fit$obj.fun, digits = 4, format = "e"))
  })

  # Estimated parameters table (model, parameter, estimate)
  output$summ <- renderTable({
    gmwm_fit <- fit()
    est <- gmwm_fit$estimate
    # Normalize estimate object to a numeric vector
    if (is.matrix(est) || is.data.frame(est)) {
      if ("Estimate" %in% colnames(est)) {
        est <- est[, "Estimate"]
      }
      est <- as.numeric(est)
    }
    model_desc <- gmwm_fit$model$desc
    param_desc <- gmwm_fit$model$obj.desc
    rows <- list()
    idx <- 1
    for (i in seq_along(model_desc)) {
      model_name <- model_desc[[i]]
      params <- param_desc[[i]]
      # If parameter names are missing or numeric, inject canonical labels
      if (is.null(params) || length(params) == 0 || is.numeric(params)) {
        if (model_name == "WN") {
          params <- "SIGMA2"
        } else if (model_name == "RW") {
          params <- "GAMMA2"
        } else if (model_name == "GM") {
          params <- c("BETA", "SIGMA2_GM")
        } else if (model_name == "QN") {
          params <- "Q2"
        } else if (model_name == "DR") {
          params <- "OMEGA"
        } else {
          params <- paste0("param", seq_len(1))
        }
      }
      for (p in seq_along(params)) {
        # One row per parameter
        rows[[length(rows) + 1]] <- data.frame(
          Model = model_name,
          Parameter = as.character(params[[p]]),
          `Estimated parameters` = est[[idx]],
          stringsAsFactors = FALSE
        )
        idx <- idx + 1
      }
    }
    df <- do.call(rbind, rows)
    # Clean column names and format numeric values
    names(df) <- gsub("\\.", " ", names(df))
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], function(x) format(x, scientific = TRUE, digits = const.nb_of_digits))
    df
  }, striped = FALSE, bordered = TRUE, spacing = "s")

  # Kalman filter transformed parameters table (uses selected frequency)
  output$summ_kf <- renderTable({
    gmwm_fit <- fit()
    freq <- freq_selected()
    df <- transform_parameters(gmwm_fit, freq)

    # Format numeric values
    num_cols <- vapply(df, is.numeric, logical(1))
    df[num_cols] <- lapply(df[num_cols], function(x) format(x, scientific = TRUE, digits = const.nb_of_digits))
    df
  }, striped = FALSE, bordered = TRUE, spacing = "s")
      

}




# run app
shinyApp(ui, server)
