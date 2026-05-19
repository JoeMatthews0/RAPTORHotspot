library(shiny)
library(MASS)
library(dplyr)
library(ggplot2)
library(DT)
library(rjags)
library(patchwork)
library(leaflet)

# Build LA choices from STATS19 Data folder at startup
stats19_files   <- list.files("STATS19 Data", pattern = "_clusters\\.csv$")
stats19_la_keys <- sub("_clusters\\.csv$", "", stats19_files)
stats19_choices <- setNames(stats19_la_keys, gsub("_", " ", stats19_la_keys))

MODELSTRING <- "
  model {
    for (i in 1:n_past) {
      y[i] ~ dnegbin(1 / c[i], lambda[i] / (c[i] - 1))
      lambda[i] <- exp(log(mu[i]) + sigma + (alpha * year[i]))
      c[i] <- exp(-year[i] * tau)
    }
    y[n_past + 1] ~ dpois(lambda[n_past + 1])
    lambda[n_past + 1] <- exp(log(mu[n_past + 1]) + sigma)
    pred ~ dnegbin(1 / exp(tau), lambda_pred / (exp(tau) - 1))
    lambda_pred <- exp(log(predmu) + sigma + alpha)
    sigma ~ dnorm(0, 0.1)
    alpha <- alpha_n * alpha_z
    alpha_n ~ dnorm(0, 1)
    alpha_z ~ dbern(p)
    p ~ dunif(0, 1)
    tau ~ dgamma(2, 20)
  }
"

# ---- UI -------------------------------------------------------------------

ui <- navbarPage(
  title = "RAPTOR Hotspot",
  id    = "mainNav",
  header = tags$head(
    tags$style(HTML("
      .navbar-default { background-color: #2c3e50; border-color: #1a252f; }
      .navbar-default .navbar-brand { color: #fff; }
      .navbar-default .navbar-nav > li > a { color: #ccc; }
      .navbar-default .navbar-nav > .active > a,
      .navbar-default .navbar-nav > .active > a:hover,
      .navbar-default .navbar-nav > .active > a:focus {
        background-color: #e1202b !important; color: #fff !important; }
      .navbar-default .navbar-nav > li > a:hover { color: #fff; background-color: #3d5166; }
      body { font-size: 14px; }
      .status-ok   { color: #27ae60; font-weight: bold; }
      .status-warn { color: #e67e22; font-weight: bold; }
      .status-err  { color: #c0392b; font-weight: bold; }
      .well { border-radius: 4px; }
      pre { background: #f8f8f8; font-size: 12px; max-height: 300px; overflow-y: auto; }
      .site-checkbox-scroll {
        max-height: 320px; overflow-y: auto;
        border: 1px solid #ddd; border-radius: 3px;
        padding: 5px 8px; background: #fff;
      }
    "))
  ),

  # Introduction ------------------------------------------------------------
  tabPanel("Introduction",
    fluidPage(fluidRow(column(8, offset = 2,
      br(),
      h2("RAPTOR Hotspot Identification"),
      p(class = "lead",
        "Traffic collision hotspot prediction using Bayesian hierarchical modelling."),
      hr(),
      h4("Overview"),
      p("RAPTOR Hotspot analyses traffic collision data and provides predicted collision
        counts in a future time period. It fits a global Safety Performance Function (SPF)
        via a negative binomial GLM across all sites, then runs a Bayesian per-site MCMC
        model (using JAGS) to capture site-specific intercept offsets, temporal trends,
        and overdispersion."),
      h4("Workflow"),
      tags$ol(
        tags$li(strong("Data Upload — "), "Upload a CSV and map columns to site ID, year,
                and collision count. All remaining numeric columns become SPF predictors."),
        tags$li(strong("Site Selection — "), "Choose which sites to include in the analysis."),
        tags$li(strong("Simulation Settings — "), "Set MCMC parameters and click Run."),
        tags$li(strong("Results — "), "View predicted counts for the next year. Click any
                row to see that site's time-series with fitted Poisson mean and SPF."),
        tags$li(strong("Site Warnings — "), "Set a threshold count to get a colour-coded
                exceedance probability list.")
      ),
      h4("Data Format"),
      p("Tidy format: one row per site per year. Required columns:"),
      tags$ul(
        tags$li("Site ID — unique location identifier"),
        tags$li("Year — numeric year or time index"),
        tags$li("Collision count — non-negative integer")
      ),
      p("All other numeric columns are treated as SPF predictor variables.")
    )))
  ),

  # Data Upload -------------------------------------------------------------
  tabPanel("Data Upload",
    fluidPage(br(), fluidRow(
      column(4, wellPanel(
        h4("Data Source"),
        radioButtons("dataSource", label = NULL,
                     choices  = c("Use example data"   = "example",
                                  "Upload my own file" = "upload",
                                  "Use STATS-19 Data"  = "stats19"),
                     selected = "example"),
        conditionalPanel("input.dataSource === 'upload'",
          fileInput("dataFile", "Choose CSV File",
                    accept = c("text/csv", ".csv"), placeholder = "No file selected"),
          checkboxInput("hasHeader", "File has header row", value = TRUE),
          selectInput("sep", "Column separator:",
                      choices = c(Comma = ",", Semicolon = ";", Tab = "\t"), selected = ",")
        ),
        conditionalPanel("input.dataSource === 'stats19'",
          selectInput("stats19LA", "Select Local Authority:",
                      choices = stats19_choices)
        ),
        hr(),
        h4("Column Mapping"),
        selectInput("idCol",    "Site ID column:",       choices = NULL),
        selectInput("yearCol",  "Year column:",          choices = NULL),
        selectInput("countCol", "Collision count column:", choices = NULL),
        actionButton("validateBtn", "Validate Selection",
                     class = "btn-primary btn-block"),
        br(),
        uiOutput("validationStatus")
      )),
      column(8,
        h4("Data Preview"),
        p(class = "text-muted", "First 20 rows"),
        DTOutput("dataPreview"),
        uiOutput("clusterMapUI")
      )
    ))
  ),

  # Site Selection ----------------------------------------------------------
  tabPanel("Site Selection",
    fluidPage(br(), fluidRow(
      column(4, wellPanel(
        h4("Select Sites to Analyse"),
        fluidRow(
          column(6, actionButton("selectAllBtn", "Select All",
                                 class = "btn-info btn-sm btn-block")),
          column(6, actionButton("clearAllBtn",  "Clear All",
                                 class = "btn-warning btn-sm btn-block"))
        ),
        br(),
        uiOutput("siteSearchUI"),
        br(),
        uiOutput("siteCheckboxes")
      )),
      column(8,
        h4("Collision Trends"),
        p(class = "text-muted", "Observed counts over time for selected sites."),
        plotOutput("collisionTrends", height = "420px")
      )
    ))
  ),

  # Simulation Settings -----------------------------------------------------
  tabPanel("Simulation Settings",
    fluidPage(br(), fluidRow(
      column(4, wellPanel(
        h4("MCMC Settings"),
        numericInput("nAdapt",  "Adaptation iterations:",       value = 1000,  min = 100,  max = 10000, step = 100),
        numericInput("nBurnin", "Burn-in iterations:",          value = 2000,  min = 100,  max = 50000, step = 500),
        numericInput("nIter",   "Sampling iterations per site:", value = 10000, min = 1000, max = 100000, step = 1000),
        numericInput("nThin",   "Thinning factor:",              value = 5,     min = 1,    max = 50,    step = 1),
        hr(),
        actionButton("runBtn", "Run Model",
                     class = "btn-success btn-lg btn-block",
                     icon  = icon("play")),
        br(), br(),
        uiOutput("runStatus")
      )),
      column(8,
        h4("Pre-flight Checks"),
        uiOutput("preflightChecks"),
        br(),
        h4("Model Log"),
        verbatimTextOutput("mcmcLog")
      )
    ))
  ),

  # Results -----------------------------------------------------------------
  tabPanel("Results",
    fluidPage(br(),
      fluidRow(column(12,
        h4("Predicted Collision Counts — Next Year"),
        p(class = "text-muted", "Click any row to add that site's time-series plot below."),
        DTOutput("resultsTable"),
        br(),
        downloadButton("downloadResults", "Download table as CSV", class = "btn-sm btn-default")
      )),
      br(),
      fluidRow(column(12,
        div(style = "display:flex; align-items:center; gap:10px;",
          h4("Time-Series Plots", style = "margin:0;"),
          actionButton("clearPlotsBtn", "Clear all",
                       class = "btn-sm btn-danger", icon = icon("times")),
          downloadButton("exportAllPlots", "Export all (PDF)",
                         class = "btn-sm btn-default")
        ),
        br(),
        uiOutput("allPlotsUI")
      ))
    )
  ),

  # Site Warnings -----------------------------------------------------------
  tabPanel("Site Warnings",
    fluidPage(br(), fluidRow(
      column(4, wellPanel(
        h4("Warning Threshold"),
        numericInput("threshold", "Collision count threshold:", value = 10, min = 0, step = 1),
        helpText("Sites are colour-coded by the posterior probability that predicted
                 collisions in the next year exceed this threshold."),
        hr(),
        div(style = "padding:8px;margin-bottom:4px;background:#e74c3c;border-radius:3px;color:white;",
            icon("exclamation-triangle"), " High risk  — P > 0.5"),
        div(style = "padding:8px;margin-bottom:4px;background:#f39c12;border-radius:3px;",
            icon("exclamation"), " Medium risk — 0.2 < P \u2264 0.5"),
        div(style = "padding:8px;background:#27ae60;border-radius:3px;color:white;",
            icon("check"), " Low risk    — P \u2264 0.2")
      )),
      column(8,
        h4("Site Warning List"),
        p(class = "text-muted", "Sorted by exceedance probability (highest first)."),
        DTOutput("warningsTable"),
        uiOutput("warningsMapUI")
      )
    ))
  )
)

# ---- Plot helpers (defined outside server so they can be called from download handlers) ----

compute_ylim <- function(site_key, mcmc_results) {
  r <- mcmc_results[[site_key]]
  if (is.null(r) || !is.null(r$error)) return(NULL)
  upper <- max(
    max(r$site_df[[r$count_col]]),
    max(apply(r$lambda_samples, 2, quantile, 0.975)),
    quantile(r$pred_samples, 0.975)
  )
  c(0, upper * 1.05)
}

make_site_plot <- function(site_key, mcmc_results, ylim = NULL) {
  r <- mcmc_results[[site_key]]
  if (is.null(r) || !is.null(r$error)) return(NULL)

  sdf      <- r$site_df
  yrc      <- r$year_col
  cnc      <- r$count_col
  lam_s    <- r$lambda_samples
  mu_vec   <- r$mu_vec
  yr_vals  <- sdf[[yrc]]
  obs_vals <- sdf[[cnc]]

  lam_mean <- colMeans(lam_s)
  lam_lo   <- apply(lam_s, 2, quantile, 0.025)
  lam_hi   <- apply(lam_s, 2, quantile, 0.975)

  hist_df <- data.frame(year = yr_vals, observed = obs_vals,
                        lam_mean = lam_mean, lam_lo = lam_lo,
                        lam_hi   = lam_hi,   spf = mu_vec)

  future_yr <- max(yr_vals) + 1
  pred_df   <- data.frame(year      = future_yr,
                          pred_mean = mean(r$pred_samples),
                          pred_lo   = quantile(r$pred_samples, 0.025),
                          pred_hi   = quantile(r$pred_samples, 0.975))

  ggplot(hist_df, aes(x = year)) +
    geom_ribbon(aes(ymin = lam_lo, ymax = lam_hi), fill = "steelblue", alpha = 0.2) +
    geom_line(aes(y = lam_mean, colour = "Fitted Poisson mean"), linewidth = 1) +
    #geom_line(aes(y = spf, colour = "SPF (global model)"), linetype = "dashed", linewidth = 0.9) +
    geom_point(aes(y = observed), colour = "black", size = 3) +
    geom_errorbar(data = pred_df,
                  aes(x = year, ymin = pred_lo, ymax = pred_hi),
                  colour = "tomato", width = 0.25, linewidth = 1) +
    geom_point(data = pred_df, aes(x = year, y = pred_mean),
               colour = "tomato", size = 4, shape = 18) +
    scale_colour_manual(values = c("Fitted Poisson mean" = "steelblue"
    #                               ,"SPF (global model)"  = "darkorange"
    )) +
    labs(x = "Year", y = "Collision Count", colour = "",
         caption = paste0("Red diamond = predicted count for year ", future_yr,
                          " (95% PI). Ribbon = 95% CI on Poisson mean.")) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank()) +
    if (!is.null(ylim)) coord_cartesian(ylim = ylim) else NULL
}

make_pred_dist_plot <- function(site_key, mcmc_results, ylim = NULL) {
  r <- mcmc_results[[site_key]]
  if (is.null(r) || !is.null(r$error)) return(NULL)

  ps        <- as.integer(round(r$pred_samples))
  pred_mean <- mean(r$pred_samples)
  pred_lo   <- quantile(r$pred_samples, 0.025)
  pred_hi   <- quantile(r$pred_samples, 0.975)

  tab    <- table(ps)
  pmf_df <- data.frame(count = as.integer(names(tab)),
                       prob  = as.numeric(tab) / length(ps))

  ggplot(pmf_df, aes(x = count, y = prob)) +
    geom_col(fill = "tomato", alpha = 0.75, width = 0.8) +
    geom_vline(xintercept = pred_mean, colour = "tomato", linewidth = 0.8,
               linetype = "dashed") +
    geom_vline(xintercept = pred_lo,   colour = "tomato", linewidth = 0.5,
               linetype = "dotted") +
    geom_vline(xintercept = pred_hi,   colour = "tomato", linewidth = 0.5,
               linetype = "dotted") +
    coord_flip(xlim = ylim) +
    labs(x = "Collision count", y = "Posterior probability",
         caption = paste0("Mean = ", round(pred_mean, 1),
                          "  \u2022  95% PI: [", round(pred_lo, 1),
                          ", ", round(pred_hi, 1), "]")) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
}

# Combine both plots into one aligned patchwork object.
# guides = "collect" pulls the TS legend out of the panel area so both panels
# are sized identically; plot_annotation adds the shared site title.
combine_site_plots <- function(p1, p2, site_key) {
  p1 + p2 +
    plot_layout(widths = c(7, 5), guides = "collect") +
    plot_annotation(title = paste("Site:", site_key)) &
    theme(legend.position = "bottom")
}

# ---- Server ---------------------------------------------------------------

server <- function(input, output, session) {

  rv <- reactiveValues(
    rawData     = NULL,
    glmModel    = NULL,
    mcmcResults = NULL,
    mcmcLog     = ""
  )

  # Data Upload -----------------------------------------------------------

  # Shared helper: populate column-mapping dropdowns from a data frame
  set_col_defaults <- function(df) {
    cols  <- colnames(df)
    id_g  <- grep("^id$|site|location|cluster_id", cols, ignore.case = TRUE, value = TRUE)
    yr_g  <- grep("^year$|time|period",  cols, ignore.case = TRUE, value = TRUE)
    cnt_g <- grep("accident|collision|crash|count", cols, ignore.case = TRUE, value = TRUE)
    updateSelectInput(session, "idCol",    choices = cols,
                      selected = if (length(id_g))  id_g[1]  else cols[1])
    updateSelectInput(session, "yearCol",  choices = cols,
                      selected = if (length(yr_g))  yr_g[1]  else cols[min(2, length(cols))])
    updateSelectInput(session, "countCol", choices = cols,
                      selected = if (length(cnt_g)) cnt_g[1] else cols[min(3, length(cols))])
  }

  load_stats19 <- function(la_key) {
    tryCatch({
      path <- file.path("STATS19 Data", paste0(la_key, "_clusters.csv"))
      df   <- read.csv(path, stringsAsFactors = FALSE)
      rv$rawData <- df
      set_col_defaults(df)
    }, error = function(e) {
      showNotification(paste("Could not load STATS-19 data:", e$message),
                       type = "error", duration = 8)
    })
  }

  # Load example data when selected
  observeEvent(input$dataSource, {
    if (input$dataSource == "example") {
      tryCatch({
        df <- read.csv("ExampleData.csv", stringsAsFactors = FALSE)
        rv$rawData <- df
        set_col_defaults(df)
      }, error = function(e) {
        showNotification(paste("Could not load example data:", e$message),
                         type = "error", duration = 8)
      })
    } else if (input$dataSource == "stats19") {
      req(input$stats19LA)
      load_stats19(input$stats19LA)
    } else {
      rv$rawData <- NULL
    }
  }, ignoreInit = FALSE)

  # Reload when user picks a different LA
  observeEvent(input$stats19LA, {
    req(input$dataSource == "stats19")
    load_stats19(input$stats19LA)
  }, ignoreInit = TRUE)

  # Load uploaded file
  observeEvent(input$dataFile, {
    req(input$dataFile, input$dataSource == "upload")
    tryCatch({
      df <- read.csv(input$dataFile$datapath,
                     header           = input$hasHeader,
                     sep              = input$sep,
                     stringsAsFactors = FALSE)
      rv$rawData <- df
      set_col_defaults(df)
    }, error = function(e) {
      showNotification(paste("Error reading file:", e$message), type = "error", duration = 8)
    })
  })

  output$dataPreview <- renderDT({
    req(rv$rawData)
    datatable(rv$rawData,
              options = list(scrollX = TRUE, pageLength = 10, dom = "tip"),
              rownames = FALSE)
  })

  # Validation ------------------------------------------------------------

  validResult <- reactiveVal(NULL)

  observeEvent(input$validateBtn, {
    req(rv$rawData, input$idCol, input$yearCol, input$countCol)
    df  <- rv$rawData
    idc <- input$idCol; yrc <- input$yearCol; cnc <- input$countCol
    errs  <- character(0)
    warns <- character(0)

    if (length(unique(c(idc, yrc, cnc))) < 3)
      errs <- c(errs, "Site ID, Year, and Count columns must all be different.")
    if (!is.numeric(df[[yrc]]))
      errs <- c(errs, paste0("Year column '", yrc, "' must be numeric."))
    if (!is.numeric(df[[cnc]]))
      errs <- c(errs, paste0("Count column '", cnc, "' must be numeric."))
    if (length(errs) == 0 && any(df[[cnc]] < 0, na.rm = TRUE))
      errs <- c(errs, "Collision counts cannot be negative.")

    min_per_site <- min(tapply(df[[yrc]], df[[idc]], function(x) length(unique(x))))
    if (length(errs) == 0 && min_per_site < 2)
      warns <- c(warns, "Some sites have only 1 time point and will be skipped.")

    pred_num <- setdiff(colnames(df), c(idc, cnc))
    pred_num <- pred_num[sapply(pred_num, function(co) is.numeric(df[[co]]))]
    if (length(pred_num) < 1)
      warns <- c(warns, "No numeric predictor columns found (besides year). SPF may be weak.")

    if (length(errs) == 0) {
      validResult(list(ok = TRUE, warns = warns,
                       n_sites = length(unique(df[[idc]])),
                       n_years = length(unique(df[[yrc]]))))
    } else {
      validResult(list(ok = FALSE, errors = errs))
    }
  })

  output$validationStatus <- renderUI({
    res <- validResult()
    if (is.null(res)) return(NULL)
    if (res$ok) {
      tagList(
        p(class = "status-ok", icon("check-circle"), " Validation passed"),
        p(paste0("Sites: ", res$n_sites, "  |  Time points: ", res$n_years)),
        lapply(res$warns, function(w) p(class = "status-warn", icon("warning"), " ", w))
      )
    } else {
      tagList(
        p(class = "status-err", icon("times-circle"), " Validation failed"),
        tags$ul(lapply(res$errors, tags$li))
      )
    }
  })

  # Site Selection --------------------------------------------------------

  availSites <- reactive({
    req(rv$rawData, input$idCol)
    sort(unique(rv$rawData[[input$idCol]]))
  })

  output$siteSearchUI <- renderUI({
    selectizeInput("siteSearch", "Search / type site IDs:",
                   choices  = availSites(),
                   selected = availSites(),
                   multiple = TRUE,
                   options  = list(placeholder = "Type to filter sites...",
                                   plugins     = list("remove_button")))
  })

  output$siteCheckboxes <- renderUI({
    div(class = "site-checkbox-scroll",
      checkboxGroupInput("selectedSites", label = NULL,
                         choices = availSites(), selected = availSites())
    )
  })

  # Keep the two selection inputs in sync using a flag to prevent loops
  syncingSites <- reactiveVal(FALSE)

  observeEvent(input$siteSearch, {
    if (isTRUE(syncingSites())) return()
    syncingSites(TRUE)
    updateCheckboxGroupInput(session, "selectedSites",
                             selected = input$siteSearch %||% character(0))
    syncingSites(FALSE)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$selectedSites, {
    if (isTRUE(syncingSites())) return()
    syncingSites(TRUE)
    updateSelectizeInput(session, "siteSearch",
                         selected = input$selectedSites %||% character(0))
    syncingSites(FALSE)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

  observeEvent(input$selectAllBtn, {
    updateCheckboxGroupInput(session, "selectedSites", selected = availSites())
    updateSelectizeInput(session, "siteSearch", selected = availSites())
  })
  observeEvent(input$clearAllBtn, {
    updateCheckboxGroupInput(session, "selectedSites", selected = character(0))
    updateSelectizeInput(session, "siteSearch", selected = character(0))
  })

  output$collisionTrends <- renderPlot({
    req(rv$rawData, input$idCol, input$yearCol, input$countCol)
    sel <- input$selectedSites
    validate(need(length(sel) > 0, "No sites selected."))
    df <- rv$rawData[rv$rawData[[input$idCol]] %in% sel, ]
    df$site_f <- as.factor(df[[input$idCol]])
    ggplot(df, aes(x = .data[[input$yearCol]], y = .data[[input$countCol]],
                   group = site_f, colour = site_f)) +
      geom_line(alpha = 0.6) +
      geom_point(size = 1.8) +
      labs(x = "Year", y = "Collision Count", colour = "Site") +
      theme_minimal(base_size = 13) +
      theme(legend.position = if (length(sel) > 15) "none" else "right",
            panel.grid.minor = element_blank())
  })

  # Pre-flight checks -----------------------------------------------------

  output$preflightChecks <- renderUI({
    c1 <- if (!is.null(rv$rawData))
      p(class = "status-ok",  icon("check"), " Data uploaded (", nrow(rv$rawData), " rows)")
    else
      p(class = "status-err", icon("times"), " No data uploaded")

    val <- validResult()
    c2 <- if (!is.null(val) && val$ok)
      p(class = "status-ok",  icon("check"), " Column mapping validated")
    else
      p(class = "status-warn", icon("warning"), " Columns not yet validated (Data Upload tab)")

    n <- length(input$selectedSites)
    c3 <- if (n > 0)
      p(class = "status-ok",  icon("check"), " ", n, " site(s) selected")
    else
      p(class = "status-err", icon("times"), " No sites selected")

    tagList(c1, c2, c3)
  })

  # Run Model -------------------------------------------------------------

  observeEvent(input$runBtn, {
    req(rv$rawData, input$idCol, input$yearCol, input$countCol)

    val <- validResult()
    if (is.null(val) || !val$ok) {
      showNotification("Validate column selection first (Data Upload tab).", type = "warning")
      return()
    }
    if (length(input$selectedSites) == 0) {
      showNotification("Select at least one site.", type = "warning")
      return()
    }

    df      <- rv$rawData
    idc     <- input$idCol
    yrc     <- input$yearCol
    cnc     <- input$countCol
    sites   <- input$selectedSites
    n_adapt  <- input$nAdapt
    n_burn   <- input$nBurnin
    n_iter   <- input$nIter
    n_thin   <- input$nThin

    # Numeric predictors (exclude ID and count; include year)
    preds <- setdiff(colnames(df), c(idc, cnc))
    preds <- preds[sapply(preds, function(co) is.numeric(df[[co]]))]

    if (length(preds) == 0) {
      showNotification("No numeric predictor columns found.", type = "error")
      return()
    }

    glm_form <- as.formula(paste(cnc, "~", paste(preds, collapse = " + ")))

    msgs <- character(0)
    append_log <- function(msg) {
      msgs <<- c(msgs, paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", msg))
      rv$mcmcLog <- paste(msgs, collapse = "\n")
    }

    withProgress(message = "RAPTOR — Running Analysis", value = 0, {

      # Fit global SPF
      incProgress(0.05, detail = "Fitting global SPF...")
      append_log(paste("Fitting SPF:", deparse(glm_form)))

      glm_fit <- tryCatch(
        glm.nb(glm_form, data = df),
        error = function(e) {
          showNotification(paste("SPF fitting failed:", e$message), type = "error", duration = 10)
          NULL
        }
      )
      if (is.null(glm_fit)) return()
      rv$glmModel <- glm_fit
      append_log(paste0("SPF fitted. AIC = ", round(AIC(glm_fit), 1)))

      df$mu_spf <- fitted(glm_fit)
      results   <- list()
      n_sites   <- length(sites)

      for (i in seq_along(sites)) {
        site <- sites[i]
        incProgress(0.9 / n_sites,
                    detail = paste0("Site ", site, " (", i, "/", n_sites, ")"))

        sdf <- df[df[[idc]] == site, ]
        sdf <- sdf[order(sdf[[yrc]]), ]
        n   <- nrow(sdf)

        if (n < 2) {
          results[[as.character(site)]] <- list(error = "< 2 observations — skipped")
          append_log(paste("SKIP", site, "— fewer than 2 obs"))
          next
        }

        y_vec    <- sdf[[cnc]]
        year_vec <- sdf[[yrc]] - max(sdf[[yrc]])
        mu_vec   <- pmax(sdf$mu_spf, 1e-6)

        predmu   <- max(mu_vec[n], 1e-6)

        jd <- list(y = y_vec, year = year_vec, mu = mu_vec,
                   predmu = predmu, n_past = n - 1)

        res <- tryCatch({
          jm <- jags.model(textConnection(MODELSTRING), data = jd,
                           n.adapt = n_adapt, quiet = TRUE)
          update(jm, n.iter = n_burn)
          samps <- coda.samples(jm, variable.names = c("lambda", "pred"),
                                n.iter = n_iter, thin = n_thin)
          smat  <- as.matrix(samps)
          lnames <- paste0("lambda[", seq_len(n), "]")
          list(pred_samples   = smat[, "pred"],
               lambda_samples = smat[, lnames, drop = FALSE],
               site_df        = sdf,
               year_vec       = year_vec,
               mu_vec         = mu_vec,
               predmu         = predmu,
               year_col       = yrc,
               count_col      = cnc,
               error          = NULL)
        }, error = function(e) list(error = paste("JAGS:", conditionMessage(e))))

        results[[as.character(site)]] <- res
        if (!is.null(res$error)) {
          append_log(paste("ERROR", site, "—", res$error))
        } else {
          append_log(paste0("OK ", site, " — mean pred = ",
                            round(mean(res$pred_samples), 2)))
        }
      }

      rv$mcmcResults <- results
      incProgress(0.05, detail = "Done.")
    })

    n_ok  <- sum(sapply(rv$mcmcResults, function(x) is.null(x$error)))
    n_err <- length(rv$mcmcResults) - n_ok
    showNotification(
      paste0("Complete. ", n_ok, " site(s) processed",
             if (n_err > 0) paste0(", ", n_err, " with errors.") else "."),
      type = if (n_err == 0) "message" else "warning", duration = 6)
  })

  output$mcmcLog  <- renderText({ rv$mcmcLog })

  output$runStatus <- renderUI({
    req(rv$mcmcResults)
    n_ok  <- sum(sapply(rv$mcmcResults, function(x) is.null(x$error)))
    n_err <- length(rv$mcmcResults) - n_ok
    tagList(
      if (n_ok  > 0) p(class = "status-ok",   icon("check"),   " ", n_ok,  " site(s) processed"),
      if (n_err > 0) p(class = "status-warn",  icon("warning"), " ", n_err, " site(s) with errors")
    )
  })

  # Results Table ---------------------------------------------------------

  resultsDF <- reactive({
    req(rv$mcmcResults)
    do.call(rbind, lapply(names(rv$mcmcResults), function(s) {
      r <- rv$mcmcResults[[s]]
      if (!is.null(r$error)) {
        data.frame(Site = s, Mean = NA_real_, Median = NA_real_,
                   Lower95 = NA_real_, Upper95 = NA_real_,
                   Status = r$error, stringsAsFactors = FALSE)
      } else {
        ps <- r$pred_samples
        data.frame(Site    = s,
                   Mean    = round(mean(ps),                 2),
                   Median  = round(median(ps),               2),
                   Lower95 = round(quantile(ps, 0.025),      2),
                   Upper95 = round(quantile(ps, 0.975),      2),
                   Status  = "OK",
                   stringsAsFactors = FALSE)
      }
    }))
  })

  output$resultsTable <- renderDT({
    req(resultsDF())
    datatable(resultsDF(), selection = "single", rownames = FALSE,
              options  = list(pageLength = 15, scrollX = TRUE, dom = "tip"),
              colnames = c("Site ID", "Mean", "Median", "Lower 95% CI", "Upper 95% CI", "Status"))
  })

  output$downloadResults <- downloadHandler(
    filename = function() paste0("RAPTOR_results_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(resultsDF(), file, row.names = FALSE)
  )

  # Time Series Plots (accumulating) --------------------------------------

  plotSites <- reactiveVal(character(0))

  observeEvent(input$resultsTable_rows_selected, {
    idx <- input$resultsTable_rows_selected
    if (length(idx) > 0) {
      site <- resultsDF()$Site[idx]
      if (!site %in% plotSites()) plotSites(c(plotSites(), site))
    }
  })

  observeEvent(input$clearPlotsBtn, { plotSites(character(0)) })

  output$allPlotsUI <- renderUI({
    sites <- plotSites()
    if (length(sites) == 0)
      return(p(class = "text-muted", "Click a row in the table above to add a plot here."))

    tagList(lapply(sites, function(site) {
      sid          <- make.names(site)
      combined_id  <- paste0("combined_", sid)
      dl_id        <- paste0("dl_plot_",  sid)
      rm_id        <- paste0("rm_plot_",  sid)
      div(style = "margin-bottom: 30px;",
        hr(),
        div(style = "display:flex; justify-content:space-between; align-items:center; margin-bottom:6px;",
          tags$span(),
          div(
            downloadButton(dl_id, "Export PNG", class = "btn-xs btn-default"),
            tags$span(" "),
            actionButton(rm_id, "", icon = icon("times"),
                         class = "btn-xs btn-danger", title = "Remove plot")
          )
        ),
        plotOutput(combined_id, height = "380px")
      )
    }))
  })

  # Dynamically create render + download + remove handlers for each plot.
  # Use session$userData to ensure each observer is registered only once.
  observe({
    req(rv$mcmcResults)
    sites <- plotSites()
    if (is.null(session$userData$plot_obs)) session$userData$plot_obs <- character(0)

    lapply(sites, function(site) {
      local({
        s           <- site
        sid         <- make.names(s)
        combined_id <- paste0("combined_", sid)
        dl_id       <- paste0("dl_plot_",  sid)
        rm_id       <- paste0("rm_plot_",  sid)

        output[[combined_id]] <- renderPlot({
          ylim <- compute_ylim(s, rv$mcmcResults)
          p1 <- make_site_plot(s, rv$mcmcResults, ylim = ylim)
          p2 <- make_pred_dist_plot(s, rv$mcmcResults, ylim = ylim)
          validate(need(!is.null(p1), paste("No results available for site", s)))
          combine_site_plots(p1, p2, s)
        })

        output[[dl_id]] <- downloadHandler(
          filename = function() paste0("RAPTOR_site_", s, "_", Sys.Date(), ".png"),
          content  = function(file) {
            ylim <- compute_ylim(s, rv$mcmcResults)
            p1 <- make_site_plot(s, rv$mcmcResults, ylim = ylim)
            p2 <- make_pred_dist_plot(s, rv$mcmcResults, ylim = ylim)
            ggsave(file, plot = combine_site_plots(p1, p2, s),
                   width = 15, height = 5, dpi = 150, device = "png")
          }
        )

        # Register remove-button observer only once per site
        if (!rm_id %in% session$userData$plot_obs) {
          session$userData$plot_obs <- c(session$userData$plot_obs, rm_id)
          observeEvent(input[[rm_id]], {
            plotSites(setdiff(plotSites(), s))
          }, ignoreInit = TRUE)
        }
      })
    })
  })

  output$exportAllPlots <- downloadHandler(
    filename = function() paste0("RAPTOR_plots_", Sys.Date(), ".pdf"),
    content  = function(file) {
      req(rv$mcmcResults)
      sites <- plotSites()
      validate(need(length(sites) > 0, "No plots to export."))
      pdf(file, width = 15, height = 5)
      for (s in sites) {
        ylim <- compute_ylim(s, rv$mcmcResults)
        p1 <- make_site_plot(s, rv$mcmcResults, ylim = ylim)
        p2 <- make_pred_dist_plot(s, rv$mcmcResults, ylim = ylim)
        if (!is.null(p1) && !is.null(p2))
          print(combine_site_plots(p1, p2, s))
      }
      dev.off()
    }
  )

  # Cluster positions helper -----------------------------------------------

  clusterPositions <- reactive({
    req(rv$rawData)
    df <- rv$rawData
    if (!all(c("centroid_lon", "centroid_lat") %in% colnames(df))) return(NULL)
    req(input$idCol %in% colnames(df))
    df[!duplicated(df[[input$idCol]]), c(input$idCol, "centroid_lon", "centroid_lat")]
  })

  # Cluster map on Data Upload tab ----------------------------------------

  output$clusterMapUI <- renderUI({
    req(clusterPositions())
    tagList(br(), h4("Cluster Locations"), leafletOutput("clusterMap", height = "420px"))
  })

  output$clusterMap <- renderLeaflet({
    pos    <- clusterPositions()
    req(pos)
    id_col <- input$idCol
    leaflet(pos) %>%
      addTiles() %>%
      addCircleMarkers(
        lng        = ~centroid_lon,
        lat        = ~centroid_lat,
        popup      = ~as.character(pos[[id_col]]),
        label      = ~as.character(pos[[id_col]]),
        radius     = 6,
        color      = "#2c3e50",
        fillColor  = "#3d5166",
        fillOpacity = 0.7,
        weight     = 1
      )
  })

  # Site Warnings ---------------------------------------------------------

  warnDF <- reactive({
    req(rv$mcmcResults, input$threshold)
    thr <- as.numeric(input$threshold)
    df  <- do.call(rbind, lapply(names(rv$mcmcResults), function(s) {
      r <- rv$mcmcResults[[s]]
      if (!is.null(r$error)) {
        data.frame(Site = s, MeanPred = NA_real_, ProbExceed = NA_real_,
                   Risk = "Error", stringsAsFactors = FALSE)
      } else {
        ps <- r$pred_samples
        pe <- mean(ps > thr)
        data.frame(Site       = s,
                   MeanPred   = round(mean(ps), 2),
                   ProbExceed = round(pe, 3),
                   Risk       = ifelse(pe > 0.5, "High",
                                       ifelse(pe > 0.2, "Medium", "Low")),
                   stringsAsFactors = FALSE)
      }
    }))
    df[order(-df$ProbExceed, na.last = TRUE), ]
  })

  output$warningsTable <- renderDT({
    warn_df <- warnDF()
    thr     <- as.numeric(input$threshold)
    datatable(warn_df, rownames = FALSE, selection = "single",
              options  = list(pageLength = 15, scrollX = TRUE, dom = "tip"),
              colnames = c("Site ID", "Mean Predicted",
                           paste0("P(count > ", thr, ")"), "Risk Level")) |>
      formatStyle("Risk",
                  backgroundColor = styleEqual(
                    c("High",    "Medium",  "Low",     "Error"),
                    c("#e74c3c", "#f39c12", "#27ae60", "#dddddd")),
                  color = styleEqual(
                    c("High",   "Low"),
                    c("white",  "white"), default = "black"))
  })

  output$warningsMapUI <- renderUI({
    req(warnDF(), clusterPositions())
    tagList(br(), h4("Risk Map"), leafletOutput("warningsMap", height = "420px"))
  })

  output$warningsMap <- renderLeaflet({
    warn_df <- warnDF()
    pos     <- clusterPositions()
    req(warn_df, pos)
    id_col  <- input$idCol
    merged  <- merge(warn_df, pos, by.x = "Site", by.y = id_col)

    risk_colors <- c("High" = "#e74c3c", "Medium" = "#f39c12",
                     "Low"  = "#27ae60", "Error"  = "#aaaaaa")
    merged$color <- risk_colors[merged$Risk]

    leaflet(merged) %>%
      addTiles() %>%
      addCircleMarkers(
        lng         = ~centroid_lon,
        lat         = ~centroid_lat,
        color       = ~color,
        fillColor   = ~color,
        fillOpacity = 0.8,
        weight      = 1,
        radius      = 8,
        popup       = ~paste0("<b>Site: ", Site, "</b><br>Risk: ", Risk,
                              "<br>Mean Predicted: ", MeanPred,
                              "<br>P(exceed): ", ProbExceed),
        label       = ~as.character(Site)
      ) %>%
      addLegend("bottomright",
                colors  = c("#e74c3c", "#f39c12", "#27ae60"),
                labels  = c("High", "Medium", "Low"),
                title   = "Risk Level",
                opacity = 0.8)
  })

  # Highlight selected row on the risk map --------------------------------

  observeEvent(input$warningsTable_rows_selected, {
    idx     <- input$warningsTable_rows_selected
    proxy   <- leafletProxy("warningsMap")
    proxy %>% clearGroup("highlight")

    if (length(idx) > 0) {
      warn_df <- warnDF()
      pos     <- clusterPositions()
      req(warn_df, pos)
      id_col       <- input$idCol
      selected_site <- warn_df$Site[idx]
      site_pos      <- pos[pos[[id_col]] == selected_site, ]

      proxy %>%
        addCircleMarkers(
          data        = site_pos,
          lng         = ~centroid_lon,
          lat         = ~centroid_lat,
          group       = "highlight",
          radius      = 14,
          color       = "white",
          fillColor   = "transparent",
          fillOpacity = 0,
          weight      = 3,
          opacity     = 1
        ) %>%
        flyTo(lng = site_pos$centroid_lon, lat = site_pos$centroid_lat, zoom = 13)
    }
  }, ignoreNULL = FALSE)
}

shinyApp(ui = ui, server = server)
