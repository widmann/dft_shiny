# dft_shiny.R - A Shiny App to illustrate DFT and demonstrate potential effects of spectral filtering
# Copyright (c) 2025 Andreas Widmann, University of Leipzig
# Author: Andreas Widmann, widmann@uni-leipzig.de

library(shiny)
library(bslib)
library(ggplot2)
library(tidyr)
library(dplyr)

# Potentially interesting signals:
# gausswin(20,6)'
# 0.0000, 0.0000, 0.0000, 0.0002, 0.0024, 0.0176, 0.0869, 0.2875, 0.6384, 0.9514, 0.9514, 0.6384, 0.2875, 0.0869, 0.0176, 0.0024, 0.0002, 0.0000, 0.0000, 0.0000
# fir_filterdcpadded(firws(50, 0.5, 'high'), 1, randn(1,20)')'
# 0.3851, -0.3447, -0.3281, 1.0845, -1.1488, 0.4373, 0.3982, -0.6941, 0.4223, -0.0680, -0.0273, 0.0881, -0.5255, 1.1685, -1.2634, 0.4673, 0.6034, -1.0597, 0.7464, -0.2347
# gausswin(20,6)' + fir_filterdcpadded(firws(50, 0.5, 'high'), 1, randn(1,20)')'
# 0.5117, -0.5396, 0.4371, -0.5070, 0.5001, -0.2053, 0.2613, -0.2892, 1.2650, 1.1429, 0.1581, 0.6580, 1.5769, -1.3640, 0.6151, -0.3489, 1.1182, -1.5572, 0.8954, -0.0129
# 1,1,1,1,1,1,1,1,1,1,1

# https://jfly.uni-koeln.de/color/
okabe<- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# UI
ui <- page_sidebar(
  
  tags$head(
    tags$style("label {font-size:80%;}")
  ),
  
  sidebar = sidebar(
    "Input data",
    textInput('fs', 'Enter sampling frequency', 100), tags$style(type='text/css', "#fs {font-size:80%;}"),
    textInput('vec1', 'Enter a data vector (comma delimited)', "0,0,0,0,0,1,1,1,1,1"), tags$style(type='text/css', "#vec1 {font-size:80%;}"),
    textInput('vec2', 'Enter a frequency weights vector (comma delimited)', "1,1,1,1,1,1"), tags$style(type='text/css', "#vec2 {font-size:80%;}")
  ),
  
  layout_columns(
    card(card_header("Time domain: Data"),  plotOutput("plot_td")),
    card(card_header("Frequency domain: Weights"),  plotOutput("plot_weights")),
    card(card_header("Frequency domain"), plotOutput("plot_fd")),
    card(card_header("Frequency domain: (Co-)sines"), plotOutput("plot_coeff")),
    card(card_header("Time domain: Synthesized data"), plotOutput("plot_synth")),
    card(card_header("Time domain: Raw - filtered data"), plotOutput("plot_diff")),

    col_widths = c(6, 6, 4, 8, 6, 6),
    row_heights = c(1.5, 3, 1.5)
  )
  
)

# Server
server <- function(input, output) {
  
  # Prepare data
  val <- reactive({
    values <- list(fs = as.numeric(input$fs), data = as.numeric(unlist(strsplit(input$vec1,","))))
    values[["n"]] <- length(values$data)
    values[["time"]] <- seq(0, values$n - 1) / values$fs
    
    f <- seq(0, values$n - 1) * values$fs / values$n
    values$f <- f[f <= values$fs / 2]
    
    weights <- rep(1, length(values$f))
    tmp <- as.numeric(unlist(strsplit(input$vec2,",")))
    tmp[tmp > 1] <- 1
    tmp[tmp < 0] <- 0
    weights[1:min(length(tmp), length(weights))] <- tmp[1:min(length(tmp), length(weights))]
    values$weights <- weights
    
    z <- fft(values$data) / values$n
    z <- z[f <= values$fs / 2]
    z[values$f > 0 & values$f < values$fs / 2] <- z[values$f > 0 & values$f < values$fs / 2] * 2
    values$z <- z
    
    z_weighted <- z * weights
    values$z_weighted <- z_weighted
    
    # print(values)
    values
  })
  
  # Prepare (co-)sines
  dat_coeff <- reactive({
    time <- val()$time
    tf_mat <- sweep(matrix(rep(time, length(val()$f)), nrow = length(time), ncol = length(val()$f), byrow = FALSE), MARGIN = 2, val()$f, `*`)
    
    cosines <- sweep(cos(2 * pi * tf_mat), MARGIN = 2, Re(val()$z_weighted), `*`)
    colnames(cosines) <- paste0("cos_", seq(1, length(val()$f)))
    
    sines <- sweep(sin(2 * pi * tf_mat), MARGIN = 2, -Im(val()$z_weighted), `*`)
    sines <- sines[, val()$f > 0 & val()$f < val()$fs / 2]
    colnames(sines) <- paste0("sin_", seq(2, sum(val()$f < val()$fs / 2)))
    
    dat_coeff <- cbind(data.frame(Time = time), cosines, sines)
    dat_coeff <- pivot_longer(dat_coeff, -c(Time), names_pattern = "(.*)_(.*)", names_to = c("type", "freq"), values_to = "Amplitude", names_transform = list(freq = as.integer))
    dat_coeff$type <- factor(dat_coeff$type)
    dat_coeff$freq <- factor(dat_coeff$freq, labels = as.character(format(val()$f, digits = 2)))
    dat_coeff
  })
  
  # Prepare (co-)sines, upsampled to increase visibility of (co-)sines
  dat_bg <- reactive({
    time_high <- seq(0, val()$n * 10 - 1) / val()$fs / 10
    tf_mat <- sweep(matrix(rep(time_high, length(val()$f)), nrow = length(time_high), ncol = length(val()$f), byrow = FALSE), MARGIN = 2, val()$f, `*`)
    
    cosines <- sweep(cos(2 * pi * tf_mat), MARGIN = 2, Re(val()$z_weighted), `*`)
    colnames(cosines) <- paste0("cos_", seq(1, length(val()$f)))
    
    sines <- sweep(sin(2 * pi * tf_mat), MARGIN = 2, -Im(val()$z_weighted), `*`)
    sines <- sines[, val()$f > 0 & val()$f < val()$fs / 2]
    colnames(sines) <- paste0("sin_", seq(2, sum(val()$f < val()$fs / 2)))
    
    dat_bg <- cbind(data.frame(Time = time_high), cosines, sines)
    dat_bg <- pivot_longer(dat_bg, -c(Time), names_pattern = "(.*)_(.*)", names_to = c("type", "freq"), values_to = "Amplitude", names_transform = list(freq = as.integer))
    dat_bg$type <- factor(dat_bg$type)
    dat_bg$freq <- factor(dat_bg$freq, labels = as.character(format(val()$f, digits = 2)))
    dat_bg
  })
  
  # Time domain: Data
  output$plot_td <- renderPlot({
    dat_td <- data.frame(Amplitude = val()$data)
    dat_td$Time <- seq(0, length(dat_td$Amplitude) - 1) / val()$fs
    
    ggplot(dat_td, aes(x = Time, y = Amplitude)) + 
      geom_line(group = 1) +
      geom_point()
  })
  
  # Frequency domain: Weights
  output$plot_weights <- renderPlot({
    dat_weights <- data.frame(Frequency = val()$f, Weight = val()$weights)

    ggplot(dat_weights, aes(x = Frequency, y = Weight)) + 
      geom_line(group = 1, linetype = "dashed", color = okabe[7]) +
      geom_point(, color = okabe[7]) +
      ylim(c(0, 1))
  })
  
  # Frequency domain
  output$plot_fd <- renderPlot({
    my_fd <- data.frame(Frequency = val()$f,
                        Rectangular_real_raw = Re(val()$z),
                        Rectangular_imag_raw = -Im(val()$z),
                        Magnitude_abs_raw = abs(val()$z),
                        Phase_abs_raw = atan2(Im(val()$z), Re(val()$z)),
                        Rectangular_real_filtered = Re(val()$z_weighted),
                        Rectangular_imag_filtered = -Im(val()$z_weighted),
                        Magnitude_abs_filtered = abs(val()$z_weighted),
                        Phase_abs_filtered = atan2(Im(val()$z_weighted), Re(val()$z_weighted))
    )
    my_fd <- my_fd %>%
      pivot_longer(-Frequency, values_to = "value", names_to = c("Response", "type", "filtered"), names_sep = "_")
    
    my_fd$Response <- factor(my_fd$Response, levels = c("Rectangular", "Magnitude", "Phase"))
    my_fd$type <- factor(my_fd$type, levels = c("real", "imag", "abs"))
    my_fd$filtered <- factor(my_fd$filtered, levels = c("raw", "filtered"))
    # print(my_fd, n = 30)

    ggplot(my_fd, aes(x = Frequency, y = value, color = type, linetype = filtered)) +
      geom_line() +
      scale_color_manual(values = okabe[c(2, 3, 1)]) +
      geom_point() +
      facet_grid(rows = vars(Response), scales="free_y") +
      guides(color = FALSE) +
      theme(axis.title.y=element_blank(), legend.title=element_blank(), legend.position = "bottom")
  })
  
  # Frequency domain: (Co-)sines
  # https://stackoverflow.com/questions/25752529/how-to-automatically-right-size-ggplot-in-shiny
  output$plot_coeff <- renderPlot({
    ggplot(dat_coeff(), aes(x = Time, y = Amplitude, col = type)) +
      geom_line(data = dat_bg()) +
      geom_point() +
      scale_color_manual(values = okabe[c(2, 3)]) +
      facet_grid(rows = vars(freq), cols = vars(type))
  })
  
  # Time domain: Synthesized data
  output$plot_synth <- renderPlot({
    dat_weighted<- dat_coeff() %>%
      group_by(Time) %>%
      summarise(Amplitude = sum(Amplitude))
    dat_weighted$Difference <- val()$data - dat_weighted$Amplitude
    
    ggplot(dat_weighted, aes(x = Time, y = Amplitude)) + 
      geom_line(group = 1) +
      geom_point() +
      ylim(min(c(dat_weighted$Amplitude, dat_weighted$Difference)), max(c(dat_weighted$Amplitude, dat_weighted$Difference)))
  })
  
  # Time domain: Raw - filtered data
  output$plot_diff <- renderPlot({
    dat_weighted <- dat_coeff() %>%
      group_by(Time) %>%
      summarise(Amplitude = sum(Amplitude))
    dat_weighted$Difference <- val()$data - dat_weighted$Amplitude
    
    ggplot(dat_weighted, aes(x = Time, y = Difference)) + 
      geom_line(group = 1) +
      geom_point() +
      ylim(min(c(dat_weighted$Amplitude, dat_weighted$Difference)), max(c(dat_weighted$Amplitude, dat_weighted$Difference)))
  })
}

shinyApp(ui = ui, server = server)
