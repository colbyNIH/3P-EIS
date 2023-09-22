# VALIDATE THE EISandPipette.R data
# this must be run after loading er and fr from the EISandPipette.R function
library(ggthemes)
library(dplyr)

# # inspect the circuits that "failed" with a maximum possible resistance
# rsol_fail <- 50
# mem_r_fail <- 4.5e5
# mem_c_fail <- 5e2
# frc_absZr_fail <- frc %>%
#   filter(
#       (pg_absZr_3 > mem_r_fail |
#       pg_absZr_4 > mem_c_fail |
#       pg_absZr_5 > mem_r_fail |
#       pg_absZr_6 > mem_c_fail |
#       pg_absZr_7 > mem_r_fail)
#   )
# 
# frc_absEr_fail <- frc %>%
#   filter(
#          (pg_absEr_3 > mem_r_fail |
#          pg_absEr_4 > mem_c_fail |
#          pg_absEr_5 > mem_r_fail |
#          pg_absEr_6 > mem_c_fail |
#          pg_absEr_7 > mem_r_fail)
#   )
# 
# # frc_absZr_pass <- frc %>%
# #   filter(
# #       (pg_absZr_3 < mem_r_fail &
# #          pg_absZr_4 < mem_c_fail &
# #          pg_absZr_5 < mem_r_fail &
# #          pg_absZr_6 < mem_c_fail &
# #          pg_absZr_7 < mem_r_fail)
# #   )

# Function to simulate the circuit ----
simulate_circuit <- function(params, f) {
  w <- 2 * pi * f
  RsolA <- params$RsolA
  RsolB <- params$RsolB
  Za <- params$Ra / (1 + 1i * w * params$Ra * params$Ca)
  Zb <- params$Rb / (1 + 1i * w * params$Rb * params$Cb)
  Rs <- params$Rs
  Zab <- Za + Zb
  Zabs <- 1 / (1 / Zab + 1 / params$Rs)
  Zabs_Rsol <- params$RsolA + params$RsolB + Zabs
  Va <- (Za * Rs + RsolA * (Rs + Za + Zb)) / (Rs + Za + Zb)
  Vb <- (Zb * Rs + RsolB * (Rs + Za + Zb)) / (Rs + Za + Zb)
  Er = abs(Va)/abs(Zabs_Rsol)
  return(list(x=Re(Zabs_Rsol), y=Im(Zabs_Rsol), Er=Er))
}

# Iterate over each unique complete_ID ----
unique_complete_IDs <- unique(fr$complete_ID)
FIG_path <- file.path(DATA_path,"Validation Data (Er)")
setwd(FIG_path)
for (complete_id in unique_complete_IDs) {
  # Extract corresponding circuit parameters
  params <- fr[fr$complete_ID == complete_id, ]
  
  # Extract measured data for the complete_ID
  measured_data <- er %>%
    filter(complete_ID %in% complete_id) %>%
    melt(
      id.vars = c("complete_ID","f"), 
      measure.vars = c("x","y","Er","xfit_absEr","yfit_absEr","Z3fit_absEr")
    ) %>%
    left_join(., stack(observation_labels), by = c(variable = "values")) %>%
    mutate(obs = ind,ind = NULL) %>% # http://127.0.0.1:31341/graphics/plot_zoom_png?width=874&height=900
    mutate(variable = recode(variable, xfit_absEr = "x", yfit_absEr = "y", Z3fit_absEr = "Er")) %>%
    dcast(complete_ID+f+variable~obs, value.var = "value")
  
  # define the frequencies to measure
  f <- unique(measured_data$f)
  
  # Simulate the circuit
  simulated_data <- simulate_circuit(params, f)
  
  # Modify the structure of the simulated data
  simulated_df <- data.frame(
    complete_ID = rep(complete_id, 3 * length(f)),
    f = rep(f, 3),
    variable = rep(c("x", "y", "Er"), each = length(f)),
    simulated = c(simulated_data$x, simulated_data$y, simulated_data$Er)
  )
  
  # Merge the simulated data with the measured data
  combined_data <- merge(measured_data, simulated_df, by=c("complete_ID","f","variable"))
  combined_data$variable <- factor(combined_data$variable, levels = c("x", "y", "Er")) # force the order of the plots
  
  # Create the plot
  p <- ggplot(combined_data, aes(x = f)) +
    geom_line(aes(y = fit, linetype = "fit"), color = 'black', size = 0.75) +
    geom_line(aes(y = simulated, color = "simulated"), size = 0.75) +
    geom_point(aes(y = measured, shape = "measured"), fill = 'white', color = 'black', size = 2, alpha = 1) +
    facet_grid(variable ~ complete_ID, scales = "free_y", labeller = labeller(variable = variable.labs, complete_ID = observation.labs)) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 5), labels = trans_format("log10", math_format(10^.x))) +
    scale_y_continuous(labels = comma) +
    scale_linetype_manual(name = element_blank(), values = c(fit = "solid")) +
    scale_shape_manual(name = element_blank(), values = c(measured = 22)) +
    scale_color_manual(name = element_blank(), values = c(simulated = '#E69F00')) +
    labs(
      y = element_blank(),
      x = "Frequency (Hz)"
    ) +
    theme_Publication() +
    guides(
      linetype = guide_legend(override.aes = list(linewidth = 1.5, size = 5)),
      color = guide_legend(override.aes = list(linewidth = 1.5, size = 5))
    )
  p
  
  file_name <- paste0("Combined data ", complete_id, ".png")
  ggsave(file_name, width = 7.1 / 2, height = 7.1, units = "in")
}

