# This script analyzes the outputs of a matlab file and performs basic
# calculations to create data visualizations for a paper about my discovery of a
# new method for directly calculating the electrical transport pathway 
# magnitudes of eptihelial tissues.
# Colby F. Lewallen 2023

# LIBRARIES ----
# load and import critical libraries for data processing and plotting
library(readxl)
library(tidyverse) # for stats
library(reshape2) # dcast
library(scales) # to access break formatting functions
library(RColorBrewer) # for color palette
library(ggdist) # for the stat_halfeye packages
library(colorspace) # for the darken and lighten scripts
library(purrr) # for stats
library(coin) # for stats
library(broom) # for stats
library(fuzzyjoin) # for merging the EIS and voltage data together
library(rstatix) # for wilcox_effsize


# library(ggthemes)
# library(cowplot)
# 
# library(svglite) # for saving SVG
# library(grid)
# #library(devtools)
# #install_github("vqv/ggbiplot")
# library(ggbiplot)

# PATHS ----
# define all of the pathways for importing and saving results
MAIN_path <- file.path("G:","My Drive", "Colby", "EIS and Pipette", fsep = .Platform$file.sep)
FITS_path <- file.path(MAIN_path,"fits", fsep = .Platform$file.sep)
SAVE_path <- file.path(MAIN_path,"figures", fsep = .Platform$file.sep)
DATA_path <- file.path(MAIN_path,"data", fsep = .Platform$file.sep)
setwd(FITS_path)

# COLOR PALETTE ---
cbPalette <- c("#000000","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
color_vector = c(cbPalette,color_vector) # add the colorblind palette if desired

# PLOT THEME ----
theme_publication <- function(base_size=11, base_family="sans") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line.x = element_line(colour="black"),
            axis.line.y = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            panel.spacing = unit(0.39*1.1, "cm"), # 0.39 cm for 11 point font width
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.4, "cm"),
            legend.spacing = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold"),
    ))
  
}

# Import ----
# load in the data frames from the csvs and xlsx files used for calculations in 
# the paper

# lookup table (LT) is the table that matches a measurement to a specific
# combination of circuit parameters for the model cell. Further, it specifies
# what chambers were used for the measurements for cross-sectional area reasons
LT <- read_excel("EIS and Pipette lookup table.xlsx")

# fit results (fr) are the resistances and capacitances that were fit to the
# data for each measurement
fr <- read.csv("fit results.csv")

# EIS results (er) are the individual data points measured for all model and
# biological cells.
er <- read.csv("EIS results.csv") 

# CALCULATE ----
# calculate new terms and values that are relevant in the remainder of the code

# COMPLETE_ID:
# create a unique id for each measurement based on the folder name and
# measurement index (the underscore after the filename in my naming convention).
# the folder name is the date-time starting in with years and working down to 
# the second the folder was created. Each unique recording was repeated up to 3 
# 3 times and each repeat corresponds to an incremental meas_idx. Therefore, the
# folder name is not sufficient to uniquely id any particular NOVA measurement. 
# call the unique identifier the complete_ID
er$complete_ID <- paste(er$meas_ID,er$meas_idx)
fr$complete_ID <- paste(fr$meas_ID,fr$meas_idx)

# TER, TEC, VDR, and complete apical and baso-lateral resistances:
# calculate the bulk parameters that allow for easy comparisons with the
# traditional methods of intracellular ephys (e.g., VDR)
# ground truth: TER TEC VDR Rapical Rbasal
fr$p_TER <- (fr$p_7*(fr$p_3+fr$p_5)/(fr$p_3+fr$p_5+fr$p_7))
fr$p_TEC <- fr$p_4*fr$p_6/(fr$p_4+fr$p_6)
fr$p_VDR <- fr$p_3/fr$p_5
fr$p_Rapical <- fr$p_1+fr$p_3
fr$p_Rbasal <- fr$p_2+fr$p_5

# EIS+Pipette Electrode ratio: fit TER TEC VDR Rapical Rbasal using Va/TEP
fr$pg_absEr_TER <- (fr$pg_absEr_7*(fr$pg_absEr_3+fr$pg_absEr_5)/(fr$pg_absEr_3+fr$pg_absEr_5+fr$pg_absEr_7))
fr$pg_absEr_TEC <- fr$pg_absEr_4*fr$pg_absEr_6/(fr$pg_absEr_4+fr$pg_absEr_6)
fr$pg_absEr_VDR <- fr$pg_absEr_3/fr$pg_absEr_5
fr$pg_absEr_Rapical <- fr$pg_absEr_1+fr$pg_absEr_3
fr$pg_absEr_Rbasal <- fr$pg_absEr_2+fr$pg_absEr_5

# EIS+Pipette Electrode ratio (max 1 kHz): fit TER TEC VDR Rapical Rbasal using Va/TEP
fr$pg_absEr1000_TER <- (fr$pg_absEr1000_7*(fr$pg_absEr1000_3+fr$pg_absEr1000_5)/(fr$pg_absEr1000_3+fr$pg_absEr1000_5+fr$pg_absEr1000_7))
fr$pg_absEr1000_TEC <- fr$pg_absEr1000_4*fr$pg_absEr1000_6/(fr$pg_absEr1000_4+fr$pg_absEr1000_6)
fr$pg_absEr1000_VDR <- fr$pg_absEr1000_3/fr$pg_absEr1000_5
fr$pg_absEr1000_Rapical <- fr$pg_absEr1000_1+fr$pg_absEr1000_3
fr$pg_absEr1000_Rbasal <- fr$pg_absEr1000_2+fr$pg_absEr1000_5

# EIS+Pipette Membrane ratio: fit TER TEC VDR Rapical Rbasal using Va/Vb
fr$pg_absZr_TER <- (fr$pg_absZr_7*(fr$pg_absZr_3+fr$pg_absZr_5)/(fr$pg_absZr_3+fr$pg_absZr_5+fr$pg_absZr_7))
fr$pg_absZr_TEC <- fr$pg_absZr_4*fr$pg_absZr_6/(fr$pg_absZr_4+fr$pg_absZr_6)
fr$pg_absZr_VDR <- fr$pg_absZr_3/fr$pg_absZr_5
fr$pg_absZr_Rapical <- fr$pg_absZr_1+fr$pg_absZr_3
fr$pg_absZr_Rbasal <- fr$pg_absZr_2+fr$pg_absZr_5

# EIS+Pipette Membrane ratio (max 1 kHz): fit TER TEC VDR Rapical Rbasal using Va/Vb
fr$pg_absZr1000_TER <- (fr$pg_absZr1000_7*(fr$pg_absZr1000_3+fr$pg_absZr1000_5)/(fr$pg_absZr1000_3+fr$pg_absZr1000_5+fr$pg_absZr1000_7))
fr$pg_absZr1000_TEC <- fr$pg_absZr1000_4*fr$pg_absZr1000_6/(fr$pg_absZr1000_4+fr$pg_absZr1000_6)
fr$pg_absZr1000_VDR <- fr$pg_absZr1000_3/fr$pg_absZr1000_5
fr$pg_absZr1000_Rapical <- fr$pg_absZr1000_1+fr$pg_absZr1000_3
fr$pg_absZr1000_Rbasal <- fr$pg_absZr1000_2+fr$pg_absZr1000_5

# EIS+no pipette (2P-EIS): fit TER TEC VDR Rapical Rbasal without extra data/measurements
fr$pg_absnone_TER <- (fr$pg_absnone_7*(fr$pg_absnone_3+fr$pg_absnone_5)/(fr$pg_absnone_3+fr$pg_absnone_5+fr$pg_absnone_7))
fr$pg_absnone_TEC <- fr$pg_absnone_4*fr$pg_absnone_6/(fr$pg_absnone_4+fr$pg_absnone_6)
fr$pg_absnone_VDR <- fr$pg_absnone_3/fr$pg_absnone_5
fr$pg_absnone_Rapical <- fr$pg_absnone_1+fr$pg_absnone_3
fr$pg_absnone_Rbasal <- fr$pg_absnone_2+fr$pg_absnone_5

# RENAME METHODS:
# Mapping of method names and short names
methods <- list(absnone = list(name = "2P-EIS", short_name = "absnone"),
                absZr = list(name = "3P-EIS", short_name = "absZr"),
                absEr = list(name = "Electrode ratio", short_name = "absEr"),
                absZr1000 = list(name = "3P-EIS (1000 Hz)", short_name = "absZr1000"),
                absEr1000 = list(name = "Electrode ratio (1000 Hz)", short_name = "absEr1000")
)

# Traditional Ephys: TER, and VDR
# we can extract what the traditional ephys measurements would have revealed
# about any circuit or epithelia. To do this, we simply need to subset the data
# and only analyze the measurements that were performed at the lowest 
# frequency measured with EIS. 
#
# for every unique measurement, get the idx of the lowest frequency (f), and set 
# TER to the real impedance only, and calculate the VDR by dividing |Va|/|Vb|
unique_IDs <- unique(fr$complete_ID)
fr$pg_TERold <- NA
fr$pg_VDRold <- NA
for (unique_ID in unique_IDs) {
  er_subset <- subset(er,complete_ID == unique_ID) # extract this specific measurement
  
  # delete Va_amp or Vb_amp greater than a threshold that seems to be caused by meas error
  er_subset <- subset(er_subset,abs(Va_amp)<0.01 | abs(Vb_amp)<0.01 )
  er_subset <- subset(er_subset,f==min(er_subset$f)) # get the minimum frequency
  
  # store old data into the frequency response data
  fr$pg_TERold[fr$complete_ID == unique_ID] <- er_subset$x
  fr$pg_VDRold[fr$complete_ID == unique_ID] <- er_subset$Va_amp/er_subset$Vb_amp
}

# Error analysis:
# calculate the median error for the data set with and without the pipette
# grab the fit and membrane parameters from the data frame. Start with a 
# function that can build a data frame that contains the relevant fit and 
# actual circuit values, per measurement, per condition (no pipette, Er, Zr, )
process_one_condition <- function(condition_data, p, method_name, method_short_name) {
  param_names <- c('RsolA', 'RsolB', 'Ra', 'Ca', 'Rb', 'Cb', 'Rs', 'TER', 'TEC', 'VDR', 'Rapical', 'Rbasal')
  param_subscript <- c('1', '2', '3', '4', '5', '6', '7', 'TER', 'TEC', 'VDR', 'Rapical', 'Rbasal')
  
  # Initialize empty data frame with an additional column for 'param_type'
  mdf <- data.frame(complete_ID = character(), cross_section = numeric(), variable = character(), actual = numeric(), guess = numeric(), diff = numeric(), error = numeric(), method = character(), param_type = factor())
  
  for (i in seq_along(param_subscript)) {
    variable <- param_subscript[i]
    actual_col_name <- paste0('p_', variable)
    guess_col_name <- paste0('pg_', method_short_name, '_', variable)
    
    # Determine the type of parameter for 'param_type' column
    scale <- 1 # used to scale capacitance values to uF
    if (param_names[i] %in% c('RsolA', 'RsolB', 'Ra', 'Rb', 'Rs', 'TER', 'VDR', 'Rapical', 'Rbasal')) {
      param_type <- "Resistance"
    } else if (param_names[i] %in% c('Ca', 'Cb', 'TEC')) {
      param_type <- "Capacitance"
      scale <- 1e6
    } else {
      param_type <- "Unknown"
    }
    
    mdf_this_param <- data.frame(
      complete_ID = p$complete_ID,
      cross_section = p$cross_section,
      variable = param_names[i],
      actual = p[[actual_col_name]]*scale,
      guess = condition_data[[guess_col_name]]*scale,
      diff = (condition_data[[guess_col_name]] - p[[actual_col_name]])*scale,
      error = (condition_data[[guess_col_name]] - p[[actual_col_name]]) / p[[actual_col_name]] * 100,
      method = method_name,
      param_type = factor(param_type, levels = c("Resistance", "Capacitance", "Unknown")) # assign the factor level
    )
    
    mdf <- rbind(mdf, mdf_this_param)
  }
  
  # add extra analysis for the error of old method
  df2 <- fr %>% filter(complete_ID %in% p$complete_ID)
  
  # TER
  mdf_this_param <- data.frame(
    complete_ID = p$complete_ID,
    cross_section = p$cross_section,
    variable = "TER",
    actual = p$p_TER,
    guess = df2$pg_TERold,
    diff = df2$pg_TERold-p$p_TER,
    error = (df2$pg_TERold-p$p_TER)/p$p_TER*100,
    method = "Traditional",
    param_type = factor("Resistance", levels = c("Resistance", "Capacitance", "Unknown"))
  )
  mdf <- rbind(mdf,mdf_this_param)
  
  # VDR
  mdf_this_param <- data.frame(
    complete_ID = p$complete_ID,
    cross_section = p$cross_section,
    variable = "VDR",
    actual = p$p_VDR,
    guess = df2$pg_VDRold,
    diff = df2$pg_VDRold-p$p_VDR,
    error = (df2$pg_VDRold-p$p_VDR)/p$p_VDR*100,
    method = "Traditional",
    param_type = factor("Resistance", levels = c("Resistance", "Capacitance", "Unknown"))
  )
  mdf <- rbind(mdf,mdf_this_param)
  
  return(mdf)
}

# merge together all unique conditions into one large data frame
process_all_conditions <- function(conditions, p, methods) {
  mdf_combined <- data.frame()
  
  for (method_short_name in names(conditions)) {
    method_name <- methods[[method_short_name]]$name
    condition_data <- conditions[[method_short_name]]
    
    result <- process_one_condition(condition_data, p, method_name, method_short_name)
    mdf_combined <- rbind(mdf_combined, result)
  }
  
  return(mdf_combined)
}

# filter out your specific measurements conditions and build the list for input 
# into process_all_conditions function that we made above this section
create_conditions_data <- function (fr) { 
  pg_absnone <- fr[grepl("pg_absnone_", colnames(fr)) | colnames(fr) == "complete_ID"]
  pg_absZr <- fr[grepl("pg_absZr_", colnames(fr)) | colnames(fr) == "complete_ID"]
  pg_absEr <- fr[grepl("pg_absEr_", colnames(fr)) | colnames(fr) == "complete_ID"]
  pg_absZr1000 <- fr[grepl("pg_absZr1000_", colnames(fr)) | colnames(fr) == "complete_ID"]
  pg_absEr1000 <- fr[grepl("pg_absEr1000_", colnames(fr)) | colnames(fr) == "complete_ID"]
  p <- fr[grepl("p_", colnames(fr)) | colnames(fr) == "complete_ID" | colnames(fr) == "cross_section"]
  conditions <- list(absnone = pg_absnone, absZr = pg_absZr, absEr = pg_absEr, absZr1000= pg_absZr1000, absEr1000=pg_absEr1000)
  return(
    list(p=p, 
         conditions=conditions)
    )
}

# Process all conditions
fit_error_conditions <- create_conditions_data(fr)
fit_error <- process_all_conditions(conditions=fit_error_conditions$conditions, p=fit_error_conditions$p, methods=methods)

# railing error conditions
# RAILING ANALYSIS
# in a large portion of the data, there was at least one condition that railed
# do analysis on the error and uniqueness of this railing condition
membrane_rail = 500000*0.99
capacitance_rail = 1e-4*0.99
solution_rail = 250*0.99
fr_rail <- fr %>%
  filter(pg_absZr_3 > membrane_rail | 
           pg_absZr_5 > membrane_rail | 
           pg_absZr_7 > membrane_rail | 
           pg_absZr_4 > capacitance_rail | 
           pg_absZr_6 > capacitance_rail|
           pg_absZr_1 > solution_rail |
           pg_absZr_2 > solution_rail)
fr_rail

# Process all conditions that railed
fit_error_conditions_rail <- create_conditions_data(fr_rail)
fit_error_rail <- process_all_conditions(conditions=fit_error_conditions_rail$conditions, p=fit_error_conditions_rail$p, methods=methods)
fit_error_rail

# count the unique cases with railing
fr_filtered_count <- fr %>%
  rowwise() %>%
  mutate(count_conditions_met = sum(c(pg_absZr_3, pg_absZr_5, pg_absZr_7) > membrane_rail)) %>%
  ungroup() %>%
  filter(chamber %in% c("model")) %>%
  filter(include_in_analysis %in% c(1) )
unique(fr_filtered_count$count_conditions_met) # did more than one ever rail at a time?
sum(fr_filtered_count$count_conditions_met) # verify how many railed?


# CLEAN ----
# remove randomly measured circuits from analysis. further, only keep the
# circuits where the ground truth parameters were analyzed. The main experiment 
# was designed around the first 982 permutations of the model circuit
LT_idx_not_orig_permute_table <- LT$meas_ID[981:nrow(LT)] 
frc <- fr %>% 
  filter(chamber %in% c("model")) %>%
  filter(include_in_analysis %in% c(1) )

# mirror the results in the EIS data
erc <- er %>%
  filter(complete_ID %in% unique(frc$complete_ID))

# PRE-PROCESS
# search the data base of measurements for representative measurements. These
# data should show low, medium, and high resistance cells with variable 
# imaginary and Er values
frc_low <- frc %>%
  filter(p_TER < 600 & p_TER > 500)

frc_medium <- frc %>%
  filter(p_TER > 1200 & p_TER < 2200)

frc_high <- frc %>%
  filter(p_TER >2200 & p_TER < 5000)

# simple function to manually calculate the TER of any sample
calculate_TER <- function(RsolA, RsolB, Ra, Rb, Rs) {
       TER = RsolA+RsolB+Rs*(Ra+Rb)/(Ra+Rb+Rs)
       return(TER)
}

# simple function to manually calculate the asymptotes of any sample
calculate_Zr <- function(RsolA, RsolB, Ra, Rb, Rs) {
  print("low frequency: ")
  print(round((Ra*Rs+RsolA*(Ra+Rb+Rs))/(Rb*Rs+RsolB*(Ra+Rb+Rs)),1))
  print("high frequency:")
  print(round(RsolA/RsolB,1))
}

# EXAMPLE DATA PLOTS | FIGURE 2 (A) and FIGURE 3 ----
# visualization of the residuals. The residuals are a measurement of the error 
# of each single frequency between the fit value and the measured value. This
# data gives insight into the accuracy of the overall fit and if the model 
# equation did a good job fitting the measured data. 

create_EIS_plot_and_dataframe <- function(EIS_example, er, fit_error, methods, theme_publication, unit="") {
  
  # Step 1: Rename observations
  observation_labels <- list(
    measured = c("x", "y", "Zr"),
    fit = c("xfit_absZr", "yfit_absZr", "Z3fit_absZr")
  )
  
  # Step 2: Rename variable labels
  variable.labs <- c("Real", "Imaginary", "Membrane ratio")
  names(variable.labs) <- c("x", "y", "Zr")
  
  # Step 3: Extract specific measurement set
  EIS_example_data <- er %>%
    filter(complete_ID %in% EIS_example) %>%
    melt(id.vars = c("complete_ID", "f","cross_section"), 
         measure.vars = c("x", "y", "Zr", "xfit_absZr", "yfit_absZr", "Z3fit_absZr")) %>%
    left_join(., stack(observation_labels), by = c(variable = "values")) %>%
    mutate(obs = ind, ind = NULL) %>%
    mutate(variable = recode(variable, xfit_absZr = "x", yfit_absZr = "y", Z3fit_absZr = "Zr")) %>%
    dcast(complete_ID + f + cross_section + variable ~ obs, value.var = "value")
  
  # Step 4: Factor and order
  EIS_example_data$complete_ID <- factor(EIS_example_data$complete_ID, levels = EIS_example)
  EIS_example_data$variable <- factor(EIS_example_data$variable, levels = c("x", "y", "Zr"))
  
  # Step 5: Rename EIS recordings
  observation.labs <- paste("Example", 1:length(EIS_example))
  names(observation.labs) <- EIS_example
  
  # optional step, scale the resistance and reactance to kOhms
  if (unit %in% c("kOhms")) {
    EIS_example_data <- EIS_example_data %>%
      mutate(
        across(
          .cols = c("fit", "measured"),
          .fns = ~ case_when(
            variable %in% c("x", "y") ~ . / 1000,
            TRUE ~ .
          )
        )
      )
  }
  
  # asymptotes code
  # min
  EIS_example_data %>% filter(complete_ID %in% unique(EIS_example_data$complete_ID)[1]) %>% filter(f %in% min(EIS_example_data$f))
  # max
  EIS_example_data %>% filter(complete_ID %in% unique(EIS_example_data$complete_ID)[1]) %>% filter(f %in% max(EIS_example_data$f))
  
  # Step 6: Plot
  ax_example_data <- ggplot(EIS_example_data, aes(x = f)) +
    geom_line(aes(y = fit*cross_section), linewidth = 0.25) +
    geom_point(aes(y = measured*cross_section), 
               shape = 22, 
               fill = 'white', 
               color = 'black', 
               size = 1, 
               alpha = 0.7) +
    facet_grid(variable ~ complete_ID, 
               scales = "free_y", 
               labeller = labeller(variable = variable.labs, complete_ID = observation.labs)) +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x, n = 5),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_continuous(labels = comma) +
    labs(
      title = paste(EIS_example, collapse = ' - '),
      y = element_blank(),
      x = "Frequency (Hz)"
    ) +
    theme_publication()
  
  # Step 7: Retrieve relevant measurements
  figure_fr <- fit_error %>% 
    filter(complete_ID %in% EIS_example) %>% 
    filter(method %in% methods[['absZr']]$name)
  figure_all_params <- figure_fr
  
  # Step 8: Filter rows
  figure_fr_filtered <- figure_fr %>%
    filter(variable %in% c("RsolA", "RsolB", "Ra", "Rb", "Rs", "Ca", "Cb"))
  
  # Step 9: Create sub columns
  figure_fr_long <- figure_fr_filtered %>%
    select(complete_ID, variable, actual, guess) %>%
    gather(key = "type", value = "value", -complete_ID, -variable) %>%
    unite("complete_ID_type", complete_ID, type, sep = "_")
  
  # Step 10: Spread columns
  table_EIS_example_actual_and_guess_vals <- figure_fr_long %>%
    spread(key = "complete_ID_type", value = "value") %>%
    select(variable, everything())
  
  # Step 11: filter and extract just the residual data for the example fits
  EIS_example_data_residual <- er %>%
    filter(complete_ID %in% EIS_example) %>%
    melt(
      id.vars = c("complete_ID","f","cross_section"), 
      measure.vars = c("residual_xfit_absZr","residual_yfit_absZr","residual_Z3fit_absZr")
    ) %>%
    mutate(variable = recode(variable, residual_xfit_absZr = "x", residual_yfit_absZr = "y", residual_Z3fit_absZr = "Zr")) %>%
    merge(EIS_example_data, by=c("complete_ID","f","cross_section","variable"))
  
  # Step 12: force the order of the plots
  EIS_example_data_residual$variable <- factor(EIS_example_data_residual$variable, levels = c("x", "y", "Zr"))
  
  # Step 13: normalize the residuals for each group
  EIS_example_data_residual <- EIS_example_data_residual %>%
    group_by(complete_ID, variable) %>%
    mutate(value = abs(value) ) %>%
    ungroup()
  
  # optional step, scale the resistance and reactance to kOhms
  if (unit %in% c("kOhms")) {
    EIS_example_data_residual <- EIS_example_data_residual %>%
      mutate(
        across(
          .cols = c("value"),
          .fns = ~ case_when(
            variable %in% c("x", "y") ~ . / 1000,
            TRUE ~ .
          )
        )
      )
  }
  
  # Step 14: calculate mean and standard deviation of the data, per real, imaginary, and impedance data
  table_EIS_example_data_residual <- EIS_example_data_residual %>%
    group_by(variable) %>%
    summarise(
      n = n(),
      mean = mean(value, na.rm=FALSE),
      sd = sd(value, na.rm=FALSE),
      error_median = quantile(value*cross_section, probs = 0.5),
      errorIQR_25 = quantile(value*cross_section, probs = 0.25), 
      errorIQR_75 = quantile(value*cross_section, probs = 0.75)
    )
  table_EIS_example_data_residual$variable = factor(table_EIS_example_data_residual$variable, levels = c("x", "y", "Zr"))
  
  # Step 15: plot the example eis residuals
  ax_example_residuals <- ggplot(EIS_example_data_residual, aes(x = 1, y = value*cross_section)) +
    geom_jitter(size = 0.2, alpha = 0.2, shape = 21, fill = "white", color = "black", width = 0.21) +
    geom_boxplot(alpha = 0.5,outlier.shape = NA) +
    #facet_grid(variable~complete_ID,
    facet_grid(variable~.,
               scales="free_y",
               labeller = labeller(variable = variable.labs,complete_ID =  observation.labs)) +
    labs(
      title = "Boxplot of Residuals",
      y = "Magnitude of Residuals",
      x = element_blank()
    ) +
    coord_flip(clip = "off") +
    scale_y_log10(
      breaks = scales::trans_breaks("log10", function(x) 10^x),
      labels = scales::trans_format("log10", scales::math_format(10^.x))
    ) +
    theme_publication() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank())
  
  return(list(
    ax_example_data = ax_example_data, 
    ax_example_residuals = ax_example_residuals, 
    table_EIS_example_actual_and_guess_vals = table_EIS_example_actual_and_guess_vals,
    table_EIS_example_data_residual = table_EIS_example_data_residual,
    table_EIS_example_actual_and_guess_vals_all = figure_all_params))
}

# use the function
# EIS_example <- unique(frc$complete_ID)
# EIS_example <- unique(frc_low$complete_ID) # use to find an low TER examples
# EIS_example <- unique(frc_medium$complete_ID) # use to find medium TER example
# EIS_example <- unique(frc_high$complete_ID) # use to find high TER examples
# EIS_example <- c("20221025_170935 2") # Example for residuals
EIS_example <- c("20221214_164937 2","20221103_175857 3","20221212_155405 3") # fig 2 a,b,c

# use the custom function to create the results data frame that has a table of 
# data and the corresponding figure
result_examples <- create_EIS_plot_and_dataframe(EIS_example, er, fit_error, methods, theme_publication, unit="kOhms")
result_examples$ax_example_data
result_examples$ax_example_residuals

# all data residual analysis
result_all_fits <- create_EIS_plot_and_dataframe(unique(frc$complete_ID), er, fit_error, methods, theme_publication, unit="kOhms")
#result_all_fits$ax_example_residuals

# sample 1
Biology_example1 <- c("20211021_183356 7", "20211014_165514 10", "20211014_165514 19")
result_biology_examples1 <- create_EIS_plot_and_dataframe(Biology_example1, er, fit_error, methods, theme_publication, unit="")
result_biology_examples1$ax_example_data
result_biology_examples1$ax_example_residuals
result_biology_examples1$table_EIS_example_actual_and_guess_vals_all

# biological examples before during and after ATP stimulation: Sample 2
Biology_example2 <- c("20211022_171622 10", "20211022_171622 13", "20211022_171622 19")
result_biology_examples2 <- create_EIS_plot_and_dataframe(Biology_example2, er, fit_error, methods, theme_publication, unit="")
result_biology_examples2$ax_example_data
result_biology_examples2$ax_example_residuals
result_biology_examples2$table_EIS_example_actual_and_guess_vals_all

# sample 3
Biology_example3 <- c("20211014_165514 5", "20211014_165514 7", "20211014_165514 17")
result_biology_examples3 <- create_EIS_plot_and_dataframe(Biology_example3, er, fit_error, methods, theme_publication, unit="")
result_biology_examples3$ax_example_data
result_biology_examples3$ax_example_residuals
result_biology_examples3$table_EIS_example_actual_and_guess_vals_all

# TABLE 1 ----
# calculation of the median errors for the different methods of calculating the 
# goodness of fit for each circuit parameter. 

# clean the fit error to only contain the fits of interest
fit_error_abs <- fit_error
fit_error_abs$diff <- abs(fit_error_abs$diff)
fit_error_abs$error <- abs(fit_error_abs$error)
fit_error_clean <- fit_error_abs %>%
  filter(complete_ID %in% unique(frc$complete_ID)) %>%
  filter(method %in% c(methods[['absnone']]$name, methods[['absZr']]$name)) 

# summarize the mean and median errors
fit_error_medianIQR <- fit_error_clean %>% 
  group_by(param_type,variable,method) %>%
  summarise(n=n(),
            error_median = quantile(error, probs = 0.5),
            errorIQR_25 = quantile(error, probs = 0.25), 
            errorIQR_75 = quantile(error, probs = 0.75)
  )
fit_error_medianIQR # display in the command window

# railing examples stats
# get complete IDs where Ra railed
fit_error_abs_rail_complete_ID <- fr %>%
  filter(pg_absZr_7 > membrane_rail) # modify term to change analysis
complete_ID_rail <- unique(fit_error_abs_rail_complete_ID$complete_ID)
fit_error_abs_rail <- fit_error_rail
fit_error_abs_rail$diff <- abs(fit_error_abs_rail$diff)
fit_error_abs_rail$error <- abs(fit_error_abs_rail$error)
fit_error_clean_abs <- fit_error_abs_rail %>%
  filter(complete_ID %in% unique(frc$complete_ID)) %>%
  filter(method %in% c(methods[['absZr']]$name)) %>%
  filter(complete_ID %in% complete_ID_rail) # in the case where Ra railed

# summarize the mean and median errors for railing conditions
fit_error_medianIQR_abs <- fit_error_clean_abs %>% 
  group_by(param_type,variable,method) %>%
  summarise(n=n(),
            error_median = quantile(error, probs = 0.5),
            errorIQR_25 = quantile(error, probs = 0.25), 
            errorIQR_75 = quantile(error, probs = 0.75)
  )
fit_error_medianIQR_abs # display in the command window

fit_error_wilcox_eff_size <- fit_error_clean %>%
  convert_as_factor(complete_ID,variable, method) %>%
  group_by(variable) %>%
  wilcox_effsize(error ~ method, paired = TRUE, alternative = "greater")

# Print the results
print(fit_error_wilcox_eff_size)

# Perform the paired Wilcoxon test and calculate the effect size
fit_error_wilcox_test <- fit_error_clean %>%
  convert_as_factor(complete_ID,variable, method) %>%
  group_by(variable) %>%
  summarise(
    wilcox_p = list(wilcox.test(error ~ method, paired = TRUE, alternative = "greater")),
    N = n()
  ) %>%
  rowwise(variable) %>%
  mutate(
    p_value = wilcox_p$p.value
  ) %>%
  ungroup() %>%
  select(-wilcox_p) %>%
  mutate(
    p_adj = p.adjust(p_value, method = "bonferroni")
  )

# Print the results
print(fit_error_wilcox_test)

# Set the parameters to be measured in the histogram plot
method_levels = c("3P-EIS", "2P-EIS")
# fit_error_subset <- fit_error_clean %>% filter(variable %in% c( "Ra", "Rb", "Rs","RsolA","RsolB")) %>%
#   filter(method %in% method_levels)
fit_error_subset <- fit_error_clean %>% filter(variable %in% c( "Ra", "Ca", "Rb", "Cb", "Rs", "RsolA", "RsolB")) %>%
  filter(method %in% method_levels)

# Define a custom labeller function
custom_labeller <- function(variable){
  return(as.character(variable))
}

# rain cloud plots of all errors
plot_levels = c("2P-EIS","3P-EIS")

# grouping the data
fit_error_subset <- fit_error_subset %>% arrange(complete_ID, method, variable)

# Convert the factor to numeric for jittering
fit_error_subset$numeric_method <- as.numeric(factor(fit_error_subset$method, levels = plot_levels))

# Create jittered positions
set.seed(1)  # for reproducibility
jitter_amount <- 0.15
fit_error_subset$jittered_x <- fit_error_subset$numeric_method +
  runif(nrow(fit_error_subset), -jitter_amount, jitter_amount)

# plot the boxplot of error
ax_error <- ggplot(fit_error_subset, aes(x=factor(method, levels = plot_levels), y = error)) + 
  facet_wrap(~variable, ncol = 2, scales = "free", dir = "h") +
  geom_boxplot(
    aes( color = factor(method, levels = plot_levels),
      color = after_scale(darken(color, .1, space = "HLS")),
      fill = after_scale(desaturate(lighten(color, .8), .4))),
    width = .55,
    outlier.shape = NA
  ) +
  # geom_line(
  #   aes(x = jittered_x,
  #       group = complete_ID,
  #       color = factor(method, levels = plot_levels)
  #   ),
  #   alpha = 0.02
  # ) +
  coord_flip(xlim = c(0.9, 2.1), clip = "off") +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))
  ) +
  geom_point(
    aes(x = jittered_x,
        color = factor(method, levels = plot_levels),
        color = after_scale(darken(color, .1, space = "HLS"))),
    fill = "white",
    shape = 21,
    stroke = 0.4,
    size = 0.25,
    alpha = .3
  ) +
  #annotation_logticks(sides = 'b') +
  scale_fill_manual(values = color_vector, guide = "none") +
  scale_color_manual(values = color_vector) +
  labs(
    fill = element_blank(),
    color = element_blank(),
    x = NULL,
    y = "Percent parameter estimation error",
    title = "Evaluation of fit error based on fitting method",
  ) +
  theme_publication() +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal",
    panel.grid.minor = element_line(colour="#f0f0f0"),
  )
ax_error

# FIGURE 4 ----
# traditional ephys data.
#sample 1 and 2 are Z8 and 3 is dominik AMDCD
biological_samples = c('20211021_183356','20211022_171622','20211014_165514')
plot_names = c(
  '20211021_183356' = "Sample 1",
  '20211022_171622' = "Sample 2",
  '20211014_165514' = "Sample 3")

fr_samples <- fr %>% filter(meas_ID %in% biological_samples)

sample_df <- function(data_ID,method_str, methods) {
  # grab the appropriate section from the calculated impedance values df
  df <- fr %>% filter(meas_ID %in% data_ID)
  
  # update all pertinent values in the data frame
  Rscale <- 1000 # convert to kOhms (1000)
  Cscale <- 1e-6 # convert to micro farads (1e-6)
  method <- paste0('pg_',method_str)
  df$method <- methods[[method_str]]$name
  df$RsolA <- df[[paste0(method,'_1')]]*df$cross_section
  df$RsolB <- df[[paste0(method,'_2')]]*df$cross_section
  df$Ra <- df[[paste0(method,'_3')]]*df$cross_section/Rscale
  df$Ca <- df[[paste0(method,'_4')]]/df$cross_section/Cscale
  df$Rb <- df[[paste0(method,'_5')]]*df$cross_section/Rscale
  df$Cb <- df[[paste0(method,'_6')]]/df$cross_section/Cscale
  df$Rs <- df[[paste0(method,'_7')]]*df$cross_section/Rscale
  df$TERold <- df$pg_TERold*df$cross_section
  df$VDRold <- df$pg_VDRold
  df$VDR <- (df$Ra*Rscale+df$RsolA)/(df$Rb*Rscale+df$RsolB)
  
  # calculate new terms
  df$TER <- df$Rs*(df$Ra+df$Rb)/(df$Ra+df$Rb+df$Rs)*Rscale
  df$TEC <- df$Ca*df$Cb/(df$Ca+df$Cb)
  
  # get ATP information
  COMMENT_path <- file.path(DATA_path,data_ID,paste0(data_ID,"_comments.txt"), fsep = .Platform$file.sep)
  comments = read.csv(COMMENT_path, sep = '\t')
  comment_idx = c(1,2) # for all files, make sure the 1st and second comment correspond to ATP
  
  # load original NOVA file for the absolute measurement times
  NOVA_path <- file.path(DATA_path,data_ID,paste0(data_ID,"_NOVAdata.txt"), fsep = .Platform$file.sep)
  NOVA <- read.csv(NOVA_path,sep='\t')
  NOVA_t = lapply(list(NOVA[,10]), function(z){ z[!is.na(z) & z != "" & z != "meas start time (s)"]})
  NOVA_t = as.numeric(NOVA_t[[1]])
  t0 = NOVA_t[1] # START TIME FOR ALL MEASUREMENTS - ABSOLUTE TIME
  
  # load the voltage data such as TEP
  VOLTAGES_path <- file.path(DATA_path,data_ID,paste0(data_ID,"_data.txt"), fsep = .Platform$file.sep)
  voltages <- read.csv(VOLTAGES_path,sep="\t")
  voltages$t <- (voltages$absolute.time.s.-t0)/60
  voltages <- voltages[voltages$t>=0, ] # remove data before the first EIS measurement
  voltages$Va <- voltages$Va..mV.
  voltages$Vb <- voltages$Vb..mV.
  voltages$TEP <- voltages$TEP..mV.
  voltages$TTL <- voltages$TTL.signal > 0.5
  voltages$meas_ID <- rep(df$meas_ID[1],nrow(voltages))
  voltages <- subset(voltages,TTL==FALSE) # replace voltages with NaN if EIS was recording
  voltages$technique <- rep("3P-Classical",nrow(voltages))
  
  # sync all meas times to start time
  df$t = (NOVA_t-t0)/60 # minutes
  comment_times = (comments[comment_idx,1] - t0)/60
  
  # Perform the fuzzy join
  df <- difference_inner_join(df, voltages, by = "t", max_dist = Inf) %>%
    arrange(abs(t.x - t.y)) %>%
    group_by(t.x) %>%
    slice_head(n = 1) %>%
    ungroup()
  df$t <- df$t.x # fix renamed variables
  df$meas_ID <- df$meas_ID.x
  
  # add the atp start and end times to the data frame
  df$ATP_start <- rep(comment_times[1],nrow(df))
  df$ATP_end <- rep(comment_times[2],nrow(df))
  
  # calculate the equivalent Isc (uA)
  df$Isc <- (df$TEP*1e-3/df$TER)*1e6 # convert to uA/cm^2
  
  # calculate the ratio of shunt conductance to transcellular conductance
  df$GpGt <- (df$Ra+df$Rb)/df$Rs
  
  # calculate the cell resistances
  df$Rcell <- df$Ra + df$Rb
  
  ###################### CROP SPECIFICALLY THE AMDCD FILE WHERE Ba2+ WAS APPLIED
  if (data_ID=='20211014_165514') {
    idx_max = which.min(abs(df$t-55))
    df = df[1:idx_max,]
    
    idx_max = which.min(abs(voltages$t-55))
    voltages = voltages[1:idx_max,]
  }
  
  # reshape the data frame to a structure suitable for plotting
  mdf <- melt(df, 
              id.vars = c("meas_ID","t","ATP_start","ATP_end"),
              measure.vars = c("RsolA","RsolB","Ra","Ca","Rb","Cb","Rs","TER","TEC","TERold","VDRold","VDR","TEP","Isc","Va","Vb","GpGt","Rcell")
              )
  
  mdf_voltages <- melt(voltages, 
              id.vars = c("meas_ID","t","technique"),
              measure.vars = c("TEP","Va","Vb")
  )
  
  # add a category property to each parameter to represent if it was an old value or new value
  technique = rep("Unknown",nrow(mdf))
  for (i in 1:nrow(mdf)) {
    if (mdf$variable[i] %in% c("Ra","Rb","Rs","TER","RsolA","RsolB","Ca","Cb","TEC","VDR","GpGt","Rcell")) {
      technique[i] = "3P-EIS"
    } else if (mdf$variable[i] %in% c("TEP","Isc")) {
      # technique[i] = "2P-Continuous"
      technique[i] = "3P-Classical"
    } else if (mdf$variable[i] %in% c("Va","Vb")) {
      # technique[i] = "3P-Continuous"
      technique[i] = "3P-Classical"
    } else {
      technique[i] = "3P-Classical"
      if (mdf$variable[i] %in% c("TERold")) {
        mdf$variable[i] = "TER"
      }
      if (mdf$variable[i] %in% c("VDRold")) {
        mdf$variable[i] = "VDR"
      }
    }
  }
  mdf$technique = technique
  
  return(list(
    mdf=mdf,
    voltages = mdf_voltages) )
}

# test the function
# "20221209_162613" <- example of a data frame that did not include ATP
for (bio_sample in biological_samples) {
  tmp <- sample_df(data_ID = bio_sample,method_str = 'absZr',methods=methods)
  if (bio_sample == biological_samples[1]){
    df_bio <- tmp$mdf
    df_bio_voltages <- tmp$voltages
  } else {
    df_bio <- rbind(df_bio,tmp$mdf)
    df_bio_voltages <- rbind(df_bio_voltages,tmp$voltages)
  }
}

# force the biological sample levels
df_bio$meas_ID <- factor(df_bio$meas_ID, levels = biological_samples)
df_bio_voltages$meas_ID <- factor(df_bio_voltages$meas_ID, levels = biological_samples)

# make a distinct data frame for the ATP rectangles
rect_data <- df_bio %>%
  select(meas_ID, ATP_start, ATP_end) %>%
  distinct()

# make the plot
plot_bio <- function (df,voltages,facet_by="",yaxis="") {
  variable_levels <- c("Ca","Cb","Ra","Rb","Rs","TER","RsolA","RsolB","VDR","TEP","Isc","Va","Vb","TEC","GpGt","Rcell")
  technique_levels <- c("3P-Classical","3P-EIS","2P-Continuous","3P-Continuous")
  
  num_vars <- length(unique(df$variable))
  num_techs <- length(unique(df$technique))
  
  unique_vars <- unique(df$variable) # get unique variables
  
  # Identify common columns
  common_cols <- intersect(names(df), names(voltages))
  
  # Filter voltages to only include common columns
  voltages <- voltages[, common_cols, drop = FALSE]
  
  # Check if unique variables are a subset of c("TEP", "Va", "Vb")
  use_voltages <- all(unique_vars %in% c("TEP", "Va", "Vb"))
  
  if (use_voltages) {
    voltages <- voltages %>% filter(variable %in% unique_vars)
    line_data <- voltages
  } else {
    line_data <- df
  }
  
  ax <- ggplot(df) +
    geom_rect(data = rect_data, 
              aes(xmin = ATP_start, xmax = ATP_end, ymin = -Inf, ymax = Inf),
              fill = color_vector[3], 
              alpha = 0.3) +
    xlim(0,69) +
    labs(title = element_blank(),
         x = "Time (min.)",
         y = "value",
         color = element_blank(),
         fill = element_blank()) +
    scale_color_manual(values = color_vector) +
    scale_fill_manual(values = color_vector) +
    theme_publication()
  
  if (facet_by %in% c("technique")) {
    ax <- ax + 
      geom_line(data = line_data, aes(x=t,y=value,color = factor(variable, levels = variable_levels)),alpha=0.99,linewidth = 0.1) +
      facet_grid(factor(technique, levels = technique_levels) ~ meas_ID, scales = "free_y")
    # change the default to white dot if only using a single facet to avoid confusion
    if (num_vars == 1) {
      ax <- ax +
        geom_point(aes(x=t,y=value,fill = factor(variable, levels = variable_levels)),
                   shape = 22,
                   fill = 'white',
                   color = 'black',
                   size = 1.25,
                   stroke = 0.4,
                   alpha = 1)
    } else {
      ax <- ax +
        geom_point(aes(x=t,y=value,fill = factor(variable, levels = variable_levels)),
                   shape = 22,
                   color = 'black',
                   size = 1.25,
                   stroke = 0.4,
                   alpha = 1)
    }
  } else {
    ax <- ax + 
      geom_line(data = line_data, aes(x=t,y=value,color = factor(technique, levels = technique_levels)),alpha=0.99,linewidth = 0.1) +
      facet_grid(factor(variable, levels = variable_levels) ~ meas_ID, scales = "free_y")
    if (num_techs == 1) {
      ax <- ax + 
        geom_point(aes(x=t,y=value,fill = factor(technique, levels = technique_levels)),
                 shape = 22,
                 fill = 'white',
                 color = 'black',
                 size = 1.25,
                 stroke = 0.4,
                 alpha = 0.7)
    } else {
      ax <- ax +
        geom_point(aes(x=t,y=value,fill = factor(technique, levels = technique_levels)),
                   shape = 22,
                   color = 'black',
                   size = 1.25,
                   stroke = 0.4,
                   alpha = 0.7)
    }
  }
  # make the y axis have log10 spacing to emphasize changes
  if (yaxis %in% "log10") {
    ax <- ax + coord_trans(y="log10")
  }
    
  return(ax)
}

# PLOT AND SAVE FIGURES ----
# in this section, we should save and plot all figures for the paper

# EXAMPLES of FITS
ax_example_data <- result_examples$ax_example_data
ax_example_data
file_name <- paste0("Example data.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = length(EIS_example)*1.5+1, height = 7.5,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 

# Biology Sample 2 example right before ATP, last data during ATP, and last data
ax_biology_examples <- result_biology_examples2$ax_example_data
ax_biology_examples
file_name <- paste0("Example biology data.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = length(EIS_example)*1.5+1, height = 7.5,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 


# TABLE OF EXAMPLE DATA AND FIT
table_EIS_example_actual_and_guess_vals <- result_examples$table_EIS_example_actual_and_guess_vals
table_EIS_example_actual_and_guess_vals
table_EIS_example_actual_and_guess_vals_all <- result_examples$table_EIS_example_actual_and_guess_vals_all
table_EIS_example_actual_and_guess_vals_all

# RESIDUALS PLOT | FIGURE 2B (supplemental)
# quantification of the residuals. The residuals are a measurement of the error 
# of each single frequency between the fit value and the measured value. This
# data gives insight into the accuracy of the overall fit and if the model 
# equation did a good job fitting the measured data. 
ax_example_residuals <- result_examples$ax_example_residuals
ax_example_residuals

file_name <- paste0("Example residual.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 1.5+1, height = 7.5,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 

# ALL FITS RESIDUALS PLOTS AND TABLES
# plot all fits of the data or show all of the residuals 
# table
table_EIS_example_data_residual = result_all_fits$table_EIS_example_data_residual
table_EIS_example_data_residual
ax_example_residuals <- result_all_fits$ax_example_residuals
ax_example_residuals
file_name <- paste0("Example residual (all).png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 5.5, height = 5,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 

# TOTAL ERROR ANALYSIS
ax_error
print(fit_error_wilcox_test)
print(fit_error_wilcox_eff_size)
fit_error_medianIQR
ggplot_build(ax_error)$data
file_name <- paste0("Fit error raincloud plot.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 5.5, height = 7.5,
       units = "in",
       dpi = 320,
       limitsize = TRUE)


# BIOLOGY PLOTS
df_voltages <- df_bio %>% filter(variable %in% c("TEP","Va"))
ax_voltages <- plot_bio(df_voltages,df_bio_voltages,facet_by = "variable",yaxis="")
ax_voltages
file_name <- paste0("Membrane voltages.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 5.5, height = 4,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 

df_traditional_values_compare <- df_bio %>% filter(variable %in% c("TER","VDR"))
ax_traditional_values_compare <- plot_bio(df_traditional_values_compare,df_bio_voltages,facet_by = "variable",yaxis="")
ax_traditional_values_compare
file_name <- paste0("3P-Classical vs 3P-EIS.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 5.5, height = 4,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 

df_3PEIS <- df_bio %>% filter(variable %in% c("Ca","Cb","Ra","Rb","Rs"))
ax_3PEIS <- plot_bio(df_3PEIS,df_bio_voltages,facet_by = "variable",yaxis="")
ax_3PEIS
file_name <- paste0("3P-EIS (membrane params).png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 5.5, height = 8,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 

df_bio_relevance <- df_bio %>% filter(variable %in% c("Ra","Rb","Rs","Rcell","TER"))
ax_bio_relevance <- plot_bio(df_bio_relevance,df_bio_voltages,facet_by = "variable",yaxis="")
ax_bio_relevance
file_name <- paste0("Bio-relevance plot.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 5.5, height = 8,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 


df_transepithelial <- df_bio %>% filter(variable %in% c("TER","TEC"))
ax_transepithelial <- plot_bio(df_transepithelial,df_bio_voltages,facet_by="variable",yaxis="")
ax_transepithelial
file_name <- paste0("Transepithelial.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 5.5, height = 4,
       units = "in",
       dpi = 320,
       limitsize = TRUE)

df_mem_resistances <- df_bio %>% filter(variable %in% c("Ra","Rb","Rs","VDR"))
ax_mem_resistances <- plot_bio(df_mem_resistances,df_bio_voltages,facet_by = "technique",yaxis="log10")
ax_mem_resistances
file_name <- paste0("Membrane resistances.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 5.5, height = 5,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 

df_mem_capacitances <- df_bio %>% filter(variable %in% c("Ca","Cb"))
ax_mem_capacitances <- plot_bio(df_mem_capacitances,df_bio_voltages,facet_by = "technique",yaxis="log10")
ax_mem_capacitances
file_name <- paste0("Membrane capacitances.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 5.5, height = 2.5,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 

df_sol_resistances <- df_bio %>% filter(variable %in% c("RsolA","RsolB"))
df_sol_resistances %>%
  group_by(variable) %>%
  summarise(
    n = n(),
    mean = mean(value, na.rm = FALSE),
    sd = sd(value, na.rm = FALSE),
    median = median(value, na.rm = FALSE),
    IQR25 = quantile(value, probs = 0.25),
    IQR75 = quantile(value, probs = 0.75)
  )
ax_sol_resistances <- plot_bio(df_sol_resistances,df_bio_voltages,facet_by = "technique",yaxis="log10")
ax_sol_resistances
file_name <- paste0("Solution resistances.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 5.5, height = 2.5,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 

df_alternative_assumptions <- df_bio %>% filter(variable %in% c("GpGt"))
df_alternative_assumptions %>%
  group_by(meas_ID) %>%
  summarise(
    n = n(),
    maxVal = max(value, na.rm = TRUE),
    minVal = min(value, na.rm = TRUE),
    diffVal = max(value, na.rm = TRUE) - min(value, na.rm = TRUE)
    )
ax_alternative_assumptions <- plot_bio(df_alternative_assumptions,df_bio_voltages,facet_by="variable",yaxis="log10")
ax_alternative_assumptions
file_name <- paste0("Prior assumptions testing.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 5.5, height = 2.5,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 

df_TEC <- df_bio %>% filter(variable %in% c("TEC"))
ax_TEC <- plot_bio(df_TEC,df_bio_voltages,facet_by = "variable",yaxis="")
ax_TEC
file_name <- paste0("TEC.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 5.5, height = 2.5,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 

df_Isc <- df_bio %>% filter(variable %in% c("Isc"))
ax_Isc <- plot_bio(df_Isc,df_bio_voltages,facet_by = "variable",yaxis="")
ax_Isc
file_name <- paste0("Isc.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 5.5, height = 2.5,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 

df_Rcell <- df_bio %>% filter(variable %in% c("Rcell"))
ax_Rcell <- plot_bio(df_Rcell,df_bio_voltages,facet_by = "variable",yaxis="")
ax_Rcell
file_name <- paste0("Rcell.png")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = 5.5, height = 2.5,
       units = "in",
       dpi = 320,
       limitsize = TRUE) 

