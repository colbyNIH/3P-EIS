# EIS+Pipette paper
# Import, processing, and analysis of all data acquired to analyze the technique
# outlined in the combined EIS and intracellular pipette paper

# Libraries ----
library(tidyverse)
library(ggthemes)
#library(devtools)
#install_github("vqv/ggbiplot")
library(ggbiplot)
library(reshape2)
library(cowplot)
library(scales) # to access break formatting functions
library(svglite) # for saving SVG
library(grid)
library(readxl)

# PATHS ----
MAIN_path <- file.path("G:","My Drive", "Colby", "EIS and Pipette", fsep = .Platform$file.sep)
DATA_path <- file.path(MAIN_path,"fits", fsep = .Platform$file.sep)
SAVE_path <- file.path(MAIN_path,"figures", fsep = .Platform$file.sep)
setwd(DATA_path)

# FIGURE SETTINGS ----
figure_width_large = 18 / 2.54 #PNAS large
figure_height_large = 22 / 2.54 #PNAS large
figure_width_medium = 11 / 2.54 #PNAS medium
figure_height_medium = 11 /2.54 #PNAS medium
figure_width_small = 9 / 2.54 #PNAS small
figure_height_small = 6 / 2.54 #PNAS small

# Import ----
LT <- read_excel("EIS and Pipette lookup table - fixed solution resistances.xlsx")
fr <- read.csv("fit results.csv")
er <- read.csv("EIS results.csv")

# NEW VARIABLES ----

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

# TER, TEC, and VDR:
# we can extract what the traditional ephys measurements would have revealed
# about any circuit or epithelia. To do this, we simply need to subset the data
# and only analyze the measurements that were performed at the lowest 
# frequency measured with EIS. 

# ground truth: TER TEC VDR Rapical Rbasal
fr$p_TER <- (fr$p_7*(fr$p_3+fr$p_5)/(fr$p_3+fr$p_5+fr$p_7))
fr$p_TEC <- fr$p_4*fr$p_6/(fr$p_4+fr$p_6)
fr$p_VDR <- fr$p_3/fr$p_5
fr$p_Rapical <- fr$p_1+fr$p_3
fr$p_Rbasal <- fr$p_2+fr$p_5

# EIS+Pipette Er: fit TER TEC VDR Rapical Rbasal
fr$pg_absEr_TER <- (fr$pg_absEr_7*(fr$pg_absEr_3+fr$pg_absEr_5)/(fr$pg_absEr_3+fr$pg_absEr_5+fr$pg_absEr_7))
fr$pg_absEr_TEC <- fr$pg_absEr_4*fr$pg_absEr_6/(fr$pg_absEr_4+fr$pg_absEr_6)
fr$pg_absEr_VDR <- fr$pg_absEr_3/fr$pg_absEr_5
fr$pg_absEr_Rapical <- fr$pg_absEr_1+fr$pg_absEr_3
fr$pg_absEr_Rbasal <- fr$pg_absEr_2+fr$pg_absEr_5

# EIS+Pipette Zr: fit TER TEC VDR Rapical Rbasal
fr$pg_absZr_TER <- (fr$pg_absZr_7*(fr$pg_absZr_3+fr$pg_absZr_5)/(fr$pg_absZr_3+fr$pg_absZr_5+fr$pg_absZr_7))
fr$pg_absZr_TEC <- fr$pg_absZr_4*fr$pg_absZr_6/(fr$pg_absZr_4+fr$pg_absZr_6)
fr$pg_absZr_VDR <- fr$pg_absZr_3/fr$pg_absZr_5
fr$pg_absZr_Rapical <- fr$pg_absZr_1+fr$pg_absZr_3
fr$pg_absZr_Rbasal <- fr$pg_absZr_2+fr$pg_absZr_5

# EIS+no pipette: fit TER TEC VDR Rapical Rbasal
fr$pg_absnone_TER <- (fr$pg_absnone_7*(fr$pg_absnone_3+fr$pg_absnone_5)/(fr$pg_absnone_3+fr$pg_absnone_5+fr$pg_absnone_7))
fr$pg_absnone_TEC <- fr$pg_absnone_4*fr$pg_absnone_6/(fr$pg_absnone_4+fr$pg_absnone_6)
fr$pg_absnone_VDR <- fr$pg_absnone_3/fr$pg_absnone_5
fr$pg_absnone_Rapical <- fr$pg_absnone_1+fr$pg_absnone_3
fr$pg_absnone_Rbasal <- fr$pg_absnone_2+fr$pg_absnone_5

# Traditional Ephys: TER, and VDR
# for every unique measurement, get the idx of the lowest frequency (f), and set 
# TER to the real impedance only, and calculate the VDR by dividing |Va|/|Vb|
unique_IDs <- unique(fr$complete_ID)
fr$pg_TERold <- NA
fr$pg_VDRold <- NA
for (i in 1:length(unique_IDs)) {
  er_subset <- subset(er,complete_ID == unique_IDs[i]) # extract this specific measurement
  
  # delete Va_amp or Vb_amp greater than a threshold that seems to be caused by meas error
  er_subset <- subset(er_subset,abs(Va_amp)<0.01 | abs(Vb_amp)<0.01 )
  er_subset <- subset(er_subset,f==min(er_subset$f)) # get the minimum frequency
  
  # store old data into the frequency response data
  fr$pg_TERold[fr$complete_ID == unique_IDs[i]] <- er_subset$x
  fr$pg_VDRold[fr$complete_ID == unique_IDs[i]] <- er_subset$Va_amp/er_subset$Vb_amp
}

# er_subset <- subset(er,complete_ID == '20211014_165514 7')
# er_subset <- subset(er_subset,abs(Va_amp)<0.01 | abs(Vb_amp)<0.01 )
# ggplot(er_subset,aes(x=log10(f),y=Va_amp)) +
#     geom_point()

# sample numbers for the paper: 20211021_183356, 20211022_171622, '20211014_165514'

fr_subset <- subset(fr,meas_ID == '20211014_165514')
mean(fr_subset$pg_absZr_3[6:7]*fr_subset$cross_section[1])
ggplot(fr_subset,aes(x=meas_idx,y=pg_absZr_3*cross_section)) +
  geom_point()

# CLEAN ----
# remove randomly measured circuits from analysis. further, only keep the
# circuits where the ground truth parameters were analyzed
meas_ID_remove <- LT$meas_ID[981:nrow(LT)] 
frc <- fr %>% 
  filter(circuitValuesKnown %in% c(1)) %>%
  filter( !( meas_ID %in% c(meas_ID_remove) ) )

# mirror the results in the EIS data
erc <- er %>%
  filter(complete_ID %in% unique(frc$complete_ID))
  
# # ERROR ----
# # grab the fit and membrane parameters from the data frame
# p <- frc[grepl("p_",colnames(frc))]
# pg_absnone <- frc[grepl("pg_absnone_",colnames(frc))]
# pg_absZr <- frc[grepl("pg_absZr_",colnames(frc))]
# 
# # delta p
# dpg_absnone <- pg_absnone-p
# dpg_absZr <- pg_absZr-p
# 
# # normalize delta p
# dpg_absnone_error <- dpg_absnone/p
# dpg_absZr_error <- dpg_absZr/p
# median((frc$pg_TERold-frc$p_TER)/(frc$p_TER))*100
# median((frc$pg_VDRold-frc$p_VDR)/(frc$p_VDR))*100
# median( ( (frc$pg_absZr_TER)/(frc$pg_absZr_7) - (frc$p_TER)/(frc$p_7) ) / (frc$p_TER)/(frc$p_7) )*100
# 
# # median errors
# library(matrixStats)
# colMedians((as.matrix(dpg_absnone_error))*100)
# colMedians((as.matrix(dpg_absZr_error))*100)
# 
# # PLOT: median errors
# data <- dpg_absZr_error[,1:7]*100
# rep_str = c('pg_absZr'='',
#             '_1'='RsolA',
#             '_2'='RsolB',
#             '_3'='Ra',
#             '_4'='Ca',
#             '_5'='Rb',
#             '_6'='Cb',
#             '_7'='Rs',
#             '_TER' = 'TER',
#             '_TEC' = 'TEC',
#             '_VDR' = 'VDR',
#             '_Rapical' = 'Rapical',
#             '_Rbasal' = 'Rbasal'
# )
# colnames(data) <- str_replace_all(colnames(data), rep_str)
# mdf <- melt(
#   data,
#   varnames = names(dimnames(data)),
#   na.rm = FALSE,
#   as.is = FALSE,
#   value.name = "value" 
# )
# ggplot(mdf) +
#   geom_histogram(aes(x=value),
#                  bins = 20) + # bins = 40
#   labs(
#     x="Percent error (%)",
#     y="Count"
#   ) +
#    scale_x_continuous(
#      limits = c(-200,200), # limits -200,200
#      breaks = seq(-150,150,by=150) # breaks -150,150 by 150
#    ) +
#   facet_grid(~variable) +
#   theme_stata()
# file_name <- paste("All fit errors (cropped).svg", sep = "")
# ggsave(file_name,
#        path = SAVE_path,
#        plot = last_plot(),
#        width = figure_width_large, height = figure_height_medium,
#        units = "in",
#        limitsize = TRUE) 

# EXAMPLE DATA ----
# select a representative recording and show the fit and predicted values. Start 
# by defining the example data to fit. Also define figure settings.
EIS_example <- c("20221117_173404 3")
EIS_example <- c("20221110_174503 3")
EIS_example <- c("20221018_191151 3")
EIS_example <- c("20221209_150838 3") # fig 3 C 
EIS_example <- c("20221103_173551 2")
EIS_example <- c("20221117_171136 2") # fig 3 B <- impedance ratio data is incorrect?
EIS_example <- c("20221221_175143 2") # fig 3 A
EIS_example <- c("20221117_171136 3")
EIS_example <- c("20221220_200140 1")
EIS_example <- c("20220923_170313 2") # <- example of a bad fit
EIS_example <- c("20221216_170401 1") # <- another example of a bad fit
EIS_example <- c("20221209_160624 3") # example of bad fit
EIS_example <- c("20221025_170935 2") # main paper figure
EIS_example <- c("20221221_175143 2","20221117_171136 2","20221209_150838 3") # fig 3 a,b,c

# rename simple variable labels in the output table  to more descriptive 
# captions
variable.labs <- c("Real","Imaginary", "Zr")
names(variable.labs) <- c("x", "y", "Zr")

# rename all observations, both fit and actual to have the same x y and Zr names
observation_labels <- list(actual=c("x","y","Zr"),
                           fit=c("xfit_absZr","yfit_absZr","Z3fit_absZr"))

# to avoid confusion in the paper, rename the EIS recording complete_ID to 
# "Example N"
observation.labs <- c("Example 1", "Example 2","Example 3")
names(observation.labs) <- EIS_example

# FILTER AND PREPROCESS EXAMPLE FREQUENCY DATA
EIS_example_data <- er %>%
  filter(complete_ID %in% EIS_example) %>%
  melt(
    id.vars = c("complete_ID","f"), 
    measure.vars = c("x","y","Zr","xfit_absZr","yfit_absZr","Z3fit_absZr")
  ) %>%
  left_join(., stack(observation_labels), by = c(variable = "values")) %>%
  mutate(obs = ind,ind = NULL) %>%
  mutate(variable = recode(variable, xfit_absZr = "x", yfit_absZr = "y", Z3fit_absZr = "Zr")) %>%
  dcast(complete_ID+f+variable~obs, value.var = "value")

fit_example_data <- fr %>%
  filter(complete_ID %in% EIS_example)
# # actual values for circuit
# fit_example_data$RsolA
# fit_example_data$RsolB
# fit_example_data$Ra
# fit_example_data$Ca
# fit_example_data$Rb
# fit_example_data$Cb
# fit_example_data$Rs
# fit_example_data$p_TER
# fit_example_data$p_totalR <- fit_example_data$p_TER + fit_example_data$RsolA + fit_example_data$RsolB
# fit_example_data$p_totalR
# fit_example_data$p_solutionR <- fit_example_data$RsolA + fit_example_data$RsolB
# fit_example_data$p_solutionR
# fit_example_data$p_tauR <- (fit_example_data$Ra*fit_example_data$Ca)/(fit_example_data$Rb*fit_example_data$Cb)
# fit_example_data$p_tauR
# fit_example_data$impedance_ratio_lowF <- (fit_example_data$Ra+fit_example_data$RsolA)/(fit_example_data$Rb+fit_example_data$RsolB)
# fit_example_data$impedance_ratio_lowF
# fit_example_data$impedance_ratio_highF <- (fit_example_data$RsolA)/(fit_example_data$RsolB)
# fit_example_data$impedance_ratio_highF
# # predicted values for circuit
# fit_example_data$pg_absZr_1
# fit_example_data$pg_absZr_2
# fit_example_data$pg_absZr_3
# fit_example_data$pg_absZr_4
# fit_example_data$pg_absZr_5
# fit_example_data$pg_absZr_6
# fit_example_data$pg_absZr_7
# fit_example_data$pg_absZr_TER
# fit_example_data$pg_absZr_totalR <- fit_example_data$pg_absZr_TER + fit_example_data$pg_absZr_1 + fit_example_data$pg_absZr_2
# fit_example_data$pg_absZr_totalR
# fit_example_data$pg_absZr_solutionR <- fit_example_data$pg_absZr_1 + fit_example_data$pg_absZr_2
# fit_example_data$pg_absZr_solutionR
# fit_example_data$pg_absZr_tauR <- (fit_example_data$pg_absZr_3*fit_example_data$pg_absZr_4)/(fit_example_data$pg_absZr_5*fit_example_data$pg_absZr_6)
# fit_example_data$pg_absZr_tauR
# fit_example_data$pg_impedance_ratio_lowF <- (fit_example_data$pg_absZr_1+fit_example_data$pg_absZr_3)/(fit_example_data$pg_absZr_2+fit_example_data$pg_absZr_5)
# fit_example_data$pg_impedance_ratio_lowF
# 
# # search for the maximum measured impedance ratio??
# max_Zr <- er %>%
#   melt(
#     id.vars = c("complete_ID","f"), 
#     measure.vars = c("x","y","Zr","xfit_absZr","yfit_absZr","Z3fit_absZr")
#   ) %>%
#   left_join(., stack(observation_labels), by = c(variable = "values")) %>%
#   mutate(obs = ind,ind = NULL) %>%
#   mutate(variable = recode(variable, xfit_absZr = "x", yfit_absZr = "y", Z3fit_absZr = "Zr")) %>%
#   dcast(complete_ID+f+variable~obs, value.var = "value") %>%
#   filter(variable %in% c("Zr")) %>%
#   filter(f %in% c(0.69556))

# EIS EXAMPLE AND FIT DATA
EIS_example_data$complete_ID <- factor(EIS_example_data$complete_ID, levels = EIS_example)
ggplot(EIS_example_data,aes(x=f)) +
  geom_line(aes(y=fit),linewidth=0.25) +
  geom_point(aes(y=actual),
             shape = 22,
             fill='white',
             color='black',
             size = 2,
             alpha = 1) +
  facet_grid(variable~complete_ID,
             scales="free_y",
             labeller = labeller(variable = variable.labs,complete_ID =  observation.labs)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x,n=5),
                labels = trans_format("log10", math_format(10^.x))) +
  scale_y_continuous(labels = comma) +
  labs(
    title = paste(EIS_example,collapse = ' - '),
    y = element_blank(),
    x = "Frequency (Hz)"
  ) +
  theme_stata(scheme = "s1color",base_size = 12) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1))
file_name <- paste0("Example data ",paste(EIS_example,collapse = ' - '),".svg")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = figure_width_large/2, height = figure_height_large,
       units = "in",
       limitsize = TRUE) 


# EXAMPLE RESIDUALS
# filter and extract just the residual data for the example fits
EIS_example_data_residual <- er %>%
  filter(complete_ID %in% EIS_example) %>%
  melt(
    id.vars = c("complete_ID","f"), 
    measure.vars = c("residual_xfit_absZr","residual_yfit_absZr","residual_Zrfit_absZr")
  ) %>%
  mutate(variable = recode(variable, residual_xfit_absZr = "x", residual_yfit_absZr = "y", residual_Zrfit_absZr = "Zr"))

# calculate mean and standard deviation of the data, per real, imaginary, and impedance data
aggregate(EIS_example_data_residual$value, list(EIS_example_data_residual$variable), FUN=length)
aggregate(EIS_example_data_residual$value, list(EIS_example_data_residual$variable), FUN=mean)
aggregate(EIS_example_data_residual$value, list(EIS_example_data_residual$variable), FUN=sd)

# normalize the residuals for each group
EIS_example_data_residual <- EIS_example_data_residual %>%
  group_by(complete_ID, variable) %>%
  mutate(value = value / max(abs(value))) %>%
  ungroup()

# plot the example eis residuals - normalized
ggplot(EIS_example_data_residual, aes(x = value)) +
  geom_histogram(binwidth = 0.05, fill=color_vector[2], color = "black", alpha = 0.7) +
  facet_grid(variable ~ complete_ID, scales = "free_y", labeller = labeller(variable = variable.labs,complete_ID =  observation.labs)) +
  labs(
    title = "Histogram of Residuals",
    y = "Frequency",
    x = "Normalized Value"
  ) +
  theme_stata(scheme = "s1color",base_size = 12) +
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1))

# ggplot(EIS_example_data_residual,aes(x=value)) +
#   facet_grid(variable~complete_ID,
#              labeller = labeller(variable = variable.labs,complete_ID =  observation.labs)) +
#   geom_histogram(
#     bins = 50,
#     color = 'black',
#     size = 0.25,
#     fill=color_vector[2],
#     alpha = 0.4
#   ) +
#   labs(
#     x= "Residual",
#     y= element_blank()
#   ) +
#   xlim(c(-5,5)) +
#   theme_stata(scheme = "s1color",base_size = 12)+
#   theme(axis.text.y = element_text(angle = 0, vjust = 0.5, hjust=1))
file_name <- paste0("Example residuals ",paste(EIS_example,collapse = ' - '),".svg")
ggsave(file_name,
       path = SAVE_path,
       width = figure_width_large/2,
       height = figure_height_large,
       units = "in")

# plot ALL eis residuals
EIS_all_data_residual <- erc %>%
  melt(
    id.vars = c("complete_ID","f"), 
    measure.vars = c("residual_xfit_absZr","residual_yfit_absZr","residual_Z3fit_absZr")
  ) %>%
  mutate(variable = recode(variable, residual_xfit_absZr = "x", residual_yfit_absZr = "y", residual_Z3fit_absZr = "Z3"))


# calculate mean and standard deviation of the data, per real, imaginary, and impedance data
aggregate(EIS_all_data_residual$value, list(EIS_all_data_residual$variable), FUN=length)
aggregate(EIS_all_data_residual$value, list(EIS_all_data_residual$variable), FUN=mean)
aggregate(EIS_all_data_residual$value, list(EIS_all_data_residual$variable), FUN=sd)

ggplot(EIS_all_data_residual,aes(x=value)) +
  facet_grid(variable~.,
             labeller = labeller(variable = variable.labs)) +
  geom_histogram(
    bins = 50,
    color = 'black',
    fill=color_vector[2],
    alpha = 0.4
  ) +
  labs(
    x= "Residual",
    y= element_blank()
  ) +
  xlim(c(-15,15)) +
  theme_lewallen()
file_name <- paste("All residuals.svg")
ggsave(file_name,
       width = figure_width,
       height = figure_height,
       units = "in")


# # ERROR bars for example
# observation_labels <- list(actual=c("Ra","Rb","Rs","Ca","Cb"),
#                            fit=c("Rag","Rbg","Rsg","Cag","Cbg"),
#                            min_bnd=c("Rag_min","Rbg_min","Rsg_min","Cag_min","Cbg_min"),
#                            max_bnd=c("Rag_max","Rbg_max","Rsg_max","Cag_max","Cbg_max")
#                            )
# variable_labels <- list(resistance=c("Ra","Rb","Rs",
#                                      "Rag","Rbg","Rsg",
#                                      "Rag_min","Rbg_min","Rsg_min",
#                                      "Rag_max","Rbg_max","Rsg_max"),
#                         capacitance=c("Ca","Cb",
#                                       "Cag","Cbg",
#                                       "Cag_min","Cbg_min",
#                                       "Cag_max","Cbg_max")
#                         )
# 
# # EIS EXAMPLE ERROR BARS
# # rename the fit and bound values for all variables to have the same column 
# # name. Then add tags that correspond to if the variable was a capacitance or 
# # resistance AND what the variable name for the actual, fit, and bounds were AND
# # if the fit was actually in range
# fit_example_data <- fr %>%
#   filter(complete_ID %in% EIS_example) %>%
#   melt(
#     id.vars = c("complete_ID"),
#     measure.vars = c("Ra","Rag","Rag_min","Rag_max",
#                      "Rb","Rbg","Rbg_min","Rbg_max",
#                      "Rs","Rsg","Rsg_min","Rsg_max",
#                      "Ca","Cag","Cag_min","Cag_max",
#                      "Cb","Cbg","Cbg_min","Cbg_max")
#   ) %>%
#   left_join(., stack(observation_labels), by = c(variable = "values")) %>%
#   mutate(tag_obs = ind,ind = NULL) %>%
#   left_join(., stack(variable_labels), by = c(variable = "values")) %>%
#   mutate(tag_var = ind,ind = NULL) %>%
#   mutate(variable = recode(variable, 
#                            Rag = "Ra",Rag_min = "Ra",Rag_max = "Ra", 
#                            Rbg = "Rb",Rbg_min = "Rb",Rbg_max = "Rb", 
#                            Rsg = "Rs",Rsg_min = "Rs",Rsg_max = "Rs",
#                            Cag = "Ca",Cag_min = "Ca",Cag_max = "Ca",
#                            Cbg = "Cb",Cbg_min = "Cb",Cbg_max = "Cb")
#          ) %>%
#   dcast(complete_ID+tag_var+variable~tag_obs, value.var = "value") %>%
#   mutate(isinCIrange = (actual>min_bnd & actual<max_bnd))
# 
# # plot the fit with confidence interval
# # resistance
# fit_example_data %>%
#   filter(tag_var %in% c('resistance')) %>%
#   ggplot(aes(x=factor(variable))) +
#   geom_errorbar(aes(ymin=min_bnd,ymax=max_bnd),
#                 width=.5,size=0.5,alpha=1) +
#   geom_point(aes(y=fit),
#              alpha = 1,
#              shape = 21,
#              size = 2,
#              fill=color_vector[2]) +
#   geom_point(aes(y=actual),
#              alpha = 1,
#              shape = 22,
#              size = 2,
#              fill='black') +
#   ylim(c(0,NA)) +
#   labs(
#     x=element_blank(),
#     y=expression(Resistance~(Omega))
#   ) +
#   theme_lewallen()
# file_name <- paste0("Example fit - ",EIS_example[1]," Resistance",".svg")
# ggsave(file_name,
#        width = figure_width,
#        height = figure_height/2,
#        units = "in")
#   
# # capacitance
# fit_example_data %>%
#   filter(tag_var %in% c('capacitance')) %>%
#   ggplot(aes(x=factor(variable))) +
#   geom_errorbar(aes(ymin=min_bnd*1e6,ymax=max_bnd*1e6),
#                 width=.5,size=0.5,alpha=1) +
#   geom_point(aes(y=fit*1e6),
#              alpha = 1,
#              shape = 21,
#              size = 2,
#              fill=color_vector[2]) +
#   geom_point(aes(y=actual*1e6),
#              alpha = 1,
#              shape = 22,
#              size = 2,
#              fill='black') +
#   ylim(c(0,NA)) +
#   labs(
#     x=element_blank(),
#     y=expression(Capacitance~(mu*F))
#   ) +
#   theme_lewallen()
# file_name <- paste0("Example fit - ",EIS_example[1]," Capacitance",".svg")
# ggsave(file_name,
#        width = figure_width,
#        height = figure_height/2,
#        units = "in")


# ANALYSIS OF ALL FITS ----
# import and analyze the goodness of fits for all circuit board measurements, 
# also analyze if the fit data was within the confidence interval

# 
# 
# observation_labels <- list(actual=c("Ra","Rb","Rs","Ca","Cb"),
#                            fit=c("Rag","Rbg","Rsg","Cag","Cbg"),
#                            min_bnd=c("Rag_min","Rbg_min","Rsg_min","Cag_min","Cbg_min"),
#                            max_bnd=c("Rag_max","Rbg_max","Rsg_max","Cag_max","Cbg_max")
# )
# variable_labels <- list(resistance=c("Ra","Rb","Rs",
#                                      "Rag","Rbg","Rsg",
#                                      "Rag_min","Rbg_min","Rsg_min",
#                                      "Rag_max","Rbg_max","Rsg_max"),
#                         capacitance=c("Ca","Cb",
#                                       "Cag","Cbg",
#                                       "Cag_min","Cbg_min",
#                                       "Cag_max","Cbg_max")
# )
# 
# fit_data <- frc %>%
#   melt(
#     id.vars = c("complete_ID","ci_alphaabs","TER_actual","TEC_actual","Rsol_r","tau_r"),
#     measure.vars = c("Ra","Rag","Rag_min","Rag_max",
#                      "Rb","Rbg","Rbg_min","Rbg_max",
#                      "Rs","Rsg","Rsg_min","Rsg_max",
#                      "Ca","Cag","Cag_min","Cag_max",
#                      "Cb","Cbg","Cbg_min","Cbg_max")
#   ) %>%
#   left_join(., stack(observation_labels), by = c(variable = "values")) %>%
#   mutate(tag_obs = ind,ind = NULL) %>%
#   left_join(., stack(variable_labels), by = c(variable = "values")) %>%
#   mutate(tag_var = ind,ind = NULL) %>%
#   mutate(variable = recode(variable, 
#                            Rag = "Ra",Rag_min = "Ra",Rag_max = "Ra", 
#                            Rbg = "Rb",Rbg_min = "Rb",Rbg_max = "Rb", 
#                            Rsg = "Rs",Rsg_min = "Rs",Rsg_max = "Rs",
#                            Cag = "Ca",Cag_min = "Ca",Cag_max = "Ca",
#                            Cbg = "Cb",Cbg_min = "Cb",Cbg_max = "Cb")
#   ) %>%
#   dcast(complete_ID+tag_var+variable+ci_alphaabs+TER_actual+TEC_actual+Rsol_r+tau_r~tag_obs, value.var = "value") %>%
#   mutate(isinCIrange = (actual>min_bnd & actual<max_bnd))
# fit_data$CI = 100*(1-fit_data$ci_alphaabs)
# fit_data$rel_span = log10(fit_data$max_bnd-fit_data$min_bnd)/log10(fit_data$actual)
# rel_span_reject_threshold = 1.5
# fit_data$rejected = abs(fit_data$rel_span) > rel_span_reject_threshold
# fit_data$pError = abs(fit_data$actual-fit_data$fit)/fit_data$fit*100






# # count the number of each unique condition
# uniqueVar = unique(fit_data$variable)
# uniqueCI = unique(fit_data$ci_alphaabs)
# 
# df_ci <-as.data.frame(matrix(nrow=length(uniqueCI),ncol=length(uniqueVar)+2))
# colnames(df_ci) <- c("n","CI",uniqueVar)
# 
# df_ci_reject = as.data.frame(matrix(nrow = 1, ncol = 7)) # multiply number of rows by measure.vars in fit_data_properties
# colnames(df_ci_reject) <- c("var","CI","param","value","nReject","nTotal","pPass")
# count = 0
# for (var in uniqueVar) {
#   fit_data_var = fit_data %>%
#     subset(variable==var) %>%
#     filter(rejected %in% c(FALSE))
#     #filter(TER_actual > 150)
#   
#   fit_data_properties = fit_data %>%
#     mutate(param=variable,variable = NULL) %>%
#     melt(id.vars = c("param","CI","rejected"),
#          measure.vars = c("actual","TER_actual","TEC_actual","Rsol_r","tau_r")) %>%
#     subset(param==var)
#   
#   idx_var = which(colnames(df_ci)==var)
#   for (alpha in uniqueCI) {
#     # set the alpha value
#     idx_alpha = which(uniqueCI==alpha)
#     df_ci[idx_alpha,2] = (1-alpha)*100
#     
#     # calculate the parameter inside confidence interval bounds
#     fit_data_var_alpha = subset(fit_data_var,ci_alphaabs==alpha)
#     
#     df_ci[idx_alpha,1] = nrow(fit_data_var_alpha)
#     nTotal = nrow(fit_data_var_alpha)
#     nInRange = sum(fit_data_var_alpha$isinCIrange)
#     df_ci[idx_alpha,idx_var] = nInRange/nTotal*100
#     
#     # evalute if the span is effected by a particular property of the circuit, 
#     # e.g., TER, TEC, etc.
#     fit_data_properties_alpha = subset(fit_data_properties,CI==(100*(1-alpha)))
#     # get the unique properties to be measured
#     uniqueVariables = unique(fit_data_properties_alpha$variable)
#     for (this_variable in uniqueVariables) {
#       fdpa = fit_data_properties_alpha %>%
#         subset(variable == this_variable)
#       uniqueValues = unique(fdpa$value)
#       for (thisUniqueValue in uniqueValues) {
#         count = count+1
#         unique_fdpa = fdpa %>%
#           subset(value == thisUniqueValue)
#         nThisUnique = nrow(unique_fdpa)
#         nReject = sum(unique_fdpa$rejected==TRUE)
#         
#         df_ci_reject[count,1] = var
#         df_ci_reject[count,2] = 100*(1-alpha)
#         df_ci_reject[count,3] = this_variable
#         df_ci_reject[count,4] = thisUniqueValue
#         df_ci_reject[count,5] = nReject
#         df_ci_reject[count,6] = nThisUnique
#         df_ci_reject[count,7] = (nThisUnique-nReject)/nThisUnique*100
#       }
#     }
#   }
# }
# 
# 
# df_ci %>% 
#   melt(id.vars = c("CI"),
#        measure.vars = uniqueVar) %>%
#   ggplot(aes(x=CI,y=value,fill=variable,color=variable)) +
#   geom_line() +
#   geom_point(shape='circle',size=2) +
#   scale_fill_manual(values=color_vector) +
#   scale_color_manual(values=color_vector) +
#   labs(
#     x="Confidence interval (%)",
#     y="Actual value within confidence interval (%)"
#   ) +
#   ylim(c(0,NaN)) +
#   theme_lewallen()
# file_name <- paste("Confidence interval convergence plot.svg")
# ggsave(file_name,
#        width = figure_width,
#        height = figure_height,
#        units = "in")



# 
# # analyze error of single fit ----
# fit_example_data <- fr %>%
#   filter(complete_ID %in% EIS_example) %>%
#   filter(ci_alphaabs %in% c(0.001))
# fit_example_data$TEC_actual
# fit_example_data$TEC12
# fit_example_data$TECabs
# 
# 
# 
# # R^2 ----
# # circuit_board_values <- c("RsolA","RsolB","Ra","Ca","Rb","Cb","Rs")
# # fit_values <- c("RsolAg","RsolBg","Rag","Cag","Rbg","Cbg","Rsg")
# # sim_values <- c("RsolAg_sim","RsolBg_sim","Rag_sim","Cag_sim","Rbg_sim","Cbg_sim","Rsg_sim")
# circuit_board_values <- c("Ra","Ca","Rb","Cb","Rs")
# fit_values <- c("Rag","Cag","Rbg","Cbg","Rsg")
# sim_values <- c("Rag_sim","Cag_sim","Rbg_sim","Cbg_sim","Rsg_sim")
# 
# circuit_board <- log10(frc[,circuit_board_values])
# fit <- log10(frc[,fit_values])
# sim <- log10(frc[,sim_values])
# 
# # calculate the SSE for each parameter
# SSE <- colSums((fit-circuit_board)^2)
# SSE_sim <- colSums((sim-circuit_board)^2)
# 
# # calculate the SST for each parameter in the circuit board. Start by finding 
# # the mean of each parameter
# circuit_board_means <- colMeans(circuit_board)
# SST <- colSums((fit-circuit_board_means)^2)
# SST_sim <- colSums((sim-circuit_board_means)^2)
# 
# # calculate the R2 for each parameter
# R2i <- 1-SSE/SST
# R2i_sim <- 1-SSE_sim/SST_sim
# 
# # determine the standard deviation of each parameter, independent of mean
# error <- (fit-circuit_board)/circuit_board
# error_mean <- colMeans(error)
# error_sd <- sqrt(colSums((error-error_mean)^2)/nrow(error))
# 
# errorp <- (10^fit-10^circuit_board)/10^circuit_board
# error_meanp <- colMeans(errorp)
# error_sd <- sqrt(colSums((errorp-error_meanp)^2)/nrow(errorp))
# 
# 
# # round the circuit board values to the nearest values
# test <- signif(circuit_board,digits=3)
# 
# 
# # R^2plots ----
# col_names <- colnames(circuit_board)
# round_digits_ax <- 1
# plot_list <- list()
# #     scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x,n=3),
# #                   labels = trans_format("log10", math_format(10^.x))) +
# #     scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=3),
# #                   labels = trans_format("log10", math_format(10^.x))) +
# for (i in col_names) {
#   this_idx <- grep(paste0("^",i,"$"),colnames(circuit_board))
#   dfi <- data.frame(cbind(circuit_board[,i],fit[,this_idx]))
#   dfi <- 10^dfi
#   if (mean(dfi$X1)<1) {
#     dfi <- dfi*1e6
#     }
#   dfi$X1 <- signif(dfi$X1,digits = round_digits_ax)
#   R2val <- paste(round(R2i[this_idx],3))
#   label_R2 =paste("R^2 ==", R2val)
#   ax <- ggplot(dfi,
#                aes(x=X1,
#                    y=X2,
#                    group=factor(X1) #,
#                    #fill=factor(X1) 
#                    ) 
#                ) +
#     geom_abline(intercept = 0, slope = 1) + 
#     geom_boxplot(
#       fill = color_vector[2],
#       alpha = 0.5,
#       outlier.shape=NA
#     ) +
#     scale_x_log10(breaks = unique(dfi$X1), labels = comma) +
#     scale_y_log10(breaks = unique(dfi$X1),labels = comma) +
#     scale_fill_manual(values=color_vector) +
#     labs(
#       title = paste(i),
#       x = paste0("Electrical element (",i,")"),
#       y = paste0("Computed (",i,")")
#     ) +
#     annotate("text", x=c(Inf),y=c(Inf),
#              label=label_R2,
#              parse=TRUE,
#              vjust = 1,
#              hjust = 1) +
#     theme_lewallen()
#   plot_list[[i]] <- ax
#   print(ax)
# }
# # paste(expression(R^2),"=",round(R2i[this_idx],3)),
# p <- plot_grid(plotlist = plot_list,
#           labels="AUTO",
#           ncol = 2,
#           axis = c("lb"))
# title <- ggdraw() + 
#   draw_label("Statistical variation of fit, per circuit element", fontface='bold') +
#   draw_label(paste("n =",nrow(fit)),fontface = 'plain', size = 12, vjust = 2.25)
# plot_grid(title, p, ncol=1, rel_heights=c(0.1, 1)) # rel_heights values control title margins
# file_name <- "R2 analysis.svg"
# ggsave(file_name,
#        width = 16,
#        height = 17,
#        units = "in"
# )
# 
# 




# EXAMPLE EIS FIGURES ----
# make an ID that is unique for each measurment
er$complete_ID <- paste(er$meas_ID,er$meas_idx)
fr$complete_ID <- paste(fr$meas_ID,fr$meas_idx)
#EIS_example <- c("20221018_191151 3","20221110_174503 3","20221117_173404 3")
#EIS_example <- c("20221025_162940 3","20221117_173404 3")
EIS_example <- c("20221209_150838 3")
fit_data <- fr %>%
  filter(complete_ID %in% EIS_example)
figure_width = 12
figure_height = 7



 
 
 
 
 
 
 # FACETED PLOT
# variable.labs <- c("Real","Imaginary", "Zr")
# names(variable.labs) <- c("x", "y", "Zr")
# observation_labels <- list(actual=c("x","y","Zr"),
#                            fit=c("xfitabs","yfitabs","Zrfitabs"))
# observation.labs <- c("Example 1", "Example 2","Example 3")
# names(observation.labs) <- EIS_example
# er %>%
#   filter(complete_ID %in% EIS_example) %>%
#   melt(
#     id.vars = c("complete_ID","f"), 
#     measure.vars = c("x","y","Zr","xfitabs","yfitabs","Zrfitabs")
#     ) %>%
#   left_join(., stack(observation_labels), by = c(variable = "values")) %>%
#   mutate(obs = ind,ind = NULL) %>%
#   mutate(variable = recode(variable, xfitabs = "x", yfitabs = "y", Zrfitabs = "Zr")) %>%
#   dcast(complete_ID+f+variable~obs, value.var = "value") %>%
#   ggplot(aes(x=f)) +
#   geom_line(aes(y=fit)) +
#   geom_point(aes(y=actual),shape=15,size = 0.9,alpha = 0.9) +
#   facet_grid(variable~complete_ID,
#              scales="free_y",
#              labeller = labeller(variable = variable.labs,complete_ID =  observation.labs)) +
#   scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x,n=5),
#                 labels = trans_format("log10", math_format(10^.x))) +
#   labs(
#     y = element_blank(),
#     x = "Frequency (Hz)"
#   ) +
#   theme_lewallen()
# file_name <- paste("Frequency component examples.svg")
# ggsave(file_name,
#        width = 15,
#        height = 18,
#        units = "in")






# CHECK HIGH RSME VALUES
high_RSME <- fr %>%
  filter(rsme_simVSmeasured > 90)
#View(high_RSME)

idx = 1
this_meas_ID <- high_RSME$meas_ID[idx]
this_meas_idx <- high_RSME$meas_idx[idx]
er %>%
  filter(meas_ID %in% high_RSME$meas_ID[idx]) %>%
  filter(meas_idx %in% high_RSME$meas_idx[idx]) %>%
  ggplot()+
  geom_point(aes(x=x, y=y),shape='square',color='black') +
  geom_line(aes(x=xfitabs,y=yfitabs)) +
  labs(
    x = "Real",
    y = "Imaginary",
    title = paste(this_meas_ID,this_meas_idx)
  )

er %>%
  filter(meas_ID %in% high_RSME$meas_ID[idx]) %>%
  filter(meas_idx %in% high_RSME$meas_idx[idx]) %>%
  ggplot()+
  geom_point(aes(x=log10(f), y=Zr),shape='square',color='black') +
  geom_line(aes(x=log10(f),y=Zrfitabs))



# FIND THE RANGE OF RPE  VALUES ----
Ussing_chamber_fit_files <- fr %>%
                              filter(chamber %in% c("Ussing"))
min_Ra <- min(na.omit(Ussing_chamber_fit_files$Rag))
max_Ra <- max(na.omit(Ussing_chamber_fit_files$Rag))
min_Rb <- min(na.omit(Ussing_chamber_fit_files$Rbg))
max_Rb <- max(na.omit(Ussing_chamber_fit_files$Rbg))
min_Rs <- min(na.omit(Ussing_chamber_fit_files$Rsg))
max_Rs <- na.omit(Ussing_chamber_fit_files$Rsg) %>%
  subset(Ussing_chamber_fit_files$Rsg<4.5e6) %>%
  max()
min_Ca <- min(na.omit(Ussing_chamber_fit_files$Cag))
max_Ca <- max(na.omit(Ussing_chamber_fit_files$Cag))
min_Cb <- min(na.omit(Ussing_chamber_fit_files$Cbg))
max_Cb <- max(na.omit(Ussing_chamber_fit_files$Cbg))
min_RsolA <- min(na.omit(Ussing_chamber_fit_files$RsolAg))
max_RsolA <- max(na.omit(Ussing_chamber_fit_files$RsolAg))
min_RsolB <- min(na.omit(Ussing_chamber_fit_files$RsolBg))
max_RsolB <- max(na.omit(Ussing_chamber_fit_files$RsolBg))






# add a row called ID that is a unique N for each measurment.
# This way we can match measred with fit/simulated data
fr_pca <- tibble::rowid_to_column(frc, "ID")
# reshape the data so that we can add the observation tags using melt
fr_pca <- melt(fr_pca,
                        id.vars= c("ID"),
                        measure.vars = c("RsolA","RsolB","Ra","Ca","Rb","Cb","Rs",
                                         "RsolAg","RsolBg","Rag","Cag","Rbg","Cbg","Rsg"))

# measure.vars = c("RsolA","RsolB","Ra","Ca","Rb","Cb","Rs",
#                  "RsolAg","RsolBg","Rag","Cag","Rbg","Cbg","Rsg"))
# measure.vars = c("Ra","Ca","Rb","Cb","Rs",
#                  "Rag","Cag","Rbg","Cbg","Rsg"))

# add an observation tag for each measurement
observation_labels <- list(actual=c("RsolA","RsolB","Ra","Ca","Rb","Cb","Rs"),
                           fit=c("RsolAg","RsolBg","Rag","Cag","Rbg","Cbg","Rsg"))
fr_pca <- fr_pca %>% 
  left_join(., stack(observation_labels), by = c(variable = "values"))

# rename the column "ind" to "obs" to reflect that it is an observation tag
names(fr_pca)[names(fr_pca) == "ind"] <- "obs"

# rename the observations so that they properly match what they were measuring
fr_pca$variable <-  gsub('RsolAg','RsolA',
                  gsub('RsolBg','RsolB',
                  gsub('Rag','Ra',
                  gsub('Cag','Ca',
                  gsub('Rbg','Rb',
                  gsub('Cbg','Cb',
                  gsub('Rsg','Rs',fr_pca$variable)))))))


fr_pca$variableID <- paste(fr_pca$variable,fr_pca$ID)
labels_to_keep <- c("variableID","variable","ID")
fr_pca_labels <- fr_pca[,labels_to_keep]

fr_pca <- dcast(fr_pca, variableID~obs, value.var = "value")
fr_pca <- left_join(fr_pca,fr_pca_labels,by="variableID")
fr_pca <- na.omit(fr_pca)

# Describe ----
# PCA analysis
reqd <- names(observation_labels)
fit_analysis_pca <- log10(fr_pca[,reqd])
fit_analysis.pca <- prcomp(fit_analysis_pca, center = TRUE,scale. = TRUE)
summary(fit_analysis.pca)

# Visualize ----
text_size = 12
ggbiplot(fit_analysis.pca,
         choice = c(1,2),
         obs.scale = 1, var.scape = 1,
         varname.size=text_size/3,
         labels = fr_pca$ID,
         groups = fr_pca$variable,
         alpha = 0.2,
         circle = TRUE) +
  theme_tufte() +
  labs(title = "PCA Analysis",
       color = "Variable") +
  theme(text = element_text(size = text_size))

# frc %>%
#   ggplot(aes(x=round(log10(Ra),1),y=log10(Rag),color=Ca))+
#   geom_abline(slope=1, intercept=0)+
#   geom_jitter() +
#   theme_wsj(color='gray') +
#   labs(title = 'Error visualization',
#        subtitle = 'How well can we guess Ra?',
#        x = expression(log[10](R_a) (Ohms*.*cm^2)))

# Analyze ----



# USSING CHAMBER EXPERIMENTS

