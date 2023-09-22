# EISandPipette.R ERROR analysis
# study the types and quality of error in the fitting algorithm
meas_ID_remove <- LT$meas_ID[981:nrow(LT)]
df1 <- fr %>% 
  filter(circuitValuesKnown %in% c(1)) %>%
  filter( !( meas_ID %in% c(meas_ID_remove) ) )

# grab the fit and membrane parameters from the data frame
p <- df1[grepl("p_",colnames(df1))]
pg_absnone <- df1[grepl("pg_absnone_",colnames(df1))]
pg_absZr <- df1[grepl("pg_absZr_",colnames(df1))]
pg_absEr <- df1[grepl("pg_absEr_",colnames(df1))]

# # calculate combined terms
# pg_absZr$pg_absZr_Rapical <- pg_absZr$pg_absZr_1+pg_absZr$pg_absZr_3
# pg_absZr$pg_absZr_Rbasal <- pg_absZr$pg_absZr_2+pg_absZr$pg_absZr_5
# pg_absEr$pg_absEr_Rapical <- pg_absEr$pg_absEr_1+pg_absEr$pg_absEr_3
# pg_absEr$pg_absEr_Rbasal <- pg_absEr$pg_absEr_2+pg_absEr$pg_absEr_5

# delta p
dpg_absnone <- pg_absnone-p
dpg_absZr <- pg_absZr-p
dpg_absEr <- pg_absEr-p

# normalize delta p
dpg_absnone_error <- dpg_absnone/p
dpg_absZr_error <- dpg_absZr/p
dpg_absEr_error <- dpg_absEr/p

# median errors for absZr (using the pipette)
data <- dpg_absZr_error*100
rep_str = c('pg_absZr'='',
            '_1'='RsolA',
            '_2'='RsolB',
            '_3'='Ra',
            '_4'='Ca',
            '_5'='Rb',
            '_6'='Cb',
            '_7'='Rs',
            '_Rapical' = 'Rapic_all',
            '_Rbasal' = 'Rbaso_all'
)
colnames(data) <- str_replace_all(colnames(data), rep_str)
mdf_absZr <- melt(
  data,
  varnames = names(dimnames(data)),
  na.rm = FALSE,
  as.is = FALSE,
  value.name = "value" 
)

# median errors for absEr (using the pipette)
data <- dpg_absEr_error*100
rep_str = c('pg_absEr'='',
            '_1'='RsolA',
            '_2'='RsolB',
            '_3'='Ra',
            '_4'='Ca',
            '_5'='Rb',
            '_6'='Cb',
            '_7'='Rs',
            '_a' = 'Rapic_all',
            '_b' = 'Rbaso_all'
)
colnames(data) <- str_replace_all(colnames(data), rep_str)
mdf_absEr <- melt(
  data,
  varnames = names(dimnames(data)),
  na.rm = FALSE,
  as.is = FALSE,
  value.name = "value" 
)

# median errors for absnone (without the pipette)
data <- dpg_absnone_error*100
rep_str = c('pg_absnone'='',
            '_1'='RsolA',
            '_2'='RsolB',
            '_3'='Ra',
            '_4'='Ca',
            '_5'='Rb',
            '_6'='Cb',
            '_7'='Rs',
            '_a' = 'Rapic_all',
            '_b' = 'Rbaso_all'
)
colnames(data) <- str_replace_all(colnames(data), rep_str)
mdf_absnone <- melt(
  data,
  varnames = names(dimnames(data)),
  na.rm = FALSE,
  as.is = FALSE,
  value.name = "value" 
)

# rename the variables for legibility
mdf_absZr <- mdf_absZr %>% filter(variable %in% c("RsolA","RsolB","Ra","Ca","Rb","Cb","Rs"))
mdf_absEr <- mdf_absEr %>% filter(variable %in% c("RsolA","RsolB","Ra","Ca","Rb","Cb","Rs"))
mdf_absnone <- mdf_absnone %>% filter(variable %in% c("RsolA","RsolB","Ra","Ca","Rb","Cb","Rs"))

# add the method name column
mdf_absZr$method <- rep("Membrane Ratio",nrow(mdf_absZr))
mdf_absEr$method <- rep("Electrode Ratio",nrow(mdf_absEr))
mdf_absnone$method <- rep("Without pipette",nrow(mdf_absnone))

# merge the data together in one large data frame
mdf2 <- rbind(mdf_absZr,mdf_absEr,mdf_absnone)

# define the color palette
library(RColorBrewer)
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#999999")
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
color_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
color_vector = c(cb_palette,color_vector) # add the colorblind palette if desired

# plot the distribution of data
ggplot(mdf2,aes(x=value,fill=factor(method,levels=c("Electrode Ratio","Membrane Ratio","Without pipette")))) +
  #geom_jitter(shape=16, position=position_jitter(0.1),size = 1,alpha = 0.5) +
  #geom_violin(alpha = 0.5) +
  geom_density(alpha = 0.7) +
  facet_grid(variable~.) +
  #coord_flip() +
  labs(
    x=element_blank(),
    y=element_blank(),
    fill=element_blank()
  ) +
  xlim(-100,250) +
  scale_fill_manual(values = color_vector) +
  theme_stata(scheme = "s1color",base_size = 12)+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_line( size=.05, color="grey" ),
        panel.grid.major.x = element_line( size=.05, color="grey" ))
file_name <- paste("All fit errors (cropped).svg", sep = "")
ggsave(file_name,
       path = SAVE_path,
       plot = last_plot(),
       width = figure_width_large, height = figure_height_large,
       units = "in",
       limitsize = TRUE) 

# show the stats in the terminal 
library(matrixStats)
tmp <- colMedians((as.matrix(dpg_absZr_error))*100)
summary(dpg_absZr_error*100)

tmp <- colMedians((as.matrix(dpg_absEr_error))*100)
summary(dpg_absEr_error*100)





# # logarithmic error
# dpg_absnone_log10 <- log10(abs(dpg_absnone_error))
# dpg_absZr_log10 <- log10(abs(dpg_absZr_error))
# 
# # mean log error
# colMeans(dpg_absnone_log10)
# colMeans(dpg_absZr_log10)
# 
# 
# # log10 error plots ----
# ylimits <- c(0,400) # force the limits on all plots to have same y axis
# xlimits <- c(-5,5) # force the limits on all plots to have the same x axis
# 
# data <- dpg_absnone_log10
# rep_str = c('pg_absnone'='',
#             '_1'='RsolA',
#             '_2'='RsolB',
#             '_3'='Ra',
#             '_4'='Ca',
#             '_5'='Rb',
#             '_6'='Cb',
#             '_7'='Rs'
#             )
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
#                  bins = 50) +
#   labs(
#     x=expression(Relative~error~(epsilon)),
#     y="Count"
#   ) +
#   #ylim(ylimits) +
#   #xlim(xlimits) +
#   facet_grid(~variable) +
#   theme_lewallen()
# file_name <- paste("mean log error absnone",".svg", sep = "")
# ggsave(file_name,
#        path = SAVE_path,
#        plot = last_plot(),
#        width = figure_width_large*1.4, 
#        height = figure_height_medium,
#        units = "in",
#        limitsize = TRUE) 
# 
# 
# data <- dpg_absZr_log10
# rep_str = c('pg_absZr'='',
#             '_1'='RsolA',
#             '_2'='RsolB',
#             '_3'='Ra',
#             '_4'='Ca',
#             '_5'='Rb',
#             '_6'='Cb',
#             '_7'='Rs'
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
#                  bins = 50) +
#   labs(
#     x=expression(Relative~error~(epsilon)),
#     y="Count"
#   ) +
#   #ylim(ylimits) +
#   #xlim(xlimits) +
#   facet_grid(~variable) +
#   theme_lewallen()
# file_name <- paste("mean log error absZr",".svg", sep = "")
# ggsave(file_name,
#        path = SAVE_path,
#        plot = last_plot(),
#        width = figure_width_large*1.4, 
#        height = figure_height_medium,
#        units = "in",
#        limitsize = TRUE) 
# 
# 
# # remove boundary fits ----
# pgmin_absZr <- df1[grepl("pgmin_absZr_",colnames(df1))]
# pgmax_absZr <- df1[grepl("pgmax_absZr_",colnames(df1))]
# tmp <- apply(pg_absZr, 1, function(x) any(x<= pgmin_absZr*1.01| x >= pgmax_absZr*0.99))
# 
# 
# 
# df_filt <- sum(rowSums(((pg_absZr <= pgmin_absZr*1.01 | pg_absZr >= pgmax_absZr*0.99)*1)[,-c(1,2)]))
# 

# 
# 
# # figure size
# figure_width = 2.5
# figure_height = 5.5
# 
# fit_results_error <- fit_results_clean
#   
# # Ra
# fit_results_error %>%
#   ggplot() +
#   geom_histogram(aes(x=Rag_epsilon),
#                  bins = 50) +
#   labs(
#     x=expression(Relative~error~(epsilon)),
#     y="Count"
#   ) +
#   theme_lewallen()
# file_name <- paste("Rag error.svg")
# ggsave(file_name,
#        width = figure_width,
#        height = figure_height,
#        units = "in")
# 
# # Rb
# fit_results_error %>%
#   ggplot() +
#   geom_histogram(aes(x=Rbg_epsilon),
#                  bins = 50) +
#   labs(
#     x=expression(Relative~error~(epsilon)),
#     y="Count"
#   ) +
#   theme_lewallen()
# file_name <- paste("Rbg error.svg")
# ggsave(file_name,
#        width = figure_width,
#        height = figure_height,
#        units = "in")
# 
# # Rs
# fit_results_error %>%
#   ggplot() +
#   geom_histogram(aes(x=Rsg_epsilon),
#                  bins = 50) +
#   labs(
#     x=expression(Relative~error~(epsilon)),
#     y="Count"
#   ) +
#   theme_lewallen()
# file_name <- paste("Rsg error.svg")
# ggsave(file_name,
#        width = figure_width,
#        height = figure_height,
#        units = "in")
# 
# # Ca
# fit_results_error %>%
#   ggplot() +
#   geom_histogram(aes(x=Cag_epsilon),
#                  bins = 50) +
#   labs(
#     x=expression(Relative~error~(epsilon)),
#     y="Count"
#   ) +
#   theme_lewallen()
# file_name <- paste("Cag error.svg")
# ggsave(file_name,
#        width = figure_width,
#        height = figure_height,
#        units = "in")
# 
# # Cb
# fit_results_error %>%
#   ggplot() +
#   geom_histogram(aes(x=Cbg_epsilon),
#                  bins = 50) +
#   labs(
#     x=expression(Relative~error~(epsilon)),
#     y="Count"
#   ) +
#   theme_lewallen()
# file_name <- paste("Cbg error.svg")
# ggsave(file_name,
#        width = figure_width,
#        height = figure_height,
#        units = "in")
# 
# # RsolA
# fit_results_error %>%
#   ggplot() +
#   geom_histogram(aes(x=RsolAg_epsilon),
#                  bins = 50) +
#   labs(
#     x=expression(Relative~error~(epsilon)),
#     y="Count"
#   ) +
#   theme_lewallen()
# file_name <- paste("RsolAg error.svg")
# ggsave(file_name,
#        width = figure_width,
#        height = figure_height,
#        units = "in")
# 
# # RsolB
# fit_results_error %>%
#   ggplot() +
#   geom_histogram(aes(x=RsolBg_epsilon),
#                  bins = 50) +
#   labs(
#     x=expression(Relative~error~(epsilon)),
#     y="Count"
#   ) +
#   theme_lewallen()
# file_name <- paste("RsolBg error.svg")
# ggsave(file_name,
#        width = figure_width,
#        height = figure_height,
#        units = "in")
# 
# 
# # build a table
# mean(fit_results_error$TER_epsilon)
# mean(fit_results_error$TEC_epsilon)
# mean(fit_results_error$GsGt_epsilon)
# mean(fit_results_error$VDR_epsilon)
# mean(fit_results_error$Rag_epsilon)
# mean(fit_results_error$Rbg_epsilon)
# mean(fit_results_error$Rsg_epsilon)
# mean(fit_results_error$Cag_epsilon)
# mean(fit_results_error$Cbg_epsilon)
# mean(fit_results_error$Aold_epsilon)
# mean(fit_results_error$TERold_epsilon)
# mean(fit_results_error$RsolAg_epsilon)
# mean(fit_results_error$RsolBg_epsilon)
# 
# sd(fit_results_error$Rag_epsilon)
# sd(fit_results_error$Rbg_epsilon)
# sd(fit_results_error$Rsg_epsilon)
# sd(fit_results_error$Cag_epsilon)
# sd(fit_results_error$Cbg_epsilon)
# 
# median(fit_results_error$RsolAg_deltaR)
# median(fit_results_error$RsolBg_deltaR)
# median(fit_results_error$Rag_deltaR)
# median(fit_results_error$Rbg_deltaR)
# median(fit_results_error$Rsg_deltaR)
# median(fit_results_error$Cag_deltaR)
# median(fit_results_error$Cbg_deltaR)
# 
#        