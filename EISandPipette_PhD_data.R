# USSING CHAMBER EXPERIMENTS ----

# define the error percentage for each parameter
# median error array: 6.415652 12.096486 -0.339046  6.186783 -7.417236  7.239384 -1.731100 -1.584547  5.496925  7.187270  1.308234 -3.977258
Rae = abs(-0.339046/100)
Rbe = abs(-7.417236/100)
Rse = abs(-1.731100/100)
Cae = abs(6.186783/100)
Cbe = abs(7.239384/100)
TERe = abs(-1.584547/100)
TECe = abs(5.496925/100)
VDRe = abs(7.187270/100)
VDRolde = abs(4.845445/100)
TERolde = abs(10.70025/100)
GsGte = abs(-1.634807e-07/100)

# sample order for paper
library(scales)
library(patchwork)
dfs <- fr %>% 
  filter(chamber %in% c('Ussing'))

# correct for the cross sectional area and rename
dfs$RsolAg<-dfs$pg_absZr_1*dfs$cross_section
dfs$RsolBg<-dfs$pg_absZr_2*dfs$cross_section
dfs$Rag<-dfs$pg_absZr_3*dfs$cross_section
dfs$Cag<-dfs$pg_absZr_4/dfs$cross_section
dfs$Rbg<-dfs$pg_absZr_5*dfs$cross_section
dfs$Cbg<-dfs$pg_absZr_6/dfs$cross_section
dfs$Rsg<-dfs$pg_absZr_7*dfs$cross_section

dfs$Rsolg<-dfs$pg_12_1*dfs$cross_section
dfs$R1g<-dfs$pg_12_2*dfs$cross_section
dfs$C1g<-dfs$pg_12_3/dfs$cross_section
dfs$R2g<-dfs$pg_12_4*dfs$cross_section
dfs$C2g<-dfs$pg_12_5/dfs$cross_section

dfs$Aold <-dfs$pg_VDRold
dfs$TERabs<-dfs$pg_absZr_TER*dfs$cross_section
dfs$TECabs<-dfs$pg_absZr_TEC/dfs$cross_section
dfs$TERold<-dfs$pg_TERold*dfs$cross_section
dfs$GsGtabs<-dfs$pg_absZr_TER/dfs$pg_absZr_7
dfs$taua<-dfs$Rag*dfs$Cag
dfs$taub<-dfs$Rbg*dfs$Cbg
dfs$tau1<-dfs$R1g*dfs$C1g
dfs$tau2<-dfs$R2g*dfs$C2g

# calculate the percentage of Rs as a function of TERabs
Gs_predicted <- 1/dfs$TERabs*0.71
dfs$Rs_predicted <- 1/Gs_predicted
dfs$tau_r_abs <- dfs$taua/dfs$taub
dfs$tau_r_12 <- dfs$tau1/dfs$tau2
for (i in 1:nrow(dfs)) {
  if (dfs$tau_r_12[i] < 1) {
    dfs$tau_r_12[i] = 1/dfs$tau_r_12[i]
    tR1 = dfs$R1g[i]
    tC1 = dfs$C1g[i]
    tR2 = dfs$R2g[i]
    tC2 = dfs$C2g[i]
    ttau1 = dfs$tau1[i]
    ttau2 = dfs$tau2[i]
    dfs$R1g[i] = tR2
    dfs$C1g[i] = tC2
    dfs$R2g[i] = tR1
    dfs$C2g[i] = tC1
    dfs$tau1[i] = ttau2
    dfs$tau2[i] = ttau1
  }
}

# preview available files
meas_in_file <- unique(dfs$meas_ID)

# SELECT DATA for plots
# selected_files = unique(dfs$meas_ID) # use this one to preview all files

# have to uncomment labeller section of plot to get files renamed
# selected_files = c('20211021_183356','20211021_160529') 
# plot_names = c(
#   '20211021_183356' = "Sample 1",
#   '20211021_160529' = "Sample 3")

# all Z8
# selected_files = c('20211021_183356','20211022_171622','20211021_160529') 
# plot_names = c(
#    '20211021_183356' = "Sample 1",
#    '20211022_171622' = "Sample 2",
#    '20211021_160529' = "Sample 3")

#sample 1 and 2 are Z8 and 3 is dominik AMDCD
selected_files = c('20211021_183356','20211022_171622','20211014_165514')
plot_names = c(
  '20211021_183356' = "Sample 1",
  '20211022_171622' = "Sample 2",
  '20211014_165514' = "Sample 3")

selected_dfs <- subset(dfs,meas_ID == selected_files[3])
ggplot(selected_dfs,aes(x=meas_idx,y=Aold)) +
  geom_line()

# all samples
# selected_files = unique(dfs$meas_ID)
# plot_names = c(
#   '20211013_173420' = "D3C 20211013_173420",
#   '20211021_160529' = "Z8 20211021_160529",
#   '20211101_150420' = "Z8 20211101_150420",
#   '20211019_173703' = "Z8 20211019_173703",
#   '20211031_174839' = "Z8 20211031_174839",
#   '20211014_165514' = "AMDCD 20211014_165514",
#   '20211021_183356' = "Z8 20211021_183356",
#   '20211022_171622' = "Z8 20211022_171622")

# all samples
# selected_files = unique(dfs$meas_ID)
# plot_names = c(
#   '20211013_173420' = "D3C 20211013_173420",
#   '20211021_160529' = "Z8 20211021_160529",
#   '20211101_150420' = "Z8 20211101_150420",
#   '20211019_173703' = "Z8 20211019_173703",
#   '20211031_174839' = "Z8 20211031_174839",
#   '20211014_165514' = "AMDCD 20211014_165514",
#   '20211021_183356' = "Z8 20211021_183356",
#   '20211022_171622' = "Z8 20211022_171622")

# plot settings
# time_limits = c(0,69)
time_limits = c(0,69)
point_size = 2
point_shape = 22
nSamples = length(selected_files)
# figure_width = 2.25*nSamples+1/8
# figure_height = 2.25
# line_color = color_vector[8]
line_color = 'black'
line_color2 = color_vector[1]
box_color = color_vector[5]
box_alpha = 0.5
ribbon_color = 'lightgrey'
ribbon_alpha = 1  
# line_alpha = 0.3
line_alpha = 0.99
# line_size <- 0.7 * (figure_width / 6)
line_size <- 0.25
vdr_scale <- 70
# init cow plot lists
plot_Ra      = list()
plot_Rb      = list()
plot_Rs      = list()
plot_Ca      = list()
plot_Cb      = list()
plot_VDRold  = list()
plot_tau_r_abs = list()
plot_tau_r_12 = list()
plot_TERold = list()
plot_TERabs  = list()
plot_Ct      = list()
plot_GsGt    = list()
plot_TEP     = list()
plot_Va      = list()
plot_Vb      = list()
plot_taua    = list()
plot_taub    = list()
plot_tau1    = list()
plot_tau2    = list()

all_dfs <- dfs %>% filter(meas_ID %in% selected_files)
mean(all_dfs$RsolAg)
sd(all_dfs$RsolAg)
length(all_dfs$RsolBg)
mean(all_dfs$RsolBg)
sd(all_dfs$RsolBg)
length(all_dfs$RsolBg)
for (i in 1:length(selected_files)) {
  # subset the EIS data for the selected file, only
  this_dfs <- dfs[dfs$meas_ID==selected_files[i],]
  
  
  # use comment file to get start and stop of ATP
  COMMENT_path <- file.path("G:","My Drive", "Colby", "EIS and Pipette","data",
                            selected_files[i],
                            paste0(selected_files[i],"_comments.txt"),
                            fsep = .Platform$file.sep)
  comments = read.csv(COMMENT_path, sep = '\t')
  comment_idx = c(1,2) # for all files, make sure the 1st and second comment correspond to atp
  
  # use NOVA file to get absolute start time of the recording
  NOVA_path <- file.path("G:","My Drive", "Colby", "EIS and Pipette","data",
                            selected_files[i],
                            paste0(selected_files[i],"_NOVAdata.txt"),
                            fsep = .Platform$file.sep)
  NOVA <- read.csv(NOVA_path,sep='\t')
  NOVA_t = lapply(list(NOVA[,10]), function(z){ z[!is.na(z) & z != "" & z != "meas start time (s)"]})
  NOVA_t = as.numeric(NOVA_t[[1]])
  # in case some nova is missing, add extra final times to the recording
  if (length(NOVA_t) != nrow(this_dfs)) {
    next
  }
  t0 = NOVA_t[1]
  
  # load the voltage data such as TEP
  voltages_path <- file.path("G:","My Drive", "Colby", "EIS and Pipette","data",
                             selected_files[i],
                             paste0(selected_files[i],"_data.txt"),
                             fsep = .Platform$file.sep)
  voltages <- read.csv(voltages_path,sep="\t")
  voltages$t <- (voltages$absolute.time.s.-t0)/60
  voltages$Va <- voltages$Va..mV.
  voltages$Vb <- voltages$Vb..mV.
  voltages$TEP <- voltages$TEP..mV.
  voltages$TTL <- voltages$TTL.signal > 0.5
  voltages$meas_ID <- rep(this_dfs$meas_ID[1],nrow(voltages))
  
  # replace voltages with NaN if EIS was recording
  # voltages$TEP[voltages$TTL==TRUE] <- NaN
  # voltages$Va[voltages$TTL==TRUE] <- NaN
  # voltages$Vb[voltages$TTL==TRUE] <- NaN
  voltages <- subset(voltages,TTL==FALSE)
  
  # synchronize all times to the absolute time files
  this_dfs$t = (NOVA_t-t0)/60
  comment_times = (comments[comment_idx,1] - t0)/60
  
  ###################### CROP SPECIFICALLY THE AMDCD FILE WHERE Ba2+ WAS APPLIED
  ############ if (selected_files[i]=='20211014_165514') {time_limits = c(0,55)}
  if (selected_files[i]=='20211014_165514') {
    idx_max = which.min(abs(this_dfs$t-55))
    this_dfs = this_dfs[1:idx_max,]
    
    idx_max = which.min(abs(voltages$t-55))
    voltages = voltages[1:idx_max,]
  }
  
  # Rag
  # y axis label
  if (i==1){
    y_text = expression(R[a]~(Omega*cm^2))
  } else {
    y_text = element_blank()
  }
  
  plot_Ra[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_ribbon(aes(x=t,
                    ymin=Rag-Rag*Rae,
                    ymax=Rag+Rag*Rae),
                alpha = ribbon_alpha,
                color = "transparent",
                fill = ribbon_color
                
    ) +
    geom_line(aes(x=t,y=Rag),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=Rag),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    xlab(element_blank()) +
    ylab(y_text) +
    labs(shape=element_blank()) + 
    scale_y_continuous(
      limits = c(0,8000),
      breaks = seq(0,8000,by=1250),
      labels = scales::comma_format(big.mark = ',',decimal.mark = '.')) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    # facet_wrap(~meas_ID,
    #            #scales = "free"#,
    #            labeller = as_labeller(plot_names[i])
    #            ) + # adds column labels to plot
    theme_lewallen() + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
  
  
  
  
  # Rbg
  # y axis label
  if (i==1){
    y_text = expression(R[b]~(Omega*cm^2))
  } else {
    y_text = element_blank()
  }
  
  plot_Rb[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_ribbon(aes(x=t,
                    ymin=Rbg-Rbg*Rbe,
                    ymax=Rbg+Rbg*Rbe),
                alpha = ribbon_alpha,
                color = "transparent",
                fill = ribbon_color
                
    ) +
    geom_line(aes(x=t,y=Rbg),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=Rbg),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    #xlab("Time (min.)") +
    xlab(element_blank())+
    ylab(y_text) +
    labs(shape=element_blank()) + 
    scale_y_continuous(
      limits = c(0,2000),
      breaks = seq(0,2000,by=400),
      labels = scales::comma_format(big.mark = ',',decimal.mark = '.')) +
    scale_x_continuous(
      limits = time_limits) +
    theme_lewallen() + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank())
  
  
  
  
  # Rsg
  # y axis label
  if (i==1){
    y_text = expression(R[s]~(Omega*cm^2))
  } else {
    y_text = element_blank()
  }
  plot_Rs[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_ribbon(aes(x=t,
                    ymin=Rsg-Rsg*Rse,
                    ymax=Rsg+Rsg*Rse),
                alpha = ribbon_alpha,
                color = "transparent",
                fill = ribbon_color
                
    ) +
    geom_line(aes(x=t,y=Rsg),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=Rsg),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    xlab("Time (min.)") +
    #xlab(element_blank()) +
    ylab(y_text) +
    labs(shape=element_blank()) + 
    scale_y_continuous(
      limits = c(0,750),
      breaks = seq(0,750,by=150),
      labels = scales::comma_format(big.mark = ',',decimal.mark = '.')) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    theme_lewallen() #+ theme(
      #axis.text.x = element_blank(),
      #axis.ticks.x = element_blank())
  
  
  # GsGt percent
  # y axis label
  if (i==1){
    y_text = expression(G[s]/G[t]~(percent))
  } else {
    y_text = element_blank()
  }
  plot_GsGt[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_ribbon(aes(x=t,
                    ymin=(GsGtabs-GsGtabs*GsGte)*100,
                    ymax=(GsGtabs+GsGtabs*GsGte)*100),
                alpha = ribbon_alpha,
                color = "transparent",
                fill = ribbon_color
                
    ) +
    geom_line(aes(x=t,y=GsGtabs*100),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=GsGtabs*100),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    # xlab("Time (min.)") +
    xlab(element_blank()) +
    ylab(y_text) +
    labs(shape=element_blank()) + 
    scale_y_continuous(
      limits = c(50, 100) #,
      #labels = label_number(accuracy = 1)
      ) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    facet_wrap(~meas_ID ,
               # scales = "free"#,
               labeller = as_labeller(plot_names[i])
    ) + # adds column labels to plot
    theme_lewallen() + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank())
  
  
  
  
  # Cag
  # y axis label
  if (i==1){
    y_text = expression(C[a]~(mu*F/cm^2))
  } else {
    y_text = element_blank()
  }
  
  plot_Ca[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_ribbon(aes(x=t,
                    ymin=(Cag-Cag*Cae)*1e6,
                    ymax=(Cag+Cag*Cae)*1e6),
                alpha = ribbon_alpha,
                color = "transparent",
                fill = ribbon_color
                
    ) +
    geom_line(aes(x=t,y=Cag*1e6),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=Cag*1e6),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    xlab(element_blank()) +
    ylab(y_text) +
    labs(shape=element_blank()) + 
    scale_y_continuous(
      limits = c(11,25),
      labels = scales::comma_format(big.mark = ',',decimal.mark = '.')) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    facet_wrap(~meas_ID ,
               # scales = "free"#,
               labeller = as_labeller(plot_names[i])
               ) + # adds column labels to plot
    theme_lewallen() + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
  
  
  
  
  # Cbg
  # y axis label
  if (i==1){
    y_text = expression(C[b]~(mu*F/cm^2))
  } else {
    y_text = element_blank()
  }
  
  plot_Cb[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_ribbon(aes(x=t,
                    ymin=(Cbg-Cbg*Cbe)*1e6,
                    ymax=(Cbg+Cbg*Cbe)*1e6),
                alpha = ribbon_alpha,
                color = "transparent",
                fill = ribbon_color
                
    ) +
    geom_line(aes(x=t,y=Cbg*1e6),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=Cbg*1e6),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    xlab("Time (min.)") +
    ylab(y_text) +
    labs(shape=element_blank()) + 
    scale_y_continuous(
      limits = c(6,11),
      labels = scales::comma_format(big.mark = ',',decimal.mark = '.')) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    # facet_wrap(~meas_ID,
    # #           scales = "free",
    #            labeller = as_labeller(plot_names[i])) + # adds column labels to plot
    theme_lewallen() #+ theme(
  #axis.text.x = element_blank(),
  #axis.ticks.x = element_blank())
  
  
  
  
  
  # Aold
  # y axis label
  if (i==1){
    y_text = expression(VDR[0.5~Hz]~(R[a]/R[b]))
  } else {
    y_text = element_blank()
  }
  
  plot_VDRold[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_ribbon(aes(x=t,
                    ymin=(Aold-Aold*VDRolde),
                    ymax=(Aold+Aold*VDRolde)),
                alpha = ribbon_alpha,
                color = "transparent",
                fill = ribbon_color
                
    ) +
    geom_line(aes(x=t,y=Aold),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=Aold),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    ylab(y_text) +
    labs(shape=element_blank()) + 
    scale_y_continuous(
      limits = c(0,7),
      breaks = seq(0,7,by=2.3),
      labels = label_number(accuracy = 0.1)
      ) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    # facet_wrap(~meas_ID,
    # #           scales = "free",
    #            labeller = as_labeller(plot_names[i])) + # adds column labels to plot
    # xlab(element_blank()) +
    xlab('Time (min).') +
    theme_lewallen() #+ theme(
      # axis.text.x = element_blank(),
      # axis.ticks.x = element_blank())
  
  
  
  
  # tau_ratio_abs
  # y axis label
  if (i==1){
    y_text = expression(Time~constant~ratio~(tau[a]/tau[b]))
  } else {
    y_text = element_blank()
  }
  
  plot_tau_r_abs[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_line(aes(x=t,y=tau_r_abs),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=tau_r_abs),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    xlab(element_blank()) +
    ylab(y_text) +
    labs(shape=element_blank()) +
    scale_y_continuous(
      limits = c(1,15),
      labels = label_number(accuracy = 0.1)) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    #ylim(c(1,9)) +
    theme_lewallen() + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
  
  
  # tau_ratio_12
  # y axis label
  if (i==1){
    y_text = expression(Time~constant~ratio~(tau[1]/tau[2]))
  } else {
    y_text = element_blank()
  }
  
  plot_tau_r_12[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_line(aes(x=t,y=tau_r_12),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=tau_r_12),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    xlab(element_blank()) +
    ylab(y_text) +
    labs(shape=element_blank()) +
    scale_y_continuous(
      limits = c(1,10),
      labels = label_number(accuracy = 0.1)) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    #ylim(c(1,9)) +
    theme_lewallen() + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
  
  
  
  
  # Ct
  # y axis label
  if (i==1){
    y_text = expression(Tissue~capacitance~(mu*F/cm^2))
  } else {
    y_text = element_blank()
  }
  
  plot_Ct[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_ribbon(aes(x=t,
                    ymin=(TECabs-TECabs*TECe)*1e6,
                    ymax=(TECabs+TECabs*TECe)*1e6),
                alpha = ribbon_alpha,
                color = "transparent",
                fill = ribbon_color
                
    ) +
    geom_line(aes(x=t,y=TECabs*1e6),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=TECabs*1e6),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    ylab(y_text) +
    labs(shape=element_blank()) +
    scale_y_continuous(
      #limits = c(4.7,6),
      labels = label_number(accuracy = 0.1)) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    xlab("Time (min).") +
    theme_lewallen() #+ theme(
  #axis.text.x = element_blank(),
  #axis.ticks.x = element_blank())
  
  
  
  # TEP
  # y axis label
  if (i==1){
    y_text = expression(TEP~(mV))
  } else {
    y_text = element_blank()
  }
  
  plot_TEP[[i]] <- ggplot(voltages) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_line(aes(x=t,y=TEP),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size+0.5) +
    ylab(y_text) +
    labs(shape=element_blank()) +
    scale_y_continuous(
      limits = c(-2,4),
      breaks = seq(-2,4,by=1),
      labels = label_number(accuracy = 1)) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    #xlab("Time (min).") +
    xlab(element_blank()) +
    facet_wrap(~meas_ID,
               #scales = "free"#,
               labeller = as_labeller(plot_names[i])
    ) + # adds column labels to plot
    theme_lewallen() + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank())
  
  
  
  
  # Va
  # y axis label
  if (i==1){
    y_text = expression(Membrane~voltages~(mV))
  } else {
    y_text = element_blank()
  }
  
  plot_Va[[i]] <- ggplot(voltages) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_line(aes(x=t,y=Va),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size+0.5) +
    geom_line(aes(x=t,y=Vb),
              color=color_vector[1],
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size+0.5) +
    ylab(y_text) +
    labs(shape=element_blank()) +
    scale_y_continuous(
      limits = c(-60,-20),
      labels = label_number(accuracy = 1)) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    #xlab("Time (min).") +
    xlab(element_blank()) +
    # facet_wrap(~meas_ID,
    #            #scales = "free"#,
    #            labeller = as_labeller(plot_names[i])
    # ) + # adds column labels to plot
    theme_lewallen() #+ theme(
      #axis.text.x = element_blank(),
      #axis.ticks.x = element_blank())
  
  
  
  # Vb
  # y axis label
  if (i==1){
    y_text = expression(V[b]~(mV))
  } else {
    y_text = element_blank()
  }
  
  plot_Vb[[i]] <- ggplot(voltages) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_line(aes(x=t,y=Vb),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size+1) +
    ylab(y_text) +
    labs(shape=element_blank()) +
    scale_y_continuous(
      limits = c(-60,-20),
      labels = label_number(accuracy = 1)) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    xlab("Time (min).") +
    # xlab(element_blank()) +
    # facet_wrap(~meas_ID,
    #            #scales = "free"#,
    #            labeller = as_labeller(plot_names[i])
    # ) + # adds column labels to plot
    theme_lewallen() #+ theme(
      # axis.text.x = element_blank(),
      # axis.ticks.x = element_blank())
  
  
  
  
  # TERold
  # y axis label
  if (i==1){
    y_text = expression(TER[0.5~Hz]~(Omega*cm^2))
  } else {
    y_text = element_blank()
  }

  plot_TERold[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_ribbon(aes(x=t,
                    ymin=(TERold-TERold*TERolde),
                    ymax=(TERold+TERold*TERolde)),
                alpha = ribbon_alpha,
                color = "transparent",
                fill = ribbon_color
                
    ) +
    geom_line(aes(x=t,y=TERold),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    # geom_line(aes(x=t,y=TERold),
    #           color=line_color,
    #           alpha=line_alpha,
    #           show.legend = FALSE,
    #           size = line_size) +
    geom_point(aes(x=t,y=TERold),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    scale_y_continuous(
      name = y_text,
      limits = c(0,720),
      breaks = seq(0,720,by=120)
    ) +
    labs(shape=element_blank()) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    facet_wrap(~meas_ID ,
               #scales = "free",
               labeller = as_labeller(plot_names[i])
               ) + # adds column labels to plot
    xlab(element_blank()) +
    # xlab('Time (min).') +
    theme_lewallen() + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
  
  
  # TERabs
  # y axis label
  if (i==1){
    y_text = expression(TER~(Omega*cm^2))
  } else {
    y_text = element_blank()
  }
  
  plot_TERabs[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    geom_ribbon(aes(x=t,
                    ymin=(TERabs-TERabs*TERe),
                    ymax=(TERabs+TERabs*TERe)),
                alpha = ribbon_alpha,
                color = "transparent",
                fill = ribbon_color
                
    ) +
    geom_line(aes(x=t,y=TERabs),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=TERabs),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    scale_y_continuous(
      name = y_text,
      limits = c(0,720),
      breaks = seq(0,720,by=120)
    ) +
    labs(shape=element_blank()) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    xlab(element_blank()) +
    # xlab('Time (min).') +
    facet_wrap(~meas_ID ,
               #scales = "free",
               labeller = as_labeller(plot_names[i])
    ) + # adds column labels to plot
    theme_lewallen() + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()) 
  
  
  
  # tau_a
  # y axis label
  if (i==1){
    y_text = expression(tau[a]~(ms))
  } else {
    y_text = element_blank()
  }
  
  plot_taua[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    # geom_ribbon(aes(x=t,
    #                 ymin=(Cbg-Cbg*Cbe)*1e6,
    #                 ymax=(Cbg+Cbg*Cbe)*1e6),
    #             alpha = ribbon_alpha,
    #             color = "transparent",
    #             fill = ribbon_color
    #             
    # ) +
    geom_line(aes(x=t,y=taua*1e3),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=taua*1e3),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    xlab("Time (min.)") +
    ylab(y_text) +
    labs(shape=element_blank()) + 
    scale_y_continuous(
      limits = c(15,120),
      labels = scales::comma_format(big.mark = ',',decimal.mark = '.')) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    # facet_wrap(~meas_ID,
    # #           scales = "free",
    #            labeller = as_labeller(plot_names[i])) + # adds column labels to plot
    theme_lewallen() + theme(
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank())
  
  
  
  # tau_b
  # y axis label
  if (i==1){
    y_text = expression(tau[b]~(ms))
  } else {
    y_text = element_blank()
  }
  
  plot_taub[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    # geom_ribbon(aes(x=t,
    #                 ymin=(Cbg-Cbg*Cbe)*1e6,
    #                 ymax=(Cbg+Cbg*Cbe)*1e6),
    #             alpha = ribbon_alpha,
    #             color = "transparent",
    #             fill = ribbon_color
    #             
    # ) +
    geom_line(aes(x=t,y=taub*1e3),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=taub*1e3),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    xlab("Time (min.)") +
    ylab(y_text) +
    labs(shape=element_blank()) + 
    scale_y_continuous(
      limits = c(3,15),
      labels = scales::comma_format(big.mark = ',',decimal.mark = '.')) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    # facet_wrap(~meas_ID,
    # #           scales = "free",
    #            labeller = as_labeller(plot_names[i])) + # adds column labels to plot
    theme_lewallen() #+ theme(
  #axis.text.x = element_blank(),
  #axis.ticks.x = element_blank())
  
  
  
  # tau_1
  # y axis label
  if (i==1){
    y_text = expression(tau[1]~(ms))
  } else {
    y_text = element_blank()
  }
  
  plot_tau1[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    # geom_ribbon(aes(x=t,
    #                 ymin=(Cbg-Cbg*Cbe)*1e6,
    #                 ymax=(Cbg+Cbg*Cbe)*1e6),
    #             alpha = ribbon_alpha,
    #             color = "transparent",
    #             fill = ribbon_color
    #             
    # ) +
    geom_line(aes(x=t,y=tau1*1e3),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=tau1*1e3),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    xlab("Time (min.)") +
    ylab(y_text) +
    labs(shape=element_blank()) + 
    scale_y_continuous(
      limits = c(0,12),
      labels = scales::comma_format(big.mark = ',',decimal.mark = '.')) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    # facet_wrap(~meas_ID,
    # #           scales = "free",
    #            labeller = as_labeller(plot_names[i])) + # adds column labels to plot
    theme_lewallen() + theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank())
  
  
  
  # tau_2
  # y axis label
  if (i==1){
    y_text = expression(tau[2]~(ms))
  } else {
    y_text = element_blank()
  }
  
  plot_tau2[[i]] <- ggplot(this_dfs) +
    geom_rect(data = data.frame(xmin = min(comment_times),
                                xmax = max(comment_times),
                                ymin = -Inf,
                                ymax = Inf),
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = box_color, alpha = box_alpha) +
    # geom_ribbon(aes(x=t,
    #                 ymin=(Cbg-Cbg*Cbe)*1e6,
    #                 ymax=(Cbg+Cbg*Cbe)*1e6),
    #             alpha = ribbon_alpha,
    #             color = "transparent",
    #             fill = ribbon_color
    #             
    # ) +
    geom_line(aes(x=t,y=tau2*1e3),
              color=line_color,
              alpha=line_alpha,
              show.legend = FALSE,
              size = line_size) +
    geom_point(aes(x=t,y=tau2*1e3),
               shape = point_shape,
               fill = 'white',
               color = 'black',
               size = point_size,
               stroke = point_size/4,
               alpha = 1) +
    xlab("Time (min.)") +
    ylab(y_text) +
    labs(shape=element_blank()) + 
    scale_y_continuous(
      limits = c(0,5),
      labels = scales::comma_format(big.mark = ',',decimal.mark = '.')) +
    scale_x_continuous(
      limits = time_limits,
      labels = label_number(accuracy = 1)) +
    # facet_wrap(~meas_ID,
    # #           scales = "free",
    #            labeller = as_labeller(plot_names[i])) + # adds column labels to plot
    theme_lewallen() #+ theme(
  #axis.text.x = element_blank(),
  #axis.ticks.x = element_blank())
  
  
  # remove the yaxis labels if it is the second or greater plot
  if (i>1) {
    plot_Ra[[i]] = plot_Ra[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_Rb[[i]] = plot_Rb[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_Rs[[i]] = plot_Rs[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_Ca[[i]] = plot_Ca[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_Cb[[i]] = plot_Cb[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_VDRold[[i]] = plot_VDRold[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_tau_r_abs[[i]] = plot_tau_r_abs[[i]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_tau_r_12[[i]] = plot_tau_r_12[[i]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_Ct[[i]] = plot_Ct[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_TERold[[i]] = plot_TERold[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_TERabs[[i]] = plot_TERabs[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_TEP[[i]] = plot_TEP[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_Va[[i]] = plot_Va[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_Vb[[i]] = plot_Vb[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_GsGt[[i]] = plot_GsGt[[i]] + 
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_taua[[i]] = plot_taua[[i]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_taub[[i]] = plot_taub[[i]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_tau1[[i]] = plot_tau1[[i]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
    plot_tau2[[i]] = plot_tau2[[i]] +
      theme(axis.text.y = element_blank(),
            axis.ticks.y = element_blank())
  } 
  
}

# function used to make columns of figures
plot_a_list <- function(master_list_with_plots, no_of_rows, no_of_cols) {
  patchwork::wrap_plots(master_list_with_plots, 
                        nrow = no_of_rows, ncol = no_of_cols)
}


# MEMBRANE RESISTANCES
ax = list()
for (i in 1:length(plot_Ra)) {
  ax[[i]] <- plot_TERabs[[i]] / plot_Ra[[i]] / plot_Rb[[i]] / plot_Rs[[i]]
  
}
plot_a_list(ax, 1, length(plot_Ra))
file_name_simple <- paste("Membrane Resistance EISandPipette",".svg", sep = "")
ggsave(file_name_simple,
       path = SAVE_path,
       plot = last_plot(),
       width = figure_width_large, height = figure_height_large,
       units = "in",
       limitsize = TRUE) 



# Old data
ax = list()
for (i in 1:length(plot_TERold)) {
  ax[[i]] <- plot_TERold[[i]] / plot_VDRold[[i]]
  
}
plot_a_list(ax, 1, length(plot_VDRold))
file_name_simple <- paste("Traditional EISandPipette",".svg", sep = "")
ggsave(file_name_simple,
       path = SAVE_path,
       plot = last_plot(),
       width = figure_width_large, height = figure_height_medium,
       units = "in",
       limitsize = TRUE) 


# VOLTAGES
ax = list()
for (i in 1:length(plot_TERold)) {
  ax[[i]] <- plot_TEP[[i]] / plot_Va[[i]]
  
}
plot_a_list(ax, 1, length(plot_TEP))
file_name_simple <- paste("Voltages EISandPipette",".svg", sep = "")
ggsave(file_name_simple,
       path = SAVE_path,
       plot = last_plot(),
       width = figure_width_large, height = figure_height_medium,
       units = "in",
       limitsize = TRUE) 


# MEMBRANE CAPACITANCES
ax = list()
for (i in 1:length(plot_Ca)) {
  ax[[i]] <- plot_Ca[[i]] / plot_Cb[[i]]
  
}
plot_a_list(ax, 1, length(plot_Ca))
file_name_simple <- paste("Membrane Capacitance EISandPipette",".svg", sep = "")
ggsave(file_name_simple,
       path = SAVE_path,
       plot = last_plot(),
       width = figure_width_large, height = figure_height_medium,
       units = "in",
       limitsize = TRUE) 


# SHUNT ANALYSIS
ax = list()
for (i in 1:length(plot_Ca)) {
  ax[[i]] <-  plot_GsGt[[i]] / plot_Rs[[i]]
  
}
plot_a_list(ax, 1, length(plot_GsGt))
file_name_simple <- paste("Gs comparison EISandPipette",".svg", sep = "")
ggsave(file_name_simple,
       path = SAVE_path,
       plot = last_plot(),
       width = figure_width_large, height = figure_height_medium,
       units = "in",
       limitsize = TRUE) 



# TAU ratio analysis abs
ax = list()
for (i in 1:length(plot_taub)) {
  ax[[i]] <-  plot_tau_r_abs[[i]] / plot_taua[[i]] / plot_taub[[i]]
  
}
plot_a_list(ax, 1, length(plot_tau_r_abs))
file_name_simple <- paste("tau ratio abs analysis EISandPipette",".svg", sep = "")
ggsave(file_name_simple,
       path = SAVE_path,
       plot = last_plot(),
       width = figure_width_large, height = figure_height_large,
       units = "in",
       limitsize = TRUE) 



# TAU ratio analysis 12
ax = list()
for (i in 1:length(plot_taub)) {
  ax[[i]] <-  plot_tau_r_12[[i]] / plot_tau1[[i]] / plot_tau2[[i]]
  
}
plot_a_list(ax, 1, length(plot_tau_r_12))
file_name_simple <- paste("tau ratio 12 analysis EISandPipette",".svg", sep = "")
ggsave(file_name_simple,
       path = SAVE_path,
       plot = last_plot(),
       width = figure_width_large, height = figure_height_large,
       units = "in",
       limitsize = TRUE) 
