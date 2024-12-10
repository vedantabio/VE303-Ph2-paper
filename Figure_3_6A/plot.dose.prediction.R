## Detection summary plots ----

## ve303.colonize.sum.alt counts detections as 1 and all else as zero in the column N_det_alt USE THIS
#### TO PLOT DETECTIONS WTHOUT FILTERING STATUS DETECTED

# Reorder visit names
ve303.colonize.sum.alt$Visit.Name.Norm = factor(ve303.colonize.sum.alt$Visit.Name.Norm, 
                                                         levels = c("Day -11","Screening", "Day 1", "Day 5", "Day 7", "Day 10", 
                                                                    "Day 12", "Day 13", "Day 14", "Day 15","Day 18", "Day 19", "Day 20", 
                                                                    "Day 21", "Day 26","Day 28", "Day 45", "Day 47","Day 56", 
                                                                    "Day 70", "Day 87", "Day 114", "Day 120", "Day 148", "Day 168"))
# Reorder treatmemt arms
ve303.colonize.sum.alt$TRT.norm = factor(ve303.colonize.sum.alt$TRT.norm,
                                                  levels = c("Placebo", "VE303 Low Dose", "VE303 High Dose"))
# 
# ggplot(ve303.colonize.sum.alt %>%
# #  filter(detection_status =="Detected") %>%
#   filter(.,Visit.Name.Norm %in% c("Screening", "Day 1", "Day 7","Day 14", "Day 28", "Day 56"  )), 
#   aes(x=Visit.Name.Norm, y=N_det_alt, fill=TRT.norm, color=TRT.norm)) + 
#   geom_boxplot(alpha=0.6, show.legend = FALSE, outlier.size = 0) +
#   geom_point(pch = 21, size = 2, alpha=0.3, position = position_jitterdodge(), show.legend = FALSE) +
#   theme_bw() + 
#   facet_wrap(~TRT.norm) + 
#   labs(title= "VE303 Colonization Summary", y="VE303 strains detected (#)", x="Timepoint", fill="Cohort") + 
#   scale_fill_manual(values = coh.cols) + 
#   scale_color_manual(values = coh.cols) + 
#   theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
#         axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
#         strip.text = element_text(size = 14)) + 
#   scale_y_continuous(breaks = c(0:8), limits = c(0,8))
# ggsave(filename = here(results, paste(Sys.Date(), taxplot,"All_VE303_colonization_detected_count_boxplot.pdf", sep = " ")),
#        width = 10, height = 5)
# 
# ####################################################################################
# ####################################################################################
# 
# days_14 <- c("Day 14")
# textsize <- 14
# #  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#  
# ggplot(ve303.colonize.sum.alt %>%
#          filter(.,Visit.Name.Norm %in% c( "Day 14" )), 
#        aes(x=TRT.norm, y=N_det_alt, fill=TRT.norm, color=TRT.norm)) + 
#   geom_boxplot(alpha=.6,  outlier.shape = NA,width = 0.5,show.legend = FALSE ) +
#   geom_point(pch = 21, size = 2.5, alpha = 0.3, position = position_jitterdodge() , show.legend = FALSE ) +
#   theme_bw() +
#   scale_fill_manual(values = coh.cols) + 
#   scale_color_manual(values = coh.cols) + 
#   #scale_y_continuous(breaks = c(0:8), limits = c(0,8))
#   
#   theme(axis.text.x = element_text(size = textsize, angle = 0, hjust = 0.5, vjust = 0.5), 
#         axis.title.x = element_text(size = textsize ),
#         plot.title = element_text(size = textsize ),
#         axis.text.y = element_text(size = textsize, hjust = 1, vjust =0.5) ,     
#         axis.title.y = element_text(size = textsize, hjust = 0.5, vjust = 1) ,
#         strip.text = element_text(size = textsize)) +
#   scale_x_discrete(labels = c("Placebo" = "Placebo" , "VE303 Low Dose" = "VE303\nLow Dose","VE303 High Dose" = "VE303\nHigh Dose" )) +
#   labs(title= "ve303 Total Detection at Day 14", y="VE303 strains detected (#)", x="", 
#        color="TRT.norm", fill="TRT.norm" )
# ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 14 only ve303_total_Proportion_boxplot_det.pdf", sep = " ")), 
#        width = 6, height = 6, dpi = 800)
# 

####################################################################################
####################################################################################
#### Plot ve303 total abundance ----
# Reorder visit names

ve303.abund.subj.total.time.det.alt$Visit.Name = factor(ve303.abund.subj.total.time.det.alt$Visit.Name, 
                                                        levels = c("Screening", "Day 1"  , "Day 7" , "Day 14" , 
                                                                   "Day 28", "Day 56" , "Day 168"))

ve303.abund.subj.total.time.det.alt$TRT.norm = factor(ve303.abund.subj.total.time.det.alt$TRT.norm,
                                                      levels = c("Placebo", "VE303 Low Dose"  , "VE303 High Dose"))


ve303.abund.subj.total.time.det.alt$collect.pre.post.abx= factor(ve303.abund.subj.total.time.det.alt$collect.pre.post.abx,
                                                                 levels = c("Pre.abx", "Post.abx"))

days <- c(  "Day 1"  , "Day 7" , "Day 14" )

#days <- c("Screening", "Day 1"  , "Day 7" , "Day 14" , "Day 28", "Day 56" )



 
# ggplot(ve303.abund.subj.total.time.det.alt  %>%
#          filter(.,Visit.Name%in% days), 
#        aes(x=Visit.Name, y=((est_relative_abundance_panel_alt_sum*100)), fill=TRT.norm, color=TRT.norm)) +
#   geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
#   geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge(0.2)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   facet_grid(~TRT.norm, scales = "free_x", space = "free") +
#   scale_fill_manual(values = coh.cols) +
#   scale_color_manual(values = coh.cols) +
#   theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
#         axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
#         strip.text = element_text(size = 14)) +
#   #scale_y_log10() +
#   labs(title= "ve303 Total Abundance - All Strains", y="Relative Abundance (Total)", x="Timepoint", color="TRT.norm", fill="TRT.norm")
# ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 1-14 ve303_total_abundance_boxplot_det.pdf", sep = " ")), 
#        width = 10, height = 5 )

####################################################################################

# 
# ggplot(ve303.abund.subj.total.time.det.alt  %>%
#          filter(.,Visit.Name%in% days), 
#        aes(x=Visit.Name, y=(log10(0.01 + est_relative_abundance_panel_alt_sum*100)), fill=TRT.norm, color=TRT.norm)) +
#   geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
#   geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge(0.2)) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   facet_grid(~TRT.norm, scales = "free_x", space = "free") +
#   scale_fill_manual(values = coh.cols) +
#   scale_color_manual(values = coh.cols) +
#   theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
#         axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
#         strip.text = element_text(size = 14)) +
#   #scale_y_log10() +
#   labs(title= "ve303 Total Abundance - All Strains", y="Relative Abundance (Total)", x="Timepoint", color="TRT.norm", fill="TRT.norm")
# ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 1-14 ve303_log10 total_abundance_boxplot_det.pdf", sep = " ")), 
#        width = 10, height = 5 )
# 
# 
# 
# 
# ####################################################################################
# ####################################################################################
# 
# days_14 <- c("Day 14")
# textsize <- 14
# ggplot(ve303.abund.subj.total.time.det.alt  %>%
#          filter(.,Visit.Name%in% days_14) %>% filter(.,TRT.norm != "Placebo"), 
#        aes(x=TRT.norm, y=((est_relative_abundance_panel_alt_sum*100) + 0.01), fill=TRT.norm, color=TRT.norm)) +
#   geom_boxplot(alpha=.6,  outlier.shape = NA,width = 0.3,show.legend = FALSE ) +
#   geom_point(pch = 21, size = 2.5, alpha = 0.3, position = position_jitterdodge() , show.legend = FALSE ) +
#   theme_bw() +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
#   #facet_grid(~TRT.norm, scales = "free_x", space = "free") +
#   scale_fill_manual(values = coh.cols) +
#   scale_color_manual(values = coh.cols) +
#   theme(axis.text.x = element_text(size = textsize, angle = 0, hjust = 0.5, vjust = 0.5), 
#         axis.title.x = element_text(size = textsize ),
#         plot.title = element_text(size = textsize ),
#         axis.text.y = element_text(size = textsize, hjust = 1, vjust =0.5) ,     
#         axis.title.y = element_text(size = textsize, hjust = 0.5, vjust = 1) ,
#         strip.text = element_text(size = textsize)) +
#   scale_x_discrete(labels = c("Placebo" = "Placebo" , "VE303 Low Dose" = "VE303\nLow Dose",
#                               "VE303 High Dose" = "VE303\nHigh Dose" )) +
#   scale_y_log10( ) +
#   labs(title= "ve303 Total Abundance at Day 14", y="Relative Abundance (Total)", x="", 
#        color="TRT.norm", fill="TRT.norm" )
# ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 14 only HD LD LOG ve303_total_abundance_boxplot_det.pdf", sep = " ")), 
#        width = 6, height = 6, dpi = 800 )

####################################################################################
####################################################################################

#### Plot ve303 total abundance and proportion as a function of input CFUs ----
# Reorder visit names

ve303.dnameta_CFUs_totalVE303$Visit.Name = factor(ve303.dnameta_CFUs_totalVE303$Visit.Name, 
                                                        levels = c("Screening", "Day 1"  , "Day 7" , "Day 14" , 
                                                                   "Day 28", "Day 56" , "Day 168"))

ve303.dnameta_CFUs_totalVE303$TRT.norm = factor(ve303.dnameta_CFUs_totalVE303$TRT.norm,
                                                      levels = c("Placebo", "VE303 Low Dose"  , "VE303 High Dose"))

ve303.dnameta_CFUs_propVE303$Visit.Name = factor(ve303.dnameta_CFUs_totalVE303$Visit.Name, 
                                                  levels = c("Screening", "Day 1"  , "Day 7" , "Day 14" , 
                                                             "Day 28", "Day 56" , "Day 168"))

ve303.dnameta_CFUs_propVE303$TRT.norm = factor(ve303.dnameta_CFUs_totalVE303$TRT.norm,
                                                levels = c("Placebo", "VE303 Low Dose"  , "VE303 High Dose"))

days <- c(    "Day 7" , "Day 14" )

####################################################################################
####################################################################################

ggplot(ve303.dnameta_CFUs_propVE303  %>%
         filter(.,Visit.Name %in% days), 
       aes(x=CFU_dosed_cumul, y=((  prop_det_alt_sum)) , group = CFU_dosed_cumul, fill=TRT.norm, color=TRT.norm)) +
  geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge(jitter.width = 20, seed=1 )) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(~Visit.Name, scales = "free_x", space = "free") +
  scale_fill_manual(values = coh.cols) +
  scale_color_manual(values = coh.cols) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
#  scale_y_log10() +
  labs(title= "ve303 Total Proportion - All Strains", y="VE303 proportion", x="Cumulative CFUs dosed", color="TRT.norm", fill="TRT.norm")
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 1-14 ve303_total_proportion_Func daily CFUS.pdf", sep = " ")), 
       width = 10, height = 5 )

####################################################################################

ggplot(ve303.dnameta_CFUs_propVE303  %>%
         filter(.,Visit.Name %in%  c(  "Day 14" )), 
       aes(x=CFU_dosed_cumul, y=((  prop_det_alt_sum)) , group = CFU_dosed_cumul, fill=TRT.norm, color=TRT.norm)) +
  #geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge(jitter.width = 20, seed=1 )) +
#  geom_smooth(method = "glm", method.args = list(family = "binomial") ) +
#  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE) +
  #  geom_line(data = fortify(model_eng_cumu), aes(x = CFU_dosed_cumul, y = .fitted)) +
  theme_bw() +
  scale_fill_manual(values = coh.cols) +
  scale_color_manual(values = coh.cols) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  #  scale_y_log10() +
  labs(title= "ve303 Total Proportion - All Strains", y="VE303 proportion", x="Cumulative CFUs dosed", color="TRT.norm", fill="TRT.norm")
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 14 ve303_total_proportion_Func daily CFUS.pdf", sep = " ")), 
       width = 10, height = 5 )

#  
# 
ggplot(ve303.dnameta_CFUs_propVE303  %>%
         filter(.,Visit.Name %in%  c(  "Day 14" )),
       aes(x=CFU_dosed_cumul, y=((  prop_det_alt_sum)) )) +
  geom_point(pch = 21, size = 2, alpha = 0.3   ) +
  #geom_smooth(method = drc::L.4()) +
  #geom_smooth(method = "glm", method.args = list(family = "") ) +
  geom_smooth(method = drm, method.args = list(fct = L.4(fixed=c(NA, 0,1, NA) )), se = FALSE) +
  #  geom_line(data = fortify(model_eng_cumu), aes(x = CFU_dosed_cumul, y = .fitted)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
 
  
  labs(title= "ve303 Total Proportion - All Strains", y="VE303 proportion", x="Cumulative CFUs dosed", color="TRT.norm", fill="TRT.norm")
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 14 ve303_total_proportion_Func daily CFUS.pdf", sep = " ")),
       width = 10, height = 5 )


####################################################################################
####################################################################################

ggplot(ve303.dnameta_CFUs_totalVE303  %>%
         filter(.,Visit.Name%in% days), 
       aes(x=CFU_dosed_cumul, y=((0.01 + est_relative_abundance_panel_alt_sum*100)) , group = CFU_dosed_cumul, fill=TRT.norm, color=TRT.norm)) +
  geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge(1.8)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(~Visit.Name, scales = "free_x", space = "free") +
  scale_fill_manual(values = coh.cols) +
  scale_color_manual(values = coh.cols) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  scale_y_log10() +
  labs(title= "ve303 Total Abundance - All Strains", y="Relative Abundance (Total)", x="Cumulative CFUs dosed", color="TRT.norm", fill="TRT.norm")
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 1-14 ve303_log10 total_abundance_Func daily CFUS.pdf", sep = " ")), 
       width = 10, height = 5 )

####################################################################################

ggplot(ve303.dnameta_CFUs_totalVE303  %>%
         filter(.,Visit.Name%in% days), 
       aes(x=CFU_dosed_cumul, y=((0.01 + est_relative_abundance_panel_alt_sum*100)) , group = CFU_dosed_cumul, fill=TRT.norm, color=TRT.norm)) +
  geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge(1.8)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(~Visit.Name, scales = "free_x", space = "free") +
  scale_fill_manual(values = coh.cols) +
  scale_color_manual(values = coh.cols) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  scale_y_log10() +
  labs(title= "ve303 Total Abundance - All Strains", y="Relative Abundance (Total)", x="Cumulative CFUs dosed", color="TRT.norm", fill="TRT.norm")
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 1-14 ve303_log10 total_abundance_Func daily CFUS.pdf", sep = " ")), 
       width = 10, height = 5 )

####################################################################################



#### Plot Response rate as a function of input CFUs ----
# Reorder visit names

dose_label <- "Cumulative CFU Dose"
pd_label <- "Recurrence Rate"
pd_data_to_plot <- recs_cohort_stats  
gg <- ggplot(data = pd_data_to_plot, aes(x = VE303.Cumulative.Dose,y = mean))
#gg <- gg + xgx_stat_ci(conf_level = .95,distribution = "binomial", geom = list("point","errorbar"))
#gg <- gg + geom_smooth(method = "glm", method.args = list(family = binomial(link = logit)), color = "black")
gg <- gg + geom_errorbar(aes(ymin= lower, ymax= upper), width=.1) +   geom_point() + geom_line()
gg <- gg + scale_y_continuous(breaks = seq(0, 1, 0.25), labels = scales::percent)
gg <- gg + coord_cartesian(ylim = c(0, 1)) 
gg <- gg + labs(x = dose_label, y = pd_label)
pp1 <- gg + ggtitle("Dose vs Response") + 
#  scale_x_log10(breaks = unique(pd_data_to_plot$VE303.Cumulative.Dose))+ ggtitle("Dose on Log Scale") +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), 
        axis.title.x = element_text(size = 14),  axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,
        axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,strip.text = element_text(size = 14)) +

ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Response Rate_Total CFUS.pdf", sep = " ")), 
       width = 7, height = 5 )


####################################################################################


#### Plot Probability of N strains engrafting as a function of input CFUs ----
# Reorder visit names


predprob_5to8 <- HD_2x_predicted[HD_2x_predicted["Proportion detected"] == "5 or more detected",]$`Probability of detection`

AllDose$`Proportion detected` = factor(AllDose$`Proportion detected`, 
                                                 levels = c("8 detected", "7 or more detected"  ,
                                                            "6 or more detected" , "5 or more detected" , 
                                                            "4 or more detected", "3 or more detected" , 
                                                            "2 or more detected", "1 or more detected"))

  
AllDose$Dose = factor(AllDose$Dose,levels = c("HD", "2x HD", "5x HD", "10x HD"))


coh_dose.cols <- c("HD" =  "#5b0076" , "2x HD" = "#08b099", "10x HD" = "#1700e3")

ggplot(AllDose %>% filter(Dose != "5x HD") , 
       aes(x=`Proportion detected`, y=(`Probability of detection`) , group = Dose, fill=Dose, color=Dose)) +
  geom_point(pch = 21, size = 2, alpha = 0.3 ) +
  geom_errorbar(aes(ymin= `Probability of detection`-`Uncertainty of detection` ,
                    ymax= `Probability of detection`+`Uncertainty of detection`), width=.1) + 
  geom_line() +
  #geom_text(x = 2, y = 0.38,  "probability = 34 %", size = 4.5) +
  geom_hline(yintercept = predprob_5to8, linetype="dotted", 
             color = "blue", size=1) +
  
  theme_bw() +
#  facet_grid(~Visit.Name, scales = "free_x", space = "free") +
  scale_fill_manual(values = coh_dose.cols) +
  scale_color_manual(values = coh_dose.cols) +
  scale_x_discrete(labels = c("8 detected" = "N = 8" , "7 or more detected" = "N >= 7",
                              "6 or more detected" = "N >= 6","5 or more detected" = "N >= 5",
                              "4 or more detected" = "N >= 4","3 or more detected" = "N >= 3",
                              "2 or more detected" = "N >= 2","1 or more detected" = "N >= 1")) +
  
  theme(axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, vjust = 0.1), 
        axis.title.x = element_text(size = 14, vjust = 0.8),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.4) ,     
        axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14), legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  labs(title= paste0( "Predicted Probability of Detection based on Logistic Model, ", model_day), y="Probability of Detection", x="# Strains", 
       color="Dose", fill="Dose")
ggsave(filename = here(results, paste(Sys.Date(), taxplot ,model_day,"Predicted Probability of Engraftment.pdf", sep = " ")), 
       width = 8, height = 5 )


####################################################################################
  


ggplot(AllDose %>% filter(Dose %in% c("HD","2x HD" ) ), 
       aes(x=`Proportion detected`, y=(`Probability of detection`) , group = Dose, fill=Dose, color=Dose)) +
  geom_point(pch = 21, size = 2, alpha = 0.3 ) +
  geom_errorbar(aes(ymin= `Probability of detection`-`Uncertainty of detection` ,
                    ymax= `Probability of detection`+`Uncertainty of detection`), width=.1) + 
  geom_line() +
  #geom_text(x = 2, y = 0.38,  "probability = 34 %", size = 4.5) +
  geom_hline(yintercept = predprob_5to8, linetype="dotted", 
             color = "blue", size=1) +
  
  theme_bw() +
  #  facet_grid(~Visit.Name, scales = "free_x", space = "free") +
  scale_fill_manual(values = coh_dose.cols) +
  scale_color_manual(values = coh_dose.cols) +
  scale_x_discrete(labels = c("8 detected" = "N = 8" , "7 or more detected" = "N >= 7",
                              "6 or more detected" = "N >= 6","5 or more detected" = "N >= 5",
                              "4 or more detected" = "N >= 4","3 or more detected" = "N >= 3",
                              "2 or more detected" = "N >= 2","1 or more detected" = "N >= 1")) +
  
  theme(axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, vjust = 0.1), 
        axis.title.x = element_text(size = 14, vjust = 0.8),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.4) ,     
        axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14), legend.text = element_text(size = 14),
        legend.title = element_text(size = 14)) +
  labs(title= paste0( "Predicted Probability of Detection based on Logistic Model, ", model_day), y="Probability of Detection", x="# Strains", 
       color="Dose", fill="Dose")
ggsave(filename = here(results, paste(Sys.Date(), taxplot ,model_day,"Predicted Probability of Engraftment HD_2xHD.pdf", sep = " ")), 
       width = 8, height = 5 )

####################################################################################
