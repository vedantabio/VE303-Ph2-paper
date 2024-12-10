## Detection summary plots ----
# Reorder visit names

rec.cols.prepost <- c("Bad" = "forestgreen", "Good" = "#561600")
rec.cols <- c("Bad" = "indianred2", "Good" = "cyan3")

ve303.colonize.sum.alt$Visit.Name = factor(ve303.colonize.sum.alt$Visit.Name, 
                                           levels = c("Screening", "Day 1"  , "Day 7" , "Day 14" , "Day 28", "Day 56" , "Day 168"))

ve303.colonize.sum.alt$TRT.norm = factor(ve303.colonize.sum.alt$TRT.norm,
                                         levels = c("Placebo", "VE303 Low Dose"  , "VE303 High Dose"))

days <- c("Screening", "Day 1"  , "Day 7" , "Day 14" , "Day 28", "Day 56" )

days_long <- c("Screening", "Day 1"  , "Day 7" , "Day 14" , "Day 28", "Day 56" , "Day 168")

####################################################################################

ggplot(ve303.colonize.sum.alt %>%
         filter(.,Visit.Name%in% days), 
       aes(x=Visit.Name, y=N_det_alt_sum, fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8)) + 
  
  geom_boxplot(alpha=0.6,  outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge(), show.legend = FALSE) +
  theme_bw() + 
  facet_grid(~TRT.norm, scales="free_x", space = "free") + 
  labs(title= "VE303 Colonization Summary", y="VE303 strains detected (#)", x="Timepoint", fill="rec.diagnosis.Wk8") + 
  scale_fill_manual(values = rec.cols) + 
  scale_color_manual(values = rec.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 13)) + 
  scale_y_continuous(limits = c(0,8), breaks = c(0,1,2,3,4,5,6,7,8,9,10,11))
ggsave(filename = here(results, paste(Sys.Date(), taxplot, 
                                      "Filtered_ve303_colonization_detected_count_boxplot_RESPONSE_TRTcol.pdf", sep = " ")),
       width = 10, height = 5)

####################################################################################

ggplot(ve303.colonize.sum.alt %>%
         filter(.,Visit.Name%in% days), 
       aes(x=Visit.Name, y=N_det_alt_sum, fill=collect.pre.post.rec, color=collect.pre.post.rec)) + 
  
  geom_boxplot(alpha=0.6,  outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge(), show.legend = FALSE) +
  theme_bw() + 
  facet_grid(~TRT.norm, scales="free_x", space = "free") + 
  labs(title= "VE303 Colonization Summary", y="VE303 strains detected (#)", x="Timepoint", fill="collect.pre.post.rec") + 
  scale_fill_manual(values = prepostrec.cols) + 
  scale_color_manual(values = prepostrec.cols) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 13)) + 
  scale_y_continuous(limits = c(0,8), breaks = c(0,1,2,3,4,5,6,7,8,9,10,11))
ggsave(filename = here(results, paste(Sys.Date(), taxplot,
                                      "Filtered_ve303_colonization_detected_count_boxplot_Pre_Post_RESPONSE.pdf", sep = " ")), width = 10, height = 5)

####################################################################################
####################################################################################

ggplot(ve303.colonize.sum.alt %>%
         filter(.,Visit.Name%in% days),  
       aes(x=Visit.Name, y=N_det_alt_sum, fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8, group=Subject.Number)) +
  
  geom_line(aes(group = Subject.Number), alpha = 0.3) +
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "median", geom = "line",  size = 0.8)+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = median_q1q3, geom = "linerange")+
  theme_bw() +
  scale_fill_manual(values = rec.cols) + 
  scale_color_manual(values = rec.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 13)) + 
  facet_grid(~TRT.norm, scales = "free") +
  labs(title= "VE303 Colonization Summary", y="VE303 strains detected (#)", x="Timepoint", fill="rec.diagnosis.Wk8") + 
  ggsave(filename = here(results, paste(Sys.Date(),taxplot, "Filtered_ve303_colonization_detected_count_lineplot_RESPONSE.pdf", sep = " ")), width = 12, height = 5)

####################################################################################

ggplot(ve303.colonize.sum.alt %>%
         filter(.,Visit.Name%in% days),  
       aes(x=Visit.Name, y=N_det_alt_sum, fill=collect.pre.post.rec, color=collect.pre.post.rec, group=Subject.Number)) +
  
  #geom_line(aes(group = Subject.Number), alpha = 0.3) +
  
  stat_summary(aes(group=collect.pre.post.rec), fun = "median", geom = "line",  size = 0.8)+
  stat_summary(aes(group=collect.pre.post.rec),fun.data = median_q1q3, geom = "linerange")+
  theme_bw() +
  scale_fill_manual(values = prepostrec.cols) + 
  scale_color_manual(values = prepostrec.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 13)) + 
  facet_grid(~TRT.norm, scales = "free") +
  labs(title= "VE303 Colonization Summary", y="VE303 strains detected (#)", x="Timepoint", fill="collect.pre.post.rec") + 
  ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Filtered_ve303_colonization_detected_count_lineplot_Pre_Post_RESPONSE.pdf", sep = " ")), width = 12, height = 5)


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


####################################################################################
####################################################################################
ggplot(ve303.abund.subj.total.time.det.alt  %>%
         filter(.,Visit.Name%in% days), 
       aes(x=Visit.Name, y=est_relative_abundance_panel_alt_sum*100, fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8)) +
  geom_boxplot(alpha=.6, outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.5, position = position_jitterdodge()) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(~TRT.norm, scales = "free_x", space = "free") +
  scale_fill_manual(values = rec.cols) +
  scale_color_manual(values = rec.cols) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 13)) + 
  
  labs(title= "VE303 Total Abundance - All Strains", y="VE303 Relative Abundance (%)",
       x="Timepoint", color="rec.diagnosis.Wk8", fill="rec.diagnosis.Wk8")
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"ve303_total_abundance_boxplot_det_response_TRTcol.pdf", sep = " ")),
       width = 10, height = 5)
####################################################################################

ggplot(ve303.abund.subj.total.time.det.alt  %>%
         filter(.,Visit.Name%in% days), 
       aes(x=Visit.Name, y=est_relative_abundance_panel_alt_sum*100, fill=collect.pre.post.rec, 
           color=collect.pre.post.rec)) +
  geom_boxplot(alpha=.6, outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.5, position = position_jitterdodge()) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(~TRT.norm, scales = "free_x", space = "free") +
  scale_fill_manual(values = prepostrec.cols) +
  scale_color_manual(values = prepostrec.cols) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 13)) + 
  
  labs(title= "VE303 Total Abundance - All Strains", y="VE303 Relative Abundance (%)",
       x="Timepoint", color="rec.diagnosis.Wk8", fill="rec.diagnosis.Wk8")
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"ve303_total_abundance_boxplot_det_response_TRTcol_Pre_Post_Rec.pdf", sep = " ")),
       width = 10, height = 5)

####################################################################################

ggplot(ve303.abund.subj.total.time.det.alt  %>%
         filter(.,Visit.Name%in% days),  
       aes(x=Visit.Name, y=(0.01+(est_relative_abundance_panel_alt_sum*100)), fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8)) +
  geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 1, position = position_jitterdodge()) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(~TRT.norm, scales = "free_x", space = "free") +
  scale_fill_manual(values = rec.cols) + 
  scale_color_manual(values = rec.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 13)) + 
  scale_y_log10() +
  labs(title= "ve303 Total Abundance - All Strains", y="VE303 Relative Abundance (%)", x="Timepoint", color="rec.diagnosis.Wk8", fill="rec.diagnosis.Wk8")
ggsave(filename = here(results, paste(Sys.Date(),taxplot, "ve303_total_abundance_boxplot_log_det_response.pdf", sep = " ")), 
       width = 10, height = 5)


####################################################################################

ggplot(ve303.abund.subj.total.time.det.alt  %>%
         filter(.,Visit.Name%in% days), 
       aes(x=Visit.Name, y=((est_relative_abundance_panel_alt_sum*100)), fill=TRT.norm, color=TRT.norm)) +
  geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(~TRT.norm, scales = "free_x", space = "free") +
  scale_fill_manual(values = coh.cols) +
  scale_color_manual(values = coh.cols) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  #scale_y_log10() +
  labs(title= "ve303 Total Abundance - All Strains", y="VE303 Relative Abundance (%)", x="Timepoint", color="TRT.norm", fill="TRT.norm")
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"ve303_total_abundance_boxplot_det.pdf", sep = " ")), 
       width = 10, height = 5 )

####################################################################################

abuntime <- ggplot(ve303.abund.subj.total.time.det.alt  %>%
         filter(.,Visit.Name%in% days_long), 
       aes(x=Visit.Name, y=((est_relative_abundance_panel_alt_sum*100)), 
           fill=TRT.norm, color=TRT.norm, show.legend = TRUE)) +
  geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  facet_grid(~TRT.norm, scales = "free_x", space = "free") +
  scale_fill_manual(values = coh.cols) +
  scale_color_manual(values = coh.cols) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  #scale_y_log10() +
  labs(title= "ve303 Total Abundance - All Strains", 
       y="VE303 Relative Abundance (%)", x="Timepoint" , fill = "Treatment", color = "Treatment") #+ guide_legend(title="Treatment")
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 168 ve303_total_abundance_boxplot_det.pdf", sep = " ")), 
       width = 10, height = 5 )

####################################################################################
####################################################################################

days_14 <- c("Day 14")
textsize <- 14
ve303.abund.subj.total.time.det.alt$TRT.norm = factor(ve303.abund.subj.total.time.det.alt$TRT.norm,
                                         levels = c("VE303 High Dose", "VE303 Low Dose"  ,"Placebo" ))

d14_abun <- ggplot(ve303.abund.subj.total.time.det.alt  %>%
         filter(.,Visit.Name%in% days_14), 
       aes(x=TRT.norm, y=((est_relative_abundance_panel_alt_sum*100)), fill=TRT.norm, color=TRT.norm)) +
  geom_boxplot(alpha=.6,  outlier.shape = NA,width = 0.5,show.legend = FALSE ) +
  geom_point(pch = 21, size = 2.5, alpha = 0.3, position = position_jitterdodge() , show.legend = FALSE ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  #facet_grid(~TRT.norm, scales = "free_x", space = "free") +
  scale_fill_manual(values = coh.cols) +
  scale_color_manual(values = coh.cols) +
  theme(axis.text.x = element_text(size = textsize, angle = 0, hjust = 0.5, vjust = 0.5), 
        axis.title.x = element_text(size = textsize ),
        plot.title = element_text(size = textsize ),
        axis.text.y = element_text(size = textsize, hjust = 1, vjust =0.5) ,     
        axis.title.y = element_text(size = textsize, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = textsize)) +
  scale_x_discrete(labels = c("Placebo" = "Placebo" , "VE303 Low Dose" = "VE303\nLow Dose",
                              "VE303 High Dose" = "VE303\nHigh Dose" )) +
  labs(title= "VE303 Total Abundance at Day 14", y="VE303 Relative Abundance (%)", x="", 
       color="TRT.norm", fill="TRT.norm" )
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 14 only ve303_total_abundance_boxplot_det.pdf", sep = " ")), 
       width = 6, height = 6, dpi = 800)




####################################################################################

days_14 <- c("Day 14")
textsize <- 14
ggplot(ve303.abund.subj.total.time.det.alt  %>%
         filter(.,Visit.Name%in% days_14) %>% filter(.,TRT.norm != "Placebo"), 
       aes(x=TRT.norm, y=((est_relative_abundance_panel_alt_sum*100)), fill=TRT.norm, color=TRT.norm)) +
  geom_boxplot(alpha=.6,  outlier.shape = NA,width = 0.3,show.legend = FALSE ) +
  geom_point(pch = 21, size = 2.5, alpha = 0.3, position = position_jitterdodge() , show.legend = FALSE ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  #facet_grid(~TRT.norm, scales = "free_x", space = "free") +
  scale_fill_manual(values = coh.cols) +
  scale_color_manual(values = coh.cols) +
  theme(axis.text.x = element_text(size = textsize, angle = 0, hjust = 0.5, vjust = 0.5), 
        axis.title.x = element_text(size = textsize ),
        plot.title = element_text(size = textsize ),
        axis.text.y = element_text(size = textsize, hjust = 1, vjust =0.5) ,     
        axis.title.y = element_text(size = textsize, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = textsize)) +
  scale_x_discrete(labels = c("Placebo" = "Placebo" , "VE303 Low Dose" = "VE303\nLow Dose",
                              "VE303 High Dose" = "VE303\nHigh Dose" )) +
  #scale_y_log10() +
  labs(title= "ve303 Total Abundance at Day 14", y="VE303 Relative Abundance (%)", x="", 
       color="TRT.norm", fill="TRT.norm" )
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 14 only HD LD ve303_total_abundance_boxplot_det.pdf", sep = " ")), 
       width = 6, height = 6, dpi = 800)

####################################################################################

days_14 <- c("Day 14")
textsize <- 14
ggplot(ve303.abund.subj.total.time.det.alt  %>%
         filter(.,Visit.Name%in% days_14) %>% filter(.,TRT.norm != "Placebo"), 
       aes(x=TRT.norm, y=((est_relative_abundance_panel_alt_sum*100) + 0.01), fill=TRT.norm, color=TRT.norm)) +
  geom_boxplot(alpha=.6,  outlier.shape = NA,width = 0.3,show.legend = FALSE ) +
  geom_point(pch = 21, size = 2.5, alpha = 0.3, position = position_jitterdodge() , show.legend = FALSE ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  #facet_grid(~TRT.norm, scales = "free_x", space = "free") +
  scale_fill_manual(values = coh.cols) +
  scale_color_manual(values = coh.cols) +
  theme(axis.text.x = element_text(size = textsize, angle = 0, hjust = 0.5, vjust = 0.5), 
        axis.title.x = element_text(size = textsize ),
        plot.title = element_text(size = textsize ),
        axis.text.y = element_text(size = textsize, hjust = 1, vjust =0.5) ,     
        axis.title.y = element_text(size = textsize, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = textsize)) +
  scale_x_discrete(labels = c("Placebo" = "Placebo" , "VE303 Low Dose" = "VE303\nLow Dose",
                              "VE303 High Dose" = "VE303\nHigh Dose" )) +
  scale_y_log10( ) +
  labs(title= "ve303 Total Abundance at Day 14", y="VE303 Relative Abundance (%)", x="", 
       color="TRT.norm", fill="TRT.norm" )
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 14 only HD LD LOG ve303_total_abundance_boxplot_det.pdf", sep = " ")), 
       width = 6, height = 6, dpi = 800 )

####################################################################################


#### Plot ve303 total abundance lineplots ----
# Reorder visit names

ve303.abund.total303.subj.det.alt$Visit.Name = factor(ve303.abund.total303.subj.det.alt$Visit.Name,
                                                      levels = c("Screening", "Day 1"  , "Day 7" ,
                                                                 "Day 14" , "Day 28", "Day 56" , "Day 168"))


ve303.abund.total303.subj.det.alt$TRT.norm = factor(ve303.abund.total303.subj.det.alt$TRT.norm,
                                                    levels = c("Placebo", "VE303 Low Dose"  , "VE303 High Dose"))
####################################################################################

ggplot(ve303.abund.total303.subj.det.alt %>% filter(Visit.Name %in% days),
       aes(x=Visit.Name, y=(100*subject_mean),fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8, group=Subject.Number)) +
  geom_line(aes(group = Subject.Number), alpha = 0.3) +
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "median", geom = "line",  size = 0.8)+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = median_q1q3, geom = "linerange")+
  theme_bw() +
  facet_grid(~TRT.norm, scales = "free") +
  scale_fill_manual(values = rec.cols) + 
  scale_color_manual(values = rec.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  labs(x="Timepoint", y="Total ve303 Relative Abundance (%)")
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"ve303_total_abundance_lineplot_det_RESPONSE.pdf", sep = " ")), width = 12, height = 5)

####################################################################################

ggplot(ve303.abund.total303.subj.det.alt %>% filter(Visit.Name %in% days),
       aes(x=Visit.Name, y=(0.01+(subject_mean*100)),fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8,group=Subject.Number)) +
  geom_line(aes(group = Subject.Number), alpha = 0.3) +
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "median", geom = "line",  size = 0.8)+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = median_q1q3, geom = "linerange")+
  scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = rec.cols) + 
  scale_color_manual(values = rec.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  facet_grid(~TRT.norm, scales = "free") +
  labs(title= "ve303 Total Abundance - All Strains",x="Timepoint", y="Log(Total ve303 Relative Abundance (%) + 0.01)")
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"ve303_total_abundance_lineplot_det_logscale_RESPONSE.pdf", sep = " ")), width = 12, height = 5)

####################################################################################

#### Plot ve303 relative abundance by Strain ----

#theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
#      axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
#      strip.text = element_text(size = 14)) +
  
# Reorder visit names
ve303.abund.sum.subj.det.alt$Visit.Name = factor(ve303.abund.sum.subj.det.alt$Visit.Name,
                                                 c("Screening", "Day 1"  , "Day 7" ,
                                                   "Day 14" , "Day 28", "Day 56" , "Day 168"))

ggplot(ve303.abund.sum.subj.det.alt %>% filter(Visit.Name %in% c("Screening", "Day 1"  , "Day 7" ,
                                                                 "Day 14" , "Day 28", "Day 56" , "Day 168")),
       aes(x=Visit.Name, y=(subject_mean*100),fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8, group=Subject.Number)) +
  geom_line(aes(group = Subject.Number), alpha = 0.3) +
  geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "mean", geom = "line",  size = 0.8)+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = mean_se, geom = "linerange")+
  
  theme_bw() +
  facet_grid(organism~TRT.norm, scales = "fixed") +
  scale_fill_manual(values = rec.cols) + 
  scale_color_manual(values = rec.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  labs(x="Timepoint", y="Relative Abundance (%)") +
  ggsave(filename = here(results, paste(Sys.Date(), taxplot,"ve303_colonization_meanabund_det_perstrain_RESPONSE.pdf", sep = " ")),
         width = 12, height = 12)

####################################################################################

## Plot ve303 relative abundance by Strain on a log scale
# Reorder visit names

ve303.abund.sum.subj.det.alt$subject_mean_pseudocount <- (0.01 + 100*ve303.abund.sum.subj.det.alt$subject_mean)

# Reorder visit names
ve303.abund.sum.subj.det.alt_Wk8R$Visit.Name = factor(ve303.abund.sum.subj.det.alt_Wk8R$Visit.Name,
                                                      c("Screening", "Day 1"  , "Day 7" ,
                                                        "Day 14" , "Day 28", "Day 56" , "Day 168"))
ve303.abund.sum.subj.det.alt_Wk8R$subject_mean_pseudocount <- (0.01 + 100*ve303.abund.sum.subj.det.alt_Wk8R$subject_mean)

####################################################################################

ggplot(ve303.abund.sum.subj.det.alt %>% filter(Visit.Name  %in% days),
       aes(x=Visit.Name , y= (subject_mean_pseudocount),fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8, group=Subject.Number)) +
  geom_line(aes(group = Subject.Number), alpha = 0.15) +
  geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "mean", geom = "line",  size = 0.8)+
  #  stat_summary(aes(group=Wk.8.Response),fun.data = median_mad, geom = "linerange")+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = mean_se , geom = "linerange")+
  scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = rec.cols) +
  scale_color_manual(values = rec.cols) +
  
  facet_grid(organism~TRT.norm, scales = "free_x") +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  
  labs(x="Timepoint", y="Log(Relative Abundance (%) + 0.01)")
ggsave(filename = here(results, paste(Sys.Date(),taxplot, "ve303_colonization_meanabund_det_perstrain_logscale_RESPONSE.pdf", sep = " ")), 
       width = 12, height = 12)

###################################################################################

ggplot(ve303.abund.sum.subj.det.alt %>% filter(Visit.Name  %in% days),
       aes(x=Visit.Name , y= (subject_mean_pseudocount),fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8, group=Subject.Number)) +
  geom_line(aes(group = Subject.Number), alpha = 0.15) +
  geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "median", geom = "line",  size = 0.8)+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = median_q1q3 , geom = "linerange")+
  scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = rec.cols) +
  scale_color_manual(values = rec.cols) +
  
  facet_grid(organism~TRT.norm, scales = "free_x") +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  
  labs(x="Timepoint", y="Log(Relative Abundance (%) + 0.01)")
ggsave(filename = here(results, paste(Sys.Date(),taxplot, "ve303_colonization_medianabund_det_perstrain_logscale_RESPONSE.pdf", sep = " ")), 
       width = 12, height = 12)

####################################################################################

ve303.abund.sum.subj.det.alt_Wk8R$Visit.Name = factor(ve303.abund.sum.subj.det.alt_Wk8R$Visit.Name,
                                                      c("Screening", "Day 1"  , "Day 7" ,
                                                        "Day 14" , "Day 28", "Day 56" , "Day 168"))
ve303.abund.sum.subj.det.alt_Wk8R$subject_mean_pseudocount <- (0.01 + 100*ve303.abund.sum.subj.det.alt_Wk8R$subject_mean)

ggplot(ve303.abund.sum.subj.det.alt_Wk8R %>% filter(Visit.Name  %in% days),
       aes(x=Visit.Name , y= (subject_mean_pseudocount),fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8, group=Subject.Number)) +
  geom_line(aes(group = Subject.Number), alpha = 0.15) +
  geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "median", geom = "line",  size = 0.4)+
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "mean", geom = "line", linetype = 'dotted',  size = 0.8)+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = mean_se, geom = "linerange")+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = median_q1q3, geom = "linerange")+
  scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = rec.cols) +
  scale_color_manual(values = rec.cols) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) , axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  
  facet_grid(organism~TRT.norm, scales = "free_x") +
  labs(x="Timepoint", y="Log(Relative Abundance (%) + 0.01)")
ggsave(filename = here(results, paste(Sys.Date(),taxplot, "ve303_colonization_meanmedianabund_det_perstrain_logscale_RESPONSE.pdf", sep = " ")), width = 11, height = 11)

####################################################################################

ve303.abund.sum.subj.det.alt_Wk8R$Visit.Name = factor(ve303.abund.sum.subj.det.alt_Wk8R$Visit.Name,
                                                      c("Screening", "Day 1"  , "Day 7" ,
                                                        "Day 14" , "Day 28", "Day 56" , "Day 168"))
ve303.abund.sum.subj.det.alt_Wk8R$subject_mean_pseudocount <- (0.01 + 100*ve303.abund.sum.subj.det.alt_Wk8R$subject_mean)

ggplot(ve303.abund.sum.subj.det.alt_Wk8R %>% filter(Visit.Name  %in% days),
       aes(x=Visit.Name , y= (subject_mean_pseudocount),fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8, group=Subject.Number)) +
  geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "median", geom = "line",  size = 0.4)+
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "mean", geom = "line", linetype = 'dotted',  size = 0.8)+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = mean_se, geom = "linerange")+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = median_q1q3, geom = "linerange")+
  scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = rec.cols) +
  scale_color_manual(values = rec.cols) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) , axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  
  facet_grid(organism~TRT.norm, scales = "free_x") +
  labs(x="Timepoint", y="Log(Relative Abundance (%) + 0.01)")
ggsave(filename = here(results, paste(Sys.Date(),taxplot, "ve303_colonization_meanmedianabund_det_perstrain_NOSUBJ_logscale_RESPONSE.pdf", sep = " ")), width = 11, height = 11)

####################################################################################

ve303.abund.sum.subj.det.alt_Wk8R$Visit.Name = factor(ve303.abund.sum.subj.det.alt_Wk8R$Visit.Name,
                                                      c("Screening", "Day 1"  , "Day 7" ,
                                                        "Day 14" , "Day 28", "Day 56" , "Day 168"))
ve303.abund.sum.subj.det.alt_Wk8R$subject_mean_pseudocount <- (0.01 + 100*ve303.abund.sum.subj.det.alt_Wk8R$subject_mean)

ggplot(ve303.abund.sum.subj.det.alt_Wk8R %>% filter(Visit.Name  %in% days & 
                                                      organism %in% c("VE303-02", "VE303-03", "VE303-05", "VE303-08")),
  aes(x=Visit.Name , y= (subject_mean_pseudocount),fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8, group=Subject.Number)) +
  #geom_line(aes(group = Subject.Number), alpha = 0.15) +
  geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "median", geom = "line",  size = 0.4)+
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "mean", geom = "line", linetype = 'dotted',  size = 0.8)+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = mean_se, geom = "linerange")+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = median_q1q3, geom = "linerange")+
  scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = rec.cols) +
  scale_color_manual(values = rec.cols) +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) , axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  
  facet_grid(organism~TRT.norm, scales = "free_x") +
  labs(x="Timepoint", y="Log(Relative Abundance (%) + 0.01)")
ggsave(filename = here(results, paste(Sys.Date(),taxplot, 
                                      "ve303_colonization_meanmedianabund_det_SIG 2358 perstrain_NOSUBJ_logscale_RESPONSE.pdf", sep = " ")), 
       width = 9, height = 7)

####################################################################################
####################################################################################

ggplot(ve303.abund.sum.subj.det.alt_Wk8R %>% filter(Visit.Name  %in% days & 
                                                      organism %in% c("VE303-02",   "VE303-08")),
       aes(x=Visit.Name , y= (subject_mean_pseudocount),fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8, group=Subject.Number)) +
  #geom_line(aes(group = Subject.Number), alpha = 0.15) +
  geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "median", geom = "line",  size = 0.4)+
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "mean", geom = "line", linetype = 'dotted',  size = 0.8)+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = mean_se, geom = "linerange")+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = median_q1q3, geom = "linerange")+
  
  scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = rec.cols) +
  scale_color_manual(values = rec.cols) +
  
  facet_grid(organism~TRT.norm, scales = "free_x") +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) , axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  labs(x="Timepoint", y="Log(Relative Abundance (%) + 0.01)")
ggsave(filename = here(results, paste(Sys.Date(),taxplot, 
                                      "ve303_colonization_meanmedianabund_det_SIG 28 perstrain_NOSUBJ_logscale_RESPONSE.pdf", sep = " ")), 
       width = 9, height = 7)


####################################################################################

ggplot(ve303.abund.sum.subj.det.alt_Wk8R %>% filter(Visit.Name  %in% days & TRT.norm != "Placebo" ),
       aes(x=Visit.Name , y= (subject_mean_pseudocount),fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8, group=Subject.Number)) +
  #geom_line(aes(group = Subject.Number), alpha = 0.15) +
  geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "median", geom = "line",  size = 0.4)+
  stat_summary(aes(group=rec.diagnosis.Wk8), fun = "mean", geom = "line", linetype = 'dotted',  size = 0.8)+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = mean_se, geom = "linerange")+
  stat_summary(aes(group=rec.diagnosis.Wk8),fun.data = median_q1q3, geom = "linerange")+
  
  scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = rec.cols) +
  scale_color_manual(values = rec.cols) +
  
  facet_grid(organism~., scales = "free_x") +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) , axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  labs(x="Timepoint", y="Log(Relative Abundance (%) + 0.01)")
ggsave(filename = here(results, paste(Sys.Date(),taxplot, "ve303_colonization_meanmedianabund_det_DOSED perstrain_NOSUBJ_logscale_RESPONSE.pdf", sep = " ")), 
       width = 9, height = 7)

####################################################################################
 
# 
# ## Plot ve303 relative abundance by Strain on a log scale
# # Reorder visit names
ve303.abund.sum.subj.det.alt$TRT.norm <- factor(ve303.abund.sum.subj.det.alt$TRT.norm,
                                                c("Placebo", "VE303 Low Dose" , "VE303 High Dose"))

ve303.abund.sum.subj.det.alt$subject_mean_pseudocount <- (0.01 + 100*ve303.abund.sum.subj.det.alt$subject_mean)

ggplot(ve303.abund.sum.subj.det.alt %>% filter(Visit.Name  %in% days),
       aes(x=Visit.Name , y= (subject_mean_pseudocount),fill=TRT.norm, color=TRT.norm, group=Subject.Number)) +
  geom_line(aes(group = Subject.Number), alpha = 0.15) +
  geom_hline(yintercept = 0.1, linetype="longdash", color="red") +

  stat_summary(aes(group=TRT.norm), fun = "median", geom = "line",  size = 0.8)+
  #  stat_summary(aes(group=TRT.norm),fun.data = mean_se, geom = "linerange")+
  stat_summary(aes(group=TRT.norm),fun.data = median_q1q3 ,geom = "linerange")+
  scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = coh.cols) +
  scale_color_manual(values = coh.cols) +
  facet_grid(organism~TRT.norm, scales = "free_x") +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) , axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  labs(x="Timepoint", y="Log(Relative Abundance (%) + 0.01)")
ggsave(filename = here(results, paste(Sys.Date(),taxplot, "ve303_colonization_abund_det_perstrain_median_logscale.pdf", sep = " ")), width = 12, height = 12)

####################################################################################
 
ve303.abund.sum.subj.det.alt$subject_mean_pseudocount <- (0.001 + 100*ve303.abund.sum.subj.det.alt$subject_mean)

ggplot(ve303.abund.sum.subj.det.alt %>% filter(Visit.Name  %in% days),
       aes(x=Visit.Name , y= (subject_mean_pseudocount),fill=TRT.norm, color=TRT.norm, group=Subject.Number)) +
  #geom_line(aes(group = Subject.Number), alpha = 0.15) +
  geom_hline(yintercept = 0.01, linetype="longdash", color="red") +
  stat_summary(aes(group=TRT.norm), fun = "mean", geom = "line", linetype = "dotted",  size = 0.8)+
  stat_summary(aes(group=TRT.norm),fun.data = mean_se, geom = "linerange")+
  stat_summary(aes(group=TRT.norm), fun = "median", geom = "line",  size = 0.4)+
  #stat_summary(aes(group=TRT.norm),fun.data = median_q1q3, geom = "linerange")+
  
  scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = coh.cols) +
  scale_color_manual(values = coh.cols) +
  facet_grid(organism~TRT.norm, scales = "free_x") +
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) , axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) +
  labs(x="Timepoint", y="Log(Relative Abundance (%) + 0.001)", fill="Treatment", color="Treatment")
ggsave(filename = here(results, paste(Sys.Date(),taxplot, "ve303_colonization_meanmedianabund_det_perstrain_NO SUBJ_logscale.pdf", sep = " ")), width = 12, height = 12)

   
#############################################################################################################


ve303.dnameta.det.alt$Visit.Name = factor(ve303.dnameta.det.alt$Visit.Name,
                                            levels = c("Screening", "Day 1"  , "Day 7" ,
                                                       "Day 14" , "Day 28", "Day 56" , "Day 168"))


ve303.dnameta.det.alt$TRT.norm = factor(ve303.dnameta.det.alt$TRT.norm,
                                       levels = c("Placebo", "VE303 Low Dose"  , "VE303 High Dose"))

#############################################################################################################
#############################################################################################################

days1 <- c("Day 7"  ,   "Day 14" )
#days1 <- c("Day 7"  ,   "Day 14"  ,  "Day 28")
 
ggplot(ve303.dnameta.det.alt %>% filter(Visit.Name  %in% days1 &
                                          organism %in% c("VE303-02","VE303-03" , "VE303-05", "VE303-08") &
                                          TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose") &
                                          collect.pre.post.rec %in% c("No.rec", "Pre.rec")),
       aes(x=collect.pre.post.rec , y= (0.01 + 100*est_relative_abundance_panel_alt),fill=collect.pre.post.rec,
           color=collect.pre.post.rec )) +
  geom_boxplot(alpha=.6, outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +

  scale_y_log10()+
  theme_bw() +
  facet_wrap(organism~., scales = "free_x") +
  scale_fill_manual(values = prepostrec.cols) +
  scale_color_manual(values = prepostrec.cols) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)) +
  theme(axis.text.y = element_text(size=12), strip.text = element_text(size=12),
        axis.text.x = element_text(size = 10), strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x="", y="Log(Relative Abundance (%) + 0.01)")
ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_abund_det_SIGperstrain_NO_TIME_logscale.pdf",
                                      sep = " ")), width = 4, height = 5)

##################################################################


ggplot(ve303.dnameta.det.alt %>% filter(Visit.Name  %in% days1 &
                                          TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose") &
                                          collect.pre.post.rec %in% c("No.rec", "Pre.rec")),
       aes(x=collect.pre.post.rec , y= (0.01 + 100*est_relative_abundance_panel_alt),fill=collect.pre.post.rec,
           color=collect.pre.post.rec )) +
  geom_boxplot(alpha=.6, outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  
  scale_y_log10()+
  theme_bw() +
  facet_wrap(organism~., scales = "free_x") +
  scale_fill_manual(values = prepostrec.cols) +
  scale_color_manual(values = prepostrec.cols) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)) +
  theme(axis.text.y = element_text(size=12), strip.text = element_text(size=12),
        axis.text.x = element_text(size = 10), strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x="", y="Log(Relative Abundance (%) + 0.01)")
ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_abund_det_ALLperstrain_NO_TIME_logscale.pdf",
                                      sep = " ")), width = 4, height = 5)

##################################################################
#################################################################

ggplot(ve303.dnameta.det.alt %>% filter(Visit.Name  %in% days1 &
                                          TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose") &
                                          collect.pre.post.rec %in% c("No.rec", "Pre.rec")),
       aes(x=collect.pre.post.rec , y= (0.01 + 100*est_relative_abundance_panel_alt_sum),fill=collect.pre.post.rec,
           color=collect.pre.post.rec )) +
  geom_boxplot(alpha=.6, outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +

  scale_y_log10()+
  theme_bw() +
  #  facet_grid(organism~., scales = "free_x") +
  scale_fill_manual(values = prepostrec.cols) +
  scale_color_manual(values = prepostrec.cols) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12)) +
  theme(axis.text.y = element_text(size=12), strip.text = element_text(size=12),
        axis.text.x = element_text(size = 10), strip.text.x = element_text(size=12), strip.text.y = element_text(size=12)) +
  labs(x="", y="Log(Relative Abundance (%) + 0.01)")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_abund_det_SIGTotal_NO_TIME_logscale.pdf",
                                      sep = " ")), width = 6, height = 7)
 
##################################################################
 
 
ggplot(ve303.dnameta.det.alt %>% filter(Visit.Name  %in% days1 &
                                          TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose") &
                                          collect.pre.post.rec.Wk8 %in% c("No.rec", "Pre.rec")),
       aes(x=rec.TRT , y= (0.01 + 100*est_relative_abundance_panel_alt_sum),fill=collect.pre.post.rec.Wk8,
           color=collect.pre.post.rec.Wk8 )) +
  geom_boxplot(alpha=.6, outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +

  scale_y_log10()+
  theme_bw() +
  #  facet_grid(organism~., scales = "free_x") +
  scale_fill_manual(values = prepostrec.cols) +
  scale_color_manual(values = prepostrec.cols) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 20)) +
  theme(axis.text.y = element_text(size=20), strip.text = element_text(size=20),
        axis.title.y = element_text(size = 20)) +
  labs(x="", y="Log(Relative Abundance (%) + 0.01)")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_interactionRecDose_SIGTotal_NO_TIME_logscale.pdf",
                                      sep = " ")), width = 9, height = 8)


##################################################################


ggplot(ve303.dnameta.det.alt %>% filter(Visit.Name  %in% days1 &
                                          TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose") &
                                          collect.pre.post.rec.Wk8 %in% c("No.rec", "Pre.rec") ),
       aes(x=rec.TRT , y= (0.01 + 100*est_relative_abundance_panel_alt), fill = collect.pre.post.rec.Wk8,
           color=collect.pre.post.rec.Wk8 )) +
  geom_boxplot(alpha=.6, outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +

  scale_y_log10()+
  theme_bw() +
  facet_wrap(organism~., scales = "free_x") +
  scale_fill_manual(values = prepostrec.cols) +
  scale_color_manual(values = prepostrec.cols) +

  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 12),
        #axis.text.x = element_blank(),
        strip.text.x = element_text(size=14), strip.text.y = element_text(size=14)) +
  scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
                            "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD")) +
  labs(x="", y="Log(Relative Abundance (%) + 0.01)")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_interactionRecDose_ALL_NO_TIME_logscale.pdf",
                                      sep = " ")), width = 12, height = 11)
 
##################################################################

days2 <- c("Day 7"  ,   "Day 14" )
days1 <- c("Day 7"  ,   "Day 14" )
#days1 <- c("Day 7"  ,   "Day 14"  ,  "Day 28")

ve303.dnameta.det.alt_early <- ve303.dnameta.det.alt %>% filter(Visit.Name  %in% days1 )
 
vvtest <- ve303.dnameta.det.alt_early %>%
  select(Subject.Number, organism, TRT.norm, Visit.Name.Norm, est_relative_abundance_panel_alt,
         est_relative_abundance_panel_alt_sum, MGS.Sample.ID, rec.TRT) %>%
  group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>%
  spread(organism, est_relative_abundance_panel_alt, fill=0) %>% # This fills in zeros for all events not detected - needed to calc mean
  gather(organism, est_relative_abundance_panel_alt, 7: ncol(.)) %>%
  group_by(Subject.Number, organism) %>%
  mutate(subject_mean_time = mean(est_relative_abundance_panel_alt)) %>%
  ungroup() %>%
  group_by(Subject.Number) %>%
  mutate(total_mean_time = mean(est_relative_abundance_panel_alt_sum)) %>%
  ungroup() %>%
  distinct(left_join(., dnameta, by = "MGS.Sample.ID"))

vvtest <- vvtest %>%
  rename_at(vars(ends_with(".x")),  ~str_replace(., "\\..$","") ) %>%
  select_at(vars(-ends_with(".y"))) %>%
  distinct(Subject.Number, organism, .keep_all = TRUE)

 
##################################################################

# THIS IS NOT QUITE RIGHT. SINCE SUBJECT ABUNDANCE HAS BEEN AVERAGED AND THERE'S ALREADY JUST ONE PER SUBJECT
# REPRESENTED IN THE DATA, FILTERING by COLLECT PRE POST REC REMOVES A COUPLE VALID SUBJECT POINTS. NEED TO NOT FILTER
# AND PLOT BAD GOOD WITH THE NO REC PRE REC COLORS
# 
# ggplot(vvtest %>% filter(TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose") &
#                            collect.pre.post.rec.Wk8 %in% c("No.rec", "Pre.rec") ),
#        aes(x=rec.TRT , y= (0.01 + 100*subject_mean_time), fill = collect.pre.post.rec.Wk8,
#            color=collect.pre.post.rec.Wk8 )) +
#   geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
#   geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
#   #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
# 
#   scale_y_log10()+
#   theme_bw() +
#   facet_wrap(organism~., scales = "free_x") +
#   scale_fill_manual(values = prepostrec.cols) +
#   scale_color_manual(values = prepostrec.cols) +
# 
#   theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
#         axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 12),
#         #axis.text.x = element_blank(),
#         strip.text.x = element_text(size=14), strip.text.y = element_text(size=14)) +
#   scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
#                             "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD")) +
#   labs(x="", y="Log(Mean Subject Abundance (%) + 0.01)")
# 
# ggsave(filename = here(results, paste(Sys.Date(),taxplot,
#                                       "TESTTESTve303_colonization_interactionRecDose_ALL_NO_TIME_MeanTime_logscale.pdf",
#                                       sep = " ")), width = 12, height = 11)
  

##################################################################
##################################################################

ggplot(vvtest %>% filter(TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose") ),
       aes(x=rec.TRT , y= (0.01 + 100*subject_mean_time), fill =  rec.diagnosis.Wk8,
           color= rec.diagnosis.Wk8 )) +
  geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  
  scale_y_log10()+
  theme_bw() +
  facet_wrap(organism~., scales = "free_x") +
  scale_fill_manual(values =  rec.cols.prepost) +
  scale_color_manual(values =  rec.cols.prepost) +
  
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 12),
        #axis.text.x = element_blank(),
        strip.text.x = element_text(size=14), strip.text.y = element_text(size=14)) +
  scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
                            "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD")) +
  labs(x="", y="Log(Mean Subject Abundance (%) + 0.01)")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "TESTTESTve303_colonization_interactionRecDose_recWk8_ALL_NO_TIME_MeanTime_logscale.pdf",
                                      sep = " ")), width = 12, height = 11)
##################################################################

# 
# ggplot(vvtest %>% filter(organism %in% c(  "VE303-02", "VE303-08") &
#                            TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose") &
#                            collect.pre.post.rec.Wk8 %in% c("No.rec", "Pre.rec")),
#        aes(x=rec.TRT , y= (0.01 + 100*subject_mean_time), fill = collect.pre.post.rec.Wk8,
#            color=collect.pre.post.rec.Wk8 )) +
#   geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
#   geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
#   #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
# 
#   theme_bw() +
#   facet_grid(.~organism , scales = "free_x") +
#   scale_fill_manual(values =  rec.cols) +
#   scale_color_manual(values = rec.cols) +
#   theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
#         axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 12),
#         #axis.text.x = element_blank(),
#         strip.text.x = element_text(size=14), strip.text.y = element_text(size=14)) +
#   scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
#                             "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD")) +
#   scale_y_log10()+
# 
#   labs(x="", y="Log(Mean Subject Abundance (%) + 0.01)")
# 
# ggsave(filename = here(results, paste(Sys.Date(),taxplot,
#                                       "ve303_colonization_interactionRecDose_SIG0208_NO_TIME_MeanTime_logscale.pdf",
#                                       sep = " ")), width = 7, height = 3)


##################################################################

ggplot(vvtest %>% filter(organism %in% c(  "VE303-02", "VE303-08") &
                           TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose")),
       aes(x=rec.TRT , y= (0.01 + 100*subject_mean_time), fill =  rec.diagnosis.Wk8,
           color= rec.diagnosis.Wk8 )) +
  geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  
  theme_bw() +
  facet_grid(.~organism , scales = "free_x") +
  scale_fill_manual(values =  rec.cols.prepost) +
  scale_color_manual(values =  rec.cols.prepost) +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 12),
        #axis.text.x = element_blank(),
        strip.text.x = element_text(size=14), strip.text.y = element_text(size=14)) +
  scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
                            "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD")) +
  scale_y_log10()+
  
  labs(x="", y="Log(Mean Subject Abundance (%) + 0.01)")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_interactionRecDose_recWk8_SIG0208_NO_TIME_MeanTime_logscale.pdf",
                                      sep = " ")), width = 7, height = 3)


# ################################################################### 


vvtest_tot <- vvtest %>% distinct(Subject.Number, .keep_all = TRUE)

ggplot(vvtest_tot %>% filter( TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose") ),
       aes(x=rec.TRT , y= (0.01 + 100*total_mean_time), fill = rec.diagnosis.Wk8,
           color= rec.diagnosis.Wk8 )) +
  geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  
  theme_bw() +
  scale_fill_manual(values =  rec.cols.prepost) +
  scale_color_manual(values =  rec.cols.prepost) +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 12),
        #axis.text.x = element_blank(),
        strip.text.x = element_text(size=14), strip.text.y = element_text(size=14)) +
  scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
                            "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD")) +
  scale_y_log10()+
  
  labs(x="", y="Log(Mean Subject Abundance (%) + 0.01)")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_interactionRecDose_recWk8_Total VE303_NO_TIME_MeanTime_logscale.pdf",
                                      sep = " ")), width = 7, height = 5)
   
# ##################################################################

ggplot(vvtest %>% filter( organism %in% c(  "VE303-01" , "VE303-02" ,"VE303-03","VE303-04",
                                            "VE303-05","VE303-06","VE303-07","VE303-08") &
                            TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose")  ),
       aes(x=rec.diagnosis.Wk8 , y= (0.01 + 100*subject_mean_time), fill =  rec.diagnosis.Wk8,
           color= rec.diagnosis.Wk8 )) +
  geom_boxplot(alpha=.6, outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +

  scale_y_log10()+
  theme_bw() +
  facet_wrap(organism~., scales = "free_x") +
  #  facet_wrap(organism~., scales = "free_x") +
  scale_fill_manual(values =  rec.cols.prepost) +
  scale_color_manual(values =  rec.cols.prepost) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 12)) +
  theme(axis.text.y = element_text(size=12), strip.text = element_text(size=12),
        axis.text.x = element_text(size = 10), strip.text.x = element_text(size=12),
        strip.text.y = element_text(size=12)) +
  #scale_x_discrete(labels=c("No.rec" = "No recurrence", "Pre.rec" = "Pre recurrence")) +

  labs(x="", y="Log(Mean Subject Abundance (%) + 0.01)")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_CombinedDose_interactionRecDose_recWk8_AllStrains_NO_TIME_MeanTime_logscale.pdf",
                                      sep = " ")), width = 9, height = 6)
  
#
 ##################################################################
  ##################################################################

vvtest3 <- ve303.dnameta.det.alt %>%
  distinct(Subject.Number, Visit.Name.Norm, .keep_all = TRUE) %>%
  filter(Day.in.treatment < 60 & Day.in.treatment > 0 )


ggplot(vvtest3 %>% filter(TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose") &
                            collect.pre.post.rec %in% c("No.rec", "Pre.rec", "Post.rec")),
       aes(x=rec.diagnosis.Wk8 , y= (shannon_diversity),fill=rec.diagnosis.Wk8 ,
           color=rec.diagnosis.Wk8  )) +
  geom_boxplot(alpha=.6, outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +

  #scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values =  rec.cols.prepost) +
  scale_color_manual(values =  rec.cols.prepost) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14)) +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14)) +
  #  scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
  #                            "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD")) +
  labs(x="", y="Alpha diversity")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_AlphaDiversity_NO_TIME_RecWk8.pdf",
                                      sep = " ")), width = 6, height = 5)
 
##################################################################


vvtest4 <- vvtest3 %>% filter(Day.in.treatment < 20 )


ggplot(vvtest4 %>% filter(TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose") &
                            collect.pre.post.rec %in% c("No.rec", "Pre.rec", "Post.rec")),
       aes(x=rec.diagnosis.Wk8 , y= (shannon_diversity),fill=rec.diagnosis.Wk8 ,
           color=rec.diagnosis.Wk8  )) +
  geom_boxplot(alpha=.6, outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  
  #scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values =  rec.cols.prepost) +
  scale_color_manual(values =  rec.cols.prepost) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14)) +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14)) +
  #  scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
  #                            "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD")) +
  labs(x="", y="Alpha diversity")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_AlphaDiversity_NO_TIME_Week 1-2 RecWk8.pdf",
                                      sep = " ")), width = 6, height = 5)

##################################################################

vvtest3 <- ve303.dnameta.det.alt %>%
  distinct(Subject.Number, Visit.Name.Norm, .keep_all = TRUE) %>%
  filter(Day.in.treatment < 60 & Day.in.treatment > 0 )


ggplot(vvtest3 %>% filter(TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose") &
                            collect.pre.post.rec %in% c("No.rec", "Pre.rec" )),
       aes(x=collect.pre.post.rec , y= (shannon_diversity),fill=collect.pre.post.rec,
           color=collect.pre.post.rec )) +
  geom_boxplot(alpha=.6, outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +

  #scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = prepostrec.cols) +
  scale_color_manual(values = prepostrec.cols) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14)) +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14)) +
  #  scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
  #                            "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD")) +
  labs(x="", y="Alpha diversity")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_AlphaDiversity_NO_TIME_prerec_norec.pdf",
                                      sep = " ")), width = 6, height = 5)
# 
# ##################################################################

vvtest_mean <- ve303.dnameta.det.alt %>%
  distinct(Subject.Number, Visit.Name.Norm, .keep_all = TRUE) %>%
  filter(Day.in.treatment < 60 & Day.in.treatment > 0 ) %>%
  group_by(Subject.Number ) %>%
  mutate(mean_shannon_subject = mean(shannon_diversity)) %>%
  ungroup() %>%
  distinct(Subject.Number, .keep_all = TRUE)

ggplot(vvtest_mean %>% filter(TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose", "Placebo")),
       aes(x=rec.TRT , y= (mean_shannon_subject),fill=rec.diagnosis.Wk8,
           color=rec.diagnosis.Wk8 )) +
  geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  
  #scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values =  rec.cols.prepost) +
  scale_color_manual(values =  rec.cols.prepost) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, size = 14)) +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(size = 14,angle = 40 ), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14)) +
  scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
                            "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD",
                            "Placebo Good" = "Placebo", "Placebo Bad" = "Placebo")) +
  labs(x="", y="Alpha diversity (Mean per Subject)")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_Interaction_AlphaDiversity_NO_TIME_MEANTIME_prerec_norec_AllTreats.pdf",
                                      sep = " ")), width = 6, height = 5)
# ##################################################################


vvtest_mean_early1 <- ve303.dnameta.det.alt %>%
  distinct(Subject.Number, Visit.Name.Norm, .keep_all = TRUE) %>%
  filter(Day.in.treatment < 20 & Day.in.treatment > 0 ) %>%
  group_by(Subject.Number ) %>%
  mutate(mean_shannon_subject = mean(shannon_diversity)) %>%
  ungroup() %>%
  mutate(mean_shannon_subject_2wk = mean_shannon_subject) %>%
  distinct(Subject.Number, .keep_all = TRUE)

## Calculate delta recurrence - non-recurrence for each treatment group
vvtest_delta_calc <- vvtest_mean_early1 %>% 
  select(Subject.Number, TRT.norm, rec.diagnosis.Wk8, mean_shannon_subject_2wk, rec.TRT) %>%
  group_by(TRT.norm, rec.diagnosis.Wk8) %>%
  mutate(mean_shannon_treatment_recurrence = mean(mean_shannon_subject_2wk)) %>%
  mutate(median_shannon_treatment_recurrence = median(mean_shannon_subject_2wk)) %>%
  ungroup() %>%
  select(rec.TRT, rec.diagnosis.Wk8, TRT.norm, mean_shannon_treatment_recurrence, median_shannon_treatment_recurrence) %>%
  distinct(rec.TRT, .keep_all = TRUE) %>%
  arrange(rec.TRT) %>%
  group_by(TRT.norm) %>%
  mutate(Delta_GoodvsBad = median_shannon_treatment_recurrence - 
           lag(median_shannon_treatment_recurrence))  %>%
  ungroup() %>%
  filter(!is.na(Delta_GoodvsBad))

Delta_labels = paste0("Delta = ", signif(vvtest_delta_calc$Delta_GoodvsBad,1))

# A data frame with labels for each facet
f_labels <- data.frame(TRT.norm = c("Placebo", "VE303 High Dose", "VE303 Low Dose"), 
                       Delta_Annotate = Delta_labels)

vvtest_mean_early <- left_join(vvtest_mean_early1, f_labels, by = "TRT.norm")

vvtest_mean_early$TRT.norm = factor(vvtest_mean_early$TRT.norm, levels = 
                                      c("Placebo", "VE303 Low Dose" , "VE303 High Dose"))
  
shannon_rec <- ggplot(vvtest_mean_early , 
                  aes(x=rec.diagnosis.Wk8 , y= (mean_shannon_subject),fill=rec.diagnosis.Wk8)) +
  geom_boxplot( aes(color = rec.diagnosis.Wk8), alpha=0.5, outlier.size = 0, outlier.shape = NA ) +
  geom_point(  aes(color = rec.diagnosis.Wk8),pch = 21, size = 2, alpha = 0.3, 
               position = position_jitterdodge(), show.legend = FALSE) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  theme_bw() +
  scale_fill_manual(name = "Diagnosis", values =  rec.cols , 
                    labels = c("Recurrent", "Non-Recurrent")) +
  scale_color_manual(name = "Diagnosis",values =  rec.cols,
                     labels = c("Recurrent", "Non-Recurrent")) +
  geom_text(aes(x = 1.03, y = 5.5, label = Delta_Annotate) , 
            size = 4.5 , color = "darkslategray" ) +
  facet_grid(~TRT.norm) + 
  theme(axis.text.x = element_blank()) +
#  theme(axis.text.x = element_text(size = 14,angle = 40 , hjust = 1, vjust = 1 )) +
  
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        strip.text.x = element_text(size=14), strip.text.y = element_text(size=14)) +
  #scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
  #                          "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD",
  #                          "Placebo Good" = "Placebo", "Placebo Bad" = "Placebo")) +
  
  labs( title = "Diversity Post-Abx and During Dosing : Recurrent vs Non-recurrent Subjects" ,x="", 
        y="Shannon Index (Mean per Subject)", fill = "Diagnosis" )


ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "AlphaDiversity_MEANTIME_Week 1-2_recWk8_AllTreats.pdf",
                                      sep = " ")), width = 9, height = 4)


vvtest_mean_early$TRT.norm = factor(vvtest_mean_early$TRT.norm, levels = 
                                      c("VE303 High Dose", "VE303 Low Dose" ,"Placebo" ))
shannon_rec_MS <- ggplot(vvtest_mean_early , 
                      aes(x=rec.diagnosis.Wk8 , y= (mean_shannon_subject),fill=rec.diagnosis.Wk8)) +
  geom_boxplot( aes(color = rec.diagnosis.Wk8), alpha=0.5, outlier.size = 0, outlier.shape = NA ) +
  geom_point(  aes(color = rec.diagnosis.Wk8),pch = 21, size = 2, alpha = 0.3, 
               position = position_jitterdodge(), show.legend = FALSE) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  theme_bw() +
  #scale_fill_manual(name = "Diagnosis", values =  rec.cols , 
  #                  labels = c("Recurrent", "Non-Recurrent")) +
  #scale_color_manual(name = "Diagnosis",values =  rec.cols,
  #                   labels = c("Recurrent", "Non-Recurrent")) +
  geom_text(aes(x = 1.23, y = 5.5, label = Delta_Annotate) , 
            size = 4.5 , color = "darkslategray" ) +
  facet_grid(~TRT.norm) + 
  #theme(axis.text.x = element_blank()) +
  #theme(axis.text.x = element_text(size = 14,angle = 40 , hjust = 1, vjust = 1 )) +
  
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        strip.text.x = element_text(size=14), strip.text.y = element_text(size=14),
        legend.position = "none") +
  scale_x_discrete(labels=c("Bad" = "Recurrent", "Good" = "Non\nRecurrent")) +

  labs( title = "Diversity During Dosing" ,x="", 
        y="Shannon Index (Mean per Subject)", fill = "Diagnosis" )



shannon_rec_MS1 <- ggplot(vvtest_mean_early , 
                         aes(x=TRT.norm, y= (mean_shannon_subject),fill=rec.diagnosis.Wk8)) +
  geom_boxplot( aes(color = rec.diagnosis.Wk8), alpha=0.5, outlier.size = 0, outlier.shape = NA ) +
  geom_point(  aes(color = rec.diagnosis.Wk8),pch = 21, size = 2, alpha = 0.3, 
               position = position_jitterdodge(), show.legend = FALSE) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  theme_bw() +
  scale_fill_manual(name = "Diagnosis", values =  rec.cols , 
                    labels = c("Recurrent", "Non-Recurrent")) +
  scale_color_manual(name = "Diagnosis",values =  rec.cols,
                     labels = c("Recurrent", "Non-Recurrent")) +
  geom_text(aes(x = 0.93, y = 5.5, label = Delta_labels[1]) , 
            size = 4.5 , color = "darkslategray" ) +
  geom_text(aes(x = 1.93, y = 5.5, label = Delta_labels[3]) , 
            size = 4.5 , color = "darkslategray" ) +
  geom_text(aes(x = 2.93, y = 5.5, label = Delta_labels[2]) , 
            size = 4.5 , color = "darkslategray" ) +
  
  # geom_text(aes(x = 1.23, y = 5.5, label = Delta_Annotate) , 
  #            size = 4.5 , color = "darkslategray" ) +
  scale_x_discrete(labels = c("Placebo" = "Placebo" , "VE303 Low Dose" = "VE303\nLow Dose",
                              "VE303 High Dose" = "VE303\nHigh Dose" )) +
  
  #theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(size=14),axis.text.x = element_text(size=12),
        strip.text = element_text(size=14),
        strip.text.x = element_text(size=14), strip.text.y = element_text(size=14),
        legend.position = "none") +
  labs( title = "Diversity During Dosing" ,x="", 
        y="Shannon Index (Mean per Subject)", fill = "Diagnosis" )

shannon_rec_MS1

# ##################################################################
# ##################################################################
vvtest3 <- ve303.dnameta.det.alt %>%
  distinct(Subject.Number, Visit.Name.Norm, .keep_all = TRUE) %>%
  filter(Day.in.treatment < 60 & Day.in.treatment > 0 )


ggplot(vvtest3 %>% filter(TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose", "Placebo") &
                            collect.pre.post.rec %in% c("No.rec", "Pre.rec" )),
       aes(x=rec.TRT , y= (shannon_diversity),fill=collect.pre.post.rec.Wk8,
           color=collect.pre.post.rec.Wk8 )) +
  geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +

  #scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = prepostrec.cols) +
  scale_color_manual(values = prepostrec.cols) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 1, size = 14)) +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(size = 14,angle = 40 ), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14)) +
  scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
                            "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD",
                            "Placebo Good" = "Placebo", "Placebo Bad" = "Placebo")) +
  labs(x="", y="Alpha diversity")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "AlphaDiversity_NO_TIME_prerec_norec_AllTreats.pdf",
                                      sep = " ")), width = 6, height = 5)

# ##################################################################

ggplot(vvtest3 %>% filter(TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose" ) &
                            collect.pre.post.rec %in% c("No.rec", "Pre.rec" )),
       aes(x=rec.TRT , y= (shannon_diversity),fill=collect.pre.post.rec.Wk8,
           color=collect.pre.post.rec.Wk8 )) +
  geom_boxplot(alpha=.6, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +

  #scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = prepostrec.cols) +
  scale_color_manual(values = prepostrec.cols) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14)) +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14)) +
  scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
                            "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD")) +
  labs(x="", y="Alpha diversity")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_Interaction_AlphaDiversity_NO_TIME_prerec_norec_Dosed.pdf",
                                      sep = " ")), width = 6, height = 5)

##################################################################
##################################################################

ggplot(vvtest3 %>% filter(TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose", "Placebo" ) &
                            collect.pre.post.rec %in% c("No.rec", "Pre.rec" , "Post.rec")),
       aes(x=rec.TRT , y= (shannon_diversity),fill=rec.diagnosis.Wk8,
           color=rec.diagnosis.Wk8 )) +
  geom_boxplot(alpha=.6, outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +

  #scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = rec.cols) +
  scale_color_manual(values = rec.cols) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 14)) +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14)) +
  scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
                            "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD",
                            "Placebo Good" = "Placebo", "Placebo Bad" = "Placebo")) +
  labs(x="", y="Alpha diversity")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_Interaction_AlphaDiversity_NO_TIME_recWk8_AllDosed.pdf",
                                      sep = " ")), width = 6, height = 5)

# ##################################################################
 
ggplot(vvtest3 %>% filter(TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose", "Placebo" ) ),
aes(x=TRT.norm , y= (shannon_diversity),fill=TRT.norm,
    color=TRT.norm )) +
  geom_boxplot(alpha=.6, outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +

  #scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = coh.cols) +
  scale_color_manual(values = coh.cols) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14)) +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14)) +
  scale_x_discrete(labels=c("VE303 High Dose" = "HD", "Placebo" = "Placebo" , "VE303 Low Dose" = "LD")) +
  labs(x="", y="Alpha diversity")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_Treatment_AlphaDiversity_NO_TIME_AllDosed.pdf",
                                      sep = " ")), width = 6, height = 5)


##################################################################

ggplot(vvtest3 %>% filter(TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose", "Placebo" ) ),
aes(x=TRT.norm , y= (species_richness),fill=TRT.norm,
    color=TRT.norm )) +
  geom_boxplot(alpha=.6, outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +

  #scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = coh.cols) +
  scale_color_manual(values = coh.cols) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 14)) +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14)) +
  scale_x_discrete(labels=c("VE303 High Dose" = "HD", "Placebo" = "Placebo" , "VE303 Low Dose" = "LD")) +
  
  ylim(0,250)+
  labs(x="", y="Species richness")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_Treatment_Richness_NO_TIME_AllDosed.pdf",
                                      sep = " ")), width = 6, height = 5)
  
# ##################################################################
# 
ggplot(vvtest3 %>% filter(TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose", "Placebo" ) &
                            #  collect.pre.post.rec %in% c("Pre.rec" )
                            rec.diagnosis.Wk8 %in% c("Good", "Bad") ),
       aes(x=rec.TRT , y= (species_richness),fill=rec.diagnosis.Wk8,
           color=rec.diagnosis.Wk8 )) +
  geom_boxplot(alpha=.6, outlier.size = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +

  #scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values = rec.cols) +
  scale_color_manual(values = rec.cols) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 14)) +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14)) +
  scale_x_discrete(labels=c("VE303 High Dose Good" = "HD", "VE303 High Dose Bad" = "HD",
                            "VE303 Low Dose Good" = "LD", "VE303 Low Dose Bad" = "LD",
                            "Placebo Good" = "Placebo", "Placebo Bad" = "Placebo")) +
  ylim(0, 250)+
  labs(x="", y="Species Richness")

ggsave(filename = here(results, paste(Sys.Date(),taxplot,
                                      "ve303_colonization_Treatment_Richness_NO_TIME_RecWk8_AllDosed.pdf",
                                      sep = " ")), width = 6, height = 5)

# ##################################################################
vvtest5 <- ve303.dnameta.det.alt %>%
  distinct(Subject.Number, Visit.Name.Norm, .keep_all = TRUE) 
days_long <- c("Screening", "Day 1"  , "Day 7" , "Day 14" , "Day 28", "Day 56" , "Day 168")


s1 <- ggplot(vvtest5 %>%
         filter(.,Visit.Name%in% days_long),  
       aes(x=Visit.Name, y=shannon_diversity, fill=TRT.norm, 
           color=TRT.norm, group=Subject.Number)) +
#  geom_line(aes(group = Subject.Number), alpha = 0.2, size = 1) +
  stat_summary(aes(group=TRT.norm), fun = "median", geom = "line", alpha = 0.6 ,size = 1.8)+
  stat_summary(aes(group=TRT.norm), fun = "median", geom = "point",  alpha = 0.6, size = 3, show.legend = FALSE)+
  stat_summary(aes(group=TRT.norm),fun.data = median_q1q3, geom = "linerange")+
  facet_grid(~TRT.norm) +
  theme_bw() +
  scale_fill_manual(values = coh.cols) + 
  scale_color_manual(values = coh.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     
        axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1.5) ,
        strip.text = element_text(size = 13)) + 
  facet_grid(~TRT.norm, scales = "free") +
  labs( title = "Alpha Diversity Recovery per Treatment Group", y="Shannon Index", 
        x="Day in Treatment", color="Treatment") + 
  ggsave(filename = here(results, paste(Sys.Date(),taxplot, "Alpha Diversity_Median lineplot_Treatment.pdf", 
                                        sep = " ")), width = 9, height = 4, dpi = 800)

s1

ggplot(vvtest5 %>%
         filter(.,Visit.Name%in% days_long),  
       aes(x=Visit.Name, y=shannon_diversity, fill=TRT.norm, 
           color=TRT.norm, group=Subject.Number)) +
  #  geom_line(aes(group = Subject.Number), alpha = 0.2, size = 1) +
  stat_summary(aes(group=TRT.norm), fun = "mean", geom = "line", alpha = 0.6 , size = 1.8)+
  stat_summary(aes(group=TRT.norm), fun = "mean", geom = "point",  alpha = 0.6, size = 3, show.legend = FALSE)+
  stat_summary(aes(group=TRT.norm),fun.data = "mean_se", geom = "linerange")+
  facet_grid(~TRT.norm) +
  theme_bw() +
  scale_fill_manual(values = coh.cols) + 
  scale_color_manual(values = coh.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     
        axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1.5) ,
        strip.text = element_text(size = 13)) + 
  facet_grid(~TRT.norm, scales = "free") +
  labs( title = "A.  Alpha Diversity Recovery per Treatment Group", y="Shannon Index", 
        x="Day in Treatment", color="Treatment") + 
  ggsave(filename = here(results, paste(Sys.Date(),taxplot, "Alpha Diversity_Mean lineplot_Treatment.pdf", 
                                        sep = " ")), width = 9, height = 4, dpi = 800)



ggplot(vvtest5 %>%
               filter(.,Visit.Name%in% days_long),  
             aes(x=Visit.Name, y=shannon_diversity, fill=TRT.norm, 
                 color=TRT.norm )) +
  geom_boxplot(alpha=.6, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge()) +
  facet_grid(~TRT.norm) +
  theme_bw() +
  scale_fill_manual(values = coh.cols) + 
  scale_color_manual(values = coh.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     
        axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1.5) ,
        strip.text = element_text(size = 13)) + 
  facet_grid(~TRT.norm, scales = "free") +
  labs( title = "Alpha Diversity Recovery per Treatment Group", y="Shannon Index", 
        x="Day in Treatment", color="Treatment", fill="Treatment") + 
  ggsave(filename = here(results, paste(Sys.Date(),taxplot, "Alpha Diversity_Boxplot_Treatment.pdf", 
                                        sep = " ")), width = 9, height = 4, dpi = 800)



###################################################################
patch <- s1 / shannon_rec + plot_annotation(tag_levels = 'A')

ggsave(filename = here(results, paste(Sys.Date(),taxplot, "Alpha Diversity_Treatment and Recurrence.pdf", 
                                      sep = " ")), width = 9, height = 8, dpi = 800)


