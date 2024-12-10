## Detection summary plots ----

## ve303.colonize.sum.alt counts detections as 1 and all else as zero in the column N_det_alt USE THIS
#### TO PLOT DETECTIONS WTHOUT FILTERING STATUS DETECTED

## ve303.colonize.sum counts all non-zero abundance as detections, if used to plot detections it would need to be filter
#### detection_status = Detcted. But this is incorrect, these non-detections should be included and treated as zeros

# Reorder visit names
ve303.colonize.sum.alt$Visit.Name.Norm = factor(ve303.colonize.sum.alt$Visit.Name.Norm, 
                                                         levels = c("Day -11","Screening", "Day 1", "Day 5", "Day 7", "Day 10", 
                                                                    "Day 12", "Day 13", "Day 14", "Day 15","Day 18", "Day 19", "Day 20", 
                                                                    "Day 21", "Day 26","Day 28", "Day 45", "Day 47","Day 56", 
                                                                    "Day 70", "Day 87", "Day 114", "Day 120", "Day 148", "Day 168"))
# Reorder treatmemt arms
ve303.colonize.sum.alt$TRT.norm = factor(ve303.colonize.sum.alt$TRT.norm,
                                                  levels = c("Placebo", "VE303 Low Dose", "VE303 High Dose"))

ggplot(ve303.colonize.sum.alt %>%
#  filter(detection_status =="Detected") %>%
  filter(.,Visit.Name.Norm %in% c("Screening", "Day 1", "Day 7","Day 14", "Day 28", "Day 56"  )), 
  aes(x=Visit.Name.Norm, y=N_det_alt, fill=TRT.norm, color=TRT.norm)) + 
  geom_boxplot(alpha=0.6, show.legend = FALSE, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha=0.3, position = position_jitterdodge(), show.legend = FALSE) +
  theme_bw() + 
  facet_wrap(~TRT.norm) + 
  labs(title= "VE303 Colonization Summary", y="VE303 strains detected (#)", x="Timepoint", fill="Cohort") + 
  scale_fill_manual(values = coh.cols) + 
  scale_color_manual(values = coh.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) + 
  scale_y_continuous(breaks = c(0:8), limits = c(0,8))
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"All_VE303_colonization_detected_count_boxplot.pdf", sep = " ")),
       width = 10, height = 5)

####################################################################################

dettime <- ggplot(ve303.colonize.sum.alt %>%
#         filter(detection_status =="Detected") %>%
         filter(.,Visit.Name.Norm %in% c("Screening", "Day 1", "Day 7","Day 14", "Day 28", "Day 56", "Day 168" )), 
       aes(x=Visit.Name.Norm, y=N_det_alt, fill=TRT.norm, color=TRT.norm , show.legend = FALSE)) + 
  geom_boxplot(alpha=0.6,  outlier.size = 0, outlier.shape = NA, show.legend = TRUE) +
  geom_point(pch = 21, size = 2, alpha=0.3, position = position_jitterdodge() , show.legend = FALSE) +
  theme_bw() + 
  facet_wrap(~TRT.norm) + 
  scale_fill_manual(values = coh.cols) + 
  scale_color_manual(values = coh.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14)) + 
  scale_y_continuous(breaks = c(0:8), limits = c(0,8)) +
  labs(title= "VE303 Colonization Summary", y="VE303 strains detected (#)", 
       x="Timepoint" , color = "Treatment", fill = "Treatment") #+ guides(color = guide_legend(show = FALSE) )

ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 168 All_VE303_colonization_detected_count_boxplot.pdf", sep = " ")),
       width = 10, height = 5)


####################################################################################

days_14 <- c("Day 14")
textsize <- 14
#  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
ve303.colonize.sum.alt$TRT.norm = factor(ve303.colonize.sum.alt$TRT.norm,
                                                      levels = c("VE303 High Dose", "VE303 Low Dose"  ,"Placebo" ))

d14_prop <- ggplot(ve303.colonize.sum.alt %>%
         filter(.,Visit.Name.Norm %in% c( "Day 14" )), 
       aes(x=TRT.norm, y=N_det_alt, fill=TRT.norm, color=TRT.norm)) + 
  geom_boxplot(alpha=.6,  outlier.shape = NA,width = 0.5,show.legend = FALSE ) +
  geom_point(pch = 21, size = 2.5, alpha = 0.3, position = position_jitterdodge() , show.legend = FALSE ) +
  theme_bw() +
  scale_fill_manual(values = coh.cols) + 
  scale_color_manual(values = coh.cols) + 
  #scale_y_continuous(breaks = c(0:8), limits = c(0,8))
  
  theme(axis.text.x = element_text(size = textsize, angle = 0, hjust = 0.5, vjust = 0.5), 
        axis.title.x = element_text(size = textsize ),
        plot.title = element_text(size = textsize ),
        axis.text.y = element_text(size = textsize, hjust = 1, vjust =0.5) ,     
        axis.title.y = element_text(size = textsize, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = textsize)) +
  scale_x_discrete(labels = c("Placebo" = "Placebo" , "VE303 Low Dose" = "VE303\nLow Dose","VE303 High Dose" = "VE303\nHigh Dose" )) +
  labs(title= "VE303 Total Detection at Day 14", y="VE303 strains detected (#)", x="", 
       color="TRT.norm", fill="TRT.norm" )
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 14 only ve303_total_Proportion_boxplot_det.pdf", sep = " ")), 
       width = 6, height = 6, dpi = 800)


####################################################################################

ggplot(ve303.colonize.sum.alt %>%
#  filter(detection_status =="Detected") %>%
  filter(.,Visit.Name.Norm %in% c("Screening", "Day 1", "Day 7","Day 14", "Day 28", "Day 56" )), 
#  aes(x=Visit.Name.Norm, y=N_det, fill=rec.diagnosis, color=rec.diagnosis)) + 
  aes(x=Visit.Name.Norm, y=N_det_alt, fill=rec.diagnosis, color=rec.diagnosis)) + 
  geom_boxplot(alpha=0.6, outlier.size = 0, outlier.shape = NA, outlier.alpha = 1) + 
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge(), show.legend = FALSE) +
  theme_bw() + 
  facet_grid(TRT.norm~.) + 
  labs(title= "VE303 Colonization Summary", y="VE303 strains detected (#)", x="Timepoint", fill="Recurrence Diagnosis") + 
  scale_fill_manual(values = rec.cols) + 
  scale_color_manual(values = rec.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 13)) + 
  scale_y_continuous(breaks = c(0:8), limits = c(0,8))
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"All_VE303_colonization_detected_count_boxplot_by_recurrenceSens2.pdf", sep = " "))
       , width = 10, height = 8)


####################################################################################

ggplot(ve303.colonize.sum.alt %>%
#         filter(detection_status =="Detected") %>%
         filter(.,Visit.Name.Norm %in% c("Screening", "Day 1", "Day 7","Day 14", "Day 28", "Day 56" )), 
       #  aes(x=Visit.Name.Norm, y=N_det, fill=rec.diagnosis, color=rec.diagnosis)) + 
       aes(x=Visit.Name.Norm, y=N_det_alt, fill=rec.diagnosis, color=rec.diagnosis)) + 
  geom_boxplot(alpha=0.6, outlier.size = 0, outlier.shape = NA, outlier.alpha = 1) + 
  geom_point(pch = 21, size = 3, position = position_jitterdodge(), show.legend = FALSE) +
  theme_bw() + 
  facet_grid(TRT.norm~.) + 
  labs(title= "VE303 Colonization Summary", y="VE303 strains detected (#)", x="Timepoint", fill="Recurrence Diagnosis") + 
  scale_fill_manual(values = rec.cols) + 
  scale_color_manual(values = rec.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 13)) + 
  scale_y_continuous(breaks = c(0:8), limits = c(0,8))
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Long All_VE303_colonization_detected_count_boxplot_by_recurrenceSens2.pdf", sep = " "))
       , width = 10, height = 8)

####################################################################################

## Area plots per subject ----

# High dose
ve303.dnameta$Visit.Name.Norm = factor(ve303.dnameta$Visit.Name.Norm, 
                                       levels = c("Screening", "Day 1", "Day 7",  
                                                  "Day 14",  
                                                  "Day 28", "Day 56",  "Day 168"))

 ggplot( ve303.dnameta %>% filter(TRT.norm == "VE303 High Dose" & Visit.Name.Norm %in% 
                                    c("Screening", "Day 1", "Day 7","Day 14", "Day 28", 
                                      "Day 56" , "Day 168" )), 
        aes(x=Visit.Name.Norm, y=freq_alt, fill=organism)) + 
   geom_bar(stat = "identity", show.legend = T) + 
   theme_bw() + 
   facet_wrap(~Subject.Number, scales = "free" ) + 
   #facet_grid(~Subject.Number, scales = "free", space ="free") + 
   scale_fill_manual(values = strain.cols) + 
   labs(title= element_blank(), y="VE303 strains detected (#)", x="Timepoint", fill="VE303 Strain") + 
   theme(axis.text.y = element_text(size=10), strip.text = element_text(size=10), 
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), 
         strip.text.x = element_text(size=10), strip.text.y = element_text(size=10)) + 
   scale_y_continuous(limits = c(0,8)) 
   #theme(strip.text.x = element_blank())
 ggsave(filename = here(results, paste(Sys.Date(), 
                                       taxplot,"HighDose_VE303_StackedBar_detection.pdf", sep = " ")), 
        width = 12, height = 12)
 
 #############
 
ggplot(ve303.dnameta %>% filter(TRT.norm == "VE303 Low Dose" & Visit.Name.Norm %in% 
                                  c("Screening", "Day 1", "Day 7","Day 14", "Day 28", 
                                    "Day 56" , "Day 168" )), aes(x=Visit.Name.Norm, y=freq_alt, fill=organism)) + 
   geom_bar(stat = "identity", show.legend = T) + 
   theme_bw() + 
   facet_wrap(~Subject.Number, scales = "free" ) +  
   scale_fill_manual(values = strain.cols) + 
   labs(title= element_blank(), y="VE303 strains detected (#)", x="Timepoint", fill="VE303 Strain") + 
   theme(axis.text.y = element_text(size=10), strip.text = element_text(size=10), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), 
        strip.text.x = element_text(size=10), strip.text.y = element_text(size=10)) + 
  scale_y_continuous(limits = c(0,8)) 
   #theme(strip.text.x = element_blank())
 ggsave(filename = here(results, paste(Sys.Date(), taxplot,"LowDose_VE303_StackedBar_detection.pdf", 
                                       sep = " ")), width = 12, height = 12)
 
 #############
 
ggplot(ve303.dnameta %>% filter(TRT.norm == "Placebo" & Visit.Name.Norm %in% 
                                  c("Screening", "Day 1", "Day 7","Day 14", "Day 28", 
                                    "Day 56" , "Day 168" )), 
       aes(x=Visit.Name.Norm, y=freq_alt, fill=organism)) +
  geom_bar(stat = "identity", show.legend = T) +
  theme_bw() +
  facet_wrap(~Subject.Number, scales = "free" ) +
  scale_fill_manual(values = strain.cols) +
  labs(title= element_blank(), y="VE303 strains detected (#)", x="Timepoint", fill="VE303 Strain") +
   theme(axis.text.y = element_text(size=10), strip.text = element_text(size=10), 
         axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7), 
         strip.text.x = element_text(size=10), strip.text.y = element_text(size=10)) + 
  scale_y_continuous(limits = c(0,8))
  #theme(strip.text.x = element_blank())
ggsave(filename = here(results, paste(Sys.Date(),taxplot, "Placebo_VE303_StackedBar_detection.pdf", 
                                      sep = " ")), width = 12, height = 12)
 


## Total plots by strain ----

############### THIS PLOT EXCLDES DETECTED STRAINS AND DOES NOT ZERO THEM #######

#Reorder visit names
ve303.colonize.sum.bystrain.tp$Visit.Name.Norm = factor(ve303.colonize.sum.bystrain.tp$Visit.Name.Norm, 
                                                      levels = c("Day -12","Screening", "Day 1", "Day 7", "Day 9", 
                                                                 "Day 12", "Day 14", "Day 18", "Day 25", "Day 28", 
                                                                 "Day 39", "Day 56", "Day 69", "Day 74", "Day 113", 
                                                                 "Day 119", "Day 147", "Day 168"))
# Reorder treatmemt arms
ve303.colonize.sum.bystrain.tp$TRT.norm = factor(ve303.colonize.sum.bystrain.tp$TRT.norm,
                                               levels = c("Placebo", "VE303 Low Dose", "VE303 High Dose"))


ggplot(ve303.colonize.sum.bystrain.tp %>%
  filter(detection_status =="Detected") %>% filter(.,Visit.Name.Norm %in% c("Screening", "Day 1", "Day 7", "Day 14", "Day 28", "Day 56")), 
  aes(x=Visit.Name.Norm, y=N_det_tot, group=TRT.norm)) +
  geom_line(show.legend=FALSE, color="slateblue", lwd = 1) +
  geom_point()+
  theme_bw() +
  facet_grid(organism~TRT.norm, scales = "free_x") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7)) +
  labs(x="Timepoint", y="VE303 Prevalence (# Subjects)")
#  scale_y_continuous(limits = c(0,9), breaks = c(0,2,4,6,8))
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"VE303_colonization_prev_det_perstrain.pdf", sep = " ")), width = 10, height = 10)
