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
    

############################# Figure 1B #######################################################

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


#### Plot ve303 total abundance lineplots ----
# Reorder visit names

ve303.abund.total303.subj.det.alt$Visit.Name = factor(ve303.abund.total303.subj.det.alt$Visit.Name,
                                                      levels = c("Screening", "Day 1"  , "Day 7" ,
                                                                 "Day 14" , "Day 28", "Day 56" , "Day 168"))


ve303.abund.total303.subj.det.alt$TRT.norm = factor(ve303.abund.total303.subj.det.alt$TRT.norm,
                                                    levels = c("Placebo", "VE303 Low Dose"  , "VE303 High Dose"))
####################################################################################
 
#### Plot ve303 relative abundance by Strain ----

#theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
#      axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
#      strip.text = element_text(size = 14)) +
  
# Reorder visit names
ve303.abund.sum.subj.det.alt$Visit.Name = factor(ve303.abund.sum.subj.det.alt$Visit.Name,
                                                 c("Screening", "Day 1"  , "Day 7" ,
                                                   "Day 14" , "Day 28", "Day 56" , "Day 168"))

ve303.abund.sum.subj.det.alt$subject_mean_pseudocount <- (0.01 + 100*ve303.abund.sum.subj.det.alt$subject_mean)

ve303.abund.sum.subj.det.alt$TRT.norm <- factor(ve303.abund.sum.subj.det.alt$TRT.norm,
                                                c("Placebo", "VE303 Low Dose" , "VE303 High Dose"))

################################ Figure 1C  ####################################################
 

abuntime_strains <- ggplot(ve303.abund.sum.subj.det.alt %>% filter(Visit.Name  %in% days),
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
#ggsave(filename = here(results, paste(Sys.Date(),taxplot, "ve303_colonization_meanmedianabund_det_perstrain_NO SUBJ_logscale.pdf", sep = " ")), width = 12, height = 12)


abuntime_strains_long <- ggplot(ve303.abund.sum.subj.det.alt %>% filter(Visit.Name  %in% days_long),
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
#ggsave(filename = here(results, paste(Sys.Date(),taxplot, "D168 ve303_colonization_meanmedianabund_det_perstrain_NO SUBJ_logscale.pdf", sep = " ")), width = 12, height = 12)

#############################################################################################################


ve303.dnameta.det.alt$Visit.Name = factor(ve303.dnameta.det.alt$Visit.Name,
                                            levels = c("Screening", "Day 1"  , "Day 7" ,
                                                       "Day 14" , "Day 28", "Day 56" , "Day 168"))


ve303.dnameta.det.alt$TRT.norm = factor(ve303.dnameta.det.alt$TRT.norm,
                                       levels = c("Placebo", "VE303 Low Dose"  , "VE303 High Dose"))

#############################################################################################################
#############################################################################################################
     

days2 <- c("Day 7"  ,   "Day 14" )
days1 <- c("Day 7"  ,   "Day 14" ) 
 
##################################################################
   
vvtest3 <- ve303.dnameta.det.alt %>%
  distinct(Subject.Number, Visit.Name.Norm, .keep_all = TRUE) %>%
  filter(Day.in.treatment < 60 & Day.in.treatment > 0 )


##################################################################


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
                                      c("VE303 High Dose", "VE303 Low Dose" ,"Placebo" ))


# ##################################################################
# ##################################################################
vvtest3 <- ve303.dnameta.det.alt %>%
  distinct(Subject.Number, Visit.Name.Norm, .keep_all = TRUE) %>%
  filter(Day.in.treatment < 60 & Day.in.treatment > 0 )

  
  
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
        x="Day in Treatment", color="Treatment") 
  ggsave(filename = here(results, paste(Sys.Date(),taxplot, "Alpha Div_lineplot_Trt_Figure5A.pdf", 
                                        sep = " ")), width = 9, height = 4, dpi = 800)

s1



