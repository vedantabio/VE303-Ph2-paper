
rec.cols.prepost <- c("Bad" = "forestgreen", "Good" = "#561600")
rec.cols <- c("Bad" = "indianred2", "Good" = "cyan3")

abx.cols <- c("Vancomycin" =  "#0487e3", "Fidaxomicin" = "#dc2d00" , 
              "Metronidazole" = "#5b0076" , "Other" = "forestgreen")

abx.cols1 <- c("On-Abx Screening" ="#dc2d00" , "True Screening" = "forestgreen"  )

abx_diversity$Visit.Name = factor(abx_diversity$Visit.Name, levels = c("Screening", "Day 1"  , "Day 7" , "Day 14" , 
                                                 "Day 28", "Day 56" , "Day 168"))
abx_diversity$TRT.norm = factor(abx_diversity$TRT.norm, levels = c("Placebo", "VE303 Low Dose"  , "VE303 High Dose"))

abx_diversity$Abx.group = factor(abx_diversity$Abx.group, levels = c("Vancomycin", "Fidaxomicin"  , 
                                                                     "Metronidazole", "Other"))


abx_indi_diversity$Visit.Name = factor(abx_indi_diversity$Visit.Name, levels = c("Screening", "Day 1"  , "Day 7" , "Day 14" , 
                                                                       "Day 28", "Day 56" , "Day 168"))
abx_indi_diversity$TRT.norm = factor(abx_indi_diversity$TRT.norm, levels = c("Placebo", "VE303 Low Dose"  , "VE303 High Dose"))

abx_indi_diversity$Abx.group = factor(abx_indi_diversity$Abx.group, levels = c("Vancomycin", "Fidaxomicin"  , 
                                                                     "Metronidazole", "Other"))

days <- c("Screening", "Day 1"  , "Day 7", "Day 14" , "Day 28", "Day 56" , "Day 168")
days_56 <- c("Screening", "Day 1"  , "Day 7", "Day 14" , "Day 28", "Day 56"  )
days_pre <- c("Screening", "Day 1"   )
 

###################################################################
# PLOT INDIVIDUAL VE303 STRAIN AVERAGE EXPOSURE WITH FIDAXOMICIN AND VANCOMYCIN
###################################################################

abx.cols2 <- c("Vancomycin" =  "turquoise4", "Fidaxomicin" = "tan4" )

av_ex <- ggplot(abx_strain_exposure %>% filter(TRT.norm != "Placebo" ), 
       aes(x=organism , y= 0.01+(100*exposure_perstrain),
           fill=Abx.group, color=Abx.group)) +
  geom_boxplot(width = 0.5,alpha=.6, outlier.size = 0, outlier.shape = NA) +
  geom_point(pch = 21, size = 2, alpha = 0.2, position = position_jitterdodge()) +
  #geom_hline(yintercept = 0.1, linetype="longdash", color="red") +
  facet_grid(TRT.norm~., scales="free_x", space = "free") + 
  scale_y_log10()+
  theme_bw() +
  scale_fill_manual(values =  abx.cols2) +
  scale_color_manual(values =  abx.cols2) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1, size = 14)) +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14)) +
  labs(x="", y="Log Mean Exposure per Subject (1 month)")

ggsave(filename = here(results_abx, paste(Sys.Date(), taxplot,"VE303_Strains_AV EXPOSURE_Fidaxo_Vanco.pdf", sep = " ")), 
       width = 10, height = 6)
 

###################################################################

 
