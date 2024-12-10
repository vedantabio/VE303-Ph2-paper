## Age and recurrence plots ----

#######################

# Reorder treatment arms
dnameta$TRT.norm = factor(dnameta$TRT.norm,
                                         levels = c("Placebo", "VE303 Low Dose", "VE303 High Dose"))


  ############# Extended data figure 5B #################################
  
  age_violin_EXTfig5b <- ggplot(dnameta %>% filter(TRT.norm != "Placebo") %>%
           distinct(Subject.Number, .keep_all = TRUE ), 
         aes(x=rec.diagnosis.Wk8, y= Age, fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8)) + 
    geom_boxplot(alpha=0.6, width = 0.45,show.legend = FALSE, fill = "white",color="black", 
                 outlier.size = 0,  outlier.shape = NA) +
    geom_point(pch = 21, size = 5,  alpha=0.5, 
               position = position_jitterdodge(), show.legend = FALSE) +
    
    geom_flat_violin(position = position_nudge(x = .25),alpha = 0.7, adjust = 0.5, show.legend = FALSE) +
    geom_hline(yintercept = 60, linetype="longdash", color="red") +
    
    theme_bw() + 
    
    scale_fill_manual(values = rec.cols) + 
    scale_color_manual(values = rec.cols) + 
    labs(title= " ", y="Age", x=" " ) + 
    theme(axis.text.x = element_text(size = 14, angle = 0, hjust = 0.5, vjust = 0.5), axis.title.x = element_text(size = 14),
          axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
          strip.text = element_text(size = 14))  +
    scale_x_discrete(labels = c("Bad" = "Recurrent" , "Good" = "Non\nRecurrent"  )) 
  ggsave(filename = here(results, paste(Sys.Date(), taxplot,"EXT Fig5B_Age_violin_plot Dosed.pdf", sep = " ")),
         width = 8.5, height = 6)
  
  #######################
  #######################
  
## Detection summary plots ----

## ve303.colonize.sum.alt counts detections as 1 and all else as zero in the column N_det_alt USE THIS

# Reorder visit names
ve303.colonize.sum.alt$Visit.Name.Norm = factor(ve303.colonize.sum.alt$Visit.Name.Norm, 
                                                         levels = c("Day -11","Screening", "Day 1", "Day 5", "Day 7", "Day 10", 
                                                                    "Day 12", "Day 13", "Day 14", "Day 15","Day 18", "Day 19", "Day 20", 
                                                                    "Day 21", "Day 26","Day 28", "Day 45", "Day 47","Day 56", 
                                                                    "Day 70", "Day 87", "Day 114", "Day 120", "Day 148", "Day 168"))
# Reorder treatmemt arms
ve303.colonize.sum.alt$TRT.norm = factor(ve303.colonize.sum.alt$TRT.norm,
                                                  levels = c("Placebo", "VE303 Low Dose", "VE303 High Dose"))
 
############################ Figure 1A ########################################################

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

####################################################################################

## VE303 strains brickplot: Extended data Figure 9 ###

ve303.dnameta$Visit.Name.Norm = factor(ve303.dnameta$Visit.Name.Norm, 
                                       levels = c("Screening", "Day 1", "Day 7",  
                                                  "Day 14",  
                                                  "Day 28", "Day 56",  "Day 168"))

brickplot_dnameta <- ve303.dnameta %>% 
  filter(TRT.norm != "Placebo" & Visit.Name.Norm %in%  c("Day 14" )) %>%
  arrange(rec.diagnosis.Wk8)

brickplot_dnameta$rec.diagnosis.Wk8= factor(brickplot_dnameta$rec.diagnosis.Wk8, 
                                            levels = c("Bad", "Good"))


ggplot( brickplot_dnameta, 
        aes(x=organism, y=freq_alt, fill=organism)) + 
  geom_bar(stat = "identity", show.legend = T) + 
  scale_y_discrete() + 
  theme_bw() + 
  facet_grid(  rec.diagnosis.Wk8+Subject.Number~organism, scales = "free", space ="free") + 
  scale_fill_manual(values = strain.cols) + 
  theme(axis.text.y = element_blank() , axis.title.x = element_text(size = 16),
        axis.text.x = element_blank(), axis.title.y = element_text(hjust = 0.5, vjust = 0.5, size = 16),
        strip.text.x = element_text(size=14) ) + 
  labs(title= element_blank(), y="Subject", x="VE303 Strain Detection", fill="VE303 Strain") 
ggsave(filename = here(results, paste(Sys.Date(), 
                                      taxplot,"Dosed_ALL_VE303_Strain_Brickplot_EXT Figure9.png", sep = " ")), 
       width = 12, height = 12, dpi = 650)

####################################################################################

## Stacked Barplots per subject ----
## Extended data figure 2A,B,C

# High dose
ve303.dnameta$Visit.Name.Norm = factor(ve303.dnameta$Visit.Name.Norm, 
                                       levels = c("Screening", "Day 1", "Day 7",  
                                                  "Day 14",  
                                                  "Day 28", "Day 56",  "Day 168"))

stack_hd_EXTfig2A <-  ggplot( ve303.dnameta %>% filter(TRT.norm == "VE303 High Dose" & Visit.Name.Norm %in% 
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
                                       taxplot,"EXT Fig2A_HD_VE303_StackedBar.pdf", sep = " ")), 
        width = 12, height = 12)
 
 #############
 
stack_ld_EXTfig2B <- ggplot(ve303.dnameta %>% filter(TRT.norm == "VE303 Low Dose" & Visit.Name.Norm %in% 
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
 ggsave(filename = here(results, paste(Sys.Date(), taxplot,"EXT Fig2B_LD_VE303_StackedBar_detection.pdf", 
                                       sep = " ")), width = 12, height = 12)
 
 #############
 
stack_placebo_EXTfig2C <- ggplot(ve303.dnameta %>% filter(TRT.norm == "Placebo" & Visit.Name.Norm %in% 
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
ggsave(filename = here(results, paste(Sys.Date(),taxplot, "EXT Fig2C_Placebo_VE303_StackedBar_detection.pdf", 
                                      sep = " ")), width = 12, height = 12)
 

