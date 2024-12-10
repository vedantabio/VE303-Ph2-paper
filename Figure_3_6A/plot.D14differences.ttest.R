
#####################################

vv <- VE303_sp1.all %>% filter(Treatment %in% c("2", "3")) #%>% distinct(Subject.Number, .keep_all = TRUE)
vv$Treats  <- plyr::mapvalues(vv$Treatment, from= c("2","3"), to=c("VE303 Low Dose","VE303 High Dose"))
vv$Recs  <- plyr::mapvalues(vv$Response_Wk8, from= c("1","2"), to=c("Bad","Good"))


p <- ggplot( vv, aes(x=Treats, y= (0.01 + 100*rel_abund), color = Treats, fill = Treats )) + 
  geom_boxplot(alpha=.6,  outlier.shape = NA,width = 0.5,show.legend = TRUE ) +
  geom_point(pch = 21, size = 2.5, alpha = 0.3, position = position_jitterdodge(0.2) , show.legend = FALSE ) +
  theme_bw() +
  facet_wrap(organism~.) + 
  scale_y_log10() +
  #  scale_y_continuous(limits = c(0.0,30)) +
  scale_color_manual(values = coh.cols) +
  scale_fill_manual(values = coh.cols) +
  labs(title= "ve303 Strain Abundance at Day 14", y="Relative abundance", x="", 
       color="Treats", fill="Treats" ) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     
        axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14),  legend.text=element_text(size=12) , legend.title = element_blank()) + 
  
  ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 14 All_VE303_Strains_abundance_boxplot.pdf", sep = " ")),
         width = 10, height = 7, dpi = 800)

#####################################


p <- ggplot( vv, aes(x=Recs, y= (0.01 + 100*rel_abund),  color = Recs , fill = Recs )) + 
  geom_boxplot(alpha=.6,  outlier.shape = NA,width = 0.5,show.legend = TRUE ) +
  geom_point(pch = 21, size = 2.5, alpha = 0.3, position = position_jitterdodge(0.2) , show.legend = FALSE ) +
  theme_bw() +
  facet_wrap(organism~.) + 
  scale_y_log10() +
  #  scale_y_continuous(limits = c(0.0,30)) +
  scale_color_manual(values = rec.cols) +
  scale_fill_manual(values = rec.cols) +
  labs(title= "ve303 Strain Abundance at Day 14", y="Relative abundance", x="", 
       color="Recs", fill="Recs" ) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.text.y = element_text(size = 14, hjust = 1, vjust =0.5) ,     
        axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 14),  legend.text=element_text(size=12) , legend.title = element_blank()) + 
  
  ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 14 All_VE303_Strains_Recurrence_abundance_boxplot.pdf", sep = " ")),
         width = 10, height = 7, dpi = 800)


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

days_14 <- c("Day 14")
textsize <- 14
ggplot(ve303.abund.subj.total.time.det.alt  %>%
         filter(.,Visit.Name%in% days_14), 
       aes(x=TRT.norm, y=((est_relative_abundance_panel_alt_sum*100)), fill=TRT.norm, color=TRT.norm)) +
  geom_boxplot(alpha=.6,  outlier.shape = NA,width = 0.5,show.legend = FALSE ) +
  geom_point(pch = 21, size = 2.5, alpha = 0.3, position = position_jitterdodge() , show.legend = FALSE ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  #facet_grid(~TRT.norm, scales = "free_x", space = "free") +
  scale_fill_manual(values = coh.cols) +
  scale_color_manual(values = coh.cols) +
  theme(axis.text.x = element_text(size = textsize, angle = 45, hjust = 1, vjust = 1), 
        axis.title.x = element_text(size = textsize ),
        plot.title = element_text(size = textsize ),
        axis.text.y = element_text(size = textsize, hjust = 1, vjust =0.5) ,     
        axis.title.y = element_text(size = textsize, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = textsize)) +
  #scale_y_log10() +
  labs(title= "ve303 Total Abundance at Day 14", y="Relative Abundance (Total)", x="", 
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
  theme(axis.text.x = element_text(size = textsize, angle = 45, hjust = 1, vjust = 1), 
        axis.title.x = element_text(size = textsize ),
        plot.title = element_text(size = textsize ),
        axis.text.y = element_text(size = textsize, hjust = 1, vjust =0.5) ,     
        axis.title.y = element_text(size = textsize, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = textsize)) +
  #scale_y_log10() +
  labs(title= "ve303 Total Abundance at Day 14", y="Relative Abundance (Total)", x="", 
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
  theme(axis.text.x = element_text(size = textsize, angle = 45, hjust = 1, vjust = 1), 
        axis.title.x = element_text(size = textsize ),
        plot.title = element_text(size = textsize ),
        axis.text.y = element_text(size = textsize, hjust = 1, vjust =0.5) ,     
        axis.title.y = element_text(size = textsize, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = textsize)) +
  scale_y_log10( ) +
  labs(title= "ve303 Total Abundance at Day 14", y="Relative Abundance (Total)", x="", 
       color="TRT.norm", fill="TRT.norm" )
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 14 only HD LD LOG ve303_total_abundance_boxplot_det.pdf", sep = " ")), 
       width = 6, height = 6, dpi = 800 )

####################################################################################
# Reorder visit names
ve303.colonize.sum.alt$Visit.Name.Norm = factor(ve303.colonize.sum.alt$Visit.Name.Norm, 
                                                levels = c("Day -11","Screening", "Day 1", "Day 5", "Day 7", "Day 10", 
                                                           "Day 12", "Day 13", "Day 14", "Day 15","Day 18", "Day 19", "Day 20", 
                                                           "Day 21", "Day 26","Day 28", "Day 45", "Day 47","Day 56", 
                                                           "Day 70", "Day 87", "Day 114", "Day 120", "Day 148", "Day 168"))
# Reorder treatmemt arms
ve303.colonize.sum.alt$TRT.norm = factor(ve303.colonize.sum.alt$TRT.norm,
                                         levels = c("Placebo", "VE303 Low Dose", "VE303 High Dose"))


days_14 <- c("Day 14")
textsize <- 14 

ggplot(ve303.colonize.sum.alt %>%
         filter(.,Visit.Name.Norm %in% c( "Day 14" )), 
       aes(x=TRT.norm, y=N_det_alt, fill=TRT.norm, color=TRT.norm)) + 
  geom_boxplot(alpha=.6,  outlier.shape = NA,width = 0.5,show.legend = FALSE ) +
  geom_point(pch = 21, size = 2.5, alpha = 0.3, position = position_jitterdodge() , show.legend = FALSE ) +
  theme_bw() +
  scale_fill_manual(values = coh.cols) + 
  scale_color_manual(values = coh.cols) + 
  #scale_y_continuous(breaks = c(0:8), limits = c(0,8))
  
  theme(axis.text.x = element_text(size = textsize, angle = 45, hjust = 1, vjust = 1), 
        axis.title.x = element_text(size = textsize ),
        plot.title = element_text(size = textsize ),
        axis.text.y = element_text(size = textsize, hjust = 1, vjust =0.5) ,     
        axis.title.y = element_text(size = textsize, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = textsize)) +
  labs(title= "ve303 Total Detection at Day 14", y="VE303 strains detected (#)", x="", 
       color="TRT.norm", fill="TRT.norm" )
ggsave(filename = here(results, paste(Sys.Date(), taxplot,"Day 14 only ve303_total_Proportion_boxplot_det.pdf", sep = " ")), 
       width = 6, height = 6, dpi = 800)
