# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#Add frequency variable for all rows - this gives every detection/non-det result a value of 1
ve303.dnameta <- data.frame(ve303.dnameta_plot, freq = 1)

# Change from character class to factor
ve303.dnameta$detection_status <- 
  factor(ve303.dnameta$detection_status, levels = c("Detected", "Insufficient data", "Not detected"))


#######################################################################

# To plot zeroes in colonization data
ve303.dnameta = ve303.dnameta %>%
  mutate(freq_alt = if_else(detection_status != "Detected", 0, freq)) %>%
  group_by(Subject.Number, Visit.Name.Norm, organism) %>%
  mutate(.,freq_sum = sum(freq_alt)) %>%
  ungroup()


#######################################################################


## Calculate detection summary stats
# Calculate number of patients in each cohort
cohort.count <- ve303.dnameta %>% 
  distinct(Subject.Number, .keep_all = TRUE) %>% 
  group_by(TRT.norm) %>% 
  summarise(N_sub = n())

# Calculate number of patients in each cohort by visit.name.norm
cohort.timepoint.count <- ve303.dnameta %>% 
  group_by(TRT.norm, Visit.Name.Norm) %>% 
  distinct(Subject.Number, .keep_all = TRUE) %>% 
  summarise(N_analyzed = n())


# Calculate number of patients in each cohort by visit.name.norm
cohort.timepoint.count.num <- ve303.dnameta %>% 
  group_by(TRT.norm, Day.in.treatment.norm,Visit.Name.Norm) %>% 
  distinct(Subject.Number, .keep_all = TRUE) %>% 
  summarise(N_analyzed = n())

# Calculate number of patients in each cohort by visit.name.norm
cohort.rec.count <- ve303.dnameta %>% 
  distinct(Subject.Number, .keep_all = TRUE) 


#######################################################################
## Plot VE303 detection and calculate summary stats  - All DETECTED events

ve303.colonize.sum <- ve303.dnameta %>% 
  group_by(MGS.Sample.ID,Visit.Name.Norm, detection_status) %>% 
  summarise(N_det = n()) %>% 
  ungroup() %>%
  select(-Visit.Name.Norm) %>%
  left_join(.,dnameta, by = "MGS.Sample.ID") %>%
  left_join(., cohort.count, by = "TRT.norm")

##To plot zeros in "All_VE303_colonization_detected_count_boxplot.pdf"
ve303.colonize.sum$N_det = as.double(ve303.colonize.sum$N_det)

ve303.colonize.sum.alt = ve303.colonize.sum %>%
  mutate(N_det_alt = if_else(detection_status != "Detected", 0, N_det)) %>%
  group_by(Subject.Number, Visit.Name, MGS.Sample.ID) %>%
  mutate(.,N_det_alt_sum = sum(N_det_alt)) %>%
  ungroup() %>%
  mutate(prop_det_alt_sum = N_det_alt_sum/8) %>%
  distinct(.,MGS.Sample.ID, .keep_all = T)


ve303.colonize.sum.editp <- ve303.colonize.sum.alt %>% 
  select(MGS.Sample.ID, Subject.Number, Visit.Name, Visit.Name.Norm,
         Day.in.treatment, prop_det_alt_sum, TRT.norm, rec.diagnosis.Wk8, rec.day.in.treatment) %>% 
  spread(Visit.Name, prop_det_alt_sum) %>% 
  select(MGS.Sample.ID, Subject.Number,  Visit.Name.Norm,TRT.norm, rec.diagnosis.Wk8, `Day 7`, `Day 14` ) %>% 
  rename("Day.14" = `Day 14` ,   "Day.7" = `Day 7` )


#s1 <- split(ve303.colonize.sum.editp,  ve303.colonize.sum.editp$Subject.Number)

ve303.colonize.sum.editp1 <- do.call(rbind , (lapply(split(ve303.colonize.sum.editp,  ve303.colonize.sum.editp$Subject.Number),
                                                     calc_prop_7or14)) )

ve303.colonize.sum.alt <- left_join(ve303.colonize.sum.alt, ve303.colonize.sum.editp1 %>% 
                                      select(MGS.Sample.ID, prop_D7orD14), by = "MGS.Sample.ID" )

#######################################################################

### VE303 Abundance - only detected strains

##To plot zeroes
ve303.dnameta.det.alt = ve303.dnameta %>%
  mutate(est_relative_abundance_panel_alt = if_else(detection_status != "Detected", 0, est_relative_abundance_panel)) %>%
  group_by(Subject.Number, Visit.Name ) %>%
  mutate(.,est_relative_abundance_panel_alt_sum = sum(est_relative_abundance_panel_alt)) %>%
  ungroup()


# Filter the data frame to only include VE303 Strain Detected by the marker panel
ve303.dnameta.det <- ve303.dnameta %>% 
  filter(detection_status == "Detected")

# Filter the data frame to include VE303 Strain Detected AND Insufficient data by the marker panel
ve303.dnameta.det_ID <- ve303.dnameta %>% 
  filter(detection_status %in% c("Detected", "Insufficient data"))


ve303.dnameta.det.tp <- ve303.dnameta.det %>% 
  arrange(Subject.Number, Visit.Name.Norm) %>% 
  group_by(Visit.Name.Norm, organism) %>% # Select a single representative sample from each Subject for each Visit.Name.Norm & result in the same df as above
  distinct(Subject.Number, .keep_all = TRUE)
 
#######################################################################

source("VE303.abundance.summary.wrangling.R")

#######################################################################


#sigtaxa <- c("VE303-02", "VE303-03", "VE303-05" ,"VE303-08")
#sigtaxa <- c("VE303-02", "VE303-03","VE303-08")
sigtaxa <- c("VE303-03",  "VE303-08")
#sigtaxa <- c("VE303-02", "VE303-05" ,"VE303-08")
#sigtaxa <- c("VE303-03", "VE303-05" ,"VE303-08")
#sigtaxa <- c("VE303-02",  "VE303-08")
#sigtaxa <- c("VE303-02","VE303-04", "VE303-05" ,"VE303-08")
#sigtaxa <- c("VE303-05", "VE303-08")


ve303.abund.subj.total.time.det.alt_sig4 <- ve303.dnameta %>% 
  mutate(est_relative_abundance_panel_alt = if_else(detection_status != "Detected", 0, est_relative_abundance_panel)) %>%
  filter(organism %in% sigtaxa) %>%
  group_by(Subject.Number, Visit.Name.Norm, organism) %>%
  mutate(.,est_relative_abundance_panel_alt_sum = sum(est_relative_abundance_panel_alt)) %>%
  ungroup()%>%
  group_by(MGS.Sample.ID, Visit.Name.Norm) %>% 
  summarise_at(vars(est_relative_abundance_panel_alt_sum), funs(sum)) %>%
  ungroup() %>% 
  select(-Visit.Name.Norm) %>% 
  left_join(., dnameta, by = "MGS.Sample.ID") %>% 
  left_join(., cohort.count, by = "TRT.norm")

#######################################################################


KM_cols <- c("TRT.norm",  "Subject.Number", "Visit.Name.Norm","Visit.Name", "Day.in.treatment",
             "est_relative_abundance_panel_alt_sum","MGS.Sample.ID", "collect.pre.post.rec.binary",
             "rec.diagnosis", "rec.diagnosis.Wk8", "collect.pre.post.abx", "rec.day.in.treatment",
             "rec.day.in.treatment.surv" )

change_name <- function(df_plot){df_plot$est_relative_abundance_panel_alt_sum <- df_plot$est_relative_abundance_panel_alt
return(df_plot)}

df_plot_KM_indiv <- split(ve303.dnameta.det.alt, ve303.dnameta.det.alt$organism)
#df_plot_KM_indiv <- split(ve303.dnameta.det.alt[KM_cols], ve303.dnameta.det.alt[KM_cols]$organism)
df_plot_KM_indiv1 <- lapply(df_plot_KM_indiv, change_name)
ve303.colonize.sum.alt1 <- ve303.colonize.sum.alt %>% rename(est_relative_abundance_panel_alt_sum = prop_det_alt_sum)
ve303.colonize.sum.alt2 <- ve303.colonize.sum.alt %>% rename(est_relative_abundance_panel_alt_sum = prop_D7orD14)

a1 <- list("Sig4 Sum" = ve303.abund.subj.total.time.det.alt_sig4[KM_cols])
a2 <- list( "Total VE303" = ve303.abund.subj.total.time.det.alt[KM_cols])

if(prop_label == "Separate Days"){a2b <- list( "Proportion VE303" = ve303.colonize.sum.alt1[KM_cols])}
if(prop_label == "D7orD14 proportion"){a2b <- list( "Proportion VE303" = ve303.colonize.sum.alt2[KM_cols])}

a3 <- df_plot_KM_indiv1


combined <- c(a1, a2, a2b, a3)

df_plot_KM_total <- combined[[ plot_tax ]][KM_cols]


#if(plot_tax == "Sig4 Sum") { df_plot_KM_total <- ve303.abund.subj.total.time.det.alt_sig4[KM_cols] } 
#if(plot_tax == "Total VE303") { df_plot_KM_total <- ve303.abund.subj.total.time.det.alt[KM_cols] } 


#######################################################################
df_plot_KM_total$rec.binary <- if_else(df_plot_KM_total$rec.diagnosis.Wk8 == "Bad", 
                                       1, 0)

df_plot_KM_total$trt.binary <- if_else(df_plot_KM_total$TRT.norm == "Placebo", 0,
                                       if_else(df_plot_KM_total$TRT.norm == "VE303 Low Dose", 1, 2))

#######################################################################

N_Plac <- length(unique(df_plot_KM_total[(df_plot_KM_total$TRT.norm == "Placebo"),]$Subject.Number))
Sub_Plac <- (unique(df_plot_KM_total[(df_plot_KM_total$TRT.norm == "Placebo"),]$Subject.Number))

N_HD <- length(unique(df_plot_KM_total[(df_plot_KM_total$TRT.norm == "VE303 High Dose"),]$Subject.Number))
Sub_HD <-  (unique(df_plot_KM_total[(df_plot_KM_total$TRT.norm == "VE303 High Dose"),]$Subject.Number))

N_LD <- length(unique(df_plot_KM_total[(df_plot_KM_total$TRT.norm == "VE303 Low Dose"),]$Subject.Number))
Sub_LD <-  (unique(df_plot_KM_total[(df_plot_KM_total$TRT.norm == "VE303 Low Dose"),]$Subject.Number))

Subs_All <- (unique(df_plot_KM_total$Subject.Number))

#######################################################################

if(plot_coh == "Dosed") {df_plot_KM_total <- df_plot_KM_total %>% filter(TRT.norm != "Placebo")}
if(plot_coh == "Low Dosed") {df_plot_KM_total <- df_plot_KM_total %>% filter(TRT.norm == "VE303 Low Dose")}
if(plot_coh == "Dosed and Placebo") {df_plot_KM_total <- df_plot_KM_total }


if(prop_label != "D7orD14 proportion" | plot_tax != "Proportion VE303") {
  df_plot_KM_total_subj <- df_plot_KM_total %>% 
  filter(Day.in.treatment < day_engraft & Day.in.treatment > day_engraft_start) %>%
  group_by(Subject.Number) %>%
  mutate(.,est_relative_abundance_panel_alt_sum_Alltime = sum(est_relative_abundance_panel_alt_sum)) %>%
  ungroup() %>%
  distinct(Subject.Number, .keep_all = TRUE) %>%
  select(TRT.norm,Subject.Number,est_relative_abundance_panel_alt_sum_Alltime, 
         rec.diagnosis, rec.diagnosis.Wk8, rec.day.in.treatment, rec.day.in.treatment.surv,
         trt.binary, rec.binary) }

if(prop_label == "D7orD14 proportion" & plot_tax == "Proportion VE303"){
  df_plot_KM_total_subj <- df_plot_KM_total %>% 
    filter(!is.na(est_relative_abundance_panel_alt_sum)) %>%
    rename(est_relative_abundance_panel_alt_sum_Alltime = (est_relative_abundance_panel_alt_sum)) %>%
    distinct(Subject.Number, .keep_all = TRUE) %>%
    select(TRT.norm,Subject.Number,est_relative_abundance_panel_alt_sum_Alltime, 
           rec.diagnosis, rec.diagnosis.Wk8, rec.day.in.treatment, rec.day.in.treatment.surv,
           trt.binary, rec.binary)  }

df_plot_KM_total_subj <- df_plot_KM_total_subj[order(-df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime),]



#if(quantile == "Median"){mthresh <- as.numeric(quantile(df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime)[3])}

if(quantile == "Zero"){mthresh <- 0.000}
if(quantile == "One"){mthresh <- 0.125}
if(quantile == "Two"){mthresh <- 0.250}
if(quantile == "Three"){mthresh <- 0.375}
if(quantile == "Median"){
  mthresh <- as.numeric(median(df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime) )}
if(quantile == "Five"){mthresh <- 0.625}
if(quantile == "Third quantile"){
  mthresh <- as.numeric(quantile(df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime)[4])}
if(quantile == "Seven"){mthresh <- 0.875}


mthresh2 <- as.numeric(quantile(df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime)[2])
mthresh3 <- as.numeric(quantile(df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime)[3])
mthresh4 <- as.numeric(quantile(df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime)[4])


as.numeric(quantile(df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime))


#######################################################################

df_plot_KM_total_subj$Engraftment <- if_else(df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime <= mthresh, 
                                             "Low Engraft", "High Engraft")

df_plot_KM_total_subj$Engraftment_Quartiles <- if_else(df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime <= mthresh2, "Q1 Engraft",
                                                       if_else( ( (mthresh2 < df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime) 
                                                                  & (df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime <= mthresh3) ), "Q2 Engraft",
                                                                if_else(((mthresh3 < df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime) 
                                                                         & (df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime<= mthresh4) ), "Q3 Engraft", "Q4 Engraft") ) )

df_plot_KM_total$trt.binary <- if_else(df_plot_KM_total$TRT.norm == "Placebo", 0,
                                       if_else(df_plot_KM_total$TRT.norm == "VE303 Low Dose", 1, 2) )

                                            
df_plot_KM_total_subj$Engraftment.binary <- if_else(df_plot_KM_total_subj$Engraftment == "High Engraft", 
                                                    1, 0)

#######################################################################

if(plot_tax != "Proportion VE303"){
ggplot(df_plot_KM_total_subj, aes(x=Engraftment , y=est_relative_abundance_panel_alt_sum_Alltime, 
                                  fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8)) + 
  geom_boxplot(alpha=0.6,  outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge(), show.legend = FALSE) +
  theme_bw() + 
  facet_grid(~Engraftment, scales="free_x", space = "free") + 
  labs(title= "VE303 Colonization Summary", y="Engraftment", x="Treatment", fill="rec.diagnosis.Wk8") + 
  scale_fill_manual(values = rec.cols) + 
  scale_color_manual(values = rec.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1), 
        strip.text = element_text(size = 13)) + 
  scale_y_log10()
  #scale_y_continuous(limits = c(0,0.5) )
  ggsave(filename = here(results, paste(Sys.Date(),  taxplot  ,plot_coh, plot_tax, day_engraft_label, "Engraftment across recurrence.pdf", sep = " ")),
       width = 10, height = 5)
}

if(plot_tax == "Proportion VE303"){
  
  ggplot(df_plot_KM_total_subj, aes(x=Engraftment , y=est_relative_abundance_panel_alt_sum_Alltime, 
                                    fill=rec.diagnosis.Wk8, color=rec.diagnosis.Wk8)) + 
    geom_boxplot(alpha=0.6,  outlier.size = 0) +
    geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge(), show.legend = FALSE) +
    theme_bw() + 
    facet_grid(~Engraftment, scales="free_x", space = "free") + 
    labs(title= "VE303 Proportion Summary", y="Engraftment", x="Treatment", fill="rec.diagnosis.Wk8") + 
    scale_fill_manual(values = rec.cols) + 
    scale_color_manual(values = rec.cols) + 
    theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
          axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1), 
          strip.text = element_text(size = 13)) 
   # scale_y_log10()
  #scale_y_continuous(limits = c(0,0.5) )
  ggsave(filename = here(results, paste(Sys.Date(),  taxplot,  plot_coh, plot_tax, day_engraft_label ,prop_label,"Engraftment across recurrence.pdf", sep = " ")),
         width = 10, height = 5)
  
}

#######################################################################

#coh1.cols <- c("Placebo" = "#66c2a5", "VE303 Low Dose" = "#fc8d62", "VE303 High Dose" = "#8da0cb")
coh1.cols <- c("Placebo" = "#1b9e77", "VE303 Low Dose" = "#7570b3", "VE303 High Dose" = "#d95f02")


xps <- cbind(rankord = as.numeric(rownames(df_plot_KM_total_subj)), (df_plot_KM_total_subj)) 
Nsamp <- length(xps$rankord)

if(plot_tax != "Proportion VE303"){
ggplot(xps, aes(x = rankord, y =  est_relative_abundance_panel_alt_sum_Alltime,fill=TRT.norm,
                color=TRT.norm )) + geom_bar(stat="identity",  width=0.6)  + ylim(0.0,0.2) +
  
  geom_vline(xintercept =  (1.5 +length(xps$est_relative_abundance_panel_alt_sum_Alltime))- whichmedian(sort(xps$est_relative_abundance_panel_alt_sum_Alltime)) , 
               linetype="longdash", color="black") +
  theme_bw() +
  scale_fill_manual(values = coh1.cols) + 
  scale_color_manual(values = coh1.cols)  + 
  #scale_y_log10() +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14)) +
  theme(axis.title = element_text(size = 14,hjust = 0.5, vjust = 0.8),
        axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14), legend.text=element_text(size=12) , legend.title = element_blank()) +
  
  labs(x="rank", y="Relative abundance")

ggsave(filename = here(results, paste(Sys.Date(),  taxplot,plot_coh, plot_tax, day_engraft_label, "Median BARS Engraftment distribution.pdf", sep = " ")), 
       width = 7, height = 4) }


if(plot_tax == "Proportion VE303"){
  
  ggplot(xps, aes(x = rankord, y = est_relative_abundance_panel_alt_sum_Alltime,fill=TRT.norm,
                  color=TRT.norm )) + geom_bar(stat="identity",  width=0.6)  +  

    geom_vline(xintercept =  (1.5 +length(xps$est_relative_abundance_panel_alt_sum_Alltime))- whichmedian(sort(xps$est_relative_abundance_panel_alt_sum_Alltime)) , 
               linetype="longdash", color="black") +
    
    theme_bw() +
    scale_fill_manual(values = coh1.cols) + 
    scale_color_manual(values = coh1.cols)  + 
    #scale_y_log10() +
    theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
          axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14)) +
    theme(axis.title = element_text(size = 14,hjust = 0.5, vjust = 0.8),
          axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
          strip.text.y = element_text(size=14), legend.text=element_text(size=12) ,
          legend.title = element_blank()) +
    
    labs(x="rank", y="VE303 Proportion")
  
  ggsave(filename = here(results, paste(Sys.Date(),  taxplot,plot_coh, plot_tax, day_engraft_label, prop_label,"Median BARS Engraftment distribution.pdf", sep = " ")), 
         width = 7, height = 4)
  
}

################################################################################################
################################################################################################

xps <- cbind(rankord = as.numeric(rownames(df_plot_KM_total_subj)), (df_plot_KM_total_subj)) 
Nsamp <- length(xps$rankord)

if(plot_tax != "Proportion VE303"){
cohdis <- ggplot(xps, aes(x = rankord, y = est_relative_abundance_panel_alt_sum_Alltime,fill=TRT.norm,
                color=TRT.norm )) + geom_point(size=2.2)  + ylim(0.0,0.4) +
    geom_vline(xintercept =  (1.5 +length(xps$est_relative_abundance_panel_alt_sum_Alltime))- whichmedian(sort(xps$est_relative_abundance_panel_alt_sum_Alltime)) , 
               linetype="longdash", color="black") +

    theme_bw() +
  scale_fill_manual(values = coh1.cols) + 
  scale_color_manual(values = coh1.cols)  + 
  scale_y_log10() +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14)) +
  theme(axis.title = element_text(size = 14,hjust = 0.5, vjust = 0.8),
        axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14), legend.text=element_text(size=12) , legend.title = element_blank()) +
  
  labs(x="rank", y="log(RA) + 1e-4")

ggsave(filename = here(results, paste(Sys.Date(),  taxplot,plot_coh, plot_tax, day_engraft_label,  "Median Engraftment curve distribution CohCol.pdf", sep = " ")), 
       width = 7, height = 4)

write.csv(xps %>% rename(strain_abun=est_relative_abundance_panel_alt_sum_Alltime) %>%
            select(rankord, strain_abun, Subject.Number, TRT.norm),
          file = here(results, paste(Sys.Date(), taxplot,plot_coh, plot_tax, day_engraft_label, "dist.csv", sep = " ")), row.names = FALSE)

}


if(plot_tax == "Proportion VE303"){
  
  cohdis <-  ggplot(xps, aes(x = rankord, y =  8*est_relative_abundance_panel_alt_sum_Alltime,fill=TRT.norm,
                  color=TRT.norm )) + geom_point(size=2.2)  +  
    geom_vline(xintercept =  (1.5 +length(xps$est_relative_abundance_panel_alt_sum_Alltime))- whichmedian(sort(xps$est_relative_abundance_panel_alt_sum_Alltime)) , 
               linetype="longdash", color="black") +  
    theme_bw() +
    scale_fill_manual(values = coh1.cols) + 
    scale_color_manual(values = coh1.cols)  + 
    theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
          axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14)) +
    theme(axis.title = element_text(size = 14,hjust = 0.5, vjust = 0.8),
          axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
          strip.text.y = element_text(size=14), legend.text=element_text(size=12) , legend.title = element_blank()) +
    
    labs(x="rank", y="VE303 Strains (N)")
    #labs(x="rank", y="VE303 Proportion")
  
  ggsave(filename = here(results, paste(Sys.Date(),  taxplot,plot_coh, plot_tax, day_engraft_label,  "Median Engraftment curve distribution CohCol.pdf", sep = " ")), 
         width = 7, height = 4)
  
}

################################################################################################
if(plot_tax != "Proportion VE303"){
  recdis <- ggplot(xps, aes(x = rankord, y =  est_relative_abundance_panel_alt_sum_Alltime,fill=rec.diagnosis.Wk8,
                color=rec.diagnosis.Wk8 )) + geom_point(size=2.2)  + ylim(0.0,0.4) +
    geom_vline(xintercept =  (1.5 +length(xps$est_relative_abundance_panel_alt_sum_Alltime))- whichmedian(sort(xps$est_relative_abundance_panel_alt_sum_Alltime)) , 
               linetype="longdash", color="black") +
  theme_bw() +
     
  scale_fill_manual(values = rec.cols,breaks=c("Bad", "Good"),labels=c("Recurrent","Non Recurrent")) + 
  scale_color_manual(values = rec.cols,breaks=c("Bad", "Good"),labels=c("Recurrent","Non Recurrent"))  + 
  scale_y_log10() +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14)) +
  theme(axis.title = element_text(size = 14,hjust = 0.5, vjust = 0.8),
        axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14), legend.text=element_text(size=12) , legend.title = element_blank()) +
  
  labs(x="rank", y="log(RA) + 1e-4")

ggsave(filename = here(results, paste(Sys.Date(),  taxplot,plot_coh, plot_tax, day_engraft_label,  "Median Engraftment curve distribution RecCol.pdf", sep = " ")), 
       width = 7, height = 4)}



if(plot_tax == "Proportion VE303"){
  
  recdis <-  ggplot(xps, aes(x = rankord, y =  8*est_relative_abundance_panel_alt_sum_Alltime,fill=rec.diagnosis.Wk8,
                  color=rec.diagnosis.Wk8 )) + geom_point(size=2.2)   +
    geom_vline(xintercept =  (1.5+length(xps$est_relative_abundance_panel_alt_sum_Alltime))- whichmedian(sort(xps$est_relative_abundance_panel_alt_sum_Alltime)) , 
               linetype="longdash", color="black") +
    
    theme_bw() +
    scale_fill_manual(values = rec.cols,breaks=c("Bad", "Good"),labels=c("Recurrent","Non Recurrent")) + 
    scale_color_manual(values = rec.cols,breaks=c("Bad", "Good"),labels=c("Recurrent","Non Recurrent"))  + 
    
    theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
          axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14)) +
    
    theme(axis.title = element_text(size = 14,hjust = 0.5, vjust = 0.8),
          axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
          strip.text.y = element_text(size=14), 
          legend.text=element_text(size=12) , legend.title = element_blank()) +
    labs(x="rank", y="VE303 Strains (N)")
    #labs(x="rank", y="VE303 Proportion")

  
  ggsave(filename = here(results, 
                         paste(Sys.Date(),  taxplot, plot_coh, plot_tax, day_engraft_label, prop_label ,
                                        "Median Eng dist RecCol.pdf", sep = " ")), 
         width = 7, height = 4)
  
}


################################################################################################
patchd_1 <- cohdis / recdis
ggsave(filename = here(results, paste(Sys.Date(),  taxplot,plot_coh, 
                                      plot_tax, day_engraft_label, prop_label , 
                                      "Coh_Rec_Engraft Distributions.pdf", sep = " ")), 
       width = 6, height = 7, dpi = 400)

################################################################################################


ggplot(xps, aes(x = rankord, y = 0.0001+est_relative_abundance_panel_alt_sum_Alltime,fill=TRT.norm,
                color=TRT.norm )) + geom_point(size=2.2)  + ylim(0.0,0.4) +
  
  geom_vline(xintercept = Nsamp*(1/4) , linetype="longdash", color="black") +
  theme_bw() +
  scale_fill_manual(values = coh1.cols) + 
  scale_color_manual(values = coh1.cols)  + 
  scale_y_log10() +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14)) +
  theme(axis.title = element_text(size = 14,hjust = 0.5, vjust = 0.8),
        axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14), legend.text=element_text(size=12) , legend.title = element_blank()) +
  
  labs(x="rank", y="log(RA) + 1e-4")

ggsave(filename = here(results, paste(Sys.Date(),  taxplot,plot_coh, plot_tax, day_engraft_label,  "Third Quantile Engraftment curve distribution.pdf", sep = " ")), 
       width = 7, height = 4)

################################################################################################



ggplot(xps, aes(x = rankord, y = 0.0001+est_relative_abundance_panel_alt_sum_Alltime,fill=TRT.norm,
                color=TRT.norm )) + geom_bar(stat="identity",  width=0.6)  + ylim(0.0,0.2) +
  
  geom_vline(xintercept = Nsamp*(1/4) , linetype="longdash", color="black") +
  theme_bw() +
  scale_fill_manual(values = coh1.cols) + 
  scale_color_manual(values = coh1.cols)  + 
  #scale_y_log10() +
  theme(axis.text.y = element_text(size=14), strip.text = element_text(size=14),
        axis.text.x = element_text(angle = 0, hjust = 1, vjust = 0.5, size = 14)) +
  theme(axis.title = element_text(size = 14,hjust = 0.5, vjust = 0.8),
        axis.text.x = element_text(size = 14), strip.text.x = element_text(size=14),
        strip.text.y = element_text(size=14), legend.text=element_text(size=12) , legend.title = element_blank()) +
  
  labs(x="rank", y="Relative abundance")

ggsave(filename = here(results, paste(Sys.Date(),  taxplot,plot_coh, plot_tax, day_engraft_label, "Third quantile BARS Engraftment distribution.pdf", sep = " ")), 
       width = 7, height = 4)


################################################################################################



df_plot_KM_total_Eng <- df_plot_KM_total %>% 
  left_join(., (df_plot_KM_total_subj[c("Subject.Number","Engraftment")]), by = "Subject.Number")  %>%
  filter(!is.na(Engraftment))

#######################################################################

meta <- "Engraftment"

inputs_trts <- c("VE303 Low Dose", "VE303 High Dose", "Placebo")
inputs_eng <- c("High Engraft", "Low Engraft")

if(meta == "Engraftment"){
  inputs <- inputs_eng}
if(meta == "Treatment"){
  inputs <- inputs_trts}

df_plot <- df_plot_KM_total_Eng

#######################################################################

#Write code using the "survival" package to generate KM curves and compute stats based on log rank test for time-to-event data


df_plot_KM_total_subj$time_rec <- df_plot_KM_total_subj$rec.day.in.treatment.surv
df_plot_KM_total_subj$status_rec <- df_plot_KM_total_subj$rec.binary
df_plot_KM_total$time_rec <- df_plot_KM_total$rec.day.in.treatment.surv
df_plot_KM_total$status_rec <- df_plot_KM_total$rec.binary


array_survival <- c("Subject.Number","trt.binary","TRT.norm","Engraftment.binary", 
                    "Engraftment_Quartiles", "Engraftment", "status_rec", "time_rec")

engrafters <-df_plot_KM_total_subj[array_survival]
 


#######################################################################



run_survival_models <- function(event_type, engrafters){
  if(event_type == "Recurrence"){ engrafters$time <- engrafters$time_rec
  engrafters$status <- engrafters$status_rec  }
  
  km <- with(engrafters, Surv(time, status))
  km_fit <- survfit(Surv(time, status) ~ 1, data=engrafters)
  km_eng_fit <- survfit(Surv(time, status) ~ Engraftment , data=engrafters)
  km_eng_quartile_fit <- survfit(Surv(time, status) ~ Engraftment_Quartiles , data=engrafters)
  km_trt_fit <- survfit(Surv(time, status) ~ TRT.norm , data=engrafters)
  cox_eng <- coxph(Surv(time, status) ~   Engraftment , data = engrafters)
  cox_eng_quartile <- coxph(Surv(time, status) ~   Engraftment_Quartiles , data = engrafters)
  cox_trt <- coxph(Surv(time, status) ~   TRT.norm , data = engrafters)
  cox_full <- coxph(Surv(time, status) ~   Engraftment + TRT.norm , data = engrafters)
  cox_full_quartile <- coxph(Surv(time, status) ~   Engraftment_Quartiles + TRT.norm , data = engrafters)
  
  Models <- list(km, km_fit, km_eng_fit, km_trt_fit, cox_eng, cox_trt, cox_full, km_eng_quartile_fit, cox_eng_quartile, cox_full_quartile)
  return(Models)}

surv_results <- run_survival_models(event_type = event_type, engrafters)

if(event_type == "Recurrence"){AE_subset <- ""}
 

#labels1 <- paste0("Below ", quantile, " colonization")
#labels2 <- paste0("Above ", quantile, " colonization")

labels1 <- paste0("0-", mthresh*8, " strains colonize")
labels2 <- paste0(1+mthresh*8,"-8", " strains colonize")

#if(plot_tax == "Proportion VE303" & plot_coh == "Dosed" & prop_label != "D7orD14 proportion"){
#  labels1 <- paste0("0-4 strains colonize")
#  labels2 <- paste0("5-8 strains colonize")}

#######################################################################

#print("THIS IS THE PLOT BERNAT HAS ASKED FOR - IN CASE YOU NEED TO MODIFY IT. It starts LINE 490 of the KM Exposure recurrence quartile loops.R")
day_end <- 100
ttt <- FALSE
linesize <- 1.5

km_eng_fit <- surv_results[[3]]

eng.cols <- c("Low Engraft" = "indianred2", "High Engraft" = "cyan3")
eng.cols1 <- c("Low Colonize" = "indianred2", "High Colonize" = "cyan3")

panel3_clin <- autoplot(km_eng_fit ,surv.size = linesize, conf.int = ttt , fun = "event") + 
#panel3_clin <- autoplot(km_eng_fit ,surv.size = linesize, conf.int = ttt) + 

  coord_cartesian(xlim=c(0, day_end)) +
  theme_bw() +
  theme_classic() +
  geom_vline(xintercept = 15 , linetype="longdash", color="black") +
  
  # ylim(0.0,1.05) +
  # annotate("rect", xmin = 0, xmax = 14, ymin = 0, ymax = 1,
  #          alpha = .1,fill = "blue") +
  # 
  # annotate('text', x = 7, y = 0.52, label = 'Dosing') +
  # annotate('text', x = 7, y = 0.48, label = 'Period') +
  
  
  ylim(0.0,0.5) +
  annotate("rect", xmin = 0, xmax = 14, ymin = 0, ymax = 0.5,
           alpha = .1,fill = "blue") +
  annotate('text', x = 7, y = 0.28, label = 'Dosing') +
  annotate('text', x = 7, y = 0.257, label = 'Period') +
  
  
  scale_colour_manual("",values=c("indianred2","cyan3"),breaks=c("Low Engraft", "High Engraft"), 
                      labels=c(labels1,labels2))+
  scale_fill_discrete(name="Engraftment",labels=c(labels1, labels2)) +
  
  theme(axis.text = element_text(size = 12,hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 12,hjust = 0.5, vjust = 0.8),
        legend.text=element_text(size=12), legend.title = element_blank(), 
        #legend.position = c(0.75, 0.43)) + 
        legend.position = c(0.75, 0.83)) + 
  
  #labs(  title= "Early Colonization", y="Recurrence-Free Probability", x="Time (Days)" ) 
  labs(  title= "Early Colonization", y="Probability of Recurrence", x="Time (Days)" )

ggsave(filename = here(results, paste(Sys.Date(),  taxplot, plot_coh, plot_tax, event_type, AE_subset, quantile 
                                      , day_engraft_label, prop_label,"KM EV_Eng.pdf", 
                                      sep = " ")), width = 8, height = 5, dpi = 800)




patchd_1 <- (cohdis / recdis | panel3_clin ) + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(1,1.4))
ggsave(filename = here(results, paste(Sys.Date(), taxplot, plot_coh, plot_tax, event_type, 
                                      AE_subset, quantile 
                                      , day_engraft_label, prop_label,
                                      "COMB CohRec_KM EV.pdf", 
                                      sep = " ")), width = 14, height = 6, dpi = 400)

if(plot_coh == "Dosed" & 
   plot_tax == "Proportion VE303" &
   event_type == "Recurrence" &
   quantile == "Median" &
   #day_engraft_label == "D14 only" &
   prop_label == "Separate Days"){
  #labels1 <- paste0("0-4 strains colonize")
  #labels2 <- paste0("5-8 strains colonize")
    
  panel3_clin_median<- autoplot(km_eng_fit ,surv.size = linesize, conf.int = ttt, fun = "event") + 
  #panel3_clin_median<- autoplot(km_eng_fit ,surv.size = linesize, conf.int = ttt) + 
      coord_cartesian(xlim=c(0, day_end)) +
    theme_bw() +
    theme_classic() +
    geom_vline(xintercept = 15 , linetype="longdash", color="black") +
    
   # ylim(0.0,1.0) +
   # annotate("rect", xmin = 0, xmax = 14, ymin = 0, ymax = 1.0,
   #          alpha = .1,fill = "blue") +
   # annotate('text', x = 7, y = 0.53, label = 'Dosing') +
   # annotate('text', x = 7, y = 0.48, label = 'Period') +
   #  
    ylim(0.0,0.5) +
    annotate("rect", xmin = 0, xmax = 14, ymin = 0, ymax = 0.5,
             alpha = .1,fill = "blue") +
    annotate('text', x = 7, y = 0.28, label = 'Dosing') +
    annotate('text', x = 7, y = 0.257, label = 'Period') +

     
    scale_colour_manual("",values=c("indianred2","cyan3"),breaks=c("Low Engraft", "High Engraft"), 
                        labels=c(labels1,labels2))+
    scale_fill_discrete(name="Engraftment",labels=c(labels1, labels2)) +
    
    theme(axis.text = element_text(size = 12,hjust = 1, vjust = 0.5),
          axis.title = element_text(size = 12,hjust = 0.5, vjust = 0.8),
          legend.text=element_text(size=12), legend.title = element_blank(), 
          #legend.position = c(0.75, 0.53)) + 
          legend.position = c(0.75, 0.83)) + 
  
    labs(  title= "Early Colonization", y="Probability of Recurrence", x="Time (Days)" )
    #labs(  title= "Early Colonization", y="Recurrence-Free Probability", x="Time (Days)" )
 
  
  ggsave(filename = here(results, paste(Sys.Date(),  taxplot, plot_coh, plot_tax, event_type, 
                                        AE_subset, quantile 
                                        , day_engraft_label, prop_label,"KM EV_Eng.pdf", 
                                        sep = " ")), width = 8, height = 5, dpi = 800) 
  
  patchd_1 <- (cohdis / recdis | panel3_clin_median ) + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(1,1.4))
  ggsave(filename = here(results, paste(Sys.Date(), taxplot, plot_coh, plot_tax, event_type, 
                                        AE_subset, quantile 
                                        , day_engraft_label, prop_label,
                                        "COMB CohRec_KM EV.pdf", 
                                        sep = " ")), width = 14, height = 6, dpi = 400)
  
  
  panel3_clin_median_eventfree <- autoplot(km_eng_fit ,surv.size = linesize, conf.int = ttt ) + 
    coord_cartesian(xlim=c(0, day_end)) +
    theme_bw() +
    theme_classic() +
    geom_vline(xintercept = 15 , linetype="longdash", color="black") +
    
    ylim(0.0,1.0) +
    annotate("rect", xmin = 0, xmax = 14, ymin = 0, ymax = 1.0,
             alpha = .1,fill = "blue") +
    annotate('text', x = 7, y = 0.53, label = 'Dosing') +
    annotate('text', x = 7, y = 0.48, label = 'Period') +
    
    
    scale_colour_manual("",values=c("indianred2","cyan3"),breaks=c("Low Engraft", "High Engraft"), 
                        labels=c(labels1,labels2))+
    scale_fill_discrete(name="Engraftment",labels=c(labels1, labels2)) +
    
    theme(axis.text = element_text(size = 12,hjust = 1, vjust = 0.5),
          axis.title = element_text(size = 12,hjust = 0.5, vjust = 0.8),
          legend.text=element_text(size=12), legend.title = element_blank(), 
          legend.position = c(0.75, 0.53)) + 
    #legend.position = c(0.75, 0.83)) + 
    
    #  labs(  title= "Early Colonization", y="Probability of Recurrence", x="Time (Days)" )
    labs(  title= "Early Colonization", y="Recurrence-Free Probability", x="Time (Days)" )
  
  
  ggsave(filename = here(results, paste(Sys.Date(),  taxplot, plot_coh, plot_tax, event_type, 
                                        AE_subset, quantile 
                                        , day_engraft_label, prop_label,"KM EVFREE_Eng.pdf", 
                                        sep = " ")), width = 8, height = 5, dpi = 800) 
  
  patchd_1 <- (cohdis / recdis | panel3_clin_median_eventfree ) + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(1,1.4))
  ggsave(filename = here(results, paste(Sys.Date(), taxplot, plot_coh, plot_tax, event_type, 
                                        AE_subset, quantile 
                                        , day_engraft_label, prop_label,
                                        "COMB CohRec_KM EVFREE.pdf", 
                                        sep = " ")), width = 14, height = 6, dpi = 400)

  }


if(plot_coh == "Dosed" & 
   plot_tax == "VE303-08" &
   event_type == "Recurrence" &
   quantile == "Median" &
   #day_engraft_label == "D14 only" &
   prop_label == "Separate Days"){
  labels1 <- paste0("Colonization = 0")
  labels2 <- paste0("Colonization > 0")
  
  panel3_clin_median_ve303_08 <- autoplot(km_eng_fit ,surv.size = linesize, conf.int = ttt, fun = "event") + 
    #panel3_clin_median<- autoplot(km_eng_fit ,surv.size = linesize, conf.int = ttt) + 
    coord_cartesian(xlim=c(0, day_end)) +
    theme_bw() +
    theme_classic() +
    geom_vline(xintercept = 15 , linetype="longdash", color="black") +
    
    # ylim(0.0,1.0) +
    # annotate("rect", xmin = 0, xmax = 14, ymin = 0, ymax = 1.0,
    #          alpha = .1,fill = "blue") +
    # annotate('text', x = 7, y = 0.53, label = 'Dosing') +
    # annotate('text', x = 7, y = 0.48, label = 'Period') +
    #  
    ylim(0.0,0.5) +
    annotate("rect", xmin = 0, xmax = 14, ymin = 0, ymax = 0.5,
             alpha = .1,fill = "blue") +
    annotate('text', x = 7, y = 0.28, label = 'Dosing') +
    annotate('text', x = 7, y = 0.257, label = 'Period') +
    
    
    scale_colour_manual("",values=c("indianred2","cyan3"),breaks=c("Low Engraft", "High Engraft"), 
                        labels=c(labels1,labels2))+
    scale_fill_discrete(name="Engraftment",labels=c(labels1, labels2)) +
    
    theme(axis.text = element_text(size = 12,hjust = 1, vjust = 0.5),
          axis.title = element_text(size = 12,hjust = 0.5, vjust = 0.8),
          legend.text=element_text(size=12), legend.title = element_blank(), 
          #legend.position = c(0.75, 0.53)) + 
          legend.position = c(0.75, 0.83)) + 
    
    labs(  title= "Early Colonization", y="Probability of Recurrence", x="Time (Days)" )
  #labs(  title= "Early Colonization", y="Recurrence-Free Probability", x="Time (Days)" )
  
  
  ggsave(filename = here(results, paste(Sys.Date(),  taxplot, plot_coh, plot_tax, event_type, 
                                        AE_subset, quantile 
                                        , day_engraft_label, prop_label,"KM EV_Eng.pdf", 
                                        sep = " ")), width = 8, height = 5, dpi = 800) 
  
  patchd_1 <- (cohdis / recdis | panel3_clin_median ) + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(1,1.4))
  ggsave(filename = here(results, paste(Sys.Date(), taxplot, plot_coh, plot_tax, event_type, 
                                        AE_subset, quantile 
                                        , day_engraft_label, prop_label,
                                        "COMB CohRec_KM EV.pdf", 
                                        sep = " ")), width = 14, height = 6, dpi = 400)

  panel3_clin_median_ve303_08_eventfree <- 
    autoplot(km_eng_fit ,surv.size = linesize, conf.int = ttt ) + 
    coord_cartesian(xlim=c(0, day_end)) +
    theme_bw() +
    theme_classic() +
    geom_vline(xintercept = 15 , linetype="longdash", color="black") +
    
    ylim(0.0,1.0) +
    annotate("rect", xmin = 0, xmax = 14, ymin = 0, ymax = 1.0,
            alpha = .1,fill = "blue") +
    annotate('text', x = 7, y = 0.53, label = 'Dosing') +
    annotate('text', x = 7, y = 0.48, label = 'Period') +
      
    
    scale_colour_manual("",values=c("indianred2","cyan3"),breaks=c("Low Engraft", "High Engraft"), 
                        labels=c(labels1,labels2))+
    scale_fill_discrete(name="Engraftment",labels=c(labels1, labels2)) +
    
    theme(axis.text = element_text(size = 12,hjust = 1, vjust = 0.5),
          axis.title = element_text(size = 12,hjust = 0.5, vjust = 0.8),
          legend.text=element_text(size=12), legend.title = element_blank(), 
          legend.position = c(0.75, 0.53)) + 
          #legend.position = c(0.75, 0.83)) + 
    
  #  labs(  title= "Early Colonization", y="Probability of Recurrence", x="Time (Days)" )
  labs(  title= "Early Colonization", y="Recurrence-Free Probability", x="Time (Days)" )
  
  
  ggsave(filename = here(results, paste(Sys.Date(),  taxplot, plot_coh, plot_tax, event_type, 
                                        AE_subset, quantile 
                                        , day_engraft_label, prop_label,"KM EVFREE_Eng.pdf", 
                                        sep = " ")), width = 8, height = 5, dpi = 800) 
  
  patchd_1 <- (cohdis / recdis | panel3_clin_median_eventfree ) + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(1,1.4))
  ggsave(filename = here(results, paste(Sys.Date(), taxplot, plot_coh, plot_tax, event_type, 
                                        AE_subset, quantile 
                                        , day_engraft_label, prop_label,
                                        "COMB CohRec_KM EVFREE.pdf", 
                                        sep = " ")), width = 14, height = 6, dpi = 400)  
}
  

#######################################################################

#km_trt_fit <- surv_results[[4]]
df_plot_KM_total_trt <- df_plot_KM_total %>% 
  select(TRT.norm, Subject.Number, rec.diagnosis.Wk8,
         rec.day.in.treatment.surv, rec.binary, time_rec, status_rec) %>%
  distinct(Subject.Number, .keep_all = TRUE)

km_trt_fit <- survfit(Surv(time_rec, status_rec) ~ TRT.norm , data=adtte_mod)
p <- autoplot(km_trt_fit ,surv.size = linesize + 0.5, conf.int = ttt, fun = "event" ,
              risk.table = TRUE) + 
      labs(x="Time (Days)", y="Recurrence free probability") +
  coord_cartesian(xlim=c(0, 75)) +
  theme_bw() +
  theme_classic() +
#  geom_vline(xintercept = 15 , linetype="longdash", color="black") +
  
#  ylim(0.0,1) +
  scale_y_continuous(limits = c(0,0.5),breaks=c(0,0.1,0.2, 0.3 ,0.4,0.5)) +
  scale_x_continuous(   breaks = c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75),
                        labels = c("0","","","","","25","","","","","50","","","","","75")) +
  #scale_x_continuous(   breaks=c(0,25,50,75)) +
  #  breaks=c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75),

    #  ylim(0.0,0.5) +
  annotate("rect", xmin = 0, xmax = 14, ymin = 0, ymax = 0.5,
           alpha = .1,fill = "blue") +
  
  annotate('text', x = 7, y = 0.28, label = '14-Day') +
  annotate('text', x = 7, y = 0.26, label = 'Dosing') +
  annotate('text', x = 7, y = 0.24, label = 'Period') +
  scale_color_manual(values = coh.cols) +  
  scale_fill_manual(values = coh.cols ) + 
  
  theme(axis.text = element_text(size = 12,hjust = 0.5, vjust = 0.5),
        axis.title = element_text(size = 12,hjust = 0.5, vjust = 0.8),
        legend.text=element_text(size=12), legend.position = c(10, 0.53),
        legend.title = element_blank()) + 
  labs(  title= " ", y="Recurrence Probability", x="Days" ) 
ggsave(filename = here(results, paste(Sys.Date(),  taxplot, plot_coh, event_type , AE_subset, "KM_Treatment.pdf", sep = " ")), 
       width = 8, height = 6, dpi = 700)

#######################################################################


cox_eng<- surv_results[[5]]
cox_trt <- surv_results[[6]]
cox_full <- surv_results[[7]]
cox_eng_quartile <- surv_results[[9]]
cox_full_quartile <- surv_results[[10]]

#capture.output(summary(cox))
#print(" NEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXT ")

#capture.output(summary(cox_trt))
#print(" NEXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXT ")

#capture.output(summary(cox_full))



#######################################################################
km_eng_quartile_fit <- surv_results[[8]]

#eng.cols <- c("Low Engraft" = "indianred2", "High Engraft" = "cyan3")


p <- autoplot(km_eng_quartile_fit ,surv.size = linesize, conf.int = ttt) + labs(x="Time (Days)", y="Recurrence free probability") +
  coord_cartesian(xlim=c(0, day_end)) +
  theme_bw() +
  theme_classic() +
  geom_vline(xintercept = 15 , linetype="longdash", color="black") +
  
  ylim(0.0,1.05) +
  annotate("rect", xmin = 0, xmax = 14, ymin = 0, ymax = 1,
           alpha = .1,fill = "blue") +
  
  annotate('text', x = 7, y = 0.52, label = 'Dosing') +
  annotate('text', x = 7, y = 0.48, label = 'Period') +
  
  scale_colour_manual("",values=c("indianred2", "violet" , "seagreen", "cyan3"),breaks=c("Q1 Engraft", "Q2 Engraft", "Q3 Engraft", "Q4 Engraft"), 
                      labels=c("First quartile engraftment","Second quartile engraftment",
                               "Third quartile engraftment","Fourth quartile engraftment") )+
  scale_fill_discrete(name="Engraftment",labels=c("First quartile engraftment","Second quartile engraftment",
                                                  "Third quartile engraftment","Fourth quartile engraftment") ) +
  
  theme(axis.text = element_text(size = 12,hjust = 1, vjust = 0.5),
        axis.title = element_text(size = 12,hjust = 0.5, vjust = 0.8),
        legend.text=element_text(size=12), legend.title = element_blank(), legend.position = c(0.8, 0.25)) + 
 
  labs(  title= "Early Engraftment <= D14", y="Recurrence Free Probability", x="Time (Days)" ) 

ggsave(filename = here(results, paste(Sys.Date(),  taxplot, plot_coh, plot_tax, event_type, AE_subset, 
                                      day_engraft_label, prop_label, "KM_Quartile Engraftment.pdf", sep = " ")), width = 7, height = 6)
p

#######################################################################
#######################################################################


save_result <- list(KM_engraft = tidy(km_eng_fit), 
                    KM_engraft_quartile = tidy(km_eng_quartile_fit), 
                    KM_trt = tidy(km_trt_fit),
                    Cox_engraft = tidy(cox_eng),
                    Cox_engraft_quartile = tidy(cox_eng_quartile),
                    Cox_trt = tidy(cox_trt),
                    Cox_trt_engraft = tidy(cox_full),
                    Cox_trt_engraft_quartile = tidy(cox_full_quartile))



## Writes to multi sheet excel file ###
writexl::write_xlsx(save_result, here(results, paste(Sys.Date(),  taxplot, plot_coh, plot_tax ,event_type , AE_subset, quantile, 
                                                     day_engraft_label, prop_label, "KM Cox summary.xlsx", sep=" ")) )

## Capture eng Cox output including hazard ratios
capture.output(summary(cox_eng), file=here(results, paste(Sys.Date(),  taxplot, plot_coh,plot_tax , event_type , AE_subset, quantile, 
                                                          day_engraft_label, prop_label,"Cox_Eng summary.doc", sep=" ")) )

## Capture eng Cox output including hazard ratios
capture.output(summary(cox_eng_quartile), file=here(results, paste(Sys.Date(),  taxplot, plot_coh,plot_tax , event_type , AE_subset, 
                                                                   day_engraft_label, prop_label, "Cox_Eng_Quartile summary.doc", sep=" ")) )

## Capture TRT Cox output including hazard ratios
capture.output(summary(cox_trt), file=here(results, paste(Sys.Date(),  taxplot, plot_coh,plot_tax , event_type , AE_subset, quantile, 
                                                          day_engraft_label, prop_label,"Cox_Trt summary.doc", sep=" ")) )

## Capture full Cox output including hazard ratios
capture.output(summary(cox_full), file=here(results, paste(Sys.Date(),  taxplot, plot_coh, plot_tax ,event_type , AE_subset, quantile,
                                                           day_engraft_label, prop_label,"Cox_EngTrt summary.doc", sep=" ")) )

## Capture full Cox output including hazard ratios
capture.output(summary(cox_full_quartile), file=here(results, paste(Sys.Date(),  taxplot, plot_coh, plot_tax ,event_type , AE_subset,
                                                                    day_engraft_label, prop_label, "Cox_EngTrt_Quartile summary.doc", sep=" ")) )



#######################################################################
 
 
  
