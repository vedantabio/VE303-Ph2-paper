
# VE303 detection and calculate summary stats
# INCLUDES All DETECTED events Pre-Abx OR Pre- and Post-Abx  depending on selection in main script


# Change from character class to factor
ve303.dnameta_plot$detection_status <- factor(ve303.dnameta_plot$detection_status, levels = c("Detected", "Insufficient data", "Not detected"))

# Change from character class to factor
ve303.dnameta_plot.alt$detection_status <- factor(ve303.dnameta_plot.alt$detection_status, levels = c("Detected", "Insufficient data", "Not detected"))


#######################################################################

# To plot zeroes in colonization data
ve303.dnameta_plot = ve303.dnameta_plot %>%
  mutate(freq_alt = if_else(detection_status != "Detected", 0, freq)) %>%
  group_by(Subject.Number, Visit.Name.Norm, organism) %>%
  mutate(.,freq_sum = sum(freq_alt)) %>%
  ungroup()

ve303.dnameta_plot.alt = ve303.dnameta_plot.alt %>%
  mutate(freq_alt = if_else(detection_status != "Detected", 0, freq)) %>%
  group_by(Subject.Number, Visit.Name.Norm, organism) %>%
  mutate(.,freq_sum = sum(freq_alt)) %>%
  ungroup()

#######################################################################


## Calculate detection summary stats
# Calculate number of patients in each cohort
cohort.count <- ve303.dnameta_plot %>% 
  distinct(Subject.Number, .keep_all = TRUE) %>% 
  group_by(TRT.norm) %>% 
  summarise(N_sub = n())

# Calculate number of patients in each cohort by visit.name.norm
cohort.timepoint.count <- ve303.dnameta_plot %>% 
  group_by(TRT.norm, Visit.Name.Norm) %>% 
  distinct(Subject.Number, .keep_all = TRUE) %>% 
  summarise(N_analyzed = n())


# Calculate number of patients in each cohort by visit.name.norm
cohort.timepoint.count.num <- ve303.dnameta_plot %>% 
  group_by(TRT.norm, Day.in.treatment.norm,Visit.Name.Norm) %>% 
  distinct(Subject.Number, .keep_all = TRUE) %>% 
  summarise(N_analyzed = n())

# Calculate number of patients in each cohort by visit.name.norm
cohort.rec.count <- ve303.dnameta_plot %>% 
  distinct(Subject.Number, .keep_all = TRUE) 


#######################################################################
## Creating dataframes for plotting

# To plot zeroes in colonization data
ve303.dnameta_plot = ve303.dnameta_plot %>%
  mutate(freq_alt = if_else(detection_status != "Detected", 0, freq)) %>%
  group_by(Subject.Number, Visit.Name.Norm, organism) %>%
  mutate(.,freq_sum = sum(freq_alt)) %>%
  ungroup()


# To plot zeroes in colonization data
ve303.dnameta_plot.alt = ve303.dnameta_plot.alt %>%
  mutate(freq_alt = if_else(detection_status != "Detected", 0, freq)) %>%
  group_by(Subject.Number, Visit.Name.Norm, organism) %>%
  mutate(.,freq_sum = sum(freq_alt)) %>%
  ungroup()

## Total colonization by sample 
# DF used for plotting "All_VE303_colonization_detected_count_boxplot.pdf"
ve303.colonize.sum <- ve303.dnameta_plot %>% 
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


# VE303 colonization by treatment arm
ve303.colonize.det.stats <-  ve303.colonize.sum.alt %>% 
#  filter(detection_status %in% c("Detected", "Not detected")) %>% 
  group_by(TRT.norm, Visit.Name.Norm) %>% 
  summarise_at(vars(N_det_alt), funs(mean, median, sd))

ve303.colonize.det.stats.mean <- ve303.colonize.det.stats %>% 
  select(TRT.norm, Visit.Name.Norm,	mean) %>%  
  spread(Visit.Name.Norm, mean) 

 
# VE303 colonization by recurrence
ve303.colonize.det.recur.stats <-  ve303.colonize.sum.alt %>% 
#  filter(detection_status %in% c("Detected", "Not detected")) %>% 
  group_by(CDISENS2,TRT.norm) %>% 
  summarise_at(vars(N_det_alt), funs(mean, median, sd))

ve303.colonize.det.recur.stats.mean <- ve303.colonize.det.recur.stats %>% 
  select(TRT.norm, CDISENS2,	mean) %>%  
  spread(CDISENS2, mean) 

 
# Total by strain
ve303.colonize.sum.bystrain <- ve303.dnameta_plot.alt %>% 
  group_by(MGS.Sample.ID, Visit.Name.Norm, organism, detection_status) %>% 
  summarise(N_det = n()) %>% 
  ungroup() %>% 
  select(-Visit.Name.Norm) %>%
  left_join(.,dnameta, by = "MGS.Sample.ID") %>%
  left_join(., cohort.count, by = "TRT.norm") %>% 
  group_by(TRT.norm, MGS.Sample.ID, Visit.Name.Norm, organism, detection_status) %>% 
  summarise(N_det_tot = sum(N_det)) %>% 
  left_join(., cohort.count, by = "TRT.norm")

ve303.colonize.sum.bystrain.tp <- ve303.dnameta_plot.alt %>%
  arrange(Subject.Number, Visit.Name.Norm) %>% 
  group_by(Visit.Name.Norm, organism) %>%      # Select a single representative sample from each Subject for each Timepoint.Calc
  distinct(Subject.Number, .keep_all = TRUE) %>% 
  group_by(MGS.Sample.ID, Visit.Name.Norm, organism, detection_status) %>% 
  summarise(N_det = n()) %>% 
  ungroup() %>% 
  select(-Visit.Name.Norm) %>% 
  left_join(., dnameta, by = "MGS.Sample.ID") %>% 
  left_join(., cohort.count, by = "TRT.norm") %>% 
  group_by(TRT.norm, Visit.Name.Norm, organism, detection_status) %>% 
  summarise(N_det_tot = sum(N_det)) %>% 
  left_join(., cohort.timepoint.count, by = c("TRT.norm", "Visit.Name.Norm")) %>% 
  mutate(perc_subj = (100 * N_det_tot/N_analyzed))

## To plot zeroes
ve303.colonize.sum.bystrain.tp.alt = ve303.colonize.sum.bystrain.tp %>%
  mutate(perc_subj_alt = if_else(detection_status != "Detected", 0, perc_subj)) %>%
  group_by(TRT.norm, Visit.Name.Norm, organism) %>%
  summarise(.,perc_subj_alt_sum = sum(perc_subj_alt))

ve303.colonize.sum.bystrain.tp.det <- ve303.colonize.sum.bystrain.tp  %>% 
  mutate(perc_subj_alt = if_else(detection_status != "Detected", 0, perc_subj)) %>%
  #  filter(detection_status %in% c("Detected", "Not detected")) %>% 
  select(-N_det_tot, -N_analyzed) %>% 
  spread(Visit.Name.Norm, perc_subj_alt, fill = NA)

 
##########################################################################
##########################################################################

##   VE303 relative abundance - Detected events only 
##& using normalized marker depth adjusted for genome size and sequencing depth


#VE303 Abundance summary by cohort and timepoint
ve303.abund.sum <- ve303.dnameta_plot.alt %>% 
  group_by(organism, TRT.norm, Visit.Name.Norm) %>% 
  summarise_at(vars(est_relative_abundance_panel), funs(mean, median, sd)) %>% 
  ungroup() %>% 
  arrange(organism, Visit.Name.Norm, TRT.norm) %>% 
  select(organism, Visit.Name.Norm, TRT.norm, everything())

ve303.abund.mean <- ve303.abund.sum %>% 
  select(organism, Visit.Name.Norm, TRT.norm, mean) %>% 
  spread(Visit.Name.Norm, mean) 
 
ve303.abund.med <- ve303.abund.sum %>% 
  select(organism, Visit.Name.Norm, TRT.norm, median) %>% 
  spread(Visit.Name.Norm, median) %>% 
  arrange(TRT.norm)
 
ve303.abund.subj.total <- ve303.dnameta_plot.alt %>% 
  group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>% 
  summarise_at(vars(est_relative_abundance_panel), funs(sum)) %>% 
  ungroup() %>% 
  arrange(TRT.norm, Subject.Number, Visit.Name.Norm) %>% 
  select(TRT.norm, Subject.Number, everything())

ve303.abund.subj.total.time <- ve303.dnameta_plot.alt %>% 
  group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>% 
  summarise_at(vars(est_relative_abundance_panel), funs(sum)) %>% 
  ungroup() %>% 
  arrange(TRT.norm, Subject.Number, Visit.Name.Norm) %>% 
  select(TRT.norm, Subject.Number, everything())

# Total abundance by sample - df used to plot ve303 abundance
ve303.abund.subj.total.time <- ve303.dnameta_plot.alt %>% 
  group_by(MGS.Sample.ID, Visit.Name.Norm) %>% 
  summarise_at(vars(est_relative_abundance_panel), funs(sum)) %>%
  ungroup() %>% 
  select(-Visit.Name.Norm) %>% 
  left_join(., dnameta, by = "MGS.Sample.ID") %>% 
  left_join(., cohort.count, by = "TRT.norm")

ve303.abund.timepoint.total.mean <- ve303.abund.subj.total.time %>% 
  group_by(TRT.norm, Visit.Name.Norm) %>% 
  summarise_at(vars(est_relative_abundance_panel), funs(mean)) %>% 
  spread(Visit.Name.Norm, est_relative_abundance_panel)
 
ve303.abund.timepoint.total.median <- ve303.abund.subj.total.time %>% 
  group_by(TRT.norm, Visit.Name.Norm) %>% 
  summarise_at(vars(est_relative_abundance_panel), funs(median)) %>%  
  spread(Visit.Name.Norm, est_relative_abundance_panel)
 
##########################################################################