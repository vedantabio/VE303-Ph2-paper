
# VE303 Abundance summary by Subject and timepoint, dfs for plotting

## Working with safety data for timelines



ve303.abund.sum.subj.det <- ve303.dnameta.det.alt %>% 
  select(Subject.Number, organism, TRT.norm, Visit.Name.Norm, est_relative_abundance_panel, MGS.Sample.ID) %>%
  group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>% 
  spread(organism, est_relative_abundance_panel, fill=0) %>% # This fills in zeros for all events not detected - needed to calc mean
  gather(organism, est_relative_abundance_panel, 5: ncol(.)) %>% 
  group_by(Subject.Number, organism, TRT.norm, Visit.Name.Norm) %>% 
  mutate(subject_mean = mean(est_relative_abundance_panel)) %>%
  ungroup() %>% 
  arrange(TRT.norm, Subject.Number, Visit.Name.Norm, organism) %>% 
  select(-est_relative_abundance_panel) %>% 
  group_by(organism, TRT.norm, Visit.Name.Norm) %>% 
  mutate(cohort_abund_mean = mean(subject_mean), cohort_abund_sd = sd(subject_mean)) %>% 
  ungroup() %>% 
  group_by(Subject.Number, organism) %>%
  distinct(Visit.Name.Norm, .keep_all = TRUE)


# VE303 Abundance summary by Subject and timepoint, df for plotting INCLUDING zeros
ve303.abund.sum.subj.det.alt <- ve303.dnameta.det.alt %>% 
  select(Subject.Number, organism, TRT.norm, Visit.Name.Norm, Visit.Name,rec.diagnosis.Wk8, 
         est_relative_abundance_panel_alt, MGS.Sample.ID) %>%
  group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>% 
  spread(organism, est_relative_abundance_panel_alt, fill=0) %>% # This fills in zeros for all events not detected - needed to calc mean
  gather(organism, est_relative_abundance_panel_alt, 7: ncol(.)) %>% 
  group_by(Subject.Number, organism, TRT.norm, Visit.Name.Norm) %>% 
  mutate(subject_mean = mean(est_relative_abundance_panel_alt)) %>%
  ungroup() %>% 
  arrange(TRT.norm, Subject.Number, Visit.Name.Norm, organism) %>% 
  select(-est_relative_abundance_panel_alt) %>% 
  group_by(organism, TRT.norm, Visit.Name.Norm) %>% 
  mutate(cohort_abund_mean = mean(subject_mean), cohort_abund_sd = sd(subject_mean)) %>% 
  ungroup() %>% 
  group_by(Subject.Number, organism) %>%
  distinct(Visit.Name.Norm, .keep_all = TRUE)


# VE303 Abundance summary by Subject and timepoint, df for plotting INCLUDING zeros
ve303.abund.sum.subj.det.alt_Wk8R <- ve303.dnameta.det.alt %>% 
  select(Subject.Number, organism, TRT.norm, Visit.Name.Norm, Visit.Name,rec.diagnosis.Wk8, est_relative_abundance_panel_alt, MGS.Sample.ID) %>%
  group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>% 
  spread(organism, est_relative_abundance_panel_alt, fill=0) %>% # This fills in zeros for all events not detected - needed to calc mean
  gather(organism, est_relative_abundance_panel_alt, 7: ncol(.)) %>% 
  group_by(Subject.Number, organism, TRT.norm, Visit.Name.Norm) %>% 
  mutate(subject_mean = mean(est_relative_abundance_panel_alt)) %>%
  ungroup() %>% 
  arrange(TRT.norm, Subject.Number, Visit.Name.Norm, organism) %>% 
  select(-est_relative_abundance_panel_alt) %>% 
  group_by(organism, TRT.norm, Visit.Name.Norm) %>% 
  mutate(cohort_abund_mean = mean(subject_mean), cohort_abund_sd = sd(subject_mean)) %>% 
  ungroup() %>% 
  group_by(Subject.Number, organism) %>%
  distinct(Visit.Name.Norm, .keep_all = TRUE)


# VE303 Abundance summary by Subject and timepoint, df for plotting
ve303.abund.sum.subj.det.log <- ve303.dnameta.det.alt %>% 
  select(Subject.Number, organism, TRT.norm, Visit.Name.Norm, est_relative_abundance_panel, MGS.Sample.ID) %>%
  group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>% 
  spread(organism, est_relative_abundance_panel, fill=0) %>% # This fills in zeros for all events not detected - needed to calc mean
  gather(organism, est_relative_abundance_panel, 5: ncol(.)) %>%
  mutate_at(vars("est_relative_abundance_panel"), function(x) x*100) %>%
  mutate_at(vars("est_relative_abundance_panel"), function(x) x+0.01) %>%
  #mutate_at(vars("est_relative_abundance_panel"), function(x) log10(x))%>%
  group_by(Subject.Number, organism, TRT.norm, Visit.Name.Norm) %>% 
  mutate(subject_mean = mean(est_relative_abundance_panel)) %>%
  ungroup() %>% 
  arrange(TRT.norm, Subject.Number, Visit.Name.Norm, organism) %>% 
  select(-est_relative_abundance_panel) %>% 
  group_by(organism, TRT.norm, Visit.Name.Norm) %>% 
  mutate(cohort_abund_mean = mean(subject_mean), cohort_abund_sd = sd(subject_mean)) %>% 
  ungroup() %>% 
  group_by(Subject.Number, organism) %>%
  distinct(Visit.Name.Norm, .keep_all = TRUE)

# # VE303 ABSOLUTE Abundance summary by Subject and timepoint, df for plotting
 ve303.ABS.abund.sum.subj.det.alt <- ve303.dnameta.det.alt %>% 
    select(Subject.Number, organism, TRT.norm, Visit.Name.Norm, absolute_abund, MGS.Sample.ID) %>%
    group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>% 
    spread(organism, absolute_abund, fill=0) %>% # This fills in zeros for all events not detected - needed to calc mean
    gather(organism, absolute_abund, 5: ncol(.)) %>% 
    group_by(Subject.Number, organism, TRT.norm, Visit.Name.Norm) %>% 
    mutate(subject_mean = mean(absolute_abund)) %>%
    ungroup() %>% 
    arrange(TRT.norm, Subject.Number, Visit.Name.Norm, organism) %>% 
    select(-absolute_abund) %>% 
    group_by(organism, TRT.norm, Visit.Name.Norm) %>% 
    mutate(cohort_abund_mean = mean(subject_mean), cohort_abund_sd = sd(subject_mean)) %>% 
    ungroup() %>% 
    group_by(Subject.Number, organism) %>%
    distinct(Visit.Name.Norm, .keep_all = TRUE)


#Total VE303 Abundance summary by Subject and timepoint
ve303.abund.total303.subj.det <- ve303.dnameta.det.alt %>% 
  select(Subject.Number, organism, TRT.norm, Visit.Name.Norm, est_relative_abundance_panel, MGS.Sample.ID) %>%
  group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>% 
  spread(organism, est_relative_abundance_panel, fill=0) %>% # This fills in zeros for all events not detected - needed to calc mean
  gather(organism, est_relative_abundance_panel, 5: ncol(.)) %>% 
  group_by(Subject.Number, Visit.Name.Norm) %>% 
  mutate(ve303_sum = sum(est_relative_abundance_panel)) %>% 
  group_by(Subject.Number, organism, TRT.norm, Visit.Name.Norm) %>%
  mutate(subject_mean = mean(ve303_sum)) %>%
  ungroup() %>% 
  arrange(TRT.norm, Subject.Number, Visit.Name.Norm, organism) %>% 
  select(-est_relative_abundance_panel) %>% 
  group_by(organism, TRT.norm, Visit.Name.Norm) %>% 
  mutate(cohort_abund_mean = mean(subject_mean), cohort_abund_sd = sd(subject_mean)) %>% 
  ungroup() %>% 
  group_by(Subject.Number, organism) %>%
  distinct(Visit.Name.Norm, .keep_all = TRUE)

# TOTAL VE303 Abundance summary by Subject and timepoint, df for plotting INCLUDING zeros

ve303.abund.total303.subj.det.alt <- ve303.dnameta.det.alt %>% 
  select(Subject.Number, organism, TRT.norm,Visit.Name,Visit.Name.Norm, rec.diagnosis.Wk8,
         est_relative_abundance_panel_alt, MGS.Sample.ID) %>%
  group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>% 
  spread(organism, est_relative_abundance_panel_alt, fill=0) %>% # This fills in zeros for all events not detected - needed to calc mean
  gather(organism, est_relative_abundance_panel_alt, 7: ncol(.)) %>% 
  group_by(Subject.Number, Visit.Name.Norm) %>% 
  mutate(ve303_sum = sum(est_relative_abundance_panel_alt)) %>% 
  group_by(Subject.Number, organism, TRT.norm, Visit.Name.Norm) %>%
  mutate(subject_mean = mean(ve303_sum)) %>%
  ungroup() %>% 
  arrange(TRT.norm, Subject.Number,rec.diagnosis.Wk8, Visit.Name.Norm, Visit.Name,organism) %>% 
  select(-est_relative_abundance_panel_alt) %>% 
  group_by(organism, TRT.norm, Visit.Name.Norm) %>% 
  mutate(cohort_abund_mean = mean(subject_mean), cohort_abund_sd = sd(subject_mean)) %>% 
  ungroup() %>% 
  group_by(Subject.Number, organism) %>%
  distinct(Visit.Name.Norm, .keep_all = TRUE) #%>% 
#  select(-Visit.Name.Norm) %>% 
# left_join(., dnameta, by = "MGS.Sample.ID") 


#VE303 Abundance summary by cohort and timepoint
ve303.abund.sum.coh.det <- ve303.dnameta.det.alt %>% 
  group_by(organism, TRT.norm, Visit.Name.Norm) %>% 
  summarise_at(vars(est_relative_abundance_panel), funs(mean, median, sd)) %>% 
  ungroup() %>% 
  arrange(organism, Visit.Name.Norm, TRT.norm) %>% 
  select(organism, Visit.Name.Norm, TRT.norm, everything())

ve303.abund.mean.det <- ve303.abund.sum.coh.det %>% 
  select(organism, Visit.Name.Norm, TRT.norm, mean) %>% 
  spread(Visit.Name.Norm, mean) %>% 
  arrange(TRT.norm)
 
ve303.abund.med.det <- ve303.abund.sum.coh.det %>% 
  select(organism, Visit.Name.Norm, TRT.norm, median) %>% 
  spread(Visit.Name.Norm, median) %>% 
  arrange(TRT.norm)
 
# Total abundance by sample to include zero events

ve303.abund.subj.total.time.det <- ve303.dnameta.det.alt %>% 
  group_by(MGS.Sample.ID, Visit.Name.Norm) %>% 
  summarise_at(vars(est_relative_abundance_panel), funs(sum)) %>%
  ungroup() %>% 
  select(-Visit.Name.Norm) %>% 
  left_join(., dnameta, by = "MGS.Sample.ID") %>% 
  left_join(., cohort.count, by = "TRT.norm")


ve303.abund.subj.total.time.det.alt <- ve303.dnameta_plot %>% 
  mutate(est_relative_abundance_panel_alt = if_else(detection_status != "Detected", 0, est_relative_abundance_panel)) %>%
  group_by(Subject.Number, Visit.Name.Norm, organism) %>%
  mutate(.,est_relative_abundance_panel_alt_sum = sum(est_relative_abundance_panel_alt)) %>%
  ungroup()%>%
  group_by(MGS.Sample.ID, Visit.Name.Norm) %>% 
  summarise_at(vars(est_relative_abundance_panel_alt_sum), funs(sum)) %>%
  ungroup() %>% 
  select(-Visit.Name.Norm) %>% 
  left_join(., dnameta, by = "MGS.Sample.ID") %>% 
  left_join(., cohort.count, by = "TRT.norm") %>%
  left_join(., (ve303.colonize.sum.alt %>% 
                  select(MGS.Sample.ID, N_det_alt_sum, prop_det_alt_sum)) , by = "MGS.Sample.ID" )


ve303.abund.timepoint.total.mean.det <- ve303.abund.subj.total.time.det %>% 
  group_by(TRT.norm, Visit.Name.Norm) %>% 
  summarise_at(vars(est_relative_abundance_panel), funs(mean)) %>% 
  spread(Visit.Name.Norm, est_relative_abundance_panel)
 
ve303.abund.timepoint.total.median.det <- ve303.abund.subj.total.time.det %>% 
  group_by(TRT.norm, Visit.Name.Norm) %>% 
  summarise_at(vars(est_relative_abundance_panel), funs(median)) %>%  
  spread(Visit.Name.Norm, est_relative_abundance_panel)
 
#print("waiting for SPR for absolute abundance")

# # Total ABSOLUTE abundance by sample, df for plotting
 ve303.ABS.abund.subj.total.time.det.alt <- ve303.dnameta.det.alt %>%
   group_by(MGS.Sample.ID, Visit.Name.Norm) %>%
   summarise_at(vars(absolute_abund), funs(sum)) %>%
   ungroup() %>%
   select(-Visit.Name.Norm) %>%
   left_join(., dnameta, by = "MGS.Sample.ID") %>%
   left_join(., cohort.count, by = "TRT.norm")
 
 ve303.ABS.abund.timepoint.total.mean.det.alt <- ve303.ABS.abund.subj.total.time.det.alt %>%
   group_by(TRT.norm, Visit.Name.Norm) %>%
   summarise_at(vars(absolute_abund), funs(mean)) %>%
   spread(Visit.Name.Norm, absolute_abund)
 
 ################################################################################

