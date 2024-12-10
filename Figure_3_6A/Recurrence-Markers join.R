
## Working with recurrence data for timelines
#Add recurrence numbers to capture recurrences/patient

#Add recurrence numbers to capture recurrences/patient
rec = rec %>% 
  filter(rec.diagnosis=="Bad") %>%
  group_by(Subject.Number) %>%
  mutate(rec.no = n())

# Create a df that has the first recurrence event for each patient 
# such that patients with multiple recurrences 
# have just the one and their first recurrence noted in this df

rec.u = rec %>%
  group_by(Subject.Number) %>%
  
  arrange(rec.day.in.treatment) %>%
  
  mutate(id=1:n()) %>%
  filter(id==1) %>%
  select(Subject.Number, Start.Date.rec, rec.day.in.treatment, Primary.CDI.Recurrence, CDISENS1, CDISENS2, rec.diagnosis,
         rec.no)

write.csv(rec.u, file = here(results, paste(Sys.Date(), "VE303_recurrence_status_by_subject.csv", sep = " ")), row.names = FALSE)

### Day in treatment schedule up to which recurrences are clinically relvant and counted as recurrences ##
rec_clin_day <- 63

########################################################################################

## Add to VE303.dnameta and dnameta
dnameta = dnameta %>%
  left_join(.,rec.u, by = "Subject.Number")

dnameta[["rec.diagnosis"]][is.na(dnameta[["rec.diagnosis"]])] <- "Good"
dnameta[["rec.no"]][is.na(dnameta[["rec.no"]])] <- 0
dnameta[["CDISENS2"]][is.na(dnameta[["CDISENS2"]])] <- "No"

dnameta$collect.pre.post.rec = if_else(dnameta$rec.diagnosis=="Bad" & dnameta$rec.day.in.treatment >= dnameta$Day.in.treatment, "Pre.rec", 
                                       if_else(dnameta$rec.diagnosis=="Bad" & dnameta$rec.day.in.treatment < dnameta$Day.in.treatment, "Post.rec", "No.rec"))

## Add metadata columns to call all timepoints as either pre or post recurrence instead of no,pre,post
dnameta$collect.pre.post.rec_joined <- as.factor(plyr::mapvalues(dnameta$collect.pre.post.rec, from= c("No.rec"), to=c("Pre.rec")))

## Add metadata columns to call all recurrence post day 60/week8 as non-recurrent or good events
dnameta$collect.pre.post.rec.Wk8 <- dnameta$collect.pre.post.rec
dnameta$rec.diagnosis.Wk8 <- dnameta$rec.diagnosis
#dnameta <- within(dnameta, rec.diagnosis.early[rec.day.in.treatment >= 60] <- "Good")

dnameta <- within(dnameta, {
  f <- rec.day.in.treatment >= rec_clin_day
  rec.diagnosis.Wk8[f] <- 'Good'
  collect.pre.post.rec.Wk8[f] <- "No.rec"
}) 

########################################################################################

# Will be needed for survival analysis ##

dnameta$TRT.norm.rec <- paste( dnameta$TRT.norm, dnameta$rec.diagnosis.Wk8, sep = " ")
dnameta$rec.day.in.treatment.surv <- dnameta$rec.day.in.treatment
dnameta_by_subject <- split(dnameta, dnameta$Subject.Number)
dnameta <- do.call(rbind, lapply(dnameta_by_subject, show.me))

dnameta$collect.pre.post.rec.binary <- if_else(dnameta$collect.pre.post.rec_joined=="Pre.rec" , 0,1)

write_csv(dnameta, path = here(results, paste(Sys.Date(), "VE303_dna_metadata_79IA_recurrence.csv", sep=" ")))

# ########################################################################################
# ########################################################################################
# ########################################################################################


ve303.dnameta = ve303.dnameta %>% 
  left_join(.,rec.u, by = "Subject.Number")

ve303.dnameta[["rec.diagnosis"]][is.na(ve303.dnameta[["rec.diagnosis"]])] <- "Good"
ve303.dnameta[["rec.no"]][is.na(ve303.dnameta[["rec.no"]])] <- 0
ve303.dnameta[["CDISENS2"]][is.na(ve303.dnameta[["CDISENS2"]])] <- "No"

ve303.dnameta$collect.pre.post.rec = if_else(ve303.dnameta$rec.diagnosis=="Bad" & ve303.dnameta$rec.day.in.treatment >= 
                                               ve303.dnameta$Day.in.treatment, "Pre.rec",
                                             if_else(ve303.dnameta$rec.diagnosis=="Bad" & ve303.dnameta$rec.day.in.treatment < 
                                                       ve303.dnameta$Day.in.treatment, "Post.rec", "No.rec"))


## Add metadata columns to call all timepoints as either pre or post recurrence instead of no,pre,post
ve303.dnameta$collect.pre.post.rec_joined <- as.factor(plyr::mapvalues(ve303.dnameta$collect.pre.post.rec, from= c("No.rec"), to=c("Pre.rec")))

## Add metadata columns to call all recurrence post day 60/week8 as non-recurrent or good events
ve303.dnameta$collect.pre.post.rec.Wk8 <- ve303.dnameta$collect.pre.post.rec
ve303.dnameta$rec.diagnosis.Wk8 <- ve303.dnameta$rec.diagnosis
ve303.dnameta$collect.pre.post.rec.binary <- if_else(ve303.dnameta$collect.pre.post.rec_joined=="Pre.rec" , 0,1)

ve303.dnameta <- within(ve303.dnameta, {
  f <- rec.day.in.treatment >= rec_clin_day
  rec.diagnosis.Wk8[f] <- 'Good'
  collect.pre.post.rec.Wk8[f] <- "No.rec"
}) 

# Will be needed for survival analysis ##

ve303.dnameta$TRT.norm.pop <- paste(ve303.dnameta$TRT.norm, ve303.dnameta$Population, sep = " ")
ve303.dnameta$TRT.norm.rec <- paste(ve303.dnameta$TRT.norm, ve303.dnameta$rec.diagnosis.Wk8, sep = " ")
ve303.dnameta$rec.day.in.treatment.surv <- ve303.dnameta$rec.day.in.treatment
ve303.dnameta_by_subject <- split(ve303.dnameta, ve303.dnameta$Subject.Number)
ve303.dnameta <- do.call(rbind, lapply(ve303.dnameta_by_subject, show.me))
 

ve303.dnameta <- ve303.dnameta %>% 
  rename_at(vars(ends_with(".x")),  ~str_replace(., "\\..$","") ) %>% 
  select_at(vars(-ends_with(".y")))

write_csv(ve303.dnameta, path = here(results, paste(Sys.Date(), "VE303_Ph2_Extended_Marker_Data_dna_metadata_recurrence_targeted_panel_RA.csv", sep=" ")))
   
 
# ########################################################################################

