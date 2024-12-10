
## Working with safety data for timelines


## Add to VE303.dnameta and dnameta
dnameta  <- dnameta %>%
  left_join(.,adverse_post, by = "Subject.Number")

dnameta$adv.day.in.treatment.surv <- dnameta$day.in.treat.AE.start

dnameta[["adv.day.in.treatment.surv"]][is.na(dnameta[["adv.day.in.treatment.surv"]])] <- 168
dnameta[["AE.first.postdosing"]][is.na(dnameta[["AE.first.postdosing"]])] <- 168
dnameta$adverse.occur <- "No"

dnameta <- within(dnameta, {
  f <- Subject.Number %in% adverse_subs
  adverse.occur[f] <- 'Yes'
}) 


dnameta <- dnameta  %>% 
  rename_at(vars(ends_with(".x")),  ~str_replace(., "\\..$","") ) %>% 
  select_at(vars(-ends_with(".y")))


write_csv(dnameta , path = here(results, paste(Sys.Date(), AE_subset,"VE303_dna_metadata_45IA_AdverseEvents.csv", sep=" ")))


########################################################################################


ve303.dnameta  <- ve303.dnameta %>% 
  left_join(.,adverse_post, by = "Subject.Number")

ve303.dnameta$adv.day.in.treatment.surv <- ve303.dnameta$day.in.treat.AE.start

ve303.dnameta[["adv.day.in.treatment.surv"]][is.na(ve303.dnameta[["adv.day.in.treatment.surv"]])] <- 168
ve303.dnameta[["AE.first.postdosing"]][is.na(ve303.dnameta[["AE.first.postdosing"]])] <- 168


ve303.dnameta$adverse.occur <- "No"

ve303.dnameta <- within(ve303.dnameta, {
  f <- Subject.Number %in% adverse_subs
  adverse.occur[f] <- 'Yes'
}) 

ve303.dnameta  <- ve303.dnameta  %>% 
  rename_at(vars(ends_with(".x")),  ~str_replace(., "\\..$","") ) %>% 
  select_at(vars(-ends_with(".y")))

write_csv(ve303.dnameta, path = here(results, paste(Sys.Date(), 
                                                    "VE303_Ph2_Extended_Marker_Data_dna_metadata_recurrence_ABX_AEs_targeted_panel_RA.csv", sep=" ")))

########################################################################################

########################################################################################
