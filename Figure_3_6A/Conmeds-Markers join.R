

relevant_metadata <- c("treat.day.len","day.in.treat.abx.start","Abx.group", 
                       "TRT.norm","CMINDC","Subject.Number")


########################################################################################


conmed_pre <- conmed[relevant_metadata] %>% 
  filter(day.in.treat.abx.start < 0 ) %>%
  group_by(Subject.Number) %>%
  mutate(abx.treat.duringdosing = max(day.in.treat.abx.start)) %>%
  mutate(treat.day.len.abx = treat.day.len ) %>%
  filter(abx.treat.duringdosing == day.in.treat.abx.start)

conmed_post <- conmed[c("day.in.treat.abx.start", "treat.day.len","Subject.Number")] %>% 
  filter((day.in.treat.abx.start > 0) ) %>%
  group_by(Subject.Number) %>%
  mutate(abx.firsttreat.postdosing = min(day.in.treat.abx.start)) %>%
  filter((abx.firsttreat.postdosing == day.in.treat.abx.start)) %>%
  mutate(treat.day.len.ftpd.abx = treat.day.len) 


conmed_excep <- conmed[c("day.in.treat.abx.start", "treat.day.len","Subject.Number")] %>% 
  filter((day.in.treat.abx.start < 0) &  (day.in.treat.abx.start + treat.day.len > 2) ) %>%
  group_by(Subject.Number) %>%
  mutate(abx.firsttreat.postdosing = 1) %>%
  mutate(abx.exceptions = abx.firsttreat.postdosing) 

relevant_metadata <- c("treat.day.len.abx","day.in.treat.abx.start","Abx.group", 
                       "TRT.norm","CMINDC","Subject.Number")

conmed_pre_post <- distinct(conmed_pre[relevant_metadata] %>%
                              left_join(conmed_post[c("abx.firsttreat.postdosing", "Subject.Number")], 
                                        by = "Subject.Number")) %>%
  left_join(conmed_excep[c("abx.exceptions", "Subject.Number")], by = "Subject.Number")

conmed_pre_post$abx.firsttreat.postdosing <- if_else(!is.na(conmed_pre_post$abx.exceptions), 
                                                     conmed_pre_post$abx.exceptions,
                                                     conmed_pre_post$abx.firsttreat.postdosing )



write.csv(conmed_pre_post, file = here(results, paste(Sys.Date(), "VE303_first_postdose_abx_by_subject.csv", sep = " ")), row.names = FALSE)

########################################################################################

conmed_pre_post$Abx.group_joined_M <- plyr::mapvalues(conmed_pre_post$Abx.group, 
                                                      from= c("Vancomycin", "Other",
                                                           "Metronidazole",  "Fidaxomicin"), 
                              to=c("Vanco_Metro","Vanco_Metro","Vanco_Metro","Fidaxo"))

## Add to VE303.dnameta and dnameta
dnameta = dnameta %>%
  left_join(.,conmed_pre_post, by = "Subject.Number")


#dnameta$collect.pre.post.abx = if_else((dnameta$abx.firsttreat.postdosing >= dnameta$Day.in.treatment) | (is.na(dnameta$abx.firsttreat.postdosing)),
#                                       "Pre.abx", "Post.abx")

dnameta$collect.pre.post.abx = if_else((dnameta$abx.firsttreat.postdosing > dnameta$Day.in.treatment) |
                                         (dnameta$abx.firsttreat.postdosing == dnameta$Day.in.treatment) |
                                         (is.na(dnameta$abx.firsttreat.postdosing)),
                                       "Pre.abx", "Post.abx")

dnameta <- dnameta %>% 
  rename_at(vars(ends_with(".x")),  ~str_replace(., "\\..$","") ) %>% 
  select_at(vars(-ends_with(".y")))

#Adding metadata variable for subjects with a pre-screening start of antibiotics
dnameta_1 <- dnameta %>% filter(Visit.Name == "Screening" | Day.in.treatment < 0) %>% 
  distinct(Subject.Number,  .keep_all = TRUE) %>%
  mutate(Screening_abx_check = if_else(day.in.treat.abx.start < Day.in.treatment-1, 
                                       "On-Abx Screening","True Screening")) %>%
  select(Subject.Number, Screening_abx_check)

dnameta <- dnameta %>%   
  left_join(dnameta_1, by = c("Subject.Number") ) 



########################################################################################


ve303.dnameta <- ve303.dnameta %>%
  left_join(.,conmed_pre_post, by = "Subject.Number")


#ve303.dnameta$collect.pre.post.abx = if_else((ve303.dnameta$abx.firsttreat.postdosing >=  ve303.dnameta$Day.in.treatment) | 
#                                               (is.na(ve303.dnameta$abx.firsttreat.postdosing)), "Pre.abx", "Post.abx")

ve303.dnameta$collect.pre.post.abx = if_else((ve303.dnameta$abx.firsttreat.postdosing >  ve303.dnameta$Day.in.treatment) | 
                                               (ve303.dnameta$abx.firsttreat.postdosing ==  ve303.dnameta$Day.in.treatment) |
                                               (is.na(ve303.dnameta$abx.firsttreat.postdosing)), "Pre.abx", "Post.abx")

ve303.dnameta <- ve303.dnameta %>% 
  rename_at(vars(ends_with(".x")),  ~str_replace(., "\\..$","") ) %>% 
  select_at(vars(-ends_with(".y")))


ve303.dnameta_1 <- ve303.dnameta %>% filter(Visit.Name == "Screening" | Day.in.treatment < 0) %>% 
  distinct(Subject.Number,  .keep_all = TRUE) %>%
  mutate(Screening_abx_check = if_else(day.in.treat.abx.start < Day.in.treatment-1, 
                                       "On-Abx Screening","True Screening")) %>%
  select(Subject.Number, Screening_abx_check)

ve303.dnameta <- ve303.dnameta %>%   
  left_join(ve303.dnameta_1, by = c("Subject.Number") ) 




########################################################################################

