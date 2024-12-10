

## Import cleaned metadata
dnameta = read_csv(here("Input files/2021-04-01 VE303_dna_metadata_45IA.csv"))

## Filter to examine unique instances of Subject per timepoint/visit.Name.Norm
dnameta = dnameta %>% filter(analysis.45=="Yes")


## Import cleaned marker panel combined with metadata file
ve303.dnameta = read_csv(here("Input files/2021-04-01 VE303_Ph2_Extended_Marker_Data_dna_metadata_targeted_panel_RA.csv"))

## Filter to examine unique instances of Subject per timepoint/visit.Name.Norm
ve303.dnameta = ve303.dnameta %>% filter(analysis.45=="Yes")

## Make all detected events with est_relative_zero as Detected 0
ve303.dnameta$detection_status_alt = if_else(ve303.dnameta$est_relative_abundance == 0 & ve303.dnameta$detection_status== "Detected", "Detected 0", ve303.dnameta$detection_status)


### CODE TO COMBINE Extended marker panel with metadata file
#### Import updated  marker panel with targeted relative abundance estimation : TO BE COMBINED with metadata file
##ve303 = read_csv(here("Input files/2021-03-25 VE303_Ph2_Extended_Marker_Data.csv"))

# Remove out empty columns,truncate MGS.Sample.ID to allow merging of metadata with ve303 detection results
#ve303 = ve303 %>%
#        select(-c(barcode_id, subject_id, cohort_id, timepoint, date_collected)) %>%
#        mutate(MGS.Sample.ID = str_remove_all(filename, pattern = ".fastq.gz"))

## Truncating after first 11 string characters - as non standard sample labels. Doesn't work just to remove ".fastq.gz"
#ve303 = ve303 %>%
#        select(-c(barcode_id, subject_id, cohort_id, timepoint, date_collected)) %>%
#        mutate(MGS.Sample.ID = substr(filename, 1,11))

#ve303.dnameta = left_join(ve303, dnameta, by = "MGS.Sample.ID")

### Filter to examine unique instances of Subject per timepoint/visit.Name.Norm
#ve303.dnameta = ve303.dnameta %>% filter(analysis.45=="Yes")


Redundant_Samps <- c("BG1869923AB","BI5711623AB","BI5712023AB","BK1494223AB","BK1494623AB","BL5001523AB","BL5001723AB","BP3746023AB")

Redundant_filename <- c("BL5001723AB.fastq.gz" ,"BG1869923AB.fastq.gz" , "BL5001523AB.fastq.gz", 
                        "BP3746023AB.fastq.gz", "BI5711623AB.fastq.gz", "BI5712023AB.fastq.gz" ,
                        "BK1494223AB.fastq.gz", "BK1494623AB.fastq.gz")

Redundant_filename_keep <- c("BP3746023AB_S157_R1_001.fastq.gz", "BK1494623AB_S156_R1_001.fastq.gz" ,"BI5711623AB_S151_R1_001.fastq.gz" ,
                             "BI5712023AB_S152_R1_001.fastq.gz", "BL5001723AB_S154_R1_001.fastq.gz" ,"BG1869923AB_S150_R1_001.fastq.gz",
                             "BK1494223AB_S155_R1_001.fastq.gz","BL5001523AB_S153_R1_001.fastq.gz")


ve303.dnameta <- ve303.dnameta %>% filter( !( (MGS.Sample.ID %in% Redundant_Samps) & (filename %in% Redundant_filename) ))
ve303.dnameta$Visit.Name.Norm1 <-  as.factor(paste("Day", ve303.dnameta$Day.in.treatment.norm, sep = " ")) 
ve303.dnameta <- ve303.dnameta %>% mutate(Visit.Name.Norm1 = str_replace(Visit.Name.Norm1, "Day Screening", "Screening"))
ve303.dnameta$Visit.Name.Norm <- ve303.dnameta$Visit.Name.Norm1
ve303.dnameta$est_relative_abundance_panel <- ve303.dnameta$targeted_panel_relative_abundance

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# ## CoreBiome SPR - DNA tab taken, cleaned, included here
# dnameta = read_csv(here("Input files/2020-06-22 VE303-002 CoreBiome Metadata-OUT.csv"))
# #BP3831023AB recorded by DNA ID belong to patient 105501 received two samples on Day 119
# dnameta = filter(dnameta, DNA.ID != "BP3831023AB")
# 
# ## CoreBiome SPR - Omnigene Gut meta tab, cleaned and included here to add omnigene gut weights
# omnimeta = read_csv(here("Input files/VE303-002 CoreBiome Stool Processing Record 2020-06-19-Omnilog-INPUT.csv"))
# names(omnimeta) = make.names(names(omnimeta), unique = T)#remove spaces
# omnimeta = omnimeta %>%
#            select(Sample.ID, OmniGene.gut.tube.full.weight..g.)
# 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


rec = read_csv(here("Input files/2021-03-19 CDI_Recurrences_IA2_treatment_assignment-OUT.csv"))
source("Recurrence-Markers join.R")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


conmed = read_csv(here("Input files/2021-03-19 VE303-002_ConMeds_treatment_assignment-OUT.csv"))
relevant_metadata <- c("treat.day.len","day.in.treat.abx.start","Abx.group", "TRT.norm","CMINDC","Subject.Number")

source("Conmeds-Markers join.R")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


print("Select the AE file to be analyzed : Diarrhea AE or all AE")
if(event_type == "Recurrence"){ ##Placeholder, not requred if Recurrence analysis but still do computations
  adverse_evs <- all_adverse_evs}
if(event_type == "AE"){adverse_evs <- all_adverse_evs}
if(event_type == "Diarrhea AE"){adverse_evs <- diarrhea_adverse_evs}

adverse_evs$EVENT_num <- 1

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


print("Select the type of AE to be analyzed, and uncomment the relevant TX.REQ, RELATIONSHIP, SEVERITY, SERIOUS, ACTION.TAKEN filters accordingly")

# # # # # # # # # # # # 

# (Not Relevant) "ALL AEs" : all filters commented
# (Not Relevant) "TEAEs" : uncomment TX.REQ == "Yes"
# # # # # # # # # # # # 
# (Relevant) "Mod-Severe TEAEs" : uncomment TX.REQ == "Yes", SEVERITY != "GRADE 1 (MILD)"
#AE : There are 14 subjects who meet these conditions
#Diarrhea AE : There are 15 subjects who meet these conditions

# (Relevant) "Serious TEAEs" (Serious AEs shows the same patients) : uncomment  TX.REQ == "Yes", SERIOUS == "Yes"
#AE : There are 5 subjects who meet these conditions (3 if considering events to week 8)
#Diarrhea AE : There are 0 subjects who meet these conditions

# (Relevant) "Drug related TEAEs" : uncomment TX.REQ == "Yes" and RELATIONSHIP == "Related". 
#AE : There are 2 subjects who meet these conditions
#Diarrhea AE : There are 3 subjects who meet these conditions

# (Relevant) "Drug related Serious TEAEs" Drug related Serious AEs or TEAEs" : uncomment TX.REQ == "Yes", RELATIONSHIP == "Related",  SERIOUS == "Yes". 
#AE : There are 0 subjects who meet these conditions
#Diarrhea AE : There are 0 subjects who meet these conditions

# (Relevant) "Drug related Mod-Severe TEAEs" : uncomment TX.REQ == "Yes", RELATIONSHIP == "Related", SEVERITY != "GRADE 1 (MILD)". 
#AE : There are 1 subjects who meet these conditions
#Diarrhea AE : There are 2 subjects who meet these conditions

# (Relevant) "TEAEs Drug Withdrawal" TEAEs leading to drug withdrawal : uncomment TX.REQ == "Yes", ACTION.TAKEN == Drug Withdrawn" 
#AE : There are 0 subjects who meet these conditions
#Diarrhea AE : There are 5 subjects who meet these conditions

# # # # # # # # # # # # 
# (Not Relevant) Drug related Mod-Severe AEs = 3 : uncomment RELATIONSHIP == "Related", SEVERITY != "GRADE 1 (MILD)". 
#AE : There are 3 subjects who meet these conditions
#Diarrhea AE : There are 0 subjects who meet these conditions
# (Not Relevant) Drug related AEs : uncomment RELATIONSHIP == "Related"
#AE : There are 9 subjects who meet these conditions
#Diarrhea AE : There are 11 subjects who meet these conditions

# # # # # # # # # # # # 


# adverse_post <- adverse_evs  %>% 
#   filter((day.in.treat.AE.start > 0) ) %>% 
#   filter(TX.REQ == "Yes") %>%
#   #filter((RELATIONSHIP == "Related") ) %>% 
#   #filter((SEVERITY != "GRADE 1 (MILD)") ) %>% 
#   #filter((SERIOUS == "Yes") ) %>% 
#   #filter(ACTION.TAKEN == "Drug Withdrawn") %>% 
#   group_by(Subject.Number) %>%
#   mutate(AE.first.postdosing = min(day.in.treat.AE.start)) %>%
#   filter((AE.first.postdosing == day.in.treat.AE.start)) 

if(AE_subset == "Mod-Severe TEAEs"){adverse_post <- adverse_evs  %>% 
  filter((day.in.treat.AE.start > 0) & (day.in.treat.AE.start < 60) ) %>% 
  filter(TX.REQ == "Yes") %>%
  #filter((RELATIONSHIP == "Related") ) %>% 
  filter((SEVERITY != "GRADE 1 (MILD)") ) %>% 
  #filter((SERIOUS == "Yes") ) %>% 
  #filter(ACTION.TAKEN == "Drug Withdrawn") %>% 
  group_by(Subject.Number) %>%
  mutate(AE.first.postdosing = min(day.in.treat.AE.start)) %>%
  ungroup() %>%
  group_by(Subject.Number, day.in.treat.AE.start) %>%
  mutate(.,events_sum_time = sum(EVENT_num)) %>%
  ungroup() %>%
  group_by(Subject.Number) %>%
  mutate(.,events_sum = sum(EVENT_num)) %>%
  ungroup()}



if(AE_subset == "Serious TEAEs"){adverse_post <- adverse_evs  %>% 
  filter((day.in.treat.AE.start > 0) & (day.in.treat.AE.start < 60) ) %>% 
  filter(TX.REQ == "Yes") %>%
  #filter((RELATIONSHIP == "Related") ) %>% 
  #filter((SEVERITY != "GRADE 1 (MILD)") ) %>% 
  filter((SERIOUS == "Yes") ) %>% 
  #filter(ACTION.TAKEN == "Drug Withdrawn") %>% 
  group_by(Subject.Number) %>%
  mutate(AE.first.postdosing = min(day.in.treat.AE.start)) %>%
  ungroup() %>%
  group_by(Subject.Number, day.in.treat.AE.start) %>%
  mutate(.,events_sum_time = sum(EVENT_num)) %>%
  ungroup() %>%
  group_by(Subject.Number) %>%
  mutate(.,events_sum = sum(EVENT_num)) %>%
  ungroup()}



if(AE_subset == "TEAEs Drug Withdrawal"){adverse_post <- adverse_evs  %>% 
  filter((day.in.treat.AE.start > 0) & (day.in.treat.AE.start < 60) ) %>% 
  filter(TX.REQ == "Yes") %>%
  #filter((RELATIONSHIP == "Related") ) %>% 
  #filter((SEVERITY != "GRADE 1 (MILD)") ) %>% 
  #filter((SERIOUS == "Yes") ) %>% 
  filter(ACTION.TAKEN == "Drug Withdrawn") %>% 
  group_by(Subject.Number) %>%
  mutate(AE.first.postdosing = min(day.in.treat.AE.start)) %>%
  ungroup() %>%
  group_by(Subject.Number, day.in.treat.AE.start) %>%
  mutate(.,events_sum_time = sum(EVENT_num)) %>%
  ungroup() %>%
  group_by(Subject.Number) %>%
  mutate(.,events_sum = sum(EVENT_num)) %>%
  ungroup()}


adverse_post <- distinct(adverse_post, Subject.Number, .keep_all = TRUE)
adverse_subs <- adverse_post$Subject.Number %>% unique()

adverse_post
adverse_post$Subject.Number %>% unique()
print(length(unique(adverse_post$Subject.Number)))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


source("Safety-Markers join.R")


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

### note: 
# f = relative abundance
# M = mass of stool 
# V = volume suspension (2 ml)
# Vh = volume homogeneizer (0.5)
# DNA = mass of DNA from the homogeneized stool
# absolute = f * DNA / M * V / Vh
# Weight of tube = 13.69 g

dnameta = dnameta %>% 
  mutate(weight.of.collected.sample..g. = OmniGene.gut.tube.full.weight..g. - 13.69)

ve303.dnameta = ve303.dnameta %>% 
  mutate(weight.of.collected.sample..g. = OmniGene.gut.tube.full.weight..g. - 13.69)

ve303.dnameta <- ve303.dnameta %>% 
  mutate(dna.norm = (DNA.yield.ug/1e9)/ weight.of.collected.sample..g.) %>% 
  mutate(absolute_abund = est_relative_abundance_panel * dna.norm * (2/0.50))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


#Add frequency variable for all rows - this gives every detection/non-det result a value of 1
ve303.dnameta <- data.frame(ve303.dnameta, freq = 1)

# Change from character class to factor
ve303.dnameta$detection_status <- factor(ve303.dnameta$detection_status, levels = c("Detected", "Insufficient data", "Not detected"))


#######################################################################

# To plot zeroes in colonization data
ve303.dnameta = ve303.dnameta %>%
  mutate(freq_alt = if_else(detection_status != "Detected", 0, freq)) %>%
  group_by(Subject.Number, Visit.Name.Norm, organism) %>%
  mutate(.,freq_sum = sum(freq_alt)) %>%
  ungroup()

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


sigtaxa <- c("VE303-02", "VE303-03", "VE303-05" ,"VE303-08")


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
ve303.dnameta.det.alt$events_sum[is.na(ve303.dnameta.det.alt$events_sum)] <- 0

events_poisson_input <- ve303.dnameta.det.alt %>% 
  distinct(Subject.Number, .keep_all = TRUE) %>% 
  select(Subject.Number, TRT.norm, events_sum)

#######################################################################

KM_cols1 <- c("MGS.Sample.ID", "organism",  "TRT.norm",  "Subject.Number", "Visit.Name", "Visit.Name.Norm","Day.in.treatment",
              "est_relative_abundance_panel_alt","collect.pre.post.rec.binary", "rec.diagnosis", "rec.diagnosis.Wk8", "collect.pre.post.abx", "rec.day.in.treatment", "rec.day.in.treatment.surv", "adv.day.in.treatment.surv", "adverse.occur", "events_sum")

df_plot_KM_indiv <- split(ve303.dnameta.det.alt[KM_cols1], ve303.dnameta.det.alt[KM_cols1]$organism)
if(plot_tax == "indiv") { 
  df_plot_KM_total <- df_plot_KM_indiv}


#######################################################################


KM_cols <- c("TRT.norm",  "Subject.Number", "Visit.Name.Norm","Visit.Name", "Day.in.treatment",
             "est_relative_abundance_panel_alt_sum","MGS.Sample.ID", "collect.pre.post.rec.binary",
             "rec.diagnosis", "rec.diagnosis.Wk8", "collect.pre.post.abx", "rec.day.in.treatment",
             "rec.day.in.treatment.surv","adv.day.in.treatment.surv", "adverse.occur", "events_sum")

if(plot_tax == "Sig4 Sum") { df_plot_KM_total <- ve303.abund.subj.total.time.det.alt_sig4[KM_cols] } 
if(plot_tax == "Total VE303") { df_plot_KM_total <- ve303.abund.subj.total.time.det.alt[KM_cols] } 


#######################################################################
df_plot_KM_total$rec.binary <- if_else(df_plot_KM_total$rec.diagnosis.Wk8 == "Bad", 
                                       1, 0)

df_plot_KM_total$adverse.binary <- if_else(df_plot_KM_total$adverse.occur == "Yes", 
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
if(plot_coh == "Dosed and Placebo") {df_plot_KM_total <- df_plot_KM_total }

df_plot_KM_total_subj <- df_plot_KM_total %>% 
  filter(collect.pre.post.abx == "Pre.abx") %>%
  filter(Day.in.treatment < 16) %>%
  group_by(Subject.Number) %>%
  mutate(.,est_relative_abundance_panel_alt_sum_Alltime = sum(est_relative_abundance_panel_alt_sum)) %>%
  ungroup() %>%
  distinct(Subject.Number, .keep_all = TRUE) %>%
  select(TRT.norm,Subject.Number,est_relative_abundance_panel_alt_sum_Alltime, rec.diagnosis, rec.diagnosis.Wk8, rec.day.in.treatment, rec.day.in.treatment.surv, adverse.occur,adv.day.in.treatment.surv,trt.binary, rec.binary, adverse.binary, events_sum)


df_plot_KM_total_subj <- df_plot_KM_total_subj[order(-df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime),]

#mthresh <- as.numeric(quantile(df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime)[4])
mthresh <- as.numeric(median(df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime))

as.numeric(quantile(df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime))

#######################################################################

df_plot_KM_total_subj$Engraftment <- if_else(df_plot_KM_total_subj$est_relative_abundance_panel_alt_sum_Alltime > mthresh, 
                                             "High Engraft", "Low Engraft")
df_plot_KM_total_subj$Engraftment.binary <- if_else(df_plot_KM_total_subj$Engraftment == "High Engraft", 
                                                    1, 0)

df_plot_KM_total_subj$events_sum[is.na(df_plot_KM_total_subj$events_sum)] <- 0


#######################################################################



ggplot(df_plot_KM_total_subj , aes(x=adverse.occur , y=0.01+100*est_relative_abundance_panel_alt_sum_Alltime, fill=adverse.occur, color=adverse.occur)) + 
  
  geom_boxplot(alpha=0.6,  outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.3, position = position_jitterdodge(), show.legend = FALSE) +
  theme_bw() + 
  #facet_grid(~Engraftment, scales="free_x", space = "free") + 
  labs(title= "Exposure-Safety Summary", y="Engraftment", x="Subjects with TEAEs", fill="adverse.occur") + 
  scale_fill_manual(values = rec.cols) + 
  scale_color_manual(values = rec.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) , strip.text = element_text(size = 13)) + 
  #  scale_y_continuous(limits = c(0,0.5) )
  #scale_y_continuous(limits = c(0,0.5) )
  scale_y_log10()
ggsave(filename = here(results, paste(Sys.Date(),  plot_coh, plot_tax, event_type,AE_subset,"Engraftment across AE groups.pdf", sep = " ")),
       width = 6, height = 5)

#######################################################################

ggplot(df_plot_KM_total_subj  , aes(x=events_sum , fill=Engraftment, color=Engraftment)) + 
  #ggplot(df_plot_KM_total_subj %>% filter(TRT.norm != "Placebo") , aes(x=events_sum , fill=Engraftment, color=Engraftment)) + 
  geom_histogram( binwidth = 0.8, position = "dodge") +
  #geom_boxplot(alpha=0.6,  outlier.size = 0) +
  #geom_point(pch = 21, size = 2, alpha = 0.5, position = position_jitterdodge(), show.legend = FALSE) +
  theme_bw() + 
  facet_grid(~Engraftment, scales="free_x", space = "free") + 
  labs(title= "Exposure-Safety Summary", x="N adverse events (Week 8)", y="N (Subjects)" ) + 
  scale_fill_manual(values = eng.cols) + 
  scale_color_manual(values = eng.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) , strip.text = element_text(size = 13)) + 
  scale_x_continuous(limits = c(-1,7) )
ggsave(filename = here(results, paste(Sys.Date(),  plot_coh, plot_tax, event_type,AE_subset,"Histogram Eng Exposure-Safety across AE groups.pdf", sep = " ")),
       width = 6, height = 5)


#######################################################################

ggplot(df_plot_KM_total_subj %>% distinct(Subject.Number, .keep_all = TRUE), aes(x=events_sum , 
                                                                                 fill=TRT.norm, color=TRT.norm)) + 
  geom_histogram(binwidth = 0.5, position = "dodge") +
  #geom_boxplot(alpha=0.6,  outlier.size = 0) +
  #geom_point(pch = 21, size = 2, alpha = 0.5, position = position_jitterdodge(), show.legend = FALSE) +
  theme_bw() + 
  #facet_grid(~TRT.norm, scales="free_x", space = "free") + 
  labs(title= "Exposure-Safety Summary", x="# Adverse Events till Week 8", y="N (Subjects)" ) + 
  scale_fill_manual(values = coh.cols) + 
  scale_color_manual(values = coh.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) ,     axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) , strip.text = element_text(size = 13) , legend.title = element_blank()) + 
  #scale_x_continuous(limits = c(-1,7) )
  ggsave(filename = here(results, paste(Sys.Date(),  plot_coh, plot_tax, event_type,AE_subset,"Histogram Cohort-Safety across AE groups.pdf", sep = " ")),
         width = 6, height = 5)

#######################################################################


d1 <- df_plot_KM_total_subj #%>% filter(TRT.norm != "Placebo")
poisson1 <- glm(events_sum ~ TRT.norm + Engraftment , family="poisson", data=d1)

summary(poisson1)

## Capture eng Cox output including hazard ratios


poisson_trt <- glm(events_sum ~ TRT.norm  , family="poisson", data=d1)
capture.output(summary(poisson_trt), file=here(results, paste(Sys.Date(), plot_coh, plot_tax, event_type , AE_subset,
                                                           "TRT Poisson Model_Exposure-Safety summary.doc", sep=" ")) )

poisson_full <- glm(events_sum ~ TRT.norm + Engraftment , family="poisson", data=d1)
capture.output(summary(poisson_full), file=here(results, paste(Sys.Date(), plot_coh, plot_tax, event_type , AE_subset, 
                                                           "TRT-Eng Poisson Model_Exposure-Safety summary.doc", sep=" ")) )

poisson_eng <- glm(events_sum ~ Engraftment , family="poisson", data=d1)
capture.output(summary(poisson_eng), file=here(results, paste(Sys.Date(), plot_coh, plot_tax, event_type , AE_subset, 
                                                           "Eng Poisson Model_Exposure-Safety summary.doc", sep=" ")) )

#######################################################################



d1 <- df_plot_KM_total_subj 

tr <- length(unique(d1$TRT.norm))
ads <- length(unique(d1$adverse.occur))

if(ads < 2){ glm1 <- glm(est_relative_abundance_panel_alt_sum_Alltime ~ TRT.norm , family="gaussian", data=d1)}
if(!(ads < 2)){ glm1 <- glm(est_relative_abundance_panel_alt_sum_Alltime ~ TRT.norm + adverse.occur , family="gaussian", data=d1)}
summary(glm1)

#glm1 <- glm(est_relative_abundance_panel_alt_sum_Alltime ~ TRT.norm + adverse.occur , family="gaussian", data=d1)
## Capture eng Cox output including hazard ratios
#glm1 <- glm(events_sum ~ TRT.norm  , family="gaussian", data=d1)

capture.output(summary(glm1), file=here(results, paste(Sys.Date(), plot_coh, plot_tax, event_type , AE_subset,
                                                       "TRT-Adverse GLM Model_Exposure-Safety summary.doc", sep=" ")) )

#######################################################################



save_result <- list(GLM_adverse_eng=tidy(glm1),
                    Poisson_engraft = tidy(poisson_eng),
                    Poisson_trt = tidy(poisson_trt),
                    Poisson_trt_engraft = tidy(poisson_full))

## Writes to multi sheet excel file ###
writexl::write_xlsx(save_result, here(results, paste(Sys.Date(), plot_coh, plot_tax, event_type , AE_subset, 
                                                     "Poisson GLM Model summary.xlsx", sep=" ")) )



#######################################################################


ggplot(df_plot_KM_total_subj %>% distinct(Subject.Number, .keep_all = TRUE), aes(x=TRT.norm , y=events_sum**(0.5), 
                                                                                 fill=TRT.norm, color=TRT.norm)) + 
  geom_boxplot(alpha=0.6,  outlier.size = 0) +
  geom_point(pch = 21, size = 2, alpha = 0.5, position = position_jitterdodge(), show.legend = FALSE) +
  theme_bw() + 
  #facet_grid(~Engraftment, scales="free_x", space = "free") + 
  labs(title= "Exposure-Safety Summary", y="N adverse events (Week 8)", x="Treatment") + 
  scale_fill_manual(values = coh.cols) + 
  scale_color_manual(values = coh.cols) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5),axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) , 
        strip.text = element_text(size = 13), legend.title = element_blank()) + 
  #  scale_y_continuous(limits = c(0,0.5) )
  ggsave(filename = here(results, paste(Sys.Date(),  plot_coh, plot_tax, event_type,AE_subset,"Cohort-Safety across AE groups.pdf", sep = " ")),
         width = 6, height = 5)

#######################################################################


#coh1.cols <- c("Placebo" = "#66c2a5", "VE303 Low Dose" = "#fc8d62", "VE303 High Dose" = "#8da0cb")
coh1.cols <- c("Placebo" = "#1b9e77", "VE303 Low Dose" = "#7570b3", "VE303 High Dose" = "#d95f02")


xps <- cbind(rankord = as.numeric(rownames(df_plot_KM_total_subj)), (df_plot_KM_total_subj)) 
Nsamp <- length(xps$rankord)


ggplot(xps, aes(x = rankord, y = 0.0001+est_relative_abundance_panel_alt_sum_Alltime,fill=TRT.norm,
                color=TRT.norm )) + geom_bar(stat="identity",  width=0.6)  + ylim(0.0,0.2) +
  
  geom_vline(xintercept = Nsamp/2 , linetype="longdash", color="black") +
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

ggsave(filename = here(results, paste(Sys.Date(),plot_coh, plot_tax, "BARS Engraftment distribution.pdf", sep = " ")), 
       width = 7, height = 4)
#######################################################################

#coh1.cols <- c("Placebo" = "#66c2a5", "VE303 Low Dose" = "#fc8d62", "VE303 High Dose" = "#8da0cb")


xps <- cbind(rankord = as.numeric(rownames(df_plot_KM_total_subj)), (df_plot_KM_total_subj)) 
Nsamp <- length(xps$rankord)


ggplot(xps, aes(x = rankord, y = 0.0001+est_relative_abundance_panel_alt_sum_Alltime,fill=TRT.norm,
                color=TRT.norm )) + geom_point(size=2.2)  + ylim(0.0,0.4) +
  
  geom_vline(xintercept = Nsamp/2 , linetype="longdash", color="black") +
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

ggsave(filename = here(results, paste(Sys.Date(),plot_coh, plot_tax, "Engraftment curve distribution.pdf", sep = " ")), 
       width = 7, height = 4)

#######################################################################
