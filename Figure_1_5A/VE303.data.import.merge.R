
## Metadata with 1. Edited sample ID for placebo patients that were accidentally excluded due to redundant label names. 
#2. pre-post recurrence information 3. stool processing record and 4. targeted marker panel information added


## Metadata with 1. Edited sample ID for placebo patients that were accidentally excluded due to redundant label names. 
#2. pre-post recurrence information 3. stool processing record and 4. targeted marker panel information added

## Import demographics data
demographics = read.csv(here("Input files/2022-10-22 VE303-002_demographics.csv")) %>% 
  mutate(Subject.Number = as.numeric(Subject)) %>%
  rename(n_prevrecs = Number.of.CDI.Recurrences.Prior.to.and.Including.Current.Episode.at.Baseline,
         BMI = BMI..kg.m2.,
         Age = Age..yrs.) %>%
  select(Subject.Number, Age, BMI, n_prevrecs)

# Merging with demographics data with age, BMI, and number of previous recurrences including baseline recurrence
dnameta <- left_join(dnameta, demographics, by = "Subject.Number")


## Add treatment days to Screening events so that it may be plotted in the cohort diagram
dnameta$Day.in.treatment = if_else(dnameta$Day.in.treatment == "Day14/EOT", "14", dnameta$Day.in.treatment)
dnameta = dnameta %>%
  rename(Day.in.treatment.norm = Day.in.treatment) %>%
  mutate(., Day.in.treatment = ifelse(Day.in.treatment.norm %in% c("Screening"), (Collection.Date-Begin.Date)+1, 
                                      Day.in.treatment.norm)) %>% filter(!is.na(Day.in.treatment))

dnameta$Day.in.treatment <- as.numeric(dnameta$Day.in.treatment)
dnameta$Visit.Name = if_else(dnameta$Visit.Name == "Day14/EOT", "Day 14", dnameta$Visit.Name)
dnameta$Visit.Name.Norm = if_else(dnameta$Visit.Name.Norm == "Day14/EOT", "Day 14", dnameta$Visit.Name.Norm)

dnameta = left_join(dnameta, study_pops, by = "Subject.Number")

### No Treatment group information available for this subject ###########
dnameta <- dnameta %>% filter(!Subject.Number == 105005)

### Exclude this subject due to data integrity issues. 
# They are in VE303 High Dose, Non Recurrent
if(analysis_78 == "yes"){dnameta <- dnameta %>% filter(!Subject.Number == 100501)}
if(analysis_78 == "no"){dnameta <- dnameta }
  
########################################################################################


### CODE TO COMBINE Extended marker panel with metadata file
#### Import updated  marker panel with targeted relative abundance estimation : TO BE COMBINED with metadata file
#ve303 = read_csv(here("Input files/2021-10-06 VE303_Ph2_Extended_Marker_Data.csv"))


# Remove out empty columns,truncate MGS.Sample.ID to allow merging of metadata with ve303 detection results
## Truncating after first 11 string characters - as non standard sample labels. Doesn't work just to remove ".fastq.gz"
ve303 = ve303 %>%
  select(-c(barcode_id, subject_id, cohort_id, timepoint, date_collected)) %>%
  mutate(MGS.Sample.ID = substr(filename, 1,11)) %>%
  #       mutate(MGS.Sample.ID = str_remove_all(filename, pattern = ".fastq.gz"))
  filter(!grepl("ATCC",MGS.Sample.ID)) ## Filter out ATCC samples

 
########################################################################################

###### check that the VE303 sequencing files from OC have the same sample IDs
## as the merged Sequencing manifest/Stool Processing Record file

meta_S <- dnameta$MGS.Sample.ID %>% unique()
seq_S <- ve303$MGS.Sample.ID %>% unique()


#diff_S <- setdiff(meta_S, seq_S)
diff_S <- setdiff(seq_S, meta_S)
both_S <- intersect(meta_S, seq_S)

 

### Removing samples from the merged sequencing sample metadata that are NOT in the OC data file
dnameta <- dnameta %>% filter(MGS.Sample.ID %in% both_S)

########################################################################################

ve303.dnameta = left_join(ve303, dnameta, by = "MGS.Sample.ID")
ve303.dnameta$Day.in.treatment <- as.numeric(ve303.dnameta$Day.in.treatment)
 
########################################################################################

## Make all detected events with est_relative_zero as Detected 0
ve303.dnameta$detection_status_alt = if_else(ve303.dnameta$est_relative_abundance == 0 & 
                                               ve303.dnameta$detection_status== "Detected", 
                                             "Detected 0", ve303.dnameta$detection_status)


########################################################################################

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
 
########################################################################################
########################################################################################

## Merging VE303 subject level data and recurrence data 
#### Uses the file rec = read_csv(here("Input files/2021-10-05 CDI_Recurrences_IA2_treatment_assignment-OUT.csv"))
source("Recurrence-Markers join.R")


########################################################################################

## Merging VE303 subject level data and ConMed data to establish time of abx dosing and interaction with recurrence
#### Uses the file conmed = read_csv(here("Input files/2021-10-05 VE303-002_ConMeds_treatment_assignment-OUT.csv"))
source("Conmeds-Markers join.R")



## Absolute abundance calculation : ADD BACK IN WHEN RECEIVE STOOL PROCESSING RECORD
### note: 
# f = relative abundance
# M = mass of stool 
# V = volume suspension (2 ml)
# Vh = volume homogeneizer (0.5 ml)
# DNA = mass of DNA from the homogeneized stool
# absolute = f * DNA / M * V / Vh
# Weight of tube = 13.69 g

print("Adding Absolute abundance calculation using SPR. Converting micro g DNA to g")

dnameta = dnameta %>% 
   mutate(weight.of.collected.sample..g. = OmniGene.gut.tube.full.weight..g. - 13.69)
 
ve303.dnameta = ve303.dnameta %>% mutate(weight.of.collected.sample..g. = OmniGene.gut.tube.full.weight..g. - 13.69)
 
ve303.dnameta <- ve303.dnameta %>% 
   mutate(dna.norm = (DNA.yield.ug/1e6)/ weight.of.collected.sample..g.) %>% 
   mutate(absolute_abund = est_relative_abundance_panel * dna.norm * (2/0.50))


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
 