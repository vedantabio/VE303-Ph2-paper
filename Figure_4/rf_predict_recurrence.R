library(here)

# Highlight the recurrence data for each subject:

# Create a directory to store results from the analysis
mainDir <- here("recurrence_model_performance")
dir.create(file.path(mainDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,sep="/")


library(gtools)
library(tidyr)
library(phyloseq)
library(vita)
library(Boruta)
library(randomForest)
library(phyloseq)
library(tidyverse)



#############################Import ############################################


# Import functions to read data from various sources:

# For example import_vanco_data

# Import vanco data
import_vanco_data <-  function(sel_time){
  
  vanco_dt  <-  read.csv(here("../Processed_data/Vancomycin.csv"))
  # Plotting stuff
  # Order the time labels
  vanco_dt$Visit <-  vanco_dt$Visit.Name
  # Remove UNSCHED
  vanco_dt <-  vanco_dt[vanco_dt$Visit %in% c("Screening","Day 1","Day 7", "Day 14", "Day 56"),]
  levels_time <- mixedsort(unique(vanco_dt$Visit))
  levels_time <- c("Screening",levels_time[1:(length(levels_time)-1)])
  vanco_dt$Visit <-  factor(vanco_dt$Visit, levels = levels_time)
  
  # Selec samples before antibiotic treatment
  vanco_dt <- vanco_dt %>% filter(collect.pre.post.abx %in% c("Pre.abx"))
  
  # Order the Treatment levels 
  vanco_dt$TRT.norm <-  factor(as.character(vanco_dt$TRT.norm),levels = c("Placebo" ,"VE303 Low Dose","VE303 High Dose"))
  # Take a log10 transform of the abundance data
  vanco_dt$Abundance <- log10(vanco_dt$concentration_corrected)
  
  vanco_dt <-  vanco_dt[!is.na(vanco_dt$concentration_corrected),]
  # Now separate the groups into groups of recurrence / Not -recurrence
  names(vanco_dt)
  unique(vanco_dt$rec.diagnosis)
  table(vanco_dt$Visit)
  library(tidyverse)
  
  ht_dt <-  vanco_dt %>%
    select(MGS.Sample.ID,Subject.Number,Visit,Analyte,Abundance,Abx.group,TRT.norm,rec.diagnosis,Screening_abx_check)
  ht_dt$Log_Abun <-  log10(ht_dt$Abundance)
  
  # Select only given time
  available_visit <-  as.character(unique(ht_dt$Visit))
  if(all(sel_time %in% available_visit)){
    ht_dt <- ht_dt %>% filter(Visit %in% c(sel_time))
    
    # Filter dosed group
    ht_dt <- ht_dt %>% filter(!TRT.norm %in% c("Placebo"))
    
    # Filter Subjects with Antibiotic treatment 
    #  ht_dt <- ht_dt %>% filter(!Screening_abx_check %in% c("Abx Pre-Screen"))
    
    # Aggregate per subject 
    library(dplyr)
    X_mean_bl <-  ht_dt %>%
      pivot_wider(names_from = c(Analyte,Visit),values_from = Log_Abun,id_cols = Subject.Number) %>% data.frame()
    # Clean up data and keep subjects that have both time points
    
    #X_mean_bl <- X_mean_bl[complete.cases(X_mean_bl),]
    return(X_mean_bl)
  }else{
    return(NULL)
  }
  
}
# Import Calprotectin data
import_calprotectin_data <- function(sel_time ){
  cal_dt  <-  read.csv(here("../Processed_data/calprotectin.csv"))
  # Plotting stuff
  # Order the time labels
  cal_dt$Visit <-  cal_dt$Visit.Name
  levels_time <- mixedsort(unique(cal_dt$Visit))
  #levels_time <- c("Screening",levels_time[1:(length(levels_time)-1)])
  cal_dt$Visit <-  factor(cal_dt$Visit, levels = levels_time)
  
  # Selec samples before antibiotic treatment
  cal_dt <- cal_dt %>% filter(collect.pre.post.abx %in% c("Pre.abx"))
  
  
  # Order the Treatment levels 
  cal_dt$TRT.norm <-  factor(as.character(cal_dt$TRT.norm),levels = c("Placebo" ,"VE303 Low Dose","VE303 High Dose"))
  
  # Throw out rows that have NA's as abundance
  #comb_dt <-  comb_dt[!is.na(comb_dt$Abun),]
  # Take a log10 transform of the abundance data
  cal_dt$Abundance <- log10(cal_dt$Result)
  cal_dt <-  cal_dt[!is.na(cal_dt$Abundance),]
  # Now separate the groups into groups of recurrence / Not -recurrence
  names(cal_dt)
  unique(cal_dt$rec.diagnosis)
  table(cal_dt$Visit)
  cal_dt$Analyte <- cal_dt$Lab.Test
  library(tidyverse)
  ht_dt <-  cal_dt %>%
    select(MGS.Sample.ID,Subject.Number,Visit,Analyte,Abundance,Abx.group,TRT.norm,rec.diagnosis,Screening_abx_check)
  ht_dt$Log_Abun <- ht_dt$Abundance
  
  # Select only given time
  available_visit <-  as.character(unique(ht_dt$Visit))
  if(all(sel_time %in% available_visit)){
    ht_dt <- ht_dt %>% filter(Visit %in% c(sel_time))
    
    # Filter dosed group
    ht_dt <- ht_dt %>% filter(!TRT.norm %in% c("Placebo"))
    
    # Filter Subjects with Antibiotic treatment 
    #  ht_dt <- ht_dt %>% filter(!Screening_abx_check %in% c("Abx Pre-Screen"))
    
    # Aggregate per subject 
    library(dplyr)
    X_mean_bl <-  ht_dt %>%
      pivot_wider(names_from = c(Analyte,Visit),values_from = Log_Abun,id_cols = Subject.Number) %>% data.frame()
    # Clean up data and keep subjects that have both time points
    
    #X_mean_bl <- X_mean_bl[complete.cases(X_mean_bl),]
    return(X_mean_bl)
  }else{
    return(NULL)
  }
}
# Import cytokines data
import_cytokines_data <-  function(sel_time){
  # Read Cytokines data
  cyto_dt  <-  read.csv(here("../Processed_data/cytokines.csv"))
  
  table(cyto_dt$Visit.Name)
  # Remove UNSCHED:
  
  cyto_dt <-  cyto_dt[cyto_dt$Visit.Name != "UNSCHED",]
  
  library(gtools)
  # Plotting stuff
  # Order the time labels
  levels_time <- mixedsort(unique(cyto_dt$Visit.Name))
  levels_time <- c("Screening",levels_time[1:(length(levels_time)-1)])
  cyto_dt$Visit <-  factor(cyto_dt$Visit.Name, levels = levels_time)
  
  # Select samples before antibiotic treatment
  cyto_dt <- cyto_dt %>% filter(collect.pre.post.abx %in% c("Pre.abx"))
  
  
  # Order the Treatment levels 
  cyto_dt$TRT.norm <-  factor(as.character(cyto_dt$TRT.norm),levels = c("Placebo" ,"VE303 Low Dose","VE303 High Dose"))
  
  # Now separate the groups into groups of recurrence / Not -recurrence
  names(cyto_dt)
  #cyto_dt$rec.diagnosis <- as.character(cyto_dt$rec.diagnosis.Wk8)
  unique(cyto_dt$rec.diagnosis)
  table(cyto_dt$Visit)
  library(tidyverse)
  ht_dt <-  cyto_dt %>%
    select(Subject.Number,Visit,cytokine,numeric_result_LOQ,Abx.group,TRT.norm,rec.diagnosis,abx_screen_status)
  ht_dt$Log_Abun <-  log10(ht_dt$numeric_result_LOQ )
  
  # Select only screening and Day 1
  
  # Select only given time
  available_visit <-  as.character(unique(ht_dt$Visit))
  if(all(sel_time %in% available_visit)){
    ht_dt <- ht_dt %>% filter(Visit %in% c(sel_time))
    #  ht_dt <- ht_dt %>% filter(Visit %in% c("Screening"))
    # Filter dosed group
    ht_dt <- ht_dt %>% filter(!TRT.norm %in% c("Placebo"))
    # Filter Subjects with Antibiotic treatment 
    #  ht_dt <- ht_dt %>% filter(!Screening_abx_check %in% c("Abx Pre-Screen"))
    
    # Aggregate per subject 
    library(dplyr)
    X_mean_bl <-  ht_dt %>%
      pivot_wider(names_from = c(cytokine,Visit),values_from = Log_Abun,id_cols = Subject.Number) %>% data.frame()
    # Clean up data and keep subjects that have both time points
    
    #  X_mean_bl <- X_mean_bl[complete.cases(X_mean_bl),]
    
    #table(is.infinite(colSums(X_mean_bl)))
    return(X_mean_bl)
  }else{
    return(NULL)
  }
}
# Import metabolomics data
# Import SCFA
import_scfa_data <-  function(sel_time){
  
  # Read SCFA data
  scfa_dt  <-  read.csv(here("../Processed_data/SCFA.csv"))
  # Plotting stuff
  # Order the time labels
  # levels_time <- mixedsort(unique(scfa_dt$Visit))
  # levels_time <- c("Screening",levels_time[1:(length(levels_time)-1)])
  # scfa_dt$Visit <-  factor(scfa_dt$Visit, levels = levels_time)
  
  
  # Select samples before antibiotic treatment
  scfa_dt <- scfa_dt %>% filter(collect.pre.post.abx %in% c("Pre.abx"))
  
  
  # Order the Treatment levels 
  scfa_dt$TRT.norm <-  factor(as.character(scfa_dt$TRT.norm),levels = c("Placebo" ,"VE303 Low Dose","VE303 High Dose"))
  
  # Throw out rows that have NA's as abundance
  #comb_dt <-  comb_dt[!is.na(comb_dt$Abun),]
  # Take a log10 transform of the abundance data
  # Now separate the groups into groups of recurrence / Not -recurrence
  names(scfa_dt)
  scfa_dt$rec.diagnosis <- as.character(scfa_dt$rec.diagnosis)
  unique(scfa_dt$rec.diagnosis)
  table(scfa_dt$Visit) 
  
  ht_dt <-  scfa_dt %>%
    select(MGS.Sample.ID,Sub_ID,Visit,Analyte,Est_Abun,Abx.group,TRT.norm,rec.diagnosis,Screening_abx_check)
  ht_dt$Log_Abun <-  log10(ht_dt$Est_Abun)
  # Select only given time
  available_visit <-  as.character(unique(ht_dt$Visit))
  if(all(sel_time %in% available_visit)){
    ht_dt <- ht_dt %>% filter(Visit %in% c(sel_time))
    
    # Select only screening and Day 1
    #  ht_dt <- ht_dt %>% filter(Visit %in% c("Screening"))
    # Filter dosed group
    ht_dt <- ht_dt %>% filter(!TRT.norm %in% c("Placebo"))
    # Filter Subjects with Antibiotic treatment 
    #ht_dt <- ht_dt %>% filter(!Screening_abx_check %in% c("Abx Pre-Screen"))
    
    # Aggregate per subject 
    library(dplyr)
    X_mean_bl <-  ht_dt %>%
      pivot_wider(names_from = c(Analyte,Visit),values_from = Log_Abun,id_cols = Sub_ID) %>% data.frame()
    # Clean up data and keep subjects that have both time points
    
    #X_mean_bl <- X_mean_bl[complete.cases(X_mean_bl),]
    #table(is.infinite(colSums(X_mean_bl)))
    names(X_mean_bl)[1] <- "Subject.Number"
    
    return(X_mean_bl)
  }else{
    return(NULL)
  }
  
}  
# Import BA
import_ba_data <-  function(sel_time){ 
  # Read BA data
  ba_dt <-   read.csv(here("../Processed_data/BA.csv"))
  # Plotting stuff
  # Order the time labels
  # levels_time <- mixedsort(unique(ba_dt$Visit))
  # levels_time <- c("Screening",levels_time[1:(length(levels_time)-1)])
  # ba_dt$Visit <-  factor(ba_dt$Visit, levels = levels_time)
  # 
  
  # Select samples before antibiotic treatment
  ba_dt <- ba_dt %>% filter(collect.pre.post.abx %in% c("Pre.abx"))
  
  
  # Order the Treatment levels 
  ba_dt$TRT.norm <-  factor(as.character(ba_dt$TRT.norm),levels = c("Placebo" ,"VE303 Low Dose","VE303 High Dose"))
  
  # Take a log10 transform of the abundance data
  ba_dt$Abundance <- log10(ba_dt$Est_Abun)
  
  ba_dt <- ba_dt[!ba_dt$Analyte %in% c("Isodeoxycholic Acid","Glycolithocholic Acid"),]
  table(ba_dt$Visit)
  
  ht_dt <-  ba_dt %>%
    select(Sub_ID,Visit,Analyte,Est_Abun,Abx.group,TRT.norm,rec.diagnosis,Screening_abx_check)
  ht_dt$Log_Abun <-  log10(ht_dt$Est_Abun)
  
  # Select only given time
  available_visit <-  as.character(unique(ht_dt$Visit))
  if(all(sel_time %in% available_visit)){
    ht_dt <- ht_dt %>% filter(Visit %in% c(sel_time))
    
    # Select only screening and Day 1
    #  ht_dt <- ht_dt %>% filter(Visit %in% c("Screening"))
    # Filter dosed group
    ht_dt <- ht_dt %>% filter(!TRT.norm %in% c("Placebo"))
    # Filter Subjects with Antibiotic treatment 
    #ht_dt <- ht_dt %>% filter(!Screening_abx_check %in% c("Abx Pre-Screen"))
    
    # Aggregate per subject 
    library(dplyr)
    X_mean_bl <-  ht_dt %>%
      pivot_wider(names_from = c(Analyte,Visit),values_from = Log_Abun,id_cols = Sub_ID) %>% data.frame()
    # Clean up data and keep subjects that have both time points
    
    #X_mean_bl <- X_mean_bl[complete.cases(X_mean_bl),]
    #table(is.infinite(colSums(X_mean_bl)))
    names(X_mean_bl)[1] <- "Subject.Number"
    
    return(X_mean_bl)
  }else{
    return(NULL)
  }
}
# Import microbiome data:

# Create a list of glommed phyloseq object
phy_glom <-  list()
tax_lev <-  "Class"
phy <- readRDS(here("../Processed_data/phy_mic.rds"))
sm_dt <-  data.frame(sample_data(phy))
# #Only Keep Screening and Baseline samples:
phy_base <-  phy
for(tax_lev in c("Class","Order","Genus")){

  phy_base <-  prune_taxa(taxa_sums(phy_base)>0,phy_base)
  phy_g <- tax_glom(phy_base,tax_lev)
  taxa_names(phy_g) <-  as.character(tax_table(phy_g)[,c(tax_lev)])
  phy_g <- prune_taxa(taxa_sums(phy_g)>0,phy_g)
  phy_glom[[tax_lev]] <- phy_g
}
phy_glom[["Species"]] <- phy_base
saveRDS(phy_glom,paste0(results_folder,"/","phy_glom_list.rds"))

phy_glom <-  readRDS(paste0(results_folder,"/","phy_glom_list.rds"))
#phy_glom[["Species"]] <- phy

import_mic_data<- function(phy_glom,tax_lev,sel_time){
  
  print(paste0("Assembling input data @ ",tax_lev))
  phy_sel <-  phy_glom[[tax_lev]]
  sm_dt <-  data.frame(sample_data(phy_sel))
  phy_melt_sel <-  psmelt(phy_sel)
  # Only HD , LD group
  #  phy_melt_sel <- phy_melt_sel %>% filter(!TRT.norm %in% c("Placebo"))
  #  phy_melt_sel <- phy_melt_sel %>% filter(!Screening_abx_check %in% c("Abx Pre-Screen"))
  
  phy_melt_sel <-  phy_melt_sel %>% filter(Visit.Name %in% c("Screening","Day 1","Day 7","Day 14","Day 28","Day 56"))
  
  # Select time  
  available_visit <-  as.character(unique(phy_melt_sel$Visit.Name))
  if(all(sel_time %in% available_visit)){
    phy_melt_sel <- phy_melt_sel %>% filter(Visit.Name %in% c(sel_time))
    phy_melt_sel <- phy_melt_sel %>% filter(!TRT.norm %in% c("Placebo"))
    phy_melt_sel <- phy_melt_sel %>% filter(collect.pre.post.abx %in% c("Pre.abx"))
    
    
    # Aggregate per subject 
    library(dplyr)
    X_mean_bl <-  phy_melt_sel %>%
      pivot_wider(names_from = c(OTU,Visit.Name),values_from = Abundance,id_cols = Subject.Number) %>% data.frame()
    # Clean up data and keep subjects that have both time points
    
    #  X_mean_bl <- X_mean_bl[complete.cases(X_mean_bl),]
    # rf_mat <-  X_mean_bl  
    # rownames(rf_mat) <-  rf_mat$Subject.Number
    # rf_mat$Subject.Number <-  NULL
    
    # diff_mat <- log2( (rf_mat[,grep("Day",names(rf_mat))] + 0.0001) /( rf_mat[,grep("Screening",names(rf_mat))]+ 0.0001))
    # names(diff_mat) <-  gsub("_Day.1","_diff",names(diff_mat))
    # 
    # rf_mat <- cbind(rf_mat, diff_mat)
    #X_mean_bl <-  data.frame(Subject.Number =  X_mean_bl$Subject.Number,rf_mat)
    return(X_mean_bl)
  }else{
    
    return(NULL)
  }
}


# Run models for each time window
# Loop across each time window
time_models <-  c("Screening","Day 1","Day 7","Day 14","Day 28","Day 56")

# # # Make a data list across all times for each data set
time_data_list <-  list()
#
for(sel_time in time_models){

  vanco_dt <-  import_vanco_data(sel_time )
  cal_dt <-  import_calprotectin_data(sel_time)
  cyto_dt <-  import_cytokines_data(sel_time )
  scfa_dt <- import_scfa_data(sel_time )
  ba_dt <- import_ba_data(sel_time )

  # Import microbiome dataset at different taxa levels:
  mic_class_dt <-  import_mic_data(phy_glom,"Class",sel_time)
  mic_order_dt <-  import_mic_data(phy_glom,"Order",sel_time)
  mic_genus_dt <-  import_mic_data(phy_glom,"Genus",sel_time)
  mic_species_dt <-  import_mic_data(phy_glom,"Species",sel_time)

  # names(mic_species_dt)
  ve303_species <- grep("VE303",names(mic_species_dt),value = T)

  library(genefilter)
  exp_mat <-t(mic_species_dt[,-1])
  rv_sp <- rowVars(exp_mat)
  summary(rv_sp)

  q75_wpn <- quantile( rowVars(exp_mat), .75)  # <= original
  q95_wpn <- quantile( rowVars(exp_mat), .9)  # <= changed to 95 quantile to reduce dataset
  expr_normalized <- exp_mat[ rv_sp > q75_wpn, ]
  rownames(expr_normalized)
  mic_species_dt <- data.frame(Subject.Number = mic_species_dt$Subject.Number, t(expr_normalized))

  data_list <-  list()
  data_list <-  list(vanco_dt,cal_dt,cyto_dt,scfa_dt,ba_dt,mic_class_dt,mic_order_dt,mic_genus_dt,mic_species_dt)
  names(data_list) <- c("Vancomycin","Calprotectin","Cytokines","SCFA","BA","Mic_Class","Mic_Order","Mic_Genus","Mic_Species")


  time_data_list[[sel_time]] <- data_list
}

saveRDS(time_data_list,paste0(results_folder,"/","time_data_list.rds"))

time_data_list <-  readRDS(paste0(results_folder,"/","time_data_list.rds"))
# Import patient metadata
pat_dt_clinical <-  read.csv(here("../Clinical_Tables/cytokine_subject_metadata.csv"))
# Things to add:
# abx_screen_status:
pat_dt_clinical$abx_screen_status
table(pat_dt_clinical$abx_screen_status)

# n_pretreatment events: Number of times a subject was prescribed ABX in a month before screening
pat_dt_clinical$n_pretreatment_events
# total_abx_days_prescreen: Number of days a subjects was exposed to any ABX within a month before screening
pat_dt_clinical$total_abx_days_prescreen
# days_on_study_ABX_pre_screen: Number of days a subjects was exposed to ABX before screening
pat_dt_clinical$days_on_study_ABX_pre_screen
# total_pre_screen_dose: Dose of ABX before screning
pat_dt_clinical$total_pre_screen_dose
# ABX type
pat_dt_clinical$Abx.group
table(pat_dt_clinical$Abx.group)

# table(pat_dt$Subject %in% pat_dt_clinical$Subject.Number)
# pat_dt$Subject %in% pat_dt_clinical$Subject.Number
# pat_dt$Subject[!pat_dt$Subject %in% pat_dt_clinical$Subject.Number]

# Filter importance variables:
# Select Subject, Age, BMI and number of previous recurrences
pat_dt_clinical_sel  <- pat_dt_clinical %>%
  select(Subject.Number,Age,BMI,total.recurrences,abx_screen_status,n_pretreatment_events,
         total_abx_days_prescreen,n_washout_days,Abx.group )

str(pat_dt_clinical_sel)  
pat_dt_clinical_sel$abx_screen_status <-  factor(pat_dt_clinical_sel$abx_screen_status)
pat_dt_clinical_sel$Abx.group <-  ifelse(pat_dt_clinical_sel$Abx.group == "Vancomycin","Vancomycin","Others")
pat_dt_clinical_sel$Abx.group <- factor(pat_dt_clinical_sel$Abx.group)
#  pat_dt[,c(1,4,11,12)]
str(pat_dt_clinical_sel)

# Add treat day len abx:
trt_abx_len_dt <-  sm_dt %>% select(Subject.Number,treat.day.len.abx)%>% unique()

pat_dt_clinical_sel <-  pat_dt_clinical_sel %>%  
  left_join(trt_abx_len_dt)

#names(pat_dt) <-  c("Subject.Number","Age","BMI","Num_prev_recurrence")


# Recurrence data
rec_dt <-  unique(sm_dt[,c("Subject.Number","rec.diagnosis.Wk8")])
names(rec_dt)[2] <- "Y_class"
rec_dt$Y_class <- factor(rec_dt$Y_class )
levels(rec_dt$Y_class) <-c("R","NR")


# Import the colonization data:
import_ve303_abun <- function(){
  # Import total VE303 abundance per subject 
  # Import marker panel data:
  m_panel_dt <-  read.csv(here("../Processed_data/VE303_abun.csv"))
  m_panel_dt$Visit <- as.character(m_panel_dt$Visit.Name)
  library(gtools)
  sorted_time <- mixedsort(unique(m_panel_dt$Visit))
  # Overall VE303 colonization:
  # Sort time points in increasing order of pre and post
  pre_labels <-  c("Screening",sorted_time[1:which(sorted_time == "Day 1")])
  # Limit the time to Day 56
  post_labels <- sorted_time[(which(sorted_time == "Day 1")+1):which(sorted_time == "Day 168")]
  
  # Subset only Screening ,Day 1, Day 7 and Day 14
  pre_labels <-  c("Screening","Day 1")
  post_labels <- c("Day 7","Day 14","Day 28","Day 56")
  
  # Subset only Screening ,Day 1, Day 7 and Day 14
  panel_dt <- m_panel_dt[m_panel_dt$Visit %in% c(pre_labels,post_labels),]
  unique(panel_dt$Visit)
  
  # Order the time labels
  levels_time <- mixedsort(unique(panel_dt$Visit))
  levels_time <- c("Screening",levels_time[1:(length(levels_time)-1)])
  panel_dt$Visit <-  factor(panel_dt$Visit, levels = levels_time)
  
  # Order the Treatment levels 
  panel_dt$TRT.norm <-  factor(as.character(panel_dt$TRT.norm),levels = c("Placebo" ,"VE303 Low Dose","VE303 High Dose"))
  panel_dt$N_det <-  ifelse(panel_dt$detection_status == "Detected", 1,0)
  #tmp <-  panel_dt[panel_dt$detection_status == "Not detected",]
  # Lets compute total number of detected strains in each sample
  num_dt <- panel_dt
  num_dt$est_rel_abundance_adj <-  num_dt$targeted_panel_relative_abundance
  num_dt$est_rel_abundance_adj[num_dt$N_det == 0] <- 0
  
  num_dt$normalized_marker_depth_adj <-  num_dt$normalized_marker_depth
  num_dt$normalized_marker_depth_adj[num_dt$N_det  == 0] <- 0
  
  # Remove pre-Abx treated group
  library(dplyr)
  #num_dt <- num_dt %>% filter(!Screening_abx_check %in% c("Abx Pre-Screen"))
  library(tidyverse)
  
  sum_num_dt <- num_dt %>% 
    group_by(Subject.Number,TRT.norm,filename,Visit) %>% 
    summarise(N_sum_det = sum(N_det)) %>%
    ungroup()
  
  
  num_strain_list <- list()
  
  total_visit <-  c(pre_labels,post_labels)
  time <- total_visit[4]
  for(time in total_visit){
    
    sel_num_dt <- sum_num_dt %>% 
      filter(Visit %in% time) %>%
      ungroup()
    
    sel_num_dt <-  sel_num_dt[!sel_num_dt$TRT.norm %in% c("Placebo"),]
    
    
    sel_num_dt <-  sel_num_dt[,c("Subject.Number", "TRT.norm","N_sum_det")]
    
    num_strain_list[[time]] <-  sel_num_dt
    
  }
  #  saveRDS(num_strain_list,paste0(results_folder,"/det_list.rds"))
  return(num_strain_list)
  
}

# Set factors for TRT.norm, Screening
#rec_dt$TRT.norm <- factor(ifelse(rec_dt$TRT.norm == "Placebo","Placebo","Dosed"))
#rec_dt$Screening_abx_check <- factor(rec_dt$Screening_abx_check)


rec_dt <-  rec_dt %>%
  left_join(pat_dt_clinical_sel)
# %>%
#  filter(!TRT.norm %in% c("Placebo"))

library(doParallel)
detectCores()
# myCluster <- makeCluster(20)
# registerDoParallel(myCluster)


library(randomForest)
library(Boruta)

# Number of undersampling
N_samp <- 100
# Number of boruta runs each undersampled data
N_boruta <- 5

library(caret)

# mcLapply version of boruta
boruta_par <-  function(dt_bor,n_bor){
  library(Boruta)
  # Set seed for Boruta
  set.seed(300 + n_bor)
  boruta_model <- Boruta(Y_class ~., data=dt_bor,
                         doTrace=2,
                         mcAdj = T, maxRuns = 200)
  boruta_signif <- names(boruta_model$finalDecision[boruta_model$finalDecision %in% c("Confirmed","Tentative")])  # collect Confirmed and Tentative variables
  
  if(length(boruta_signif) == 0){
    var_stat <-  attStats(boruta_model)
    var_stat <-  var_stat %>% arrange(desc(meanImp))
    # If there are no species
    # Get top 10 taxa interms of meanImp
    
    if(ntax < 10){
      boruta_signif <- rownames(var_stat)[1:ntax]
    }else {
      boruta_signif <- rownames(var_stat)[1:10]
    }
  }
  
  return(boruta_signif)
  
}

################################################################################

# Create undersampled datasets for LOOCV
# Loop through each sampled time:
time_slots <-  names(time_data_list)
sel_time <-  time_slots[1]

#test_data_list <- list()
#train_data_list <- list()
under_sampled_data_list <- list()
boruta_list <-  list()
pred_time_list <- list()
for(sel_time in time_slots){
  print(paste0("Evaluation at time: ",sel_time))
  data_list <- time_data_list[[sel_time]]
  # Remove NULL elements
  data_list <-  purrr::compact(data_list)
  # Now loop over each data set:
  dt_set <-  names(data_list)[9]
  
  pred_set <-  list()
  for(dt_set in names(data_list)){
    
    print(paste0("Evaluation using dataset: ",dt_set))
    rf_mat <-  data_list[[dt_set]]
    #Import strain list
    num_st_dt <- import_ve303_abun()[[sel_time]]
    # Join the number of strain detected
    mic_datasets <-grep("Mic",names(data_list),value = T) 
    if(dt_set %in% mic_datasets){
      rf_mat <- rf_mat %>%
        left_join(num_st_dt %>% select(Subject.Number,N_sum_det))
    }
    rf_mat <- rf_mat %>%
      left_join(rec_dt[,c(names(rec_dt)[-2],"Y_class")])
    rf_mat <- rf_mat[complete.cases(rf_mat),]
    rownames(rf_mat) <-  rf_mat$Subject.Number
    rf_mat$Subject.Number <-  NULL
    
    t_data <-  rf_mat
    ntax <- ncol(t_data) -1
    num_R <-  as.numeric(table(t_data$Y_class)[1])
    # Only use datasets that have atleast n >= 5 instances
    # and recurring subjects > 3
    
    # Save prediction in each loo cases
    prediction_list <- list()
    if(nrow(t_data) >= 5 & num_R >4) {
      
      # Create stratified Kfolds
      # Set the number of folds for cross validation
      num_folds <- 5
      
      # Create the folds using stratified sampling
      set.seed(1057)
      folds <- createFolds(t_data$Y_class, k = num_folds, list = TRUE, returnTrain = FALSE)
      
      # Create undersampled data 
      for(loo in seq_along(folds)){
        print(paste0("StratKFold :" ,loo))
        # Get the training and testing data for this fold
        train_set_dt <- t_data[-folds[[loo]], ]
        test_set_dt <- t_data[folds[[loo]], ]
        
        # Create undersampled data for training set:
        # Save undersampled datalist for training
        #print(paste0("Undersampling loo: ",loo ))
        
        # Save predictions using each sampled dataset
        prediction_sampled <-  c()
        for(samp_n in 1:N_samp){
          
          pred_test <-  c()
          seed_sample <- 300 + samp_n
          #set.seed(seed_sample)
          # under_sampled_data <- ovun.sample(Y_class~., data = train_set_dt,
          #                                   p=0.5, seed=seed_sample, method="under")$data
          x_mat <- train_set_dt[,1:(ncol(train_set_dt) -1)]
          y_mat <- train_set_dt$Y_class
          # Use downsampling
          set.seed(seed_sample)
          under_sampled_data <- downSample(x = x_mat, y = y_mat, yname = "Y_class")
          
          under_sampled_data_list[[sel_time]][[dt_set]][[as.character(samp_n)]] <-   under_sampled_data
          
          #print(paste0("Running boruta: ",samp_n))
          
          # Run boruta parallely
          sig_sp_list <- mclapply(1:N_boruta,boruta_par,dt_bor = under_sampled_data,mc.cores = 15)
          imp_vars <-  unique(unlist(sig_sp_list))
          boruta_list[[sel_time]][[dt_set]][[as.character(samp_n)]] <- imp_vars
          
          
          # Run random forest model
          if(length(imp_vars) > 0){
            
            set.seed(300 + loo + samp_n*1000)
            rf_model_fmt <- randomForest(Y_class~.,data=under_sampled_data[,c(imp_vars,"Y_class")],
                                         importance = T,
                                         ntree = 5000)
            
            ntax <-  length(imp_vars)
            #print(paste0("Predicting iteration: ",samp_n))
            # Now perform predictions:
            if(ntax == 1){
              pred_test <-   predict(rf_model_fmt,
                                              newdata=test_set_dt[,1:(ncol(test_set_dt)-1),drop = F],
                                              type = "prob")[,1]
            }else {
              pred_test <-   predict(rf_model_fmt,
                                              newdata=test_set_dt[,1:(ncol(test_set_dt)-1)],
                                              type = "prob")[,1]
              }
            
          }else{
            
            pred_test <- NA
          }
          
          prediction_sampled[[as.character(samp_n)]] <-  as.numeric(pred_test)
          
          
        }# End of the undersampling loop
        
        prediction_list[[as.character(loo)]] <- prediction_sampled
        
      }# End of the loo loop 
      
      pred_loo_list <-  list()
      # Convert list into dataframe:
      for(len_list in 1:length(prediction_list)){
        
        pred_samp_list <- prediction_list[[len_list]]
        for(sam_length in 1:length(pred_samp_list)){
          dt_pred <- data.frame(fold_num = len_list,fold_id = folds[[len_list]],samp_n =  sam_length,pred = prediction_list[[len_list]][[sam_length]]) 
          pred_loo_list[[paste0(len_list,"_",sam_length)]] <- dt_pred
        }
      }
      dt_sample_pred <- do.call("rbind",pred_loo_list)
      obs_dt  <-  data.frame(fold_id = 1:nrow(t_data),obs = t_data$Y_class)
      obs_dt$dataset <- dt_set
      obs_dt$Time <- sel_time
      obs_dt <- dt_sample_pred %>% left_join(obs_dt)
      
      pred_set[[dt_set]] <- obs_dt 
      
    }# if case only passes through when the data has enough instances of recurrence and rows 
    
    
    
  }# Loop over datasets 
  
  pred_time_list[[sel_time]] <- do.call("rbind",pred_set)
  
}# Loop over time

# Save pred list
saveRDS(pred_time_list,paste0(results_folder,"/","pred_time_set_list_all.rds"))

# Save under sampled data and boruta list:
saveRDS(boruta_list,paste0(results_folder,"/boruta_undersampled_list.rds"))  
saveRDS(under_sampled_data_list,paste0(results_folder,"/under_sampled_data_list.rds"))  


