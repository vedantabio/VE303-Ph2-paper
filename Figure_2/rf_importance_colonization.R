# Create a directories to store results from the analysis
library(here)
mainDir <- here("colonization_model_importance")
dir.create(file.path(mainDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,sep="/")

# Performance:
library(randomForest)
library(Boruta)
library(parallel)

# Import prediction results
pred_sp_list <-  readRDS(here("colonization_model_performance/pred_sp_list.rds"))

final_pred_dt <-  do.call("rbind",pred_sp_list)
f_pred_dt <-  final_pred_dt 
names(f_pred_dt)[4] <- "pred_prob"
f_pred_dt$obs <- ifelse(f_pred_dt$obs == "C",1,0)



compute_metrics <- function(pred_dt){
  # Returns a list of variables
  # F1 score
  # Precision
  # Recall
  # AUC
  # ROC dataframe 
  # Create a list 
  # names(pred_dt)[1] <- "pred_prob"
  # pred_dt$obs <- ifelse(pred_dt$obs == "Yes",1,0)
  # 
  print(pred_dt$dataset[1])
  print(pred_dt$Time[1])
  print(table(pred_dt$obs))
  metric <- list()
  
  # ROC curve
  library(caret)
  #create confusion matrix and calculate metrics related to confusion matrix
  cm <-  confusionMatrix(factor(ifelse(pred_dt$pred_prob > 0.5,1,0)), factor(pred_dt$obs), mode = "everything", positive="1")
  
  library(pROC)
  roc_mod  <- roc(pred_dt$obs, pred_dt$pred_prob,smooth = F, plot=TRUE,legacy.axes=TRUE,
                  xlab="False Positive Percentage",ylab="True Postive Percentage",col="#377eb8", lwd=4,
                  direction = "<",
                  print.auc=TRUE, print.thres=TRUE)
  
  metric[["AUC"]] <- as.numeric(gsub("*: ",'',roc_mod$auc))
  
  # library(MLmetrics)
  # auc <- MLmetrics::AUC(pred_prob, ifelse(t_data$Y_class == "R",1,0))
  roc_dt <-  data.frame(TPP = roc_mod$sensitivities, FPP = 1- roc_mod$specificities,
                        AUC = as.numeric(gsub("*: ",'',roc_mod$auc)))
  metric[["ROC"]] <-  roc_dt
  
  
  # Precision 
  precision <-  cm[["byClass"]][["Precision"]]
  metric[["precision"]] <- precision
  # Precision 
  recall <-  cm[["byClass"]][["Recall"]]
  metric[["recall"]] <- recall
  # F1 score
  f1_score <- cm[["byClass"]][["F1"]] 
  metric[["f1_score"]] <- f1_score
  
  return(metric)
  
  # pred_prob <- prediction_fmt_list[[as.character(seed)]]  
  # obs <-  ifelse(t_data$Y_class == "R",1,0)
  # table(obs)
  # 
}

ve303_species <- names(pred_sp_list)
data_names <-  unique(f_pred_dt$dataset)

metric_sp_list <-  list()

N_samp <-  length(unique(f_pred_dt$samp_n))
library(tidyverse)
#options(warn = 0)    

freq_list <- list()

ve_sp <- ve303_species[1]
library(tidyverse)
for(ve_sp in ve303_species){
  
  dt_pred_dt <-  f_pred_dt %>% filter(Species == ve_sp)
  data_names <- unique(f_pred_dt$dataset)
  for(dt_set in data_names){
    
    for(sampN in 1:N_samp){
      pred_dt <- dt_pred_dt %>%
        filter(dataset == dt_set) %>%
        filter(samp_n == sampN) 
      
      
      
      num_C <-  as.numeric(table(pred_dt$obs)[1])
      num_NC <- as.numeric(table(pred_dt$obs)[2])
      num_Case <- min(num_C,num_NC)
      
      # Only use datasets that have atleast n >= 5 instances
      # and recurring subjects > 3
      if(nrow(pred_dt) >= 5 & num_Case >4) {
        
        m_list <-  compute_metrics(pred_dt)
        auc <-  m_list$AUC
        f1_score <- m_list$f1_score
        precision <-  m_list$precision
        recall <-  m_list$recall
        
        metric_dt <-  data.frame(auc ,f1_score,
                                 #precision,
                                 #recall ,
                                 samp_n = sampN,
                                 Species = ve_sp,
                                 dataset = dt_set)
        
        metric_sp_list[[paste0(ve_sp,dt_set,sampN)]] <- metric_dt
      }
    } 
    pred_freq_dt <- dt_pred_dt %>%
      filter(dataset == dt_set) %>%
      filter(samp_n == 1) 
    # Compute frequency
    freq_dt <-  data.frame(t(as.matrix(table(pred_dt$obs))))
    freq_dt$dataset <-  dt_set
    freq_dt$Species  <-  ve_sp
    
    freq_list[[paste0(ve_sp,dt_set)]] <-  freq_dt
    
    
    
  }
}  

# Frequency:
freq_dt <-  do.call("rbind",freq_list)

freq_dt$N  <-  freq_dt$X0 + freq_dt$X1
freq_dt$N <- paste0(freq_dt$N,"(", freq_dt$X1,"/", freq_dt$X0,")")

freq_dt_w <- freq_dt %>%
  select(dataset,Species,N)%>%
  pivot_wider(names_from = dataset,values_from = N)

# Subset and order dataset

freq_dt_w <- freq_dt_w %>%
  select( Vancomycin,Cytokines,SCFA,BA,Mic_Species)
names(freq_dt_w)[c(5)] <- c("Microbiome")


m_dt <-  do.call("rbind", metric_sp_list)

metric_dt_long <- m_dt %>%
  pivot_longer(cols = c(auc,f1_score
                        #precision,recall
  ),names_to = "Metric",values_to = "value")

head(metric_dt)
metric_dt_long <- na.omit(metric_dt_long)

# Barplots with se for each model and taxa:
# Calculates mean, sd, se and IC
sum_dt <- metric_dt_long %>%
  dplyr::group_by(Species,dataset,Metric) %>%
  summarise( 
    n=n(),
    mean=mean(value),
    sd=sd(value)
  ) %>%
  mutate( se=sd/sqrt(n))  

sum_dt$dataset <- factor(sum_dt$dataset , levels = c("Mic_Class","Mic_Order","Mic_Genus",
                                                     "Mic_Species",
                                                     "Vancomycin","SCFA","BA","Cytokines"))
library(hues)
set.seed(100)
col_pred <-  iwanthue(length(unique(sum_dt$dataset)), plot=TRUE)
names(col_pred) <- levels(sum_dt$dataset)
#sum_dt$tax_level <- factor(sum_dt$tax_level , levels = c("Class","Order","Genus","Species"))
# Standard Error

library(wesanderson)

sum_dt <- sum_dt[sum_dt$Metric != "f1_score",]
facet_lab <-  c("auc" = "AUC", "f1_score" = "F1-score")

perf_dt <-  sum_dt %>%
  left_join(freq_dt %>% select(dataset,Species, N))
perf_dt$dataset <- factor(perf_dt$dataset , levels = c("Mic_Class","Mic_Order","Mic_Genus",
                                                       "Mic_Species","Vancomycin","SCFA","BA","Cytokines"))

library(lemon)
# bar_p <- ggplot(perf_dt,aes(x = dataset,y = mean ,fill = dataset)) +
#   geom_bar( position = position_dodge(preserve = "single"), stat="identity",color = "black", alpha=0.8) +
#   geom_errorbar( aes(x=dataset, ymin=mean-se, ymax=mean+se), width=.2,
#                  position=position_dodge(width = .9,preserve = "single"), colour="black", alpha=0.9, size=0.5) +
#   # facet_rep_grid(Time ~ Metric,repeat.tick.labels = T )+
#   geom_text(aes(label=N,group = dataset, y= 0.95),position=position_dodge(width = 1),
#             hjust = 'left',angle = 90,size=3) +
#   facet_rep_grid( ~ Species,repeat.tick.labels = T )+
#   geom_hline(aes(yintercept = 0.5),color = "red")+
#   scale_fill_manual(name = "Dataset",values = col_pred)+
#   scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.1),
#                      labels = c("0","0.2","0.4","0.6","0.8","1",""),
#                      limits = c(0,1.1))+
#   ylab("AUC")+
#   xlab("Dataset")+
#   theme_bw()+
#   theme(axis.text.x=element_text(size=12,face="bold",angle = 60,hjust=1),
#         #text.y=element_text(size=10,face="bold"),
#         axis.title.y = element_text(size=10,face="bold"),
#         axis.title.x = element_text(size=10,face="bold"),
#         strip.text.x = element_text(size=10,face="bold"),
#         legend.title=element_text(size=10), 
#         legend.text=element_text(size=10),
#         axis.text.y = element_text(size=12,face="bold"))+
#   ggtitle("Model performance")
# 
# pdf(paste0(results_folder,"/RF_Model_performance.pdf"),width = 13,height = 6)
# print(bar_p)
# dev.off()



library(lemon)
bar_p <- ggplot(perf_dt,aes(x = dataset,y = mean ,fill = dataset)) +
  
  geom_hline(aes(yintercept = 0.5),color = "red",linetype = "longdash",linewidth = 1)+
  geom_linerange(aes(x=dataset,color = dataset,ymin = mean-sd, ymax = mean+sd),
                 alpha = 0.8,linewidth = 1.5)+
  geom_point(alpha=0.8,shape = 23,stroke = 1,color = "black",size = 3)+
  geom_text(aes(label=N,group = dataset, y= 0.95),position=position_dodge(width = 1),
            hjust = 'left',angle = 90,size=3) +
  facet_rep_grid( ~ Species,repeat.tick.labels = T )+
  scale_fill_manual(name = "Dataset",values = col_pred)+
  scale_color_manual(name = "Dataset",values = col_pred)+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1,1.1),
                     labels = c("0","0.2","0.4","0.6","0.8","1",""),
                     limits = c(0,1.1))+
  ylab("AUC")+
  xlab("Dataset")+
  theme_bw()+
  theme(axis.text.x=element_text(size=12,face="bold",angle = 60,hjust=1),
        #text.y=element_text(size=10,face="bold"),
        axis.title.y = element_text(size=10,face="bold"),
        axis.title.x = element_text(size=10,face="bold"),
        strip.text.x = element_text(size=10,face="bold"),
        legend.title=element_text(size=10), 
        legend.text=element_text(size=10),
        axis.text.y = element_text(size=12,face="bold"))+
  ggtitle("Model performance")

pdf(paste0(results_folder,"/RF_Colonization_Model_performance.pdf"),width = 16,height = 6)
print(bar_p)
dev.off()


write.csv(perf_dt,paste0(results_folder,"/Performance_colonization.csv")) 

library(gtools)
library(tidyr)
library(phyloseq)
library(vita)
library(Boruta)
library(randomForest)
library(phyloseq)
library(tidyverse)
library(lemon)


# Import functions to read data from various sources:

time_data_list <-  readRDS(here("colonization_model_performance/time_data_list.rds"))
# Import patient metadata
#pat_dt <-  read.csv("../../data/Clinical_Tables/2022_10_22_VE303_002_demographics.csv")

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

#names(pat_dt) <-  c("Subject.Number","Age","BMI","Num_prev_recurrence")
phy <- readRDS(here("../Processed_data/phy_mic.rds"))
sm_dt <-  data.frame(sample_data(phy))
# Add treat day len abx:
trt_abx_len_dt <-  sm_dt %>% select(Subject.Number,treat.day.len.abx)%>% unique()

pat_dt_clinical_sel <-  pat_dt_clinical_sel %>%  
  left_join(trt_abx_len_dt)

# VE303 Colonization data:
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
  post_labels <- c("Day 7","Day 14","Day 28")
  
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
  det_list <-  list()
  ve303_list <-  list()
  org <-  unique(num_dt$organism)[1]
  for(org in unique(num_dt$organism)){
    
    sum_num_dt <- num_dt %>% 
      filter(organism %in% org) %>%
      group_by(Subject.Number,TRT.norm,filename,Visit) %>% 
      summarise(N_sum_det = sum(N_det),
                sum_abun = sum(est_rel_abundance_adj)) %>%
      ungroup()
    
    # Select either Day 7 or Day 14
    
    sum_num_dt <-  sum_num_dt[sum_num_dt$Visit %in% c("Day 7","Day 14","Day 28"),]
    sum_num_dt <-  sum_num_dt[!sum_num_dt$TRT.norm %in% c("Placebo"),]
    
    # Pivot Wider:
    # sum_num_dt_w <-  sum_num_dt %>%
    #   dplyr::select(Subject.Number,TRT.norm,Visit,sum_abun)%>%
    #   pivot_wider(names_from = Visit,values_from =sum_abun) %>%
    #   data.frame()
    
    # Detection
    sum_num_dt_w <-  sum_num_dt %>%
      dplyr::select(Subject.Number,TRT.norm,Visit,N_sum_det)%>%
      pivot_wider(names_from = Visit,values_from =N_sum_det) %>%
      data.frame() %>%
      rowwise()%>%
      mutate(sum_det = sum(c_across(Day.7:Day.28), na.rm = TRUE)) 
    
    sum_num_dt_w <-  sum_num_dt_w[,c("Subject.Number", "TRT.norm","Day.7","Day.14","Day.28", "sum_det" )]
    
    det_list[[org]] <-  sum_num_dt_w
    # sum_num_dt_w$Day.14[is.na(sum_num_dt_w$Day.14)] <-  sum_num_dt_w$Day.7[is.na(sum_num_dt_w$Day.14)]
    sum_num_dt_w <- sum_num_dt_w[,c("Subject.Number","sum_det")]
    names(sum_num_dt_w)[2] <- "Y_class"
    # Convert abundances to presence and absence 
    sum_num_dt_w$Y_class <- ifelse(sum_num_dt_w$Y_class > 1,"C","NC")
    sum_num_dt_w$Y_class <-  factor(sum_num_dt_w$Y_class)
    ve303_list[[org]] <- sum_num_dt_w
    
  }
  saveRDS(det_list,paste0(results_folder,"/det_list.rds"))
  return(ve303_list)
  
}


# Set factors for  Screening
#rec_dt$Screening_abx_check <- factor(rec_dt$Screening_abx_check)

# %>%
#  filter(!TRT.norm %in% c("Placebo"))

library(doParallel)
detectCores()
library(randomForest)
library(Boruta)

# Number of boruta runs each undersampled data
N_boruta <- 10

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




# Save boruta list for LOOCV for each seed
library(randomForest)
library(Boruta)
# Lets use Boruta for feature selection:
# Boruta to select important species.
# Implementation of Boruta 
# First step to filter out species that relates to the genotype
# Loop over multi seed


ve303_col_list <-  import_ve303_abun()
ve303_sp <-  names(ve303_col_list)[c(1:5,7:8)]
# ve303_sp <-  names(ve303_col_list)

# Loop through each dataset:
data_list <- time_data_list[["Screening"]]
# Remove NULL elements
data_list <-  purrr::compact(data_list)
names(data_list)


# Save predictions for each species
imp_dt_list <- list()
#Loop over species
ve_sp <-  ve303_sp[1]
for( ve_sp in ve303_sp){
  
  # Looping over VE303 species
  print(paste0("Looping over: ",ve_sp))
  # Now loop over each data set:
  # Save prediction for each dataset
  for(dt_set in names(data_list)){
    
    print(paste0("Looping over: ",dt_set))
    rf_mat <-  data_list[[dt_set]]
    
    colo_dt <- ve303_col_list[[ve_sp]] 
    colo_dt <-  colo_dt %>%
      left_join(pat_dt_clinical_sel)
    
    rf_mat <- rf_mat %>%
      left_join(colo_dt[,c(names(colo_dt)[-2],"Y_class")])
    
    rf_mat <- rf_mat[complete.cases(rf_mat),]
    rownames(rf_mat) <-  rf_mat$Subject.Number
    rf_mat$Subject.Number <-  NULL
    
    t_data <-  rf_mat
    ntax <- ncol(t_data) -1
    
    num_C <-  as.numeric(table(t_data$Y_class)[1])
    num_NC <- as.numeric(table(t_data$Y_class)[2])
    num_Case <- min(num_C,num_NC)
    
    imp_list <- list()
    
    # Only use datasets that have atleast n >= 5 instances
    # and recurring subjects > 3
    # Save prediction in each loo cases
    if(nrow(t_data) >= 5 & num_Case >4) {
      
      train_set_dt <- t_data
      # Run boruta parallely
      N_boruta <-  10
      # Run boruta
      sig_sp_list <- mclapply(1:N_boruta,boruta_par,dt_bor = train_set_dt,mc.cores = 15)
      
      for(j in 1:length(sig_sp_list)){
        
        imp_vars <- sig_sp_list[[j]]
        
        if(length(imp_vars) > 0){
          
          set.seed(300  + j*1000)
          # Run random forest
          rf_model_fmt <- randomForest(Y_class~.,data=train_set_dt[,c(imp_vars,"Y_class")],
                                       importance = T,
                                       ntree = 5000)
          # Extract variable importance using Boruta picked variables:
          library(vita)
          set.seed(340)
          pimp <- PIMP(train_set_dt[,c(imp_vars),drop = F], train_set_dt$Y_class,rf_model_fmt,parallel=TRUE, ncores=10 ,seed = 340)
          
          p_test <- NULL
          pimp_all <-  NULL
          if(length(imp_vars) == 1){
            
            # Doesn't work for only one variable 
            # Empirical cumulative null distribution
            Fn0 =  ecdf(pimp$PerVarImp)
            p.val <-  1  -  Fn0(pimp$VarImp)
            pimp_all <- data.frame(VarImp =  pimp$VarImp, p.value = p.val)
            
          }else{
            p_test <- PimpTest(pimp)
            pimp_all <- data.frame(orig =  p_test$VarImp, p_test$pvalue)
          }
          
          imp_dt <- pimp_all
          imp_dt$variable <- rownames(imp_dt)
          imp_dt$dataset <-  dt_set
          imp_dt$Time <- "Screening"
          imp_dt$Species <- ve_sp
          imp_dt$rep <-  j
          imp_list[[j]] <-  imp_dt
        }
        
      }
      
      # if statement ends here for enough instances of colonization
      # Variable importance
      # Create a summary of variable importance alongside frequency
      imp_dt_sum <-  do.call("rbind",imp_list)
      # Filter by pval
      imp_dt_sum <-  imp_dt_sum %>% 
        filter(VarImp > 0) %>%
        filter(p.value <= 0.1)
      
      imp_dt_list[[ve_sp]][[dt_set]] <- imp_dt_sum
      
    }else {
      imp_dt_list[[ve_sp]][[dt_set]] <- NA
    } 
    
  } # Loop over datasets
  
  
} #  End Looping over VE303 species   


saveRDS(imp_dt_list,paste0(results_folder,"/imp_dt_list.rds"))  



imp_dt_list <-  readRDS(paste0(results_folder,"/imp_dt_list.rds"))


# Species naming 
sp_name_dt <-  readxl::read_excel(here("../Processed_data/Feb 2021 Vedanta database.xlsx"))%>%
  select(GTDB_Species,NCBI_classification,StrainID,GTDB_closest_placement_taxonomy)
sp_name_dt$GTDB_closest_species <- gsub(".*s__","",sp_name_dt$GTDB_closest_placement_taxonomy)
sp_name_dt$NCBI_species <- gsub(".*s__","",sp_name_dt$NCBI_classification)

sp_name_dt <- sp_name_dt[grep("MP|HF",sp_name_dt$StrainID),]


imp_comp_list <-  list()
for(ve_sp in names(imp_dt_list)){
  sp_list <-  purrr::compact(imp_dt_list[[ve_sp]])
  imp_comp_list[[ve_sp]] <-  sp_list[["Mic_Species"]]
}
imp_comp_dt <-  do.call("rbind",imp_comp_list)
imp_comp_dt <-  imp_comp_dt[grep("MP|HF",imp_comp_dt$variable),]

mp_names <-  gsub("_Screening","",unique(imp_comp_dt$variable))
mp_dt <-  data.frame(variable= unique(unique(imp_comp_dt$variable)),MP_ID = mp_names)
mp_dt$MP_ID <-  gsub(".*\\.\\.MP","",mp_dt$MP_ID)
mp_dt$MP_ID <-  gsub(".*\\.\\.HF","",mp_dt$MP_ID)
mp_dt$MP_ID <-  gsub("_.*","",mp_dt$MP_ID)
mp_dt$MP_ID[grep("MP",mp_dt$variable)] <-  paste0("MP",mp_dt$MP_ID[grep("MP",mp_dt$variable)])
mp_dt$MP_ID[grep("HF",mp_dt$variable)] <-  paste0("HF",mp_dt$MP_ID[grep("HF",mp_dt$variable)])

length(mp_names)

gtdb_sp_name <-  c()
ncbi_sp_name <- c()
for(mp_id in mp_dt$MP_ID){
  gtdb_sp <- unique( sp_name_dt[grep(mp_id,sp_name_dt$StrainID),]$GTDB_closest_species)
  ncbi_sp <- unique( sp_name_dt[grep(mp_id,sp_name_dt$StrainID),]$NCBI_species)
  gtdb_sp_name <-  c(gtdb_sp_name,gtdb_sp)
  ncbi_sp_name <-  c(ncbi_sp_name,ncbi_sp)
}
mp_dt$gtdb_sp_name <- gtdb_sp_name
mp_dt$ncbi_sp_name <- ncbi_sp_name
mp_dt$ncbi_sp_name <- gsub(";.*","",mp_dt$ncbi_sp_name)
#mp_dt$ncbi_sp_name[is.na(mp_dt$ncbi_sp_name)] <- "Bacteroides sp. OF02-3LB"

mp_dt$ncbi_sp_name<- paste0(mp_dt$ncbi_sp_name, "_",mp_dt$MP_ID)
# For duplicate names, add MP id at the end:
# idx_dup <- duplicated(mp_dt$ncbi_sp_name) | duplicated(mp_dt$ncbi_sp_name, fromLast = TRUE)
# mp_dt$ncbi_sp_name[idx_dup] <- paste0(mp_dt$ncbi_sp_name[idx_dup], "_",mp_dt$MP_ID[idx_dup])


# Plot the importance plots and barplot

# Import prediction results
time_slots <- c("Screening")

p_imp_list  <- list()
p_abun_list <- list()

N_samp <-  10
fil_imp_dt_list <-  list()
heatmap_imp_dt_list <- list()
for(ve_sp in ve303_sp){ 
  print(paste0("Running for ",ve_sp))
  # Now loop over each data set:
  dt_set <-  names(data_list)[9]
  
  # pdf(paste0(results_folder,"/",sel_time,"_Variable_imp_RF.pdf"), width = 25, height = 25)
  for(dt_set in names(data_list)){
    # Train data
    print(paste0("Looping over: ",dt_set))
    rf_mat <-  data_list[[dt_set]]
    
    colo_dt <- ve303_col_list[[ve_sp]] 
    colo_dt <-  colo_dt %>%
      left_join(pat_dt_clinical_sel)
    
    rf_mat <- rf_mat %>%
      left_join(colo_dt[,c(names(colo_dt)[-2],"Y_class")])
    
    rf_mat <- rf_mat[complete.cases(rf_mat),]
    rownames(rf_mat) <-  rf_mat$Subject.Number
    rf_mat$Subject.Number <-  NULL
    
    t_data <-  rf_mat
    ntax <- ncol(t_data) -1
    
    num_C <-  as.numeric(table(t_data$Y_class)[1])
    num_NC <- as.numeric(table(t_data$Y_class)[2])
    num_Case <- min(num_C,num_NC)
    
    if(nrow(t_data) >= 5 & num_Case >4) {
      
      train_set_dt <- t_data
      
      # Variable importance
      imp_dt_fil <-  imp_dt_list[[ve_sp]][[dt_set]] 
      
      if(nrow(imp_dt_fil)> 0){
        
        sum_imp_dt <-  imp_dt_fil %>%
          select(variable, VarImp,dataset) %>%
          group_by(variable)%>%
          summarise(Freq = n()/N_samp,
                    Mean_imp = mean(VarImp),
                    Median_imp = median(VarImp),
                    se =  sd(VarImp)/sqrt(n()))
        
        sum_imp_dt <-  sum_imp_dt %>%
          left_join(mp_dt)
        
        
        fil_dt <- sum_imp_dt
        fil_dt$dataset <-  dt_set
        fil_dt$Species <-  ve_sp
        fil_imp_dt_list[[paste0(ve_sp,dt_set)]] <- fil_dt
        
        library(ggprism)
        # Directionality:
        sel_dt <-  t_data[,c(unique(sum_imp_dt$variable),"Y_class")]
        
        #  Only if it contains factored columns:
        check_factor <-   names(Filter(is.factor,sel_dt[,1:(ncol(sel_dt) -1),drop = FALSE]))
        
        library(caret)
        if(!is.null(check_factor)){
          # dummify the data
          dmy <- dummyVars(" ~ .", data = sel_dt[,1:(ncol(sel_dt) -1),drop = FALSE])
          ohe_dt <- data.frame(predict(dmy, newdata = sel_dt[,1:(ncol(sel_dt) -1),drop = FALSE]))
          
          # Remove columns after OHE
          rem_cols <-  c("Abx.group.Others","abx_screen_status.ON.ABX.screen","TRT.norm.Placebo")
          ohe_dt <-  ohe_dt[, !names(ohe_dt) %in% rem_cols,drop=FALSE]
          
          sel_dt <-  cbind(ohe_dt,Y_class = sel_dt$Y_class)
        }
        # Summary
        sel_dt_long <-  sel_dt %>%
          pivot_longer(cols = !c(Y_class))
        summary_dt <-  sel_dt_long %>%
          group_by(Y_class,name)%>%
          summarise(mean = mean(value))%>%
          pivot_wider(id_cols = name,names_from = Y_class,values_from = mean) %>%
          mutate(mean_diff = C - NC)
        summary_dt$dir <-  ifelse(summary_dt$mean_diff > 0, "C", "NC")
        
        # Remove strings:
        summary_dt$name <-  gsub("\\.Vancomycin","",summary_dt$name)
        summary_dt$name <-  gsub("\\.True\\.screen","",summary_dt$name)
        
        
        imp_dt_final <- sum_imp_dt %>%
          left_join(summary_dt %>% select(name,mean_diff,dir), by = c("variable" = "name"))
        
        # Clean up the names of the variables
        imp_dt_final$clean_names <- imp_dt_final$variable
        # Clean up the predictors names
        imp_dt_final$clean_names <- gsub("_Screening","",imp_dt_final$clean_names)
        imp_dt_final$clean_names <- gsub("\\.\\."," ", imp_dt_final$clean_names)
        imp_dt_final$clean_names <- gsub("\\."," ", imp_dt_final$clean_names)
        imp_dt_final$clean_names <- gsub("^X","", imp_dt_final$clean_names)
        imp_dt_final$clean_names <- gsub("total_abx_days_prescreen","Total Days of Abx", imp_dt_final$clean_names)
        imp_dt_final$clean_names <- gsub("n_pretreatment_events","N Rounds of Abx", imp_dt_final$clean_names)
        
        # Flip names from MP to ncbi
        idx_mp <-  which(!is.na(imp_dt_final$ncbi_sp_name))
        if(length(idx_mp) > 0){
          imp_dt_final$clean_names[idx_mp] <-  imp_dt_final$ncbi_sp_name[idx_mp]
        }
        
        imp_dt_final$clean_names <- gsub("\\[","", imp_dt_final$clean_names)
        imp_dt_final$clean_names <- gsub("\\]","", imp_dt_final$clean_names)
        
        
        # NO filter on freq
        imp_dt_final <-  imp_dt_final %>% filter(Freq >= 0.1)
        
        heatmap_fil_dt <- imp_dt_final
        heatmap_fil_dt$dataset <-  dt_set
        heatmap_fil_dt$Species <-  ve_sp
        heatmap_imp_dt_list[[paste0(ve_sp,dt_set)]] <- heatmap_fil_dt
        
        
        p_imp <- NULL
        bar_p <- NULL
        
        if(nrow(imp_dt_final) > 0){
          imp_dt_final <-  imp_dt_final %>% arrange(Mean_imp)
          
          
          imp_dt_final$clean_names <- factor(imp_dt_final$clean_names, levels = imp_dt_final$clean_names)
          imp_dt_final$dir <- factor(imp_dt_final$dir , levels = c("C", "NC"))
          library(ggside)
          library(ggstar)
          library(wesanderson)
          col_R <- c("NC"= "#EE6363","C" ="#00CDCD")
          
          plot_imp_dt <- imp_dt_final
          
          p_imp <- ggplot( plot_imp_dt, aes(y= clean_names, x=Mean_imp))+
            geom_bar(aes(x = Mean_imp, y = clean_names, alpha = Freq ),
                     stat="identity",fill = "darkblue",color = "black",
                     data = plot_imp_dt)+
            geom_star(data = plot_imp_dt,aes(fill = dir,y = clean_names),starshape = 15,size = 8, alpha = 0.8)+
            geom_errorbar(aes(y = clean_names, xmin=Mean_imp-se, xmax=Mean_imp+se),
                          width=.2,colour="black", alpha=0.9) +
            scale_alpha_continuous(limits = c(0,1),range = c(0,1))+
            scale_fill_manual(values = col_R) +
            xlab("Variable Importance")+
            ylab("Features")+
            ggtitle(paste0(ve_sp," (",dt_set,")"))+
            theme_prism(base_size = 22)
          
          
          #print(p_imp)
          library(wesanderson)
          set.seed(100)
          
          bar_dt <-  sel_dt_long 
          # Remove strings:
          bar_dt$name <-  gsub("\\.Vancomycin","",bar_dt$name)
          bar_dt$name <-  gsub("\\.True\\.screen","",bar_dt$name)
          bar_dt <-  bar_dt %>% 
            left_join(imp_dt_final %>% select(variable,clean_names), by = c("name" = "variable"))
          
          # bar_dt <-  bar_dt %>% filter(name %in% levels(imp_dt_final$))
          bar_dt$clean_names <-  factor(bar_dt$clean_names, levels = rev(levels(imp_dt_final$clean_names)))
          
          bar_p <- ggplot(bar_dt,aes(y = value,x = Y_class )) +
            #geom_violin(aes(fill = Y_class),width=1.4, alpha = 0.3) +
            geom_boxplot(aes(y=value, x= Y_class,fill = Y_class),outlier.colour = NA,alpha = 0.3)+
            geom_point( aes(y=value, x= Y_class,fill = Y_class), colour="black", alpha=0.5, size=2,
                        position = position_jitter(w = 0.1, h = 0),
                        shape = 21,stroke = 1) +
            
            facet_wrap(~clean_names,scales = "free", ncol = 4)+
            #facet_wrap(~name~TRT.norm,scales = "free")+
            scale_fill_manual(name = "",values = col_R)+
            ylab("Abundance")+
            xlab("")+
            theme_bw()+
            theme(axis.text.x=element_text(size=10,face="bold",angle = 60,hjust=1),
                  #axis.text.y=element_text(size=10,face="bold"),
                  axis.title.y = element_text(size=10,face="bold"),
                  axis.title.x = element_text(size=10,face="bold"),
                  strip.text.x = element_text(size=10,face="bold"),
                  legend.title=element_text(size=10),
                  legend.text=element_text(size=10),
                  axis.text.y = element_text(size=10,face="bold"))+
            ggtitle(dt_set)
          # Print variable importance
          #  write.csv(imp_dt,paste0(results_folder,"/VIMP_",sp,".csv"))
          p_imp_list[[ve_sp]][[dt_set]] <- p_imp
          p_abun_list[[ve_sp]][[dt_set]] <- bar_p
          
          library(patchwork)
          # print(p_imp / bar_p)
          
        }  
        
      }
      
    }   
    
  }
  # dev.off()
}

# sum_imp_final_dt <-  do.call("rbind",fil_imp_dt_list)
# write.csv(sum_imp_final_dt,paste0(results_folder,"/Imp_vars_models.csv"),row.names = F)

sum_imp_final_dt <-  do.call("rbind",heatmap_imp_dt_list)
write.csv(sum_imp_final_dt,paste0(results_folder,"/Imp_vars_models.csv"),row.names = F)


saveRDS(p_imp_list,paste0(results_folder,"/Imp_plot_list.rds"))
saveRDS(p_abun_list,paste0(results_folder,"/Abun_plot_list.rds"))  
saveRDS(heatmap_imp_dt_list,paste0(results_folder,"/Heatmap_imp_list.rds"))


# # Just print the importance plots
# for(ve_sp in ve303_species){ 
#   pdf(paste0(results_folder,"/",ve_sp,"_Variable_imp_RF.pdf"), width = 10, height = 8)
#   for(dt_set in names(data_list)){
#     p_imp <-  p_imp_list[[ve_sp]][[dt_set]]
#     #p_abun_list[[sel_time]][[dt_set]] <- bar_p
#     print(p_imp )
#     
#   }   
#   
#   dev.off()
# }
# 
# 
# # Just print the box plots
# for(ve_sp in ve303_species){ 
#   pdf(paste0(results_folder,"/",ve_sp,"_Boxplot_imp_RF.pdf"), width = 10, height = 10)
#   for(dt_set in names(data_list)){
#     p_abun <-  p_abun_list[[ve_sp]][[dt_set]]
#     #p_abun_list[[sel_time]][[dt_set]] <- bar_p
#     print(p_abun  )
#     
#   }   
#   
#   dev.off()
# }



