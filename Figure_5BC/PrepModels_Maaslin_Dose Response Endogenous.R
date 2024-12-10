
print("Preparing model types for running through MaasLin - different combinations of excluding post recurrence, post abx and post early (D14/D28) samples")

if(treats == "Dosed and Placebo"){
  
  ### MODEL 0 ####
  # ARRAY FOR ALL TREATMENTS,  KEEPING ALL TIMEPOINTS
  VE303_sp.all_Maaslin_All <- VE303_sp.all_Maaslin 
  
  ### MODEL 1 ####
  # ARRAY FOR ALL TREATMENTS,  KEEPING ALL TIMEPOINTS PRE-RECURRENCE 
  VE303_sp.all_Maaslin_postrecfiltered <- VE303_sp.all_Maaslin %>% filter(collect.pre.post.rec != "Post.rec")
  
  ### MODEL 2 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING TIMEPOINTS THAT ARE BOTH PRE-RECURRENCE +++ EARLY TIMEPOINTS UP TO DAY 28 ONLY
  timefilt <- c( "Day 1","Day 7","Day 14","Day 28")
#  timefilt <- c("Screening","Day 1","Day 7","Day 14","Day 28")
  VE303_sp.all_Maaslin_time_postrecfiltered_D28 <- VE303_sp.all_Maaslin %>% filter((collect.pre.post.rec != "Post.rec") & (Visit.Name %in% timefilt))
  
  ### MODEL 3 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING TIMEPOINTS THAT ARE BOTH PRE-RECURRENCE +++ EARLY TIMEPOINTS UP TO DAY 14 ONLY
  timefilt1 <- c( "Day 1","Day 7","Day 14")
#  timefilt1 <- c("Screening","Day 1","Day 7","Day 14")
  VE303_sp.all_Maaslin_time_postrecfiltered_D14 <- VE303_sp.all_Maaslin %>% filter((collect.pre.post.rec != "Post.rec") & (Visit.Name %in% timefilt1))
  
  ### MODEL 4 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING EARLY UP TO DAY 28 REGARDLESS OF RECURRENCE STATUS
  VE303_sp.all_Maaslin_time_filtered <- VE303_sp.all_Maaslin %>% filter((Visit.Name %in% timefilt))
  
  
  ### MODEL 5 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING EARLY UP TO DAY 14 REGARDLESS OF RECURRENCE STATUS
  VE303_sp.all_Maaslin_time_filtered_D14 <- VE303_sp.all_Maaslin %>% filter((Visit.Name %in% timefilt1))
  
  
  ### MODEL 6 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING DAY 14 REGARDLESS OF RECURRENCE STATUS
  timefilt2 <- c("Day 14")
  VE303_sp.all_Maaslin_time_filtered_D14_only <- VE303_sp.all_Maaslin %>% filter((Visit.Name %in% timefilt2))
  
  ### MODEL 7 ####
  # ARRAY FOR ALL TREATMENTS, KEEPING DAY 7 and later timepoints + Pre-recurrence  
  timefilt2 <- c("Day 14")
  VE303_sp.all_Maaslin_time_postrecfiltered_D1Onward <- VE303_sp.all_Maaslin %>% filter((collect.pre.post.rec != "Post.rec") & (Visit.Name != "Day 1"))
  
  models <- c("Model 0","Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6", "Model 7" )
  model_description <- c("AllTime","Exclude post-recurrence", "Exclude post-recurrence_D28", "Exclude post-recurrence_D14",
                         "Exclude post-D28",  "Exclude post-D14",  "Include D14 Only",  "Exclude Day1 and post-recurrence" )
  Arrays <- list(VE303_sp.all_Maaslin_All,VE303_sp.all_Maaslin_postrecfiltered, VE303_sp.all_Maaslin_time_postrecfiltered_D28, VE303_sp.all_Maaslin_time_postrecfiltered_D14, VE303_sp.all_Maaslin_time_filtered,VE303_sp.all_Maaslin_time_filtered_D14, VE303_sp.all_Maaslin_time_filtered_D14_only, VE303_sp.all_Maaslin_time_postrecfiltered_D1Onward)
  
}


if(treats == "Dosed"){
  
  ### MODEL 0 ####
  # ARRAY FOR ALL TREATMENTS,  KEEPING ALL TIMEPOINTS
  VE303_sp.all_Maaslin_All <- VE303_sp.all_Maaslin_filtered
  
  ### MODEL 1 ####
  # ARRAY FOR ALL TREATMENTS,  KEEPING ALL TIMEPOINTS PRE-RECURRENCE 
  VE303_sp.all_Maaslin_postrecfiltered <- VE303_sp.all_Maaslin_filtered %>% filter(collect.pre.post.rec != "Post.rec")
  
  ### MODEL 2 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING TIMEPOINTS THAT ARE BOTH PRE-RECURRENCE +++ EARLY TIMEPOINTS UP TO DAY 28 ONLY
  timefilt <- c("Screening","Day 1","Day 7","Day 14","Day 28")
  VE303_sp.all_Maaslin_time_postrecfiltered_D28 <- VE303_sp.all_Maaslin_filtered %>% filter((collect.pre.post.rec != "Post.rec") & (Visit.Name %in% timefilt))
  
  ### MODEL 3 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING TIMEPOINTS THAT ARE BOTH PRE-RECURRENCE +++ EARLY TIMEPOINTS UP TO DAY 14 ONLY
  timefilt1 <- c("Screening","Day 1","Day 7","Day 14")
  VE303_sp.all_Maaslin_time_postrecfiltered_D14 <- VE303_sp.all_Maaslin_filtered %>% filter((collect.pre.post.rec != "Post.rec") & (Visit.Name %in% timefilt))
  
  ### MODEL 4 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING EARLY UP TO DAY 28 REGARDLESS OF RECURRENCE STATUS
  VE303_sp.all_Maaslin_time_filtered <- VE303_sp.all_Maaslin_filtered %>% filter((Visit.Name %in% timefilt))
  
  
  ### MODEL 5 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING EARLY UP TO DAY 14 REGARDLESS OF RECURRENCE STATUS
  VE303_sp.all_Maaslin_time_filtered_D14 <- VE303_sp.all_Maaslin_filtered %>% filter((Visit.Name %in% timefilt1))
  
  ### MODEL 6 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING DAY 14 REGARDLESS OF RECURRENCE STATUS
  timefilt2 <- c("Day 14")
  VE303_sp.all_Maaslin_time_filtered_D14_only <- VE303_sp.all_Maaslin_filtered %>% filter((Visit.Name %in% timefilt2))
  
  ### MODEL 7 ####
  # ARRAY FOR ALL TREATMENTS, KEEPING DAY 7 and later timepoints + Pre-recurrence  
  timefilt3 <- c("Day 7","Day 14")
  VE303_sp.all_Maaslin_time_postrecfiltered_D1Onward <- VE303_sp.all_Maaslin_filtered %>% filter((collect.pre.post.rec != "Post.rec") & (Visit.Name %in% timefilt3))
  
  models <- c("Model 0","Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6", "Model 7" )
  model_description <- c("AllTime","Exclude post-recurrence", "Exclude post-recurrence_D28", "Exclude post-recurrence_D14",
                         "Exclude post-D28",  "Exclude post-D14",  "Include D14 Only",  "Exclude Day1 and post-recurrence" )
  Arrays <- list(VE303_sp.all_Maaslin_All,VE303_sp.all_Maaslin_postrecfiltered, VE303_sp.all_Maaslin_time_postrecfiltered_D28, VE303_sp.all_Maaslin_time_postrecfiltered_D14, VE303_sp.all_Maaslin_time_filtered,VE303_sp.all_Maaslin_time_filtered_D14, VE303_sp.all_Maaslin_time_filtered_D14_only, VE303_sp.all_Maaslin_time_postrecfiltered_D1Onward)
  
  
}



if(treats == "Placebo"){
  
  ### MODEL 0 ####
  # ARRAY FOR ALL TREATMENTS,  KEEPING ALL TIMEPOINTS
  VE303_sp.all_Maaslin_All <- VE303_sp.all_Maaslin_filtered_placebo
  
  ### MODEL 1 ####
  # ARRAY FOR ALL TREATMENTS,  KEEPING ALL TIMEPOINTS PRE-RECURRENCE 
  VE303_sp.all_Maaslin_postrecfiltered <- VE303_sp.all_Maaslin_filtered_placebo %>% 
    filter(collect.pre.post.rec != "Post.rec")
  
  ### MODEL 2 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING TIMEPOINTS THAT ARE BOTH PRE-RECURRENCE +++ EARLY TIMEPOINTS UP TO DAY 28 ONLY
  timefilt <- c("Screening","Day 1","Day 7","Day 14","Day 28")
  VE303_sp.all_Maaslin_time_postrecfiltered_D28 <- VE303_sp.all_Maaslin_filtered_placebo %>% 
    filter((collect.pre.post.rec != "Post.rec") & (Visit.Name %in% timefilt))
  
  ### MODEL 3 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING TIMEPOINTS THAT ARE BOTH PRE-RECURRENCE +++ EARLY TIMEPOINTS UP TO DAY 14 ONLY
  timefilt1 <- c("Screening","Day 1","Day 7","Day 14")
  VE303_sp.all_Maaslin_time_postrecfiltered_D14 <- VE303_sp.all_Maaslin_filtered_placebo %>% 
    filter((collect.pre.post.rec != "Post.rec") & (Visit.Name %in% timefilt))
  
  ### MODEL 4 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING EARLY UP TO DAY 28 REGARDLESS OF RECURRENCE STATUS
  VE303_sp.all_Maaslin_time_filtered <- VE303_sp.all_Maaslin_filtered_placebo %>% 
    filter((Visit.Name %in% timefilt))
  
  
  ### MODEL 5 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING EARLY UP TO DAY 14 REGARDLESS OF RECURRENCE STATUS
  VE303_sp.all_Maaslin_time_filtered_D14 <- VE303_sp.all_Maaslin_filtered_placebo %>% 
    filter((Visit.Name %in% timefilt1))
  
  ### MODEL 6 ####
  # ARRAY FOR ALL TREATMENTS, ONLY KEEPING DAY 14 REGARDLESS OF RECURRENCE STATUS
  timefilt2 <- c("Day 14")
  VE303_sp.all_Maaslin_time_filtered_D14_only <- VE303_sp.all_Maaslin_filtered_placebo %>%
    filter((Visit.Name %in% timefilt2))
  
  ### MODEL 7 ####
  # ARRAY FOR ALL TREATMENTS, KEEPING DAY 7 and later timepoints + Pre-recurrence  
  timefilt3 <- c("Day 7","Day 14")
  VE303_sp.all_Maaslin_time_postrecfiltered_D1Onward <- VE303_sp.all_Maaslin_filtered_placebo %>% 
    filter((collect.pre.post.rec != "Post.rec") & (Visit.Name %in% timefilt3))
  
  models <- c("Model 0","Model 1", "Model 2", "Model 3", "Model 4", "Model 5", "Model 6", "Model 7" )
  model_description <- c("AllTime","Exclude post-recurrence", "Exclude post-recurrence_D28", 
                         "Exclude post-recurrence_D14",
                         "Exclude post-D28",  "Exclude post-D14",  "Include D14 Only",  
                         "Exclude Day1 and post-recurrence" )
  Arrays <- list(VE303_sp.all_Maaslin_All,VE303_sp.all_Maaslin_postrecfiltered,
                 VE303_sp.all_Maaslin_time_postrecfiltered_D28,
                 VE303_sp.all_Maaslin_time_postrecfiltered_D14, VE303_sp.all_Maaslin_time_filtered,
                 VE303_sp.all_Maaslin_time_filtered_D14, VE303_sp.all_Maaslin_time_filtered_D14_only, 
                 VE303_sp.all_Maaslin_time_postrecfiltered_D1Onward)
}