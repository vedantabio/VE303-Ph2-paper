## Create Table for PD results across all taxonomic levels

print("Function called _cleanup_summary_tables_ for loading older analysis CSV filed for Volcano Plots")
print("The combined analysis results csv loaded to make Volcano plots for VE303 Dosed Bad vs VE303 Dosed Good 
      is derived from a Dosed + Placebo data cut. The model used for fixed effects was Species ~ rec.TRT.combined, TRT.norm, Batch. 
      If wanting to load the results on ONLY Dosed group, modify --Dosed and Placebo-- in the filename below to read --Dosed combined--. 
      The date remains the same.")

cleanup_summary_tables <- function(folder_name, tax_l){
  labelload <- paste("all_results", tax_l, ".tsv", sep = "_")
  model_result = read_tsv(here(paste(folder_name, labelload, sep = "/")))
  models_list <- list()
  
  for (x1 in (model_result$metadata %>% unique())) {
    model_result1 <- model_result %>% 
      filter(metadata == x1)
    
    model_result2 <- model_result1[
      with(model_result1, order(value, pval)),] %>%
      filter(pval < 1.0)
    
    models_list <-  append(models_list, list(model_result2)) }
  
  models_1 <- lapply(models_list, function(x) cbind(" "=rownames(x), x))
  ## Writes to multi sheet excel file ###
  writexl::write_xlsx(models_1, here(results, paste("model result pfiltered" , ".xlsx", sep=" ")) )
  
  for (x in 1:length(model_result$metadata %>% unique()) ) {
    model_test <- models_1 [[x]]
    covariate <- model_test$metadata %>% unique()
    
    if(covariate == "rec.TRT.combined"){
      pm1 <- model_test %>% filter(value == reference_level) %>%
        select( feature, coef, stderr, pval, qval)
      pm1$Covariate <- "VE303 Dosed Bad vs Good"
      pm1 <- pm1[c( "feature", "Covariate", "coef", "stderr", "pval", "qval" )]
    }
    
    if(covariate == "VE303"){
      pm2  <- model_test %>% 
        select( feature, coef, stderr, pval, qval)
      pm2$Covariate <- "VE303 abundance"
      pm2 <- pm2[c( "feature", "Covariate", "coef", "stderr", "pval", "qval" )]
    }
  }
  return(list(pm1, pm2))}

