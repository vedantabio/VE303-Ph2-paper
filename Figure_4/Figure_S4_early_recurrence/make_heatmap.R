# Highlight the recurrence data for each subject:
# Create a directories to store results from the analysis
library(here)
mainDir <- here("Heatmaps")
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
library(lemon)

imp_dt_list <-  readRDS(here("early_recurrence_model_importance/Heatmap_imp_list.rds"))

sum_imp_final_dt <-  do.call("rbind",imp_dt_list)

# Import Performance of the models:
perf_dt <-  read.csv(here("early_recurrence_model_importance/Performance_colonization.csv"))
head(perf_dt)
head(sum_imp_final_dt)

# Filter the importance:
imp_fil_dt <-  sum_imp_final_dt %>%
  left_join(perf_dt %>% select(dataset,Time, Metric,mean))%>%
  select(clean_names,dataset,dir,Time,Freq,Metric,mean) %>%
  filter(mean > 0.5)

# Remove Mic_ from the group
imp_fil_dt$dataset <-  gsub("Mic_","",imp_fil_dt$dataset)
unique(imp_fil_dt$dataset)

clinical_factors <-  c("Age", "BMI","total recurrences","abx_screen_status","n_washout_days","N Rounds of Abx","Total Days of Abx",
                       "Abx group","treat day len abx")
imp_fil_dt$Categories <- ifelse(imp_fil_dt$clean_names %in% clinical_factors,"Clinical factors","Others")
# Replace others with datasets:
imp_fil_dt$Categories[imp_fil_dt$Categories == "Others"] <- imp_fil_dt$dataset[imp_fil_dt$Categories == "Others"]


# Replace Total days of Abx to "N days of non-CDI Abx"
imp_fil_dt$clean_names[imp_fil_dt$clean_names == "Total Days of Abx"] <- "N Days of non-CDI Abx"
# Replace N Rounds of Abx to "N Rounds of non-CDI Abx"
imp_fil_dt$clean_names[imp_fil_dt$clean_names == "N Rounds of Abx"] <- "N Rounds of non-CDI Abx"




# Remove "TGFa" ,"TNFb","EGF"
imp_fil_dt <- imp_fil_dt[!imp_fil_dt$clean_names %in% c("TGFa" ,"TNFb","EGF"),]


# Replace names of Cyokines :
cyt_idx <- imp_fil_dt$Categories == "Cytokines"
imp_fil_dt$clean_names[cyt_idx] <-  gsub(" ","-",imp_fil_dt$clean_names[cyt_idx])
unique(imp_fil_dt$clean_names[cyt_idx])
imp_fil_dt$clean_names[imp_fil_dt$clean_names == "MIG-CXCL9"] <- "CXCL9"
imp_fil_dt$clean_names[imp_fil_dt$clean_names == "MIP-1b"] <- "MIP-1\u03B2"
imp_fil_dt$clean_names[imp_fil_dt$clean_names == "SDF-1a-b"] <- "SDF-1\u03B1+\u03B2"
imp_fil_dt$clean_names[imp_fil_dt$clean_names == "MIP-1d"] <- "MIP-1\u03B4"
imp_fil_dt$clean_names[imp_fil_dt$clean_names == "TGFb2"] <- "TGF-\u03B22"

# Remove Saccharomyces:
imp_fil_dt <-  imp_fil_dt[!grepl("Sacc",imp_fil_dt$clean_names),]

#write.csv(imp_fil_dt,paste0(results_folder,"/","imp_fil_dt.csv"))

# Now pivot tables to create a matrix for a heatmap:
imp_dt_w <- imp_fil_dt %>%
  select(clean_names,Freq,Time,dir)%>%
  group_by(clean_names,Time,dir) %>%
   summarise(Freq = max(abs(Freq)))
imp_dt_w$Freq <-  imp_dt_w$Freq* (ifelse(imp_dt_w$dir == "R",1,-1))
  
imp_dt_w$Time <-  factor(imp_dt_w$Time,levels = c("Screening","Day 1","Day 7","Day 14"))

  imp_dt_w <- imp_dt_w %>%
  pivot_wider(id_cols = clean_names, names_from = Time,values_from = Freq,values_fill = 0)%>%
  data.frame()

# Now Make a Heatmap of these data:
mat_logic <-  imp_dt_w
rownames(mat_logic) <-  imp_dt_w$clean_names
mat_logic$clean_names <- NULL

# Remove X
rownames(mat_logic) <- gsub("^X","",rownames(mat_logic))

library(ComplexHeatmap)

colnames(mat_logic) <- gsub("\\."," ",colnames(mat_logic))

mat_logic <- mat_logic[,c("Screening","Day 1","Day 7")]

meta_cat <-  imp_fil_dt %>%
  select(clean_names, Categories)%>%
  unique()%>%
  data.frame()
meta_cat <- meta_cat[match(imp_dt_w$clean_names,meta_cat$clean_names),]

split_rows <-  meta_cat$Categories
split_rows <-  factor(split_rows,levels = rev(c("Clinical factors","Vancomycin","SCFA","BA","Cytokines","Order","Class","Genus","Species")))

rec.cols <- c("No" = "#00CDCD", "Yes" = "#EE6363", "Good" = "#00CDCD", "Bad" = "#EE6363")

library(circlize)
col_mat = colorRamp2(c(-1, 0, 1), c("#00CDCD", "white", "#EE6363"))

#ve303_col <- ifelse(!all_species %in% ve303_sp,"black","darkred")

rnames <-  rownames(mat_logic)
exp_rnames <- paste0(paste0("'", rnames, "'"), collapse = ",")

ht1 = Heatmap(as.matrix(mat_logic), name = "Coeff", column_title = NA, 
              row_split = split_rows,
              row_gap = unit(2, "mm"),border = TRUE,
              row_title_rot = 0,
              column_title_rot = 60,
              column_names_rot = 60,
              row_title_gp = gpar(fontface = "bold"),
              row_names_gp = gpar(fontface = "bold",fontsize = 9),
              #column_names_max_height =  unit(8, "cm"),
              #left_annotation = ha1,
              #clustering_method_rows = "complete",
              cluster_rows = T,
              cluster_row_slices =   F,
              show_row_dend = F,
              row_names_side = "left", color_space = "LAB",
              col=col_mat,
              cluster_columns = F,
              width=2, 
              #column_title_gp = gpar(fontsize = 8,face = "bold"),
              #row_names_max_width = unit(8, "cm"),show_column_names= T,
              row_names_max_width = max_text_width(rownames(mat_logic),
                                                   gp = gpar(fontsize = 12)),
              column_names_max_height = max_text_width(colnames(mat_logic),
                                                       gp = gpar(fontsize = 12)),
              show_heatmap_legend = F,
              #heatmap_height = nrow(mat_logic)*unit(8, "mm"),
              na_col="white")
#draw(ht1)

lgd <- Legend(at = c(-1,0,1), title = "Diagnosis", labels = c("Non Recurrent","Not Significant","Recurrent"),col_fun = col_mat )

# pdf(paste(results_folder,"/Heatmap_Summary_recurrence_AUC_0.5.pdf",sep=""),height =12, width = 7)
# draw(ht1,annotation_legend_list = list(lgd))
# dev.off()
library(Cairo)
cairo_pdf(paste(results_folder,"/Heatmap_Summary_recurrence_AUC_0.5.pdf",sep=""),height =13, width = 7)
draw(ht1,annotation_legend_list = list(lgd))
dev.off()



