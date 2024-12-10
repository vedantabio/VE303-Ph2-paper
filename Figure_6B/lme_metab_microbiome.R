# Load libraries
library(gtools)
library(tidyverse)
# Create a folder for LME analysis
library(here)
mainDir <- here("LME_Microbiome_Metabolites")
dir.create(file.path(mainDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,sep="/")

######################### VE303 ##############################################
# Import microbiome data  :
library(phyloseq)
phy <- readRDS(here("../Processed_data/phy_mic.rds"))
sm_dt <-  data.frame(sample_data(phy))
# Only HD , LD group
#  phy_melt_sel <- phy_melt_sel %>% filter(!TRT.norm %in% c("Placebo"))
#  phy_melt_sel <- phy_melt_sel %>% filter(!Screening_abx_check %in% c("Abx Pre-Screen"))
# Look at only dosed group:
phy_dose <- subset_samples(phy,!TRT.norm %in% c("Placebo") )
phy_dose <- subset_samples(phy_dose,Visit.Name %in% c("Screening","Day 1","Day 7","Day 14","Day 28","Day 56") )
phy_dose <-  prune_taxa(taxa_sums(phy_dose)>0,phy_dose)

# Only select prevalent taxa:
perc_samples <-  0.15
#perc_samples <-  0.05
phy_fil = filter_taxa(phy_dose, function(x) sum(x > 0) > (perc_samples*length(x)), TRUE)

otu_dt <-  data.frame(otu_table(phy_fil)) %>% t() %>% data.frame()
otu_dt$SampleID <-  rownames(otu_dt)
sm_dt_fil <-  data.frame(sample_data(phy_fil))
otu_dt$SampleID == sm_dt_fil$MGS.Sample.ID
otu_dt$SampleID <- paste0(sm_dt_fil$Subject.Number,"_",sm_dt_fil$Visit.Name)

######################################### SCFA #########################################
# Import SCFA :
# Read SCFA abundance data
scfa_dt  <-  read.csv(here("../Processed_data/SCFA.csv"))

scfa_dt$Visit <- as.character(scfa_dt$Visit)
# Limit the time to Day 56
sorted_time <- mixedsort(unique(scfa_dt$Visit))

pre_labels <-  c("Screening",sorted_time[1:which(sorted_time == "Day 1")])
post_labels <- sorted_time[(which(sorted_time == "Day 1")+1):which(sorted_time == "Day 56")]

#pre_labels <-  c("Screening","Day 1")
pre_labels <-  c("Screening","Day 1")
post_labels <- c("Day 7","Day 14", "Day 56")

scfa_dt <- scfa_dt[scfa_dt$Visit %in% c(pre_labels,post_labels),]

scfa_dt$Time <- scfa_dt$Visit 

#scfa_dt$Time <- factor(scfa_dt$Time,levels = c("Screening","Day 1","Day 7","Day 14","Day 56"))
scfa_dt$Time <- factor(scfa_dt$Time,levels = c("Screening","Day 1","Day 7","Day 14","Day 56"))
# Treatment
unique(scfa_dt$TRT.norm)
scfa_dt$Treatment <- factor(as.character(scfa_dt$TRT.norm), levels = c("Placebo" ,"VE303 Low Dose", "VE303 High Dose"))

# Analyte
scfa_dt$Analyte <- as.character(scfa_dt$Analyte)
scfa_dt$SampleID <-  paste0(scfa_dt$Subject.Number,"_",scfa_dt$Visit)
 
# Match it with the MGS sample ID
library(tidyverse)
comb_dt <-  scfa_dt  %>% 
  inner_join(otu_dt, by = c("SampleID"))

levels(scfa_dt$Time)


table(comb_dt$Visit)
# Now Look for association of VE303 abundance across all the SCFA 
# Subset timepoints from Day 1 to Day 56
lme_dt <-  comb_dt[comb_dt$Visit %in% c("Day 7","Day 14"),]


# Loop over three conditions:

# Highly engrafted subjects
# All dosed subjects
condns <-  c("dosed")
condn <-  "dosed"
for(condn in condns){
  
  # Subset VE303 dose group
  lme_dt <- lme_dt[lme_dt$TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose"),]
  #lme_dt <- lme_dt[lme_dt$rec.diagnosis %in% c("Good"),]
  
  lme_dt <-  data.frame(lme_dt)
  
  # Only keep the high engrafted subjects
  if(condn == "high_engraft"){
    lme_dt <-  lme_dt[lme_dt$Sub_ID %in% eng_dt$Subject.Number,]
  }
  
  # lme_dt <- lme_dt[lme_dt$rec.diagnosis %in% c("Good"),]
  length(unique(lme_dt$Sub_ID))
  # For each metabolite, check for its association with VE303 strains:
  unique(lme_dt$Analyte)
  #mod_var <-  "Screening_abx_check"
  result_lme <-  list()
  scfa <- unique(lme_dt$Analyte)
  scfa_var <- scfa[3]
  for (scfa_var in scfa){
    
    scfa_sel <- lme_dt[lme_dt$Analyte %in% scfa_var,]
    strains  <-  names(scfa_sel)[30:ncol(scfa_sel)]
    list_strain <-  list()
    strain <- strains[5]
    for(strain in  strains){
      print(strain)    
      sum_mod_dt <- tryCatch({
        
        # Remove the species with NA absolute abundance
        scfa_sel_st <- scfa_sel[complete.cases(scfa_sel[, "Est_Abun"]),]
        scfa_sel_st <- scfa_sel[complete.cases(scfa_sel[, c(strain)]),]
        #scfa_sel_st <- scfa_sel[scfa_sel[, c(strain)] != 0 , ]
        
        
        scfa_sel_st$Abundance <- log10(scfa_sel_st$Est_Abun)
        scfa_sel_st$St_Abundance <-  asin(sqrt(scfa_sel_st[, c(strain)]))
        
        scfa_sel_st$Subject.Number <- factor(scfa_sel_st$Subject.Number)
        scfa_sel_st$Time <- scale(as.numeric(gsub("Day ","",scfa_sel_st$Time)),center = F)
        
        
        #scfa_sel_time$rec.diagnosis <-  factor(scfa_sel_time$rec.diagnosis, c("Bad","Good"))
        library("nlme")
        fml <- as.formula( paste( "Abundance", "~", paste("St_Abundance + Time") ) )
        mod_bac <- lme(fml,random = ~ 1|Subject.Number, scfa_sel_st)
        sum_mod <-  summary(mod_bac)
        #print(anova.lme(mod_bac))
        sum_mod_dt <- data.frame(sum_mod$tTable)
        sum_mod_dt$scfa_var <- scfa_var
        sum_mod_dt$Var <-  rownames(sum_mod_dt)
        sum_mod_dt$Strain <- strain
        #sum_mod_dt$trt <-  trt
        
        
        list_strain[[strain]] <-  sum_mod_dt
        
      }, error = function(err){
        sum_mod_dt <- NA
        sum_mod_dt
      }
      )
      list_strain[[strain]] <-  sum_mod_dt
    }
    
    result_st_dt <- do.call("rbind", list_strain)
    names(result_st_dt)[5] <- "pval"
    result_st_dt <- result_st_dt[!is.na(result_st_dt$Var),]
    #result_time_dt <-  result_time_dt[result_time_dt$Var != "(Intercept)",]
    #result_time_dt$pval <-  p.adjust(result_time_dt$pval,method = "BH")
    result_lme[[scfa_var]] <- result_st_dt
    
  }
  
  result_lme_dt <- do.call("rbind", result_lme)
  
  result_lme_dt$adj.pval <- p.adjust(result_lme_dt$pval,method = "BH")
  # result_lme_dt <- result_lme_dt %>%
  #   group_by(Strain) %>%
  #   mutate(adj.pval =  p.adjust(pval,method = "BH"))
  result_lme_dt <-  result_lme_dt[result_lme_dt$Var != "(Intercept)",]
  result_lme_dt <-  result_lme_dt[result_lme_dt$Var != "Time",]
  
  write.csv(result_lme_dt,paste0(results_folder,"/",paste0("SCFA_lme_VE303_",condn,".csv")))
  
  
}
######################### Bile Acid #################################################################### 
# Import BA:
# Read SCFA abundance data
scfa_dt  <-  read.csv(here("../Processed_data/BA.csv"))

scfa_dt$Visit <- as.character(scfa_dt$Visit)
# Limit the time to Day 56
sorted_time <- mixedsort(unique(scfa_dt$Visit))

pre_labels <-  c("Screening",sorted_time[1:which(sorted_time == "Day 1")])
post_labels <- sorted_time[(which(sorted_time == "Day 1")+1):which(sorted_time == "Day 56")]

#pre_labels <-  c("Screening","Day 1")
pre_labels <-  c("Screening","Day 1")
post_labels <- c("Day 7","Day 14", "Day 56")

scfa_dt <- scfa_dt[scfa_dt$Visit %in% c(pre_labels,post_labels),]

scfa_dt$Time <- scfa_dt$Visit 

#scfa_dt$Time <- factor(scfa_dt$Time,levels = c("Screening","Day 1","Day 7","Day 14","Day 56"))
scfa_dt$Time <- factor(scfa_dt$Time,levels = c("Screening","Day 1","Day 7","Day 14","Day 56"))
# Treatment
unique(scfa_dt$TRT.norm)
scfa_dt$Treatment <- factor(as.character(scfa_dt$TRT.norm), levels = c("Placebo" ,"VE303 Low Dose", "VE303 High Dose"))

# Analyte
scfa_dt$Analyte <- as.character(scfa_dt$Analyte)

# Remove 
scfa_dt <- scfa_dt[!scfa_dt$Analyte %in% c("Isodeoxycholic Acid","Glycolithocholic Acid"),]
scfa_dt$SampleID <-  paste0(scfa_dt$Subject.Number,"_",scfa_dt$Visit)
# Match it with the MGS sample ID
library(tidyverse)
comb_dt <-  scfa_dt  %>% 
  inner_join(otu_dt, by = c("SampleID"))

levels(scfa_dt$Time)


table(comb_dt$Visit)
# Now Look for association of VE303 abundance across all the SCFA 
# Subset timepoints from Day 1 to Day 56
lme_dt <-  comb_dt[comb_dt$Visit %in% c("Day 7","Day 14"),]

# Loop over three conditions:

# Highly engrafted subjects
# All dosed subjects
condns <-  c("dosed")
condn <-  "dosed"
for(condn in condns){
  
  
  # Subset VE303 dose group
  lme_dt <- lme_dt[lme_dt$TRT.norm %in% c("VE303 Low Dose", "VE303 High Dose"),]
  #lme_dt <- lme_dt[lme_dt$rec.diagnosis %in% c("Good"),]
  
  lme_dt <-  data.frame(lme_dt)
  
  # Only keep the high engrafted subjects
  if(condn == "high_engraft"){
    lme_dt <-  lme_dt[lme_dt$Sub_ID %in% eng_dt$Subject.Number,]
  }
  
  # For each metabolite, check for its association with VE303 strains:
  
  unique(lme_dt$Analyte)
  #mod_var <-  "Screening_abx_check"
  result_lme <-  list()
  scfa <- unique(lme_dt$Analyte)
  scfa_var <- scfa[3]
  for (scfa_var in scfa){
    
    scfa_sel <- lme_dt[lme_dt$Analyte %in% scfa_var,] 
    strains  <-  names(scfa_sel)[30:ncol(scfa_sel)]
    list_strain <-  list()
    strain <- strains[5]
    for(strain in  strains){
      
      print(strain)
      sum_mod_dt <- tryCatch({
        
        # Remove the species with NA absolute abundance
        scfa_sel_st <- scfa_sel[complete.cases(scfa_sel[, "Est_Abun"]),]
        scfa_sel_st <- scfa_sel[complete.cases(scfa_sel[, c(strain)]),]
        #scfa_sel_st <- scfa_sel[scfa_sel[, c(strain)] != 0 , ]
        
        
        scfa_sel_st$Abundance <- log10(scfa_sel_st$Est_Abun)
        scfa_sel_st$St_Abundance <-  asin(sqrt(scfa_sel_st[, c(strain)]))
        
        scfa_sel_st$Subject.Number <- factor(scfa_sel_st$Subject.Number)
        scfa_sel_st$Time <- scale(as.numeric(gsub("Day ","",scfa_sel_st$Time)),center = F)
        
        
        #scfa_sel_time$rec.diagnosis <-  factor(scfa_sel_time$rec.diagnosis, c("Bad","Good"))
        library("nlme")
        fml <- as.formula( paste( "Abundance", "~", paste("St_Abundance + Time") ) )
        mod_bac <- lme(fml,random = ~ 1|Subject.Number, scfa_sel_st)
        sum_mod <-  summary(mod_bac)
        #print(anova.lme(mod_bac))
        sum_mod_dt <- data.frame(sum_mod$tTable)
        sum_mod_dt$scfa_var <- scfa_var
        sum_mod_dt$Var <-  rownames(sum_mod_dt)
        sum_mod_dt$Strain <- strain
        #sum_mod_dt$trt <-  trt
        
        
        list_strain[[strain]] <-  sum_mod_dt
        
      }, error = function(err){
        sum_mod_dt <- NA
        sum_mod_dt
      }
      )
      list_strain[[strain]] <-  sum_mod_dt
    }
    
    result_st_dt <- do.call("rbind", list_strain)
    names(result_st_dt)[5] <- "pval"
    result_st_dt <- result_st_dt[!is.na(result_st_dt$Var),]
    #result_time_dt <-  result_time_dt[result_time_dt$Var != "(Intercept)",]
    #result_time_dt$pval <-  p.adjust(result_time_dt$pval,method = "BH")
    result_lme[[scfa_var]] <- result_st_dt
    
  }
  
  result_lme_dt <- do.call("rbind", result_lme)
  
  result_lme_dt$adj.pval <- p.adjust(result_lme_dt$pval,method = "BH")
  
  # result_lme_dt <- result_lme_dt %>%
  #   group_by(Strain) %>%
  #   mutate(adj.pval =  p.adjust(pval,method = "BH"))
  result_lme_dt <-  result_lme_dt[result_lme_dt$Var != "(Intercept)",]
  result_lme_dt <-  result_lme_dt[result_lme_dt$Var != "Time",]
  
  result_lme_dt <-  result_lme_dt %>%
    mutate(.,Class = case_when(scfa_var %in% c("Cholic Acid", "Chenodeoxycholic Acid") ~ "Primary BA",
                               
                               scfa_var %in% c("Glycocholic Acid", "Glycochenodeoxycholic Acid", "Taurocholic Acid", "Taurochenodeoxycholic Acid") ~ "Primary BA",
                               
                               scfa_var %in% c("Lithocholic Acid", "Deoxycholic Acid", "Isodeoxycholic Acid", "Ursodeoxycholic Acid") ~ "Secondary BA",
                               
                               scfa_var %in% c("Alloiso_Isolithocholic Acid", "Dehydrolithocholic Acid", "Glycodeoxycholic Acid", "Glycolithocholic Acid", "Glycoursodeoxycholic Acid", "Taurodeoxycholic Acid", "Taurolithocholic Acid", "Tauroursodeoxycholic Acid") ~ "Secondary BA"))%>%
    data.frame()
  
  
  
  write.csv(result_lme_dt,paste0(results_folder,"/",paste0("BA_lme_VE303_",condn,".csv")))
  
  
}

# Combine both results into a heatmap
# For all dosed
result_lme_dt_scfa  <- read.csv(here("LME_Microbiome_Metabolites/SCFA_lme_VE303_dosed.csv"))
result_lme_dt_scfa$Class <-  "SCFA"
result_lme_dt_ba  <- read.csv(here("LME_Microbiome_Metabolites/BA_lme_VE303_dosed.csv"))
result_lme_dt_ba <-  result_lme_dt_ba %>%
  mutate(.,Class = case_when(scfa_var %in% c("Cholic Acid", "Chenodeoxycholic Acid") ~ "Primary BA",
                             
                             scfa_var %in% c("Glycocholic Acid", "Glycochenodeoxycholic Acid", "Taurocholic Acid", "Taurochenodeoxycholic Acid") ~ "Primary BA",
                             
                             scfa_var %in% c("Lithocholic Acid", "Deoxycholic Acid", "Isodeoxycholic Acid", "Ursodeoxycholic Acid") ~ "Secondary BA",
                             
                             scfa_var %in% c("Alloiso_Isolithocholic Acid", "Dehydrolithocholic Acid", "Glycodeoxycholic Acid", "Glycolithocholic Acid", "Glycoursodeoxycholic Acid", "Taurodeoxycholic Acid", "Taurolithocholic Acid", "Tauroursodeoxycholic Acid") ~ "Secondary BA"))%>%
  data.frame()

result_lme_dt <-  rbind(result_lme_dt_scfa,result_lme_dt_ba)
#result_lme_dt <-  result_lme_dt_scfa

# Display the lme results for time:
sig_lme_dt <-  result_lme_dt[result_lme_dt$adj.pval < 0.05,]
unique(sig_lme_dt$scfa_var)
# Remove Sach
sig_lme_dt <-  sig_lme_dt[!grepl("Saccharomyces.cerevisiae",sig_lme_dt$Strain),]

# Clean up species name and replace MPs with NCBI names
# Species naming 
sp_name_dt <-  readxl::read_excel("../Processed_data/Feb 2021 Vedanta database.xlsx")%>%
  select(GTDB_Species,NCBI_classification,StrainID,GTDB_closest_placement_taxonomy)
sp_name_dt$GTDB_closest_species <- gsub(".*s__","",sp_name_dt$GTDB_closest_placement_taxonomy)
sp_name_dt$NCBI_species <- gsub(".*s__","",sp_name_dt$NCBI_classification)

sp_name_dt <- sp_name_dt[grep("MP|HF",sp_name_dt$StrainID),]

# From analysis
species_name <- unique(sig_lme_dt$Strain)
mp_names <-  species_name[grep("MP|HF",species_name)]

mp_dt <-  data.frame(variable= mp_names ,MP_ID = mp_names)
mp_dt$MP_ID <-  gsub(".*\\.\\.MP","",mp_dt$MP_ID)
mp_dt$MP_ID <-  gsub(".*\\.\\.HF","",mp_dt$MP_ID)
#mp_dt$MP_ID <-  gsub("_.*","",mp_dt$MP_ID)
mp_dt$MP_ID[grep("MP",mp_dt$variable)] <-  paste0("MP",mp_dt$MP_ID[grep("MP",mp_dt$variable)])
mp_dt$MP_ID[grep("HF",mp_dt$variable)] <-  paste0("HF",mp_dt$MP_ID[grep("HF",mp_dt$variable)])

gtdb_sp_name <-  c()
ncbi_sp_name <- c()
for(mp_id in mp_dt$MP_ID){
  gtdb_sp <- unique( sp_name_dt[grep(mp_id,sp_name_dt$StrainID),]$GTDB_closest_species)
  ncbi_sp <- unique( sp_name_dt[grep(mp_id,sp_name_dt$StrainID),]$NCBI_species)
  # print(gtdb_sp)
  #print(ncbi_sp)
  gtdb_sp_name <-  c(gtdb_sp_name,gtdb_sp)
  ncbi_sp_name <-  c(ncbi_sp_name,ncbi_sp)
}
mp_dt$gtdb_sp_name <- gtdb_sp_name
mp_dt$ncbi_sp_name <- ncbi_sp_name
mp_dt$ncbi_sp_name <- gsub(";.*","",mp_dt$ncbi_sp_name)
mp_dt$ncbi_sp_name<- paste0(mp_dt$ncbi_sp_name, "_",mp_dt$MP_ID)

sig_lme_dt <-  sig_lme_dt %>%
  left_join(mp_dt,by = c("Strain" = "variable"))

sig_lme_dt$clean_names <-  sig_lme_dt$Strain

# Flip names from MP to ncbi
idx_mp <-  which(!is.na(sig_lme_dt$ncbi_sp_name))
if(length(idx_mp) > 0){
  sig_lme_dt$clean_names[idx_mp] <-  sig_lme_dt$ncbi_sp_name[idx_mp]
}
sig_lme_dt$clean_names <-  make.names(sig_lme_dt$clean_names)
sig_lme_dt$clean_names <- gsub("\\.\\."," ", sig_lme_dt$clean_names)
sig_lme_dt$clean_names <- gsub("\\."," ", sig_lme_dt$clean_names)
sig_lme_dt$clean_names <- gsub("^X","", sig_lme_dt$clean_names)



# Display the heatmap of VE303 Vs SCFA:
ht_dt <-  sig_lme_dt %>%
  select(scfa_var,Strain,Value)%>%
  pivot_wider(names_from = Strain,values_from = Value,values_fill = 0)%>%
  data.frame() 

rownames(ht_dt) <- ht_dt$scfa_var
ht_dt$scfa_var <-  NULL
ht_dt[ht_dt > 0] <- 1
ht_dt[ht_dt < 0] <- -1


# Now Make a Heatmap of these data:
mat_logic <-  ht_dt

mat_logic[mat_logic == 1] <- "Positive"
mat_logic[mat_logic == -1] <- "Negative"
mat_logic[mat_logic == 0] <- "NS"

cols <- structure(c("#3b5998","#990000", "white"), names = c("Positive","Negative","NS"))

cols <-  c("Positive" = "#3b5998","Negative" = "#990000","NS" = "white")
library(ComplexHeatmap)

colnames(mat_logic) <- gsub("\\.","-",colnames(mat_logic))


#unique(result_lme_dt_scfa$scfa_var)

ann_dt <- sig_lme_dt %>% select(scfa_var,Class) %>% unique()
ann_dt <-  ann_dt[match(ann_dt$scfa_var,rownames(mat_logic)),]
split_rows <-  factor(ann_dt$Class,levels = c("SCFA","Primary BA","Secondary BA"))
levels(split_rows) <- c("SCFA","PBA","SBA")

# split_rows <-  ifelse(rownames(mat_logic) %in% unique(result_lme_dt_scfa$scfa_var),"SCFA","BA")
# split_rows <-  factor(split_rows, levels = c("SCFA","BA"))

# Import taxa 
phy <- readRDS(here("../Processed_data/phy_mic.rds"))
tax_dt <-  data.frame(tax_table(phy)@.Data)
tax_dt$Species_t <-  make.names(tax_dt$Species)

colnames(mat_logic) <- make.names(colnames(mat_logic))

tax_dt <-  tax_dt[tax_dt$Species_t %in% colnames(mat_logic),]
tax_dt <- tax_dt[match(colnames(mat_logic),tax_dt$Species_t),]

# Merge mp_dt
tax_dt <-  tax_dt %>%
  left_join(sig_lme_dt %>% select(Strain, clean_names) %>% unique(), by = c("Species_t"= "Strain"))


split_cols <-  factor(tax_dt$Class)

colnames(mat_logic) <- tax_dt$clean_names

all_species <-  colnames(mat_logic)
ve303_sp <-  grep("VE303",all_species,value = T)

ve303_col <- ifelse(!all_species %in% ve303_sp,"black","darkred")

ht1 = Heatmap(as.matrix(mat_logic), name = "Coeff", column_title = NA, 
              #clustering_distance_rows = "euclidean",
              column_split = split_cols,
              column_names_gp = gpar(col = ve303_col),
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
              row_names_side = "left", color_space = "LAB",
              col=cols,
              width=2, 
              #column_title_gp = gpar(fontsize = 8,face = "bold"),
              #row_names_max_width = unit(8, "cm"),show_column_names= T,
              row_names_max_width = max_text_width(rownames(mat_logic),
                                                   gp = gpar(fontsize = 12)),
              column_names_max_height = max_text_width(colnames(mat_logic),
                                                       gp = gpar(fontsize = 12)),
              show_heatmap_legend = F,
              heatmap_height = nrow(mat_logic)*unit(8, "mm"),
              na_col="white")

lgd <- Legend(at = c("Positive", "Negative"), title = "Relation", legend_gp = gpar(fill = cols))

pdf(paste(results_folder,"/Heatmap_Comb_Species_logical_dosed.pdf",sep=""),height =10, width = 16)
draw(ht1,annotation_legend_list = list(lgd))
dev.off()

