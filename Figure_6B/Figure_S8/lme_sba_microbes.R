library(tidyverse)
library(Biostrings)
library(nlme)
library(gtools)
library(circlize)
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(phyloseq)
library(here )

# Create a directory to store results from the analysis
mainDir <- here("SBA_Microbes_Analysis")
dir.create(file.path(mainDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,sep="/")

# Read Processed data:
phy_sba <- readRDS(here("../Processed_data/shortbred/phy_sba.rds"))
phy_sba <- subset_samples(phy_sba,! TRT.norm %in% c("Placebo") )
phy_sba <- subset_samples(phy_sba,Visit.Name %in% c("Day 7","Day 14","Day 28","Day 56") )

perc_samples <-  0.1
perc_samples*nsamples(phy_sba)
phy_sba_fil = filter_taxa(phy_sba, function(x) sum(x > 0) > (perc_samples*length(x)), TRUE)

sba_dt <-  data.frame(otu_table(phy_sba_fil)) %>% t() %>% data.frame()
sba_dt$SampleID <-  rownames(sba_dt)
sm_dt_fil <-  data.frame(sample_data(phy_sba_fil))
sba_dt$SampleID == sm_dt_fil$MGS.Sample.ID
sba_dt$SampleID <- paste0(sm_dt_fil$Subject.Number,"_",sm_dt_fil$Visit.Name)

melt_sba_dt <-  psmelt(phy_sba_fil)
melt_sba_dt$SampleID <-  paste0(melt_sba_dt$Subject.Number,"_",melt_sba_dt$Visit.Name)
melt_sba_dt <-  melt_sba_dt %>% select(OTU,SampleID,Subject.Number,Visit.Name,Abundance)
names(melt_sba_dt)[1] <- "SBA"
# Import microbiome data
phy_mic <- readRDS(here("../Processed_data/phy_mic.rds"))
sm_dt <-  data.frame(sample_data(phy_mic))
# Only HD , LD group
#  phy_melt_sel <- phy_melt_sel %>% filter(!TRT.norm %in% c("Placebo"))
#  phy_melt_sel <- phy_melt_sel %>% filter(!Screening_abx_check %in% c("Abx Pre-Screen"))
# Look at only dosed group:
phy_dose <- subset_samples(phy_mic,!TRT.norm %in% c("Placebo") )
phy_dose <- subset_samples(phy_dose,Visit.Name %in% c("Day 7","Day 14","Day 28","Day 56") )
phy_dose <-  prune_taxa(taxa_sums(phy_dose)>0,phy_dose)

# Only select prevalent taxa:
perc_samples <-  0.1
perc_samples*nsamples(phy_dose)
phy_mic_fil = filter_taxa(phy_dose, function(x) sum(x > 0) > (perc_samples*length(x)), TRUE)

otu_dt <-  data.frame(otu_table(phy_mic_fil)) %>% t() %>% data.frame()
otu_dt$SampleID <-  rownames(otu_dt)
sm_dt_fil <-  data.frame(sample_data(phy_mic_fil))
otu_dt$SampleID == sm_dt_fil$MGS.Sample.ID
otu_dt$SampleID <- paste0(sm_dt_fil$Subject.Number,"_",sm_dt_fil$Visit.Name)

# Combine two data
# Analyte
# Match it with the MGS sample ID
library(tidyverse)
comb_dt <-  melt_sba_dt  %>% 
  inner_join(otu_dt, by = c("SampleID"))

lme_dt <-  comb_dt

# tr|A0A2N0UR24|A0A2N0UR24_9FIRM
# Dorea.sp..VE303.06

result_lme <-  list()
scfa <- unique(lme_dt$SBA)
scfa_var <- scfa[3]
for (scfa_var in scfa){
  scfa_sel <- lme_dt[lme_dt$SBA %in% scfa_var,]
  strains  <-  names(scfa_sel)[6:ncol(scfa_sel)]
  list_strain <-  list()
  strain <- strains[5]
  for(strain in  strains){
    print(strain)    
    sum_mod_dt <- tryCatch({
      
      # Remove the species with NA absolute abundance
      scfa_sel_st <- scfa_sel[complete.cases(scfa_sel[, "Abundance"]),]
      scfa_sel_st <- scfa_sel[complete.cases(scfa_sel[, c(strain)]),]
      #scfa_sel_st <- scfa_sel[scfa_sel[, c(strain)] != 0 , ]
      scfa_sel_st$Abundance <- log(scfa_sel_st$Abundance + 0.001)
      scfa_sel_st$St_Abundance <-  asin(sqrt(scfa_sel_st[, c(strain)]))
      
      scfa_sel_st$Subject.Number <- factor(scfa_sel_st$Subject.Number)
      scfa_sel_st$Time <- scale(as.numeric(gsub("Day ","",scfa_sel_st$Visit.Name)),center = F)
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
  result_lme[[scfa_var]] <- result_st_dt
  
} 

saveRDS(result_lme,paste0(results_folder,"/",paste0("LME_Microbes_SBA.rds")))
result_lme_dt <- do.call("rbind", result_lme)
result_lme_dt <-  result_lme_dt[grep("St_Abun",result_lme_dt$Var),]
result_lme_dt <- result_lme_dt %>%
  group_by(scfa_var)%>%
  mutate(adj.pval = p.adjust(pval, method='BH'))
#result_lme_dt$adj.pval <- p.adjust(result_lme_dt$pval,method = "BH")

write.csv(result_lme_dt,paste0(results_folder,"/",paste0("SBA_lme_Microbes.csv")))

# Display the lme results for time:
sig_lme_dt <-  result_lme_dt[result_lme_dt$adj.pval < 0.05,]
unique(sig_lme_dt$scfa_var)

# Display the heatmap of VE303 Vs SCFA:
ht_dt <-  sig_lme_dt %>%
  select(scfa_var,Strain,Value)%>%
  pivot_wider(names_from = Strain,values_from = Value,values_fill = 0)%>%
  data.frame() 

rownames(ht_dt) <- ht_dt$scfa_var
ht_dt$scfa_var <-  NULL

# Remove ht_dt
table(colSums(ht_dt != 0) > 2)
idx_col <-  colSums(ht_dt != 0) > 2
table(rowSums(ht_dt != 0)>2)
idx_row <- rowSums(ht_dt != 0)> 2

ht_dt <-  ht_dt[idx_row,idx_col]

v_clust <- hclust(dist(ht_dt %>% as.matrix() %>% t()))
v_clust$order
ordered_labels <- names(ht_dt)[v_clust$order]

row_clust <- hclust(dist(ht_dt %>% as.matrix() ))
row_clust$order
row_labels <- rownames(ht_dt)[row_clust$order]



ht_dt[ht_dt > 0] <- 1
ht_dt[ht_dt < 0] <- -1

# Now Make a Heatmap of these data:
mat_logic <-  ht_dt[row_labels,ordered_labels]

mat_logic[mat_logic == 1] <- "Positive"
mat_logic[mat_logic == -1] <- "Negative"
mat_logic[mat_logic == 0] <- "NS"

cols <- structure(c("#3b5998","#990000", "white"), names = c("Positive","Negative","NS"))
library(ComplexHeatmap)

#colnames(mat_logic) <- gsub("\\.","-",colnames(mat_logic))

tax_dt <-  data.frame(tax_table(phy_sba_fil)@.Data)
tax_dt <-  tax_dt[match(rownames(mat_logic),tax_dt$fasta_id),]

split_rows <-  tax_dt$gene_desc
split_rows <-  rep("Other",length(split_rows))
split_rows[grep("Bai|Bile acid",tax_dt$gene_desc)] <- "Bai_genes"
split_rows[grep("3a|3-alpha|3alpha",tax_dt$gene_desc)] <- "3-a-HSDH"
split_rows[grep("3b|3-beta|3beta",tax_dt$gene_desc)] <- "3-b-HSDH"
split_rows[grep("7a|7-alpha|7alpha",tax_dt$gene_desc)] <- "7-a-HSDH"
split_rows[grep("7b|7-beta|7beta",tax_dt$gene_desc)] <- "7-b-HSDH"
split_rows[grep("Choloylglycine hydrolase|choloylglycine hydrolase|Bile salt",tax_dt$gene_desc)] <- "BSH"

tax_dt$gene_desc[grep("Other",split_rows)]


#max(mat_logic)
# col_grad <- colorRampPalette(c("blue", "white", "red"))(5)
# col_coef <- colorRamp2(c(-100,-50, 0, 50,100),col_grad) 

cols <-  c("Positive" = "#3b5998","Negative" = "#990000","NS" = "white")

sp_dt <-  data.frame(tax_table(phy_mic_fil)@.Data)
sp_dt$Name <-  rownames(sp_dt)

sp_dt <-  sp_dt[match(colnames(mat_logic),sp_dt$Name),]
sp_dt$Name == colnames(mat_logic)

split_cols <-  sp_dt$Phylum

colnames(mat_logic) <- sp_dt$Species

all_species <-  colnames(mat_logic)
ve303_sp <-  grep("VE303",all_species,value = T)
ve303_col <- ifelse(!all_species %in% ve303_sp,"black","darkred")

ht1 = Heatmap(as.matrix(mat_logic), name = "Coeff", column_title = NA, 
              #clustering_distance_rows = "euclidean",
              #clustering_distance_columns = "euclidean",
              cluster_column_slices = F,
              cluster_columns = F,
              cluster_rows = F,
              row_split = split_rows,
              column_split = split_cols,
              column_names_gp = gpar(col = ve303_col),
              column_title_rot = 60,
              column_names_rot = 60,
              row_gap = unit(2, "mm"),border = TRUE,
              row_title_gp = gpar(fontsize = 10,fontface = "bold"),
              row_title_rot = 0,
              #left_annotation = ha1,
              clustering_method_rows = "complete",row_names_side = "left", km=1, color_space = "LAB",
              col=cols,
              width=2, 
              row_names_gp = gpar(fontsize = 9),
              show_row_dend = F,
              show_column_dend = F,
              row_names_max_width = max_text_width(rownames(mat_logic),
                                                   gp = gpar(fontsize = 12)),
              column_names_max_height = max_text_width(colnames(mat_logic),
                                                       gp = gpar(fontsize = 12)),
              na_col="white")



pdf(paste(results_folder,"/",paste0("Heatmap_SBA_Microbes.pdf"),sep=""),height =12, width = 25)
draw(ht1)
dev.off()
