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
mainDir <- here("SBA_Analysis")
dir.create(file.path(mainDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,sep="/")

# Read Processed data:
phy_sba <-  readRDS(here("../Processed_data/shortbred/phy_sba.rds"))
phy_sba <- subset_samples(phy_sba,Visit.Name %in% c("Screening","Day 1","Day 7","Day 14","Day 28","Day 56") )


# Compare Placebo Vs LD Vs HD

#phy_lme <- subset_samples(phy_sba, Visit.Name %in% c("Day 7","Day 14","Day 28") )
phy_lme <- subset_samples(phy_sba, Visit.Name %in% c("Day 7","Day 14","Day 28", "Day 56") )
#phy_lme <- subset_samples(phy_lme, TRT.norm %in% c("Placebo","VE303 Low Dose") )
#phy_lme <- subset_samples(phy_lme, collect.pre.post.abx %in% c("Pre.abx") )

phy_lme <- prune_taxa(taxa_sums(phy_lme) >0 , phy_lme)

# Keep prevalent species
library(microbiome)
prev_taxa <-  stack(prevalence(phy_lme, detection=0, sort=TRUE, count=TRUE))
names(prev_taxa) <- c("Prevalence","Taxa")

avg_abun <- stack(apply(data.frame(otu_table(phy_lme)),1,mean))
names(avg_abun)  <- c("Avg_abun","Taxa")

prevdf <-  prev_taxa %>% inner_join(avg_abun)

prev_frac<-0.1
prev_cutoff <- prev_frac*nsamples(phy_lme) # Cut-off
#ab_cutoff <- 10^-3 # Cut-off

# Prevalence
prevdf_fil <- prevdf[prevdf$Prevalence >= prev_cutoff, ]
# Abundance
#prevdf_fil <- prevdf_fil[prevdf_fil$Avg_abun >= ab_cutoff, ]


phy_lme_fil <- prune_taxa(as.character(prevdf_fil$Taxa),phy_lme)

# Transform
sm_dt <- data.frame(sample_data(phy_lme_fil))

pat_dt <-  unique(sm_dt[,c("Subject.Number","TRT.norm","rec.diagnosis.Wk8")])
table(pat_dt$rec.diagnosis.Wk8,pat_dt$TRT.norm)

phy_melt <-  psmelt(phy_lme_fil)
names(phy_melt)

phy_melt$Subject.Number <-  factor(phy_melt$Subject.Number)
unique(phy_melt$Subject.Number)

phy_melt$Timepoint  <-  factor(phy_melt$Visit.Name)

names(phy_melt)

unique(phy_melt$Subject.Number)
unique(phy_melt$Timepoint)

result_lme <-  list()
all_taxa <- unique(phy_melt$OTU)
taxa <- all_taxa[1]

trt <- unique(phy_melt$TRT.norm)[1]

list_trt <- list()

taxa_trt <- phy_melt
#[phy_melt$TRT.norm %in% trt,]
result_lme <-  list()
for (taxa in all_taxa){
  print(taxa)
  sum_mod_dt <- tryCatch({
    taxa_sel <- taxa_trt[taxa_trt$OTU %in% taxa,]
    taxa_sel$t_Abun <- log(taxa_sel$Abundance + 0.01)
    taxa_sel$Time <- as.numeric(gsub("Day ","",taxa_sel$Timepoint))
    taxa_sel$TRT.norm <- factor(taxa_sel$TRT.norm)
    library("nlme")
    mod_var <-  "TRT.norm + Time "
    fml <- as.formula( paste( "t_Abun", "~", paste(mod_var, collapse="+") ) )
    mod_bac <- lme(fml,random = ~ 1|Subject.Number, taxa_sel)
    sum_mod <-  summary(mod_bac)
    sum_mod$tTable
    #print(anova.lme(mod_bac))
    library(emmeans)
    #emmeans(mod_bac)
    em_mod_bac <- emmeans(mod_bac,revpairwise ~ TRT.norm ,adjust = "none")
    #summary(em_mod_bac)
    #pairs(em_mod_bac)
    sum_mod_dt <- data.frame(summary(em_mod_bac)$contrasts)
    sum_mod_dt$taxa <- taxa
    # #sum_mod_dt$trt <- trt
    #sum_mod_dt$Var <-  rownames(sum_mod_dt)
    # sum_mod_dt <- data.frame(sum_mod$tTable)
    # sum_mod_dt$taxa <- taxa
    # #sum_mod_dt$trt <- trt
    # sum_mod_dt$Var <-  rownames(sum_mod_dt)
    #sum_mod_dt
    
    result_lme[[taxa]] <-  sum_mod_dt
    
    
    
    
    
  }, error = function(err){
    sum_mod_dt <- NA
    sum_mod_dt
  }
  )
  result_lme[[taxa]] <-  sum_mod_dt
}

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[6] <- "pval"

result_lme_dt <- result_lme_dt %>%
               # filter(contrast %in% c("Placebo - VE303 High Dose"))%>%
                 group_by(contrast)%>%
               mutate(pval.adj = p.adjust(pval, method='BH'))

final_lme_dt <-  result_lme_dt
write.csv(final_lme_dt,paste0(results_folder,"/LME_Dosed_Group_Comparison.csv"))

sig_lme_dt <-  final_lme_dt[final_lme_dt$pval.adj  < 0.1,]

unique(sig_lme_dt$taxa)

library(tidyverse)

phy_sel <-  prune_taxa( unique(sig_lme_dt$taxa), phy_sba) 

phy_ht <-  phy_sel
# Make a heatmap of SBA genes in Phase 2:
# Relative abundance:
mat <- data.frame(otu_table(phy_ht))
met <- data.frame(sample_data(phy_ht))

met$Visit.Name <- factor(met$Visit.Name, levels = c("Screening","Day 1","Day 7","Day 14","Day 28","Day 56"))
met <-  met[order(met$TRT.norm,met$Visit.Name,met$Subject.Number),]
mat <-  mat[,match(rownames(met),colnames(mat))]
colnames(mat) == rownames(met)


library(gtools)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(yingtools2)
library(data.table)
library(hues)
set.seed(1057)

library(scico)
library(ComplexHeatmap)

col_mat <- c("white",viridis::viridis(10))
splitcols <- factor(paste0(met$TRT.norm))
#"#440154FF""#482878FF"

coh.cols <- c("Placebo" = "#0487e3", "VE303 High Dose" = "#dc2d00", "VE303 Low Dose" = "#5b0076")


ha_column = HeatmapAnnotation(Cohort = met$TRT.norm,
                              Time = met$Visit.Name,
                              show_annotation_name = F, 
                              col=list(Cohort = coh.cols,
                                       Time = c("Screening" = "white","Day 1" = "grey50","Day 7" = "lightblue","Day 14" = "blue", "Day 28" = "darkblue",
                                                "Day 56" = "purple")))


library(RColorBrewer)
# status_col <-  brewer.pal(length(unique(met$Visit)), "Set1")
# names(status_col) <- unique(unique(met$Visit))

library(ComplexHeatmap)
library(tidyverse)
library(gtools)
library(RColorBrewer)
library(circlize)
col_bar = colorRamp2(c(0, 10), c("white", "blue"))
# Taxa annotations
library(data.table)
library(hues)
set.seed(1057)
col_b <- c("0" = "grey","1"= "blue")

# Order columns based on CFU
#col_time <- scico(2,palette = "devon")
#col_trt <- c("0 Day" = "grey50","2 Wk"= "#cb2314","2 Mo" = "#fad510","6 Mo" = "purple")
#splitcols <-factor(as.character(met$Visit),levels = names(col_trt))
#splitcols <-factor(as.character(met$PID))
#splitcols <- factor(splitcols,levels = c("preTreatment","Treatment"))
tax_dt <-  data.frame(tax_table(phy_ht)@.Data)
tax_dt <-  tax_dt[match(rownames(mat),tax_dt$fasta_id),]

#rownames(mat) <- paste0(tax_dt$fasta_id,"_",tax_dt$gene_desc)

#mat_sel <-  mat[grep("stx|aaiC|aagR|eae|aatA|bfpA|elt|estI",rownames(mat)),]
mat_sel <- mat
min(mat_sel[mat_sel>0])
split_rows <-  factor(tax_dt$gene_desc)


#Make row annotation:
coeff_dt <-  sig_lme_dt[,c("taxa","contrast","estimate")]
coeff_dt_w <-  coeff_dt %>%
  pivot_wider(names_from = contrast,values_from = estimate)%>%
  data.frame()
coeff_dt_w <-  coeff_dt_w[match(rownames(mat),coeff_dt_w$taxa),]
#names(coeff_dt_w) <- c("taxa","HD_Placebo","LD_Placebo","HD_LD")
names(coeff_dt_w) <- c("taxa","HD_Placebo","LD_Placebo")

col_grad <- colorRampPalette(c("blue", "white", "red"))(5)
col_coef <- colorRamp2(c(-2,-1, 0, 1,2),col_grad) 

ha_left =  rowAnnotation(HD_Vs_placebo = coeff_dt_w$HD_Placebo,
                         LD_Vs_placebo = coeff_dt_w$LD_Placebo,
                         col = list(HD_Vs_placebo = col_coef,
                                    LD_Vs_placebo = col_coef),
                         show_legend = c(TRUE,F,F),
                         show_annotation_name = c(TRUE,TRUE,T),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         annotation_label = c("HD Vs P","LD Vs P"),
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 10,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 10,fontface = "bold" ),
                                                        at =c(-2,-1, 0, 1,2),
                                                        legend_height = unit(3, "cm")))

ht  =  Heatmap(log(mat_sel+ 0.001),name = "Log(Abun + 0.01)",
               column_split = splitcols,
               row_split =  split_rows,
               top_annotation = ha_column,
               left_annotation = ha_left,
               col = col_mat,
               row_names_side = "left",
               row_gap = unit(2, "mm"),
               show_column_names = F,
               row_title_gp = gpar(fontsize = 10,fontface = "bold"),
               row_title_rot = 0,
               column_title_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 0,
               cluster_rows = T,
               cluster_row_slices = T,
               show_parent_dend_line = F,
               show_row_dend = F,
               cluster_columns = F,
               border = T,
               heatmap_legend_param = list(title = "Log(RPKM + 0.01)", 
                                           legend_width = unit(10, "cm"),
                                           legend_direction = "horizontal",
                                           title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                           labels_gp = gpar(fontsize = 12,fontface = "bold" )),
               row_names_max_width = max_text_width(rownames(mat_sel),
                                                    gp = gpar(fontsize = 12)),
               column_names_max_height = max_text_width(colnames(mat_sel),
                                                        gp = gpar(fontsize = 12)))
pdf(paste0(results_folder,"/Heatmap_SBA_Sig.pdf"),height = 5, width = 20,useDingbats = F)
draw(ht,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
#draw(ht)
dev.off()


# Now compare R vs NR in dosed subjects:

phy_lme <- subset_samples(phy_sba,! TRT.norm %in% c("Placebo") )
#phy_lme <- subset_samples(phy_lme, Visit.Name %in% c("Day 7","Day 14","Day 28", "Day 56") )
phy_lme <- subset_samples(phy_lme, Visit.Name %in% c("Day 7","Day 14","Day 28") )
#phy_lme <- subset_samples(phy_lme, TRT.norm %in% c("Placebo","VE303 Low Dose") )
#phy_lme <- subset_samples(phy_lme, collect.pre.post.abx %in% c("Pre.abx") )

phy_lme <- prune_taxa(taxa_sums(phy_lme) >0 , phy_lme)

# Keep prevalent species
library(microbiome)
prev_taxa <-  stack(prevalence(phy_lme, detection=0, sort=TRUE, count=TRUE))
names(prev_taxa) <- c("Prevalence","Taxa")

avg_abun <- stack(apply(data.frame(otu_table(phy_lme)),1,mean))
names(avg_abun)  <- c("Avg_abun","Taxa")

prevdf <-  prev_taxa %>% inner_join(avg_abun)

prev_frac<-0.2
prev_cutoff <- prev_frac*nsamples(phy_lme) # Cut-off
#ab_cutoff <- 10^-3 # Cut-off

# Prevalence
prevdf_fil <- prevdf[prevdf$Prevalence >= prev_cutoff, ]
# Abundance
#prevdf_fil <- prevdf_fil[prevdf_fil$Avg_abun >= ab_cutoff, ]


phy_lme_fil <- prune_taxa(as.character(prevdf_fil$Taxa),phy_lme)

# Transform
sm_dt <- data.frame(sample_data(phy_lme_fil))
#sm_dt <- data.frame(sample_data(phy_sba))
pat_dt <-  unique(sm_dt[,c("Subject.Number","TRT.norm","rec.diagnosis.Wk8")])
table(pat_dt$rec.diagnosis.Wk8,pat_dt$TRT.norm)

phy_melt <-  psmelt(phy_lme_fil)
names(phy_melt)

phy_melt$Subject.Number <-  factor(phy_melt$Subject.Number)
unique(phy_melt$Subject.Number)

phy_melt$Timepoint  <-  factor(phy_melt$Visit.Name)

names(phy_melt)

unique(phy_melt$Subject.Number)
unique(phy_melt$Timepoint)

result_lme <-  list()
all_taxa <- unique(phy_melt$OTU)
taxa <- all_taxa[1]

trt <- unique(phy_melt$TRT.norm)[1]

list_trt <- list()

taxa_trt <- phy_melt
#[phy_melt$TRT.norm %in% trt,]
result_lme <-  list()
for (taxa in all_taxa){
  print(taxa)
  sum_mod_dt <- tryCatch({
    taxa_sel <- taxa_trt[taxa_trt$OTU %in% taxa,]
    taxa_sel$t_Abun <- log(taxa_sel$Abundance + 0.01)
    taxa_sel$Time <- as.numeric(gsub("Day ","",taxa_sel$Timepoint))
    taxa_sel$TRT.norm <- factor(taxa_sel$TRT.norm)
    taxa_sel$Status <- factor(taxa_sel$rec.diagnosis.Wk8, levels = c("Bad","Good"))
    levels(taxa_sel$Status) <-  c("NR","R")
    library("nlme")
    mod_var <-  "Status + Time "
    fml <- as.formula( paste( "t_Abun", "~", paste(mod_var, collapse="+") ) )
    mod_bac <- lme(fml,random = ~ 1|Subject.Number, taxa_sel)
    sum_mod <-  summary(mod_bac)
    #sum_mod$tTable
    #print(anova.lme(mod_bac))
    library(emmeans)
    #emmeans(mod_bac)
    em_mod_bac <- emmeans(mod_bac,revpairwise ~ Status ,adjust = "none")
    #summary(em_mod_bac)
    #pairs(em_mod_bac)
    sum_mod_dt <- data.frame(summary(em_mod_bac)$contrasts)
    sum_mod_dt$taxa <- taxa
    result_lme[[taxa]] <-  sum_mod_dt
    
    
    
    
    
  }, error = function(err){
    sum_mod_dt <- NA
    sum_mod_dt
  }
  )
  result_lme[[taxa]] <-  sum_mod_dt
}

result_lme_dt <- do.call("rbind", result_lme)
names(result_lme_dt)[6] <- "pval"

result_lme_dt <- result_lme_dt %>%
  # filter(contrast %in% c("Placebo - VE303 High Dose"))%>%
  group_by(contrast)%>%
  mutate(pval.adj = p.adjust(pval, method='BH'))

final_lme_dt <-  result_lme_dt
write.csv(final_lme_dt,paste0(results_folder,"/LME_Recurrence_Comparison.csv"))
sig_lme_dt <-  final_lme_dt[final_lme_dt$pval.adj  < 0.2,]

library(tidyverse)
phy_sel <- subset_samples(phy_sba,! TRT.norm %in% c("Placebo") )
phy_sel <- subset_samples(phy_sel, Visit.Name %in% c("Day 7","Day 14","Day 28") )

phy_sel <-  prune_taxa( unique(sig_lme_dt$taxa), phy_sel) 

phy_ht <-  phy_sel
# Make a heatmap of SBA genes in Phase 2:
# Relative abundance:
mat <- data.frame(otu_table(phy_ht))
met <- data.frame(sample_data(phy_ht))

met$Visit.Name <- factor(met$Visit.Name, levels = c("Screening","Day 1","Day 7","Day 14","Day 28","Day 56"))
met <-  met[order(met$rec.diagnosis.Wk8,met$Visit.Name,met$Subject.Number),]
mat <-  mat[,match(rownames(met),colnames(mat))]
colnames(mat) == rownames(met)


library(gtools)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(yingtools2)
library(data.table)
library(hues)
set.seed(1057)

library(scico)
library(ComplexHeatmap)

col_mat <- c("white",viridis::viridis(10))
splitcols <- factor(paste0(met$rec.diagnosis.Wk8))
levels(splitcols) <- c("R","NR")
#"#440154FF""#482878FF"

coh.cols <- c("Placebo" = "#0487e3", "VE303 High Dose" = "#dc2d00", "VE303 Low Dose" = "#5b0076")


ha_column = HeatmapAnnotation(Cohort = met$TRT.norm,
                              Time = met$Visit.Name,
                              show_annotation_name = F, 
                              col=list(Cohort = coh.cols,
                                       Time = c("Screening" = "white","Day 1" = "grey50","Day 7" = "lightblue","Day 14" = "blue", "Day 28" = "darkblue",
                                                "Day 56" = "purple")))


library(RColorBrewer)
# status_col <-  brewer.pal(length(unique(met$Visit)), "Set1")
# names(status_col) <- unique(unique(met$Visit))

library(ComplexHeatmap)
library(tidyverse)
library(gtools)
library(RColorBrewer)
library(circlize)
col_bar = colorRamp2(c(0, 10), c("white", "blue"))
# Taxa annotations
library(data.table)
library(hues)
set.seed(1057)
col_b <- c("0" = "grey","1"= "blue")

# Order columns based on CFU
#col_time <- scico(2,palette = "devon")
#col_trt <- c("0 Day" = "grey50","2 Wk"= "#cb2314","2 Mo" = "#fad510","6 Mo" = "purple")
#splitcols <-factor(as.character(met$Visit),levels = names(col_trt))
#splitcols <-factor(as.character(met$PID))
#splitcols <- factor(splitcols,levels = c("preTreatment","Treatment"))
tax_dt <-  data.frame(tax_table(phy_ht)@.Data)
tax_dt <-  tax_dt[match(rownames(mat),tax_dt$fasta_id),]

#rownames(mat) <- paste0(tax_dt$fasta_id,"_",tax_dt$gene_desc)

#mat_sel <-  mat[grep("stx|aaiC|aagR|eae|aatA|bfpA|elt|estI",rownames(mat)),]
mat_sel <- mat
min(mat_sel[mat_sel>0])
split_rows <-  factor(tax_dt$gene_desc)


#Make row annotation:
coeff_dt <-  sig_lme_dt[,c("taxa","contrast","estimate")]
coeff_dt_w <-  coeff_dt %>%
  pivot_wider(names_from = contrast,values_from = estimate)%>%
  data.frame()
coeff_dt_w <-  coeff_dt_w[match(rownames(mat),coeff_dt_w$taxa),] %>% as.data.frame()
#names(coeff_dt_w) <- c("taxa","HD_Placebo","LD_Placebo","HD_LD")
names(coeff_dt_w) <- c("taxa","R_NR")

col_grad <- colorRampPalette(c("blue", "white", "red"))(5)
col_coef <- colorRamp2(c(-2,-1, 0, 1,2),col_grad) 

ha_left =  rowAnnotation(R_Vs_NR = coeff_dt_w$R_NR,
                         col = list( R_Vs_NR = col_coef),
                         show_legend = c(TRUE),
                         show_annotation_name = c(TRUE),
                         annotation_name_gp = gpar(fontface = "bold" ),
                         annotation_name_rot = 45,
                         annotation_label = c("R Vs NR"),
                         annotation_legend_param = list(title = "Coeff",
                                                        title_gp = gpar(fontsize = 10,fontface = "bold" ),
                                                        labels_gp = gpar(fontsize = 10,fontface = "bold" ),
                                                        at =c(-2,-1, 0, 1,2),
                                                        legend_height = unit(3, "cm")))
rownames(mat_sel) <- gsub(".*\\|","",rownames(mat_sel))
colnames(mat_sel) <- met$Subject.Number
ht  =  Heatmap(log(mat_sel+ 0.001),name = "Log(Abun + 0.01)",
               column_split = splitcols,
               row_split =  split_rows,
               top_annotation = ha_column,
               left_annotation = ha_left,
               col = col_mat,
               row_names_side = "left",
               row_gap = unit(2, "mm"),
               show_column_names = T,
               row_title_gp = gpar(fontsize = 10,fontface = "bold"),
               row_title_rot = 0,
               column_title_gp = gpar(fontsize = 10,fontface = "bold"),
               column_title_rot = 0,
               cluster_rows = T,
               cluster_row_slices = T,
               show_parent_dend_line = F,
               show_row_dend = F,
               cluster_columns = F,
               border = T,
               heatmap_legend_param = list(title = "Log(RPKM + 0.01)", 
                                           legend_width = unit(10, "cm"),
                                           legend_direction = "horizontal",
                                           title_gp = gpar(fontsize = 14,fontface = "bold" ),
                                           labels_gp = gpar(fontsize = 12,fontface = "bold" )),
               row_names_max_width = max_text_width(rownames(mat_sel),
                                                    gp = gpar(fontsize = 12)),
               column_names_max_height = max_text_width(colnames(mat_sel),
                                                        gp = gpar(fontsize = 12)))
pdf(paste0(results_folder,"/Heatmap_SBA_Sig_R_vs_NR.pdf"),height = 4, width = 25,useDingbats = F)
draw(ht,heatmap_legend_side = c("top"),annotation_legend_side="right",legend_grouping = "original")
#draw(ht)
dev.off()

