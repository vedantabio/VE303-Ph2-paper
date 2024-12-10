# Don't forget to comment on all the code
# Only reason is for you to know when you come back to this code after few months

#Load the  libraries
library(readxl)
library(tidyr)
library(dplyr)
library(gtools)
library(ggplot2)
library(ggthemes)
library(lemon)
library(here)
# Create a directories to store results from the analysis
mainDir <- here("Comb_SCFA_BA_plot")
dir.create(file.path(mainDir), showWarnings = TRUE,recursive = TRUE)
results_folder <- paste(mainDir,sep="/")

# Create a function that takes two input  : a) Metadata with treatment b) SCFA/BA abundance data
# And Outputs a combined and clean SCFA/BA abudance along with Treatment arms

#######################################BA##############
grp = "BA"
# Read SCFA abundance data
scfa_dt  <-  read.csv("../Processed_data/BA.csv")
unique(scfa_dt$Time.Point)

comb_dt <- scfa_dt
# Plotting stuffs
unique(comb_dt$Visit)

# Order the time labels
levels_time <- mixedsort(unique(comb_dt$Visit))
#levels_time <- c("Screening",levels_time[1:(length(levels_time)-1)])
comb_dt$Visit <-  factor(comb_dt$Visit, levels = levels_time)

# Order the Treatment levels 
comb_dt$TRT.norm <-  factor(as.character(comb_dt$TRT.norm),levels = c("Placebo" ,"VE303 Low Dose","VE303 High Dose"))


library(dplyr)
comb_dt  = comb_dt %>%
  mutate(.,BA.class = case_when(Analyte %in% c("Cholic Acid", "Chenodeoxycholic Acid") ~ "Primary BA",
                                
                                Analyte %in% c("Glycocholic Acid", "Glycochenodeoxycholic Acid", "Taurocholic Acid", "Taurochenodeoxycholic Acid") ~ "Primary BA",
                                
                                Analyte %in% c("Lithocholic Acid", "Deoxycholic Acid", "Isodeoxycholic Acid", "Ursodeoxycholic Acid") ~ "Secondary BA",
                                
                                Analyte %in% c("Alloiso_Isolithocholic Acid", "Dehydrolithocholic Acid", "Glycodeoxycholic Acid", "Glycolithocholic Acid", "Glycoursodeoxycholic Acid", "Taurodeoxycholic Acid", "Taurolithocholic Acid", "Tauroursodeoxycholic Acid") ~ "Secondary BA"))%>%
  data.frame()


comb_dt  <- comb_dt %>%
  dplyr::select(Sample.ID.metabolite, Subject.Number,BA.class,Visit.Name,TRT.norm,Est_Abun,rec.diagnosis,collect.pre.post.abx)%>%
  dplyr::group_by(Sample.ID.metabolite,Subject.Number,BA.class,Visit.Name,TRT.norm)%>%
  dplyr::mutate(Tot_abun = sum(Est_Abun))%>%
  dplyr::select(-c(Est_Abun))%>%
  unique()
comb_dt$log10_Abun <- log10(comb_dt$Tot_abun)

# Remove UNSCHEDULE points
#comb_dt <- comb_dt[!grepl("UNSCHED",comb_dt$Visit),]

# Process data for SCFA/BA 
# Keep necessary time points
# Prep for plotting
# Returns  clean dataframe
comb_dt$Visit <- as.character(comb_dt$Visit.Name)
sorted_time <- mixedsort(unique(comb_dt$Visit))
table(comb_dt$Visit)



# Sort time points in increasing order of pre and post
pre_labels <-  c("Screening",sorted_time[1:which(sorted_time == "Day 1")])

# Limit the time to Day 56
post_labels <- sorted_time[(which(sorted_time == "Day 1")+1):which(sorted_time == "Day 56")]

# Subset only Screening ,Day 1, Day 7 and Day 14
pre_labels <-  c("Screening","Day 1")
post_labels <- c("Day 7","Day 14","Day 56")

# Subset only Screening ,Day 1, Day 7 and Day 14
comb_dt <- comb_dt[comb_dt$Visit %in% c(pre_labels,post_labels),]
unique(comb_dt$Visit)

# Order the time labels
levels_time <- mixedsort(unique(comb_dt$Visit))
levels_time <- c("Screening",levels_time[1:(length(levels_time)-1)])
comb_dt$Visit <-  factor(comb_dt$Visit, levels = levels_time)

# Order the Treatment levels 
comb_dt$TRT.norm <-  factor(as.character(comb_dt$TRT.norm),levels = c("Placebo" ,"VE303 Low Dose","VE303 High Dose"))



# Now separate the groups into groups of recurrence / Not -recurrence
names(comb_dt)
comb_dt$rec.diagnosis <- as.character(comb_dt$rec.diagnosis)
# Replace NA with Good
#comb_dt$rec.diagnosis[is.na(comb_dt$rec.diagnosis)] <- "Good"
unique(comb_dt$rec.diagnosis)

tot_ba_dt <- comb_dt


# Now import scfa data:
#################################
grp = "SCFA"
# Read SCFA abundance data
scfa_dt  <-  read.csv("../Processed_data/SCFA.csv")
unique(scfa_dt$Time.Point)

comb_dt <- scfa_dt
# Plotting stuffs
unique(comb_dt$Visit)

# Order the time labels
levels_time <- mixedsort(unique(comb_dt$Visit))
#levels_time <- c("Screening",levels_time[1:(length(levels_time)-1)])
comb_dt$Visit <-  factor(comb_dt$Visit, levels = levels_time)

# Order the Treatment levels 
comb_dt$TRT.norm <-  factor(as.character(comb_dt$TRT.norm),levels = c("Placebo" ,"VE303 Low Dose","VE303 High Dose"))


library(dplyr)

comb_dt  <- comb_dt %>%
  dplyr::select(Sample.ID.metabolite, Subject.Number,Visit.Name,TRT.norm,Est_Abun,rec.diagnosis,collect.pre.post.abx)%>%
  dplyr::group_by(Sample.ID.metabolite,Subject.Number,Visit.Name,TRT.norm)%>%
  dplyr::mutate(Tot_abun = sum(Est_Abun))%>%
  dplyr::select(-c(Est_Abun))%>%
  unique()
comb_dt$log10_Abun <- log10(comb_dt$Tot_abun)

# Remove UNSCHEDULE points
#comb_dt <- comb_dt[!grepl("UNSCHED",comb_dt$Visit),]

# Process data for SCFA/BA 
# Keep necessary time points
# Prep for plotting
# Returns  clean dataframe
comb_dt$Visit <- as.character(comb_dt$Visit.Name)
sorted_time <- mixedsort(unique(comb_dt$Visit))
table(comb_dt$Visit)



# Sort time points in increasing order of pre and post
pre_labels <-  c("Screening",sorted_time[1:which(sorted_time == "Day 1")])

# Limit the time to Day 56
post_labels <- sorted_time[(which(sorted_time == "Day 1")+1):which(sorted_time == "Day 56")]

# Subset only Screening ,Day 1, Day 7 and Day 14
pre_labels <-  c("Screening","Day 1")
post_labels <- c("Day 7","Day 14","Day 56")

# Subset only Screening ,Day 1, Day 7 and Day 14
comb_dt <- comb_dt[comb_dt$Visit %in% c(pre_labels,post_labels),]
unique(comb_dt$Visit)

# Order the time labels
levels_time <- mixedsort(unique(comb_dt$Visit))
levels_time <- c("Screening",levels_time[1:(length(levels_time)-1)])
comb_dt$Visit <-  factor(comb_dt$Visit, levels = levels_time)

# Order the Treatment levels 
comb_dt$TRT.norm <-  factor(as.character(comb_dt$TRT.norm),levels = c("Placebo" ,"VE303 Low Dose","VE303 High Dose"))



# Now separate the groups into groups of recurrence / Not -recurrence
names(comb_dt)
comb_dt$rec.diagnosis <- as.character(comb_dt$rec.diagnosis)
# Replace NA with Good
#comb_dt$rec.diagnosis[is.na(comb_dt$rec.diagnosis)] <- "Good"
unique(comb_dt$rec.diagnosis)


tot_scfa_dt <-  comb_dt
tot_scfa_dt$BA.class = "SCFA"
# Now combine the two 
comb_final_dt <- rbind(tot_ba_dt,tot_scfa_dt)

names(tot_ba_dt)
names(tot_scfa_dt)
head(comb_final_dt)
# Plot points on top of the boxplot
library(ggplot2)
library(ggthemes)
library(lemon)
plot_dt_sel <-  comb_final_dt
coh.cols <- c("Placebo" = "#0487e3", "VE303 High Dose" = "#dc2d00", "VE303 Low Dose" = "#5b0076")
# Good Vs Bad by removing Post Abx samples
sample_cols <-  c('#377eb8','#e41a1c')
names(sample_cols) <- c(  "Good","Bad")
# Order groups for Recurrence
plot_dt_sel$rec.diagnosis <-factor(as.character(plot_dt_sel$rec.diagnosis),levels = rev(names(sample_cols)))

plot_dt_sub <-  plot_dt_sel[plot_dt_sel$collect.pre.post.abx == "Pre.abx",]
plot_dt_sub$BA.class <-  factor(plot_dt_sub$BA.class,levels = c("Primary BA","Secondary BA","SCFA"))
head(plot_dt_sub)
p_scfa <- ggplot(plot_dt_sub, aes(x = Visit, y = log10_Abun)) +
  geom_boxplot(aes(y = log10_Abun, fill = rec.diagnosis, color = rec.diagnosis), alpha = 0.3, outlier.colour = NA,
               position=position_dodge2(width=0.75,preserve = "single"))+
  geom_point(aes(y = log10_Abun, color = rec.diagnosis),
             position=position_dodge(width=0.75,preserve = "total"),
             alpha = 0.3,size = 2)+
  scale_fill_manual(name = "Recurrence",values = sample_cols )+
  scale_color_manual(name = "Recurrence",values = sample_cols )+
  facet_rep_wrap(BA.class~TRT.norm , scales = "free",ncol = 3)+
  theme_base()+
  theme(legend.position = "top")+
  theme(axis.text.x=element_text(angle=45, hjust=1))+
  xlab("Day")+
  ylab("Log10(Abundance)")

pdf(paste(results_folder,paste0("Total_SCFA_BA_Pre_Abx",'.pdf'),sep="/"),height =14 , width = 10, useDingbats = FALSE)
print(p_scfa)
dev.off()


