

# Input files
# for VE303 -- ve303.dnameta
# for Vanco -- ve303.vanco.meta
# for calpro -- ve303_cal. recode Study_Day_numeric as Day.in.treatment
# for cytokines -- cytokine.meta.sampling

plot_timeline <- function(X, data_type, result_folder){
  X$TRT.norm = factor(X$TRT.norm, levels = c("Placebo", "VE303 Low Dose", "VE303 High Dose"))
  
  ggplot(X, aes(x=as.numeric(Day.in.treatment), y=as.character(Subject.Number), color=TRT.norm)) + 
    scale_y_discrete() + 
    geom_rect(data = filter(X, TRT.norm == "Placebo"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
    geom_rect(data = filter(X, TRT.norm == "VE303 High Dose"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
    geom_rect(data = filter(X, TRT.norm == "VE303 Low Dose"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
    geom_point(size=2.5, show.legend = FALSE) + 
    theme_bw() + 
    labs(title= element_blank(), x="Day in Study", y="Subject", colour = "TRT.norm", shape=1) + 
    theme(axis.title = element_text(size=16), axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12), panel.grid.major.y = element_line(size = 1, color = "darkgrey"), panel.grid.major.x = element_line(size = .2, color = "darkgrey"), strip.text.y = element_text(angle = 0, size = 12)) + 
    facet_grid(TRT.norm ~ ., scales = "free_y", space = "free_y") + 
    scale_color_manual(values = coh.cols) + 
    scale_x_continuous(breaks = c( seq(-28, 168, by=14)))
  ggsave(filename = here(result_folder, paste(Sys.Date(), data_type,"Cohort_timeline-VE303-002-Unblind.png", sep = " ")), width = 14, height = 10)
  
  ##################################################################################################

  ggplot( X, aes(as.numeric(Day.in.treatment), y=as.character(Subject.Number), color= variable)) + 
    scale_y_discrete() + 
    geom_rect(data = filter(X, TRT.norm == "Placebo"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
    geom_rect(data = filter(X, TRT.norm == "VE303 High Dose" ), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
    geom_rect(data = filter(X, TRT.norm == "VE303 Low Dose" ), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
    
    geom_point(aes(x=as.numeric(Day.in.treatment), y=as.character(Subject.Number),color=TRT.norm), size=2.5, show.legend = FALSE) + 
    geom_point(aes(x=as.numeric(rec.day.in.treatment), y=as.character(Subject.Number),color=TRT.norm), shape = 13, size=3, show.legend = FALSE) + 
    geom_point(aes(x=as.numeric(abx.firsttreat.postdosing), y=as.character(Subject.Number),color=TRT.norm), shape = 2, size=2, show.legend = FALSE) + 
    
    theme_bw() + 
    labs(title= element_blank(), x="Day in Study", y="Subject", colour = "TRT.norm", shape=1) + 
    theme(axis.title = element_text(size=16), axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12), panel.grid.major.y = element_line(size = 1, color = "darkgrey"), panel.grid.major.x = element_line(size = .2, color = "darkgrey"), strip.text.y = element_text(angle = 0, size = 12)) + 
    facet_grid(TRT.norm.rec ~ ., scales = "free_y", space = "free_y") + 
    scale_color_manual(values = coh.cols) + 
    scale_x_continuous(breaks = c( seq(-28, 168, by=14)))
  ggsave(filename = here(result_folder, paste(Sys.Date(), data_type,"Cohort_abx_rec_timeline-VE303-002-Unblind.png", sep = " ")), width = 14, height = 10)

  ##################################################################################################
  X$day.in.treat.abx.end <- X$day.in.treat.abx.start + X$treat.day.len.abx
  X <- X %>% filter(Day.in.treatment < 16 & Screening_abx_check == "On-Abx Screening")
  
  
  ggplot( X, aes(as.numeric(Day.in.treatment), y=as.character(Subject.Number), color= variable)) + 
    scale_y_discrete() + 
    geom_rect(data = filter(X, TRT.norm == "Placebo"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
    geom_rect(data = filter(X, TRT.norm == "VE303 High Dose" ), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
    geom_rect(data = filter(X, TRT.norm == "VE303 Low Dose" ), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
    
    geom_point(aes(x=as.numeric(Day.in.treatment), y=as.character(Subject.Number),color=TRT.norm), size=2.5, show.legend = FALSE) + 
    geom_point(aes(x=as.numeric(day.in.treat.abx.start), y=as.character(Subject.Number),
                   color=TRT.norm), shape = 2, size=3, show.legend = FALSE) + 
    geom_point(aes(x=as.numeric(day.in.treat.abx.end), y=as.character(Subject.Number),
                   color=TRT.norm), shape = 2, size=3, show.legend = FALSE) + 
    
    theme_bw() + 
    labs(title= element_blank(), x="Day in Study", y="Subject", colour = "TRT.norm", shape=1) + 
    theme(axis.title = element_text(size=16), axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12), panel.grid.major.y = element_line(size = 1, color = "darkgrey"), panel.grid.major.x = element_line(size = .2, color = "darkgrey"), strip.text.y = element_text(angle = 0, size = 12)) + 
    facet_grid(TRT.norm.rec ~ ., scales = "free_y", space = "free_y") + 
    scale_color_manual(values = coh.cols) + 
    scale_x_continuous(breaks = c( seq(-28, 16, by=14)), limits = c(-28, 16))
  ggsave(filename = here(result_folder, paste(Sys.Date(), data_type,"Cohort_abx_rec_SCREEN_ABX-VE303-002-Unblind.png", sep = " ")), width = 14, height = 10)
  
    
} 
   
 