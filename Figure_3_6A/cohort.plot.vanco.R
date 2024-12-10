

ve303.dnameta$TRT.norm = factor(ve303.dnameta$TRT.norm, levels = c("Placebo", "VE303 Low Dose", "VE303 High Dose"))
ve303.vanco.meta$TRT.norm = factor(ve303.vanco.meta$TRT.norm, levels = c("Placebo", "VE303 Low Dose", "VE303 High Dose"))
c1$TRT.norm = factor(c1$TRT.norm, levels = c("Placebo", "VE303 Low Dose", "VE303 High Dose"))


## Vanco sampling timeline
##################################################################################################
##################################################################################################

ggplot(ve303.vanco.meta, aes(x=as.numeric(Day.in.treatment), y=as.character(Subject.Number), color=TRT.norm)) + 
  scale_y_discrete() + 
  geom_rect(data = filter(ve303.vanco.meta, TRT.norm == "Placebo"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
  geom_rect(data = filter(ve303.vanco.meta, TRT.norm == "VE303 High Dose"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
  geom_rect(data = filter(ve303.vanco.meta, TRT.norm == "VE303 Low Dose"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
  geom_point(size=2.5, show.legend = FALSE) + 
  theme_bw() + 
  labs(title= element_blank(), x="Day in Study", y="Subject", colour = "TRT.norm", shape=1) + 
  theme(axis.title = element_text(size=16), axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12), panel.grid.major.y = element_line(size = 1, color = "darkgrey"), panel.grid.major.x = element_line(size = .2, color = "darkgrey"), strip.text.y = element_text(angle = 0, size = 12)) + 
  facet_grid(TRT.norm ~ ., scales = "free_y", space = "free_y") + 
  scale_color_manual(values = coh.cols) + 
  scale_x_continuous(breaks = c( seq(-28, 168, by=14)))
ggsave(filename = here(results, paste(Sys.Date(), "Cohort_timeline-Vanco-VE303-002-Unblind.png", sep = " ")), width = 14, height = 10)
 
##################################################################################################

ggplot( ve303.vanco.meta, aes(as.numeric(Day.in.treatment), y=as.character(Subject.Number), color= variable)) + 
  scale_y_discrete() + 
  geom_rect(data = filter(ve303.vanco.meta, TRT.norm == "Placebo"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
  geom_rect(data = filter(ve303.vanco.meta, TRT.norm == "VE303 High Dose" ), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
  geom_rect(data = filter(ve303.vanco.meta, TRT.norm == "VE303 Low Dose" ), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
  
  geom_point(aes(x=as.numeric(Day.in.treatment), y=as.character(Subject.Number),color=TRT.norm), size=2.5, show.legend = FALSE) + 
  geom_point(aes(x=as.numeric(rec.day.in.treatment), y=as.character(Subject.Number),color=TRT.norm), shape = 13, size=3, show.legend = FALSE) + 
  geom_point(aes(x=as.numeric(abx.firsttreat.postdosing), y=as.character(Subject.Number),color=TRT.norm), shape = 2, size=2, show.legend = FALSE) + 
  
  theme_bw() + 
  labs(title= element_blank(), x="Day in Study", y="Subject", colour = "TRT.norm", shape=1) + 
  theme(axis.title = element_text(size=16), axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12), panel.grid.major.y = element_line(size = 1, color = "darkgrey"), panel.grid.major.x = element_line(size = .2, color = "darkgrey"), strip.text.y = element_text(angle = 0, size = 12)) + 
  facet_grid(TRT.norm.rec ~ ., scales = "free_y", space = "free_y") + 
  scale_color_manual(values = coh.cols) + 
  scale_x_continuous(breaks = c( seq(-28, 168, by=14)))
ggsave(filename = here(results, paste(Sys.Date(), "Cohort_recurrence_abxrec_timeline-Vanco-VE303-002-Unblind.png", sep = " ")), width = 14, height = 10)

##################################################################################################

#Calprotectin timeline
##################################################################################################

ggplot( c1, aes(as.numeric(Study_day_numeric), y=as.character(Subject.Number), color= variable)) + 
  scale_y_discrete() + 
  geom_rect(data = filter(c1, TRT.norm == "Placebo"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
  geom_rect(data = filter(c1, TRT.norm == "VE303 High Dose" ), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
  geom_rect(data = filter(c1, TRT.norm == "VE303 Low Dose" ), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
  
  geom_point(aes(x=as.numeric(Study_day_numeric), y=as.character(Subject.Number),color=TRT.norm), size=2.5, show.legend = FALSE) + 
  geom_point(aes(x=as.numeric(rec.day.in.treatment), y=as.character(Subject.Number),color=TRT.norm), shape = 13, size=3, show.legend = FALSE) + 
  geom_point(aes(x=as.numeric(abx.firsttreat.postdosing), y=as.character(Subject.Number),color=TRT.norm), shape = 2, size=2, show.legend = FALSE) + 
  
  theme_bw() + 
  labs(title= element_blank(), x="Day in Study", y="Subject", colour = "TRT.norm", shape=1) + 
  theme(axis.title = element_text(size=16), axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12), panel.grid.major.y = element_line(size = 1, color = "darkgrey"), panel.grid.major.x = element_line(size = .2, color = "darkgrey"), strip.text.y = element_text(angle = 0, size = 12)) + 
  facet_grid(TRT.norm.rec ~ ., scales = "free_y", space = "free_y") + 
  scale_color_manual(values = coh.cols) + 
  scale_x_continuous(breaks = c( seq(-28, 168, by=14)))
ggsave(filename = here(results, paste(Sys.Date(), "Cohort_recurrence_abxrec_timeline-Calpro-VE303-002-Unblind.png", sep = " ")), width = 14, height = 10)

##################################################################################################
##################################################################################################

## Metagenomics sequencing sampling timeline

##################################################################################################
##################################################################################################


ggplot(ve303.dnameta, aes(as.numeric(Day.in.treatment), y=as.character(Subject.Number), color= variable)) + 
  scale_y_discrete() + 
  geom_rect(data = filter(ve303.dnameta, TRT.norm == "Placebo"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
  geom_rect(data = filter(ve303.dnameta, TRT.norm == "VE303 High Dose" ), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
  geom_rect(data = filter(ve303.dnameta, TRT.norm == "VE303 Low Dose" ), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
  
  geom_point(aes(x=as.numeric(Day.in.treatment), y=as.character(Subject.Number),color=TRT.norm), size=2.5, show.legend = FALSE) + 
  geom_point(aes(x=as.numeric(rec.day.in.treatment), y=as.character(Subject.Number),color=TRT.norm), shape = 13, size=3, show.legend = FALSE) + 
  geom_point(aes(x=as.numeric(abx.firsttreat.postdosing), y=as.character(Subject.Number),color=TRT.norm), shape = 2, size=2, show.legend = FALSE) + 
  
  theme_bw() + 
  labs(title= element_blank(), x="Day in Study", y="Subject", colour = "TRT.norm", shape=1) + 
  theme(axis.title = element_text(size=16), axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12), panel.grid.major.y = element_line(size = 1, color = "darkgrey"), panel.grid.major.x = element_line(size = .2, color = "darkgrey"), strip.text.y = element_text(angle = 0, size = 12)) + 
  facet_grid(TRT.norm.rec ~ ., scales = "free_y", space = "free_y") + 
  scale_color_manual(values = coh.cols) + 
  scale_x_continuous(breaks = c( seq(-28, 168, by=14)))
ggsave(filename = here(results, paste(Sys.Date(), "Cohort_recurrence_abxrec_timeline-VE303-002-Unblind.png", sep = " ")), width = 14, height = 10)

# 
# ##################################################################################################
# 
# ggplot(ve303.dnameta, aes(x=as.numeric(Day.in.treatment), y=as.character(Subject.Number), color=TRT.norm)) + 
#   scale_y_discrete() + 
#   geom_rect(data = filter(ve303.dnameta, TRT.norm == "Placebo"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
#   geom_rect(data = filter(ve303.dnameta, TRT.norm == "VE303 High Dose"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
#   geom_rect(data = filter(ve303.dnameta, TRT.norm == "VE303 Low Dose"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
#   geom_point(size=2.5, show.legend = FALSE) + 
#   theme_bw() + 
#   labs(title= element_blank(), x="Day in Study", y="Subject", colour = "TRT.norm", shape=1) + 
#   theme(axis.title = element_text(size=16), axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12), panel.grid.major.y = element_line(size = 1, color = "darkgrey"), panel.grid.major.x = element_line(size = .2, color = "darkgrey"), strip.text.y = element_text(angle = 0, size = 12)) + 
#   facet_grid(TRT.norm ~ ., scales = "free_y", space = "free_y") + 
#   scale_color_manual(values = coh.cols) + 
#   scale_x_continuous(breaks = c( seq(-28, 168, by=14)))
# ggsave(filename = here(results, paste(Sys.Date(), "Cohort_timeline-VE303-002-Unblind.png", sep = " ")), width = 14, height = 10)
# 
# ##################################################################################################
# 
# ggplot(ve303.dnameta, aes(x=as.numeric(Day.in.treatment), y=as.character(Subject.Number), color=TRT.norm)) + 
#   scale_y_discrete() + 
#   geom_rect(data = filter(ve303.dnameta, TRT.norm == "Placebo"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
#   geom_rect(data = filter(ve303.dnameta, TRT.norm == "VE303 High Dose" ), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
#   geom_rect(data = filter(ve303.dnameta, TRT.norm == "VE303 Low Dose" ), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), alpha = 0.05, fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
#   geom_point(size=2.5, show.legend = FALSE) + 
#   theme_bw() + 
#   labs(title= element_blank(), x="Day in Study", y="Subject", colour = "TRT.norm", shape=1) + 
#   theme(axis.title = element_text(size=16), axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12), panel.grid.major.y = element_line(size = 1, color = "darkgrey"), panel.grid.major.x = element_line(size = .2, color = "darkgrey"), strip.text.y = element_text(angle = 0, size = 12)) + 
#   facet_grid(TRT.norm.rec ~ ., scales = "free_y", space = "free_y") + 
#   scale_color_manual(values = coh.cols) + 
#   scale_x_continuous(breaks = c( seq(-28, 168, by=14)))
# ggsave(filename = here(results, paste(Sys.Date(), "Cohort_recurrence_timeline-VE303-002-Unblind.png", sep = " ")), width = 14, height = 10)
# 
# ##################################################################################################


# ggplot(ve303.dnameta, aes(x=as.numeric(Day.in.treatment), y=as.character(Subject.Number))) + 
#   scale_y_discrete() + 
#   geom_rect(data = filter(ve303.dnameta, TRT.norm == "Placebo"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
#   geom_rect(data = filter(ve303.dnameta, TRT.norm == "VE303 High Dose"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
#   geom_rect(data = filter(ve303.dnameta, TRT.norm == "VE303 Low Dose"), aes(xmin=1, xmax=14, ymin=0, ymax=Inf), fill="lightskyblue", color="lightskyblue", show.legend = FALSE) + 
#   geom_point(size=2.5, show.legend = FALSE) + 
#   theme_bw() + 
#   labs(title= element_blank(), x="Day in Study", y="Subject", shape=1) + 
#   theme(axis.title = element_text(size=16), axis.text.y = element_text(size = 12), axis.text.x = element_text(size=12), panel.grid.major.y = element_line(size = 1, color = "darkgrey"), panel.grid.major.x = element_line(size = .2, color = "black"), strip.text.y = element_text(angle = 0, size = 12)) + 
#   #facet_grid(TRT.norm ~ ., scales = "free_y", space = "free_y") + 
#   #scale_color_manual(values = coh.cols) + 
#   scale_x_continuous(breaks = c( seq(-26, 168, by=10)))
# #theme(axis.text.y = element_blank())
# 
# ggsave(filename = here(results, paste(Sys.Date(), "Cohort_timeline-VE303-002-Blind2.pdf", sep = " ")), width = 12, height = 8)
