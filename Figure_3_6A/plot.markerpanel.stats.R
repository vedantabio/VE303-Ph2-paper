
## Marker panel investigations 

########################################################################################

ggplot(data = ve303.dnameta, aes(x=0.01+normalized_marker_depth, y = 0.01+(100*est_relative_abundance_panel),  fill=detection_status, color=detection_status)) +
  geom_point(pch = 21, size = 1, show.legend = TRUE) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(~organism, scale = 'fixed')
ggsave(filename = here(results, paste(Sys.Date(), "MarkerPanelTest_normMarkerDepth_EstRA.pdf", sep = " ")), width = 12, height = 5)

########################################################################################

ggplot(data = ve303.dnameta, aes(x=0.01+normalized_marker_depth, y = 0.01+(100*proportion_classified_reads),  fill=detection_status, color=detection_status)) +
  geom_point(pch = 21, size = 1, show.legend = TRUE) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(~organism, scale = 'fixed')
ggsave(filename = here(results, paste(Sys.Date(), "MarkerPanelTest_normMarkerDepth_propClassifiedReads.pdf", sep = " ")), width = 12, height = 5)

########################################################################################

ggplot(data = ve303.dnameta, aes(x= 0.01+(100*proportion_classified_reads), y = 0.01+100*est_relative_abundance_panel,  fill=detection_status, color=detection_status)) +
  geom_point(pch = 21, size = 1, show.legend = TRUE) +
  scale_y_log10() +
  scale_x_log10() +
  facet_wrap(~organism)
ggsave(filename = here(results, paste(Sys.Date(), "MarkerPanelTest_propClassifiedReads_estRA.pdf", sep = " ")), width = 12, height = 5)

########################################################################################

ve303.dnameta$detection_status_figure <- 
  if_else((ve303.dnameta$detection_status == "Insufficient data"),"Low Confidence", as.character(ve303.dnameta$detection_status))

det.cols.fig <- c("Detected" = "#4daf4a", "Not detected" = "#e41a1c", "Low Confidence" = "#984ea3", "Probable" = "#377eb8")


ggplot(data = ve303.dnameta, aes(x= mean_marker_depth, y =  ( 100*proportion_markers_detected),  fill=detection_status_figure, color=detection_status_figure)) +
  geom_point(pch = 21, size = 1, show.legend = TRUE) +
  #  scale_y_log10() +
  
  scale_x_log10() +
  scale_color_manual(values = det.cols.fig) + 
  #  scale_fill_manual(values = None) + 
  scale_fill_manual(values = det.cols.fig ) + 
  theme(axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5), axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, hjust = 1, vjust =0.5) , axis.title.y = element_text(size = 14, hjust = 0.5, vjust = 1) ,
        strip.text = element_text(size = 13), legend.text=element_text(size=12), legend.title = element_blank())+
  facet_wrap(~TRT.norm) +
  labs(x="VE303 Marker Depth", y="VE303 Markers Detected (%)")
ggsave(filename = here(results, paste(Sys.Date(), "MarkerPanelTest_Coverage vs Depth.pdf", sep = " ")), width = 10, height = 4)


########################################################################################
