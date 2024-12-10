
print("Preparing array for inputting into MaasLin model")
print(tax_level)

if (tax_level == 'phylum'){ 
  d1 <-  oc.ph.long.rel.meta
  d1 <- d1 %>%  rename( organism = phylum_name)
  input_depvar <- filtphylum
}
if (tax_level == 'class'){ 
  d1 <-  oc.cl.long.rel.meta
  d1 <- d1 %>%  rename( organism = class_name)
  input_depvar <- filtclass}

if (tax_level == 'order'){ 
  d1 <-  oc.order.long.rel.meta
  d1 <- d1 %>%  rename( organism = order_name)
  input_depvar <- filtorder}

if (tax_level == 'family'){ 
  d1 <-  oc.family.long.rel.meta
  d1 <- d1 %>%  rename( organism = family_name)
  input_depvar <- filtfamily}

if (tax_level == 'genus'){ 
  d1 <-  oc.genus.long.rel.meta
  d1 <- d1 %>%  rename( organism = genus_name)
  input_depvar <- filtgenus}

if (tax_level == 'species'){ 
  d1 <-  oc.sp.long.rel.meta
  d1 <- d1 %>%  rename( organism = Species.corrected)
  input_depvar <- filtspecies}


print("Filtering out timepoints post-abx if taxplot = Pre.abx, otherwise keeping post-abx timepoints")
if(taxplot == "Pre Abx"){d11.total <- d1 %>% filter(collect.pre.post.abx == "Pre.abx")
print("excluded")}
if(taxplot == "All Abx"){d11.total <- d1 }



d11.total$collect.pre.post.rec_joined <- 
  as.factor(plyr::mapvalues(d11.total$collect.pre.post.rec, from= c("No.rec"), to=c("Pre.rec")))
days_numeric <- as.numeric(d11.total$Day.in.treatment)
days_numeric_adjusted <- 
  as.numeric(plyr::mapvalues(d11.total$Day.in.treatment, 
                             from= c("Screening","-11" ,"10",  "13" ), to=c(-11,-11 ,14 ,14)))


abx_factor <- plyr::mapvalues(d11.total$Abx.group, from= c("Vancomycin", "Other",
                                                           "Metronidazole",  "Fidaxomicin"), 
                              to=c("Vanco_Metro","Vanco_Metro","Vanco_Metro","Fidaxo"))
batch_factor <- plyr::mapvalues(d11.total$Seq.batch, from= c("Batch 1", "Batch 2", "Batch 3", "Batch 4"),
                                to=c(1, 2,3, 4))
response_factor <- plyr::mapvalues(d11.total$rec.diagnosis, from= c("Bad", "Good"), to=c(1, 2))
cohort_factor <- plyr::mapvalues(d11.total$TRT.norm, from= c("Placebo","VE303 Low Dose","VE303 High Dose"), to=c(1, 2,3))
response_factor_joined <- plyr::mapvalues(d11.total$collect.pre.post.rec_joined, from= c("Post.rec", "Pre.rec"), to=c(1, 2))
response_pre_post_no <- plyr::mapvalues(d11.total$collect.pre.post.rec, from= c("Post.rec", "Pre.rec", "No.rec"), to=c(1, 2, 3))


d11.total$Timepoint <- as.numeric(days_numeric)
d11.total$Timepoint_adj <- as.numeric(days_numeric_adjusted)
d11.total$Response <- as.numeric(response_factor)
d11.total$Response_joined_prepostnone <- as.numeric(response_pre_post_no)
d11.total$Response_joined <- as.numeric(response_factor_joined)
d11.total$Treatment <- as.numeric(cohort_factor)
d11.total$Batch <- as.factor(batch_factor)
d11.total$Abx.group_joined <- as.factor(abx_factor)

d11.total$rec.TRT1 =  paste((d11.total$TRT.norm ),d11.total$rec.diagnosis.Wk8)
d11.total$rec.TRT <- factor(d11.total$rec.TRT1,levels = c("Placebo Bad" , "Placebo Good","VE303 High Dose Bad" ,"VE303 High Dose Good",  "VE303 Low Dose Bad","VE303 Low Dose Good"))


d11.total$rec.TRT.combined =  if_else(d11.total$rec.TRT1 %in% c("VE303 High Dose Good", "VE303 Low Dose Good"), "VE303 Dosed Good",
                                       if_else(d11.total$rec.TRT1 %in% c("VE303 High Dose Bad", "VE303 Low Dose Bad"), "VE303 Dosed Bad",d11.total$rec.TRT1))
#d11.total$rec.TRT.combined <- factor(d11.total$rec.TRT.combined1,levels = c("Placebo Bad" , "Placebo Good","VE303 Dosed Bad" ,"VE303 Dosed Good"))




