
if(taxonomy_system == "GTDB") { 
  
  if(abun_choice =="reads_w_children"){ 
    
      if(desired_tax_rank == "species") {
    #Species level reads summary IF USING READCOUNT WITH CHILDREN
    oc.sp.wide <- oc.abund.meta %>% 
      filter(superkingdom_name == "Bacteria", tax_rank == "species") %>%
      mutate(abun_readcount_w_children = readcount_w_children/n_reads_classified) %>%
      group_by(MGS.Sample.ID, Species.corrected) %>% 
      summarise_at(vars(abun_readcount_w_children), funs(sum)) %>%
      ungroup() %>% 
      spread(Species.corrected, abun_readcount_w_children, fill = 0)
  
    MGS.Sample.ID <- oc.sp.wide$MGS.Sample.ID
    oc.sp.wide <- oc.sp.wide %>% 
      select(-MGS.Sample.ID)
    row.names(oc.sp.wide) <- MGS.Sample.ID
    
    ### DECOSTAND NOT REQUIRED AS NORMALIZING TO N CLASSIFIED READS  #####
    #oc.sp.wide.rel <- decostand(oc.sp.wide, method = "total")
    oc.sp.wide.rel <- oc.sp.wide
    oc.sp.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.sp.wide.rel), oc.sp.wide.rel)
    
    oc.sp.long.rel.meta <- oc.sp.wide.rel %>% 
      gather(Species.corrected, rel_abund, 2:ncol(.)) %>% 
      left_join(., dnameta, by = "MGS.Sample.ID") %>% 
      group_by(Species.corrected) %>% 
      mutate(max_abund = max(rel_abund)) %>% 
      mutate(min_abund = min(rel_abund)) %>% 
      arrange(Species.corrected) %>% 
      filter(max_abund >= 0.001) %>% 
      ungroup()
    
    oc.sp.long.rel.meta <- oc.sp.long.rel.meta %>% 
      group_by(Subject.Number, Visit.Name.Norm, Species.corrected) %>% 
      mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
      ungroup()
    
    oc.sp.long.rel.meta <- oc.sp.long.rel.meta %>% 
      mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
      mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
      mutate(absolute_abund_log = (absolute_abund)+0.01)
    
    oc.sp.long.rel.meta <- oc.sp.long.rel.meta %>% 
      filter(!is.na(TRT.norm))
    
    oc.sp.long.rel.meta.ve303 <- oc.sp.long.rel.meta %>% 
      filter(grepl("VE303", Species.corrected)) 
    
    ve303.sp.sum <- oc.sp.long.rel.meta.ve303 %>% 
      group_by(MGS.Sample.ID) %>% 
      summarise(ve303_tot_perc = (sum(rel_abund) *100))
    
    oc.sp.long.rel.meta <- left_join(oc.sp.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
    oc.sp.long.rel.meta.ve303 <- left_join(oc.sp.long.rel.meta.ve303, ve303.sp.sum, by = "MGS.Sample.ID")
    
    oc.sp.long.rel.meta1 <- oc.sp.wide.rel %>% 
      gather(Species.corrected, rel_abund, 2:ncol(.)) %>% 
      left_join(., dnameta, by = "MGS.Sample.ID") %>% 
      group_by(Species.corrected, Subject.Number) %>% 
      mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
      mutate(max_abund = max(rel_abund)) %>% 
      mutate(mean_abund = mean(rel_abund)) %>% 
      arrange(Species.corrected) %>% 
      #filter(max_abund >= 0.005 & binary_abund >= 0.5 & mean_abund >= 0.001) %>% 
      filter(max_abund >= threshabun & binary_abund >= threshprev & mean_abund >= threshmean) %>% 
      ungroup()
    
    filtspecies <- (oc.sp.long.rel.meta1$Species.corrected %>% unique())
    allspecies <- (oc.sp.long.rel.meta$Species.corrected %>% unique()) }
  
  ## Generate genus-level relative abundance tables
  if(desired_tax_rank == "genus") {
    #genus level - Based on Read count - GTDB OR NCBI genus assignment
    oc.genus.wide <- oc.abund.meta %>% 
      filter(superkingdom_name == "Bacteria") %>% 
      filter(tax_rank != "no rank") %>%
      filter(!is.na(genus_name)) %>%
      mutate(abun_readcount_w_children = readcount_w_children/n_reads_classified) %>%
      select(MGS.Sample.ID, genus_name, abun_readcount_w_children) %>% 
      group_by(MGS.Sample.ID, genus_name) %>%
      mutate(reads_total = sum(abun_readcount_w_children)) %>%
      select(-abun_readcount_w_children) %>% 
      group_by(MGS.Sample.ID) %>% 
      distinct(genus_name, .keep_all = TRUE) %>% 
      spread(genus_name, reads_total, fill = 0)
    
    MGS.Sample.ID <- oc.genus.wide$MGS.Sample.ID
    oc.genus.wide <- oc.genus.wide[-1]
    row.names(oc.genus.wide) <- MGS.Sample.ID
    
    ### DECOSTAND NOT REQUIRED AS NORMALIZING TO N CLASSIFIED READS  #####
#    oc.genus.wide.rel <- decostand(oc.genus.wide, method = "total")
    oc.genus.wide.rel <- oc.genus.wide 
    oc.genus.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.genus.wide.rel), oc.genus.wide.rel)
    
    oc.genus.long.rel.meta <- oc.genus.wide.rel %>% 
      gather(genus_name, rel_abund, 2:ncol(.)) %>% 
      left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
      group_by(genus_name) %>% 
      mutate(max_abund = max(rel_abund)) %>% 
      arrange(genus_name) %>% 
      filter(max_abund >= 0.001) %>% 
      ungroup()
    
    oc.genus.long.rel.meta <- oc.genus.long.rel.meta %>% 
      group_by(Subject.Number, Visit.Name.Norm, genus_name) %>% 
      mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
      ungroup()
    
    oc.genus.long.rel.meta <- oc.genus.long.rel.meta %>% 
      mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
      mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
      mutate(absolute_abund_log = (absolute_abund)+0.01)
    
    oc.genus.long.rel.meta <- oc.genus.long.rel.meta %>% 
      filter(!is.na(TRT.norm))
    
    # Add higher taxa assignments to the table
    df.taxa <- oc.abund.meta %>% 
      select(phylum_name, genus_name, order_name, family_name, class_name, Species.corrected) %>% 
      filter(!is.na(Species.corrected)) %>% 
      distinct(Species.corrected, .keep_all = TRUE)
    
    df.taxa.genus <- df.taxa %>% 
      select(phylum_name, genus_name) %>% 
      distinct(genus_name, .keep_all = TRUE)
    
    oc.genus.long.rel.meta <- left_join(oc.genus.long.rel.meta, df.taxa.genus, by = "genus_name")
    oc.genus.long.rel.meta <- left_join(oc.genus.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
    
    oc.genus.long.rel.meta1 <- oc.genus.wide.rel %>% 
      gather(genus_name, rel_abund, 2:ncol(.)) %>% 
      left_join(., dnameta, by = "MGS.Sample.ID") %>% 
      group_by(genus_name, Subject.Number) %>% 
      mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
      mutate(max_abund = max(rel_abund)) %>% 
      arrange(genus_name) %>% 
      filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
      ungroup()
    filtgenus <- (oc.genus.long.rel.meta1$genus_name %>% unique()) }

  ## Generate family-level relative abundance tables
  if(desired_tax_rank == "family") {
    #family level - Based on Read count - GTDB family assignment
    oc.family.wide <- oc.abund.meta %>% 
      filter(superkingdom_name == "Bacteria") %>% 
      filter(tax_rank != "no rank") %>%
      filter(!is.na(family_name)) %>% 
      mutate(abun_readcount_w_children = readcount_w_children/n_reads_classified) %>%
      select(MGS.Sample.ID, family_name, abun_readcount_w_children) %>% 
      group_by(MGS.Sample.ID, family_name) %>%
      mutate(reads_total = sum(abun_readcount_w_children)) %>%
      select(-abun_readcount_w_children) %>% 
      group_by(MGS.Sample.ID) %>% 
      distinct(family_name, .keep_all = TRUE) %>% 
      spread(family_name, reads_total, fill = 0)
    
    MGS.Sample.ID <- oc.family.wide$MGS.Sample.ID
    oc.family.wide <- oc.family.wide[-1]
    row.names(oc.family.wide) <- MGS.Sample.ID
 
    ### DECOSTAND NOT REQUIRED AS NORMALIZING TO N CLASSIFIED READS  #####
    oc.family.wide.rel <-  oc.family.wide 
    oc.family.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.family.wide.rel), oc.family.wide.rel)
    
    oc.family.long.rel.meta <- oc.family.wide.rel %>% 
      gather(family_name, rel_abund, 2:ncol(.)) %>% 
      left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
      group_by(family_name) %>% 
      mutate(max_abund = max(rel_abund)) %>% 
      arrange(family_name) %>% 
      filter(max_abund >= 0.001) %>% 
      ungroup()
    
    oc.family.long.rel.meta <- oc.family.long.rel.meta %>% 
      group_by(Subject.Number, Visit.Name.Norm, family_name) %>% 
      mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
      ungroup()
    
    oc.family.long.rel.meta <- oc.family.long.rel.meta %>% 
      mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
      mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
      mutate(absolute_abund_log = (absolute_abund)+0.01)
    
    oc.family.long.rel.meta <- oc.family.long.rel.meta %>% 
      filter(!is.na(TRT.norm))
    
    # Add higher taxa assignments to the table
    df.taxa <- oc.abund.meta %>% 
      select(phylum_name, genus_name, order_name, family_name, class_name, Species.corrected) %>% 
      filter(!is.na(Species.corrected)) %>% 
      distinct(Species.corrected, .keep_all = TRUE)
    
    df.taxa.family <- df.taxa %>% 
      select(phylum_name, family_name) %>% 
      distinct(family_name, .keep_all = TRUE)
    
    oc.family.long.rel.meta <- left_join(oc.family.long.rel.meta, df.taxa.family, by = "family_name")
    oc.family.long.rel.meta <- left_join(oc.family.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
    
    
    oc.family.long.rel.meta1 <- oc.family.wide.rel %>% 
      gather(family_name, rel_abund, 2:ncol(.)) %>% 
      left_join(., dnameta, by = "MGS.Sample.ID") %>% 
      group_by(family_name, Subject.Number) %>% 
      mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
      mutate(max_abund = max(rel_abund)) %>% 
      arrange(family_name) %>% 
      filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
      ungroup()
    
    filtfamily <- (oc.family.long.rel.meta1$family_name %>% unique())  }

  ## Generate class-level relative abundance tables
  if(desired_tax_rank == "class") {
    #Class level - Based on Read count - GTDB class assignment
    oc.cl.wide <- oc.abund.meta %>% 
      filter(superkingdom_name == "Bacteria") %>% 
      filter(tax_rank != "no rank") %>% 
      filter(!is.na(class_name)) %>% 
      mutate(abun_readcount_w_children = readcount_w_children/n_reads_classified) %>%
      select(MGS.Sample.ID, class_name, abun_readcount_w_children) %>% 
      group_by(MGS.Sample.ID, class_name) %>%
      mutate(reads_total = sum(abun_readcount_w_children)) %>%
      select(-abun_readcount_w_children) %>% 
      group_by(MGS.Sample.ID) %>% 
      distinct(class_name, .keep_all = TRUE) %>% 
      spread(class_name, reads_total, fill = 0)
    
    MGS.Sample.ID <- oc.cl.wide$MGS.Sample.ID
    
    oc.cl.wide <- oc.cl.wide[-1]
    row.names(oc.cl.wide) <- MGS.Sample.ID
    
    ### DECOSTAND NOT REQUIRED AS NORMALIZING TO N CLASSIFIED READS  #####
    oc.cl.wide.rel <-  oc.cl.wide  
    oc.cl.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.cl.wide.rel), oc.cl.wide.rel)
    
    oc.cl.long.rel.meta <- oc.cl.wide.rel %>% 
      gather(class_name, rel_abund, 2:ncol(.)) %>% 
      left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
      group_by(class_name) %>% 
      mutate(max_abund = max(rel_abund)) %>% 
      arrange(class_name) %>% 
      filter(max_abund >= 0.001) %>% 
      ungroup()
    
    oc.cl.long.rel.meta <- oc.cl.long.rel.meta %>% 
      group_by(Subject.Number, Visit.Name.Norm, class_name) %>% 
      mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
      ungroup()
    
    oc.cl.long.rel.meta <- oc.cl.long.rel.meta %>% 
      mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
      mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
      mutate(absolute_abund_log = (absolute_abund)+0.01)
    
    oc.cl.long.rel.meta <- oc.cl.long.rel.meta %>% 
      filter(!is.na(TRT.norm))
    
    # Add higher taxa assignments to the table
    df.taxa <- oc.abund.meta %>% 
      select(phylum_name, class_name, order_name, family_name, genus_name, Species.corrected) %>% 
      filter(!is.na(Species.corrected)) %>% 
      distinct(Species.corrected, .keep_all = TRUE)
    
    df.taxa.cl <- df.taxa %>% 
      select(phylum_name, class_name) %>% 
      distinct(class_name, .keep_all = TRUE)
    
    oc.cl.long.rel.meta <- left_join(oc.cl.long.rel.meta, df.taxa.cl, by = "class_name")
    oc.cl.long.rel.meta <- left_join(oc.cl.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
    
    oc.cl.long.rel.meta1 <- oc.cl.wide.rel %>% 
      gather(class_name, rel_abund, 2:ncol(.)) %>% 
      left_join(., dnameta, by = "MGS.Sample.ID") %>% 
      group_by(class_name, Subject.Number) %>% 
      mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
      mutate(max_abund = max(rel_abund)) %>% 
      arrange(class_name) %>% 
      filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
      ungroup()
    
    filtclass <- (oc.cl.long.rel.meta1$class_name %>% unique())
    length(filtclass)  }
  
  ## Generate order-level relative abundance tables
  if(desired_tax_rank == "order") {
    #order level - Based on Read count - GTDB order assignment
    oc.order.wide <- oc.abund.meta %>% 
      filter(superkingdom_name == "Bacteria") %>% 
      filter(tax_rank != "no rank") %>%
      filter(!is.na(order_name)) %>% 
      mutate(abun_readcount_w_children = readcount_w_children/n_reads_classified) %>%
      select(MGS.Sample.ID, order_name, abun_readcount_w_children) %>% 
      group_by(MGS.Sample.ID, order_name) %>%
      mutate(reads_total = sum(abun_readcount_w_children)) %>%
      select(-abun_readcount_w_children) %>% 
      group_by(MGS.Sample.ID) %>% 
      distinct(order_name, .keep_all = TRUE) %>% 
      spread(order_name, reads_total, fill = 0)
    
    MGS.Sample.ID <- oc.order.wide$MGS.Sample.ID
    
    oc.order.wide <- oc.order.wide[-1]
    row.names(oc.order.wide) <- MGS.Sample.ID
    
    ### DECOSTAND NOT REQUIRED AS NORMALIZING TO N CLASSIFIED READS  ####
    oc.order.wide.rel <-  oc.order.wide 
    oc.order.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.order.wide.rel), oc.order.wide.rel)
    
    oc.order.long.rel.meta <- oc.order.wide.rel %>% 
      gather(order_name, rel_abund, 2:ncol(.)) %>% 
      left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
      group_by(order_name) %>% 
      mutate(max_abund = max(rel_abund)) %>% 
      arrange(order_name) %>% 
      filter(max_abund >= 0.001) %>% 
      ungroup()
    
    oc.order.long.rel.meta <- oc.order.long.rel.meta %>% 
      group_by(Subject.Number, Visit.Name.Norm, order_name) %>% 
      mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
      ungroup()
    
    oc.order.long.rel.meta <- oc.order.long.rel.meta %>% 
      mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
      mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
      mutate(absolute_abund_log = (absolute_abund)+0.01)
    
    oc.order.long.rel.meta <- oc.order.long.rel.meta %>% 
      filter(!is.na(TRT.norm))
    
    # Add higher taxa assignments to the table
    df.taxa.order <- df.taxa %>% 
      select(phylum_name, order_name) %>% 
      distinct(order_name, .keep_all = TRUE)
    
    oc.order.long.rel.meta <- left_join(oc.order.long.rel.meta, df.taxa.order, by = "order_name")
    oc.order.long.rel.meta <- left_join(oc.order.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
    
    oc.order.long.rel.meta1 <- oc.order.wide.rel %>% 
      gather(order_name, rel_abund, 2:ncol(.)) %>% 
      left_join(., dnameta, by = "MGS.Sample.ID") %>% 
      group_by(order_name, Subject.Number) %>% 
      mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
      mutate(max_abund = max(rel_abund)) %>% 
      arrange(order_name) %>% 
      filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
      ungroup()
    
    filtorder <- (oc.order.long.rel.meta1$order_name %>% unique()) }

  ## Generate Phylum-level relative abundance tables
  if(desired_tax_rank == "phylum") {
     oc.ph.wide <- oc.abund.meta %>% 
      filter(superkingdom_name == "Bacteria") %>% 
      filter(tax_rank != "no rank") %>%
      filter(!is.na(phylum_name)) %>% 
      mutate(abun_readcount_w_children = readcount_w_children/n_reads_classified) %>%
      select(MGS.Sample.ID, phylum_name , abun_readcount_w_children) %>%
      group_by(MGS.Sample.ID, phylum_name ) %>% 
      summarise_at(vars(abun_readcount_w_children), funs(sum)) %>%
      ungroup() %>% 
      spread(phylum_name, abun_readcount_w_children, fill = 0)
    
    MGS.Sample.ID <- oc.ph.wide$MGS.Sample.ID
    oc.ph.wide <- oc.ph.wide %>% select(-MGS.Sample.ID)
    row.names(oc.ph.wide) <- MGS.Sample.ID
    ### DECOSTAND NOT REQUIRED AS NORMALIZING TO N CLASSIFIED READS  #####
    oc.ph.wide.rel <-  oc.ph.wide 
    oc.ph.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.ph.wide.rel), oc.ph.wide.rel)
    
    oc.ph.long.rel.meta <- oc.ph.wide.rel %>% 
      gather(phylum_name, rel_abund, 2:ncol(.)) %>% 
      left_join(., dnameta, by = "MGS.Sample.ID") %>% 
      group_by(phylum_name) %>% 
      mutate(max_abund = max(rel_abund)) %>% 
      arrange(phylum_name) %>% 
      filter(max_abund >= 0.001) %>% 
      ungroup()
    
    oc.ph.long.rel.meta <- oc.ph.long.rel.meta %>% 
      mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
      mutate(absolute_abund = rel_abund * dna.norm * (2/0.5))
    
    oc.ph.long.rel.meta <- oc.ph.long.rel.meta %>% 
      group_by(Subject.Number, Visit.Name.Norm, phylum_name) %>% 
      mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
      mutate(abs.abund.avg.timepoint = mean(absolute_abund)) %>% 
      ungroup()
    
    oc.ph.long.rel.meta <- oc.ph.long.rel.meta %>% 
      filter(!is.na(TRT.norm))
    
    # Generate a wide format table of average abundance values to calc significance
    oc.ph.wide.avg.meta <- oc.ph.long.rel.meta %>% 
      select(phylum_name, abund.avg.timepoint, Visit.Name.Norm, TRT.norm, Subject.Number) %>% 
      group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>%
      distinct(phylum_name, .keep_all = TRUE) %>% 
      spread(phylum_name, abund.avg.timepoint, fill = 0)
    
    oc.ph.long.avg.meta <- oc.ph.long.rel.meta %>% 
      select(phylum_name, abund.avg.timepoint, Visit.Name.Norm, TRT.norm, Subject.Number, MGS.Sample.ID, rec.diagnosis, Day.in.treatment) %>% 
      group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>%
      distinct(phylum_name, .keep_all = TRUE)
    
    oc.ph.long.avg.meta <- left_join(oc.ph.long.avg.meta, ve303.sp.sum, by = "MGS.Sample.ID")
    filtphylum <- (oc.ph.long.rel.meta$phylum_name %>% unique())
    
    # Calculate total biomass per subject over time
    biomass.per.subj <- oc.ph.long.rel.meta %>%
      left_join(.,select(oc.abund.meta, MGS.Sample.ID, n_reads), by="MGS.Sample.ID") %>%
      group_by(Subject.Number, Visit.Name.Norm) %>% 
      mutate(tot.biomass = sum(absolute_abund)) %>% 
      group_by(Subject.Number, Visit.Name.Norm) %>% 
      mutate(biomass.tp.mean = mean(tot.biomass), biomass.tp.sd = sd(tot.biomass), biomass.tp.median = median(tot.biomass)) %>% 
      group_by(Subject.Number) %>% 
      distinct(Visit.Name.Norm, .keep_all = TRUE) %>% 
      select(Subject.Number, TRT.norm, Visit.Name.Norm, Day.in.treatment, tot.biomass, biomass.tp.mean, biomass.tp.sd, biomass.tp.median, n_reads, weight.of.collected.sample..g.) %>% 
      arrange(Subject.Number, Visit.Name.Norm) %>% 
      group_by(TRT.norm, Visit.Name.Norm) %>% 
      mutate(cohort.med.tp = median(biomass.tp.median)) %>% 
      group_by(Visit.Name.Norm) %>% 
      mutate(overall.med.tp = median(biomass.tp.median)) %>% 
      ungroup()  }
  
  ## Generate Clostridia, Bacilli, Gammaproteo relative abundance tables
  if(desired_tax_rank == "Clos_Bac_Gam"){
    oc.abund.meta$abundance_w_children[is.na(oc.abund.meta$abundance_w_children)] <- 0
    
    oc.abund.long.sp <- oc.abund.meta %>% 
      filter(superkingdom_name == "Bacteria", tax_rank == "species") %>% 
      group_by(sample_name) %>%
      distinct(Species.corrected, .keep_all = TRUE) %>% 
      group_by(Species.corrected) %>% 
      mutate(max_abund = max(abundance_w_children), est.abund.per = (abundance_w_children)*100) %>%
      ungroup()
    
    clost.top.sp.list <- oc.abund.long.sp %>% 
      filter(abundance_w_children >0) %>%
      filter(grepl("Clostridia", class_name)) %>%
      group_by(sample_id) %>% 
      top_n(n=2, wt = abundance_w_children) %>% 
      ungroup() %>% 
      distinct(Species.corrected)
    
    clost.top.sp <- oc.abund.long.sp %>% 
      filter(Species.corrected %in% clost.top.sp.list$Species.corrected)
    
    gtdb.genus.sum <- oc.abund.long.sp %>% 
      group_by(MGS.Sample.ID, phylum_name, genus_name, Visit.Name.Norm) %>% 
      summarise(total_abund = (sum(abundance_w_children)*100)) %>%
      select(-Visit.Name.Norm) %>%
      left_join(., dnameta, by = "MGS.Sample.ID")
    
    gtdb.gen.sum.clost <- gtdb.genus.sum %>% 
      filter(genus_name %in% clost.top.sp$genus_name) %>%
      mutate(total_abund_log = total_abund + 0.01)
    
    #######################################################
    
    bac.top.sp.list <- oc.abund.long.sp %>% 
      filter(abundance_w_children >0) %>%
      filter(grepl("Bacilli", class_name)) %>%
      group_by(sample_id) %>% 
      top_n(n=2, wt = abundance_w_children) %>% 
      ungroup() %>% 
      distinct(Species.corrected)
    
    bac.top.sp <- oc.abund.long.sp %>% 
      filter(Species.corrected %in% bac.top.sp.list$Species.corrected)
    
    gtdb.gen.sum.bac <- gtdb.genus.sum %>% 
      filter(genus_name %in% bac.top.sp$genus_name) %>%
      mutate(total_abund_log = total_abund + 0.01)
    ######################
    
    gprot.top.sp.list <- oc.abund.long.sp %>% 
      filter(abundance_w_children >0) %>%
      filter(grepl("Gammaproteobacteria", class_name)) %>%
      group_by(sample_id) %>% 
      top_n(n=2, wt = abundance_w_children) %>% 
      ungroup() %>% 
      distinct(Species.corrected)
    
    gprot.top.sp <- oc.abund.long.sp %>% 
      filter(Species.corrected %in% gprot.top.sp.list$Species.corrected)
    
    gtdb.gen.sum.gprot <- gtdb.genus.sum %>% 
      filter(genus_name %in% gprot.top.sp$genus_name) %>%
      mutate(total_abund_log = total_abund + 0.01) }
    ################################## 
    }
  
    if(abun_choice =="abun_w_children"){
        if(desired_tax_rank == "species") {
        #Species level abundance summary
        oc.sp.wide <- oc.abund.meta %>% 
          filter(superkingdom_name == "Bacteria", tax_rank == "species") %>% 
          select(MGS.Sample.ID, Species.corrected, abundance_w_children) %>%
          group_by(MGS.Sample.ID, Species.corrected) %>% 
          summarise_at(vars(abundance_w_children), funs(sum)) %>%
          ungroup() %>% 
          spread(Species.corrected, abundance_w_children, fill = 0)
        
        
        MGS.Sample.ID <- oc.sp.wide$MGS.Sample.ID
        
        oc.sp.wide <- oc.sp.wide %>% 
          select(-MGS.Sample.ID)
        
        row.names(oc.sp.wide) <- MGS.Sample.ID

        oc.sp.wide.rel <- oc.sp.wide
        oc.sp.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.sp.wide.rel), oc.sp.wide.rel)
        
        oc.sp.long.rel.meta <- oc.sp.wide.rel %>% 
          gather(Species.corrected, rel_abund, 2:ncol(.)) %>% 
          left_join(., dnameta, by = "MGS.Sample.ID") %>% 
          group_by(Species.corrected) %>% 
          mutate(max_abund = max(rel_abund)) %>% 
          mutate(min_abund = min(rel_abund)) %>% 
          arrange(Species.corrected) %>% 
          filter(max_abund >= 0.001) %>% 
          ungroup()
        
        oc.sp.long.rel.meta <- oc.sp.long.rel.meta %>% 
          group_by(Subject.Number, Visit.Name.Norm, Species.corrected) %>% 
          mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
          ungroup()
        
        oc.sp.long.rel.meta <- oc.sp.long.rel.meta %>% 
          mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
          mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
          mutate(absolute_abund_log = (absolute_abund)+0.01)
        
        oc.sp.long.rel.meta <- oc.sp.long.rel.meta %>% 
          filter(!is.na(TRT.norm))
        
        oc.sp.long.rel.meta.ve303 <- oc.sp.long.rel.meta %>% 
          filter(grepl("VE303", Species.corrected)) 
        
        ve303.sp.sum <- oc.sp.long.rel.meta.ve303 %>% 
          group_by(MGS.Sample.ID) %>% 
          summarise(ve303_tot_perc = (sum(rel_abund) *100))
        
        oc.sp.long.rel.meta <- left_join(oc.sp.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
        oc.sp.long.rel.meta.ve303 <- left_join(oc.sp.long.rel.meta.ve303, ve303.sp.sum, by = "MGS.Sample.ID")
        
        oc.sp.long.rel.meta1 <- oc.sp.wide.rel %>% 
          gather(Species.corrected, rel_abund, 2:ncol(.)) %>% 
          left_join(., dnameta, by = "MGS.Sample.ID") %>% 
          group_by(Species.corrected, Subject.Number) %>% 
          mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
          mutate(max_abund = max(rel_abund)) %>% 
          mutate(mean_abund = mean(rel_abund)) %>% 
          arrange(Species.corrected) %>% 
          #filter(max_abund >= 0.005 & binary_abund >= 0.5 & mean_abund >= 0.001) %>% 
          filter(max_abund >= threshabun & binary_abund >= threshprev & mean_abund >= threshmean) %>% 
          ungroup()
        
        filtspecies <- (oc.sp.long.rel.meta1$Species.corrected %>% unique())
        allspecies <- (oc.sp.long.rel.meta$Species.corrected %>% unique()) }
    

      ## Generate genus-level relative abundance tables
      if(desired_tax_rank == "genus") {
        #genus level - Based on Read count - GTDB OR NCBI genus assignment
        oc.genus.wide <- oc.abund.meta %>%
          filter(superkingdom_name == "Bacteria") %>% 
          filter(tax_rank != "no rank") %>%
          filter(!is.na(genus_name)) %>%
          select(MGS.Sample.ID, genus_name, abundance_w_children) %>% 
          group_by(MGS.Sample.ID, genus_name) %>%
          mutate(reads_total = sum(abundance_w_children)) %>%
          select(-abundance_w_children) %>% 
          group_by(MGS.Sample.ID) %>% 
          distinct(genus_name, .keep_all = TRUE) %>% 
          spread(genus_name, reads_total, fill = 0)
        
        MGS.Sample.ID <- oc.genus.wide$MGS.Sample.ID
        oc.genus.wide <- oc.genus.wide[-1]
        row.names(oc.genus.wide) <- MGS.Sample.ID
        
        oc.genus.wide.rel <- oc.genus.wide 
        oc.genus.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.genus.wide.rel), oc.genus.wide.rel)
        
        oc.genus.long.rel.meta <- oc.genus.wide.rel %>% 
          gather(genus_name, rel_abund, 2:ncol(.)) %>% 
          left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
          group_by(genus_name) %>% 
          mutate(max_abund = max(rel_abund)) %>% 
          arrange(genus_name) %>% 
          filter(max_abund >= 0.001) %>% 
          ungroup()
        
        oc.genus.long.rel.meta <- oc.genus.long.rel.meta %>% 
          group_by(Subject.Number, Visit.Name.Norm, genus_name) %>% 
          mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
          ungroup()
        
        oc.genus.long.rel.meta <- oc.genus.long.rel.meta %>% 
          mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
          mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
          mutate(absolute_abund_log = (absolute_abund)+0.01)
        
        oc.genus.long.rel.meta <- oc.genus.long.rel.meta %>% 
          filter(!is.na(TRT.norm))
        
        # Add higher taxa assignments to the table
        df.taxa <- oc.abund.meta %>% 
          select(phylum_name, genus_name, order_name, family_name, class_name, Species.corrected) %>% 
          filter(!is.na(Species.corrected)) %>% 
          distinct(Species.corrected, .keep_all = TRUE)
        
        df.taxa.genus <- df.taxa %>% 
          select(phylum_name, genus_name) %>% 
          distinct(genus_name, .keep_all = TRUE)
        
        oc.genus.long.rel.meta <- left_join(oc.genus.long.rel.meta, df.taxa.genus, by = "genus_name")
        oc.genus.long.rel.meta <- left_join(oc.genus.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
        
        oc.genus.long.rel.meta1 <- oc.genus.wide.rel %>% 
          gather(genus_name, rel_abund, 2:ncol(.)) %>% 
          left_join(., dnameta, by = "MGS.Sample.ID") %>% 
          group_by(genus_name, Subject.Number) %>% 
          mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
          mutate(max_abund = max(rel_abund)) %>% 
          arrange(genus_name) %>% 
          filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
          ungroup()
        
        filtgenus <- (oc.genus.long.rel.meta1$genus_name %>% unique()) }
    
      ## Generate family-level relative abundance tables
      if(desired_tax_rank == "family") {
        #family level - Based on abundance - GTDB family assignment
        oc.family.wide <- oc.abund.meta %>% 
          filter(superkingdom_name == "Bacteria") %>% 
          filter(tax_rank != "no rank") %>%
          filter(!is.na(family_name)) %>% 
          select(MGS.Sample.ID, family_name, abundance_w_children) %>% 
          group_by(MGS.Sample.ID, family_name) %>%
          mutate(reads_total = sum(abundance_w_children)) %>%
          select(-abundance_w_children) %>% 
          group_by(MGS.Sample.ID) %>% 
          distinct(family_name, .keep_all = TRUE) %>% 
          spread(family_name, reads_total, fill = 0)
        
        MGS.Sample.ID <- oc.family.wide$MGS.Sample.ID
        oc.family.wide <- oc.family.wide[-1]
        row.names(oc.family.wide) <- MGS.Sample.ID
        
        oc.family.wide.rel <- oc.family.wide 
        oc.family.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.family.wide.rel), oc.family.wide.rel)
        
        oc.family.long.rel.meta <- oc.family.wide.rel %>% 
          gather(family_name, rel_abund, 2:ncol(.)) %>% 
          left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
          group_by(family_name) %>% 
          mutate(max_abund = max(rel_abund)) %>% 
          arrange(family_name) %>% 
          filter(max_abund >= 0.001) %>% 
          ungroup()
        
        oc.family.long.rel.meta <- oc.family.long.rel.meta %>% 
          group_by(Subject.Number, Visit.Name.Norm, family_name) %>% 
          mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
          ungroup()
        
        oc.family.long.rel.meta <- oc.family.long.rel.meta %>% 
          mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
          mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
          mutate(absolute_abund_log = (absolute_abund)+0.01)
        
        oc.family.long.rel.meta <- oc.family.long.rel.meta %>% 
          filter(!is.na(TRT.norm))
        
        # Add higher taxa assignments to the table
        df.taxa <- oc.abund.meta %>% 
          select(phylum_name, genus_name, order_name, family_name, class_name, Species.corrected) %>% 
          filter(!is.na(Species.corrected)) %>% 
          distinct(Species.corrected, .keep_all = TRUE)
        
        df.taxa.family <- df.taxa %>% 
          select(phylum_name, family_name) %>% 
          distinct(family_name, .keep_all = TRUE)
        
        oc.family.long.rel.meta <- left_join(oc.family.long.rel.meta, df.taxa.family, by = "family_name")
        oc.family.long.rel.meta <- left_join(oc.family.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
        
        oc.family.long.rel.meta1 <- oc.family.wide.rel %>% 
          gather(family_name, rel_abund, 2:ncol(.)) %>% 
          left_join(., dnameta, by = "MGS.Sample.ID") %>% 
          group_by(family_name, Subject.Number) %>% 
          mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
          mutate(max_abund = max(rel_abund)) %>% 
          arrange(family_name) %>% 
          filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
          ungroup()
        
        filtfamily <- (oc.family.long.rel.meta1$family_name %>% unique()) }
    
      ## Generate class-level relative abundance tables
      if(desired_tax_rank == "class") {
        #Class level - Based on Read count - GTDB class assignment
        oc.cl.wide <- oc.abund.meta %>% 
          filter(superkingdom_name == "Bacteria") %>% 
          filter(tax_rank != "no rank") %>%
          filter(!is.na(class_name)) %>% 
          select(MGS.Sample.ID, class_name, abundance_w_children) %>% 
          group_by(MGS.Sample.ID, class_name) %>%
          mutate(reads_total = sum(abundance_w_children)) %>%
          select(-abundance_w_children) %>% 
          group_by(MGS.Sample.ID) %>% 
          distinct(class_name, .keep_all = TRUE) %>% 
          spread(class_name, reads_total, fill = 0)
        
        MGS.Sample.ID <- oc.cl.wide$MGS.Sample.ID
        
        oc.cl.wide <- oc.cl.wide[-1]
        row.names(oc.cl.wide) <- MGS.Sample.ID
        
        oc.cl.wide.rel <-  oc.cl.wide 
        oc.cl.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.cl.wide.rel), oc.cl.wide.rel)
        
        oc.cl.long.rel.meta <- oc.cl.wide.rel %>% 
          gather(class_name, rel_abund, 2:ncol(.)) %>% 
          left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
          group_by(class_name) %>% 
          mutate(max_abund = max(rel_abund)) %>% 
          arrange(class_name) %>% 
          filter(max_abund >= 0.001) %>% 
          ungroup()
        
        oc.cl.long.rel.meta <- oc.cl.long.rel.meta %>% 
          group_by(Subject.Number, Visit.Name.Norm, class_name) %>% 
          mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
          ungroup()
        
        oc.cl.long.rel.meta <- oc.cl.long.rel.meta %>% 
          mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
          mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
          mutate(absolute_abund_log = (absolute_abund)+0.01)
        
        oc.cl.long.rel.meta <- oc.cl.long.rel.meta %>% 
          filter(!is.na(TRT.norm))
        
        # Add higher taxa assignments to the table
        df.taxa <- oc.abund.meta %>% 
          select(phylum_name, class_name, order_name, family_name, genus_name, Species.corrected) %>% 
          filter(!is.na(Species.corrected)) %>% 
          distinct(Species.corrected, .keep_all = TRUE)
        
        df.taxa.cl <- df.taxa %>% 
          select(phylum_name, class_name) %>% 
          distinct(class_name, .keep_all = TRUE)
        
        oc.cl.long.rel.meta <- left_join(oc.cl.long.rel.meta, df.taxa.cl, by = "class_name")
        oc.cl.long.rel.meta <- left_join(oc.cl.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
        
        
        oc.cl.long.rel.meta1 <- oc.cl.wide.rel %>% 
          gather(class_name, rel_abund, 2:ncol(.)) %>% 
          left_join(., dnameta, by = "MGS.Sample.ID") %>% 
          group_by(class_name, Subject.Number) %>% 
          mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
          mutate(max_abund = max(rel_abund)) %>% 
          arrange(class_name) %>% 
          filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
          ungroup()
        
        filtclass <- (oc.cl.long.rel.meta1$class_name %>% unique()) }
      
      ## Generate order-level relative abundance tables
      if(desired_tax_rank == "order") {
        
        #order level - Based on abun - GTDB order assignment
        oc.order.wide <- oc.abund.meta %>% 
          filter(superkingdom_name == "Bacteria") %>% 
          filter(tax_rank != "no rank") %>%
          filter(!is.na(order_name)) %>% 
          select(MGS.Sample.ID, order_name,abundance_w_children) %>% 
          group_by(MGS.Sample.ID, order_name) %>%
          mutate(reads_total = sum(abundance_w_children)) %>%
          select(-abundance_w_children) %>% 
          group_by(MGS.Sample.ID) %>% 
          distinct(order_name, .keep_all = TRUE) %>% 
          spread(order_name, reads_total, fill = 0)
        
        MGS.Sample.ID <- oc.order.wide$MGS.Sample.ID
        oc.order.wide <- oc.order.wide[-1]
        row.names(oc.order.wide) <- MGS.Sample.ID
        
        oc.order.wide.rel <-  oc.order.wide
        oc.order.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.order.wide.rel), oc.order.wide.rel)
        
        oc.order.long.rel.meta <- oc.order.wide.rel %>% 
          gather(order_name, rel_abund, 2:ncol(.)) %>% 
          left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
          group_by(order_name) %>% 
          mutate(max_abund = max(rel_abund)) %>% 
          arrange(order_name) %>% 
          filter(max_abund >= 0.001) %>% 
          ungroup()
        
        oc.order.long.rel.meta <- oc.order.long.rel.meta %>% 
          group_by(Subject.Number, Visit.Name.Norm, order_name) %>% 
          mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
          ungroup()
        
        oc.order.long.rel.meta <- oc.order.long.rel.meta %>% 
          mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
          mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
          mutate(absolute_abund_log = (absolute_abund)+0.01)
        
        oc.order.long.rel.meta <- oc.order.long.rel.meta %>% 
          filter(!is.na(TRT.norm))
        
        # Add higher taxa assignments to the table
        
        df.taxa.order <- df.taxa %>% 
          select(phylum_name, order_name) %>% 
          distinct(order_name, .keep_all = TRUE)
        
        oc.order.long.rel.meta <- left_join(oc.order.long.rel.meta, df.taxa.order, by = "order_name")
        oc.order.long.rel.meta <- left_join(oc.order.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
        
        
        oc.order.long.rel.meta1 <- oc.order.wide.rel %>% 
          gather(order_name, rel_abund, 2:ncol(.)) %>% 
          left_join(., dnameta, by = "MGS.Sample.ID") %>% 
          group_by(order_name, Subject.Number) %>% 
          mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
          mutate(max_abund = max(rel_abund)) %>% 
          arrange(order_name) %>% 
          filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
          ungroup()
        
        filtorder <- (oc.order.long.rel.meta1$order_name %>% unique()) }
      
      
      ## Generate Phylum-level relative abundance tables
      if(desired_tax_rank == "phylum") {
        
        oc.ph.wide <- oc.abund.meta %>% 
          filter(!is.na(phylum_name)) %>% 
          filter(superkingdom_name == "Bacteria", tax_rank == "species") %>% 
          select(MGS.Sample.ID, phylum_name, abundance_w_children) %>%
          group_by(MGS.Sample.ID, phylum_name) %>% 
          summarise_at(vars(abundance_w_children), funs(sum)) %>%
          ungroup() %>% 
          spread(phylum_name, abundance_w_children, fill = 0) 
         
        MGS.Sample.ID <- oc.ph.wide$MGS.Sample.ID
        oc.ph.wide <- oc.ph.wide %>% 
          select(-MGS.Sample.ID)
        
        row.names(oc.ph.wide) <- MGS.Sample.ID
        
        oc.ph.wide.rel <-  oc.ph.wide  
        oc.ph.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.ph.wide.rel), oc.ph.wide.rel)
        
        oc.ph.long.rel.meta <- oc.ph.wide.rel %>% 
          gather(phylum_name, rel_abund, 2:ncol(.)) %>% 
          left_join(., dnameta, by = "MGS.Sample.ID") %>% 
          group_by(phylum_name) %>% 
          mutate(max_abund = max(rel_abund)) %>% 
          arrange(phylum_name) %>% 
          filter(max_abund >= 0.001) %>% 
          ungroup()
        
        oc.ph.long.rel.meta <- oc.ph.long.rel.meta %>% 
          mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
          mutate(absolute_abund = rel_abund * dna.norm * (2/0.5))
        
        oc.ph.long.rel.meta <- oc.ph.long.rel.meta %>% 
          group_by(Subject.Number, Visit.Name.Norm, phylum_name) %>% 
          mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
          mutate(abs.abund.avg.timepoint = mean(absolute_abund)) %>% 
          ungroup()
        
        oc.ph.long.rel.meta <- oc.ph.long.rel.meta %>% 
          filter(!is.na(TRT.norm))
        
        # Generate a wide format table of average abundance values to calc significance
        oc.ph.wide.avg.meta <- oc.ph.long.rel.meta %>% 
          select(phylum_name, abund.avg.timepoint, Visit.Name.Norm, TRT.norm, Subject.Number) %>% 
          group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>%
          distinct(phylum_name, .keep_all = TRUE) %>% 
          spread(phylum_name, abund.avg.timepoint, fill = 0)
        
        oc.ph.long.avg.meta <- oc.ph.long.rel.meta %>% 
          select(phylum_name, abund.avg.timepoint, Visit.Name.Norm, TRT.norm, Subject.Number, MGS.Sample.ID, rec.diagnosis, Day.in.treatment) %>% 
          group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>%
          distinct(phylum_name, .keep_all = TRUE)
        
        oc.ph.long.avg.meta <- left_join(oc.ph.long.avg.meta, ve303.sp.sum, by = "MGS.Sample.ID")
        filtphylum <- (oc.ph.long.rel.meta$phylum_name %>% unique())
        
        # Calculate total biomass per subject over time
        biomass.per.subj <- oc.ph.long.rel.meta %>%
          left_join(.,select(oc.abund.meta, MGS.Sample.ID, n_reads), by="MGS.Sample.ID") %>%
          group_by(Subject.Number, Visit.Name.Norm) %>% 
          mutate(tot.biomass = sum(absolute_abund)) %>% 
          group_by(Subject.Number, Visit.Name.Norm) %>% 
          mutate(biomass.tp.mean = mean(tot.biomass), biomass.tp.sd = sd(tot.biomass), biomass.tp.median = median(tot.biomass)) %>% 
          group_by(Subject.Number) %>% 
          distinct(Visit.Name.Norm, .keep_all = TRUE) %>% 
          select(Subject.Number, TRT.norm, Visit.Name.Norm, Day.in.treatment, tot.biomass, biomass.tp.mean, biomass.tp.sd, biomass.tp.median, n_reads, weight.of.collected.sample..g.) %>% 
          arrange(Subject.Number, Visit.Name.Norm) %>% 
          group_by(TRT.norm, Visit.Name.Norm) %>% 
          mutate(cohort.med.tp = median(biomass.tp.median)) %>% 
          group_by(Visit.Name.Norm) %>% 
          mutate(overall.med.tp = median(biomass.tp.median)) %>% 
          ungroup() }
      
      ## Generate Clostridia, Bacilli, Gammaproteo relative abundance tables
      if(desired_tax_rank == "Clos_Bac_Gam"){
        oc.abund.meta$abundance_w_children[is.na(oc.abund.meta$abundance_w_children)] <- 0
          oc.abund.long.sp <- oc.abund.meta %>% 
          filter(superkingdom_name == "Bacteria", tax_rank == "species") %>% 
          group_by(sample_name) %>%
          distinct(Species.corrected, .keep_all = TRUE) %>% 
          group_by(Species.corrected) %>% 
          mutate(max_abund = max(abundance_w_children), est.abund.per = (abundance_w_children)*100) %>%
          ungroup()
        
        clost.top.sp.list <- oc.abund.long.sp %>% 
          filter(abundance_w_children >0) %>%
          filter(grepl("Clostridia", class_name)) %>%
          group_by(sample_id) %>% 
          top_n(n=2, wt = abundance_w_children) %>% 
          ungroup() %>% 
          distinct(Species.corrected)
        
        clost.top.sp <- oc.abund.long.sp %>% 
          filter(Species.corrected %in% clost.top.sp.list$Species.corrected)
        
        gtdb.genus.sum <- oc.abund.long.sp %>% 
          group_by(MGS.Sample.ID, phylum_name, genus_name, Visit.Name.Norm) %>% 
          summarise(total_abund = (sum(abundance_w_children)*100)) %>%
          select(-Visit.Name.Norm) %>%
          left_join(., dnameta, by = "MGS.Sample.ID")
        
        gtdb.gen.sum.clost <- gtdb.genus.sum %>% 
          filter(genus_name %in% clost.top.sp$genus_name) %>%
          mutate(total_abund_log = total_abund + 0.01)
        
        #######################################################
        
        bac.top.sp.list <- oc.abund.long.sp %>% 
          filter(abundance_w_children >0) %>%
          filter(grepl("Bacilli", class_name)) %>%
          group_by(sample_id) %>% 
          top_n(n=2, wt = abundance_w_children) %>% 
          ungroup() %>% 
          distinct(Species.corrected)
        
        bac.top.sp <- oc.abund.long.sp %>% 
          filter(Species.corrected %in% bac.top.sp.list$Species.corrected)
        
        gtdb.gen.sum.bac <- gtdb.genus.sum %>% 
          filter(genus_name %in% bac.top.sp$genus_name) %>%
          mutate(total_abund_log = total_abund + 0.01)
        ######################
        
        gprot.top.sp.list <- oc.abund.long.sp %>% 
          filter(abundance_w_children >0) %>%
          filter(grepl("Gammaproteobacteria", class_name)) %>%
          group_by(sample_id) %>% 
          top_n(n=2, wt = abundance_w_children) %>% 
          ungroup() %>% 
          distinct(Species.corrected)
        
        gprot.top.sp <- oc.abund.long.sp %>% 
          filter(Species.corrected %in% gprot.top.sp.list$Species.corrected)
        
        gtdb.gen.sum.gprot <- gtdb.genus.sum %>% 
          filter(genus_name %in% gprot.top.sp$genus_name) %>%
          mutate(total_abund_log = total_abund + 0.01) }
        ##################################
    }
}




if(taxonomy_system %in% c("NCBI", "NCBI OC")) { 
  
  if(abun_choice =="reads_w_children"){ 
    if(desired_tax_rank == "species") {
      #Species level reads summary IF USING READCOUNT WITH CHILDREN
      oc.sp.wide <- oc.abund.meta %>% 
        filter(superkingdom_name == "Bacteria", tax_rank == "species") %>%
        mutate(abun_readcount_w_children = readcount_w_children/n_reads_classified) %>%
        group_by(MGS.Sample.ID, Species.corrected) %>% 
        summarise_at(vars(abun_readcount_w_children), funs(sum)) %>%
        ungroup() %>% 
        spread(Species.corrected, abun_readcount_w_children, fill = 0)
      
      MGS.Sample.ID <- oc.sp.wide$MGS.Sample.ID
      oc.sp.wide <- oc.sp.wide %>% 
        select(-MGS.Sample.ID)
      row.names(oc.sp.wide) <- MGS.Sample.ID
      
      ### DECOSTAND NOT REQUIRED AS NORMALIZING TO N CLASSIFIED READS  #####
      #oc.sp.wide.rel <- decostand(oc.sp.wide, method = "total")
      oc.sp.wide.rel <- oc.sp.wide
      oc.sp.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.sp.wide.rel), oc.sp.wide.rel)
      
      oc.sp.long.rel.meta <- oc.sp.wide.rel %>% 
        gather(Species.corrected, rel_abund, 2:ncol(.)) %>% 
        left_join(., dnameta, by = "MGS.Sample.ID") %>% 
        group_by(Species.corrected) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        mutate(min_abund = min(rel_abund)) %>% 
        arrange(Species.corrected) %>% 
        filter(max_abund >= 0.001) %>% 
        ungroup()
      
      oc.sp.long.rel.meta <- oc.sp.long.rel.meta %>% 
        group_by(Subject.Number, Visit.Name.Norm, Species.corrected) %>% 
        mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
        ungroup()
      
      oc.sp.long.rel.meta <- oc.sp.long.rel.meta %>% 
        mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
        mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
        mutate(absolute_abund_log = (absolute_abund)+0.01)
      
      oc.sp.long.rel.meta <- oc.sp.long.rel.meta %>% 
        filter(!is.na(TRT.norm))
      
      oc.sp.long.rel.meta.ve303 <- oc.sp.long.rel.meta %>% 
        filter(grepl("VE303", Species.corrected)) 
      
      ve303.sp.sum <- oc.sp.long.rel.meta.ve303 %>% 
        group_by(MGS.Sample.ID) %>% 
        summarise(ve303_tot_perc = (sum(rel_abund) *100))
      
      oc.sp.long.rel.meta <- left_join(oc.sp.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
      oc.sp.long.rel.meta.ve303 <- left_join(oc.sp.long.rel.meta.ve303, ve303.sp.sum, by = "MGS.Sample.ID")
      
      oc.sp.long.rel.meta1 <- oc.sp.wide.rel %>% 
        gather(Species.corrected, rel_abund, 2:ncol(.)) %>% 
        left_join(., dnameta, by = "MGS.Sample.ID") %>% 
        group_by(Species.corrected, Subject.Number) %>% 
        mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        mutate(mean_abund = mean(rel_abund)) %>% 
        arrange(Species.corrected) %>% 
        #filter(max_abund >= 0.005 & binary_abund >= 0.5 & mean_abund >= 0.001) %>% 
        filter(max_abund >= threshabun & binary_abund >= threshprev & mean_abund >= threshmean) %>% 
        ungroup()
      
      filtspecies <- (oc.sp.long.rel.meta1$Species.corrected %>% unique())
      allspecies <- (oc.sp.long.rel.meta$Species.corrected %>% unique()) }
    
    ## Generate genus-level relative abundance tables
    if(desired_tax_rank == "genus") {
      #genus level - Based on Read count -  OR NCBI genus assignment
      oc.genus.wide <- oc.abund.meta %>% 
        filter(superkingdom_name == "Bacteria", tax_rank == "genus") %>% 
        mutate(abun_readcount_w_children = readcount_w_children/n_reads_classified) %>%
        select(MGS.Sample.ID, genus_name, abun_readcount_w_children) %>% 
        spread(genus_name, abun_readcount_w_children, fill = 0)
      
      MGS.Sample.ID <- oc.genus.wide$MGS.Sample.ID
      oc.genus.wide <- oc.genus.wide[-1]
      row.names(oc.genus.wide) <- MGS.Sample.ID
      
      ### DECOSTAND NOT REQUIRED AS NORMALIZING TO N CLASSIFIED READS  #####
      #    oc.genus.wide.rel <- decostand(oc.genus.wide, method = "total")
      oc.genus.wide.rel <- oc.genus.wide 
      oc.genus.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.genus.wide.rel), oc.genus.wide.rel)
      
      oc.genus.long.rel.meta <- oc.genus.wide.rel %>% 
        gather(genus_name, rel_abund, 2:ncol(.)) %>% 
        left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
        group_by(genus_name) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(genus_name) %>% 
        filter(max_abund >= 0.001) %>% 
        ungroup()
      
      oc.genus.long.rel.meta <- oc.genus.long.rel.meta %>% 
        group_by(Subject.Number, Visit.Name.Norm, genus_name) %>% 
        mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
        ungroup()
      
      oc.genus.long.rel.meta <- oc.genus.long.rel.meta %>% 
        mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
        mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
        mutate(absolute_abund_log = (absolute_abund)+0.01)
      
      oc.genus.long.rel.meta <- oc.genus.long.rel.meta %>% 
        filter(!is.na(TRT.norm))
      
      # Add higher taxa assignments to the table
      df.taxa <- oc.abund.meta %>% 
        select(phylum_name, genus_name, order_name, family_name, class_name, Species.corrected) %>% 
        filter(!is.na(Species.corrected)) %>% 
        distinct(Species.corrected, .keep_all = TRUE)
      
      df.taxa.genus <- df.taxa %>% 
        select(phylum_name, genus_name) %>% 
        distinct(genus_name, .keep_all = TRUE)
      
      oc.genus.long.rel.meta <- left_join(oc.genus.long.rel.meta, df.taxa.genus, by = "genus_name")
      oc.genus.long.rel.meta <- left_join(oc.genus.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
      
      oc.genus.long.rel.meta1 <- oc.genus.wide.rel %>% 
        gather(genus_name, rel_abund, 2:ncol(.)) %>% 
        left_join(., dnameta, by = "MGS.Sample.ID") %>% 
        group_by(genus_name, Subject.Number) %>% 
        mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(genus_name) %>% 
        filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
        ungroup()
      filtgenus <- (oc.genus.long.rel.meta1$genus_name %>% unique()) }
    
    ## Generate family-level relative abundance tables
    if(desired_tax_rank == "family") {
      #family level - Based on Read count -  family assignment
      oc.family.wide <- oc.abund.meta %>% 
        filter(superkingdom_name == "Bacteria", tax_rank == "family") %>% 
        mutate(abun_readcount_w_children = readcount_w_children/n_reads_classified) %>%
        select(MGS.Sample.ID, tax_name, abun_readcount_w_children) %>% 
        mutate(family_name = tax_name) %>%
        select(-tax_name) %>%
        spread(family_name, abun_readcount_w_children, fill = 0)
      
      MGS.Sample.ID <- oc.family.wide$MGS.Sample.ID
      oc.family.wide <- oc.family.wide[-1]
      row.names(oc.family.wide) <- MGS.Sample.ID
      
      ### DECOSTAND NOT REQUIRED AS NORMALIZING TO N CLASSIFIED READS  #####
      oc.family.wide.rel <-  oc.family.wide 
      oc.family.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.family.wide.rel), oc.family.wide.rel)
      
      oc.family.long.rel.meta <- oc.family.wide.rel %>% 
        gather(family_name, rel_abund, 2:ncol(.)) %>% 
        left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
        group_by(family_name) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(family_name) %>% 
        filter(max_abund >= 0.001) %>% 
        ungroup()
      
      oc.family.long.rel.meta <- oc.family.long.rel.meta %>% 
        group_by(Subject.Number, Visit.Name.Norm, family_name) %>% 
        mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
        ungroup()
      
      oc.family.long.rel.meta <- oc.family.long.rel.meta %>% 
        mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
        mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
        mutate(absolute_abund_log = (absolute_abund)+0.01)
      
      oc.family.long.rel.meta <- oc.family.long.rel.meta %>% 
        filter(!is.na(TRT.norm))
      
      # Add higher taxa assignments to the table
      df.taxa <- oc.abund.meta %>% 
        select(phylum_name, genus_name, order_name, family_name, class_name, Species.corrected) %>% 
        filter(!is.na(Species.corrected)) %>% 
        distinct(Species.corrected, .keep_all = TRUE)
      
      df.taxa.family <- df.taxa %>% 
        select(phylum_name, family_name) %>% 
        distinct(family_name, .keep_all = TRUE)
      
      oc.family.long.rel.meta <- left_join(oc.family.long.rel.meta, df.taxa.family, by = "family_name")
      oc.family.long.rel.meta <- left_join(oc.family.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
      
      
      oc.family.long.rel.meta1 <- oc.family.wide.rel %>% 
        gather(family_name, rel_abund, 2:ncol(.)) %>% 
        left_join(., dnameta, by = "MGS.Sample.ID") %>% 
        group_by(family_name, Subject.Number) %>% 
        mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(family_name) %>% 
        filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
        ungroup()
      
      filtfamily <- (oc.family.long.rel.meta1$family_name %>% unique())  }
    
    ## Generate class-level relative abundance tables
    if(desired_tax_rank == "class") {
      #Class level - Based on Read count -  class assignment
      oc.cl.wide <- oc.abund.meta %>% 
        filter(superkingdom_name == "Bacteria", tax_rank == "class") %>% 
        mutate(abun_readcount_w_children = readcount_w_children/n_reads_classified) %>%
        select(MGS.Sample.ID, class_name, abun_readcount_w_children) %>% 
        spread(class_name,  abun_readcount_w_children, fill = 0)
      
      MGS.Sample.ID <- oc.cl.wide$MGS.Sample.ID
      
      oc.cl.wide <- oc.cl.wide[-1]
      row.names(oc.cl.wide) <- MGS.Sample.ID
      
      ### DECOSTAND NOT REQUIRED AS NORMALIZING TO N CLASSIFIED READS  #####
      oc.cl.wide.rel <-  oc.cl.wide  
      oc.cl.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.cl.wide.rel), oc.cl.wide.rel)
      
      oc.cl.long.rel.meta <- oc.cl.wide.rel %>% 
        gather(class_name, rel_abund, 2:ncol(.)) %>% 
        left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
        group_by(class_name) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(class_name) %>% 
        filter(max_abund >= 0.001) %>% 
        ungroup()
      
      oc.cl.long.rel.meta <- oc.cl.long.rel.meta %>% 
        group_by(Subject.Number, Visit.Name.Norm, class_name) %>% 
        mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
        ungroup()
      
      oc.cl.long.rel.meta <- oc.cl.long.rel.meta %>% 
        mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
        mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
        mutate(absolute_abund_log = (absolute_abund)+0.01)
      
      oc.cl.long.rel.meta <- oc.cl.long.rel.meta %>% 
        filter(!is.na(TRT.norm))
      
      # Add higher taxa assignments to the table
      df.taxa <- oc.abund.meta %>% 
        select(phylum_name, class_name, order_name, family_name, genus_name, Species.corrected) %>% 
        filter(!is.na(Species.corrected)) %>% 
        distinct(Species.corrected, .keep_all = TRUE)
      
      df.taxa.cl <- df.taxa %>% 
        select(phylum_name, class_name) %>% 
        distinct(class_name, .keep_all = TRUE)
      
      oc.cl.long.rel.meta <- left_join(oc.cl.long.rel.meta, df.taxa.cl, by = "class_name")
      oc.cl.long.rel.meta <- left_join(oc.cl.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
      
      oc.cl.long.rel.meta1 <- oc.cl.wide.rel %>% 
        gather(class_name, rel_abund, 2:ncol(.)) %>% 
        left_join(., dnameta, by = "MGS.Sample.ID") %>% 
        group_by(class_name, Subject.Number) %>% 
        mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(class_name) %>% 
        filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
        ungroup()
      
      filtclass <- (oc.cl.long.rel.meta1$class_name %>% unique())
      length(filtclass)  }
    
    ## Generate order-level relative abundance tables
    if(desired_tax_rank == "order") {
      #order level - Based on Read count -  order assignment
      oc.order.wide <- oc.abund.meta %>% 
        filter(superkingdom_name == "Bacteria", tax_rank == "order") %>% 
        mutate(abun_readcount_w_children = readcount_w_children/n_reads_classified) %>%
        select(MGS.Sample.ID, order_name, abun_readcount_w_children) %>% 
        spread(order_name,  abun_readcount_w_children, fill = 0)
      
      MGS.Sample.ID <- oc.order.wide$MGS.Sample.ID
      
      oc.order.wide <- oc.order.wide[-1]
      row.names(oc.order.wide) <- MGS.Sample.ID
      
      ### DECOSTAND NOT REQUIRED AS NORMALIZING TO N CLASSIFIED READS  ####
      oc.order.wide.rel <-  oc.order.wide 
      oc.order.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.order.wide.rel), oc.order.wide.rel)
      
      oc.order.long.rel.meta <- oc.order.wide.rel %>% 
        gather(order_name, rel_abund, 2:ncol(.)) %>% 
        left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
        group_by(order_name) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(order_name) %>% 
        filter(max_abund >= 0.001) %>% 
        ungroup()
      
      oc.order.long.rel.meta <- oc.order.long.rel.meta %>% 
        group_by(Subject.Number, Visit.Name.Norm, order_name) %>% 
        mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
        ungroup()
      
      oc.order.long.rel.meta <- oc.order.long.rel.meta %>% 
        mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
        mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
        mutate(absolute_abund_log = (absolute_abund)+0.01)
      
      oc.order.long.rel.meta <- oc.order.long.rel.meta %>% 
        filter(!is.na(TRT.norm))
      
      # Add higher taxa assignments to the table
      df.taxa.order <- df.taxa %>% 
        select(phylum_name, order_name) %>% 
        distinct(order_name, .keep_all = TRUE)
      
      oc.order.long.rel.meta <- left_join(oc.order.long.rel.meta, df.taxa.order, by = "order_name")
      oc.order.long.rel.meta <- left_join(oc.order.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
      
      oc.order.long.rel.meta1 <- oc.order.wide.rel %>% 
        gather(order_name, rel_abund, 2:ncol(.)) %>% 
        left_join(., dnameta, by = "MGS.Sample.ID") %>% 
        group_by(order_name, Subject.Number) %>% 
        mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(order_name) %>% 
        filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
        ungroup()
      
      filtorder <- (oc.order.long.rel.meta1$order_name %>% unique()) }
    
    ## Generate Phylum-level relative abundance tables
    if(desired_tax_rank == "phylum") {
      
      oc.ph.wide <- oc.abund.meta %>% 
        filter(superkingdom_name == "Bacteria", tax_rank == "phylum") %>% 
        mutate(abun_readcount_w_children = readcount_w_children/n_reads_classified) %>%
        select(MGS.Sample.ID, phylum_name , abun_readcount_w_children) %>%
        spread(phylum_name, abun_readcount_w_children, fill = 0)
      
      MGS.Sample.ID <- oc.ph.wide$MGS.Sample.ID
      oc.ph.wide <- oc.ph.wide %>% select(-MGS.Sample.ID)
      row.names(oc.ph.wide) <- MGS.Sample.ID
      ### DECOSTAND NOT REQUIRED AS NORMALIZING TO N CLASSIFIED READS  #####
      oc.ph.wide.rel <-  oc.ph.wide 
      oc.ph.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.ph.wide.rel), oc.ph.wide.rel)
      
      oc.ph.long.rel.meta <- oc.ph.wide.rel %>% 
        gather(phylum_name, rel_abund, 2:ncol(.)) %>% 
        left_join(., dnameta, by = "MGS.Sample.ID") %>% 
        group_by(phylum_name) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(phylum_name) %>% 
        filter(max_abund >= 0.001) %>% 
        ungroup()
      
      oc.ph.long.rel.meta <- oc.ph.long.rel.meta %>% 
        mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
        mutate(absolute_abund = rel_abund * dna.norm * (2/0.5))
      
      oc.ph.long.rel.meta <- oc.ph.long.rel.meta %>% 
        group_by(Subject.Number, Visit.Name.Norm, phylum_name) %>% 
        mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
        mutate(abs.abund.avg.timepoint = mean(absolute_abund)) %>% 
        ungroup()
      
      oc.ph.long.rel.meta <- oc.ph.long.rel.meta %>% 
        filter(!is.na(TRT.norm))
      
      # Generate a wide format table of average abundance values to calc significance
      oc.ph.wide.avg.meta <- oc.ph.long.rel.meta %>% 
        select(phylum_name, abund.avg.timepoint, Visit.Name.Norm, TRT.norm, Subject.Number) %>% 
        group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>%
        distinct(phylum_name, .keep_all = TRUE) %>% 
        spread(phylum_name, abund.avg.timepoint, fill = 0)
      
      oc.ph.long.avg.meta <- oc.ph.long.rel.meta %>% 
        select(phylum_name, abund.avg.timepoint, Visit.Name.Norm, TRT.norm, Subject.Number, MGS.Sample.ID, rec.diagnosis, Day.in.treatment) %>% 
        group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>%
        distinct(phylum_name, .keep_all = TRUE)
      
      oc.ph.long.avg.meta <- left_join(oc.ph.long.avg.meta, ve303.sp.sum, by = "MGS.Sample.ID")
      filtphylum <- (oc.ph.long.rel.meta$phylum_name %>% unique())
      
      # Calculate total biomass per subject over time
      biomass.per.subj <- oc.ph.long.rel.meta %>%
        left_join(.,select(oc.abund.meta, MGS.Sample.ID, n_reads), by="MGS.Sample.ID") %>%
        group_by(Subject.Number, Visit.Name.Norm) %>% 
        mutate(tot.biomass = sum(absolute_abund)) %>% 
        group_by(Subject.Number, Visit.Name.Norm) %>% 
        mutate(biomass.tp.mean = mean(tot.biomass), biomass.tp.sd = sd(tot.biomass), biomass.tp.median = median(tot.biomass)) %>% 
        group_by(Subject.Number) %>% 
        distinct(Visit.Name.Norm, .keep_all = TRUE) %>% 
        select(Subject.Number, TRT.norm, Visit.Name.Norm, Day.in.treatment, tot.biomass, biomass.tp.mean, biomass.tp.sd, biomass.tp.median, n_reads, weight.of.collected.sample..g.) %>% 
        arrange(Subject.Number, Visit.Name.Norm) %>% 
        group_by(TRT.norm, Visit.Name.Norm) %>% 
        mutate(cohort.med.tp = median(biomass.tp.median)) %>% 
        group_by(Visit.Name.Norm) %>% 
        mutate(overall.med.tp = median(biomass.tp.median)) %>% 
        ungroup()  }
    
    ## Generate Clostridia, Bacilli, Gammaproteo relative abundance tables
    if(desired_tax_rank == "Clos_Bac_Gam"){
      oc.abund.meta$abundance_w_children[is.na(oc.abund.meta$abundance_w_children)] <- 0
      
      oc.abund.long.sp <- oc.abund.meta %>% 
        filter(superkingdom_name == "Bacteria", tax_rank == "species") %>% 
        group_by(sample_name) %>%
        distinct(Species.corrected, .keep_all = TRUE) %>% 
        group_by(Species.corrected) %>% 
        mutate(max_abund = max(abundance_w_children), est.abund.per = (abundance_w_children)*100) %>%
        ungroup()
      
      clost.top.sp.list <- oc.abund.long.sp %>% 
        filter(abundance_w_children >0) %>%
        filter(grepl("Clostridia", class_name)) %>%
        group_by(sample_id) %>% 
        top_n(n=2, wt = abundance_w_children) %>% 
        ungroup() %>% 
        distinct(Species.corrected)
      
      clost.top.sp <- oc.abund.long.sp %>% 
        filter(Species.corrected %in% clost.top.sp.list$Species.corrected)
      
      gtdb.genus.sum <- oc.abund.long.sp %>% 
        group_by(MGS.Sample.ID, phylum_name, genus_name, Visit.Name.Norm) %>% 
        summarise(total_abund = (sum(abundance_w_children)*100)) %>%
        select(-Visit.Name.Norm) %>%
        left_join(., dnameta, by = "MGS.Sample.ID")
      
      gtdb.gen.sum.clost <- gtdb.genus.sum %>% 
        filter(genus_name %in% clost.top.sp$genus_name) %>%
        mutate(total_abund_log = total_abund + 0.01)
      
      #######################################################
      
      bac.top.sp.list <- oc.abund.long.sp %>% 
        filter(abundance_w_children >0) %>%
        filter(grepl("Bacilli", class_name)) %>%
        group_by(sample_id) %>% 
        top_n(n=2, wt = abundance_w_children) %>% 
        ungroup() %>% 
        distinct(Species.corrected)
      
      bac.top.sp <- oc.abund.long.sp %>% 
        filter(Species.corrected %in% bac.top.sp.list$Species.corrected)
      
      gtdb.gen.sum.bac <- gtdb.genus.sum %>% 
        filter(genus_name %in% bac.top.sp$genus_name) %>%
        mutate(total_abund_log = total_abund + 0.01)
      ######################
      
      gprot.top.sp.list <- oc.abund.long.sp %>% 
        filter(abundance_w_children >0) %>%
        filter(grepl("Gammaproteobacteria", class_name)) %>%
        group_by(sample_id) %>% 
        top_n(n=2, wt = abundance_w_children) %>% 
        ungroup() %>% 
        distinct(Species.corrected)
      
      gprot.top.sp <- oc.abund.long.sp %>% 
        filter(Species.corrected %in% gprot.top.sp.list$Species.corrected)
      
      gtdb.gen.sum.gprot <- gtdb.genus.sum %>% 
        filter(genus_name %in% gprot.top.sp$genus_name) %>%
        mutate(total_abund_log = total_abund + 0.01) }
    ################################## 
  }
  
  if(abun_choice =="abun_w_children"){
    if(desired_tax_rank == "species") {
      #Species level abundance summary
      oc.sp.wide <- oc.abund.meta %>% 
        filter(superkingdom_name == "Bacteria", tax_rank == "species") %>% 
        select(MGS.Sample.ID, Species.corrected, abundance_w_children) %>%
        group_by(MGS.Sample.ID, Species.corrected) %>% 
        summarise_at(vars(abundance_w_children), funs(sum)) %>%
        ungroup() %>% 
        spread(Species.corrected, abundance_w_children, fill = 0)
      
      MGS.Sample.ID <- oc.sp.wide$MGS.Sample.ID
      oc.sp.wide <- oc.sp.wide %>% 
        select(-MGS.Sample.ID)
      
      row.names(oc.sp.wide) <- MGS.Sample.ID
      
      oc.sp.wide.rel <- oc.sp.wide
      oc.sp.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.sp.wide.rel), oc.sp.wide.rel)
      
      oc.sp.long.rel.meta <- oc.sp.wide.rel %>% 
        gather(Species.corrected, rel_abund, 2:ncol(.)) %>% 
        left_join(., dnameta, by = "MGS.Sample.ID") %>% 
        group_by(Species.corrected) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        mutate(min_abund = min(rel_abund)) %>% 
        arrange(Species.corrected) %>% 
        filter(max_abund >= 0.001) %>% 
        ungroup()
      
      oc.sp.long.rel.meta <- oc.sp.long.rel.meta %>% 
        group_by(Subject.Number, Visit.Name.Norm, Species.corrected) %>% 
        mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
        ungroup()
      
      oc.sp.long.rel.meta <- oc.sp.long.rel.meta %>% 
        mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
        mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
        mutate(absolute_abund_log = (absolute_abund)+0.01)
      
      oc.sp.long.rel.meta <- oc.sp.long.rel.meta %>% 
        filter(!is.na(TRT.norm))
      
      oc.sp.long.rel.meta.ve303 <- oc.sp.long.rel.meta %>% 
        filter(grepl("VE303", Species.corrected)) 
      
      ve303.sp.sum <- oc.sp.long.rel.meta.ve303 %>% 
        group_by(MGS.Sample.ID) %>% 
        summarise(ve303_tot_perc = (sum(rel_abund) *100))
      
      oc.sp.long.rel.meta <- left_join(oc.sp.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
      oc.sp.long.rel.meta.ve303 <- left_join(oc.sp.long.rel.meta.ve303, ve303.sp.sum, by = "MGS.Sample.ID")
      
      oc.sp.long.rel.meta1 <- oc.sp.wide.rel %>% 
        gather(Species.corrected, rel_abund, 2:ncol(.)) %>% 
        left_join(., dnameta, by = "MGS.Sample.ID") %>% 
        group_by(Species.corrected, Subject.Number) %>% 
        mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        mutate(mean_abund = mean(rel_abund)) %>% 
        arrange(Species.corrected) %>% 
        #filter(max_abund >= 0.005 & binary_abund >= 0.5 & mean_abund >= 0.001) %>% 
        filter(max_abund >= threshabun & binary_abund >= threshprev & mean_abund >= threshmean) %>% 
        ungroup()
      
      filtspecies <- (oc.sp.long.rel.meta1$Species.corrected %>% unique())
      allspecies <- (oc.sp.long.rel.meta$Species.corrected %>% unique()) }
    
    
    ## Generate genus-level relative abundance tables
    if(desired_tax_rank == "genus") {
      
      #genus level - Based on Read count -  OR NCBI genus assignment
      oc.genus.wide <- oc.abund.meta %>%
        filter(superkingdom_name == "Bacteria", tax_rank == "genus") %>% 
        select(MGS.Sample.ID, genus_name, abundance_w_children) %>% 
        spread(genus_name, abundance_w_children, fill = 0)
      
      MGS.Sample.ID <- oc.genus.wide$MGS.Sample.ID
      oc.genus.wide <- oc.genus.wide[-1]
      row.names(oc.genus.wide) <- MGS.Sample.ID
      
      oc.genus.wide.rel <- oc.genus.wide 
      oc.genus.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.genus.wide.rel), oc.genus.wide.rel)
      
      oc.genus.long.rel.meta <- oc.genus.wide.rel %>% 
        gather(genus_name, rel_abund, 2:ncol(.)) %>% 
        left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
        group_by(genus_name) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(genus_name) %>% 
        filter(max_abund >= 0.001) %>% 
        ungroup()
      
      oc.genus.long.rel.meta <- oc.genus.long.rel.meta %>% 
        group_by(Subject.Number, Visit.Name.Norm, genus_name) %>% 
        mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
        ungroup()
      
      oc.genus.long.rel.meta <- oc.genus.long.rel.meta %>% 
        mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
        mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
        mutate(absolute_abund_log = (absolute_abund)+0.01)
      
      oc.genus.long.rel.meta <- oc.genus.long.rel.meta %>% 
        filter(!is.na(TRT.norm))
      
      # Add higher taxa assignments to the table
      df.taxa <- oc.abund.meta %>% 
        select(phylum_name, genus_name, order_name, family_name, class_name, Species.corrected) %>% 
        filter(!is.na(Species.corrected)) %>% 
        distinct(Species.corrected, .keep_all = TRUE)
      
      df.taxa.genus <- df.taxa %>% 
        select(phylum_name, genus_name) %>% 
        distinct(genus_name, .keep_all = TRUE)
      
      oc.genus.long.rel.meta <- left_join(oc.genus.long.rel.meta, df.taxa.genus, by = "genus_name")
      oc.genus.long.rel.meta <- left_join(oc.genus.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
      
      oc.genus.long.rel.meta1 <- oc.genus.wide.rel %>% 
        gather(genus_name, rel_abund, 2:ncol(.)) %>% 
        left_join(., dnameta, by = "MGS.Sample.ID") %>% 
        group_by(genus_name, Subject.Number) %>% 
        mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(genus_name) %>% 
        filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
        ungroup()
      
      filtgenus <- (oc.genus.long.rel.meta1$genus_name %>% unique()) }
    
    ## Generate family-level relative abundance tables
    if(desired_tax_rank == "family") {
      #family level - Based on abundance -  family assignment
      oc.family.wide <- oc.abund.meta %>%
        filter(superkingdom_name == "Bacteria", tax_rank == "family") %>% 
        select(MGS.Sample.ID, tax_name, abundance_w_children) %>% 
        mutate(family_name = tax_name) %>%
        select(-tax_name) %>%
        spread(family_name, abundance_w_children, fill = 0)
      
      MGS.Sample.ID <- oc.family.wide$MGS.Sample.ID
      oc.family.wide <- oc.family.wide[-1]
      row.names(oc.family.wide) <- MGS.Sample.ID
      
      oc.family.wide.rel <- oc.family.wide 
      oc.family.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.family.wide.rel), oc.family.wide.rel)
      
      oc.family.long.rel.meta <- oc.family.wide.rel %>% 
        gather(family_name, rel_abund, 2:ncol(.)) %>% 
        left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
        group_by(family_name) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(family_name) %>% 
        filter(max_abund >= 0.001) %>% 
        ungroup()
      
      oc.family.long.rel.meta <- oc.family.long.rel.meta %>% 
        group_by(Subject.Number, Visit.Name.Norm, family_name) %>% 
        mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
        ungroup()
      
      oc.family.long.rel.meta <- oc.family.long.rel.meta %>% 
        mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
        mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
        mutate(absolute_abund_log = (absolute_abund)+0.01)
      
      oc.family.long.rel.meta <- oc.family.long.rel.meta %>% 
        filter(!is.na(TRT.norm))
      
      # Add higher taxa assignments to the table
      df.taxa <- oc.abund.meta %>% 
        select(phylum_name, genus_name, order_name, family_name, class_name, Species.corrected) %>% 
        filter(!is.na(Species.corrected)) %>% 
        distinct(Species.corrected, .keep_all = TRUE)
      
      df.taxa.family <- df.taxa %>% 
        select(phylum_name, family_name) %>% 
        distinct(family_name, .keep_all = TRUE)
      
      oc.family.long.rel.meta <- left_join(oc.family.long.rel.meta, df.taxa.family, by = "family_name")
      oc.family.long.rel.meta <- left_join(oc.family.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
      
      oc.family.long.rel.meta1 <- oc.family.wide.rel %>% 
        gather(family_name, rel_abund, 2:ncol(.)) %>% 
        left_join(., dnameta, by = "MGS.Sample.ID") %>% 
        group_by(family_name, Subject.Number) %>% 
        mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(family_name) %>% 
        filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
        ungroup()
      
      filtfamily <- (oc.family.long.rel.meta1$family_name %>% unique()) }
    
    ## Generate class-level relative abundance tables
    if(desired_tax_rank == "class") {
      #Class level - Based on abun -  class assignment
      oc.cl.wide <- oc.abund.meta %>% 
        filter(superkingdom_name == "Bacteria", tax_rank == "class") %>% 
        select(MGS.Sample.ID, class_name, abundance_w_children) %>% 
        spread(class_name, abundance_w_children, fill = 0)
      
      MGS.Sample.ID <- oc.cl.wide$MGS.Sample.ID
      
      oc.cl.wide <- oc.cl.wide[-1]
      row.names(oc.cl.wide) <- MGS.Sample.ID
      
      oc.cl.wide.rel <-  oc.cl.wide 
      oc.cl.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.cl.wide.rel), oc.cl.wide.rel)
      
      oc.cl.long.rel.meta <- oc.cl.wide.rel %>% 
        gather(class_name, rel_abund, 2:ncol(.)) %>% 
        left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
        group_by(class_name) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(class_name) %>% 
        filter(max_abund >= 0.001) %>% 
        ungroup()
      
      oc.cl.long.rel.meta <- oc.cl.long.rel.meta %>% 
        group_by(Subject.Number, Visit.Name.Norm, class_name) %>% 
        mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
        ungroup()
      
      oc.cl.long.rel.meta <- oc.cl.long.rel.meta %>% 
        mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
        mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
        mutate(absolute_abund_log = (absolute_abund)+0.01)
      
      oc.cl.long.rel.meta <- oc.cl.long.rel.meta %>% 
        filter(!is.na(TRT.norm))
      
      # Add higher taxa assignments to the table
      df.taxa <- oc.abund.meta %>% 
        select(phylum_name, class_name, order_name, family_name, genus_name, Species.corrected) %>% 
        filter(!is.na(Species.corrected)) %>% 
        distinct(Species.corrected, .keep_all = TRUE)
      
      df.taxa.cl <- df.taxa %>% 
        select(phylum_name, class_name) %>% 
        distinct(class_name, .keep_all = TRUE)
      
      oc.cl.long.rel.meta <- left_join(oc.cl.long.rel.meta, df.taxa.cl, by = "class_name")
      oc.cl.long.rel.meta <- left_join(oc.cl.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
      
      
      oc.cl.long.rel.meta1 <- oc.cl.wide.rel %>% 
        gather(class_name, rel_abund, 2:ncol(.)) %>% 
        left_join(., dnameta, by = "MGS.Sample.ID") %>% 
        group_by(class_name, Subject.Number) %>% 
        mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(class_name) %>% 
        filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
        ungroup()
      
      filtclass <- (oc.cl.long.rel.meta1$class_name %>% unique()) }
    
    ## Generate order-level relative abundance tables
    if(desired_tax_rank == "order") {
      
      #order level - Based on abun - order assignment
      oc.order.wide <- oc.abund.meta %>% 
        filter(superkingdom_name == "Bacteria", tax_rank =="order") %>% 
        select(MGS.Sample.ID, order_name,abundance_w_children) %>% 
        spread(order_name, abundance_w_children, fill = 0)
      
      MGS.Sample.ID <- oc.order.wide$MGS.Sample.ID
      oc.order.wide <- oc.order.wide[-1]
      row.names(oc.order.wide) <- MGS.Sample.ID
      
      oc.order.wide.rel <-  oc.order.wide
      oc.order.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.order.wide.rel), oc.order.wide.rel)
      
      oc.order.long.rel.meta <- oc.order.wide.rel %>% 
        gather(order_name, rel_abund, 2:ncol(.)) %>% 
        left_join(.,dnameta, by = "MGS.Sample.ID") %>% 
        group_by(order_name) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(order_name) %>% 
        filter(max_abund >= 0.001) %>% 
        ungroup()
      
      oc.order.long.rel.meta <- oc.order.long.rel.meta %>% 
        group_by(Subject.Number, Visit.Name.Norm, order_name) %>% 
        mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
        ungroup()
      
      oc.order.long.rel.meta <- oc.order.long.rel.meta %>% 
        mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
        mutate(absolute_abund = rel_abund * dna.norm * (2/0.5)) %>%
        mutate(absolute_abund_log = (absolute_abund)+0.01)
      
      oc.order.long.rel.meta <- oc.order.long.rel.meta %>% 
        filter(!is.na(TRT.norm))
      
      # Add higher taxa assignments to the table
      
      df.taxa.order <- df.taxa %>% 
        select(phylum_name, order_name) %>% 
        distinct(order_name, .keep_all = TRUE)
      
      oc.order.long.rel.meta <- left_join(oc.order.long.rel.meta, df.taxa.order, by = "order_name")
      oc.order.long.rel.meta <- left_join(oc.order.long.rel.meta, ve303.sp.sum, by = "MGS.Sample.ID")
      
      
      oc.order.long.rel.meta1 <- oc.order.wide.rel %>% 
        gather(order_name, rel_abund, 2:ncol(.)) %>% 
        left_join(., dnameta, by = "MGS.Sample.ID") %>% 
        group_by(order_name, Subject.Number) %>% 
        mutate(binary_abund = sum(sign(rel_abund))/sum(sign(rel_abund+0.000001))) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(order_name) %>% 
        filter(max_abund >= threshabun & binary_abund >= threshprev) %>% 
        ungroup()
      
      filtorder <- (oc.order.long.rel.meta1$order_name %>% unique()) }
    
    
    ## Generate Phylum-level relative abundance tables
    if(desired_tax_rank == "phylum") {
      
      oc.ph.wide <- oc.abund.meta %>% 
        filter(superkingdom_name == "Bacteria", tax_rank == "phylum") %>% 
        select(MGS.Sample.ID, phylum_name, abundance_w_children) %>%
        spread(phylum_name, abundance_w_children, fill = 0) 
      
      MGS.Sample.ID <- oc.ph.wide$MGS.Sample.ID
      oc.ph.wide <- oc.ph.wide %>% 
        select(-MGS.Sample.ID)
      
      row.names(oc.ph.wide) <- MGS.Sample.ID
      
      oc.ph.wide.rel <-  oc.ph.wide  
      oc.ph.wide.rel <- data.frame(MGS.Sample.ID = row.names(oc.ph.wide.rel), oc.ph.wide.rel)
      
      oc.ph.long.rel.meta <- oc.ph.wide.rel %>% 
        gather(phylum_name, rel_abund, 2:ncol(.)) %>% 
        left_join(., dnameta, by = "MGS.Sample.ID") %>% 
        group_by(phylum_name) %>% 
        mutate(max_abund = max(rel_abund)) %>% 
        arrange(phylum_name) %>% 
        filter(max_abund >= 0.001) %>% 
        ungroup()
      
      oc.ph.long.rel.meta <- oc.ph.long.rel.meta %>% 
        mutate(dna.norm = (DNA.yield.ug/1e6) / weight.of.collected.sample..g.) %>% 
        mutate(absolute_abund = rel_abund * dna.norm * (2/0.5))
      
      oc.ph.long.rel.meta <- oc.ph.long.rel.meta %>% 
        group_by(Subject.Number, Visit.Name.Norm, phylum_name) %>% 
        mutate(abund.avg.timepoint = (mean(rel_abund)*100)) %>% 
        mutate(abs.abund.avg.timepoint = mean(absolute_abund)) %>% 
        ungroup()
      
      oc.ph.long.rel.meta <- oc.ph.long.rel.meta %>% 
        filter(!is.na(TRT.norm))
      
      # Generate a wide format table of average abundance values to calc significance
      oc.ph.wide.avg.meta <- oc.ph.long.rel.meta %>% 
        select(phylum_name, abund.avg.timepoint, Visit.Name.Norm, TRT.norm, Subject.Number) %>% 
        group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>%
        distinct(phylum_name, .keep_all = TRUE) %>% 
        spread(phylum_name, abund.avg.timepoint, fill = 0)
      
      oc.ph.long.avg.meta <- oc.ph.long.rel.meta %>% 
        select(phylum_name, abund.avg.timepoint, Visit.Name.Norm, TRT.norm, Subject.Number, MGS.Sample.ID, rec.diagnosis, Day.in.treatment) %>% 
        group_by(Subject.Number, TRT.norm, Visit.Name.Norm) %>%
        distinct(phylum_name, .keep_all = TRUE)
      
      oc.ph.long.avg.meta <- left_join(oc.ph.long.avg.meta, ve303.sp.sum, by = "MGS.Sample.ID")
      filtphylum <- (oc.ph.long.rel.meta$phylum_name %>% unique())
      
      # Calculate total biomass per subject over time
      biomass.per.subj <- oc.ph.long.rel.meta %>%
        left_join(.,select(oc.abund.meta, MGS.Sample.ID, n_reads), by="MGS.Sample.ID") %>%
        group_by(Subject.Number, Visit.Name.Norm) %>% 
        mutate(tot.biomass = sum(absolute_abund)) %>% 
        group_by(Subject.Number, Visit.Name.Norm) %>% 
        mutate(biomass.tp.mean = mean(tot.biomass), biomass.tp.sd = sd(tot.biomass), biomass.tp.median = median(tot.biomass)) %>% 
        group_by(Subject.Number) %>% 
        distinct(Visit.Name.Norm, .keep_all = TRUE) %>% 
        select(Subject.Number, TRT.norm, Visit.Name.Norm, Day.in.treatment, tot.biomass, biomass.tp.mean, biomass.tp.sd, biomass.tp.median, n_reads, weight.of.collected.sample..g.) %>% 
        arrange(Subject.Number, Visit.Name.Norm) %>% 
        group_by(TRT.norm, Visit.Name.Norm) %>% 
        mutate(cohort.med.tp = median(biomass.tp.median)) %>% 
        group_by(Visit.Name.Norm) %>% 
        mutate(overall.med.tp = median(biomass.tp.median)) %>% 
        ungroup() }
    
    ## Generate Clostridia, Bacilli, Gammaproteo relative abundance tables
    if(desired_tax_rank == "Clos_Bac_Gam"){
      oc.abund.meta$abundance_w_children[is.na(oc.abund.meta$abundance_w_children)] <- 0
      oc.abund.long.sp <- oc.abund.meta %>% 
        filter(superkingdom_name == "Bacteria", tax_rank == "species") %>% 
        group_by(sample_name) %>%
        distinct(Species.corrected, .keep_all = TRUE) %>% 
        group_by(Species.corrected) %>% 
        mutate(max_abund = max(abundance_w_children), est.abund.per = (abundance_w_children)*100) %>%
        ungroup()
      
      clost.top.sp.list <- oc.abund.long.sp %>% 
        filter(abundance_w_children >0) %>%
        filter(grepl("Clostridia", class_name)) %>%
        group_by(sample_id) %>% 
        top_n(n=2, wt = abundance_w_children) %>% 
        ungroup() %>% 
        distinct(Species.corrected)
      
      clost.top.sp <- oc.abund.long.sp %>% 
        filter(Species.corrected %in% clost.top.sp.list$Species.corrected)
      
      gtdb.genus.sum <- oc.abund.long.sp %>% 
        group_by(MGS.Sample.ID, phylum_name, genus_name, Visit.Name.Norm) %>% 
        summarise(total_abund = (sum(abundance_w_children)*100)) %>%
        select(-Visit.Name.Norm) %>%
        left_join(., dnameta, by = "MGS.Sample.ID")
      
      gtdb.gen.sum.clost <- gtdb.genus.sum %>% 
        filter(genus_name %in% clost.top.sp$genus_name) %>%
        mutate(total_abund_log = total_abund + 0.01)
      
      #######################################################
      
      bac.top.sp.list <- oc.abund.long.sp %>% 
        filter(abundance_w_children >0) %>%
        filter(grepl("Bacilli", class_name)) %>%
        group_by(sample_id) %>% 
        top_n(n=2, wt = abundance_w_children) %>% 
        ungroup() %>% 
        distinct(Species.corrected)
      
      bac.top.sp <- oc.abund.long.sp %>% 
        filter(Species.corrected %in% bac.top.sp.list$Species.corrected)
      
      gtdb.gen.sum.bac <- gtdb.genus.sum %>% 
        filter(genus_name %in% bac.top.sp$genus_name) %>%
        mutate(total_abund_log = total_abund + 0.01)
      ######################
      
      gprot.top.sp.list <- oc.abund.long.sp %>% 
        filter(abundance_w_children >0) %>%
        filter(grepl("Gammaproteobacteria", class_name)) %>%
        group_by(sample_id) %>% 
        top_n(n=2, wt = abundance_w_children) %>% 
        ungroup() %>% 
        distinct(Species.corrected)
      
      gprot.top.sp <- oc.abund.long.sp %>% 
        filter(Species.corrected %in% gprot.top.sp.list$Species.corrected)
      
      gtdb.gen.sum.gprot <- gtdb.genus.sum %>% 
        filter(genus_name %in% gprot.top.sp$genus_name) %>%
        mutate(total_abund_log = total_abund + 0.01) }
    ##################################
  }
}


