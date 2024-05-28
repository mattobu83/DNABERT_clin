check_same_sign <- function(column) {
  # Remove NA values
  non_na_values <- column[!is.na(column)]
  
  # Get the sign of non-NA values (+1 for positive, -1 for negative, 0 for zero)
  signs <- sign(non_na_values)
  
  # Check if all signs are the same
  all_same_sign <- all(signs == signs[1])
  
  return(all_same_sign)
}

check_same_value <- function(column) {
  # Remove NA values
  non_na_values <- column[!is.na(column)]
  if (length(non_na_values) == 0) {
    all_same_value = TRUE #since all of the values were NA 
  } else {
    all_same_value <- length(unique(non_na_values)) == 1
  }
  return(all_same_value)
}

prefilter_log_diff <- function(log_diff_path,variant_caller,percentile = FALSE){
  df <- read.delim(log_diff_path,row.names = 1)
  df <- round(df,2)
  muts <- df[grepl(variant_caller,rownames(df)),] # determine the variant caller
  
  #drop any columns with all NAs 
  muts <- muts[,colSums(is.na(muts))<nrow(muts)]
  
  if (variant_caller == "pindel"){
    rownames(muts) <- gsub("_sanger_raw_pindel","",rownames(muts))
  }else{
    rownames(muts) <- gsub("_CaVEMan","",rownames(muts))
  }
  
  if (percentile == TRUE){
    #zscore and filter for 95th percentile 
    mut_info <- data.frame(coord = colnames(muts), ld = rowMeans(t(muts),na.rm = T))
    mut_info$zscore <- (mut_info$ld - mean(mut_info$ld) / sd(mut_info$ld))
    z95 <- mut_info[,mut_info$zscore >= 1.645]
    filt_muts <- muts[,z95]
  } else {
    filt_muts <- muts
  }

  
  return(filt_muts)
}

unicox_mut <- function(surv_mut){
  #remove patients with no survival data or mutation data 
  surv_mut <- surv_mut %>% drop_na()
  
  #do univariate CPH to find significant mutations
  covariates <-colnames(surv_mut)[4:length(colnames(surv_mut))]
  print("last feature column:")
  print(colnames(surv_mut)[3])
  print("first mutation column:")
  print(colnames(surv_mut)[4])
  univ_formulas <- sapply(covariates,
                          function(x) as.formula(paste('Surv(days, status)~', (x))))
  
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = surv_mut)})
  # Extract data 
  univ_results <- lapply(univ_models,
                         function(x){ 
                           x <- summary(x)
                           p.value<-signif(x$wald["pvalue"], digits=2)
                           wald.test<-signif(x$wald["test"], digits=2)
                           beta<-signif(x$coef[1], digits=2);#coeficient beta
                           HR <-signif(x$coef[2], digits=2);#exp(beta)
                           HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                           HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                           #CONF <- paste0("(", 
                           #               HR.confint.lower, "-", HR.confint.upper, ")")
                           res<-c(beta, HR, HR.confint.lower,HR.confint.upper,wald.test, p.value)
                           names(res)<-c("beta", "HR", "95% CI low","95% CI high", "wald.test", 
                                         "p.value")
                           return(res)
                           #return(exp(cbind(coef(x),confint(x))))
                         })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  res_df <- as.data.frame(res)
  res_df$padj <- p.adjust(res_df$p.value,method = "BH")
  print("Completed Univariate CoxPH. Dimensions of results table:")
  print(dim(res_df))
  return(res_df)
}


run_cph <- function(merged_df,surv_mut,mut_count,t2g,GDC_project, variant_caller,n_mut_filter,percent_patients){ #need merged df as well to get the log difference
  cre <- unlist(strsplit(getwd(),"/"))[length(unlist(strsplit(getwd(),"/")))]
  print(paste0("Running: ", cre))
  #run cph and filter for significant only 
  res_df <- unicox_mut(surv_mut)
  res_df <- res_df[!is.na(res_df$p.value),]
  res_df_counted <- merge(mut_count,res_df, by = 0 ) #add patient counts 
  
  print(paste0("Filtering for muts present in at least ",n_mut_filter ," patients, with a p-value < 0.05" ))
  sig_df <- res_df_counted[res_df_counted$p.value < 0.05 & res_df_counted$count >= n_mut_filter,]
  
  if(dim(sig_df)[1] == 0){
    print("no patients passed filtering, returning mutations CPH results")
    sig_df_counted <- res_df_counted
    output_filename <- paste0("univariate_coxph_results_", cre,"_" ,variant_caller,"_",GDC_project,"_nonNA_pval.tsv")
  } else {
    sig_df_counted <- sig_df
    output_filename <- paste0("univariate_coxph_results_", cre,"_" ,variant_caller,"_",GDC_project,"_p05_mut",percent_patients*100, "p.tsv")
  }
  
  
  sig_df_counted_ord <-  sig_df_counted[order(sig_df_counted$HR,decreasing = T),]
  sig_mut_info <- sig_df_counted_ord[,c('mut','count',"HR","p.value","padj")]
  sig_mut_ld <- merged_df[,sig_mut_info$mut,drop=FALSE] #the corresponding log differences
  sig_mut_info$same_value <- sapply(sig_mut_ld, check_same_value)
  sig_mut_info$same_sign <- sapply(sig_mut_ld, check_same_sign)
  sig_mut_info$mean_ld <- rowMeans(t(sig_mut_ld),na.rm = T)
  sig_mut_info$ensg <- apply(sig_mut_info,MARGIN = 1,get_gene_name,t2g)
  sig_mut_info$hgnc <- sapply(sig_mut_info$ensg,get_hugo_name,t2g)
  print("Significant mutations: ")
  print(dim(sig_mut_info))
  print(head(sig_mut_info))
  
  write.table(sig_mut_info,output_filename,sep = "\t",quote = F, row.names = F)
  
  return(sig_mut_info)
}

plot_km <- function(surv_mut_filename,sig_mut_df,mut_of_interest){
  surv_mut <- read.delim(surv_mut_filename,row.names = 1)
  hugo <- sig_mut_df[sig_mut_df$mut == mut_of_interest,"hgnc"]
  univ_formula <- as.formula(paste('Surv(days, status)~', (mut_of_interest)))
  x <- coxph(univ_formula, data = surv_mut)
  coxph_sum <- summary(x)
  fit <- survminer::surv_fit(univ_formula, data = surv_mut)
  p.value<-signif(coxph_sum$wald["pvalue"], digits=2)
  kmplot <- ggsurvplot(fit, risk.table = TRUE,pval = p.value,title = paste("KM Plot: ",hugo))
  print(kmplot)
}

process_log_diff_surv <- function(mutation_file,
                             clinical_file,
                             mapping_file,
                             variant_caller, t2g = t2g, #mart
                             GDC_project = "all",
                             percent_patients = 0.1, percentile = FALSE){


  print("Reading in and preparing Files")
  mapping = read.delim(mapping_file)
  muts <- prefilter_log_diff(mutation_file,variant_caller,percentile)
  
  merged_df <- merge(mapping[,c('prefix','patient')],muts,by.x = 1, by.y = 0)
  
  #process survival data 
  clin <- read.delim(clinical_file)
  table(clin$'case_submitter_id' %in% mapping$patient)
  surv <- clin[,c("case_submitter_id",'vital_status', 'days_to_death', 'days_to_last_follow_up','project_id')]
  surv <- surv[!duplicated(surv$case_submitter_id),]
  surv$days <- ifelse(surv$days_to_death == "'--", as.numeric(surv$days_to_last_follow_up), as.numeric(surv$days_to_death))
  surv$status <- ifelse(surv$vital_status == "Alive",0,1)
  
  #remove patients that have -ve numbers and NAs
  surv <- surv[!is.na(surv$days),]
  surv <- surv[!surv$days < 0, c("case_submitter_id", 'project_id','days','status')]
  print("Dimensions of survival dataframe")
  print(dim(surv))
  
  merged_df <- merged_df[!merged_df$patient %in% c("1105273","C3N-02788"),]
  rownames(merged_df) <- merged_df$patient
  merged_df <- merged_df[,-c(1,2)]
  print("Dimensions of mapped mutation dataframe")
  print(dim(merged_df))
  
  
  #get binary matrix
  merged_df_bin <- merged_df %>%
    mutate(across(where(is.numeric), ~ifelse(is.na(.), 0, 1))  )
  
  surv_mut <- merge(surv,merged_df_bin,by.x = 1, by.y = 0) #merge with survival

  rownames(surv_mut) <- surv_mut[,1]
  surv_mut <- surv_mut[-1]
  print("Dimensions of merged surival and mutation dataframe.")
  print(dim(surv_mut))
  

  
  #filter for project  
  if(GDC_project != "all"){
    print("Filtering for GDC project:")
    if(GDC_project == "GBM-only"){
      surv_mut_filtered <- surv_mut[surv_mut$project_id == "CPTAC-3" | surv_mut$project_id == "TCGA-GBM",]
      merged_df_bin_filtered <- merged_df_bin[rownames(surv_mut_filtered),]
      merged_df_filtered <- merged_df[rownames(surv_mut_filtered),]
      print(table(surv_mut_filtered$project_id))
    } else {
      surv_mut_filtered <- surv_mut[surv_mut$project_id == GDC_project,]
      merged_df_bin_filtered <- merged_df_bin[rownames(surv_mut_filtered),]
      merged_df_filtered <- merged_df[rownames(surv_mut_filtered),]
      print(table(surv_mut_filtered$project_id))
    }
    

  }else{
    print("No project selected, processing all patients in dataframe ")
    surv_mut_filtered <- surv_mut
    merged_df_bin_filtered <- merged_df_bin
    merged_df_filtered <- merged_df 
  }

  print("Dimensions after project selection:")
  print(dim(merged_df_bin_filtered))
  
  
  #output survival mutation , do not overwrite if the file exists already
  surv_mut_filename <- paste0("surv_mut_",variant_caller,"_",GDC_project,".tsv")
  if (!file.exists(surv_mut_filename)){
    print("Writing merged survival and mutation dataframe to file")
    write.table(surv_mut_filtered,surv_mut_filename,sep = "\t",row.names = T,quote = F)
  } else {
    print("surv mut for current combo already written. Skipping.")
  }
  
  #get number of mutants
  mut_count <- data.frame('mut' = colnames(merged_df_bin_filtered),'count' = rowSums(t(merged_df_bin_filtered)))
  n_mut_filter <- round(length(rownames(merged_df_bin_filtered))*percent_patients)
  print("Running CoxPH, returning significant mutatnts")
  sig_mut <- run_cph(merged_df_filtered,surv_mut_filtered,mut_count,t2g,GDC_project,variant_caller,n_mut_filter,percent_patients)
  print("Summary of Significant Mutants:")
  print("Are all log diffs in column the same value?")
  print(table(sig_mut$same_value))
  print("Are all log diffs in column the same sign?")
  print(table(sig_mut$same_sign))
  print("Disruption type?")
  print(table(sign(sig_mut$mean_ld)))
  print("Significant after FDR adjustment?")
  print(table(sig_mut$padj < 0.05))
  return(sig_mut)
  
}



