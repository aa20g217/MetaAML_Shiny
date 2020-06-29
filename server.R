# metaAML server.R
#
# Brooks Benard
# bbenard@stanford.edu
#
#

### This script compiles the mutation calls across four different cohorts of AML patients (Stanford, TCGA, BeatAML, and MultiStage) and analyses the co-occurence, mutual exclusive, and survival patterns for different mutation pairs. It allows the user to perform subset analyses on de novo, secondary, and repalse patient cohorts, in addition to subsetting analysis to desired gene sets.

#### load required packages and data ####
if (!require('shiny')) install.packages('shiny'); library('shiny')
if (!require('scales')) install.packages('scales'); library('scales')
if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('readxl')) install.packages('readxl'); library('readxl')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('reshape2')) install.packages('reshape2'); library('reshape2')
if (!require('plyr')) install.packages('plyr'); library('plyr')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('UpSetR')) install.packages('UpSetR'); library('UpSetR')
if (!require('muhaz')) install.packages('muhaz'); library('muhaz')
if (!require('data.table')) install.packages('data.table'); library('data.table')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('ggrepel')) install.packages('ggrepel'); library('ggrepel')
if (!require('ggpubr')) install.packages('ggpubr'); library('ggpubr')
if (!require('RCurl')) install.packages('RCurl'); library('RCurl')
if (!require('muhaz')) install.packages('muhaz'); library('muhaz')
if (!require('survival')) install.packages('survival'); library('survival')
if (!require('survMisc')) install.packages('survMisc'); library('survMisc')
if (!require('survminer')) install.packages('survminer'); library('survminer')
#
#
## source data compile script
# source("~/Desktop/MetaAML_results/MetaAML_data_compile_script.R")
# load R data file from data compile script
load("~/Desktop/MetaAML_results/final_data_matrix.RData")

    
#### Mutation plot function ####
mutation_plot_function <- function(pt_subset, n_muts, vaf_or_variant, include, gene_list){
      
      # subset to desired cohorts
      if(pt_subset == "All"){
        final_data_matrix_sub <- final_data_matrix
        final_data_matrix_sub <-  final_data_matrix_sub[!duplicated(final_data_matrix_sub[1:2]),]
        label1 <- "All"
      }
      if(pt_subset == "De novo"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
        final_data_matrix_sub <-  final_data_matrix_sub[!duplicated(final_data_matrix_sub[1:2]),]
         label1 <- "De novo"
      }
      if(pt_subset == "Secondary"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "secondary")
        final_data_matrix_sub <-  final_data_matrix_sub[!duplicated(final_data_matrix_sub[1:2]),]
        label1 <- "Secondary"
      }
      if(pt_subset == "Relapse"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "relapse")
        final_data_matrix_sub <-  final_data_matrix_sub[!duplicated(final_data_matrix_sub[1:2]),]
        label1 <- "Relapse"
      }
      
      if(pt_subset == "Therapy related"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "therapy")
        final_data_matrix_sub <-  final_data_matrix_sub[!duplicated(final_data_matrix_sub[1:2]),]
        label1 <- "Therapy related"
      }
      if(pt_subset == "Other"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "other")
        final_data_matrix_sub <-  final_data_matrix_sub[!duplicated(final_data_matrix_sub[1:2]),]
        label1 <- "Other"
      }
      
      
      if(gene_list == "All Genes"){
        final_data_matrix_sub <- final_data_matrix_sub
      }
      
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      # subset plot to show individual gene list relationships
      if(include == "Default (all genes)"){
        final_data_matrix_sub <- final_data_matrix_sub
      }
      if(include == "Include (all rows)"){
        gene_list_sub <- as.data.frame(gene_list)
        colnames(gene_list_sub)[1] <- "genes_sub"
        gene_list_sub$genes_sub <- as.character(gene_list_sub$genes_sub)

        final_data_matrix_sub_1 <- setDT(final_data_matrix_sub)[Gene %chin% gene_list_sub$genes_sub] 
        final_data_matrix_sub <- setDT(final_data_matrix_sub)[Sample %chin% final_data_matrix_sub_1$Sample]
      }
      if(include == "Include (select rows)"){
        gene_list_sub <- as.data.frame(gene_list)
        colnames(gene_list_sub)[1] <- "genes_sub"
        gene_list_sub$genes_sub <- as.character(gene_list_sub$genes_sub)
        
        final_data_matrix_sub <- setDT(final_data_matrix_sub)[Gene %chin% gene_list_sub$genes_sub] 

      }
      if(include == "Exclude"){
        
        gene_list_sub <- as.data.frame(gene_list)
        colnames(gene_list_sub)[1] <- "genes_sub"
        gene_list_sub$genes_sub <- as.character(gene_list_sub$genes_sub)
        
        final_data_matrix_sub <- setDT(final_data_matrix_sub)[ ! Gene %chin% gene_list_sub$genes_sub, ] 
      }
      
      
      ### recalculate the frequency of mutations per gene and number of mutations per patient sample
      # find the frequency of mutations per patient
      combined_freq <- aggregate(data.frame(count = final_data_matrix_sub), list(value = final_data_matrix_sub$Sample), length)
      
      combined_pt_freq <- select(combined_freq, "value", "count.Sample")
      colnames(combined_pt_freq)[1] <- "Sample"
      
      combined_mid <- dplyr::left_join(final_data_matrix_sub, combined_pt_freq, by = "Sample")
      
      # find the frequency of mutations per gene
      combined_mut_freq <- aggregate(data.frame(count = final_data_matrix_sub), list(value = final_data_matrix_sub$Gene), length)
      
      combined_mut_freq <- select(combined_mut_freq, "value", "count.Gene")
      colnames(combined_mut_freq)[1] <- "Gene"
      
      combined_final <- dplyr::left_join(combined_mid, combined_mut_freq, by = "Gene")
      
      
      ## for aesthetic purposes, filter to only genes mutated at least "n" times
      n <- as.numeric(n_muts)
      
      combined_final_sub <- subset(combined_final, combined_final$count.Gene >= n)
      
      
      
      # make list of recurrent mutations
      mutations <- as.data.frame(unique(as.character(combined_final_sub$Gene)))
      colnames(mutations)[1] <- "Gene"
      mutations$Gene <- as.character(mutations$Gene)
      
      n_muts_1 <- as.numeric(nrow(mutations))
      n_muts_2 <- (n_muts_1 + 1)
      
      # make list of patient ids
      patients <- as.data.frame(unique(combined_final_sub$Sample))
      colnames(patients)[1] <- "Sample"
      
      n_pts <- as.numeric(nrow(patients))
      
      # make new data frame to populate with results
      temp_dat2 <- data.frame(matrix(NA, nrow = n_pts, ncol = n_muts_2))
      colnames(temp_dat2)[2:n_muts_2] <- c(mutations$Gene)
      
      colnames(temp_dat2)[1] <- "Sample"
      temp_dat2$Sample <- as.character(patients$Sample)
      
      n_pts <- as.numeric(nrow(temp_dat2))
      
      for(i in 1:nrow(combined_final_sub)){
        # i <- 2
        id <- as.character(combined_final_sub$Sample[i])
        mut <- as.character(combined_final_sub$Gene[i])
        
        r_num <- as.numeric(which(temp_dat2$Sample == id))
        c_num <- as.numeric(which(colnames(temp_dat2) == mut))
        
        temp_dat2[r_num,c_num] <- 1
      }
      
      temp_dat2[is.na(temp_dat2)] <- 0
      
      temp_dat2 <- as.data.frame(temp_dat2)
      
      combined_final_sub <- combined_final_sub[order(combined_final_sub$count.Gene, decreasing = T),]
      
      
      sort_desc <- c(as.character(combined_final_sub$Gene))
      
      
      # order the plot
      temp_dat2 <- temp_dat2 %>% 
        arrange_at(sort_desc, desc)
    
      pt_order <- as.data.frame(temp_dat2$Sample)
      colnames(pt_order)[1] <- "Sample"
      pt_order$order <- seq.int(nrow(pt_order))
      

      mut_table_final <- left_join(pt_order, combined_final_sub, by = "Sample")
      
      mut_table_final$variant_type <- as.character(mut_table_final$variant_type)

        
        # manually annotate 'insertion' and 'tandem duplications' variants as ITD in FLT3 cases
        for (i in 1:nrow(mut_table_final)) {
          if(!is.na(mut_table_final$variant_type[i])){
          if(mut_table_final$Gene[i] == "FLT3" & mut_table_final$variant_type[i] == "insertion"){
            mut_table_final$variant_type[i] <- "ITD"
          }
          if(mut_table_final$variant_type[i] == "tandem_duplication"){
            mut_table_final$variant_type[i] <- "ITD"
          }
          }
        }
      
      pt_order_by_muts <- select(mut_table_final, Sample, order, mut_freq_pt)
      pt_order_by_muts <- unique(pt_order_by_muts)
      pt_ord <- as.vector(pt_order_by_muts$Sample)
      
     p <-  ggbarplot(pt_order_by_muts, x = "Sample", y = "mut_freq_pt", order = pt_ord, ylab = "# muts") +
                # color = "black",
                # palette = c("#66c2a4", "#66c2a4","#66c2a4","#66c2a4", "#a6bddb","#a6bddb","#a6bddb"), alpha = .5) +
        theme(legend.position="none",axis.text.x=element_blank(), axis.title.x = element_blank()) +
        scale_y_continuous(expand = c(0,0), limits = c(0, max(mut_table_final$mut_freq_pt) + 5)) +
        theme(plot.title = element_text(hjust = 0.5))
      print(p)
      ggsave(filename = "~/Desktop/muts_per_patient.png", dpi = 300, width = 20, height = 1, units = "in")
      
      
        # plot data in CoMut plot
        if(vaf_or_variant == "Variants"){
          a <- ggplot(mut_table_final, aes( x = reorder(Sample, order), reorder(Gene, count.Gene), height = 0.8, width = 0.8)) +
            geom_tile(aes(fill = variant_type, color = variant_type)) +
            xlab(paste("n = ", n_pts, sep = "")) +
            ggtitle(paste("Mutation landscape in", label1, "AML", sep = " ")) +
            scale_fill_manual(values = c("Deletion" = "#374E55FF", "INDEL" = "#DF8F44FF", "Insertion" = "#00A1D5FF", "ITD" = "#79AF97FF", "SNV" = "#B24745FF", "Splicing" = "#6A6599FF", "Unknown" = "#80796BFF")) +      
            scale_color_manual(values = c("Deletion" = "#374E55FF", "INDEL" = "#DF8F44FF", "Insertion" = "#00A1D5FF", "ITD" = "#79AF97FF", "SNV" = "#B24745FF", "Splicing" = "#6A6599FF", "Unknown" = "#80796BFF")) + 
            theme(
              legend.key = element_rect(fill='NA'),
              legend.key.size = unit(.75, 'cm'),
              legend.title = element_blank(),
              legend.position = "right",
              legend.text = element_text(size=15, face="bold"),
              axis.ticks.x = element_blank(),
              axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(colour = "#252525", size = 15),
              axis.title.x = element_text(face = "bold"),
              axis.title.y = element_blank(),
              panel.grid.major.x = element_blank(),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.background = element_blank())
          print(a)
          ggsave(filename = "~/Desktop/oncoprint_variants.png", dpi = 300, width = 7.5, height = 3, units = "in")
          
        }

        if(vaf_or_variant == "VAFs"){
          a <- ggplot(mut_table_final, aes( x = reorder(Sample, order), reorder(Gene, count.Gene), height = 0.8, width = 0.8)) +
            geom_tile(aes(fill = VAF)) + 
            scale_fill_gradient2(low="#08519c", mid = "grey", high="#a50f15", midpoint = 0.5, na.value = "#252525") +
            scale_color_gradient2(low="#08519c", mid = "grey", high="#a50f15", midpoint = 0.5, na.value = "#252525") +
            xlab(paste("n = ", n_pts, sep = "")) +
            ggtitle(paste("Mutation co-occurence and VAF in", label1, "AML", sep = " ")) +
            theme( 
              legend.key = element_rect(fill='NA'),
              legend.key.size = unit(1, 'cm'),
              legend.title = element_blank(),
              legend.position = "right",
              legend.text = element_text(size=15, face="bold"),
              axis.ticks.x = element_blank(), 
              axis.ticks.y = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_text(colour = "#252525"), 
              axis.title.x = element_text(face = "bold"), 
              axis.title.y = element_blank(),
              panel.grid.major.x = element_blank(), 
              panel.grid.major.y = element_blank(),
              panel.grid.minor.x = element_blank(), 
              panel.grid.minor.y = element_blank(), 
              panel.background = element_blank())
          print(a)
          ggsave(filename = "~/Desktop/oncoprint_vafs.png", dpi = 300, width = 20, height = 10, units = "in")
        }
        
      }

# mutation_plot_function(pt_subset = "All", n_muts = 25, vaf_or_variant = "variant", include = "Include (all rows)", gene_list = "All Genes")
    
#### run co-mutation analysis ####
co_mutation_analysis_function <- function(pt_subset, include, gene_list){
      
      # subset to desired cohorts
      if(pt_subset == "All"){
        final_data_matrix_sub <- final_data_matrix
        label1 <- "All"
      }
      if(pt_subset == "De novo"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
        label1 <- "De novo"
      }
      if(pt_subset == "Secondary"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "secondary")
        label1 <- "Secondary"
      }
      if(pt_subset == "Relapse"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "relapse")
        label1 <- "Relapse"
      }
      if(pt_subset == "Therapy related"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "therapy")
        label1 <- "Therapy related"
      }
      if(pt_subset == "Other"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "other")
        label1 <- "Other"
      }
      
      if(gene_list == "All Genes"){
        final_data_matrix_sub <- final_data_matrix_sub
      }
      
      
      ### recalculate the frequency of mutations per gene and number of mutations per patient sample
      # find the frequency of mutations per patient
      combined_freq <- aggregate(data.frame(count = final_data_matrix_sub), list(value = final_data_matrix_sub$Sample), length)
      
      combined_pt_freq <- select(combined_freq, "value", "count.Sample")
      colnames(combined_pt_freq)[1] <- "Sample"
      
      combined_mid <- dplyr::left_join(final_data_matrix_sub, combined_pt_freq, by = "Sample")
      
      # find the frequency of mutations per gene
      combined_mut_freq <- aggregate(data.frame(count = final_data_matrix_sub), list(value = final_data_matrix_sub$Gene), length)
      
      combined_mut_freq <- select(combined_mut_freq, "value", "count.Gene")
      colnames(combined_mut_freq)[1] <- "Gene"
      
      combined_final <- dplyr::left_join(combined_mid, combined_mut_freq, by = "Gene")
      
      
      ## for aesthetic purposes, filter to only genes mutated at least "n" times
      n <- 20
      
      combined_final_sub <- subset(combined_final, combined_final$count.Gene >= n)

      # make list of recurrent mutations
      mutations <- as.data.frame(unique(as.character(combined_final_sub$Gene)))
      colnames(mutations)[1] <- "Gene"
      mutations$Gene <- as.character(mutations$Gene)
      
      n_muts_1 <- as.numeric(nrow(mutations))
      n_muts_2 <- (n_muts_1 + 1)
      
      # make list of patient ids
      patients <- as.data.frame(unique(combined_final_sub$Sample))
      colnames(patients)[1] <- "Sample"
      
      n_pts <- as.numeric(nrow(patients))
      
      # make new data frame to populate with results
      temp_dat2 <- data.frame(matrix(NA, nrow = n_pts, ncol = n_muts_2))
      colnames(temp_dat2)[2:n_muts_2] <- c(mutations$Gene)
      
      colnames(temp_dat2)[1] <- "Sample"
      temp_dat2$Sample <- as.character(patients$Sample)
      
      n_pts <- as.numeric(nrow(temp_dat2))
      
      for(i in 1:nrow(combined_final_sub)){
        # i <- 2
        id <- as.character(combined_final_sub$Sample[i])
        mut <- as.character(combined_final_sub$Gene[i])
        
        r_num <- as.numeric(which(temp_dat2$Sample == id))
        c_num <- as.numeric(which(colnames(temp_dat2) == mut))
        
        temp_dat2[r_num,c_num] <- 1
      }
      
      temp_dat2[is.na(temp_dat2)] <- 0
      
      temp_dat2 <- as.data.frame(temp_dat2)
      
      combined_final_sub <- combined_final_sub[order(combined_final_sub$count.Gene, decreasing = T),]
      
      
      sort_desc <- c(as.character(combined_final_sub$Gene))
      
      
      # order the plot
      temp_dat2 <- temp_dat2 %>% 
        arrange_at(sort_desc, desc)
      
      
      pt_order <- as.data.frame(temp_dat2$Sample)
      colnames(pt_order)[1] <- "Sample"
      pt_order$order <- seq.int(nrow(pt_order))
      
      mut_table_final <- left_join(pt_order, combined_final_sub, by = "Sample")
      
      mut_table_final$variant_type <- as.character(mut_table_final$variant_type)
      
      # manually annotate 'insertion' and 'tandem duplications' variants as ITD in FLT3 cases
      for (i in 1:nrow(mut_table_final)) {
        if(mut_table_final$Gene[i] == "FLT3" & mut_table_final$variant_type[i] == "insertion"){
          mut_table_final$variant_type[i] <- "ITD"
        }
        if(mut_table_final$variant_type[i] == "tandem_duplication"){
          mut_table_final$variant_type[i] <- "ITD"
        }
      }
      
      
      #### perform Fisher's exact test on mutation occurences in the combined data set ####
      # deliniate between FLT3-ITD and FLT3-TKD because they will have different patterns with mutations (...I think)
      
      # create fresh dataframe
      final_data_matrix <- mut_table_final
      
      # manually annotate FLT as TKD or ITD
      for (i in 1:nrow(final_data_matrix)) {
        if(final_data_matrix$Gene[i] == "FLT3" & final_data_matrix$variant_type[i] == "ITD"){
          final_data_matrix$Gene[i] <- "FLT3_ITD"
        }
        if(final_data_matrix$Gene[i] == "FLT3" & final_data_matrix$variant_type[i] == "deletion"){
          final_data_matrix$Gene[i] <- "FLT3_ITD"
        }
        if(final_data_matrix$Gene[i] == "FLT3" & final_data_matrix$variant_type[i] == "SNV"){
          final_data_matrix$Gene[i] <- "FLT3"
        }
      }
      
      # make list of recurrent mutations
      mutations_2 <- as.data.frame(unique(as.character(final_data_matrix$Gene)))
      colnames(mutations_2)[1] <- "Gene"
      mutations_2$Gene <- as.character(mutations_2$Gene)
      
      n_muts_2 <- as.numeric(nrow(mutations_2))
      
      n_muts_3 <- (n_muts_2 + 1)
      
      # make new data frame to populate
      temp_dat3 <- data.frame(matrix(NA, nrow = n_pts, ncol = n_muts_3))
      colnames(temp_dat3)[2:n_muts_3] <- c(mutations_2$Gene)
      
      temp_dat3[,order(colnames(temp_dat3))]
      
      colnames(temp_dat3)[1] <- "Sample"
      temp_dat3$Sample <- as.character(patients$Sample)
      
      temp_dat3 <- temp_dat3[,order(colnames(temp_dat3))]
      
      n_pts <- as.numeric(nrow(temp_dat3))
      
      for(i in 1:nrow(final_data_matrix)){
        # i <- 2
        id <- as.character(final_data_matrix$Sample[i])
        mut <- as.character(final_data_matrix$Gene[i])
        
        r_num <- as.numeric(which(temp_dat3$Sample == id))
        c_num <- as.numeric(which(colnames(temp_dat3) == mut))
        
        temp_dat3[r_num,c_num] <- 1
      }
      
      temp_dat3[is.na(temp_dat3)] <- 0
      
      temp_dat3 <- as.data.frame(temp_dat3)
      
      
      
      # create new dataframe for test results
      f_test_matrix <- data.frame(matrix(NA, nrow = n_muts_2, ncol = n_muts_2))
      temp_dat3$Sample <- NULL
      c_names <- colnames(temp_dat3)
      
      colnames(f_test_matrix) <- c(c_names)
      rownames(f_test_matrix) <- c(c_names)
      
      # create dataframe for p-values from fisher's exact test
      f_test_matrix_p <- f_test_matrix
      
      for(i in 1:nrow(f_test_matrix)){
# i <- 1
        r_name <- as.character(rownames(f_test_matrix[i,]))
        for(j in 1:ncol(f_test_matrix)){
          if(i < j | i > j){
            # j <- 10
            c_name <- as.character(colnames(f_test_matrix)[j])
            
            # select the columns in the binarized mutation file
            c2_num <- as.numeric(grep(r_name, colnames(temp_dat3)))
            
            if(r_name == "FLT3" | c_name == "FLT3"){
              c2_num <- c2_num[1]
            }
            
            c1_num <- as.numeric(grep(c_name, names(temp_dat3)))
            
            if(r_name == "FLT3" | c_name == "FLT3"){
              c1_num <- c1_num[1]
            }
            
            f_test <- as.data.frame(temp_dat3[,c(c1_num,c2_num)])
            f_test$sum <- (f_test[,1] + f_test[,2])
            
            f_test_sub <- data.frame(matrix(NA, nrow = 2, ncol = 2))
            
            colnames(f_test_sub)[1:2] <- c(r_name, c_name)
            rownames(f_test_sub)[1:2] <- c(c_name, r_name)
            
            f_test_sub[1,1] <- as.numeric(nrow(subset(f_test, f_test[,1] == 0 & f_test[,2] == 0)))
            f_test_sub[1,2] <- as.numeric(nrow(subset(f_test, f_test[,1] == 1)))
            f_test_sub[2,1] <- as.numeric(nrow(subset(f_test, f_test[,2] == 1)))
            f_test_sub[2,2] <- as.numeric(nrow(subset(f_test, f_test[,1] == 1 & f_test[,2] == 1)))
  
            # odds ratio
            odds <- (f_test_sub[1,1]*f_test_sub[2,2])/(f_test_sub[1,2]*f_test_sub[2,1])
            
            p_val <- fisher.test(f_test_sub)$p.value
            
            f_test_matrix[i,j] <- odds
            f_test_matrix_p[i,j] <- p_val 
          }
        }
      }
      

    
      # Get lower triangle of the correlation matrix
      get_lower_tri<-function(f_test_matrix){
        f_test_matrix[lower.tri(f_test_matrix)] <- NA
        return(f_test_matrix)
      }

      final_cor_frame_d <- get_lower_tri(f_test_matrix)

      final_cor_frame_d$genes <- rownames(final_cor_frame_d)

      # Melt the correlation matrix
      melted_cormat_d <- melt(final_cor_frame_d, na.rm = TRUE)
      colnames(melted_cormat_d)[3] <- "Delta"
      
      
      # Get lower triangle of the correlation matrix
      get_lower_tri<-function(f_test_matrix_p){
        f_test_matrix_p[lower.tri(f_test_matrix_p)] <- NA
        return(f_test_matrix_p)
      }
      
      final_cor_frame_p <- get_lower_tri(f_test_matrix_p)
      
      final_cor_frame_p$genes <- rownames(final_cor_frame_p)
      
      # Melt the correlation matrix
      melted_cormat_p <- melt(final_cor_frame_p, na.rm = TRUE)
      colnames(melted_cormat_p)[3] <- "P_value"
      
      # create the final dataframe with all variables used for plotting
      melted_cormat <- left_join(melted_cormat_p, melted_cormat_d, by = c("genes", "variable"))
      colnames(melted_cormat)[1:4] <- c("Gene_1", "Gene_2", "P_value", "Delta")

      melted_cormat$P_value <- as.numeric(as.character(melted_cormat$P_value))
      
      melted_cormat$fdr_q_value <- p.adjust(melted_cormat$P_value, method = "fdr")
      
      # remove all p values that are not significant
      for(i in 1:nrow(melted_cormat)){
        if(melted_cormat$fdr_q_value[i] > 0.1){
          melted_cormat$fdr_q_value[i] <- NA
        }
      }

      
      # subset to genes of interest
      if(include == "Include (select rows)"){
        gene_list_sub <- as.data.frame(gene_list)
        colnames(gene_list_sub)[1] <- "genes_sub"
        gene_list_sub$genes_sub <- as.character(gene_list_sub$genes_sub)
        
        melted_cormat <- setDT(melted_cormat)[Gene_1 %chin% gene_list_sub$genes_sub]
        # melted_cormat <- setDT(melted_cormat)[Gene_2 %chin% gene_list_sub$genes_sub]
        
      }

      
      # plot the data
      # define color bins
      melted_cormat$odds_bin <- cut(melted_cormat$Delta, breaks = c(-0.001, 0, 0.01, 0.1, 1, 10, 30))
      for(i in 1:nrow(melted_cormat)){
        if(melted_cormat$P_value[i] > 0.05){
          melted_cormat$odds_bin[i] <- NA
        }
      }
    
      a <-  ggplot(data = melted_cormat, aes(x=Gene_1, y = Gene_2), fill = odds_bin) +
        geom_tile(aes(fill = odds_bin) , color = "#252525")+
        scale_fill_manual(values=c("darkred", '#d73027', '#fcbba1', "#91bfdb", "#4575b4"),
                          na.value = "white",
                          # breaks = c("0.001", "0.01", "0.1", "1", "10", "100", "1000"),
                          labels = c("0", "0.01-0.1", "0.1-0.09", "1.1-10", "ns"),
                          name= "Odds ratio \n(p<0.05)", guide = guide_legend(reverse = F)) +
        geom_point(aes(size = fdr_q_value)) +
        scale_size(range = c(2,0), name = "q value") +
        # guides(
        #   color= "black", 
        #   size=guide_legend(override.aes = list(size = fdr_q_value))
        # ) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5,
                                         size = 8, hjust = 0),
              # panel.grid.major = element_line(size = 0.5, linetype = 'solid',
              #                                 colour = "#f0f0f0"),
              axis.text.y = element_text(size = 8),
              legend.position = c(0.65,0.35))+
        theme(axis.text.x.top = element_text(vjust = 0.5)) +
        scale_x_discrete(position = "top") +
        xlab(label= "") +
        ylab(label="") +
        labs(title = "Correlation of Mutations (Fisher's Exact)")
      

 print(a)
 ggsave(filename = "~/Desktop/fishers_test.png", dpi = 300, width = 7, height = 7, units = "in")
 }

    
##### Survival analysis and function ######
survival_function <- function(pt_subset, background, include, survival_mutations){
      
      # subset to desired cohorts
      if(pt_subset == "All"){
        final_data_matrix_sub <- final_data_matrix
        label1 <- "All"
      }
      if(pt_subset == "De novo"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
        label1 <- "De novo"
      }
      if(pt_subset == "Secondary"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "secondary")
        label1 <- "Secondary"
      }
      if(pt_subset == "Relapse"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "relapse")
        label1 <- "Relapse"
      }
      if(pt_subset == "Therapy related"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "therapy")
        label1 <- "Therapy related"
      }
      if(pt_subset == "Other"){
        final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "other")
        label1 <- "Other"
      }
      
      if(background != "None"){
        if(background != "DNMT3A" & background != "FLT3" & background != "NPM1"){
          final_data_matrix_sub <- subset(final_data_matrix_sub, subset=(Sex == background | Cytogenetics == background | Risk == background))
        }
        if(background == "DNMT3A" | background == "FLT3" | background == "NPM1"){
          pts <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == background)
          final_data_matrix_sub <- setDT(final_data_matrix_sub)[(Sample) %chin% pts$Sample] 
        }
      }

      
      gene1 <- survival_mutations[1]
      gene2 <- survival_mutations[2]
      
      # find patients with mutations in the genes
      final_data_matrix_sub1 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == gene1)
      final_data_matrix_sub2 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == gene2)
      final_data_matrix_sub3 <- setDT(final_data_matrix_sub1)[(Sample) %chin% final_data_matrix_sub2$Sample]
      
      # create strata columns
      final_data_matrix_sub1$STRATA1 <- 1
      final_data_matrix_sub2$STRATA2 <- 2
      final_data_matrix_sub3$STRATA3 <- 3
      
      # subset to informative columns
      final_data_matrix_sub1 <- select(final_data_matrix_sub1, Sample, STRATA1)
      final_data_matrix_sub2 <- select(final_data_matrix_sub2, Sample, STRATA2)
      final_data_matrix_sub3 <- select(final_data_matrix_sub3, Sample, STRATA3)
      
      # add strata to the larger data matrix
      final_data_matrix_sub <- left_join(final_data_matrix_sub, final_data_matrix_sub1, by = "Sample")
      final_data_matrix_sub <- left_join(final_data_matrix_sub, final_data_matrix_sub2, by = "Sample")
      final_data_matrix_sub <- left_join(final_data_matrix_sub, final_data_matrix_sub3, by = "Sample")
      
      # create consensus strate column
      final_data_matrix_sub$STRATA <- rowSums(final_data_matrix_sub[,c("STRATA1", "STRATA2", "STRATA3")], na.rm=TRUE)
      
      final_data_matrix_sub$STRATA <- as.numeric(final_data_matrix_sub$STRATA)
      
      final_data_matrix_sub[,14:16] <- NULL
      
      # filter to unique patients
      final_data_matrix_sub_final <- final_data_matrix_sub[!duplicated(final_data_matrix_sub[1]),]
      # final_data_matrix_sub_final <- na.omit(final_data_matrix_sub_final)
      
      nstrata <- as.numeric(length(unique(final_data_matrix_sub_final$STRATA)))
      
      # 
      if(nstrata != 3){
      c <- c("#374E55FF", "#DF8F44FF", "#00A1D5FF", "#B24745FF")
      l_labs <- c("neither", gene1, gene2, "both")
      }

      if(nstrata == 3){
        c <- c("#374E55FF", "#DF8F44FF", "#00A1D5FF", "#B24745FF")
        l_labs <- c("neither", gene1, gene2)
      }
      
      # change the survival time to years format
      final_data_matrix_sub_final$Time_to_OS <- (final_data_matrix_sub_final$Time_to_OS/365)
      
      # create the survival data 
      final_data_matrix_sub_final$OS <- with(final_data_matrix_sub_final, Surv(Time_to_OS, Censor == 1))
      
      # create the survival objects used to plot kaplan-meyer curves
      OS <- survfit(OS ~ STRATA, data = final_data_matrix_sub_final, conf.type = "log-log")
  
      
      # find the different p-values for the different comparisons
      final_data_matrix_sub_final$Censor <- as.numeric(final_data_matrix_sub_final$Censor)  
      res <- pairwise_survdiff(Surv(Time_to_OS, Censor) ~ STRATA,
                               data = final_data_matrix_sub_final)
      print(res)

        
        # plots the survival
        surv_plot <- ggsurvplot(OS,
                                data = final_data_matrix_sub_final,
                                log = (OS),
                                log.rank.weights = c("survdiff"),
                                pval = T,
                                test.for.trend = T,
                                pval.method.size = 3,
                                pval.coord = c(0, 0),
                                conf.int = F,
                                censor = T,
                                surv.median.line = "none",
                                risk.table = F,
                                risk.table.title = "",
                                risk.table.fontsize = 4,
                                risk.table.height = .25,
                                risk.table.y.text = T,
                                break.time.by = 5,
                                risk.table.pos = c("out"),
                                palette = c,
                                title = paste("Survival by ", gene1, " and ", gene2, " mutations ", "(",label1,")", sep = ""),
                                 xlab = "Years",
                                ylim = c(0, 1.0),
                                ylab =  "Survival Probability",
                                font.main = c(20, "plain", "#252525"),
                                pval.size = 5,
                                font.x = c(20, "plain", "#252525"),
                                font.y =  c(20, "plain", "#252525"),
                                font.legend = c(15, "plain"),
                                font.tickslab = c(15, "plain", "#252525"),
                                legend.labs = l_labs,
                                legend.title = "Mutation status",
                                legend = "right",
                                ggtheme = theme(plot.title = element_text(hjust = 0.5)))
        
        print(surv_plot)

        png(filename = "~/Desktop/meta_aml_survival.png", res = 300, width = 8, height = 7.5, units = "in")

        surv_plot
        print(surv_plot)
        dev.off()
        
    }

 
### Scatterplot function to plot the VAF for two user-defined mutation inputs ####

vaf_scatterplot_function <- function(pt_subset, gene_1_2){
  
  # subset to desired cohorts
  if(pt_subset == "All"){
    final_data_matrix_sub <- final_data_matrix
    label1 <- "All"
  }
  if(pt_subset == "De novo"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
    label1 <- "De novo"
  }
  if(pt_subset == "Secondary"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "transformed")
    label1 <- "Secondary"
  }
  if(pt_subset == "Relapse"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "relapse")
    label1 <- "Relapse"
  }
  if(pt_subset == "Therapy related"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "therapy")
    label1 <- "Therapy related"
  }
  if(pt_subset == "Other"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "other")
    label1 <- "Other"
  }
  
  # select mutations of interest
  gene_x <- gene_1_2[1]
  gene_y <- gene_1_2[2]
  
  sub_gene_1 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == gene_x)
  sub_gene_2 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == gene_y)
  
  # select patients with both mutations
  gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
  
  # add point color columns for visualizing clonal/subclonal trends
  gene_1_and_2$vaf_ratio <- (gene_1_and_2$VAF.x - gene_1_and_2$VAF.y)
  
  gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$vaf_ratio), ] 
  
  # define point color
  gene_1_and_2$Clonality <- NA
  
  for(i in 1:nrow(gene_1_and_2)){
    if(gene_1_and_2$vaf_ratio[i] <= .05 & gene_1_and_2$vaf_ratio[i] >= -.05){
      gene_1_and_2$Clonality[i] <- "Ambiguous"
    }
    if(gene_1_and_2$vaf_ratio[i] > .05){
      gene_1_and_2$Clonality[i] <- paste(gene_x, "first", sep = " ")
    }
    if(gene_1_and_2$vaf_ratio[i] < -.05){
      gene_1_and_2$Clonality[i] <- paste(gene_y, "first", sep = " ")
    }
  }
  
  gene_1_and_2 <- gene_1_and_2[order(gene_1_and_2$Sample, -gene_1_and_2$VAF.x),]
  gene_1_and_2= gene_1_and_2[!duplicated(gene_1_and_2$Sample),]
  
  if(nrow(gene_1_and_2) > 4){
    # make the scatterplot
    
    scatter_plot = ggplot(gene_1_and_2, aes(x = gene_1_and_2$VAF.y, y = gene_1_and_2$VAF.x)) +
      # geom_rect(mapping=aes(xmin=-Inf, xmax=.5, ymin=-Inf, ymax=.5), fill = "lightgrey", color=NA, alpha=.5) +
      geom_point(aes(color = Clonality), size = 4, alpha = 0.75) +
      xlim(0,1)+
      ylim(0,1)+
      xlab(paste(gene_y, " VAF", sep ="")) +
      ylab(paste(gene_x, " VAF", sep ="")) +  
      scale_color_manual(values = c("#374E55FF","#8c510a","#01665e")) +
      geom_abline(intercept = .05, slope = (1), color="#969696",
                  linetype="dashed", size=.5)+
      geom_abline(intercept = -.05, slope = (1), color="#969696",
                  linetype="dashed", size=.5)+
      geom_point(shape = 1, size =  4,colour = "black") +
      theme(plot.title = element_text(hjust = 0.5, paste(gene_y, "vs.", gene_x, sep = "")), legend.position = "right", legend.title = element_blank(), legend.text = element_text(size=15, face="bold")) 
    
    print(scatter_plot)

      ggsave(filename =paste("~/Desktop/MetaAML_",gene_x, "_",gene_y,"_scatterplot.png", sep = ""), dpi = 300, width = 6.5, height = 4.5, units = "in")
  }
}


### Scatterplot survival function to analyse differences in survival between clonal and subclonal patterns in gene relationships ####
vaf_scatterplot_survival_function <- function(pt_subset, gene_1_2){
  
  # subset to desired cohorts
  if(pt_subset == "All"){
    final_data_matrix_sub <- final_data_matrix
    label1 <- "All"
  }
  if(pt_subset == "De novo" | pt_subset == "De_novo" | pt_subset == "de_novo"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
    label1 <- "De novo"
  }
  if(pt_subset == "Secondary"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "transformed")
    label1 <- "Secondary"
  }
  if(pt_subset == "Relapse"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "relapse")
    label1 <- "Relapse"
  }
  if(pt_subset == "Therapy related"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "therapy")
    label1 <- "Therapy related"
  }
  if(pt_subset == "Other"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "other")
    label1 <- "Other"
  }
  
  # select mutations of interest
  g1 <- gene_1_2[1]
  g2 <- gene_1_2[2]
  
  sub_gene_1 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == g1)
  sub_gene_2 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == g2)
  
  # select patients with both mutations
  gene_1_and_2 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
  
  # add point color columns for visualizing clonal/subclonal trends
  gene_1_and_2$vaf_ratio <- (gene_1_and_2$VAF.x - gene_1_and_2$VAF.y)
  
  gene_1_and_2 <- gene_1_and_2[complete.cases(gene_1_and_2$vaf_ratio), ] 
  
  # define point color
  gene_1_and_2$Clonality <- NA
  
  # assign order
  for(i in 1:nrow(gene_1_and_2)){
    if(gene_1_and_2$vaf_ratio[i] <= .05 & gene_1_and_2$vaf_ratio[i] >= -.05){
      gene_1_and_2$Clonality[i] <- "Ambiguous"
    }
    if(gene_1_and_2$vaf_ratio[i] > .05){
      gene_1_and_2$Clonality[i] <- paste(g1, "first", sep = " ")
    }
    if(gene_1_and_2$vaf_ratio[i] < -.05){
      gene_1_and_2$Clonality[i] <- paste(g2, "first", sep = " ")
    }
  }
  
  # if(ambiguous_curve == F){
    gene_1_and_2 = subset(gene_1_and_2, gene_1_and_2$Clonality != "Ambiguous")
    p_val = T
    pal = c("#8c510a","#01665e")
    lab = c(paste(g2, "first",sep = " "), paste(g1, "first",sep = " "))
  # }
  
  # make sure to only include one row per patient if there are two mutations in the same gene
  gene_1_and_2 <- gene_1_and_2[order(gene_1_and_2$Sample, -gene_1_and_2$VAF.x),]
  gene_1_and_2= gene_1_and_2[!duplicated(gene_1_and_2$Sample),]
  
  # create the survival data object
  gene_1_and_2=as.data.frame(gene_1_and_2) %>% distinct(Sample, Censor.x, Time_to_OS.x, Clonality, .keep_all = F)
  gene_1_and_2=na.omit(gene_1_and_2)
  gene_1_and_2$Censor.x=as.numeric(gene_1_and_2$Censor.x)
  gene_1_and_2$Time_to_OS.x <- (gene_1_and_2$Time_to_OS.x/365)
  gene_1_and_2$Time_to_OS.x=as.numeric(gene_1_and_2$Time_to_OS.x)
  
  # create survival object
  gene_1_and_2$OS <- with(gene_1_and_2, Surv(Time_to_OS.x, Censor.x == 1))
  OS <- survfit(OS ~ Clonality, data = gene_1_and_2, conf.type = "log-log")
  
  # find all the p-values for the different comparisons
  res <- pairwise_survdiff(Surv(Time_to_OS.x, Censor.x) ~ Clonality,
                           data = gene_1_and_2)
  print(res)
  
  # plots the survival
  surv_plot <- ggsurvplot(OS,
                          data = gene_1_and_2,
                          log = (OS),
                          log.rank.weights = c("survdiff"),
                          pval = p_val,
                          test.for.trend = F,
                          pval.method.size = 3,
                          pval.coord = c(0, 0),
                          conf.int = F,
                          censor = T,
                          surv.median.line = "none",
                          risk.table = T,
                          risk.table.title = "",
                          risk.table.fontsize = 4,
                          risk.table.height = .3,
                          risk.table.y.text = T,
                          break.time.by = 5,
                          risk.table.pos = c("out"),
                          palette = pal,
                          xlab = "Years",
                          ylim = c(0, 1.0),
                          ylab =  "Survival Probability",
                          font.main = c(15, "plain", "#252525"),
                          pval.size = 4,
                          font.x = c(12, "plain", "#252525"),
                          font.y =  c(12, "plain", "#252525"),
                          font.legend = c(12, "plain"),
                          font.tickslab = c(12, "plain", "#252525"),
                          legend.labs = lab,
                          legend.title = 'Mutation order',
                          legend = "right",
                          ggtheme = theme(plot.title = element_text(hjust = 0.5)))
  
  print(surv_plot)
  png(filename = paste("~/Desktop/",g1,"_",g2,".png", sep = ""), res = 300, width = 5, height = 5, units = "in")
  surv_plot
  print(surv_plot)
  dev.off()
}
    
co_occuring_clonality_function <- function(pt_subset, gene_1_2_3, order_gene){
  
  # subset to desired cohorts
  if(pt_subset == "All"){
    final_data_matrix_sub <- final_data_matrix
    label1 <- "All"
  }
  if(pt_subset == "De novo"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "de_novo")
    label1 <- "De novo"
  }
  if(pt_subset == "Secondary"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "secondary")
    label1 <- "Secondary"
  }
  if(pt_subset == "Relapse"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "relapse")
    label1 <- "Relapse"
  }
  if(pt_subset == "Therapy related"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "therapy")
    label1 <- "Therapy related"
  }
  if(pt_subset == "Other"){
    final_data_matrix_sub <- subset(final_data_matrix, final_data_matrix$Subset == "other")
    label1 <- "Other"
  }
  
  # select mutations of interest
  gene_x <- gene_1_2_3[1]
  gene_y <- gene_1_2_3[2]
  gene_z <- gene_1_2_3[3]
  
  sub_gene_1 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == gene_x)
  sub_gene_2 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == gene_y)
  sub_gene_3 <- subset(final_data_matrix_sub, final_data_matrix_sub$Gene == gene_z)
  
  # select patients with both mutations
  gene_12 <- inner_join(sub_gene_1, sub_gene_2, by = "Sample")
  gene_123 <- inner_join(gene_12, sub_gene_3, by = "Sample")
  
  gene_123 <- na.omit(select(gene_123, Sample, Gene, VAF, VAF.x, VAF.y))
  gene_123 <- gene_123[!duplicated(gene_123[,c("Sample")]),]
  
  n_sample <- as.numeric(nrow(gene_123))
  
  gene_123_sub <- select(gene_123, Sample, VAF.x, VAF.y, VAF)
  colnames(gene_123_sub) <- c("Sample", paste(gene_x), paste(gene_y), paste(gene_z))
  
  
  gene_123_sub_m <- melt(gene_123_sub)
  
  gene <- order_gene
  
  gene_123_sub <- subset(gene_123_sub_m, gene_123_sub_m$variable == gene)
  gene_123_sub <- gene_123_sub[order(gene_123_sub$value, decreasing = T),]
  order <- as.data.frame(gene_123_sub$Sample)
  order$order_sample <- 1:nrow(order)
  colnames(order)[1] <- "Sample"
  
  gene_123_sub_m <- left_join(gene_123_sub_m, order, by = "Sample")
  
  three_gene_vaf_plot <-  ggplot(gene_123_sub_m, aes(x = reorder(Sample, order_sample), value)) +
    geom_point(aes(color = as.factor(variable)), size = 3.5, alpha = .75) +
    scale_color_manual(values = c("#374E55FF", "#DF8F44FF", "#00A1D5FF")) +
    ggtitle(paste("Clonal Patterns of Co-occuring\n",gene_x, ", ", gene_y, ", & ", gene_z, " Mutations", sep = "")) +
    xlab(paste("Samples (n = ", n_sample, ")", sep = "")) +
    ylab("VAF") +
    theme(axis.text.x=element_blank(),
          legend.key = element_rect(fill='NA'),
          legend.key.size = unit(.5, 'cm'),
          legend.title = element_blank(),
          legend.position = "right",
          legend.text = element_text(size=15, face="bold"))
  
  print(three_gene_vaf_plot)
}



#### Main server function ####
shinyServer(function(input, output, session){
  
  # Listens for the "analysis_start" variable which is the clickable button
  # When the button is clicked, these functions run
  run_comut <- eventReactive(input$analysis_start, {
    req(input$cohort_choices)
    req(input$inclusion_cutoff)
    req(input$subset_to_variant_of_vaf)
    req(input$gene_subsets)
    req(input$subset_to_genes)
    
    # functions
    mutation_plot_function(pt_subset = input$cohort_choices,  n_muts = input$inclusion_cutoff, vaf_or_variant = input$subset_to_variant_of_vaf, include = input$subset_to_genes, gene_list = input$gene_subsets)


    
  })
  
  run_cooccurrance <- eventReactive(input$analysis_start, {
    req(input$cohort_choices)
    req(input$inclusion_cutoff)
    req(input$subset_to_genes)
    req(input$gene_subsets)

    # functions
    co_mutation_analysis_function(pt_subset = input$cohort_choices,include = input$subset_to_genes, gene_list = input$gene_subsets)
    
  })
  
  run_scatter <- eventReactive(input$analysis_start, {
    req(input$cohort_choices)
    req(input$gene_scatter)
    
    # function
    vaf_scatterplot_function(pt_subset = input$cohort_choices, gene_1_2 = input$gene_scatter)
    
  })
  
  run_survscatter <- eventReactive(input$analysis_start, {
    # req(input$cohort_choices)
    # req(input$gene_scatter)
    
    # function
    vaf_scatterplot_survival_function(pt_subset = input$cohort_choices, gene_1_2 = input$gene_scatter)    
    
  })
  
  run_vaf_ordering <- eventReactive(input$analysis_start, {
    req(input$cohort_choices)
    req(input$genes_vaf)
    
    # function
    co_occuring_clonality_function(pt_subset = input$cohort_choices, gene_1_2_3 = input$genes_vaf, order_gene = input$genes_vaf[1])    
    
  })
  
  run_survival <- eventReactive(input$analysis_start, {
    req(input$cohort_choices)
    req(input$subset_to_genes)
    
    # print(input$subset_to_genes)
    
    # function
    survival_function(pt_subset = input$cohort_choices, background = input$genetic_background, include = input$subset_to_genes, survival_mutations = input$gene_survival)

    
  })
  
  # Renders the plots that are output by the plotting functions
  output$comut_plot <- renderPlot({
    run_comut()
  })
  
  output$cooccurrance_plot <- renderPlot({
    run_cooccurrance()
  })
  
  output$scatter_plot <- renderPlot({
    run_scatter()
  })
  
  output$surv_scatter_plot <- renderPlot({
    run_survscatter()
  })
  
  output$three_gene_vaf_plot <- renderPlot({
    run_vaf_ordering()
  })
  
  output$survival_plot <- renderPlot({
    run_survival()
  })
})
