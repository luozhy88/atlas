
library(dplyr)
library(tibble)
library(DescTools)
library(purrr)
library(ggplot2)
library(patchwork)


pathway_pdf=function(meta_info=meta_info,pathway_df=pathway_df,meta_col=meta_col,out.name=out.name){
  
  ####args
  # meta_info<-meta.348 # the column of meta_info  is c("sampleID", "Group", "Gender")
  # meta_col<-c("sampleID", "Group", "Gender")
  # pathway_df<-kegg.pathway2 # the column of pathway_df is SampleID; The rowname of pathway_df is pathway name.
  # out.name="KEGG_"
  ###
  
  DF <- pathway_df %>%column_to_rownames(var = "name_id") %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column(var = "sampleID") %>%
    inner_join(meta_info) %>%
    dplyr::select(meta_col, everything())
  DF_features_col<-colnames(DF)[-c(1:length(meta_col))]
  ##
  table(DF$Group)
  
  
  ### create all combinations
  pair <- CombPairs(unique(DF$Group))
  
  for (j in 1:nrow(pair)) {
    # j=1
    DF2 <- subset(DF, Group == pair[j, 1] | Group == pair[j, 2])
    
    comparison <- paste0(pair[j, 1], "_", pair[j, 2])
    
    
    ## or ignore the subgroups
    DF.tss <- DF2
    # 计算每列的方差
    DF.tss[,-c(1:2)] <- DF.tss[,-c(1:2)]      %>% 
                          select_if(is.numeric)       %>%
                          select(where(~var(.x) != 0))
    # DF.tss <- DF2 %>%janitor::adorn_percentages("row") %>% data.frame()
    DF.tss <- DF.tss[, colMeans(is.na(DF.tss)) <= 0.3]
    DF.tss_features_col =colnames(DF.tss)
    
    
    ## real mean difference
    diff.real <- DF.tss %>%
      select_if(is.numeric) %>%
      map_df(~ MeanDiffCI(. ~ Group, data = DF.tss), .id = "var")
    ## model based mean difference
    diff.model <- DF.tss %>%
      select_if(is.numeric) %>%
      map_df(~ broom::tidy(wilcox.test(. ~ Group, data = DF.tss, conf.int = T, conf.level = .95)), .id = "var")
    diff.model$p.adj <- p.adjust(diff.model$p.value, "BH")
    
    ### merge
    diff <- left_join(diff.model, diff.real)
    #    diff.sig <- subset(diff, p.value < 0.01)
    diff.sig <- subset(diff, p.adj < 0.05)
    diff.sig$col <- ifelse(diff.sig$meandiff > 0, "01plus", "02minus")
    
    
    
    DF.tss2 <- reshape2::melt(  DF.tss[,  c("Group", DF.tss_features_col) ]  )
    #  DF.tss2$Group <- factor(DF.tss2$Group, levels = c("young", str_replace(comparison, pattern = "young-vs-", replacement = "")))
    DF.tss2$variable <- as.character(DF.tss2$variable)
    ##
    DF.tss3 <- DF.tss2[DF.tss2$variable %in% diff.sig$var, ]
    
    ## base plot
    p1 <- ggplot(data = DF.tss3, aes(x = variable, y = value, fill = Group))
    
    ### add alternated background
    for (i in 1:(nrow(diff.sig))) {
      p1 <- p1 + annotate("rect",
                          xmin = i - 0.5, xmax = i + 0.5, ymin = 0, ymax = Inf,
                          fill = ifelse(i %% 2 == 0, "white", "gray95")
      ) + scale_x_discrete()
    }
    
    p1 <- p1 +
      geom_bar(color = "black", position = "dodge", stat = "summary", fun = "mean") +
      xlab("") +
      ylab("Proportion (%)") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(color = "black"),
        axis.ticks = element_line(color = "black", linewidth = 0.5)
      ) +
      coord_flip() +
      theme(legend.position = "top")
    # Use custom colors
    p1 <- p1 + scale_fill_manual(values = c("darkorange", "deepskyblue"))
    
    
    ### base plot
    p2 <- ggplot(diff.sig, aes(x = var, y = meandiff, fill = col))
    ### add alternated background
    for (i in 1:(nrow(diff.sig))) {
      p2 <- p2 + annotate("rect",
                          xmin = i - 0.5, xmax = i + 0.5, ymin = -Inf, ymax = Inf,
                          fill = ifelse(i %% 2 == 0, "white", "gray95")
      ) + scale_x_discrete()
    }
    ## add the other components
    p2 <- p2 +
      geom_point(stat = "identity", size = 3, aes(col = col), position = position_dodge(width = 1)) +
      geom_errorbar(aes(ymin = lwr.ci, ymax = upr.ci),
                    width = .2, size = 0.5,
                    position = "dodge"
      ) +
      xlab("") +
      ylab("Difference between proportions (%)") +
      geom_hline(yintercept = 0, lty = 2) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(color = "black"),
        axis.ticks.x = element_line(color = "black", linewidth = 0.5)
      ) +
      coord_flip() +
      theme(legend.position = "none") +
      theme(axis.text.y = element_blank())
    p2 <- p2 + scale_color_manual(values = c("darkorange", "deepskyblue"))
    
    p1 + p2
    
    ###
    ##
    dir.create("output",recursive = T)
    ggsave(filename = glue::glue("output/", out.name, comparison, "_metagenomics_smallSize.pdf"), width = 12, height = 5.5)
    ggsave(filename = glue::glue("output/", out.name, comparison, "_metagenomics_mediumSize.pdf"), width = 12.5, height = 7.0)
    ggsave(filename = glue::glue("output/", out.name, comparison, "_metagenomics_bigSize.pdf"), width = 12.5, height = 8.0)
    
    #  }
  }
  
}


############################################################################################################
## meta
metadata<-read.csv("../../16s/input/ONC-2019-LZ-023.csv")
metadata$sample <- paste0("X",metadata$sample)
metadata <- metadata[,c("sample","Groups")]
rownames(metadata)<-metadata$sample

## kegg
# SamleID_lab<- read.csv("input/samples.tsv",sep = "\t") 
# SamleID_lab.ID= SamleID_lab$Reads_raw_R1
# SamleID_lab.ID= stringr::word( basename(SamleID_lab.ID),2,2,stringr::fixed("_") )

kegg.table.df <- readr::read_csv("input/genomes_KEGG_Module_count_grouped_step_coverage0.8_tss.csv") 
colnames(kegg.table.df)[1:2]=c("module","module_name")
kegg.table.df$name_id<- paste0(kegg.table.df$module,"_",kegg.table.df$module_name)
kegg.table.df =kegg.table.df %>% dplyr::select(name_id,everything()) %>% dplyr::select(-module,-module_name)
#metadata_ <-metadata[colnames(kegg.table.df)[-1],]
#commonid <- intersect(rownames(metadata),colnames(kegg.table.df)[-1])
select_rows <- rownames(metadata)
#kegg.table.df[-1] <- dplyr::select(kegg.table.df[-1],colnames(kegg.table.df[-1]) %in% select_rows)
selected_cols <- kegg.table.df[, colnames(kegg.table.df) %in% select_rows]
name_id  <- data.frame(kegg.table.df$name_id)
colnames(name_id)<-"name_id"
kegg.table.df <- cbind(name_id,selected_cols)
#顺序不同时，核对两个数据框的行名和列名是否一致
identical(sort(colnames(kegg.table.df)[-1]), sort(rownames(metadata)))
#colnames(kegg.table.df)[-1]==rownames(metadata)
colnames(metadata)<-c("sampleID", "Group")

# kegg.table.df$name_id <-kegg.table.df$name_id %>% make.names()
pathway_pdf(meta_info=metadata,pathway_df=kegg.table.df,meta_col=c("sampleID", "Group"),out.name="KEGG_" )









