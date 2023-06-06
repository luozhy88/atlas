library(dplyr)
library(tibble)
library(tidyr)
library(foreach)
library(doParallel)
library(arrow)
library(rhdf5)
library(readr)


# genome module

## read
atlas_wd_folder=""
taxonomy_file = paste0(atlas_wd_folder,"genomes/taxonomy/gtdb_taxonomy.tsv")
tree_file = paste0(atlas_wd_folder,"genomes/tree/gtdbtk.bac120.nwk")
quality_file= paste0(atlas_wd_folder,"genomes/checkm/completeness.tsv")
counts_file= paste0(atlas_wd_folder,"genomes/counts/counts_genomes.parquet")
abundance_file = paste0(atlas_wd_folder,"genomes/counts/median_coverage_genomes.parquet")
readstats_file= paste0(atlas_wd_folder,"stats/read_counts.tsv")
keggmodules_file = paste0(atlas_wd_folder,"genomes/annotations/dram/kegg_modules.tsv")
product_file_name = paste0(atlas_wd_folder,"genomes/annotations/dram/distil/product.tsv")


## function
kegg_modules= readr::read_tsv(keggmodules_file,col_select = -1,show_col_types = FALSE)
module_names= kegg_modules %>% dplyr::select(c('module','module_name')) %>% distinct()
# convert to dataframe with row names
module_names=  data.frame(row.names= module_names$module, name= module_names$module_name)


### Relative abundance

#For the relative abundance, we take the coverage over the genome, not the raw counts. This implicitly normalizes for genome size. The coverage is calculated as the median of the coverage values calculated in 1kb blocks.
D <- arrow::read_parquet(abundance_file) %>%
  column_to_rownames(var = "index") %>%
  as.data.frame()
# calculate relative abundance
rel_ab <- sweep(D, 1, rowSums(D),`/`)




### coverage_threshold
step_coverage_threshold= 0.8

module_step_coverage_matrix = pivot_wider(kegg_modules,  
                                          id_cols = genome,
                                          names_from = module,
                                          values_from = step_coverage
) %>%
  column_to_rownames("genome") %>% as.matrix()

module_step_coverage_matrix = module_step_coverage_matrix[, colSums(module_step_coverage_matrix) > 0]


module_presence_matrix = 1 * (module_step_coverage_matrix>step_coverage_threshold)

module_presence_matrix = module_presence_matrix[, colSums(module_presence_matrix) > 0]

#Sum of rel_ab for all species where a module is presence is equel to the matrix multiplication
stopifnot(dim(rel_ab)[2] == dim(module_presence_matrix)[1]  )
module_rel_ab <- as.matrix(rel_ab) %*% module_presence_matrix
module_rel_ab <- data.frame(t(module_rel_ab))

module_rel_ab_names=merge(module_names,module_rel_ab,by=0) 

write.csv(module_rel_ab_names,file = "genomes/annotations/genomes_KEGG_Module_count_grouped_step_coverage0.8_tss.csv",row.names = F,quote = T)



# genomes metabolism

metabolism_file_name = paste0(atlas_wd_folder,"genomes/annotations/dram/distil/metabolism_summary.xlsx")
list_sheetname=c("MISC","carbon utilization","Transporters","Energy","Organic Nitrogen","carbon utilization (Woodcroft)")
tables_all=list()
for (sheetname in list_sheetname){
    print(sheetname)
    # sheetname="MISC"
    metabolism_file=readxl::read_excel(metabolism_file_name,sheetname) 
    basis_info=metabolism_file[,1:4] %>% dplyr::distinct(gene_id,.keep_all = T)%>% column_to_rownames("gene_id")
    basis_info$class=sheetname
    metabolism_count=metabolism_file[,-c(2:5)]
    metabolism_count[,-1]<-apply(metabolism_count[,-1], 2, as.numeric)  %>% data.frame()
    metabolism_count=metabolism_count %>% dplyr::group_by(gene_id) %>% summarise_all(.funs = sum) %>% column_to_rownames("gene_id") %>% t() %>% data.frame()
    
    rel_ab_metab=rel_ab[,rownames(metabolism_count)]
    # Sum of rel_ab for all species where a module is presence is equel to the matrix multiplication
    stopifnot(dim(rel_ab_metab)[2] == dim(metabolism_count)[1]  )
    metab_rel_ab <- as.matrix(rel_ab_metab) %*% as.matrix(metabolism_count)#the rownames of metabolism_count is MAG
    metab_rel_ab <- data.frame(t(metab_rel_ab))
    metab_rel_ab_basis=merge(basis_info,metab_rel_ab,by=0) 
    tables_all[[sheetname]]=metab_rel_ab_basis
}
tables_all_metabolism<-plyr::rbind.fill(tables_all)

# write.csv(tables_all_metabolism,file = "genomes/annotations/metabolism_summary_merge_all.csv",row.names = F,quote = T)

# product
product_file= readr::read_tsv(product_file_name,show_col_types = FALSE) %>% column_to_rownames("genome")

product_file_presence =1*product_file
product_file_presence = 1 * (product_file_presence>step_coverage_threshold)
product_file_presence = product_file_presence[, colSums(product_file_presence) > 0]

#Sum of rel_ab for all species where a module is presence is equel to the matrix multiplication
stopifnot(dim(rel_ab)[2] == dim(product_file_presence)[1]  )
product_file_presence_rel_ab <- as.matrix(rel_ab) %*% product_file_presence
product_file_presence_rel_ab <- data.frame(t(product_file_presence_rel_ab))

write.csv(product_file_presence_rel_ab,file = "genomes/annotations/product_file_presence_rel_ab_tss.csv",row.names = T,quote = T)
