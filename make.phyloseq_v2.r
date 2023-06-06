library(phyloseq)
library(dplyr)
atlas_wd_folder=""
abundance_file = paste0(atlas_wd_folder,"genomes/counts/median_coverage_genomes.parquet")

raw_tax<-read.table("genomes/taxonomy/gtdb_taxonomy.tsv",sep = "\t",row.names = 1,header = T)
raw_count <- arrow::read_parquet(abundance_file) %>% column_to_rownames(var = "index") %>% t() %>% as.data.frame()
raw_tax<-raw_tax %>% as.matrix()

phy0 <- phyloseq::phyloseq(otu_table( as.matrix(raw_count) ,taxa_are_rows = TRUE),tax_table(raw_tax))

saveRDS(phy0, "genomes/atlas_binning_phyloseq.rds")
