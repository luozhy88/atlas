library(phyloseq)
library(dplyr)
raw_count<-read.table("genomes/counts/raw_counts_genomes.tsv",row.names = 1,header = T)
raw_tax<-read.table("genomes/taxonomy/gtdb_taxonomy.tsv",sep = "\t",row.names = 1,header = T)
# colnames(raw_tax)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
raw_tax<-raw_tax %>% as.matrix()
phy0 <- phyloseq::phyloseq(otu_table( as.matrix(raw_count) ,taxa_are_rows = TRUE),tax_table(raw_tax))


###############
library(phyloseq)
library(Biostrings)
library(RDPutils)
rep.seqs <- Biostrings::readDNAStringSet(paste0(workpath,"/representative_sequences/filtered/sequences.fasta"), format = "fasta")
phyloseq_new <- phyloseq::phyloseq(otu.table, rep.seqs)
phyloseq_new
tax_table(phyloseq_new)<-tax_table(physeq)
#phy_tree(phyloseq_new)<-phy_tree(physeq)
phyloseq_new



Biostrings::writeXStringSet(rep.seqs,file ="rep.seqs.fasta",   format = "fasta")
phyloseq_new
print(phyloseq_new)
saveRDS(phyloseq_new, paste0(workpath,"/atlas_binning_ampliseq_phyloseq.rds"))
