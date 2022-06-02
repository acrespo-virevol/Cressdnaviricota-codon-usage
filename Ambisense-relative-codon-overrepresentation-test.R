#This script performs hypergeometric tests of relative codon overrepresentation in
#ambisense coding sequence pairs from annotated Genbank accessions. Requires text file
#with list of Genbank accession numbers and FASTA file with coding sequence pairs as input.
#One of the sequences in each analyzed pair must contain the word 'complement' in the FASTA header.
#The example files used here are available at: https://github.com/acrespo-virevol/Cressdnaviricota-codon-usage

#Required packages

library(Biostrings)
library(openxlsx)
library(dplyr)
library(Rcpp)

#Input files

acc <- file("Cyclovirus_RefSeq_Accessions_n42.txt", open="r") #File with list of accession numbers
RefSeq_Accessions<-readLines(acc)
fas <-readDNAStringSet("Cyclovirus_RefSeqs_Rep-and-CP_n42.fasta") #FASTA file with ambisense coding sequence pairs

#Empty data frames to store phyper hypergeometric test results

nnt_vir_sense_p_val<- data.frame()
nnt_anti_sense_p_val <- data.frame()

nna_vir_sense_p_val<- data.frame()
nna_anti_sense_p_val <- data.frame()

nnc_vir_sense_p_val<- data.frame()
nnc_anti_sense_p_val <- data.frame()

nng_vir_sense_p_val<- data.frame()
nng_anti_sense_p_val <- data.frame()

#Obtaining parameters for codon overrepresentation phyper tests
#Run '?phyper' to learn more about command arguments

for (i in 1:length(RefSeq_Accessions)){
  
  genome_seq<-fas[grep(RefSeq_Accessions[i],names(fas))]
  
  vir_sense_seq<-genome_seq[-grep("complement",names(genome_seq))]
  
  vir_sense_codon_table<- trinucleotideFrequency(vir_sense_seq[grep(RefSeq_Accessions[i],names(vir_sense_seq))], step=3) # vir sense codon frequencies
  vir_sense_codon_table2<- subset(vir_sense_codon_table, select = -c(ATG,TGG,TAA,TAG,TGA)) # remove Met, Trp and STOP  
  
  vir_sense_all_codon<-rowSums(as.data.frame(vir_sense_codon_table2)) # k arg for all vir sense phyper tests
  
  anti_sense_seq<-genome_seq[grep("complement",names(genome_seq))]
  
  anti_sense_codon_table<- trinucleotideFrequency(anti_sense_seq[grep(RefSeq_Accessions[i],names(anti_sense_seq))], step=3) # anti sense codon frequencies
  anti_sense_codon_table2<- subset(anti_sense_codon_table, select = -c(ATG,TGG,TAA,TAG,TGA)) #remove Met, Trp and STOP
  
  anti_sense_all_codon<-rowSums(as.data.frame(anti_sense_codon_table2)) # k arg for all anti sense phyper tests
  
  #genome nnt counts
  
  vir_sense_nnt<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)=="T"]) #q arg for vir sense nnt phyper test
  anti_sense_nnt<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)=="T"]) 
  
  genome_nnt<-sum(vir_sense_nnt, anti_sense_nnt) # m for nnt phyper test
  
  #genome nna counts
  
  vir_sense_nna<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)=="A"])
  anti_sense_nna<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)=="A"]) 
  
  genome_nna<-sum(vir_sense_nna, anti_sense_nna) # m for nna phyper test
  
  #genome nna counts
  
  vir_sense_nnc<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)=="C"])
  anti_sense_nnc<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)=="C"]) 
  
  genome_nnc<-sum(vir_sense_nnc, anti_sense_nnc) # m for nnc phyper test
  
  #genome nng counts
  
  vir_sense_nng<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)=="G"])
  anti_sense_nng<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)=="G"]) 
  
  genome_nng<-sum(vir_sense_nng, anti_sense_nng) # m for nng phyper test
  
  #genome non-nnt ending counts
  
  vir_sense_non_nnt<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)!="T"])
  anti_sense_non_nnt<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)!="T"]) 
  
  genome_non_nnt<-sum(vir_sense_non_nnt, anti_sense_non_nnt) # n for nnt phyper test
  
  #genome non-nna ending counts
  
  vir_sense_non_nna<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)!="A"])
  anti_sense_non_nna<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)!="A"]) 
  
  genome_non_nna<-sum(vir_sense_non_nna, anti_sense_non_nna) # n for nna phyper test
  
  #genome non-nnc ending counts
  
  vir_sense_non_nnc<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)!="C"])
  anti_sense_non_nnc<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)!="C"]) 
  
  genome_non_nnc<-sum(vir_sense_non_nnc, anti_sense_non_nnc) # n for nnc phyper test
  
  
  #genome non-nng ending counts
  
  vir_sense_non_nng<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)!="G"])
  anti_sense_non_nng<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)!="G"]) 
  
  genome_non_nng<-sum(vir_sense_non_nng, anti_sense_non_nng) # n for nng phyper test
  
  
  #phyper tests

for (i in 1:length(vir_sense_seq)){
  nnt_vir_sense_p_test<-phyper(vir_sense_nnt-1,genome_nnt,
                                        genome_non_nnt,vir_sense_all_codon, lower.tail = FALSE)
    
  nnt_vir_sense_result_df<- data.frame('nnt_vir_sense_p_val'=nnt_vir_sense_p_test)
  nnt_vir_sense_p_val<- bind_rows(nnt_vir_sense_p_val, nnt_vir_sense_result_df)
  
  nna_vir_sense_p_test<-phyper(vir_sense_nna-1,genome_nna,
                               genome_non_nna,vir_sense_all_codon, lower.tail = FALSE)
  
  nna_vir_sense_result_df<- data.frame('nna_vir_sense_p_val'=nna_vir_sense_p_test)
  nna_vir_sense_p_val<- bind_rows(nna_vir_sense_p_val, nna_vir_sense_result_df)
  
  nnc_vir_sense_p_test<-phyper(vir_sense_nnc-1,genome_nnc,
                               genome_non_nnc,vir_sense_all_codon, lower.tail = FALSE)
  
  nnc_vir_sense_result_df<- data.frame('nnc_vir_sense_p_val'=nnc_vir_sense_p_test)
  nnc_vir_sense_p_val<- bind_rows(nnc_vir_sense_p_val, nnc_vir_sense_result_df)
  
  nng_vir_sense_p_test<-phyper(vir_sense_nng-1,genome_nng,
                               genome_non_nng,vir_sense_all_codon, lower.tail = FALSE)
  
  nng_vir_sense_result_df<- data.frame('nng_vir_sense_p_val'=nng_vir_sense_p_test)
  nng_vir_sense_p_val<- bind_rows(nng_vir_sense_p_val, nng_vir_sense_result_df)
 }
  
for (i in 1:length(anti_sense_seq)){
  nnt_anti_sense_p_test<-phyper(anti_sense_nnt-1,genome_nnt,
                                genome_non_nnt,anti_sense_all_codon,lower.tail = FALSE)
  
  nnt_anti_sense_result_df<- data.frame('nnt_anti_sense_p_val'=nnt_anti_sense_p_test)
  nnt_anti_sense_p_val<- bind_rows(nnt_anti_sense_p_val, nnt_anti_sense_result_df)
  
  nna_anti_sense_p_test<-phyper(anti_sense_nna-1,genome_nna,
                                genome_non_nna,anti_sense_all_codon,lower.tail = FALSE)
  
  nna_anti_sense_result_df<- data.frame('nna_anti_sense_p_val'=nna_anti_sense_p_test)
  nna_anti_sense_p_val<- bind_rows(nna_anti_sense_p_val, nna_anti_sense_result_df)
  
  nnc_anti_sense_p_test<-phyper(anti_sense_nnc-1,genome_nnc,
                                genome_non_nnc,anti_sense_all_codon,lower.tail = FALSE)
  
  nnc_anti_sense_result_df<- data.frame('nnc_anti_sense_p_val'=nnc_anti_sense_p_test)
  nnc_anti_sense_p_val<- bind_rows(nnc_anti_sense_p_val, nnc_anti_sense_result_df)
  
  nng_anti_sense_p_test<-phyper(anti_sense_nng-1,genome_nng,
                                genome_non_nng,anti_sense_all_codon,lower.tail = FALSE)
  
  nng_anti_sense_result_df<- data.frame('nng_anti_sense_p_val'=nng_anti_sense_p_test)
  nng_anti_sense_p_val<- bind_rows(nng_anti_sense_p_val, nng_anti_sense_result_df)
 }
}

nnt_results<- cbind(RefSeq_Accessions, nnt_vir_sense_p_val, nnt_anti_sense_p_val)
nna_results<- cbind(RefSeq_Accessions, nna_vir_sense_p_val, nna_anti_sense_p_val)
nnc_results<- cbind(RefSeq_Accessions, nnc_vir_sense_p_val, nnc_anti_sense_p_val)
nng_results<- cbind(RefSeq_Accessions, nng_vir_sense_p_val, nng_anti_sense_p_val)

#Create empty workbook and add sheets

result_file<- createWorkbook()

addWorksheet(result_file, "nnt results")
addWorksheet(result_file, "nna results")
addWorksheet(result_file, "nng results")
addWorksheet(result_file, "nnc results")

#Write data to Excel sheets 

writeData(result_file, sheet ="nnt results", nnt_results)
writeData(result_file, sheet ="nna results", nna_results)
writeData(result_file, sheet ="nng results", nng_results)
writeData(result_file, sheet ="nnc results", nnc_results)

#Export file

saveWorkbook(result_file, file = "output-file.xlsx")
