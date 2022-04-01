#This script performs hypergeometric tests of relative codon overrepresentation in
#ambisense coding sequence pairs from annotated Genbank accessions using the 'phyper' function in the R stats package. 
#Requires text file with list of Genbank accession numbers and FASTA file with coding sequence pairs as input.
#One of the sequences in each analyzed pair must contain the word 'complement' in the FASTA header.

#Required packages

library(Biostrings)
library(openxlsx)
library(dplyr)
library(Rcpp)

#Input files

acc <- file("accessions-list.txt", open="r") #File with list of accession numbers
RefSeq_Accessions<-readLines(acc)
fas <-readDNAStringSet("sequence-file.fasta") #FASTA file with ambisense coding sequence pairs

#Empty data frames to store phyper hypergeometric test results

NNT_vir_sense_p_val<- data.frame()
NNT_anti_sense_p_val <- data.frame()

NNA_vir_sense_p_val<- data.frame()
NNA_anti_sense_p_val <- data.frame()

NNC_vir_sense_p_val<- data.frame()
NNC_anti_sense_p_val <- data.frame()

NNG_vir_sense_p_val<- data.frame()
NNG_anti_sense_p_val <- data.frame()

#Obtaining parameters for codon overrepresentation phyper tests
#Run '?phyper' to learn more about command arguments

for (i in 1:length(RefSeq_Accessions)){
  
  genome_seq<-fas[grep(RefSeq_Accessions[i],names(fas))]
  
  vir_sense_seq<-genome_seq[-grep("comp",names(genome_seq))]
  
  vir_sense_codon_table<- trinucleotideFrequency(vir_sense_seq[grep(RefSeq_Accessions[i],names(vir_sense_seq))], step=3) # vir sense codon frequencies
  vir_sense_codon_table2<- subset(vir_sense_codon_table, select = -c(ATG,TGG,TAA,TAG,TGA)) # remove Met, Trp and STOP  
  
  vir_sense_all_codon<-rowSums(as.data.frame(vir_sense_codon_table2)) # k arg for all vir sense phyper tests
  
  anti_sense_seq<-genome_seq[grep("comp",names(genome_seq))]
  
  anti_sense_codon_table<- trinucleotideFrequency(anti_sense_seq[grep(RefSeq_Accessions[i],names(anti_sense_seq))], step=3) # anti sense codon frequencies
  anti_sense_codon_table2<- subset(anti_sense_codon_table, select = -c(ATG,TGG,TAA,TAG,TGA)) #remove Met, Trp and STOP
  
  anti_sense_all_codon<-rowSums(as.data.frame(anti_sense_codon_table2)) # k arg for all anti sense phyper tests
  
  #genome NNT counts
  
  vir_sense_NNT<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)=="T"]) #q arg for vir sense NNT phyper test
  anti_sense_NNT<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)=="T"]) 
  
  genome_NNT<-sum(vir_sense_NNT, anti_sense_NNT) # m for NNT phyper test
  
  #genome NNA counts
  
  vir_sense_NNA<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)=="A"])
  anti_sense_NNA<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)=="A"]) 
  
  genome_NNA<-sum(vir_sense_NNA, anti_sense_NNA) # m for NNA phyper test
  
  #genome NNA counts
  
  vir_sense_NNC<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)=="C"])
  anti_sense_NNC<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)=="C"]) 
  
  genome_NNC<-sum(vir_sense_NNC, anti_sense_NNC) # m for NNC phyper test
  
  #genome NNG counts
  
  vir_sense_NNG<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)=="G"])
  anti_sense_NNG<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)=="G"]) 
  
  genome_NNG<-sum(vir_sense_NNG, anti_sense_NNG) # m for NNG phyper test
  
  #genome non-NNT ending counts
  
  vir_sense_non_NNT<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)!="T"])
  anti_sense_non_NNT<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)!="T"]) 
  
  genome_non_NNT<-sum(vir_sense_non_NNT, anti_sense_non_NNT) # n for NNT phyper test
  
  #genome non-NNA ending counts
  
  vir_sense_non_NNA<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)!="A"])
  anti_sense_non_NNA<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)!="A"]) 
  
  genome_non_NNA<-sum(vir_sense_non_NNA, anti_sense_non_NNA) # n for NNA phyper test
  
  #genome non-NNC ending counts
  
  vir_sense_non_NNC<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)!="C"])
  anti_sense_non_NNC<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)!="C"]) 
  
  genome_non_NNC<-sum(vir_sense_non_NNC, anti_sense_non_NNC) # n for NNC phyper test
  
  
  #genome non-NNG ending counts
  
  vir_sense_non_NNG<-rowSums(as.data.frame(vir_sense_codon_table2)[substr(names(as.data.frame(vir_sense_codon_table2)),3,3)!="G"])
  anti_sense_non_NNG<-rowSums(as.data.frame(anti_sense_codon_table2)[substr(names(as.data.frame(anti_sense_codon_table2)),3,3)!="G"]) 
  
  genome_non_NNG<-sum(vir_sense_non_NNG, anti_sense_non_NNG) # n for NNG phyper test
  
  
  #phyper tests

for (i in 1:length(vir_sense_seq)){
  NNT_vir_sense_p_test<-phyper(vir_sense_NNT-1,genome_NNT,
                                        genome_non_NNT,vir_sense_all_codon, lower.tail = FALSE)
    
  NNT_vir_sense_result_df<- data.frame('NNT_vir_sense_p_val'=NNT_vir_sense_p_test)
  NNT_vir_sense_p_val<- bind_rows(NNT_vir_sense_p_val, NNT_vir_sense_result_df)
  
  NNA_vir_sense_p_test<-phyper(vir_sense_NNA-1,genome_NNA,
                               genome_non_NNA,vir_sense_all_codon, lower.tail = FALSE)
  
  NNA_vir_sense_result_df<- data.frame('NNA_vir_sense_p_val'=NNA_vir_sense_p_test)
  NNA_vir_sense_p_val<- bind_rows(NNA_vir_sense_p_val, NNA_vir_sense_result_df)
  
  NNC_vir_sense_p_test<-phyper(vir_sense_NNC-1,genome_NNC,
                               genome_non_NNC,vir_sense_all_codon, lower.tail = FALSE)
  
  NNC_vir_sense_result_df<- data.frame('NNC_vir_sense_p_val'=NNC_vir_sense_p_test)
  NNC_vir_sense_p_val<- bind_rows(NNC_vir_sense_p_val, NNC_vir_sense_result_df)
  
  NNG_vir_sense_p_test<-phyper(vir_sense_NNG-1,genome_NNG,
                               genome_non_NNG,vir_sense_all_codon, lower.tail = FALSE)
  
  NNG_vir_sense_result_df<- data.frame('NNG_vir_sense_p_val'=NNG_vir_sense_p_test)
  NNG_vir_sense_p_val<- bind_rows(NNG_vir_sense_p_val, NNG_vir_sense_result_df)
 }
  
for (i in 1:length(anti_sense_seq)){
  NNT_anti_sense_p_test<-phyper(anti_sense_NNT-1,genome_NNT,
                                genome_non_NNT,anti_sense_all_codon,lower.tail = FALSE)
  
  NNT_anti_sense_result_df<- data.frame('NNT_anti_sense_p_val'=NNT_anti_sense_p_test)
  NNT_anti_sense_p_val<- bind_rows(NNT_anti_sense_p_val, NNT_anti_sense_result_df)
  
  NNA_anti_sense_p_test<-phyper(anti_sense_NNA-1,genome_NNA,
                                genome_non_NNA,anti_sense_all_codon,lower.tail = FALSE)
  
  NNA_anti_sense_result_df<- data.frame('NNA_anti_sense_p_val'=NNA_anti_sense_p_test)
  NNA_anti_sense_p_val<- bind_rows(NNA_anti_sense_p_val, NNA_anti_sense_result_df)
  
  NNC_anti_sense_p_test<-phyper(anti_sense_NNC-1,genome_NNC,
                                genome_non_NNC,anti_sense_all_codon,lower.tail = FALSE)
  
  NNC_anti_sense_result_df<- data.frame('NNC_anti_sense_p_val'=NNC_anti_sense_p_test)
  NNC_anti_sense_p_val<- bind_rows(NNC_anti_sense_p_val, NNC_anti_sense_result_df)
  
  NNG_anti_sense_p_test<-phyper(anti_sense_NNG-1,genome_NNG,
                                genome_non_NNG,anti_sense_all_codon,lower.tail = FALSE)
  
  NNG_anti_sense_result_df<- data.frame('NNG_anti_sense_p_val'=NNG_anti_sense_p_test)
  NNG_anti_sense_p_val<- bind_rows(NNG_anti_sense_p_val, NNG_anti_sense_result_df)
 }
}

NNT_results<- cbind(RefSeq_Accessions, NNT_vir_sense_p_val, NNT_anti_sense_p_val)
NNA_results<- cbind(RefSeq_Accessions, NNA_vir_sense_p_val, NNA_anti_sense_p_val)
NNC_results<- cbind(RefSeq_Accessions, NNC_vir_sense_p_val, NNC_anti_sense_p_val)
NNG_results<- cbind(RefSeq_Accessions, NNG_vir_sense_p_val, NNG_anti_sense_p_val)

#Create empty wprkbook and add sheets

result_file<- createWorkbook()

addWorksheet(result_file, "NNT results")
addWorksheet(result_file, "NNA results")
addWorksheet(result_file, "NNG results")
addWorksheet(result_file, "NNC results")

#Write data to Excel sheets 

writeData(result_file, sheet ="NNT results", NNT_results)
writeData(result_file, sheet ="NNA results", NNA_results)
writeData(result_file, sheet ="NNG results", NNG_results)
writeData(result_file, sheet ="NNC results", NNC_results)

#Export file

saveWorkbook(result_file, file = "output-file.xlsx")
