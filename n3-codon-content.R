#This script calculates nnt, nna, nng, and nnc codon content
#(where n= any nucleotide, t=thymine, a= adenine, g= guanine and c=cytosine)
#for each coding sequence in an input FASTA file.

#Required packages

library(Biostrings)
library(openxlsx)
library(dplyr)

#Input FASTA file

fas<-readDNAStringSet("sequence-file.fasta") #input FASTA file with coding sequences

Sequence<- names(fas) #Sequence names

#Empty data frames for n3 content

nnt<- data.frame()
nna<- data.frame()
nng<- data.frame()
nnc<- data.frame()


for (i in 1:length(fas)){
  
  codon_table<- trinucleotideFrequency(fas[i], step=3, with.labels= TRUE) # codon frequencies
  all_codon<-rowSums(as.data.frame(codon_table)) #sequence length in codons
  
  t3_count <-rowSums(as.data.frame(codon_table)[substr(names(as.data.frame(codon_table)),3,3)=="T"])
  a3_count<- rowSums(as.data.frame(codon_table)[substr(names(as.data.frame(codon_table)),3,3)=="A"])
  g3_count<- rowSums(as.data.frame(codon_table)[substr(names(as.data.frame(codon_table)),3,3)=="G"])
  c3_count<- rowSums(as.data.frame(codon_table)[substr(names(as.data.frame(codon_table)),3,3)=="C"])
  
  
  t3_content<- (t3_count/all_codon)*100
  a3_content<- (a3_count/all_codon)*100
  g3_content<- (g3_count/all_codon)*100
  c3_content<- (c3_count/all_codon)*100
  
 
  t3_cont_df<- data.frame('nnt'=t3_content)
  nnt<- bind_rows(nnt, t3_cont_df)
  
  a3_cont_df<- data.frame('nna'=a3_content)
  nna<- bind_rows(nna, a3_cont_df)
  
  g3_cont_df<- data.frame('nng'=g3_content)
  nng<- bind_rows(nng, g3_cont_df)
  
  c3_cont_df<- data.frame('nnc'=c3_content)
  nnc<- bind_rows(nnc, c3_cont_df)
  
  
}
 
n3_content<- cbind(Sequence, nnt, nna, nng, nnc)
 
 result_file<- createWorkbook()
 addWorksheet(result_file, "n3 content")
 
 writeData(wb= result_file, sheet= "n3 content", x= n3_content)
 
 #Export file
 
 saveWorkbook(result_file, file = "output-file-name.xlsx")