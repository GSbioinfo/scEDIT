library(tidyverse)
library(patchwork)
library(stringr)
library(patchwork)
library(DBI)
library(RSQLite)
library(here)
library(extrafont)
library(Biostrings)
source("ploting_functions.R")
ampgRNA_details<- readxl::read_excel("/scratch/PRJNA889818/gRNA_AMP_details.xlsx") %>%
  mutate(AMPID=gsub("\\$.*","",AMP_gRNA_ID))
csv_file_path= list(c(  "gatarep1",
                       "/scratch/PRJNA889818/allSampRun001/GATA_rep1/gata1_rep1/result_out/final_count/For_cellBC_gRNA_AMP_seqcount_freq.csv"),
                    c("gatarep2",
                      "/scratch/PRJNA889818/allSampRun001/GATA_rep2/gata1_rep2/result_out/final_count/For_cellBC_gRNA_AMP_seqcount_freq.csv"),
                    c("gatarep3",
                      "/scratch/PRJNA889818/allSampRun001/GATA_rep3/gata1_rep3/result_out/final_count/For_cellBC_gRNA_AMP_seqcount_freq.csv"),
                    c("gatarep4",
                      "/scratch/PRJNA889818/allSampRun001/GATA_rep4/gata1_rep4/result_out/final_count/For_cellBC_gRNA_AMP_seqcount_freq.csv"),
                    c("hbfrep1",
                      "/scratch/PRJNA889818/allSampRun001/HbF_rep1/HbF_rep1/result_out/final_count/For_cellBC_gRNA_AMP_seqcount_freq.csv"),
                    c("hbfrep2",
                      "/scratch/PRJNA889818/allSampRun001/HbF_rep2/HbF_rep2/result_out/final_count/For_cellBC_gRNA_AMP_seqcount_freq.csv"))
FColNames<- c("AMPNAME","CELLID","gRNASeq","gRNASeqRC", "SEQ", "seqcount", "freq", "namefreq")
RColNames<- c("AMPNAME","CELLID","gRNASeq","gRNASeqRC", "SEQ", "seqcount", "freq", "namefreq")

all_gRNAseq_qulify <- tibble()
tidymsa_all <- tibble()
all_edits_stats<- tibble()
all_cell_qulify<- tibble()
all_full_data_edit<- tibble()
for(i in seq(1,4)){
for_data<- read.delim(csv_file_path[[i]][2], sep = ',', header = F) %>%
  set_names(FColNames) %>%
  mutate(AMP_CELL_ID= paste0(gsub("!Edite.*","",AMPNAME),"_",CELLID), READN=c('FOR'))%>%
  select( READN,CELLID, AMP_CELL_ID,gRNASeq,gRNASeqRC,SEQ,seqcount,freq,namefreq) 

full_data<- read.delim(gsub("For_cellBC_", "Rev_cellBC_", csv_file_path[[i]][2]), sep = ',', header = F) %>%
   set_names(RColNames) %>%
  mutate(AMP_CELL_ID= paste0(gsub("!Edite.*","",AMPNAME),"_",CELLID), READN=c('REV'))%>%
  select( READN, CELLID, AMP_CELL_ID,gRNASeq,gRNASeqRC,SEQ,seqcount,freq,namefreq) %>% bind_rows(for_data)


fuul_data_grp<- full_data %>% group_by(READN,AMP_CELL_ID,CELLID) %>%
  group_split()
full_data_edit<- tibble()
for(grp in fuul_data_grp){
  plttitle<- paste0(unique(grp$CELLID),":",unique(grp$gRNASeq),":",unique(grp$gRNASeqRC),":",unique(grp$READN))
  if(nrow(grp)<2){
    next
  }
  tmpref<- grp %>% filter(grepl("Aref_",namefreq)) %>%
    group_by(READN, gRNASeq, SEQ) %>% 
    summarise(seqcount=sum(seqcount)) %>%
    ungroup() %>%
    mutate(freq = 0) %>%
    mutate(namefreq = paste0('Aref_',READN))
  if(i > 4){
  tmpedits<- grp %>% filter(!grepl("Aref_",namefreq)) %>%
    filter(SEQ!=tmpref$SEQ & !str_ends(SEQ, "-")) %>%
    filter(sum(freq)>=25 & sum(seqcount) > 20)
  }else{
    tmpedits<- grp %>% filter(!grepl("Aref_",namefreq)) %>%
      filter(SEQ!=tmpref$SEQ & !str_ends(SEQ, "-")) %>%
      filter(sum(freq)>=25 & sum(seqcount) > 2)
  }
  editStat = 'N'
  if(nrow(tmpedits)<1){
    editStat = 'N'
    tmpgrpedit<- grp %>% mutate(EDTISTATUS=c(editStat))
    full_data_edit <- bind_rows(tmpgrpedit,full_data_edit)
    next
  }
  
  editStat = 'Y'
  tmpgrpedit<- grp %>% mutate(EDTISTATUS=c(editStat))
  full_data_edit <- bind_rows(tmpgrpedit,full_data_edit)
  
}

all_full_data_edit<- full_data_edit %>% 
  mutate(sample=c(csv_file_path[[i]][1])) %>%
  bind_rows(all_full_data_edit)

 edited_cell_qulify<- full_data_edit %>% filter(EDTISTATUS=="Y") %>%
  group_by(READN,AMP_CELL_ID,CELLID) %>%
  filter(sum(seqcount)>5) %>% 
  ungroup()

  cell_qulify<- full_data_edit %>% 
    group_by(READN,AMP_CELL_ID,CELLID) %>%
    filter(sum(seqcount)>5) %>% 
    ungroup()
  
  all_cell_qulify <- cell_qulify %>% mutate(sample=c(csv_file_path[[i]][1])) %>%
    bind_rows(all_cell_qulify)
  
  edits_stats_list<- cell_qulify %>% 
   group_by(CELLID) %>% 
   slice_max(seqcount,n=1,with_ties = F) %>%
   ungroup()%>%
   mutate(HAPLO=case_when((EDTISTATUS == 'Y' & freq <= 60 ) ~ "HETERO",
                          (EDTISTATUS == 'Y' & freq > 60) ~ "HOMO",
                          EDTISTATUS == 'N' ~ "UNEDIT")) %>%
  group_by(gRNASeq, READN) %>%
  group_split()
  
 edits_stats <- tibble()
 for(edt_stat in edits_stats_list){
   edits_stats <- edt_stat %>% group_by(gRNASeq, READN) %>%
     count(HAPLO, EDTISTATUS, sort = T) %>%
     group_by(EDTISTATUS) %>%
     mutate(EDCOUNT= sum(n)) %>%
     ungroup()%>%
     group_by(gRNASeq,READN) %>%
     mutate(CELLCOUNT = sum(n)) %>%
     ungroup()%>%
     bind_rows(edits_stats)
 } 
 edits_stats <- edits_stats %>% arrange(CELLCOUNT) %>%
   mutate(gRNASeq = factor(gRNASeq, levels = unique(gRNASeq)),
          sample = c(csv_file_path[[i]][1]))
 
 all_edits_stats<- edits_stats %>% bind_rows(all_edits_stats)
 
analysi_cell_qulify <- edited_cell_qulify %>% group_by(CELLID) %>%
  slice_max(seqcount, n = 1, with_ties = F) %>%
  distinct(CELLID,.keep_all = T) %>%
  ungroup() %>%
  mutate(EDITS = case_when(freq >= 60 ~ "Homozygous",
                           .default = "Heterozygous"))


uniq_gRNAseqs<- full_data_edit %>% group_by(gRNASeq) %>% filter(EDTISTATUS=='Y') %>%
  slice_max(seqcount, n=1, with_ties = F) %>%
  distinct(gRNASeq,.keep_all = T) %>%
  ungroup() %>%
  mutate(EDITS = case_when(freq >= 60 ~ "Homozygous",
                           .default = "Heterozygous"))

uniq_cellIDs<- full_data_edit %>% group_by(CELLID) %>%
  distinct(CELLID,.keep_all = T) %>%
  ungroup() 


analysi_gRNAseq_qulify <- analysi_cell_qulify %>% filter(EDTISTATUS=="Y") %>%
  group_by(gRNASeq) %>%
  ungroup() %>%
  count(gRNASeq)
top10_gRNAseq_qulify <- analysi_cell_qulify %>% filter(EDTISTATUS=="Y") %>%
  group_by(gRNASeq) %>%
  ungroup() %>%
  count(gRNASeq) %>%
  slice_max(n, n=10,with_ties = F)
  
  
dataplt_top_gRNA_grp <- full_data_edit %>% filter(EDTISTATUS=='Y') %>%
  filter(gRNASeq %in% as.character(top10_gRNAseq_qulify$gRNASeq)) %>%
  group_by(gRNASeq, READN, CELLID ) %>%
  mutate(seqsum = sum(seqcount), sumfreq = sum(freq)) %>%
  ungroup()%>%
  group_by(gRNASeq)%>%
  group_split()
temp_tidymsa_all<- tibble()

edigRNA_plt_list<- ggplot()

for(grp in dataplt_top_gRNA_grp){
  if(nrow(grp)<2){
    next
  }
  
  tmgrp <- grp %>% filter(sum(seqcount)>10)
  if(nrow(tmgrp)<1){
    next
  }
  tmgrp <- grp %>% group_by(READN,CELLID) %>%
    mutate(sumseq = sum(seqcount)) %>%
    ungroup()%>%
    filter(sumseq>min(sumseq)) %>%
    select(names(grp)) %>%
    group_by(CELLID) %>%
    group_split()
  
  
  
  #editplt<- plt_editseqs(temp_tidymsa, plttitle, T)
  #edigRNA_plt_list<- edigRNA_plt_list+editplt
  grna_cell_collect<- tibble()
  
  for(newgrp in tmgrp){
    plttitle<- paste0(unique(newgrp$CELLID),":",unique(newgrp$gRNASeq),":",unique(newgrp$gRNASeqRC),":",unique(newgrp$READN))
    if(length(plttitle)>1){
      ref_tmp<-  newgrp %>% filter(READN =="FOR") %>% filter(grepl("ref_",namefreq)) 
      
    }else{
      temp_tidymsa <- transform_data_fn ( newgrp, plttitle )
    }
  }
  for(newgrp in tmgrp){
    
    plttitle<- paste0(unique(newgrp$CELLID),":",unique(newgrp$gRNASeq),":",unique(newgrp$gRNASeqRC),":",unique(newgrp$READN))
    if(length(plttitle)>1){
      temp_tidymsa <- transform_data_fn ( newgrp %>% filter(READN =="FOR"), plttitle[1] )
    }else{
      temp_tidymsa <- transform_data_fn ( newgrp, plttitle )
    }
    ref_tidy<- temp_tidymsa %>% filter(grepl("ref_",freqX)) %>%
      mutate(gRNASeq= c(unique(newgrp$gRNASeq)), cellid = c(unique(newgrp$CELLID)),sample = c(csv_file_path[[i]][1])) %>%
      mutate(freqX = paste0("000",freqX)) 
    temp_tidymsa_all<- temp_tidymsa %>% 
      mutate(gRNASeq= c(unique(newgrp$gRNASeq)), cellid = c(unique(newgrp$CELLID)),sample = c(csv_file_path[[i]][1])) %>%
      mutate(freqX = paste0(cellid,":",freqX)) %>%
     bind_rows(temp_tidymsa_all)
    #temp_tidymsa_all<- temp_tidymsa_all %>% filter(!grepl("ref_",freqX)) 
    
    #temp_tidymsa_all<- bind_rows(temp_tidymsa_all, ref_tidy)
  }
    
  
  #print(editplt)
}
tidymsa_all <- temp_tidymsa_all %>% #filter(!grepl("ref_",freqX)) %>%
  #mutate(freqX = paste0(cellid,":",freqX)) %>%
  #bind_rows(temp_tidymsa %>% filter(grepl("ref_",freqX)) %>%
  #            mutate(freqX = c("0ref_000"), cellid = c("CELLID"),sample = c(csv_file_path[[i]][1]))) %>%
  bind_rows(tidymsa_all)

analysi_gRNAseq_qulify$sample<- c(csv_file_path[[i]][1])
all_gRNAseq_qulify <- bind_rows(all_gRNAseq_qulify,analysi_gRNAseq_qulify)
}

 ggplot((all_gRNAseq_qulify %>% filter(n>1) %>%
          mutate(sample=gsub("hbf","",sample))),aes(gRNASeq,n)) + 
  geom_bar(stat = "identity",aes(color=sample,fill = sample)) +
 
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

 Most_prevelent_gRNA<- all_edits_stats %>% 
   group_by(gRNASeq) %>%
   slice_max(CELLCOUNT,n=1, with_ties = F) %>%
   ungroup() %>%
   group_by(gRNASeq) %>%
   mutate(sum(CELLCOUNT))

 fullseq_ammps<- getFullAmps()
 grna_ampid_table<- data.frame(matrix(ncol = 3,nrow = 0))
 names(grna_ampid_table)<- c("gRNASeq", "AMPID", "STRAND")
for(irow in 1:nrow(Most_prevelent_gRNA)){
  ampstest = F
  sgrna<- as.character(Most_prevelent_gRNA[irow,"gRNASeq"][[1]][1])
  for(xni in 1:nrow(fullseq_ammps)){
    ampid<- fullseq_ammps[xni,"AMPID"]
    seqamp<- fullseq_ammps[xni,"seq"]
    seqampRC <- fullseq_ammps[xni,"seqRC"]
    seqamplen <- fullseq_ammps[xni,"seq_len"]
    if(grepl(sgrna,seqamp) || grepl(sgrna,seqampRC) ){
      ampstest = T
      if(grepl(sgrna,seqamp)){
        tdf<- data.frame( sgrna, ampid, "+")
        names(tdf)<- names(grna_ampid_table)
        grna_ampid_table <-  rbind.data.frame(grna_ampid_table,tdf) 
      }
      if(grepl(sgrna,seqampRC)){
        tdf<- data.frame( sgrna, ampid, "-")
        names(tdf)<- names(grna_ampid_table)
        grna_ampid_table <-  rbind.data.frame(grna_ampid_table,tdf) 
        
      }
      
    }
  }
  if(!ampstest){
    tdf<- data.frame( sgrna, "NT", ".")
    names(tdf)<- names(grna_ampid_table)
    grna_ampid_table <-  rbind.data.frame(grna_ampid_table,tdf) 
  }
}
 grna_ampid_table<- grna_ampid_table %>% separate(AMPID,c('Target',"LOCA"), sep = "_GS:",remove = F)
 grna_ampid_table<- grna_ampid_table %>% 
   left_join(ampgRNA_details %>% select(Chromosome,AmpStart,AmpEnd,Target),by = 'Target') %>%
   group_by(gRNASeq) %>%
   distinct(gRNASeq, .keep_all = T) %>%
   mutate(sgRNA_id= paste0(gRNASeq,"_guide"))

#write_csv(grna_ampid_table%>%select(gRNASeq), "GATA1_gRNASeq.txt")
 if(i < 5){
   write_csv(grna_ampid_table, "GATA1_gRNASeq_chrom_table.csv")
   }else{
write_csv(grna_ampid_table, "HBF1_gRNASeq_chrom_table.csv")
}
gRNA_levels<- all_edits_stats %>% 
  group_by(gRNASeq)%>%
  slice_max(CELLCOUNT,n=1,with_ties = F) %>%
  arrange(CELLCOUNT)
 all_edits_stats$gRNASeq<- factor(all_edits_stats$gRNASeq, levels = gRNA_levels$gRNASeq)
 all_edits_stats$HAPLO<- factor(all_edits_stats$HAPLO, levels = c('UNEDIT','HOMO','HETERO'))
grna_cell_edit_plt<- ggplot(all_edits_stats, aes(gRNASeq, n)) +
  geom_bar(stat = "identity",aes(color = HAPLO, fill = HAPLO)) +
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) +
  theme(panel.spacing.x = unit(20, "pt")) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank()) +
  coord_flip()+
  xlab("")+ylab("")+
  facet_wrap(.~sample,  nrow =1)
if(i < 5){
ggsave("gRNA_editing_cellNumbers_GATA1_repeats.pdf", grna_cell_edit_plt, device = "pdf", height = 12, width = 10)
}else{
  ggsave("gRNA_editing_cellNumbers_HBF1_repeats.pdf", grna_cell_edit_plt, device = "pdf", height = 12, width = 10) 
}
topgRNAs<- all_edits_stats %>% filter(EDTISTATUS=='Y') %>%
  group_by(sample)%>%
  top_n(n = 20, wt = EDCOUNT) %>%
  ungroup() %>%
  arrange(desc(EDCOUNT))
top_gRNA_levels<- topgRNAs %>% 
  group_by(gRNASeq)%>%
  slice_max(EDCOUNT,n=1,with_ties = F) %>%
  arrange(EDCOUNT)
topgRNAs$gRNASeq<- factor(topgRNAs$gRNASeq, levels = top_gRNA_levels$gRNASeq)
topgrna_cell_edit_plt <- ggplot(topgRNAs, aes(gRNASeq, n)) +
  geom_bar(stat = "identity",aes(color = HAPLO, fill = HAPLO)) +
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  theme(axis.text.y = element_text(size = 8, colour = 'black')) +
  theme(panel.spacing.x = unit(10, "pt")) +
  theme(strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'top') +
  coord_flip()+
  xlab("")+ylab("")+
  facet_wrap(.~sample,  nrow =1, scale='free_x')
if(i < 5){
ggsave("Top_gRNA_editing_cellNumbers_GATA1_repeats.pdf", topgrna_cell_edit_plt, device = "pdf", height = 8, width = 8)
}else{
  ggsave("Top_gRNA_editing_cellNumbers_HBF1_repeats.pdf", topgrna_cell_edit_plt, device = "pdf", height = 8, width = 8)
  }
#cechkx<- tidymsa_all %>%filter(grepl("NA",freqnew))



topgRNA_full_data<- all_full_data_edit %>% filter(EDTISTATUS =="Y") %>% 
  filter(gRNASeq %in% as.character(top_gRNA_levels$gRNASeq)) %>%
  group_by(gRNASeq,READN) %>%
  group_split()
for(toptim in topgRNA_full_data){
  reftop <- toptim %>% filter(grepl("ref_", namefreq)) %>%
    mutate(namefreq = paste0("000",gRNASeq,":",namefreq)) %>%
    group_by(CELLID) %>% group_split()
  tmptop <- toptim %>% filter(!grepl("ref_", namefreq)) %>%
    mutate(namefreq = paste0(CELLID,":",namefreq)) %>%
    bind_rows(reftop[[1]])
  plttitle= paste0(reftop[[1]]$gRNASeq,":",reftop[[1]]$gRNASeqRC)
  temp_tidymsa<- transform_data_fn ( tmptop, plttitle )
  editplt<- plt_editseqs(temp_tidymsa , plttitle, T)+
    theme(axis.text.y = element_text(size = 6))
  #print(editplt)
}

all_grnas_fullData<- all_full_data_edit %>% #filter(EDTISTATUS =="Y") %>% 
  filter(gRNASeq %in% as.character(Most_prevelent_gRNA$gRNASeq)) %>%
  group_by(gRNASeq,READN) %>%
  group_split()

aiml_data_table<- tibble()
for(ttimx in all_grnas_fullData){
  ttim_grp<- ttimx %>% group_by(CELLID) %>% group_split()
  for(ttim in ttim_grp){
  reftm <- ttim %>% filter(grepl("ref_", namefreq)) %>%
    mutate(namefreq = paste0("000",gRNASeq,":",namefreq)) %>%
    group_by(CELLID) %>% group_split()
  tmpref <- reftm[[1]]
  
  if(grepl(tmpref$gRNASeqRC, tmpref$SEQ)){
    
    tmpref<- tmpref%>% mutate(SEQ=as.character(reverseComplement(DNAString(as.character(SEQ)))))
    refseq = tmpref$SEQ
    tmpt <- ttim %>% filter(!grepl("ref_", namefreq)) %>%
      mutate(SEQ=reveseComp(as.character(SEQ))) %>%
    mutate(namefreq = paste0(CELLID,":",namefreq)) %>%
    bind_rows(tmpref ) %>%
      mutate(mismatch= seq_diff_count(SEQ,refseq), REFSEQ=refseq)
  }else{
    refseq = tmpref$SEQ
    tmpt <- ttim %>% filter(!grepl("ref_", namefreq)) %>%
      mutate(namefreq = paste0(CELLID,":",namefreq)) %>%
      bind_rows(tmpref ) %>%
      mutate(mismatch= seq_diff_count(SEQ,refseq), REFSEQ=refseq)
  }
  
  plttitle= paste0(tmpref$gRNASeq,":",tmpref$gRNASeqRC)
  #temp_tidymsa<- transform_data_fn( tmpt, plttitle )
  #editplt<- plt_editseqs(temp_tidymsa , plttitle, T)+
  #  theme(axis.text.y = element_text(size = 6))
  #print(editplt)
  
  aiml_data_table<- tmpt %>% mutate(PAMSeq= gsub(gRNASeq,"|",SEQ)) %>%
    bind_rows(aiml_data_table)
  }
}

#tmpt %>% mutate(mismatch= pairwiseAlignment(refseq,refseq))
#pwa <- pairwiseAlignment("CCTTTTTTGAATTTCTCACCACAG",refseq,)
#str_diff(tmpt$SEQ, refseq)
if(i<5){
write_csv(aiml_data_table %>% filter((freq > 10 & seqcount >5)), "GATA1_aiml_data_set_from_all_replicate.csv") 
}else{
  write_csv(aiml_data_table %>% filter((freq > 10 & seqcount >5)), "HbF1_aiml_data_set_from_all_replicate.csv") 
}

if(i<5){
  write_csv(aiml_data_table, "GATA1_aiml_complete_data_set_from_all_replicate.csv") 
}else{
  write_csv(aiml_data_table, "HbF1_aiml_complete_data_set_from_all_replicate.csv") 
}


#calculate the editing freq at each base
fre_edit_cal<- aiml_data_table %>% filter((freq > 10 & seqcount >2)) %>%
  group_by(gRNASeq) %>% group_split()
fre_edit_cal[[2]]$gRNASeq

fre_edit_cal_keys<- aiml_data_table %>% filter((freq > 10 & seqcount >2)) %>%
  group_by(gRNASeq) 

ediProfil<- tibble()
all_edit_summary<- tibble()
for(freEdiTemp in fre_edit_cal){
  if('CCAGGGGCCGGCGGCTGGCT' == unique(freEdiTemp$gRNASeq)){
    ofint_grna<- freEdiTemp
  }
  tmpprofX<- freEdiTemp %>% select(gRNASeq,CELLID,SEQ,REFSEQ,namefreq, EDTISTATUS) %>% 
    mutate(X= str_split(SEQ,'',simplify = F)) %>% 
    select(contains('X')) %>% unnest_wider(col = X, simplify = T,names_sep = '') %>%
    
    mutate(across(everything(), ~replace_na(.,'N')))%>% 
    mutate(across(everything(), ~ if_else(. == '-', 'D', .)))%>% 
    pivot_longer(cols = everything(), names_to = "column", values_to = "value") %>%
    count(column, value) %>%
    mutate(location= as.numeric(gsub("X","",column))) %>%
    arrange(location) %>%
    pivot_wider( id_cols = column, names_from = value, values_from = n, values_fill = 0) %>%
    rowwise() %>%
    mutate(across(where(is.numeric), ~ .x / sum(c_across(where(is.numeric))))) %>%
    column_to_rownames(var = 'column')
  
  tmpprofY<- freEdiTemp %>% select(gRNASeq,CELLID,SEQ,REFSEQ,namefreq, EDTISTATUS) %>% 
    mutate(Y= str_split(REFSEQ,'',simplify = F)) %>% 
    select(contains('Y')) %>% unnest_wider(col = Y, simplify = T,names_sep = '') %>%
    
    mutate(across(everything(), ~replace_na(.,'N')))%>% 
    
    pivot_longer(cols = everything(), names_to = "column", values_to = "value") %>%
    count(column, value) %>%
    mutate(location= as.numeric(gsub("Y","",column))) %>%
    arrange(location) %>%
    pivot_wider( id_cols = column, names_from = value, values_from = n, values_fill = 0) %>%
    rowwise() %>%
    mutate(across(where(is.numeric), ~ .x / sum(c_across(where(is.numeric))))) %>%
    mutate(across(where(is.numeric), ~ if_else(. == 1, 0, 1))) %>%
    mutate(D=c(1)) %>%
    column_to_rownames(var = 'column')
  
  ediProfil<- (tmpprofX * tmpprofY[names(tmpprofX)]) %>%
    rownames_to_column(var = "column") %>%
    rowwise() %>%
    mutate(lowestP = sum(c_across(where(is.numeric)), na.rm = TRUE)) %>%
    ungroup() %>%
    select(column,lowestP) %>%
    pivot_wider( names_from = column, values_from = lowestP, values_fill = 0) %>%
    mutate(rowid= unique(freEdiTemp$gRNASeq)) %>%
    bind_rows(ediProfil)
  
  dfx<- freEdiTemp %>% select(gRNASeq,CELLID,SEQ,REFSEQ,namefreq, EDTISTATUS) 
  # --- Step 1: Expand sequences into base-by-base comparison ---
  df_expanded <- dfx %>% 
    select(SEQ,REFSEQ) %>%
    rowwise() %>%
    mutate(
      ref = list(split_chars(REFSEQ)),
      edit = list(split_chars(SEQ))
    ) %>%
    ungroup() %>%
    mutate(row_id = row_number()) %>%
    select(row_id, ref, edit) %>%
    unnest_wider(ref, names_sep = "_") %>%
    unnest_wider(edit, names_sep = "_")
  
  # --- Step 2: Convert wide to long format ---
  n_bases <- nchar(dfx$REFSEQ[1])
  
  base_long <- df_expanded %>%
    pivot_longer(
      cols = -row_id,
      names_to = c("type", "position"),
      names_sep = "_",
      values_to = "base"
    ) %>%
    pivot_wider(names_from = type, values_from = base) %>%
    mutate(position = as.integer(position)) %>%
    select(row_id, position, ref, edit)
  
  # --- Step 3: Identify edits ---
  all_edit_summary <- base_long %>%
    mutate(edit_type = if_else(ref != edit, paste0(ref, ">", edit), NA)) %>%
    filter(!is.na(edit_type)) %>%
    count(position, ref, edit, edit_type, name = "edit_count") %>%
    group_by(position) %>%
    mutate(edit_frequency = edit_count / sum(edit_count)) %>%
    ungroup() %>%
    mutate(overall_frequency = edit_count / dim(dfx)[1],
           gRNAseq=unique(dfx$gRNASeq)) %>%
    bind_rows(all_edit_summary)
  
} 

if(i<5){ 
  
  uniquTopgRNA<- all_edits_stats %>% 
    group_by(sample) %>%
    top_n(n = 50, wt = CELLCOUNT) %>%
    ungroup() %>%
    filter(HAPLO=="UNEDIT") %>%
    select(gRNASeq,HAPLO,CELLCOUNT,sample) %>%
    group_by(gRNASeq) %>%
    mutate(TCOUNT = sum(CELLCOUNT)) %>%
    ungroup() %>%
    distinct(gRNASeq,.keep_all = T) %>%
    arrange(desc(TCOUNT))
  
all_topgRNAs<- all_edits_stats %>% 
  group_by(sample) %>%
  top_n(n = 50, wt = CELLCOUNT) %>%
  ungroup() %>%
  filter(HAPLO !="HETERO") %>%
  select(gRNASeq,HAPLO,EDCOUNT,CELLCOUNT,sample) %>%
  group_by(gRNASeq,HAPLO)%>%
  mutate(TCOUNT=sum(EDCOUNT)) %>%
  ungroup() %>%
  arrange(desc(TCOUNT)) %>%
  mutate(gRNA_HAP= paste0(gRNASeq,"_",HAPLO)) %>%
  distinct(gRNA_HAP,.keep_all = T) %>%
  mutate(EDITSTAT=case_when(HAPLO =="HOMO" ~'EDITED',
                            HAPLO =="UNEDIT" ~ 'UNEDITED')) %>%
  select(gRNASeq,EDITSTAT,TCOUNT) %>%
  mutate(gRNASeq = factor(gRNASeq,levels = as.character(uniquTopgRNA$gRNASeq)))
  
  
}else{
  uniquTopgRNA<- all_edits_stats %>% 
    group_by(sample) %>%
    top_n(n = 75, wt = CELLCOUNT) %>%
    ungroup() %>%
    filter(HAPLO=="UNEDIT") %>%
    select(gRNASeq,HAPLO,CELLCOUNT,sample) %>%
    group_by(gRNASeq) %>%
    mutate(TCOUNT = sum(CELLCOUNT)) %>%
    ungroup() %>%
    distinct(gRNASeq,.keep_all = T) %>%
    arrange(desc(TCOUNT))
  
  all_topgRNAs<- all_edits_stats %>% 
    group_by(sample) %>%
    top_n(n = 75, wt = CELLCOUNT) %>%
    ungroup() %>%
    filter(HAPLO !="HETERO") %>%
    select(gRNASeq,HAPLO,EDCOUNT,CELLCOUNT,sample) %>%
    group_by(gRNASeq,HAPLO)%>%
    mutate(TCOUNT=sum(EDCOUNT)) %>%
    ungroup() %>%
    arrange(desc(TCOUNT)) %>%
    mutate(gRNA_HAP= paste0(gRNASeq,"_",HAPLO)) %>%
    distinct(gRNA_HAP,.keep_all = T) %>%
    mutate(EDITSTAT=case_when(HAPLO =="HOMO" ~'EDITED',
                              HAPLO =="UNEDIT" ~ 'UNEDITED')) %>%
    select(gRNASeq,EDITSTAT,TCOUNT) %>%
    mutate(gRNASeq = factor(gRNASeq,levels = as.character(uniquTopgRNA$gRNASeq)))
  
  
}


data_plot_ediProfil <- ediProfil %>% filter(rowid %in% as.character(uniquTopgRNA$gRNASeq)) %>%
  mutate(across(where(is.numeric),~ if_else(. == 1, 0, .))) %>%
  pivot_longer(cols = contains('X'), names_to = "column") %>%
  mutate(rowid=factor(rowid, levels=as.character(uniquTopgRNA$gRNASeq)),
         column=factor(column, levels=names(ediProfil)),
         location=as.numeric(gsub("X","",column)),
         value = case_when(is.na(value) ~ 0,
                           .default = value),
         frequency= value)
p_tile<- ggplot(data_plot_ediProfil, aes(location,rowid)) +
  geom_tile(aes(fill = frequency),color='gray',stat = 'identity')+
  
  scale_fill_gradientn(colors = c("white", "darkblue"))+
  scale_x_continuous(expand = c(0,0),breaks = seq(1,30))+
  scale_y_discrete(expand = c(0,0))+
  theme_classic()+
  theme(panel.border = element_rect(colour="black",fill = NA,size = 1))+
  theme(axis.line = element_blank())+
  theme(axis.text.y = element_text(size = 10,colour = 'black'))+
  theme(axis.text.x = element_blank(),# element_text(size = 12,colour = 'black'),
        plot.margin = margin(0, 0, 0, 0))+
  theme(legend.position = 'top') +
  xlab("")+ylab("")

p_point<- data_plot_ediProfil %>%  mutate(edit_site = factor(as.character(location),levels = as.character(seq(1,30))) ) %>%
  ggplot( aes(edit_site,value)) +
  #stat_summary(fun = mean, geom = "bar", aes(fill = column))+
  geom_point()+#aes(colour = column))+
  #scale_x_continuous(breaks = seq(1,30))+
  theme_classic()+
  theme(panel.border = element_rect(colour="black",fill = NA,size = 1))+
  theme(axis.text.y = element_text(size = 12,colour = 'black'))+
  theme(axis.line = element_blank())+
  theme(axis.text.x = element_blank(),#element_text(size = 12,colour = 'black'),
        plot.margin = margin(0, 0, 0, 0))+
  theme(legend.position = 'none') +
  xlab("")+ylab("Edit frequency")

p_bar<- all_topgRNAs %>%
  ggplot(aes(gRNASeq,TCOUNT))+geom_bar(aes(fill=EDITSTAT,color=EDITSTAT),stat = 'identity') +
  scale_y_continuous(expand = c(0.0,0.0))+
  scale_fill_manual(values = c('#1b9e77','#F8766D')) +
  scale_color_manual(values = c('#1b9e77','#F8766D')) +
  coord_flip()+
  theme_classic()+
  theme(panel.border = element_rect(colour="black",fill = NA,size = 1))+
  theme(axis.line = element_blank())+
  theme(axis.text.x = element_text(size = 12,colour = 'black'))+
  theme(axis.text.y = element_blank(),
        plot.margin = margin(0, 0, 0, 0))+
  theme(legend.position = 'top') +
  xlab("")+ylab("# of cell")

p_edit_point<- all_edit_summary %>% filter(gRNAseq %in% as.character(uniquTopgRNA$gRNASeq)) %>%
  mutate(edit_site = factor(as.character(position),levels = as.character(seq(1,30))) )%>%
  mutate(across(where(is.numeric),~ if_else(. == 1, 0, .))) %>%
  mutate(Edit_shape = case_when(edit_type == "A>G" ~'A>G',
                                edit_type == "T>C" ~'T>C',
                                grepl(">-",edit_type) ~'Del',
                                .default = c('Sub')) ) %>% 
  mutate(Edit_shape = factor(Edit_shape,levels = c('A>G','T>C','Del','Sub')) )%>%
  ggplot(aes(edit_site,overall_frequency)) +
  geom_point(aes(colour = Edit_shape),size = 1.5 )+
  #scale_x_continuous(breaks = seq(1,30))+
  theme_classic()+
  theme(panel.border = element_rect(colour="black",fill = NA,size = 1))+
  theme(axis.text.y = element_text(size = 12,colour = 'black'))+
  theme(axis.line = element_blank())+
  theme(axis.text.x = element_text(size = 11,colour = 'black'),
        plot.margin = margin(0, 0, 0, 0))+
  theme(legend.position = 'bottom') +
  xlab("")+ylab("Edit frequency")

bottom_row <- p_point + plot_spacer() + plot_layout(widths = c(3, 0.95))
Add_bottom_row <- p_edit_point + plot_spacer() + plot_layout(widths = c(3, 0.95))
bottom_row <- bottom_row/Add_bottom_row+ plot_layout(heights = c(1, 1))
top_row <- p_tile + p_bar + plot_layout(widths = c(3.4, 1))


# Final layout: point plot on top, tile + bar below
final_plot <- top_row / bottom_row + plot_layout(heights = c(5, 1.5))

# Display
final_plot <- final_plot + coord_equal()
final_plot
if(i<5){
  write_csv(all_topgRNAs,"Gata1_topRNA_cellcount_data.csv")
  ggsave("GATA1_topgRNA_location_edit_freq_cellcout_patchwork_plot.pdf",final_plot,device = 'pdf',height = 14,width = 11) 
}else{
  write_csv(all_topgRNAs,"Hbf1_topRNA_cellcount_data.csv")
  ggsave("Hbf1_topgRNA_location_edit_freq_cellcout_patchwork_plot.pdf",final_plot,device = 'pdf',height = 12,width = 10)
}

test_freq<-all_edit_summary %>% mutate(edit_site = factor(as.character(position),levels = as.character(seq(1,30))) )%>%
  mutate(across(where(is.numeric),~ if_else(. == 1, 0, .))) %>%
  mutate(Edit_shape = case_when(edit_type == "A>G" ~'A>G',
                                edit_type == "T>C" ~'T>C',
                                grepl(">-",edit_type) ~'Del',
                                .default = c('Sub')) ) %>% 
  mutate(Edit_shape = factor(Edit_shape,levels = c('A>G','T>C','Del','Sub')) ) %>%
  group_by(edit_site,ref,edit_type) %>%
  reframe(site_freq = edit_count / sum(edit_count)) 
  

top_point_plt<- ediProfil %>%
  mutate(across(where(is.numeric),~ if_else(. == 1, 0, .))) %>%
  pivot_longer(cols = contains('X'), names_to = "column") %>%
  mutate(rowid=factor(rowid, levels=as.character(uniquTopgRNA$gRNASeq)),
         column=factor(column, levels=names(ediProfil)),
         location=as.numeric(gsub("X","",column)),
         value = case_when(is.na(value) ~ 0,
                           .default = value),
         frequency= value) %>%
  mutate(edit_site = factor(as.character(location),levels = as.character(seq(1,30))) ) %>%
  ggplot( aes(edit_site,value)) +
  #stat_summary(fun = mean, geom = "bar", aes(fill = column))+
  geom_point()+#aes(colour = column))+
  scale_y_continuous(expand = c(0,0))+
  theme_classic()+
  theme(panel.border = element_rect(colour="black",fill = NA,size = 1))+
  theme(axis.text.y = element_text(size = 12,colour = 'black'))+
  theme(axis.line = element_blank())+
  theme(axis.text.x = element_blank(),#element_text(size = 12,colour = 'black'),
        plot.margin = margin(0, 0, 0, 0))+
  theme(legend.position = 'none') +
  xlab("")+ylab("Edit frequency")


bottom_edit_point_plt<- all_edit_summary %>% 
  mutate(edit_site = factor(as.character(position),levels = as.character(seq(1,30))) )%>%
  mutate(across(where(is.numeric),~ if_else(. == 1, 0, .))) %>%
  mutate(Edit_shape = case_when(edit_type == "A>G" ~'A>G',
                                edit_type == "T>C" ~'T>C',
                                grepl(">-",edit_type) ~'Del',
                                .default = c('Sub')) ) %>% 
  mutate(Edit_shape = factor(Edit_shape,levels = c('A>G','T>C','Del','Sub')) )%>%
  ggplot(aes(edit_site,overall_frequency)) +
  geom_point(aes(colour = Edit_shape),size = 1.5 )+
  #geom_text(aes(label = ref, x= edit_site),y = 0.0,  vjust = 2,check_overlap = T)+
  scale_y_continuous(expand = c(0,0))+
  coord_cartesian(clip = "off") + 
  theme_classic()+
  theme(panel.border = element_rect(colour="black",fill = NA,size = 1))+
  theme(axis.text.y = element_text(size = 12,colour = 'black'))+
  theme(axis.line = element_blank())+
  theme(axis.text.x = element_text(size = 11,colour = 'black'),
        plot.margin = margin(0, 0, 0, 0))+
  theme(legend.position = 'bottom') +
  xlab("")+ylab("Edit frequency") +
  facet_grid(ref~.)

edit_freq_plt<- top_point_plt / bottom_edit_point_plt + plot_layout(heights = c(1, 4))
edit_freq_plt
if(i<5){
  ggsave("GATA1_allRNA_location_edit_freq_plot.pdf",edit_freq_plt,device = 'pdf',height = 12,width = 8) 
}else{
  ggsave("Hbf1_allgRNA_location_edit_freq_plot.pdf",edit_freq_plt,device = 'pdf',height = 12,width = 8)
}

if(i<5){
plt_grf_abstract<- all_edit_summary %>% 
  mutate(edit_site = factor(as.character(position),levels = as.character(seq(1,30))) )%>%
  mutate(across(where(is.numeric),~ if_else(. == 1, 0, .))) %>%
  mutate(Edit_shape = case_when(edit_type == "A>G" ~'A>G',
                                edit_type == "T>C" ~'T>C',
                                grepl(">-",edit_type) ~'Del',
                                .default = c('Sub')) ) %>% 
  mutate(Edit_shape = factor(Edit_shape,levels = c('A>G','T>C','Del','Sub')) )%>%
  ggplot(aes(edit_site,overall_frequency)) +
  geom_point(aes(colour = Edit_shape),size = 1.5 )+
  geom_text(aes(label = ref, x= edit_site),y = 0.0,  vjust = 1.5,check_overlap = T)+
  scale_y_continuous(expand = c(0,0))+
  scale_x_discrete(position = "top")+
  coord_cartesian(clip = "off") + 
  theme_classic()+
  theme(panel.spacing.y = unit(1.5, "lines"))+
  theme(panel.border = element_rect(colour="black",fill = NA,size = 1))+
  theme(axis.text.y = element_text(size = 12,colour = 'black'))+
  theme(axis.line = element_blank())+
  theme(axis.text.x = element_text(size = 12,colour = 'black'),
        plot.margin = margin(15, 15, 15, 15))+
  theme(
    strip.text = element_blank(),       # Remove text
    strip.background = element_blank()  # Remove strip box
  )+
  theme(legend.position = 'bottom') +
  xlab("")+ylab("Edit frequency") +
  facet_wrap(.~ref, ncol=1)

plt_grf_abstract
ggsave("GATA1_Grap_abstract_freq_plot.pdf",plt_grf_abstract,device = 'pdf',height = 10,width = 10) 
}





blastsgRNA<- list.files(path = "./",pattern = "^mmc*")
blastsgRNA <- blastsgRNA[grepl("mm[2-5]_crophg38_blast.bsl",blastsgRNA)]
all_gRNAs<- tibble()
for(blsfile in blastsgRNA){
  sgRNA_mapping<- read_tsv(blsfile, col_names = F) %>%
    filter(X8==X15) %>%
    mutate(Strand=(X4-X3)/abs(X4-X3)) %>%
    mutate(gRNAID= X1) %>%
    mutate(Chromosome=X2, 
           gRNASeq=X16) %>%
    mutate(`Start Location`=if_else(Strand>0,X3,X4),	`End Location` = if_else(Strand>0,X4,X3)) %>%
    mutate(, Midpoint = `Start Location`+((X3-X4)/2)) %>%
    select(Chromosome, `Start Location`,	`End Location`, gRNAID, Strand, Midpoint, gRNASeq)
  #write.xlsx(as.data.frame(sgRNA_mapping),"All_sgRNA_crophg38_blast.xlsx",sheetName = gsub("_crophg38_blast.bsl","",blsfile), row.names=F, append = T)
  all_gRNAs <- bind_rows(all_gRNAs, sgRNA_mapping)
}
all_gRNAs_csv <- all_gRNAs %>% filter(gRNASeq %in% Most_prevelent_gRNA$gRNASeq)


pltGrps<- tidymsa_all %>% group_by(sample) %>%group_split()


for(pltgr in pltGrps){
  
  for(grplt in pltgr %>% group_by(gRNASeq) %>% group_split()){
    ref_plt<- grplt %>% filter(grepl("ref",freqX)) %>% group_by(gRNASeq,cellid) %>%
      group_split()
    tplt<- grplt %>% filter(!grepl("ref",freqX)) %>%
      bind_rows(ref_plt[[1]])
    editplt<- plt_editseqs(grplt %>% mutate(fas=freqX) , as.character(unique(grplt$gRNASeq)), T)+
    coord_equal()+
    theme(axis.text.y = element_text(size = 8))
    print(editplt)
  }
}




Total_cells<- all_cell_qulify %>% 
  group_by(CELLID,sample) %>% 
  slice_max(seqcount,n=1,with_ties = F) %>%
  ungroup()%>%
  mutate(HAPLO=case_when((EDTISTATUS == 'Y' & freq <= 60 ) ~ "HETERO",
                         (EDTISTATUS == 'Y' & freq > 60) ~ "HOMO",
                         EDTISTATUS == 'N' ~ "UNEDIT")) %>%
  group_by(CELLID) %>%
  ungroup() %>%
  count(HAPLO, sample) %>% 
  ungroup() 
Total_cells<- Total_cells %>% group_by(sample) %>%
 mutate(totalCells=sum(n)) %>%
  mutate(percentage = n/totalCells*100)
  
Plt_cell_ed_per<- Total_cells %>% 
  gather(key = "Keys", value = "value",-HAPLO, -sample, -totalCells, -percentage)
Plt_cell_ed_per$HAPLO<- factor(Plt_cell_ed_per$HAPLO, levels = c('UNEDIT', 'HOMO',  'HETERO')) 
ggplt_cell_per<- ggplot(Plt_cell_ed_per, aes(sample,percentage)) +
  geom_bar(stat = "identity", aes(color=HAPLO, fill=HAPLO), width = 1)+
  coord_polar(theta = "y")+
  geom_text(aes(label = paste0(sprintf("%.1f",percentage),"%")), size = 5, vjust="outward", angle = 45)+
  geom_text(aes(label = sprintf("%d",totalCells)), size = 6,nudge_x = -0.5, nudge_y = 0.5, check_overlap = T)+
  theme_void()+
  facet_wrap(Keys~sample, scale = 'free_x', nrow = 2, ncol=2)
if(i < 5){
ggsave("Editing_per_GATA_all_repeats.pdf",ggplt_cell_per,device = "pdf",width = 10, height = 10)
}else{
  ggsave("Editing_per_HBF_all_repeats.pdf",ggplt_cell_per,device = "pdf",width = 5, height = 5)
}

#######################################

stop("Done excutions")

#######################################

######################################### gRNA speicific editing ###############################################
gRNA_data_grp<- full_data %>% 
  group_by( gRNASeq, READN) %>%
  group_split()
  
  
  
gRNA_edit_data<- tibble()

for(grp in gRNA_data_grp){
  if(nrow(grp)<2){
    next
  }
  tmpref<- grp %>% filter(grepl("Aref_",namefreq)) %>%
    group_by(gRNASeq,READN,SEQ) %>% 
    summarise(seqcount=sum(seqcount)) %>%
    ungroup() %>% mutate(freq = 0) %>%
    mutate(namefreq = paste0('Aref_',READN))
  
  tmpedits<- grp %>% filter(!grepl("Aref_",namefreq))
  if(nrow(tmpedits)==0 ){
    editStat = 'N'
    gRNA_edit_data <- grp %>% filter(!grepl("Aref_",namefreq)) %>%
      group_by(gRNASeq,READN,SEQ) %>% 
      summarise(seqcount=sum(seqcount)) %>%
      mutate(freq=seqcount/sum(seqcount)*100)%>%
      bind_rows(tmpref) %>%
      mutate(EDTISTATUS=c(editStat)) %>% bind_rows(gRNA_edit_data)
    #gRNA_edit_data <- bind_rows(tmpgrpedit,gRNA_edit_data)
    next
  }
  tmpedits<- grp %>% filter(!grepl("Aref_",namefreq)) %>%
    group_by(gRNASeq,READN,SEQ) %>% 
    summarise(seqcount=sum(seqcount)) %>%
    ungroup() %>%
    mutate(freq=seqcount/sum(seqcount)*100)%>%
    filter(SEQ==tmpref$SEQ) %>%
    filter(sum(freq)>80)
    
  editStat = 'N'
  if(nrow(tmpedits)==1 ){
    editStat = 'N'
    gRNA_edit_data <- grp %>% filter(!grepl("Aref_",namefreq)) %>%
      group_by(gRNASeq,READN,SEQ) %>% 
      summarise(seqcount=sum(seqcount)) %>%
      ungroup() %>%
      mutate(freq=seqcount/sum(seqcount)*100)%>%
      arrange(desc(freq))%>%
      filter(freq>1) %>%
      mutate(namefreq = paste0("seq",seq(1:length(freq))," ",seqcount," (",sprintf("%.1f",freq),"%)")) %>%
      ungroup() %>%
      select(names(tmpref)) %>% bind_rows(tmpref) %>%
      mutate(EDTISTATUS=c(editStat)) %>% bind_rows(gRNA_edit_data)
    #gRNA_edit_data <- bind_rows(tmpgrpeditx,gRNA_edit_data)
    next
  }
  
  editStat = 'Y'
  gRNA_edit_data <- grp %>% filter(!grepl("Aref_",namefreq)) %>%
    group_by(gRNASeq,READN,SEQ) %>% 
    summarise(seqcount=sum(seqcount)) %>%
    ungroup() %>%
    mutate(freq=seqcount/sum(seqcount)*100)%>%
    arrange(desc(freq))%>%
    filter(freq>1) %>%
    mutate(namefreq = paste0("seq",seq(1:length(freq))," ",seqcount," (",sprintf("%.1f",freq),"%)")) %>%
    select(names(tmpref)) %>% bind_rows(tmpref) %>%
    mutate(EDTISTATUS=c(editStat)) %>% bind_rows(gRNA_edit_data)
  #gRNA_edit_data <- bind_rows(tmpgrpedit,gRNA_edit_data)
  
}



grans_uni<- gRNA_edit_data %>% filter(!grepl("Aref_",namefreq)) %>% 
  filter(EDTISTATUS == 'Y') %>%
 distinct(gRNASeq, .keep_all = TRUE) 

edited_grna_df<- tibble()
for(grp in gRNA_data_grp){
   
    if(nrow(grp)<2){
      next
    }
   
   tmp<- grp %>% filter(!grepl("Aref_",namefreq)) %>%
    group_by(READN, gRNASeq, SEQ) %>% 
    summarise(seqcount=sum(seqcount)) %>% 
    ungroup()
   if(nrow(tmp)<2){
     next
   }
    tmp<- tmp %>% mutate(freq = seqcount/sum(seqcount)*100) %>%
    arrange(desc(freq))%>%
    filter(freq>1) %>%
    mutate(namefreq = paste0("seq",sprintf("%03d",seq(1:length(freq)))," ",seqcount," (",sprintf("%.1f",freq),"%)")) %>%
    bind_rows(
      grp %>% filter(grepl("Aref_",namefreq)) %>%
        group_by(READN, gRNASeq, SEQ) %>% 
        summarise(seqcount=sum(seqcount)) %>%
        ungroup() %>%
        mutate(freq = 0) %>%
        mutate(namefreq = paste0('Aref_',READN))
    )
    reseq<- tmp %>% filter(grepl("Aref_",namefreq)) %>%
      distinct(SEQ)
    edseq<- tmp %>% filter(SEQ!=as.character(reseq)) %>%
      arrange(desc(freq)) %>% top_n(1,freq) %>%
     filter(sum(freq)>50)
    plttitle<- paste0(unique(grp$gRNASeq),":",unique(grp$gRNASeqRC),":",unique(grp$READN))
    if(nrow(edseq)<1){
      next
    }
    edited_grna_df<- tmp %>% bind_rows(tmp,edited_grna_df)
    temp_tidymsa <- transform_data_fn ( tmp, plttitle )
    editplt<- plt_editseqs(temp_tidymsa, plttitle, T)
    print(editplt)
}
  
  qualified_grna<- edited_grna_df %>% group_by(gRNASeq) %>%
  summarise(sumcount = sum(seqcount)) %>%
  ungroup() %>%
  filter(sumcount>20) 
  edited_grna_ploting<- edited_grna_df %>% 
    filter(gRNASeq %in% qualified_grna$gRNASeq) %>%
    group_by(gRNASeq,READN) %>%
    group_split()
  
  for(pltgrp in edited_grna_ploting){
    plttitle <- paste0(unique(pltgrp$AMPNAME),":",unique(pltgrp$gRNASeq),":",unique(pltgrp$READN))
    temp_tidymsa <- transform_data_fn ( pltgrp, plttitle )
    editplt<- plt_editseqs(temp_tidymsa, plttitle, T)
    print(editplt)
  }
  

##################### Plotting bulk analysis #####################################

bulk_file_path= list(c(  "gatarep1",
                          "/scratch/PRJNA889818/allSampRun001/GATA_rep1/gata1_rep1/result_out/final_count/For_gRNA_AMP_seqcount_freq.csv"))

FColNames1<- c("AMPNAME","gRNASeq","gRNASeqRC", "SEQ", "seqcount", "freq", "namefreq")
for_amp_data<- read.delim(bulk_file_path[[1]][2], sep = ' ', header = F) %>%
  unite("V9", V7, V8, V9, sep = " ", na.rm = TRUE) %>% set_names(FColNames1) %>%
  mutate( READN=c('FOR'))%>%
  select( READN,AMPNAME,gRNASeq,gRNASeqRC,SEQ,seqcount,freq,namefreq)
full_amp_data<- read.delim(gsub("For_", "Rev_", bulk_file_path[[1]][2]), sep = ' ', header = F) %>%
  unite("V9", V7, V8, V9, sep = " ", na.rm = TRUE) %>% set_names(FColNames1) %>%
  mutate( READN=c('REV'))%>%
  select( READN,AMPNAME,gRNASeq,gRNASeqRC,SEQ,seqcount,freq,namefreq) %>% bind_rows(for_amp_data)

full_amp_data <- full_amp_data %>% mutate(AMPNAME = gsub("\\!Edite.*","",AMPNAME))

edited_amp_ploting <- full_amp_data %>% group_by(AMPNAME,READN,gRNASeq) %>% 
  group_split()

for(pltgrp in edited_amp_ploting){
  testplt <- pltgrp %>% mutate(totalsum= sum(seqcount)) %>%
    filter(totalsum >100)
  if(nrow(testplt)<3){
    next
  }
  plttitle <- paste0(unique(pltgrp$gRNASeq),":",unique(pltgrp$READN))
  temp_tidymsa <- transform_data_fn ( pltgrp, plttitle )
  editplt<- plt_editseqs(temp_tidymsa, plttitle, T)
  print(editplt)
}
