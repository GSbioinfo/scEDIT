args <- commandArgs(trailingOnly = TRUE)
make_align_plt = FALSE #TRUE
write_missing_BC = FALSE
library(tidyverse)
library(stringr)

library(DBI)
library(RSQLite)
library(here)
library(extrafont)
library(Biostrings)
create_r_packages_file <- function(file_name = "R-packages.R") {
  info <- sessionInfo()
  packages <- c(
    info$otherPkgs, # Attached packages
    info$loadedOnly # Loaded via namespace
  )
  
  lines <- sapply(packages, function(pkg) {
    sprintf('install.packages("%s", version = "%s", repos = "https://cloud.r-project.org/")', 
            pkg$Package, pkg$Version)
  })
  
  writeLines(lines, file_name)
}
sample_num = args[1]
# Generate the R-packages.R file
create_r_packages_file("R-packages.R")
extrafont::loadfonts(quiet = TRUE)
amp_full_seqs<- readDNAStringSet("/scratch/PRJNA889818/allSampRun001/GATA_rep1/gata1_rep1/runinput_dir/AMP_fullseq_GS.fasta")
#amp_full_seqs<- readDNAStringSet("temp_dir/Context_AMP_fullseq_GS.fasta")
amp_full_seqs<- amp_full_seqs[grepl("Edite",names(amp_full_seqs)),]
amp_full_seqs <- data.frame(ampnames=names(amp_full_seqs),seq=as.character(amp_full_seqs))
amp_full_seqs$AMPID<- gsub("!.*","",amp_full_seqs$ampnames)
f2<- function(x) as.character(reverseComplement(DNAString(x)))
amp_full_seqs$seqRC <- unlist(lapply(amp_full_seqs$seq, f2))
amp_full_seqs<- amp_full_seqs%>% mutate(seq_len=str_length(seq)) %>% separate(ampnames,into = c("Target","Nothing"), sep = "\\|", remove = F)

ampgRNA_details<- readxl::read_excel("/scratch/PRJNA889818/gRNA_AMP_details.xlsx") %>%
  mutate(AMPID=gsub("\\$.*","",AMP_gRNA_ID))
csv_file_path= list(c( "/scratch/PRJNA889818/allSampRun001/GATA_rep1/gata1_rep1/result_out/CellFastqs/FastaReads_R1_sqluploads.csv", "gatarep1",
                     "/scratch/PRJNA889818/allSampRun001/GATA_rep1/gata1_rep1/result_out/final_count/sgRNA_BC_GATA1Rep1.csv"),
c("/scratch/PRJNA889818/allSampRun001/GATA_rep2/gata1_rep2/result_out/CellFastqs/FastaReads_R1_sqluploads.csv","gatarep2",
  "/scratch/PRJNA889818/allSampRun001/GATA_rep2/gata1_rep2/result_out/final_count/sgRNA_BC_GATA1Rep2.csv"),
c("/scratch/PRJNA889818/allSampRun001/GATA_rep3/gata1_rep3/result_out/CellFastqs/FastaReads_R1_sqluploads.csv","gatarep3",
  "/scratch/PRJNA889818/allSampRun001/GATA_rep3/gata1_rep3/result_out/final_count/sgRNA_BC_GATA1Rep3.csv"),
c("/scratch/PRJNA889818/allSampRun001/GATA_rep4/gata1_rep4/result_out/CellFastqs/FastaReads_R1_sqluploads.csv","gatarep4",
  "/scratch/PRJNA889818/allSampRun001/GATA_rep4/gata1_rep4/result_out/final_count/sgRNA_BC_GATA1Rep4.csv"),
c("/scratch/PRJNA889818/allSampRun001/HbF_rep1/HbF_rep1/result_out/CellFastqs/FastaReads_R1_sqluploads.csv","hbfrep1",
  "/scratch/PRJNA889818/allSampRun001/HbF_rep1/HbF_rep1/result_out/final_count/sgRNA_BC_HbF1.csv"),
c("/scratch/PRJNA889818/allSampRun001/HbF_rep2/HbF_rep2/result_out/CellFastqs/FastaReads_R1_sqluploads.csv","hbfrep2",
  "/scratch/PRJNA889818/allSampRun001/HbF_rep2/HbF_rep2/result_out/final_count/sgRNA_BC_HbF2.csv"))


mydb <- dbConnect(SQLite(), ":memory:")
for(i in seq(1,6)){
  missingBC<- paste0(dirname(csv_file_path[[i]][1]),"/Missing_gRNA_FastaReads_R2_sqluploads.csv")
  misBCTable<- paste0(csv_file_path[[i]][2],"missBC")
  dbExecute(mydb, paste0("DROP TABLE IF EXISTS ",misBCTable))
  dbExecute(mydb, paste0("CREATE TABLE ",misBCTable," (READID VARCHAR(255)  NOT NULL, 
                       CELLID VARCHAR(255) NOT NULL, AMPNAME VARCHAR(255) NOT NULL,
                       EDIT VARCHAR(255) NOT NULL, SEQUENCE VARCHAR(255) NOT NULL,
                       CIGAR VARCHAR(255) NOT NULL)"))
  read_csv_chunked( paste0(dirname(csv_file_path[[i]][1]),"/Missing_gRNA_FastaReads_R2_sqluploads.csv"), 
                    callback = function(chunk, dummy){
                      dbWriteTable(mydb, misBCTable, chunk, append = T)}, 
                    chunk_size = 10000, col_types = "cccccc", 
                    col_names = c("READID","CELLID","AMPNAME","EDIT", "SEQUENCE", "CIGAR"))
  dbExecute(mydb, paste0("CREATE UNIQUE INDEX ",csv_file_path[[i]][2],"_misBC ON ",misBCTable,"(READID)"))
}

for(i in seq(1,6)){
  misBCTable<- paste0(csv_file_path[[i]][2],"missBC")
  grnamis_bc<- unlist(unique(dbGetQuery(mydb, paste0("SELECT CELLID FROM ", misBCTable))))
  grnamis_bc<- paste(shQuote(grnamis_bc), collapse = ",")
  temp_misBC<- dbGetQuery(mydb, paste0("SELECT CELLID, SEQUENCE FROM ", misBCTable," WHERE CELLID IN(",grnamis_bc,")" )) %>%
    group_by(CELLID, SEQUENCE) %>% count(SEQUENCE) %>% 
    group_by(CELLID) %>%
    slice_max(n, n=1, with_ties = FALSE) %>%
    filter(n>1)
  misgrnaCellIDamp<- data.frame(matrix(ncol = 3,nrow = 0))
  names(misgrnaCellIDamp)<- c("CELLID","AMPNAME","gRNASeq")
  for(inx in 1:nrow(temp_misBC)){
    cid <- as.character(temp_misBC[inx,'CELLID'])
    grn <- as.character(temp_misBC[inx,'SEQUENCE'])
    for(xni in 1:nrow(amp_full_seqs)){
      ampname<- amp_full_seqs[xni,"ampnames"]
      seqamp<- amp_full_seqs[xni,"seq"]
      seqampRC <- amp_full_seqs[xni,"seqRC"]
      #print(seqamp)
      if(grepl(grn,seqamp) || grepl(grn,seqampRC)){
        tdf<- data.frame(cid,ampname,grn)
        names(tdf)<- names(misgrnaCellIDamp)
        misgrnaCellIDamp<- rbind.data.frame(misgrnaCellIDamp,tdf)
      }
    }
  }
  missingBCTable<- paste0(csv_file_path[[i]][2],"BCgSeq")
  dbExecute(mydb, paste0("DROP TABLE IF EXISTS ",missingBCTable))
  dbWriteTable(mydb, missingBCTable, temp_misBC)         
  FoundBCTable<- paste0(csv_file_path[[i]][2],"FoundAMPgSeq")
  dbExecute(mydb, paste0("DROP TABLE IF EXISTS ",FoundBCTable))
  dbWriteTable(mydb, FoundBCTable, misgrnaCellIDamp )
  }

  
 
  


 
for(i in seq(sample_num,sample_num)){
  fastaTable1<- paste0(csv_file_path[[i]][2],"fastaR1")
  dbExecute(mydb, paste0("DROP TABLE IF EXISTS ",fastaTable1))
  dbExecute(mydb, paste0("CREATE TABLE ",fastaTable1," (READID VARCHAR(255)  NOT NULL, 
                       CELLID VARCHAR(255) NOT NULL, AMPNAME VARCHAR(255) NOT NULL,
                       EDIT VARCHAR(255) NOT NULL, SEQUENCE VARCHAR(255) NOT NULL,
                       CIGAR VARCHAR(255) NOT NULL)"))
  print(csv_file_path[[i]][2])
  read_csv_chunked( csv_file_path[[i]][1], 
                    callback = function(chunk, dummy){
                      dbWriteTable(mydb, fastaTable1, chunk, append = T)}, 
                    chunk_size = 10000, col_types = "cccccc", 
                    col_names = c("READID","CELLID","AMPNAME","EDIT", "SEQUENCE", "CIGAR"))
  
  
  #dbExecute(mydb, paste0("UPDATE ",csv_file_path[[i]][2]," SET (READID, CELLID, AMPNAME, EDIT, SEQUENCE, CIGAR) = (SELECT * FROM ",csv_file_path[[i]][2]," ORDER BY CELLID,  AMPNAME)"))
  dbExecute(mydb, paste0("CREATE UNIQUE INDEX ",csv_file_path[[i]][2],"_1Idx ON ",fastaTable1,"(READID)"))
  
  
  
  missingBC<- paste0(dirname(csv_file_path[[i]][1]),"/Missing_gRNA_FastaReads_R2_sqluploads.csv")
  misBCTable<- paste0(csv_file_path[[i]][2],"missBC")
  dbExecute(mydb, paste0("DROP TABLE IF EXISTS ",misBCTable))
  dbExecute(mydb, paste0("CREATE TABLE ",misBCTable," (READID VARCHAR(255)  NOT NULL, 
                       CELLID VARCHAR(255) NOT NULL, AMPNAME VARCHAR(255) NOT NULL,
                       EDIT VARCHAR(255) NOT NULL, SEQUENCE VARCHAR(255) NOT NULL,
                       CIGAR VARCHAR(255) NOT NULL)"))
  read_csv_chunked( paste0(dirname(csv_file_path[[i]][1]),"/Missing_gRNA_FastaReads_R2_sqluploads.csv"), 
                    callback = function(chunk, dummy){
                      dbWriteTable(mydb, misBCTable, chunk, append = T)}, 
                    chunk_size = 10000, col_types = "cccccc", 
                    col_names = c("READID","CELLID","AMPNAME","EDIT", "SEQUENCE", "CIGAR"))
  dbExecute(mydb, paste0("CREATE UNIQUE INDEX ",csv_file_path[[i]][2],"_misBC ON ",misBCTable,"(READID)"))
  csvTable<- paste0(csv_file_path[[i]][2],"table")
  dbExecute(mydb, paste0("DROP TABLE IF EXISTS ",csvTable))
  tabdf<- read_csv(csv_file_path[[i]][3])
  dbWriteTable(mydb, csvTable, tabdf)
  dbExecute(mydb, paste0("CREATE UNIQUE INDEX ",csv_file_path[[i]][2],"_iBC ON ",csvTable,"(CellID)"))
  dbListTables(mydb)
  
#}
#for(i in seq(1,1)){
    grnaTable<-  read_csv(csv_file_path[[i]][3]) %>%
    column_to_rownames('CellID') %>% 
    mutate(rowsum=rowSums(.)) %>% 
    filter(rowsum>0) %>% 
    select(-rowsum) %>%
    rownames_to_column('CELLID') %>%
    pivot_longer(-CELLID) %>%
    #pivot_wider(names_from = CellID, values_from = value) %>% 
    #pivot_longer(cols = everything()) %>%      # Make it long: column name + value
    group_by(CELLID) %>%                         # Group by the original column
    slice_max(value, n = 1, with_ties = FALSE) %>% # Pick top 5 largest values
    filter(value >0 ) %>%
    separate(col = name, c("AMP","gRNAID"), sep = '\\$', fill = 'right') %>%
    mutate(gRNAID=gsub(":.*","", gRNAID))

    knownGRNAs<- grnaTable %>% filter(!is.na(gRNAID))
    missingBCTable<- paste0(csv_file_path[[i]][2],"BCgSeq")
    FoundBCTable<- paste0(csv_file_path[[i]][2],"FoundAMPgSeq")
    misCellBC<- dbGetQuery(mydb, paste0("SELECT * FROM ", missingBCTable ))
    cellBCFound<-  dbGetQuery(mydb, paste0("SELECT * FROM ",  FoundBCTable)) %>%
      mutate(AMPID=gsub("!.*","",AMPNAME))
    temgRNAs<- ampgRNA_details %>% filter(gRNAID %in% knownGRNAs$gRNAID) %>%
      select(AMPID,gRNASeq,gRNAID)
    knownGRNAs <- knownGRNAs %>% inner_join(temgRNAs,by = 'gRNAID', multiple='any')
    allKnowgRNAs<- knownGRNAs %>% select(CELLID,gRNASeq,AMPID)
    allKnowgRNAs <- bind_rows(allKnowgRNAs, cellBCFound %>% select(CELLID,gRNASeq,AMPID)) %>%
      mutate(gRNASeqRC=as.character(reverseComplement(DNAString(gRNASeq))))
    grnaCellIDampSeq<- data.frame(matrix(ncol = 10,nrow = 0))
    names(grnaCellIDampSeq)<- c("CELLID","AMPID","gRNASeq","gRNASeqRC","FSTART","FEND","RSTART","REND", "F_REF", "R_REF")
    for(inx in 1:nrow(allKnowgRNAs)){
      cid <- as.character(allKnowgRNAs[inx,'CELLID'])
      grn <- as.character(allKnowgRNAs[inx,'gRNASeq'])
      grnRC <- as.character(allKnowgRNAs[inx,'gRNASeqRC'])
      for(xni in 1:nrow(amp_full_seqs)){
        ampid<- amp_full_seqs[xni,"AMPID"]
        seqamp<- amp_full_seqs[xni,"seq"]
        seqampRC <- amp_full_seqs[xni,"seqRC"]
        seqamplen <- amp_full_seqs[xni,"seq_len"]
        if(grepl(grn,seqamp) || grepl(grn,seqampRC) ){
          if(grepl(grn,seqamp)){
            stlocF<- regexpr(grn,seqamp)
            enlocF<- stlocF+str_length(grn)-2
            stlocR<- seqamplen-enlocF
            enlocR<- seqamplen-stlocF  
          }
          if(grepl(grn,seqampRC)){
            stlocR<- regexpr(grn,seqampRC)
            enlocR<- stlocR+str_length(grn)-2
            stlocF<- seqamplen-enlocR
            enlocF<- seqamplen-stlocR
          }
          final_stlocF=stlocF-5
          final_enlocF=enlocF+6
          final_stlocR=stlocR-5
          final_enlocR=enlocR+6
          if(final_stlocF<0)
            final_stlocF=0
          if(final_stlocF>seqamplen)
            final_stlocF=seqamplen
          if(final_stlocR < 0)
            final_stlocR = 0
          if(final_enlocR > seqamplen)
            final_enlocR = seqamplen
          tdf<- data.frame(cid,ampid,grn,grnRC,final_stlocF,final_enlocF,final_stlocR,final_enlocR,str_sub(seqamp,final_stlocF,final_enlocF),str_sub(seqampRC,final_stlocR,final_enlocR))
          names(tdf)<- names(grnaCellIDampSeq)
          grnaCellIDampSeq<- rbind.data.frame(grnaCellIDampSeq,tdf)
        }
      }
    }
    grnaCellIDampSeq<- grnaCellIDampSeq %>% mutate(across(where(is.numeric), ~ pmax(., 0)))
    grnaCellIDampSeq$SAMPLE<- csv_file_path[[i]][2]
    
#}
  allForBC_tibble<- tibble()
  
  for(indx in 1:nrow(grnaCellIDampSeq)){
    cbc<- grnaCellIDampSeq[indx,'CELLID']
    ampid<- grnaCellIDampSeq[indx,'AMPID']
    fend<- grnaCellIDampSeq[indx,'FEND']
    rend<- grnaCellIDampSeq[indx,'REND']
    fstart<- grnaCellIDampSeq[indx,'FSTART']
    rstart<- grnaCellIDampSeq[indx,'RSTART']
    grnaSeq<- grnaCellIDampSeq[indx,'gRNASeq']
    grnaSeqRC<- grnaCellIDampSeq[indx,'gRNASeqRC']
    F_ref<- grnaCellIDampSeq[indx,'F_REF']
    R_ref<- grnaCellIDampSeq[indx,'R_REF']
    if(F_ref == ''){
      print(cbc, fstart,  fend, grnaSeq)
      stop()
    }
    amp_bcR1 <- dbGetQuery(mydb, paste0("SELECT * FROM ", fastaTable1," WHERE CELLID = '",cbc,"' AND AMPNAME LIKE '", ampid, "%'" )) %>%
      mutate(SEQID=gsub("\\|.*","",READID)) %>%
      filter(str_length(SEQUENCE)>=fend) %>%
      mutate(FSEQ = str_sub(SEQUENCE, fstart, fend), gRNASeq=c(grnaSeq), gRNASeqRC=c(grnaSeqRC)) %>%
      group_by(AMPNAME,CELLID,gRNASeq,gRNASeqRC) %>% count(FSEQ, name = 'seqcount') %>%
      mutate(freq = seqcount / sum(seqcount)*100) %>% filter(freq > 1) %>%
      arrange(desc(freq)) 
      if(nrow(amp_bcR1)> 0){
        amp_bcR1 <- amp_bcR1 %>% mutate(nameFreq = paste0("seq",sprintf("%03d",seq(1:length(freq)))," ",seqcount," (",sprintf("%.1f",freq),"%)")) %>% 
        ungroup() %>% 
        add_row(AMPNAME=ampid, CELLID=cbc, gRNASeq = grnaSeq, gRNASeqRC=grnaSeqRC,FSEQ=F_ref,seqcount=0,freq = 0, nameFreq = 'Aref_F')
       }else{
         amp_bcR1 <- amp_bcR1 %>% ungroup() %>%
           mutate(nameFreq= NA_character_) %>%
           add_row(AMPNAME=ampid, CELLID=cbc, gRNASeq = grnaSeq, gRNASeqRC=grnaSeqRC,FSEQ=F_ref,seqcount=0,freq = 0, nameFreq = 'Aref_F')
       }
      allForBC_tibble <- bind_rows(allForBC_tibble,amp_bcR1)
      
  }
  write.table(allForBC_tibble, paste0(dirname(csv_file_path[[i]][3]),"/For_cellBC_gRNA_AMP_seqcount_freq.csv"),row.names = F, col.names = F, quote = F, sep = ",")
  dbExecute(mydb, paste0("DROP TABLE IF EXISTS ",fastaTable1))
}

#dbExecute(mydb, paste0("DROP TABLE IF EXISTS ",fastaTable1))
for(i in seq(sample_num,sample_num)){
  fastaTable2<- paste0(csv_file_path[[i]][2],"fastaR2")
  dbExecute(mydb, paste0("DROP TABLE IF EXISTS ",fastaTable2))
  dbExecute(mydb, paste0("CREATE TABLE ",fastaTable2," (READID VARCHAR(255)  NOT NULL, 
                         CELLID VARCHAR(255) NOT NULL, AMPNAME VARCHAR(255) NOT NULL,
                         EDIT VARCHAR(255) NOT NULL, SEQUENCE VARCHAR(255) NOT NULL,
                         CIGAR VARCHAR(255) NOT NULL)"))
  read_csv_chunked( gsub("_R1_sqluploads","_R2_sqluploads",csv_file_path[[i]][1]), 
                    callback = function(chunk, dummy){
                      dbWriteTable(mydb, fastaTable2, chunk, append = T)}, 
                    chunk_size = 10000, col_types = "cccccc", 
                    col_names = c("READID","CELLID","AMPNAME","EDIT", "SEQUENCE", "CIGAR"))
  dbExecute(mydb, paste0("CREATE UNIQUE INDEX ",csv_file_path[[i]][2],"_2Idx ON ",fastaTable2,"(READID)"))
  


  allRevBC_tibble<- tibble()
  for(indx in 1:nrow(grnaCellIDampSeq)){
    cbc<- grnaCellIDampSeq[indx,'CELLID']
    ampid<- grnaCellIDampSeq[indx,'AMPID']
    fend<- grnaCellIDampSeq[indx,'FEND']
    rend<- grnaCellIDampSeq[indx,'REND']
    fstart<- grnaCellIDampSeq[indx,'FSTART']
    rstart<- grnaCellIDampSeq[indx,'RSTART']
    grnaSeq<- grnaCellIDampSeq[indx,'gRNASeq']
    grnaSeqRC<- grnaCellIDampSeq[indx,'gRNASeqRC']
    F_ref<- grnaCellIDampSeq[indx,'F_REF']
    R_ref<- grnaCellIDampSeq[indx,'R_REF']
  amp_bcR2 <- dbGetQuery(mydb, paste0("SELECT * FROM ", fastaTable2," WHERE CELLID = '",cbc,"' AND AMPNAME LIKE '", ampid, "%'" )) %>%
    mutate(SEQID=gsub("\\|.*","",READID)) %>%
    filter(str_length(SEQUENCE) >= rend) %>%
    mutate(RSEQ = str_sub(SEQUENCE, rstart, rend), gRNASeq=c(grnaSeq), gRNASeqRC=c(grnaSeqRC)) %>%
    group_by(AMPNAME,CELLID,gRNASeq,gRNASeqRC) %>% count(RSEQ, name = 'seqcount')  %>%
    mutate(freq = seqcount / sum(seqcount)*100) %>% filter(freq > 1) %>%
    arrange(desc(freq)) 
  if(nrow(amp_bcR2)> 0){
    amp_bcR2 <- amp_bcR2 %>% mutate(nameFreq = paste0("seq",sprintf("%03d",seq(1:length(freq)))," ",seqcount," (",sprintf("%.1f",freq),"%)")) %>% 
      ungroup() %>%
      add_row(AMPNAME=ampid, CELLID=cbc, gRNASeq = grnaSeq, gRNASeqRC=grnaSeqRC,RSEQ=R_ref,seqcount=0,freq = 0, nameFreq = 'Aref_R')
  }else{
    amp_bcR2 <- amp_bcR2 %>% 
      ungroup() %>%
      mutate(nameFreq= NA_character_) %>%
      add_row(AMPNAME=ampid, CELLID=cbc, gRNASeq = grnaSeq, gRNASeqRC=grnaSeqRC,RSEQ=R_ref,seqcount=0,freq = 0, nameFreq = 'Aref_R')
  }
  allRevBC_tibble <- bind_rows(allRevBC_tibble,amp_bcR2)
  }
  write.table(allRevBC_tibble, paste0(dirname(csv_file_path[[i]][3]),"/Rev_cellBC_gRNA_AMP_seqcount_freq.csv"),row.names = F, col.names = F, quote = F, sep = ",")
  dbExecute(mydb, paste0("DROP TABLE IF EXISTS ",fastaTable2))
}
dbDisconnect(mydb)


