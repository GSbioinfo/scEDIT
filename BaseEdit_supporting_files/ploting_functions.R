library(Biostrings)
# --- Helper: Split string into character vector ---
split_chars <- function(s) unlist(strsplit(s, ""))

reveseComp<- function(x){
  re_revcom=c()
  for (seq in x){
    re_revcom =  c(re_revcom, as.character(reverseComplement(DNAString(seq))))
  }
  return(re_revcom)
  
} 

seq_diff_count <- function(listx,y){
  misatchs <- c()
  for(sti in listx){
    a= str_split_1(sti,"")
    b= str_split_1(y,"")
    msum=0
    for(j in 1:length(a)){
      if(a[j]!=b[j]) msum=msum+1
    }
    misatchs<- c(misatchs, msum)
  }
  return(misatchs)
}

transform_data_fn<- function(rawdat, plttitle){
  rawdat %>% mutate(X= str_split(SEQ,'',simplify = F)) %>% 
  select(namefreq,contains('X')) %>% unnest_wider(col = X, simplify = T,names_sep = '') %>%
  mutate(across(everything(), ~replace_na(.,'N'))) %>%
  gather(key = 'Xloc',value = 'Char',-namefreq)  %>% 
  mutate(value = case_when(Char=='A' ~ 1, 
                           Char == 'T' ~ 2, 
                           Char=='G' ~ 3, 
                           Char == 'C' ~ 4,
                           Char == 'N' ~ 5,
                           !grepl('[A-Za-z]',Char) ~ 0))%>% 
  mutate(value = as.numeric(value))%>%
  group_by(Xloc) %>% 
  mutate(DIFF = ifelse(row_number()==n(), F, ifelse(value - last(value) < 0, F, T))) %>%
  ungroup() %>%
  mutate(colorX = case_when(Char=='A' ~ "#1B9E77", 
                            Char == 'T' ~ "#D95F02", 
                            Char=='G' ~ "#7570B3", 
                            Char == 'C' ~ "#E7298A",
                            Char == 'N' ~ "#eeeeee",
                            !grepl('[A-Za-z]',Char) ~ "white"))%>% 
  mutate(posit = as.numeric(gsub("X","", Xloc)), freqX = factor(namefreq), fas = plttitle) %>%
  mutate(value = case_when(Char=='A' ~ 1, 
                           Char == 'T' ~ 2, 
                           Char=='G' ~ 3, 
                           Char == 'C' ~ 4,
                           Char == 'N' ~ 5,
                           !grepl('[A-Za-z]',Char) ~ 0))%>% 
  mutate(value = as.numeric(value))%>%
  group_by(Xloc) %>% 
  mutate(useDot = ifelse(row_number() == n(), Char, ifelse(value - last(value) != 0, Char, '.'))) %>%
  #mutate(useDot = ifelse(!DIFF, Char, '.')) %>%
  mutate(useDot = ifelse(freqX == "Ref", Char, useDot) ) %>%
  mutate(colordot = case_when(grepl('[A-Za-z]', useDot) ~ "black",
                              useDot == '.' ~ "lightgray",
                              useDot == '-' ~ "black"))%>%
  mutate(color_bg = case_when(grepl('[A-Za-z]', useDot) ~ "black",
                              useDot == '.' ~ "white",
                              useDot == '-' ~ "gray"))%>%
  mutate(color_bg = ifelse(grepl('[A-Za-z]', useDot),colorX,color_bg))%>%
  ungroup() -> trans_data
  return(trans_data)
}
plt_editseqs<- function(temp_tidymsa, pltname, usedot){
  if(usedot){
    tempplt <- ggplot(temp_tidymsa, aes(posit,freqX)) + geom_tile(fill = temp_tidymsa$color_bg, colour = 'black', alpha = 0.35) +
      coord_fixed()  + geom_text(aes(label = useDot),colour = temp_tidymsa$colordot,show.legend = T,fontface = 'bold')+
      scale_x_continuous(expand = c(0,0))+
      scale_y_discrete(position = "right", expand = c(0,0),limit=rev)+
      #scale_color_manual(values = temp_tidymsa$colorX) +
      theme_minimal()+
      theme(axis.text.y = element_text(size = 12, colour = 'black', face = 'bold'))+
      theme(axis.text.x = element_blank())+
      xlab('')+ylab('')+ 
      theme(axis.line = element_blank())+
      theme(panel.background = element_rect(fill = 'white'))+
      #theme(legend.position="none") +
      theme(axis.ticks = element_blank())+
      theme(axis.ticks = element_blank())+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      ggtitle(pltname)
    
  }else{
    tempplt <- ggplot(temp_tidymsa, aes(posit,freqX)) + geom_tile(fill = temp_tidymsa$colorX, colour = 'gray', alpha = 0.55) +
      coord_fixed()  + geom_text(aes(label = Char),colour = 'black',show.legend = T, size = 2,fontface = 'bold')+
      scale_x_continuous(expand = c(0,0))+
      scale_y_discrete(position = "right", expand = c(0,0),limit=rev)+
      #scale_color_manual(values = temp_tidymsa$colorX) +
      theme_minimal()+
      theme(axis.text.y = element_text(size = 12, colour = 'black', face = 'bold'))+
      theme(axis.text.x = element_blank())+
      xlab('')+ylab('')+ 
      theme(axis.line = element_blank())+
      theme(panel.background = element_rect(fill = 'white'))+
      #theme(legend.position="none") +
      theme(axis.ticks = element_blank())+
      theme(axis.ticks = element_blank())+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
      ggtitle(pltname)
  }
  return(tempplt)
}


haplo_plt_do <- function(pltdo, tabelX, Fileabbri){
  genex = paste0(unique(pltdo$GeneName))
  print(genex)
  tabid<-  tabelX[grepl(genex,tabelX$ampId),'tableN'][[1]]
  if(grepl("R1", tabid)){
    winstart= win_amp_merge[grep(genex,win_amp_merge$GeneName), "Win_start"]
    winend= win_amp_merge[grep(genex,win_amp_merge$GeneName), "Win_end"]
    winseq = Win_full_seqs[grep(genex,Win_full_seqs$GeneName), "seq"]
    
  }else{
    winseq = Win_full_seqs[grepl(genex,Win_full_seqs$GeneName), "seqRC"]
    winstart= win_amp_merge[grepl(genex,win_amp_merge$GeneName), "WinR_start"]
    winend= win_amp_merge[grepl(genex,Win_full_seqs$GeneName), "WinR_end"]
  }
  temp_tidymsa<- pltdo %>% select(EDITSEQ, nameFreq) %>% add_row(EDITSEQ = winseq, nameFreq ="ARef") %>%
    mutate(X= str_split(EDITSEQ,'',simplify = F)) %>% 
    select(nameFreq,contains('X')) %>% unnest_wider(col = X, simplify = T,names_sep = '') %>%
    mutate(across(everything(), ~replace_na(.,'N'))) %>%
    gather(key = 'Xloc',value = 'Char',-nameFreq)  %>% 
    mutate(value = case_when(Char=='A' ~ 1, 
                             Char == 'T' ~ 2, 
                             Char=='G' ~ 3, 
                             Char == 'C' ~ 4,
                             Char == 'N' ~ 5,
                             !grepl('[A-Za-z]',Char) ~ 0))%>% 
    mutate(value = as.numeric(value))%>%
    group_by(Xloc) %>% 
    mutate(DIFF = ifelse(row_number()==n(), F, ifelse(value - last(value) < 0, F, T))) %>%
    ungroup() %>%
    mutate(colorX = case_when(Char=='A' ~ "#1B9E77", 
                              Char == 'T' ~ "#D95F02", 
                              Char=='G' ~ "#7570B3", 
                              Char == 'C' ~ "#E7298A",
                              Char == 'N' ~ "#eeeeee",
                              !grepl('[A-Za-z]',Char) ~ "white"))%>% 
    mutate(posit = as.numeric(gsub("X","", Xloc)), freqX = factor(nameFreq), fas = pltname) %>%
    mutate(value = case_when(Char=='A' ~ 1, 
                             Char == 'T' ~ 2, 
                             Char=='G' ~ 3, 
                             Char == 'C' ~ 4,
                             Char == 'N' ~ 5,
                             !grepl('[A-Za-z]',Char) ~ 0))%>% 
    mutate(value = as.numeric(value))%>%
    group_by(Xloc) %>% 
    mutate(useDot = ifelse(row_number() == n(), Char, ifelse(value - last(value) != 0, Char, '.'))) %>%
    #mutate(useDot = ifelse(!DIFF, Char, '.')) %>%
    mutate(useDot = ifelse(freqX == "Ref", Char, useDot) ) %>%
    mutate(colordot = case_when(grepl('[A-Za-z]', useDot) ~ "black",
                                useDot == '.' ~ "lightgray",
                                useDot == '-' ~ "black"))%>%
    mutate(color_bg = case_when(grepl('[A-Za-z]', useDot) ~ "black",
                                useDot == '.' ~ "white",
                                useDot == '-' ~ "gray"))%>%
    mutate(color_bg = ifelse(grepl('[A-Za-z]', useDot),colorX,color_bg))%>%
    ungroup() 
  ycolor= c('black', pltdo$seqCol)
  editplt<- plt_editseqs(temp_tidymsa, genex, T)
  ggsave(filename = paste0(genex,"_",Fileabbri,"_alinment_Petal_plt.pdf"),plot = editplt, width = 8, height = 6)
  print(editplt)
  return(temp_tidymsa)
}

getFullAmps<- function(){
amp_full_seqs<- readDNAStringSet("/scratch/PRJNA889818/allSampRun001/GATA_rep1/gata1_rep1/runinput_dir/AMP_fullseq_GS.fasta")
amp_full_seqs<- amp_full_seqs[grepl("Edite",names(amp_full_seqs)),]
amp_full_seqs <- data.frame(ampnames=names(amp_full_seqs),seq=as.character(amp_full_seqs))
amp_full_seqs$AMPID<- gsub("!.*","",amp_full_seqs$ampnames)
f2<- function(x) as.character(reverseComplement(DNAString(x)))
amp_full_seqs$seqRC <- unlist(lapply(amp_full_seqs$seq, f2))
amp_full_seqs<- amp_full_seqs%>% mutate(seq_len=str_length(seq)) %>% separate(ampnames,into = c("Target","Nothing"), sep = "\\|", remove = F)
return(amp_full_seqs)
}
