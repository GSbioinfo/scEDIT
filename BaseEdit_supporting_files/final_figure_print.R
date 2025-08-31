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

aiml_data_table <- read_csv("GATA1_aiml_complete_data_set_from_all_replicate.csv") %>%
  mutate(gene=paste0(gRNASeq,"_GATA1"))
aiml_data_table <- read_csv("HbF1_aiml_complete_data_set_from_all_replicate.csv") %>%
  mutate(gene=paste0(gRNASeq,"_HbF")) %>%
  bind_rows(aiml_data_table)



#calculate the editing freq at each base
fre_edit_cal<- aiml_data_table %>% filter((freq > 10 & seqcount >2)) %>%
  group_by(gene) %>% group_split()
fre_edit_cal[[2]]$gene

fre_edit_cal_keys<- aiml_data_table %>% filter((freq > 10 & seqcount >2)) %>%
  group_by(gene) 

ediProfil<- tibble()
all_edit_summary<- tibble()
for(freEdiTemp in fre_edit_cal){
  tmpprofX<- freEdiTemp %>% select(gene,CELLID,SEQ,REFSEQ,namefreq, EDTISTATUS) %>% 
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
  
  tmpprofY<- freEdiTemp %>% select(gene,CELLID,SEQ,REFSEQ,namefreq, EDTISTATUS) %>% 
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
    mutate(rowid= unique(freEdiTemp$gene)) %>%
    bind_rows(ediProfil)
  
  dfx<- freEdiTemp %>% select(gene,CELLID,SEQ,REFSEQ,namefreq, EDTISTATUS) 
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
    mutate(edit_type = if_else(ref != edit, paste0(ref, ">", edit), paste0(ref, "=", edit))) %>%
    filter(!is.na(edit_type)) %>%
    count(position, ref, edit, edit_type, name = "edit_count") %>%
    group_by(position) %>%
    mutate(edit_frequency = edit_count / sum(edit_count)) %>%
    ungroup() %>%
    mutate(overall_frequency = edit_count / dim(dfx)[1],
           gene=unique(dfx$gene)) %>%
    bind_rows(all_edit_summary)
  
} 
data_plot_ediProfil <- ediProfil %>% #filter(rowid %in% as.character(uniquTopgRNA$gRNASeq)) %>%
  mutate(across(where(is.numeric),~ if_else(. == 1, 0, .))) %>%
  pivot_longer(cols = contains('X'), names_to = "column") %>%
  mutate(#rowid=factor(rowid, levels=as.character(uniquTopgRNA$gRNASeq)),
         column=factor(column, levels=names(ediProfil)),
         location=as.numeric(gsub("X","",column)),
         value = case_when(is.na(value) ~ 0,
                           .default = value),
         frequency= value)
top_point_plt<- ediProfil %>%
  mutate(across(where(is.numeric),~ if_else(. == 1, 0, .))) %>%
  pivot_longer(cols = contains('X'), names_to = "column") %>%
  mutate(#rowid=factor(rowid, levels=as.character(uniquTopgRNA$gRNASeq)),
         column=factor(column, levels=names(ediProfil)),
         location=as.numeric(gsub("X","",column)),
         value = case_when(is.na(value) ~ 0,
                           .default = value),
         frequency= value) %>%
  mutate(edit_site = factor(as.character(location),levels = as.character(seq(1,30))) ) %>%
  separate(rowid,c('gRNASeq','Gene'), sep = "_",remove = F) %>%
  ggplot( aes(edit_site,value)) +
  #stat_summary(fun = mean, geom = "bar", aes(fill = column))+
  geom_point(aes(shape = Gene))+#aes(colour = column))+
  scale_y_continuous(expand = c(0.01,0.01),limits = c(0,1))+
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
                                  grepl("=",edit_type) ~NA,
                                  .default = c('Sub')) ) %>% 
    mutate(edit_stat = case_when(grepl(">",edit_type) ~T,
                                  grepl("=",edit_type) ~F,
                                  .default = c(F)) ) %>% 
    filter(!is.na(Edit_shape)) %>%
    mutate(Edit_shape = factor(Edit_shape,levels = c('A>G','T>C','Del','Sub')) )%>%
    separate(gene,c('gRNASeq','Gene'), sep = "_",remove = F) %>%
    ggplot(aes(edit_site,overall_frequency)) +
    geom_point(aes(colour = Edit_shape, shape = Gene),size = 1.5 )+
    #geom_text(aes(label = ref, x= edit_site),y = 0.0,  vjust = 1.5,check_overlap = T)+
    scale_y_continuous(expand = c(0.01,0.01),limits = c(0,1))+
    #scale_x_discrete(position = "top")+
    coord_cartesian(clip = "off") + 
    theme_classic()+
    theme(panel.spacing.y = unit(1, "lines"))+
    theme(panel.border = element_rect(colour="black",fill = NA,size = 0.5))+
    theme(axis.text.y = element_text(size = 12,colour = 'black'))+
    theme(axis.line = element_blank())+
    theme(axis.text.x = element_text(size = 12,colour = 'black'),
          plot.margin = margin(15, 15, 15, 15))+
    #theme(
    #  strip.text = element_blank(),       # Remove text
    #  strip.background = element_blank()  # Remove strip box
    #)+
    theme(legend.position = 'bottom') +
    xlab("")+ylab("Edit frequency") +
    facet_grid(ref~.)
  
bottom_edit_point_plt
  edit_freq_plt<- top_point_plt / bottom_edit_point_plt + plot_layout(heights = c(1, 4))
  edit_freq_plt
  ggsave("Final_figure_freq_plot.pdf",edit_freq_plt,device = 'pdf',height = 12,width = 10) 

  edit_cal<- all_edit_summary %>% 
    mutate(edit_site = factor(as.character(position),levels = as.character(seq(1,30))) )%>%
    mutate(across(where(is.numeric),~ if_else(. == 1, 0, .))) %>%
    mutate(Edit_shape = case_when(edit_type == "A>G" ~'A>G',
                                  edit_type == "T>C" ~'T>C',
                                  grepl(">-",edit_type) ~'Del',
                                  grepl("=",edit_type) ~NA,
                                  .default = c('Sub')) ) %>% 
    filter(!is.na(Edit_shape)) %>%
    mutate(Edit_shape = factor(Edit_shape,levels = c('A>G','T>C','Del','Sub')) )%>%
    group_by(edit_site) %>%
    count(ref) 

  length(unique(edit_cal$gRNAseq))



stop("Done excutions")