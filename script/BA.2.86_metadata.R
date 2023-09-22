#!/usr/bin/env R

library(stringr)
library(readr)
library(ggtree)
library(treeio)
library(ggplot2)
library(ape)
library(dplyr)
library(tidyr)
library(tidytree)
library(openxlsx)
library(data.table)
library(magrittr)

#Extraction of BA.2.86 and BA.2.86.1 data from GISAID metadata
metadata <- fread('metadata.tsv')
colnames(metadata) <- gsub("[ _-]", ".", colnames(metadata))
metadata.BA.2.86 <- metadata %>% filter(Pango.lineage %in% c("BA.2.86", "BA.2.86.1"))
BA.2.86_name_list <- metadata.BA.2.86$Virus.name
metadata.filtered <- metadata %>% distinct(Accession.ID,.keep_all=T) %>%
  filter(Host == "Human",
         !N.Content > 0.01 | is.na(N.Content),
         str_length(Collection.date) == 10,
         Pango.lineage != "",
         Pango.lineage != "None",
         Pango.lineage != "Unassigned",
         !str_detect (Additional.location.information,"[Qq]uarantine")
  )
metadata.filtered <- metadata.filtered %>% filter(!(Virus.name %in% BA.2.86_name_list))

write_tsv(x=metadata.filtered, file='$metadata.filtered.tsv')
write_tsv(x=metadata.BA.2.86, file='$metadata.BA.2.86.tsv')

##COMMENTOUT
#cut the nextclade.tsv, maintaing only rows 1-8
#See BA.2.86_nextclade_metadata_prep.txt

#read all the needed files to merge GISAID and Nextclade metadata
metadata.filtered <- fread('metadata.tsv')
metadata.BA.2.86 <- fread('metadata.BA.2.86.tsv')
nextclade <- fread('nextclade.cut.tsv')
nextclade$seqName <- str_split(nextclade$seqName, "\\|", simplify=T)[,1]

#keep-both-names
combo.df <- nextclade.cut %>% inner_join(metadata.filtered, join_by(seqName==Virus.name), relationship = "many-to-many")
combo.df.1 <- nextclade.cut %>% inner_join(metadata.BA.2.86, join_by(seqName==Virus.name))
combo.df <- combo.df %>% distinct(seqName,.keep_all=T)
combo.df.1 <- combo.df.1 %>% distinct(seqName,.keep_all=T)

#sampled 20 per clade based on clade_nextstrain
num.samples <- 20
combo.df.sampled <- combo.df %>%
  filter(
    !is.na(clade_nextstrain), 
    clade_nextstrain != "", 
    clade_nextstrain != "recombinant", #removing any recombinant group
    clade_nextstrain != "21M" #removing 21M (Omicron, B.1.1.529)
  ) %>%
  group_by(clade_nextstrain) %>%
  slice_sample(n = num.samples, replace = TRUE) %>%
  ungroup()

combo.df.1.filtered <- combo.df.1 %>% distinct(Accession.ID,.keep_all=T) %>%
  filter(Host == "Human",
         !N.Content > 0.03 | is.na(N.Content),
         str_length(Collection.date) == 10,
         Pango.lineage != "",
         Pango.lineage != "None",
         Pango.lineage != "Unassigned",
         #!str_detect (Additional.location.information,"[Qq]uarantine")
  )

final.df <- rbind(combo.df.sampled, combo.df.1.filtered)
final.df <- final.df %>% distinct(seqName,.keep_all=T)

write_tsv(x=final.df, file='$merged.final.tsv')

