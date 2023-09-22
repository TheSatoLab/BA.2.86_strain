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

write_tsv(x=metadata.filtered, file='metadata.filtered.tsv')
write_tsv(x=metadata.BA.2.86, file='metadata.BA.2.86.tsv')

#done
