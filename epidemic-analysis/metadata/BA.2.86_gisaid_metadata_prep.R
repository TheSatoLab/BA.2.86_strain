#!/usr/bin/env R

library(stringr)
library(readr)
library(data.table)
library(dplyr)
library(stringr)
library(magrittr)

#Extraction of BA.2.86 and BA.2.86.1 data from GISAID metadata
metadata <- fread('metadata.tsv')
colnames(metadata) <- gsub("[ _-]", ".", colnames(metadata))
metadata.filtered <- metadata %>% distinct(Accession.ID,.keep_all=T) %>%
  filter(Host == "Human",
         #!N.Content > 0.01 | is.na(N.Content),
         str_length(Collection.date) == 10,
         Pango.lineage != "",
         Pango.lineage != "None",
         Pango.lineage != "Unassigned",
         !str_detect (Additional.location.information,"[Qq]uarantine")
  )

write_tsv(x=metadata.filtered, file='metadata.filtered.tsv')

#done
