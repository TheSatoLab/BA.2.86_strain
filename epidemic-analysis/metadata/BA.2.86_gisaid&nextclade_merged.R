#!/usr/bin/env R

library(stringr)
library(readr)
library(ape)
library(dplyr)
library(data.table)
library(magrittr)

#read all the needed files to merge GISAID and Nextclade metadata from the previous scripts
metadata.filtered <- fread('metadata.tsv')
nextclade <- fread('nextclade.cut.tsv')
nextclade$seqName <- str_split(nextclade$seqName, "\\|", simplify=T)[,1]

#keep-both-names
merged.df <- nextclade.cut %>% inner_join(metadata.filtered, join_by(seqName==Virus.name), relationship = "many-to-many")
merged.df <- merged.df %>% distinct(seqName,.keep_all=T)

final.metadata <- merged.df %>% distinct(seqName,.keep_all=T)
final.metadata <- final.metadata %>% filter(Collection.date >= date.start)

final.metadata <- final.metadata %>%
  mutate(Collection.date = as.Date(Collection.date),
         region = str_split(Location,"/",simplify = T)[,1],
         country = str_split(Location,"/",simplify = T)[,2],
         state = str_split(Location,"/",simplify = T)[,3])

final.metadata$clade_nextstrain <- ifelse(grepl("BA\\.2\\.86", final.metadata$Pango.lineage) | grepl("BA\\.2\\.86\\.1", final.metadata$Pango.lineage), "BA.2.86", final.metadata$clade_nextstrain) #rename-BA.2.86 (and BA.2.86.1) as BA.2.86
final.metadata.filtered <- final.metadata %>% 
  mutate(group = case_when(
    clade_nextstrain %in% c('19A', '19B', '20A', '20B', '20C', '20D', '20E', '20F', '20G', '21B', '21C', '21D', '21E', '21G', '21H', '21F') ~ 'outgroup',
    clade_nextstrain %in% c('20I') ~ 'Alpha',
    clade_nextstrain %in% c('20H') ~ 'Beta',
    clade_nextstrain %in% c('20J') ~ 'Gamma',
    clade_nextstrain %in% c('21A', '21I', '21J') ~ 'Delta',
    clade_nextstrain %in% c('21K') ~ 'BA.1',
    clade_nextstrain %in% c('22A') ~ 'BA.4',
    clade_nextstrain %in% c('21L') ~ 'BA.2',
    clade_nextstrain %in% c('22B') ~ 'BA.5',
    clade_nextstrain %in% c('22D') ~ 'BA.2.75',
    clade_nextstrain %in% c('22E') ~ 'BQ.1',
    clade_nextstrain %in% c('22F') ~ 'XBB.1',
    clade_nextstrain %in% c('23A') ~ 'XBB.1.5',
    clade_nextstrain %in% c('23B') ~ 'XBB.1.16',
    clade_nextstrain %in% c('23F') ~ 'EG.5.1',
    clade_nextstrain %in% c('BA.2.86') ~ 'BA.2.86',
    TRUE ~ NA_character_ 
  )) %>% filter(!is.na(group))

#write-final.metadata.filtered
write_tsv(x=final.metadata.filtered, file='final.metadata.filtered.tsv')

#done
