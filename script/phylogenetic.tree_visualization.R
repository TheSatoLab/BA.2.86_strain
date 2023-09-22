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

metadata <- fread("merged.final.tsv") #this tsv file was from the BA.2.86_gisaid&nextclade_merged.R script
nwk <- system.file("extdata/SARS-CoV-2", "${name}.nwk", package="treeio")
tree <- read.tree(nwk)
modify.tree <- function(tree) {
  new.tip.labels <- gsub(".*\\|(EPI_ISL_\\d+)\\|.*", "\\1", tree$tip.label)
  tree$tip.label <- new.tip.labels
  return(tree)
}
tree <- modify.tree(tree)
plot(tree)

tip.labels <- as.data.frame(tree$tip.label) %>% set_colnames("tip.labels")
updated.metadata <- metadata %>% inner_join(tip.labels, join_by(Accession.ID == tip.labels))
updated.metadata <- updated.metadata %>% distinct(Accession.ID,.keep_all=T)
updated.metadata$clade_nextstrain <- ifelse(grepl("BA\\.2\\.86", updated.metadata$Pango.lineage) | grepl("BA\\.2\\.86\\.1", updated.metadata$Pango.lineage), "21.L.1", updated.metadata$clade_nextstrain)

updated.metadata <- updated.metadata %>% 
  mutate(group = case_when(
    clade_nextstrain %in% c('19A', '19B', '20A', '20B', '20C', '20D', '20E', '20F', '20G', '21B', '21C', '21D', '21E', '21G', '21H', '21F') ~ 'outgroup',
    clade_nextstrain %in% c('20I') ~ 'Alpha',
    clade_nextstrain %in% c('20H') ~ 'Beta',
    clade_nextstrain %in% c('20J') ~ 'Gamma',
    clade_nextstrain %in% c('21A', '21I', '21J') ~ 'Delta',
    clade_nextstrain %in% c('21K') ~ 'BA.1',
    clade_nextstrain %in% c('22A') ~ 'BA.4',
    clade_nextstrain %in% c('21L','22C') ~ 'BA.2',
    clade_nextstrain %in% c('22B') ~ 'BA.5',
    clade_nextstrain %in% c('22D', '23C') ~ 'BA.2.75',
    clade_nextstrain %in% c('22E') ~ 'BQ.1',
    clade_nextstrain %in% c('22F', '23E') ~ 'XBB.1',
    clade_nextstrain %in% c('23A') ~ 'XBB.1.5',
    clade_nextstrain %in% c('23B') ~ 'XBB.1.16',
    clade_nextstrain %in% c('23D', '23F') ~ 'EG.5.1',
    clade_nextstrain %in% c('21.L.1') ~ 'BA.2.86',
    TRUE ~ NA_character_ 
  ))

groupList <- as.list(NULL)
for(i in 1:nrow(updated.metadata)){
  tip.name <- as.character(updated.metadata[i,]$Accession.ID)
  tip.class <- as.character(updated.metadata[i,]$group)
  if(! tip.class %in% names(groupList)){
    groupList[[tip.class]] <- as.vector(NULL)
  }
  groupList[[tip.class]] <- c(groupList[[tip.class]],tip.name)
}

tree <- groupOTU(tree,groupList)
color.matching <- c(
  "outgroup" = "grey",
  "Alpha" = "#1C67A8",
  "Beta" = "#39A04A",
  "Gamma" = "#D7C046",
  "Delta" = "#E77A25",
  "BA.1" = "#339680",
  "BA.2" = "#02686B",
  "BA.5" = "#803345",
  "BA.2.75" = "#898731",
  "BQ.1" = "#67468C",
  "XBB.1" = "#736357",
  "XBB.1.5" = "#981D2A",
  "XBB.1.16" = "#D95A24",
  "EG.5.1" = "#DE0303",
  "BA.2.86" = "#1D9FF0",
  "BA.4" = "#8C2788"
)

g <- ggtree(tree, layout="equal_angle") + 
  geom_tippoint(aes(color = group)) +
  scale_color_manual(values = color.matching)
g

#done
