#!/usr/bin/env R

library(tidyverse)
library(ape)
library(ggtree)
library(Biostrings)
library(outliers)

args = commandArgs(trailingOnly=T)

pval.threash <- 1.0E-5


data.fasta.name <- args[1]

fastaFile <- readDNAStringSet(data.fasta.name, "fasta")

seq_name.v <- names(fastaFile)
seq.v <- paste(fastaFile)
data.fasta <- data.frame(label = seq_name.v, seq = seq.v)
data.fasta$label <- gsub(" .+","",data.fasta$label)

tree.name <- args[2]
tree <- read.tree(tree.name)
tree.tip.df <- data.frame(tip_Id = 1:length(tree$tip.label), species = tree$tip.label)
tree.info.df <- ggtree(tree)$data
tree.info.df <- tree.info.df[,c('label','branch.length')]

tree.info.df.OTU <- tree.info.df %>% filter(label %in% tree.tip.df$species)
tree.info.df.OTU$branch.length.log <- log(tree.info.df.OTU$branch.length+1,10)

hist(tree.info.df.OTU$branch.length.log)


while(nrow(tree.info.df.OTU)>0){
  max <- max(tree.info.df.OTU$branch.length.log)
  fit <- grubbs.test(tree.info.df.OTU$branch.length.log)
  pval <- fit$p.value
  if(pval > pval.threash){
    break
  }
  tree.info.df.OTU <- tree.info.df.OTU %>% filter(branch.length.log < max)
}


hist(tree.info.df.OTU$branch.length.log)


data.fasta.filtered <- data.fasta %>% filter(label %in% as.character(tree.info.df.OTU$label))



data.fasta.filtered$label <- paste('>',data.fasta.filtered$label,sep="")

out.name <- args[3]
write.table(data.fasta.filtered,out.name,col.names=F,row.names=F,sep="\n",quote=F)
