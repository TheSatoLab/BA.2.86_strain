#!/usr/bin/env R

#read all the needed files to merge GISAID and Nextclade metadata from the previous scripts
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

write_tsv(x=final.df, file='merged.final.tsv')

#done
