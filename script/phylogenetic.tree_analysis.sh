#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 18
#$ -l exclusive,s_vmem=10G,mem_req=10G
##$ -o /dev/null
##$ -e /dev/null
##$ -t 1-2:1

#align using minimap2
minimap2 -a -x asm20 --sam-hit-only --secondary=no --score-N=0 -t 10 \
NC_045512.fasta \
${name}.fasta \
-o ${name}.sam

#convert sam file to fasta file
gofasta-linux-amd64 sam toMultiAlign \
-s ${name}.sam \
-t 10 \
--reference NC_045512.fasta \
--trimstart 265 \
--trimend 29674 \
--trim \
--pad > ${name}.aligned.fasta

#trim uncertain sites
trimal \
-in ${name}.aligned.fasta \
-out ${name}.aligned.trimmed.fasta \
-htmlout /dev/null \
-gt 0.1

##three-step-protocol
#construction of the first tree
iqtree2 \
  -s ${name}.aligned.trimmed.fasta \
  -nt AUTO -mem 8G -m GTR -bb 1000

#tips with longer external branch to be removed
R --vanilla --no-echo --args \
${name}.aligned.trimmed.fasta \
${name}.treefile \
${name}.aligned.trimmed.2nd.fasta  \
< ~/cov-phylo/omicron_tree/BA2/script/detect_outlier_OTU.R

#construction of the second tree
iqtree2 \
-s ${name}.aligned.trimmed.2nd.fasta \
-nt AUTO -mem 8G -m GTR -bb 1000

#done
#convert the tree file into a newick file for phylogenetic.tree_visualization.R script
