#!/bin/sh
#$ -S /bin/bash
#$ -pe def_slot 18
#$ -l exclusive,s_vmem=10G,mem_req=10G
##$ -o /dev/null
##$ -e /dev/null
##$ -t 1-2:1

cut -f 1-8 nextclade.tsv > nextclade.cut.tsv
