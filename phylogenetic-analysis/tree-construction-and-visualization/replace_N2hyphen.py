#!/usr/bin/env python

import sys, re
argvs = sys.argv

fasta_f = open(argvs[1])

seq_d = {}
for line in fasta_f:
  line = line.strip()
  if '>' in line:
    name = line.replace('>','')
    seq_d[name] = ""
  else:
    seq = line
    seq_d[name] += seq

for name in seq_d:
  seq = seq_d[name]
  seq = seq.replace('N','-')
  res = ">" + name + "\n" + seq
  print(res)
