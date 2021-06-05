#!/data/SW/anaconda3/envs/myenv/bin/python

import os,sys,operator, math, glob, pickle
import gzip
import numpy as np
from Bio import AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna, generic_protein
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from scipy.stats import hypergeom

fs="\t"

offset={"3sn6.cif.gz":0,
        "6cmo.cif.gz":0,
        "6d9h.cif.gz":0,
        "6g79.cif.gz":3,
        "6gdg.cif.gz":5,
        "6k42.cif.gz":0,
        "6vms.cif.gz":0,
        "6oya.cif.gz":-9
       }


###Script to generate consensus columns from an input MSA 

#print (ac2gn.keys())
###Loading the MSA through Biopython AlignIO function
track={}
aln = AlignIO.read(sys.argv[1], "fasta")
ref_flag=0
    
    
#print ("REFSEQ", ref_id)
#print ("Seq_id\tDeletions\tInsertions\tDivergent_positions")
gaps=[".", "-"]
poscc={}
for a in aln:
  aid=a.id.split(",")[0]
  track[aid]={}
  ###Looping now through the sequences to spot differences
  ###Using three counters, one for the canonical (ref)sequence and one for the target (tar) sequence and one for the MSA
  ref_pos=0
  tar_pos=0 
  ii=0
  jj=0
  #print (aid, a.seq)
  for aa in a.seq:
    if aa not in gaps:
      jj+=1
      track[aid][ii]=jj
      #print (aa,jj)
      if ii not in poscc:
        poscc[ii]=1
      else:
        poscc[ii]+=1
    ii+=1

conpos=[]
for k in poscc.keys():
  if poscc[k] == len(aln):
    conpos.append(k)
  
for aid in track.keys():
  outstr=aid+fs
  for p in conpos:
    outstr+=str(track[aid][p]+offset[aid])+","
    #outstr+=str(p)+","
  print (outstr)
    
    


