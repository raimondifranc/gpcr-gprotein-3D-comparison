#!/data/SW/anaconda3/envs/myenv/bin/python

import os,sys
import pandas as pd

fs="\t"

Gprot_r={}
Gpcr_r={}
for l in open("../Gprot_pdbs_consensus.txt", "rt"):
  pdb=l.split(fs)[0].split("_")[0]
  chain=l.split(fs)[0].split("_")[1]
  res=l.split(fs)[1].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gprot_r[pdb]=[chain,res]

for l in open("../GPCR_pdbs_consensus.txt", "rt"):
  pdb=l.split(fs)[0].split("_")[0]
  chain=l.split(fs)[0].split("_")[1]
  res=l.split(fs)[1].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gpcr_r[pdb]=[chain,res]
  
#print (list(Gpcr_r.keys()))
print ("#!/bin/bash")
for jj in range(len(Gpcr_r.keys())):
  pdb1="../3sn6.cif"
  rc1=Gpcr_r[pdb1][0]
  r_res1=Gpcr_r[pdb1][1]

  pdb2=list(Gpcr_r.keys())[jj]
  rc2=Gpcr_r[pdb2][0]
  r_res2=Gpcr_r[pdb2][1]

  mc1=Gprot_r[pdb1][0]
  m_res1=Gprot_r[pdb1][1]
  mc2=Gprot_r[pdb2][0]
  m_res2=Gprot_r[pdb2][1]

  print ("/data/Users/francesco/code/bio_align_fit_n_rmsd.py -refe ../%s -mobi ../%s -r_fit_chain %s -m_fit_chain %s -r_fit_sele %s -m_fit_sele %s -r_rms_chain %s -m_rms_chain %s -r_rms_sele %s -m_rms_sele %s > %s_%s_fit.txt" % (pdb1, pdb2, rc1, rc2, r_res1, r_res2, mc1, mc2, m_res1, m_res2, pdb1.split("/")[1].split(".cif")[0], pdb2.split("/")[1].split(".cif")[0]))
