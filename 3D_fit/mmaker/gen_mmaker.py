#!/data/SW/anaconda3/envs/myenv/bin/python

import os,sys
import pandas as pd

fs="\t"

Gprot_r={}
Gpcr_r={}
for l in open("Gprot_pdbs_consensus.txt", "rt"):
  pdb=l.split(fs)[0]
  res=l.split(fs)[1].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gprot_r[pdb]=res

for l in open("GPCR_pdbs_consensus.txt", "rt"):
  pdb=l.split(fs)[0]
  res=l.split(fs)[1].strip("\n").strip(",")
  #print (pdb, len(res.split(",")))
  Gpcr_r[pdb]=res
  

infiles=open("../GPCR-Gprotein_pair_chains.tsv","rt").readlines()
for ii in range(1, len(infiles)):
  for jj in range(ii+1, len(infiles)):
    pdb1=infiles[ii].split(fs)[0]
    pdb2=infiles[jj].split(fs)[0]
    pdb1_n=pdb1+".cif.gz"
    pdb2_n=pdb2+".cif.gz"
    outfile=open("%s_%s_ct_mmaker_nogui" % (pdb1, pdb2), "wt")
    outfile.write("open ../%s\n" % (pdb1_n))
    outfile.write("open ../%s\n" % (pdb2_n))

    r_chain1=infiles[ii].split(fs)[4]
    r_chain2=infiles[jj].split(fs)[4]
    g_chain1=infiles[ii].split(fs)[8].strip("\n")
    g_chain2=infiles[jj].split(fs)[8].strip("\n")
    g_uid1=infiles[ii].split(fs)[7]
    g_uid2=infiles[jj].split(fs)[7]
    r_pc1=pdb1+"_"+r_chain1
    r_pc2=pdb2+"_"+r_chain2
    g_pc1=pdb1+"_"+g_chain1
    g_pc2=pdb2+"_"+g_chain2
    
    r_sele1=Gpcr_r[r_pc1]
    r_sele2=Gpcr_r[r_pc2]
    g_sele1=Gprot_r[g_pc1]
    g_sele2=Gprot_r[g_pc2]
    
    
###First fitting the receptor chains
    #outfile.write("mmaker #0:%s.%s #1:%s.%s pair ss\n" % (r_sele1, r_chain1, r_sele2, r_chain2))
    outfile.write("mmaker #0:.%s #1:.%s pair ss\n" % (r_chain1, r_chain2))

###Then calculating the RMSD on G-protein selections
    #outfile.write("rmsd #0:360-384.D@CA #1:370-394.A@CA\n")
    print (len(g_sele1.split(",")), len(g_sele2.split(",")))
    outfile.write("rmsd #0:%s.%s@CA #1:%s.%s@CA\n" % (g_sele1, g_chain1, g_sele2, g_chain2))
