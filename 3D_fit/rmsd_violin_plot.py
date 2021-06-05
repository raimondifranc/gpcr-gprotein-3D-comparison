#!/data/SW/anaconda3/envs/myenv/bin/python

import sys, os, operator,gzip,glob, math
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits import *
from pylab import *
import pylab
from mpl_toolkits.axes_grid1 import host_subplot
import mpl_toolkits.axisartist as AA
import matplotlib.pyplot as plt
from matplotlib import ticker
#from parsing_tools import GetID
from scipy import stats
import scipy.spatial.distance as ssd
from scipy.cluster.hierarchy import linkage,dendrogram,ward
import seaborn as sns
import pandas as pd


fs="\t"

SS=["N-term","H1","ICL1","H2","ECL1","H3","ICL2","H4","ECL2","H5","ICL3","H6","ECL3","H7","H8","C-term"]

infiles=open("../GPCR-Gprotein_pair_chains.tsv","rt").readlines()
gprotein_lab={}
for ii in range(1, len(infiles)):
  pdb1=infiles[ii].split(fs)[0]
  pdb1_n=pdb1+".cif.gz"
  g_chain1=infiles[ii].split(fs)[8].strip("\n")
  g_uid1=infiles[ii].split(fs)[7]
  r_gn=infiles[ii].split(fs)[2]
  g_gn=infiles[ii].split(fs)[6]
  if g_gn.find("GNAI") != -1 or g_gn.find("GNAO1") != -1 or g_gn.find("GNAT1") != -1 or g_gn.find("Gnai") != -1:
    g_lab="Gi/o"
  else:
    g_lab="Gs"
  #g_lab=r_gn+"|"+g_gn
  gprotein_lab[pdb1]=g_lab


complex_lab=[]
labs={}
rmsd={}
pdbs=[]
stat={"Gs":[],"Gi/o":[]}
for infile in glob.glob("*fit.txt"):
  pdb1=infile.split("/")[-1].split("_")[0]
  pdb2=infile.split("/")[-1].split("_")[1]
  glab1=gprotein_lab[pdb1]
  glab2=gprotein_lab[pdb2] 
  if glab1 == glab2:
    for l in open(infile, "r"):
      if  l.find("Selection RMSD") != -1:
        rms=float(l.split()[3].strip("\n"))
        stat[glab1].append(rms)
 

#print (stat)

#d = dict( oncogene = np.array(role_stat1["oncogene"]), oncogene_TSG = np.array(role_stat1["oncogene_TSG"]), TSG = np.array(role_stat1["TSG"]), fusion = np.array(role_stat1["fusion"]), oncogene_TSG_fusion = np.array(role_stat1["oncogene_TSG_fusion"]), oncogene_fusion = np.array(role_stat1["oncogene_fusion"]))


df = pd.DataFrame(dict(Gs=pd.Series(stat['Gs']),\
                       Gio=pd.Series(stat['Gi/o'])))

print (df)
sns.set_style('whitegrid')
"#ID_SAMPLE\tSAMPLE_NAME\tGENE_NAME\tREGULATION\tZ_SCORE\tID_STUDY"



ax = sns.violinplot(data=[stat["Gs"], stat["Gi/o"]])
ax.set_xticklabels(['Gs', 'Gi/o'])

print (np.median(stat["Gs"]), np.std(stat["Gs"]), stats.sem(stat["Gs"]))
print (np.median(stat["Gi/o"]), np.std(stat["Gi/o"]), stats.sem(stat["Gi/o"]))
print (stats.ranksums(stat["Gs"], stat["Gi/o"]))

x1, x2 = 0, 1   # columns 'Sat' and 'Sun' (first column: 0, see plt.xticks())
y, h, col = max(np.array(stat["Gi/o"])) + 4, 0.5, 'gray'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h+0.15, "P-val=5e-17", ha='center', va='bottom', color='black')

ax.set_xticklabels(["Gs", "Gi/o"], size = 15)

#sns.axlabel('', 'Nr. of pathways', fontsize=15)
#plt.title('Cosmic Census Oncodrivers', fontsize=15)
plt.ylim([-4,27])

plt.savefig('Gs_vs_Gio_rmsd.png')


sys.exit()
