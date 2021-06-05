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
  #if g_gn.find("GNAI") != -1 or g_gn.find("GNAO1") != -1 or g_gn.find("GNAT1") != -1 or g_gn.find("Gnai") != -1:
  #  g_lab="Gi/o"
  #else:
  #  g_lab="Gs"
  g_lab=r_gn+"|"+g_gn
  gprotein_lab[pdb1]=g_lab


complex_lab=[]
labs={}
rmsd={}
pdbs=[]
for infile in glob.glob("*fit.txt"):
  pdb1=infile.split("/")[-1].split("_")[0]
  pdb2=infile.split("/")[-1].split("_")[1] 
  if pdb1 not in pdbs:
    pdbs.append(pdb1)
  if pdb2 not in pdbs:
    pdbs.append(pdb2)
  for l in open(infile, "r"):
    if  l.find("Selection RMSD") != -1:
        rms=float(l.split()[3].strip("\n"))
        if pdb1 not in rmsd and pdb2 not in rmsd:
          rmsd[pdb1]={}
          rmsd[pdb1][pdb2]=rms
        elif pdb1 in rmsd and pdb2 not in rmsd[pdb1]:
          rmsd[pdb1][pdb2]=rms
        elif pdb2 in rmsd and pdb1 not in rmsd[pdb2]:
          rmsd[pdb2][pdb1]=rms
print (pdbs)
#print (rmsd)
 

mat=[]
y=0
x=0
x_ref=[]
y_ref=[]
for ii in range(len(pdbs)):
  pdb1=pdbs[ii]
  row=[]
  x=0
  for jj in range(len(pdbs)):
    pdb2=pdbs[jj]
    if pdb1 in rmsd and pdb2 in rmsd[pdb1]:
      row.append(rmsd[pdb1][pdb2])
    elif pdb2 in rmsd and pdb1 in rmsd[pdb2]:
      row.append(rmsd[pdb2][pdb1])
    else:
      row.append(0)

  mat.append(row)
#for r in mat:
#  print (r)
#print mat
mat=np.array(mat, dtype=float)

###Plotting here
rc('font', weight='bold')
fig = plt.figure( figsize=(20,20))
#f, (ax0, ax1) = plt.subplots(2, sharex=True)
ax0 = fig.add_axes([0.075, 0.2,0.8,0.6])


#Reading the alignment for Y axes protein
dend1ax = fig.add_axes([0.175, 0.81,0.6,0.075])
DistMatrix=ssd.pdist(mat)
#X=ssd.squareform(DistMatrix)
#Y = linkage(DistMatrix, method='complete')
Y = linkage(DistMatrix, method='ward')
dend1 = dendrogram(Y, orientation='top')
dend1ax.set_xticks([])
dend1ax.set_yticks([])

#Reading the alignment for Y axes protein
dend2ax = fig.add_axes([0.8, 0.2,0.1,0.6])
DistMatrix=ssd.pdist(mat.T)
#X=ssd.squareform(DistMatrix)
#Y = linkage(DistMatrix, method='complete')
Y = linkage(DistMatrix, method='ward')
dend2 = dendrogram(Y, orientation='right')
dend2ax.set_xticks([])
dend2ax.set_yticks([])


index1 = dend1['leaves']
index2 = dend2['leaves']
mat_r = mat[index1,:]
mat_r2 = mat_r[:,index2]

#complex_lab_reorder=[]
#for i in index1:
#  newlab=complex_lab[i]+"("
#  for lab in labs[complex_lab[i]]:
#    newlab+=lab+":"
#  newlab=(newlab.strip(":"))+")"
#  complex_lab_reorder.append(newlab)
labs_r=[pdbs[i]+"("+gprotein_lab[pdbs[i]]+")" for i in index1]
####Plotting the contact matrix

#cmap = matplotlib.colors.ListedColormap(['white', 'dodgerblue', 'lightgreen', 'tomato'])
#bounds=[0, 1, 2, 3, 4]
#norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
img=ax0.imshow(mat_r2, cmap='Reds', interpolation='nearest', aspect=1, origin='lower')
#img=ax0.imshow(mat, cmap='Reds', interpolation='nearest', aspect=1)

minorLocator = IndexLocator(1,0.5)
majorLocator= IndexLocator(1,0)
minorLocatorY = IndexLocator(1,0)
#print minorLocatorY


ax0.xaxis.set_minor_locator(minorLocator)
ax0.xaxis.set_major_locator(majorLocator)
ax0.yaxis.set_minor_locator(minorLocatorY)
ax0.set_xticklabels(labs_r, rotation=90, minor=True)
major_ticksy = np.arange(0, len(labs_r), 1)
ax0.set_yticks(major_ticksy)
ax0.set_yticklabels(labs_r, rotation=0, ha="right")
ax0.tick_params(axis='x', which='minor', labelsize=17)
ax0.tick_params(axis='y', which='major', labelsize=17)
#ax1.tick_params(axis='y', which='major', labelsize=15)
ax0.xaxis.set_major_formatter(ticker.NullFormatter())
ax0.grid(axis="x", which='major', color="gray", linestyle="-")
ax0.grid(axis="y", which='minor', color="gray", linestyle="-")

axcolor = fig.add_axes([0.95, 0.45,0.01,0.15])
cbar=pylab.colorbar(img, cax=axcolor)#, ticks=[0,1,2,3,4,5,6])
cbar.ax.tick_params(labelsize=15)
#cbar.ax.set_yticklabels(['0', '1', '2', '3', '4', '5', '6']) 

plt.draw()
plt.savefig('rmsd_clusters.eps')
plt.savefig('rmsd_clusters.png')


sys.exit()
