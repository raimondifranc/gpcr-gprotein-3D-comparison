#!/data/SW/anaconda3/envs/myenv/bin/python

import sys, os, operator,gzip,glob, math, pickle
sys.path.insert(0, "./")
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
from parsing_tools import GetID

fs="\t"

AA_color={'A':'sienna','F':'sandybrown','G':'peru','I':'bisque','L':'tan','M':'olive','P':'y','V':'yellow','Y':'greenyellow','W':'lime','D':'red','E':'crimson','K':'b','R':'cyan','Q':'violet','N':'orchid','S':'fuchsia','T':'salmon','H':'skyblue','C':'pink'}
pfam_aa="gNllvllviltkkslrtstnyfilsLaisDlllgllvlpfaiiyeltnewlfgsvlCklwlsldvvlstasilnllaisiDRYvaivkplkykrlkrrakavillvWvlslllavpllvletvnaqdtnptksClittllpsvyiifsavlsFilPllvivilyarifrtvrktsrtmtvkskngslakgsetkkkvkrrslrkekkaartllivvvvFivcwlPftiiyllqalskqtvaaeniptslftiavwlgysnsavNPiiY"
pfam_ss="HHHHHHHHHHH-CCC-SHHHHHHHHHHHHHHHHHHHCHHHHHHHHHHSSCHSHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHSHHHHHCCHHHHHHHHHHHHHHHHHHHHHHHHTCCEHHCCCTTGGGEESSSTCCHHHHHHHHHHHCHHHHHHHHHHHHHHHHCHHCCC............XXXXXXXXXXX-CCCHCCHHHHHHHHHHHHHHHHHHHHHHHHHHHHCHCTTSTHHHHHHCHHHHHHHHHHHHHHHHHHHHH"
BW=["1.49","1.50","1.51","1.52","1.53","1.54","1.55","1.56","1.57","1.58","1.59","ICL1(12)","ICL1(13)","ICL1(14)","ICL1(15)","ICL1(16)","ICL1(17)","2.38","2.39","2.40","2.41","2.42","2.43","2.44","2.45","2.46","2.47","2.48","2.49","2.50","2.51","2.52","2.53","2.54","2.55","2.56","2.57","2.58","2.59","2.60","2.61","2.62","2.63","2.64","2.65","2.65","2.66","2.67","ECL1(49)","ECL1(50)","ECL1(51)","ECL1(52)","3.21","3.22","3.23","3.24","3.25","3.26","3.27","3.28","3.29","3.30","3.31","3.32","3.33","3.34","3.35","3.36","3.37","3.38","3.39","3.40","3.41","3.42","3.43","3.44","3.45","3.46","3.47","3.48","3.49","3.50","3.51","3.52","3.53","3.54","3.55","3.56","3.57","3.58","3.59","3.60","3.61","ICL2(94)","ICL2(95)","4.39","4.40","4.41","4.42","4.43","4.44","4.45","4.46","4.47","4.48","4.49","4.50","4.51","4.52","4.53","4.54","4.55","4.56","4.57","4.58","4.59","4.60","4.61","4.62","ECL2(120)","ECL2(121)","ECL2(122)","ECL2(123)","ECL2(124)","ECL2(125)","ECL2(126)","ECL2(127)","ECL2(128)","ECL2(129)","ECL2(130)","ECL2(131)","ECL2(132)","ECL2(133)","ECL2(134)","ECL2(135)","ECL2(136)","ECL2(137)","ECL2(138)","ECL2(139)","ECL2(140)","ECL2(141)","5.36","5.37","5.38","5.39","5.40","5.41","5.42","5.43","5.44","5.45","5.46","5.47","5.48","5.49","5.50","5.51","5.52","5.53","5.54","5.55","5.56","5.57","5.58","5.59","5.60","5.61","5.62","5.63","5.64","5.65","5.66","5.67","ICL3(174)","ICL3(175)","ICL3(176)","ICL3(177)","ICL3(178)","ICL3(179)","ICL3(180)","ICL3(181)","ICL3(182)","ICL3(183)","ICL3(184)","ICL3(185)","ICL3(186)","ICL3(187)","ICL3(188)","ICL3(189)","ICL3(190)","ICL3(191)","ICL3(192)","ICL3(193)","ICL3(194)","ICL3(195)","ICL3(196)","ICL3(197)","ICL3(198)","ICL3(199)","6.25","6.26","6.27","6.28","6.29","6.30","6.31","6.32","6.33","6.34","6.35","6.36","6.37","6.38","6.39","6.40","6.41","6.42","6.43","6.44","6.45","6.46","6.47","6.48","6.49","6.50","6.51","6.52","6.53","6.54","6.55","6.56","6.57","6.58","6.59","6.60","ECL3(236)","ECL3(237)","ECL3(238)","ECL3(239)","ECL3(240)","7.26","7.27","7.28","7.29","7.30","7.31","7.32","7.33","7.34","7.35","7.36","7.37","7.38","7.39","7.40","7.41","7.42","7.43","7.44","7.45","7.46","7.47","7.48","7.49","7.50","7.51","7.52","7.53"]

SS=["N-term","H1","ICL1","H2","ECL1","H3","ICL2","H4","ECL2","H5","ICL3","H6","ECL3","H7","H8","C-term"]

def GetGpcrdb(GPCR,GPCRDB):
  uac=GPCR[0]
  p=GPCR[1]
  if p in GPCRDB[uac]:
    lab=GPCRDB[uac][p]
  else:
    lab="-"  
  return lab  

def GetCgn(GPROT,CGN):
  uac=GPROT[0]
  p=GPROT[1]
  if p in CGN[uac]:
    lab=CGN[uac][p]
  else:
    lab="-"
  return lab

unipid_uac={}

#with open('/data/DB/uniprot/HUMAN_9606_ID2AC.pickle', 'rb') as handle:
#  unipid_uac = pickle.load(handle)

unipac_gn={}
gn_unipid={}
gn_unipac={}
for line in gzip.open("./Uniprot_Genenames.tsv.gz", "rt"):
  unipid=line.split("\t")[1]
  unipac=line.split("\t")[0]
  if len(line.split("\t")) > 2:# and unipid.find("_HUMAN") != -1:
    gn=line.split("\t")[2]
    gn=gn.strip("\n")
    if unipid not in unipid_uac:
      unipid_uac[unipid]=[unipac]
    elif unipac not in unipid_uac[unipid]:
      unipid_uac[unipid].append(unipac)
    unipac_gn[unipac]=gn
    if gn not in gn_unipid:
      gn_unipid[gn]=[]
      gn_unipid[gn].append(unipid)
    elif gn in gn_unipid:
      gn_unipid[gn].append(unipid)
    if gn not in gn_unipac:
      gn_unipac[gn]=[]
      gn_unipac[gn].append(unipac)
    elif gn in gn_unipac:
      gn_unipac[gn].append(unipac)


gpcrdb={}
for infile in glob.glob("*gpcrdb.tsv.txt"):
  uac=infile.split("/")[-1].split("_")[1]
  gpcrdb[uac]={}
  for l in open(infile, "rt"):
    num=l.split(fs)[0]
    pos=l.split(fs)[1][1:].strip("\n").strip("\r")
    gpcrdb[uac][pos]=num


cgns={}
gprot=[]
for l in open("CGN_human_sequence_table.txt","rt"):
  if l.split(fs)[0] == "Sort number":
    for g in l.split(fs)[1:]:
      g=g.strip("\n")
      gprot.append(g)
      for gg in g.split("/"):
        if gg in unipid_uac:
          for uac in unipid_uac[gg]:
            cgns[uac]={}
  else:
    cgn=l.split(fs)[0]
    cgn_id=l.split(fs)[1]
    for ii in range(1, len(l.split(fs))):
      pos=(l.split(fs)[ii][1:]).strip("\n")
      g=gprot[ii-1]
      for gg in g.split("/"):
        if gg in unipid_uac:
          for uac in unipid_uac[gg]:
            cgns[uac][pos]=[cgn, cgn_id]

pdb_list=[]
for l in open(sys.argv[1], "rt"):
#for l in open("../Gs_pdbs.txt", "rt"):
  pdb=l.split(fs)[0].strip("\n")
  pdb_list.append(pdb)
  
#print (pdb_list)
#print gpcrdb.keys()
#print cgns.keys()
int_stat={}
aa_stat={}
total_prot_pairs=[]
total_pos_pairs=[]
gpcrpos_list=[]
gprotpos_list=[]
ref_stat={}
for infile in glob.glob("./*_c*_c*_d*_contacts.txt"):
  #print infile
  for l in open(infile, "rt"):
    if l[0] != "#":
      s1_uac=l.split(fs)[0]
      s2_uac=l.split(fs)[1]
      s1_pdb=l.split(fs)[2]
      s2_pdb=l.split(fs)[3]
      #print (s1_pdb)
      if s1_pdb.split("/")[0] in pdb_list:
        #print (l)
        uac1=s1_uac.split("/")[0]
        uac2=s2_uac.split("/")[0] 
        pos1=GetID(s1_uac.split("/")[1])
        pos2=GetID(s2_uac.split("/")[1])
        num1=""
        num2=""
        flag1=""
        flag2=""
        if  uac1 in gpcrdb:
          site=[uac1,pos1]
          num1=GetGpcrdb(site,gpcrdb)
          flag1="GPCR"
        elif uac1 in cgns:
          site=[uac1,pos1]
          num1=GetCgn(site,cgns)[0]
          flag1="GPROT"
        if  uac2 in gpcrdb:
          site=[uac2,pos2]
          num2=GetGpcrdb(site,gpcrdb)
          flag2="GPCR"
        elif uac2 in cgns:
          site=[uac2,pos2]
          num2=GetCgn(site,cgns)[0]
          flag2="GPROT"
        #print uac1, pos1, s1_uac.split("/")[1][0], num1, flag1, uac2, pos2,  s2_uac.split("/")[1][0], num2, flag2
        flag=0
        if flag1 == "GPROT" and flag2 == "GPCR":
          gprotpos=num1
          gpcrpos=num2
          flag=1
        elif flag1 == "GPCR" and flag2 == "GPROT":
          gprotpos=num2
          gpcrpos=num1
          flag=1
        if num1 != "-" and num2 != "-" and flag == 1:
          if gpcrpos not in gpcrpos_list:
            gpcrpos_list.append(gpcrpos)
          if gprotpos not in gprotpos_list:
            gprotpos_list.append(gprotpos)
          pair=num1+"|"+num2
          invpair=num2+"|"+num1
          protpair=uac1+"|"+uac2
          invprotpair=uac2+"|"+uac1
          if pair not in total_pos_pairs and invpair not in total_pos_pairs:
            total_pos_pairs.append(pair) 
          if protpair not in total_prot_pairs and invprotpair not in total_prot_pairs:
            total_prot_pairs.append(protpair)
          if gpcrpos not in int_stat:
            int_stat[gpcrpos]={}
            int_stat[gpcrpos][gprotpos]=[protpair]
          elif gpcrpos in int_stat and gprotpos not in int_stat[gpcrpos]:
            int_stat[gpcrpos][gprotpos]=[protpair]
          elif protpair not in int_stat[gpcrpos][gprotpos] and invprotpair not in int_stat[gpcrpos][gprotpos]:
            int_stat[gpcrpos][gprotpos].append(protpair)
          if s1_pdb.find("mt1-gi1") != -1:
            if gpcrpos not in ref_stat:
              ref_stat[gpcrpos]={}
              ref_stat[gpcrpos][gprotpos]=[protpair]
            elif gpcrpos in ref_stat and gprotpos not in ref_stat[gpcrpos]:
              ref_stat[gpcrpos][gprotpos]=[protpair]
            elif protpair not in ref_stat[gpcrpos][gprotpos] and invprotpair not in ref_stat[gpcrpos][gprotpos]:
              ref_stat[gpcrpos][gprotpos].append(protpair)        
        print (uac1+fs+pos1+fs+s1_uac.split("/")[1][0]+fs+num1+fs+flag1+fs+s1_pdb+fs+uac2+fs+pos2+fs+s2_uac.split("/")[1][0]+fs+num2+fs+flag2+fs+s2_pdb.strip("\n"))

gprotpos_list=['G.HN.48', 'G.HN.51', 'G.HN.52', 'G.hns1.2', 'G.hns1.3', 'G.S1.2', 'G.S3.1', 'G.S3.3', 'G.h3s5.2', 'G.H4.3', 'G.H4.5', 'G.H4.22', 'G.H4.26', 'G.h4s6.3', 'G.h4s6.8', 'G.h4s6.9', 'G.h4s6.10', 'G.h4s6.11', 'G.h4s6.12', 'G.h4s6.20', 'G.S6.1', 'G.S6.2', 'G.S6.3', 'G.H5.2', 'G.H5.5', 'G.H5.6', 'G.H5.8', 'G.H5.9', 'G.H5.10', 'G.H5.11', 'G.H5.12', 'G.H5.13', 'G.H5.15', 'G.H5.16', 'G.H5.17', 'G.H5.19', 'G.H5.20', 'G.H5.21', 'G.H5.22', 'G.H5.23', 'G.H5.24', 'G.H5.25', 'G.H5.26']
gpcrpos_list=['1x27', '1x60', '2x36', '2x37', '2x38', '2x39', '2x40', '2x41', '2x42', '2x43', '2x60', '3x49', '3x50', '3x52', '3x53', '3x54', '3x55', '3x56', '4x38', '4x39', '4x40', '4x41', '5x47', '5x48', '5x58', '5x61', '5x62', '5x64', '5x65', '5x66', '5x67', '5x68', '5x69', '5x70', '5x71', '5x72', '5x73', '5x74', '5x75', '5x76', '6x23', '6x24', '6x25', '6x26', '6x27', '6x28', '6x29', '6x30', '6x31', '6x32', '6x33', '6x34', '6x35', '6x36', '6x37', '6x38', '6x39', '6x40', '7x53', '7x55', '7x56', '8x47', '8x48', '8x49', '8x50', '8x51', '8x52']

	
#gpcrpos_list=sorted(gpcrpos_list)
#print int_stat 
outstr="GPCRpos\t"
for gp in gprotpos_list:
  outstr+=gp+fs
#print outstr

mat=[]
y=0
x=0
x_ref=[]
y_ref=[]
for rec in gpcrpos_list:
  outstr=rec+fs
  row=[]
  x=0
  for gp in gprotpos_list:
    if rec in int_stat and gp in int_stat[rec]:
      fraction="%2.1f" % (float(len(int_stat[rec][gp]))/float(len(total_prot_pairs)))
      row.append(float(len(int_stat[rec][gp]))/float(len(total_prot_pairs)))
    else:
      fraction="0.0"
      row.append(0.0)
    outstr+=fraction+fs
    if rec in ref_stat and gp in ref_stat[rec]:
      x_ref.append(x)
      y_ref.append(y)
    x+=1
  #print outstr
  mat.append(row)
  y+=1


mat=np.array(mat, dtype=float)
#print x_ref
#print y_ref

###Plotting here
rc('font', weight='bold')
fig = plt.figure( figsize=(20,20))
#f, (ax0, ax1) = plt.subplots(2, sharex=True)
ax0 = fig.add_axes([0.175, 0.082,0.8,0.6])




####Plotting the contact matrix

#cmap = matplotlib.colors.ListedColormap(['white', 'dodgerblue', 'lightgreen', 'tomato'])
#bounds=[0, 1, 2, 3, 4]
#norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
img=ax0.imshow(mat, cmap='coolwarm', interpolation='nearest', aspect=1)

#sc=ax0.scatter(x_ref,y_ref, c="lime", edgecolors="lime")

minorLocator = IndexLocator(1,0.5)
majorLocator= IndexLocator(1,0)
minorLocatorY = IndexLocator(1,0)
#print minorLocatorY

#print (gpcrpos_list)
#print (gprotpos_list)

ax0.xaxis.set_minor_locator(minorLocator)
ax0.xaxis.set_major_locator(majorLocator)
ax0.yaxis.set_minor_locator(minorLocatorY)
ax0.set_xticklabels(gprotpos_list, rotation=90, minor=True)
major_ticksy = np.arange(0, len(gpcrpos_list), 1)
ax0.set_yticks(major_ticksy)
ax0.set_yticklabels(gpcrpos_list, rotation=0, ha="right")
ax0.tick_params(axis='x', which='minor', labelsize=12.5)
ax0.tick_params(axis='y', which='major', labelsize=15)
#ax1.tick_params(axis='y', which='major', labelsize=15)
ax0.xaxis.set_major_formatter(ticker.NullFormatter())
ax0.grid(axis="x", which='major', color="gray", linestyle="-")
ax0.grid(axis="y", which='minor', color="gray", linestyle="-")

axcolor = fig.add_axes([0.8, 0.3,0.01,0.15])
cbar=pylab.colorbar(img, cax=axcolor)
cbar.ax.tick_params(labelsize=15)

plt.draw()
plt.savefig('GPCR-Gs-contact_stat_comp.eps')
plt.savefig('GPCR-Gs-contact_stat_comp.png')
ax0.set_xlim([0,len(gprotpos_list)])
ax0.set_ylim([0,len(gpcrpos_list)])

sys.exit()
