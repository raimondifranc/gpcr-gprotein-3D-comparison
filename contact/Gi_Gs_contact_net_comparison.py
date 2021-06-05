#!/data/SW/anaconda3/envs/myenv/bin/python

import sys, os, operator,gzip,glob, math
sys.path.insert(0, "/data/Users/francesco/code/")
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

def GetNet(inlist, clist):
  tot=[]
  outlist={}
  for l in open(inlist, "rt"):
    #print (l)
    uac1=l.split(fs)[0]
    aa1=l.split(fs)[2]
    pos1=l.split(fs)[3]
    class1=l.split(fs)[4]
    pdb1=l.split(fs)[5].split("/")[0]
    uac2=l.split(fs)[6]
    aa2=l.split(fs)[8]
    pos2=l.split(fs)[9]
    class2=l.split(fs)[10]
    pdb2=l.split(fs)[11].split("/")[0]
    #print (class1, class2)
    if (class1 == "GPROT" and class2 == "GPCR") or (class1 == "GPCR" and class2 == "GPROT"):
      if class1 == "GPROT":
        pair=pos1+"|"+pos2
        contact=aa1+"|"+aa2
        comp=uac1+"|"+uac2
      elif class1 == "GPCR":
        pair=pos2+"|"+pos1
        contact=aa2+"|"+aa1
        comp=uac2+"|"+uac1
      if pair not in outlist:
        outlist[pair]=[{},[comp]]
        outlist[pair][0][contact]=[comp]
      elif pair in outlist:
        if comp not in outlist[pair][1]:
          outlist[pair][1].append(comp)
        if contact not in outlist[pair][0]:
          outlist[pair][0][contact]=[comp]
        else:
          outlist[pair][0][contact].append(comp)
      if pair not in clist:
        clist.append(pair) 
      if comp not in tot:
        tot.append(comp) 
  return outlist,tot

consensus_net=[]
net1,tot1=GetNet("Gi_contacts.log",consensus_net)
net2,tot2=GetNet("Gs_contacts.log",consensus_net)

#print consensus_net
outlinks=open("links.txt", "w")
outnodes=open("nodes.txt", "w")
link_stat={}
node_stat={}
for pair in consensus_net:
  if pair in net1:
    contact_fraction1=float(len(set(net1[pair][1])))/float(len(set(tot1)))
    contact_sorted=sorted(net1[pair][0].items(), key=lambda a: len(set(a[1])), reverse=True)
    top_contact1=contact_sorted[0][0]
    top_contact_fraction1=float(len(set(contact_sorted[0][1])))/float(len(set(tot1)))
  else:
    contact_fraction1=0
    top_contact1="-"
    top_contact_fraction1=0.0
  if pair in net2:
    contact_fraction2=float(len(set(net2[pair][1])))/float(len(set(tot2)))
    contact_sorted=sorted(net2[pair][0].items(), key=lambda a: len(set(a[1])), reverse=True)
    top_contact2=contact_sorted[0][0]
    top_contact_fraction2=float(len(set(contact_sorted[0][1])))/float(len(set(tot2)))
  else:
    contact_fraction2=0
    top_contact2="-"
    top_contact_fraction2=0.0
  
  print (pair+fs+str(contact_fraction1)+fs+str(top_contact1)+fs+str(top_contact_fraction1)+fs+str(contact_fraction2)+fs+str(top_contact2)+fs+str(top_contact_fraction2))
  if pair.split("|")[0].find(".") != -1:
    ss1=pair.split("|")[0].rsplit(".",1)[0]
  elif pair.split("|")[1].find(".") != -1:
    ss1=pair.split("|")[1].rsplit(".",1)[0]
  if pair.split("|")[1].find("x") != -1:
    if pair.split("|")[1].find("34") != -1:
      ss2="ICL2"
    else:
      ss2="TM"+pair.split("|")[1].split("x")[0]
  elif pair.split("|")[0].find("x") != -1:
    if pair.split("|")[0].find("34") != -1:
      ss2="ICL2"
    else:
      ss2="TM"+pair.split("|")[0].split("x")[0]

  sspair=ss1+"|"+ss2
  if sspair not in link_stat:
    link_stat[sspair]=[[],[]]
  if contact_fraction1 != 0:
    link_stat[sspair][0].append(contact_fraction1)
  if contact_fraction2 != 0:
    link_stat[sspair][1].append(contact_fraction2)
  if ss1 not in node_stat:
    node_stat[ss1]=[[],[]]
  if ss2 not in node_stat:
    node_stat[ss2]=[[],[]]
  if contact_fraction1 != 0:
    node_stat[ss1][0].append(pair)
    node_stat[ss2][0].append(pair)
  if contact_fraction2 != 0:
    node_stat[ss1][1].append(pair)
    node_stat[ss2][1].append(pair)

outnodes.write("SS\tNr_contacts_Gi\tNr_contacts_Gs\n")
for node in node_stat.keys():
  outnodes.write(node+fs+str(len(set(node_stat[node][0])))+fs+str(len(set(node_stat[node][1])))+"\n")

outlinks.write("SS1\tSS2\tSS_contacts_Gi\tAvg_fraction_Gi\tSS_contacts_Gs\tAvg_fraction_Gs\n")
for link in link_stat.keys():
  ss1=link.split("|")[0]
  ss2=link.split("|")[1]
  avg1=0.0
  avg2=0.0
  avg1=np.average(link_stat[link][0])
  avg2=np.average(link_stat[link][1])
  if np.isnan(avg1):
    avg1 = 0.0
  if np.isnan(avg2):
    avg2 = 0.0
    
  outlinks.write(ss1+fs+ss2+fs+str(len(set(link_stat[link][0])))+fs+str(avg1)+fs+str(len(set(link_stat[link][1])))+fs+str(avg2)+"\n")
sys.exit()
