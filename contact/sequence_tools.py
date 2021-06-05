#!/usr/bin/env python


import os, sys, operator, math, gzip, sqlite3, glob, itertools
import numpy as np

################################################################################################################################
####### A collection of python functions to parse various formats of text alignment files                                ######
#################################Written by Francesco Raimondi##################################################################
################################################################################################################################

f="\t"


def ExtractFasta(db, idslist, fastaout):
  alias1={}
  alias2={}
  alias3={}
  fastaoutfile=open(fastaout,"w")

  conn = sqlite3.connect(db)

  c = conn.cursor()
  ###We need to retrieve all the pdb2uniprot matches for the input XL
  sql_query=""
  #sql_query = 'select * from swissprot where ' + ' or '.join(("id_ac = '" + str(n)+"'"+" or id_id = '" + str(n)+"'"+" or id_gn = '" + str(n)+"'" for n in idslist)) 
  for n in idslist:
    ###This is a workaround I implemented to process as much data as possible from XlinkAnalyzer
    ###There are a lot of identifiers screwed up, i.e. they are lower case with no correspondence in the Uniprot
    ###while instead Uniprot cotains the corresponding uppercase identifier.
    sql_query = 'select * from uniprot where id_ac = '+"'"+ str(n)+"'"
    c.execute(sql_query)
    data=c.fetchall()
    if len(data) == 0:
      N=n.upper()
      sql_query = 'select * from uniprot where id_ac = '+"'"+ str(N)+"'"
      #c.execute(sql_query):
      #data=c.fetchall()
    print (sql_query)
    
    for line in c.execute(sql_query):
      uniac=line[0]
      unid=line[1]
      gn=line[2]
      fasta=line[3]
      print (gn, fasta)
      fastaoutfile.write(fasta)
      alias1[uniac]=unid
      alias2[uniac]=gn
      alias3[n]=uniac

  conn.close()
  
  return alias1, alias2, alias3


def GetBlastOut(infile, evalcutoff,nriter,overlap_flag):
  if infile.find(".gz") != -1:
    blastout=gzip.open(infile,"rt")
  else:
    blastout=open(infile,"rt")
    
  pdb_matches={}
  flag1=0
  flag2=0
  flag3=0
  start=0
  nrround=0
  protein_matches=0
  buffer_list={}
  for line in blastout:
    if line.find("Results from round")!= -1:

      start=1
      nrround=int((line.split()[3]).strip("\n"))
      continue
      
    if start == 1 and len(line.split()) > 0 and line.split()[0] == "Query=":
      flag1=1
      uniac=line.split()[1]
      unid=uniac.split("|")[2]
      uniac=uniac.split("|")[1]
      matches=0

    if start == 1 and line[0] == ">":
      chain=line.rsplit("_",1)[1]
      pdbid=line.rsplit("_",1)[0]
      chain=chain.split()[0]
      pdbid=(pdbid.lstrip(">")).lstrip()
      flag2=1
      flag3=0
      cont=0
      query_start=0
      query_seq=""
      query_stop=0
      sbjc_start=0
      sbjc_seq=""
      sbjc_stop=0
      f1=0
      f2=0
      matches+=1
      continue

    if start == 1 and flag1 == 1 and flag2 == 1:
      if len(line.split()) > 0 and line.split()[0] == "Score" and line.split()[1] == "=":
        Evalue=(line.split()[7]).strip(",")
        if Evalue[0] == "e":
          Evalue="1"+Evalue
        Evalue=float(Evalue)          
      if len(line.split()) > 0 and line.split()[0] == "Identities" and line.split()[1] == "=":
        identity=(line.split()[3]).strip(",")
        identity=identity.strip("(")
        identity=identity.strip(")")
        identity=identity.strip("%")
        if Evalue < evalcutoff:
          print (uniac, unid, line)
          flag3=1
          #flag2=0
          #continue
        #continue
      #continue

    if start == 1 and  len(line.split()) > 0 and line.split()[0] == "Query" and flag3 == 1:
      #print line
      query_start=line.split()[1]
      query_seq=line.split()[2]
      query_stop=line.split()[3]
      f1=1

    if start == 1 and len(line.split()) > 0 and line.split()[0] == "Sbjct" and flag3 == 1:
      #print line
      sbjc_start=line.split()[1]
      sbjc_seq=line.split()[2]
      sbjc_stop=line.split()[3]
      f2=1

    if start == 1 and flag3 == 1 and f1 == 1 and f2 == 1:
      ###Using targ_range and ref_range to dynamically monitor whether a given sequence alignment is overlapping with previously aligned regions
      ###If so, discard them to avoid confusing sequence-structure assignments which more often happen with repeat motives
      f1=0
      f2=0
      #print ("Save stuff!")
      if pdbid not in buffer_list:
        buffer_list[pdbid]={}
        buffer_list[pdbid][uniac]=[]
        buffer_list[pdbid][uniac].append([Evalue,int(query_start), int(query_stop),query_seq,int(sbjc_start), int(sbjc_stop), sbjc_seq, chain, identity])
        targ_range=[]
        for ii in range(int(query_start), int(query_stop)):
          targ_range.append(ii)
        ref_range=[]
        for ii in range(int(sbjc_start), int(sbjc_stop)):
          ref_range.append(ii)
      elif pdbid in buffer_list and uniac not in buffer_list[pdbid]:
        buffer_list[pdbid][uniac]=[]
        buffer_list[pdbid][uniac].append([Evalue,int(query_start), int(query_stop),query_seq,int(sbjc_start), int(sbjc_stop), sbjc_seq, chain, identity])
        targ_range=[]
        for ii in range(int(query_start), int(query_stop)):
          targ_range.append(ii)
        ref_range=[]
        for ii in range(int(sbjc_start), int(sbjc_stop)):
          ref_range.append(ii)
      else:
        if overlap_flag == 0:
          if int(query_start) not in targ_range and int(query_stop) not in targ_range and int(sbjc_start) not in ref_range and int(sbjc_stop) not in ref_range:
            buffer_list[pdbid][uniac].append([Evalue,int(query_start), int(query_stop),query_seq,int(sbjc_start), int(sbjc_stop), sbjc_seq, chain, identity])
            for ii in range(int(query_start), int(query_stop)):
              targ_range.append(ii)
            for ii in range(int(sbjc_start), int(sbjc_stop)):
              ref_range.append(ii)
        elif overlap_flag == 1:
          buffer_list[pdbid][uniac].append([Evalue,int(query_start), int(query_stop),query_seq,int(sbjc_start), int(sbjc_stop), sbjc_seq, chain, identity])
                
    if (start == 1 and line.find("Search has CONVERGED!") != -1) or (start == 1 and line.find("Effective search space used:") != -1 and nrround == nriter):
      for pid in buffer_list.keys():
        if pid not in pdb_matches:
          pdb_matches[pid]={}
        for uac in buffer_list[pid].keys():
          #print pdbid, uniac
          pdb_matches[pid][uac]=buffer_list[pid][uac]
      start=0
      #flag2=0
      buffer_list={}
      protein_matches+=1

  return pdb_matches

def GetBlastOutX(infile, evalcutoff,nriter,overlap_flag):
  blastout=open(infile,"r")
  pdb_matches={}
  flag1=0
  flag2=0
  flag3=0
  start=0
  nrround=0
  protein_matches=0
  buffer_list={}
  for line in blastout:
    if line.find("Results from round")!= -1:

      start=1
      nrround=int((line.split()[3]).strip("\n"))
      continue
      
    if start == 1 and len(line.split()) > 0 and line.split()[0] == "Query=":
      flag1=1
      uniac=line.split()[1]
      unid=uniac.split("|")[2]
      uniac=uniac.split("|")[1]
      matches=0

    if start == 1 and line[0] == ">":
      chain=line.rsplit("_",1)[1]
      pdbid=line.rsplit("_",1)[0]
      chain=chain.split()[0]
      pdbid=(pdbid.lstrip(">")).lstrip()
      flag2=1
      flag3=0
      cont=0
      query_start=0
      query_seq=""
      query_stop=0
      sbjc_start=0
      sbjc_seq=""
      sbjc_stop=0
      f1=0
      f2=0
      matches+=1
      continue

    if start == 1 and flag1 == 1 and flag2 == 1:
      if len(line.split()) > 0 and line.split()[0] == "Score" and line.split()[1] == "=":
        Evalue=(line.split()[7]).strip(",")
        if Evalue[0] == "e":
          Evalue="1"+Evalue
        Evalue=float(Evalue)          
      if len(line.split()) > 0 and line.split()[0] == "Identities" and line.split()[1] == "=":
        identity=(line.split()[3]).strip(",")
        identity=identity.strip("(")
        identity=identity.strip(")")
        identity=identity.strip("%")
        if Evalue < evalcutoff:
          print (uniac, unid, line)
          flag3=1
          Xcount=0
          #flag2=0
          #continue
        #continue
      #continue

    if start == 1 and  len(line.split()) > 0 and line.split()[0] == "Query" and flag3 == 1:
      #print line
      query_start=line.split()[1]
      query_seq=line.split()[2]
      query_stop=line.split()[3]
      f1=1

    if start == 1 and len(line.split()) > 0 and line.split()[0] == "Sbjct" and flag3 == 1:
      #print line
      sbjc_start=line.split()[1]
      sbjc_start2=int(line.split()[1])-Xcount
      sbjc_seq=line.split()[2]
      Xcount+=sbjc_seq.count("X")
      sbjc_stop=line.split()[3]
      sbjc_stop2=int(line.split()[3])-Xcount
      f2=1

    if start == 1 and flag3 == 1 and f1 == 1 and f2 == 1:
      ###Using targ_range and ref_range to dynamically monitor whether a given sequence alignment is overlapping with previously aligned regions
      ###If so, discard them to avoid confusing sequence-structure assignments which more often happen with repeat motives
      f1=0
      f2=0
      #print ("Save stuff!")
      if pdbid not in buffer_list:
        buffer_list[pdbid]={}
        buffer_list[pdbid][uniac]=[]
        buffer_list[pdbid][uniac].append([Evalue,int(query_start), int(query_stop),query_seq,int(sbjc_start), int(sbjc_stop), sbjc_seq, chain, identity, sbjc_start2, sbjc_stop2])
        targ_range=[]
        for ii in range(int(query_start), int(query_stop)):
          targ_range.append(ii)
        ref_range=[]
        for ii in range(int(sbjc_start), int(sbjc_stop)):
          ref_range.append(ii)
      elif pdbid in buffer_list and uniac not in buffer_list[pdbid]:
        buffer_list[pdbid][uniac]=[]
        buffer_list[pdbid][uniac].append([Evalue,int(query_start), int(query_stop),query_seq,int(sbjc_start), int(sbjc_stop), sbjc_seq, chain, identity, sbjc_start2, sbjc_stop2])
        targ_range=[]
        for ii in range(int(query_start), int(query_stop)):
          targ_range.append(ii)
        ref_range=[]
        for ii in range(int(sbjc_start), int(sbjc_stop)):
          ref_range.append(ii)
      else:
        if overlap_flag == 0:
          if int(query_start) not in targ_range and int(query_stop) not in targ_range and int(sbjc_start) not in ref_range and int(sbjc_stop) not in ref_range:
            buffer_list[pdbid][uniac].append([Evalue,int(query_start), int(query_stop),query_seq,int(sbjc_start), int(sbjc_stop), sbjc_seq, chain, identity, sbjc_start2, sbjc_stop2])
            for ii in range(int(query_start), int(query_stop)):
              targ_range.append(ii)
            for ii in range(int(sbjc_start), int(sbjc_stop)):
              ref_range.append(ii)
        elif overlap_flag == 1:
          buffer_list[pdbid][uniac].append([Evalue,int(query_start), int(query_stop),query_seq,int(sbjc_start), int(sbjc_stop), sbjc_seq, chain, identity, sbjc_start2, sbjc_stop2])
                
    if (start == 1 and line.find("Search has CONVERGED!") != -1) or (start == 1 and line.find("Effective search space used:") != -1 and nrround == nriter):
      for pid in buffer_list.keys():
        if pid not in pdb_matches:
          pdb_matches[pid]={}
        for uac in buffer_list[pid].keys():
          #print pdbid, uniac
          pdb_matches[pid][uac]=buffer_list[pid][uac]
      start=0
      #flag2=0
      buffer_list={}
      protein_matches+=1

  return pdb_matches
