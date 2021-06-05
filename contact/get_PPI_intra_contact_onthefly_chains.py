#!/data/SW/anaconda3/envs/myenv/bin/python

import os,sys,operator, math
import gzip, sqlite3
import numpy as np
import glob, itertools
from Bio import PDB
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC, generic_dna
from Bio.Align.Applications import ClustalwCommandline
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.PDB import PDBIO
from Bio.PDB.Polypeptide import *
from difflib import SequenceMatcher
sys.path.insert(0, "./")
from sequence_tools import ExtractFasta, GetBlastOut, GetBlastOutX

AA=['K','R','H','E','D','N','Q','Y','S','T','C','A','G','P','F','I','L','V','M','W']

f="\t"

class Element_type:

  def __init__(self):
    self.c_list=[]
    self.n_list=[]
    self.o_list=[]
    self.all_list=[]
    
  def assign_element(self, residue):
    for atom in residue:
      self.all_list.append(atom)
      if atom.element == "C":
        self.c_list.append(atom)
      elif atom.element == "O":
        self.o_list.append(atom)
      elif atom.element == "N":
        self.n_list.append(atom)
    return


def GetPdbRes(struct, ch):
  resmaps={}
  for m in struct:
    for c in m:
      if c.id == ch:
        ii=1
        for r in c:
          resmaps[ii]=r.id[1]
          ii+=1
  
  return resmaps
  
def GetID(aa):
  labels=""
  resid=""
  for c in aa:
    if c.isdigit() == False:
      labels+=c
    elif c.isdigit() == True:
      resid+=c

  if int(resid)>0:
    return int(resid)


def FindPos(aa, start, seq):
  idx=0
  #print (aa, start, seq)
  for s in seq:
    if int(start) == int(aa):
      #print "Hello"
      return idx
    #print s, start,idx,aa
    #if s != "-" and s != ".":
    if s != "-" and s != "." and s != "X":
      start+=1
    idx+=1


def GetPos(idx,start,seq):
  pos=0
  #print (idx, start, seq)
  for s in seq:
    #print (pos, idx)
    if int(pos) == int(idx):
      return start,s
    #if s != "-" and s != ".":
    if s != "-" and s != "." and s != "X":
      start+=1
    pos+=1
    
def CheckXRes(pos, refseq):
  aac=1
  Xc=0
  outstr=""
  for aa in refseq.split("\n")[1]:
    if aa == "X":
        Xc += 1
    outstr+=aa
    if int(aac) == int(pos):
      #print (refseq)
      #print (pos,aac, Xc,outstr)
      return pos-Xc
    aac += 1

####Maybe replace this with some similarity based on substitution matrices
def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

def MapVariants(qstart,qseq,sstart,sseq,var):
  simcutoff=0.5
  offset=5
  varidx=FindPos(var,sstart,sseq)
  #print (var, sstart, sseq, varidx)
  try:
    fampos,res=GetPos(varidx,qstart,qseq)
  except:
    return
  if similar(sseq[varidx:varidx+offset], qseq[varidx:varidx+offset]) >= simcutoff and sseq[varidx:varidx+offset] != "" and qseq[varidx:varidx+offset] != "":# or sseq[varidx:varidx-offset] == qseq[varidx:varidx-offset]:
    #print sseq[varidx:varidx+offset], qseq[varidx:varidx+offset] 
    return fampos, res
  elif  similar(sseq[varidx:varidx-offset], qseq[varidx:varidx-offset]) >= simcutoff  and sseq[varidx:varidx-offset] != "" and qseq[varidx:varidx-offset] != "":
    #print sseq[varidx:varidx-offset], qseq[varidx:varidx-offset] 
    return fampos, res
  else:
    #print sseq[varidx:varidx+offset], qseq[varidx:varidx+offset]
    return 0,""







def calc_residue_contact(r1, r2):
  """Returns the C-alpha distance between two residues"""
  contact=0
  for c1 in r1.c_list:
    for c2 in r2.c_list:
      diff_VdW  = c1.coord - c2.coord
      ###C-C distance(VdW) ) <= 5 A
      if np.sqrt(np.sum(diff_VdW * diff_VdW)) <= 5:
        contact=1
        return contact
 
  for o1 in r1.o_list:
    for n2 in r2.n_list:
      diff_HB = o1.coord - n2.coord
      diff_SB = o1.coord - n2.coord
      ###N-O distance (hydrogen bond) <= 3.5 A
      if np.sqrt(np.sum(diff_HB * diff_HB)) <= 3.5:
        contact=1
        return contact
         ###N-O distance (salt bridge) <= 5.5 A
      elif np.sqrt(np.sum(diff_SB * diff_SB)) <= 5.5:
        contact=1
        return contact
 
  ###Checking the other way round
  for n1 in r1.n_list:
    for o2 in r2.o_list:
      diff_HB = n1.coord - o2.coord
      diff_SB = n1.coord - o2.coord
      if np.sqrt(np.sum(diff_HB * diff_HB)) <= 3.5:
        contact=1
        return contact
      elif np.sqrt(np.sum(diff_SB * diff_SB)) <= 5.5:
        contact=1
        return contact
  
  return contact

def calc_residue_contact_coarse(r1, r2, dcutoff):
  contact=0
  for c1 in r1.all_list:
    for c2 in r2.all_list:
      diff_VdW  = c1.coord - c2.coord
      ###C-C distance(VdW) ) <= 5 A
      #print r1, r2
      if np.sqrt(np.sum(diff_VdW * diff_VdW)) <= dcutoff:
        
        contact=1
        break
  return contact
  

def calc_dist_matrix(chainA, chainB, dcutoff):
  ###Returns an array with n residues length with 1 for each residues in contact with the partner chain
  contact = np.zeros(len(chainA),np.float)
  clist={}
  for c, res1 in enumerate(chainA):
    res1_atm=Element_type()
    res1_atm.assign_element(res1)
    clist[res1]=[]
    for cc, res2 in enumerate(chainB):
      if (res2.id[1]-res1.id[1]) >5 or (res1.id[1]-res2.id[1]) >5:
        res2_atm=Element_type()
        res2_atm.assign_element(res2)
        #cont=calc_residue_contact(res1_atm, res2_atm)
        cont=calc_residue_contact_coarse(res1_atm, res2_atm, dcutoff)
        contact[c]+=cont
        if cont > 0 and res2 not in clist[res1]:
          clist[res1].append(res2)
  return contact, clist


def aaLabel(resname):
  if resname in standard_aa_names:
    aaout=three_to_one(resname)
  else:
    aaout=resname

  return aaout

def GetResIdx(chain, residue):
  cc=0
  for rr in chain:
    if rr.id[1] == residue.id[1]:
      return cc  
    cc+=1


def GetFasta(infile, outseqs):
  seq=""
  flag=0
  tracker=[]
  if infile[:-3] == ".gz":
    infile_stream=gzip.open(infile,"r")
  else:
    infile_stream=open(infile, "r")
  for l in infile_stream:
    if l[0] == ">":
      flag=0
      seq=""
      #gn=l.split()[0]
      seqid=l.split()[0].lstrip(">")
      outseqs[seqid]=l
    if l[0] != ">":
      outseqs[seqid]+=l.strip("\n")


#####################################################################################################################################################################################
if len(sys.argv) != 5:
  print ("Calculate contacts between specified chains of a given pdb\nUsage: pdb, chain 1, chain 2, distance cutoff\nIt checks all the contacts between chain-pairs in a given pdb")
  print ("It extract the sequences from the PDB, generates a database for blast searches, and then align uniprot sequences on it")
  print ("Remember to remove from the PDB insertions (as in engineered GPCRs,e.g. 3SN6, otherwise the Biopython ca-atom extraction will just fail-improve this part)")
  sys.exit()

pdb=sys.argv[1]
pdbid=pdb.split(".cif.gz")[0]
pdbid=pdbid.lower()
c1_ref=sys.argv[2]
c2_ref=sys.argv[3]
discutoff=float(sys.argv[4])

###Getting uniprot accession corresponding to each PDB chain
unip_list=[]
pdb2uac={}
for l in gzip.open("/data/DB/SIFTS/pdb_chain_uniprot.tsv.gz", "rt"):
  if l[0] != "#":
    pid=l.split(f)[0]
    chain=l.split(f)[1]
    uac=l.split(f)[2]
    if pid == pdbid:
      if pid not in pdb2uac:
        pdb2uac[pid]=[[chain,uac]]
      elif [chain, uac] not in pdb2uac[pid]:
        pdb2uac[pid].append([chain, uac])
      unip_list.append(uac)


####Retrieving corresponding fasta sequences
fastaname=("uniprot_query.fasta")

####Running the script for any use case requires using an internal SQLite version of Uniprot
####Plase contact francesco.raimondi@sns.it if you require any help 
#pac2pid, pac2gn,origid2pac = ExtractFasta("/data/DB/uniprot/uniprot_seq_ids_new.db", unip_list, fastaname)


####I use here the same workaround employed for XL_interpreter when applied to user input structures (i.e. not taken from the PDB)
### I create a db for blast searches with the sequences from the input PDB and I blast on hit SPROT sequences 
#pdbdir="/net/home.isilon/ds-russell/pdb/%s/pdb%s.ent.gz" % (pdb[1:3],pdb)
pdbdir="%s" % (pdb)
pdbpath=open(pdbdir, "rt")

seq_list=[]
for record in SeqIO.parse(pdbpath, "cif-atom"):
  pdb=pdb.strip()
  print (pdbid+"_"+record.annotations["chain"], record.seq)
  seq=SeqRecord(record.seq, id=pdbid+"_"+record.annotations["chain"], description="Sequence from PDB")
  seq_list.append(seq)

SeqIO.write(seq_list, "%s.fa" % (pdbid), "fasta")
os.system("makeblastdb -dbtype prot -in %s.fa" % (pdbid))

###Aligning uniprot sequences with PDB ones: this is to translate back to uniprot sequences the structural information

blastout_file="%s_blast.out" % (pdb)
Eval=10
psiblast_nriter=1

if os.path.exists(blastout_file):
  pass
else:
  os.system("psiblast -evalue %f -db %s.fa -query %s -num_iterations %d > %s -num_threads 3" % (Eval, pdbid, fastaname, psiblast_nriter, blastout_file))

###I can improve this part by considering only the best sequence match to a structure based on E-value

pdb2uniprot=GetBlastOut(blastout_file,Eval,psiblast_nriter,1)
#pdb2uniprot=GetBlastOutX(blastout_file,Eval,psiblast_nriter,1)
#print pdb2uac[pdbid]

print (pdb2uniprot)
outfile=open("%s_c%s_c%s_d%s_contacts.txt" % (pdbid, c1_ref, c2_ref, discutoff),"w")
CASEQ={}
fastas=pdbid+".fa"
GetFasta(fastas, CASEQ)

if os.path.isfile(pdbdir):
  #print pdb
  Pdb=pdb.upper()
  pdbpath=gzip.open(pdbdir, "rt")
  parser = PDB.MMCIFParser(QUIET=True)
  struct = parser.get_structure('XXX',pdbpath)
  model = struct[0]
  
  for rr in range(len(pdb2uac[pdbid])):
    for ss in range(rr, len(pdb2uac[pdbid])):
      #print (pdb2uac[pdbid][rr], pdb2uac[pdbid][ss])
      c1=pdb2uac[pdbid][rr][0]
      uac1=pdb2uac[pdbid][rr][1]
      c2=pdb2uac[pdbid][ss][0]
      uac2=pdb2uac[pdbid][ss][1]
      #contact_list = calc_dist_matrix(model[c1], model[c2])
      outstr=""

      if (c1 == c1_ref or c1 == c2_ref) and (c2 == c2_ref or c2 == c1_ref):
        if uac1 in pdb2uniprot[pdbid] and uac2 in pdb2uniprot[pdbid]:
          #print "YOOOO!!"
          contact_list, contact_res = calc_dist_matrix(model[c1], model[c2], discutoff)
          ###Define a function when chain residues are enumerated once for all and the correspondences stored in a dictionary to be called afterwards.
          ###This should speed up the analysis
          #residx=GetResIdx([model[c1], model[c2]])
          if (sum(contact_list) > 0):
            for ii, res in enumerate(model[c1]):
              aa=""
              r1=""
              outstr=""
              idx=0
              aa1=""
              r2_string=""
              r2_str=""
              aa2=""
              if contact_list[ii] > 0:
                aa=aaLabel(res.get_resname())
                #print ii, aa, res
                if ii == 0 and aa == "ACE":
                  kk=0
                elif ii >= 0:
                  kk=1
                r1=aa+str(res.id[1])
                
                #print (pdbid, uac1, pdb2uniprot[pdbid][uac1])
                for pseq in pdb2uniprot[pdbid][uac1]:
                  seq_chain1=CASEQ[pdbid+"_"+c1]
                  start_x=CheckXRes(pseq[4], seq_chain1)
                  stop_x=CheckXRes(pseq[5], seq_chain1)
                  if int(ii+kk) in range(start_x, stop_x+1):
#def MapVariants(tstart,tseq,rstart,rseq,start_count,simcutoff,seqchain):
                    #idx, aa1=MapVariants(pseq[1],pseq[3],pseq[4],pseq[6],ii+kk)
                    try:
                      idx, aa1=MapVariants(pseq[1],pseq[3],start_x,pseq[6],ii+kk)
                    except:
                      continue
                    #print aa1,idx
                    if aa1 != None and aa1 != "":
                      #print (c1, ii+kk, res, aa, r1, pseq[1],pseq[3], pseq[4],pseq[5]+1, pseq[6], aa1, idx)
                      break
                contacts=[]
                for r in set(contact_res[res]):
                  aa2=aaLabel(r.get_resname())
                  r2_string+=aa2+str(r.id[1])+","
                  jj=GetResIdx(model[c2],r)
                  for pseq2 in pdb2uniprot[pdbid][uac2]:
                    #print pdbid, uac2, pdb2uniprot[pdbid][uac2]
                    idx2=0
                    aa2_unip=""
                    #print (ii, res, pseq)
                    #print (jj, r, pseq2)
                    seq_chain2=CASEQ[pdbid+"_"+c2]
                    start_x2=CheckXRes(pseq2[4], seq_chain2)
                    stop_x2=CheckXRes(pseq2[5], seq_chain2)
                    if aa1 in AA and int(jj+kk) in range(start_x2, stop_x2+1):
                      #idx2, aa2_unip=MapVariants(pseq2[1],pseq2[3],pseq2[4],pseq2[6],jj+kk)
                      try:
                        idx2, aa2_unip=MapVariants(pseq2[1],pseq2[3],start_x2,pseq2[6],jj+kk)
                      except:
                        continue
                      if aa2_unip != "" and aa2_unip in AA:
                        #print ("MATCH: ",c2, jj+kk, r, aa2_unip, r2_string, pseq2[1],pseq2[3], pseq2[4],pseq2[5]+1, pseq2[6], aa2_unip, idx2)
                        r2=aa2_unip+str(idx2)+aa2_unip
                        r2_str+=r2+","
                        contacts.append(r2)
                        if aa in AA:
                          outfile.write(uac1+"/"+aa1+str(idx)+aa1+"\t"+uac2+"/"+r2+"\t"+pdbid+"/"+c1+"/"+r1+"\t"+pdbid+"/"+c2+"/"+aa2+str(r.id[1])+"\n")
                          
                      else:
                        continue
                    else:
                      continue
                #outfile1.write(Pdb+"/"+c1+"/"+r1+"\t"+str(len(set(contact_res[res])))+"\t"+Pdb+"/"+c1+"\t"+r2_string+"\n")
                #if aa1 == "" or len(aa) > 1:
                #  continue
                #elif aa1 in AA:
                #  outstr+=uniac+"/"+aa1+str(idx)+aa1+"\t"+str(len(set(contacts)))+"\t"+Pdb+"/"+c1+"\t"+r2_str+"\n"
                #outfile2.write(outstr)



#print dist_matrix

sys.exit()


