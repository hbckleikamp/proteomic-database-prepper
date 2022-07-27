

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 17 12:01:51 2021

@author: hugokleikamp
"""

#%% change directory to script directory (should work on windows and mac)
import os
from pathlib import Path
from inspect import getsourcefile
os.chdir(str(Path(os.path.abspath(getsourcefile(lambda:0))).parents[0]))
script_dir=os.getcwd()
print(os.getcwd())

basedir=os.getcwd()

#%% import 

import Bio
from Bio import SeqIO
import pandas as pd
import itertools
import random



#%% Options

Bacterial_only=True   # retain only Bacteria (and Archaea(!)) in database  
Equate_IL=True        # change I and J into L 
Remove_ambiguous=True # remove ambiguous amino acids "B","O","U","X","Z","[","(" , and J in case IL is not equated
Add_decoy=True        # append decoy of reversed peptides

#%% Variables


Path_to_db="path/to/db"          # path to protein database that uses NCBI taxonomy


#for Bacterial_only
Path_to_taxonomy="path/to/taxonomy" # path to parsed NCBI taxonomy (as created by https://github.com/hbckleikamp/NCBI2Lineage)
Taxid_delimiter="TaxID=" #Important variable! OX= for UniProtKB, TaxID= for Uniref50, check how this is done

#for Remove_ambiguous
Ambiguous_AAs=["B","O","U","X","Z","[","("]

#Add decoy
decoy_delimiter="decoy_"
decoy_method="reverse" #or "scramble"


#%% Functions

def chunk_gen(it,size=10**6):
    c=itertools.count()
    for _,g in itertools.groupby(it,lambda _:next(c)//size):
        yield g
        
#%%

#parse output_path
Output_path=Path_to_db
if Bacterial_only:   Output_path=Output_path.replace(".fasta","_BacArch.fasta")
if Remove_ambiguous: Output_path=Output_path.replace(".fasta","_NoAmb.fasta")
if Equate_IL:        Output_path=Output_path.replace(".fasta","_IJeqL.fasta")
if Add_decoy:        Output_path=Output_path.replace(".fasta","_Decoy.fasta")

# tax database and files
if Bacterial_only:
    ranks=["superkingdom","phylum","class","order","family","genus","species"] 
    taxdf=pd.read_csv(Path_to_taxonomy,sep="\t")
    taxdf.columns=['OX']+ranks+["OS"]
    taxdf=taxdf[(taxdf["superkingdom"]=="Bacteria") | (taxdf["superkingdom"]=="Archaea")].astype(str)

if not Equate_IL: Ambiguous_AAs+=["J"]

#read    
recs=SeqIO.parse(Path_to_db,format="fasta")
chunks=chunk_gen(recs)

#write IL datbase
print("writing "+Path(Output_path).stem)
with open(Output_path,"w+") as f:

    for ic,c in enumerate(chunks):
        print("chunk "+ str(ic))
        
        taxs=[]
        rs=[]
        for r in c:
            
            s=str(r.seq)
            d=r.description
            
            if Bacterial_only:
                if Taxid_delimiter in d:
                    tax=r.description.split(Taxid_delimiter)[1].split()[0] #get taxonomic identifiers
                    taxs.append(tax)
                else: print("taxonomy delimiter not found!")
            
            if Equate_IL:
                s=s.replace("I","L").replace("J","L")
            
            if Remove_ambiguous:
                if sum([i in s for i in Ambiguous_AAs])!=0: 
                    continue
                
            if Add_decoy:
                if decoy_method=="reverse":
                    rs.append([decoy_delimiter+d,s[::-1]])
                if decoy_method=="scramble":
                    rs.append([decoy_delimiter+d,''.join(random.sample(s, len(s)))])
                if Bacterial_only:
                    taxs.append(tax)
            
            rs.append([d,s])
            
        
        if Bacterial_only:
            rs=pd.Series(rs)[pd.Series(taxs).isin(taxdf.OX)].tolist() #select only those that have superkingdom Archaea or Bacteria

        s="\n".join([">"+"\n".join(r) for r in rs])
        f.write(s+"\n")
    
      
#%%
