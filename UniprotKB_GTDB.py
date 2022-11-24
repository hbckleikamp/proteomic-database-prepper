# -*- coding: utf-8 -*-
"""
Created on Wed Nov 23 11:28:47 2022

@author: ZR48SA
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
import requests
import ftputil, urllib, gzip, zipfile, shutil, tarfile
import Bio
from Bio import SeqIO

import urllib
import ftputil
import requests
import zipfile
from openpyxl import Workbook, load_workbook 
import datetime
import pandas as pd
import time


#Variables
Path_to_db="H:/Databases/UniprotKB/UniprotKB.fasta"         # path to protein database that uses NCBI taxonomy


#for Bacterial_only
Path_to_taxonomy="" #unused
Taxid_delimiter="OX=" #unused

#for Remove_ambiguous
Ambiguous_AAs=["B","X","Z","[","("] #"O","U","*"

#Add decoy
decoy_delimiter="" #unused
decoy_method="" #unused


import Bio
from Bio import SeqIO
import pandas as pd
import itertools
import random
import subprocess

# Options
GTDB_only=True  # retain only Bacteria (and Archaea(!)) in database  
Equate_IL=True        # change I and J into L 
Remove_ambiguous=True # remove ambiguous amino acids "B","O","U","X","Z","[","(" , and J in case IL is not equated
No_Fragments=True     # remove incomplete sequences from UniprotKB contain (Fragment)/(Fragments) in header
Add_decoy=False        # append decoy of reversed or scrambled peptides
decoy_method="reverse" #or "scrambled

# Functions

def chunk_gen(it,size=10**6):
    c=itertools.count()
    for _,g in itertools.groupby(it,lambda _:next(c)//size):
        yield g

#parse output_path
Output_path=Path_to_db
if GTDB_only:   Output_path=Output_path.replace(".fasta","_GTDB.fasta")
if Remove_ambiguous: Output_path=Output_path.replace(".fasta","_NoAmb.fasta")
if No_Fragments:     Output_path=Output_path.replace(".fasta","_NoFrag.fasta")
if Equate_IL:        Output_path=Output_path.replace(".fasta","_IJeqL.fasta")
if Add_decoy:        Output_path=Output_path.replace(".fasta","_Decoy.fasta")

#add minimum length?

# GTDB tax database and files
ranks=["superkingdom","phylum","class","order","family","genus","species"] 
metadata=[
"H:/Databases/GTDB/Metadata/ar53_metadata_r207.tsv",
"H:/Databases/GTDB/Metadata/bac120_metadata_r207.tsv"]
taxdf=pd.concat([pd.read_csv(i,sep="\t",usecols=["ncbi_taxid","gtdb_taxonomy"]) for i in metadata]).drop_duplicates()
taxdf["ncbi_taxid"]=taxdf["ncbi_taxid"].astype(str)
taxdf[ranks]=taxdf["gtdb_taxonomy"].str.rsplit(";",expand=True)

sizes=taxdf.groupby("ncbi_taxid").size()
taxdf=taxdf[~taxdf["ncbi_taxid"].isin(sizes[sizes>1].index)] #remove dump taxa (ncbi taxa that map to multiple gtdb taxa)


#%%

if not Equate_IL: Ambiguous_AAs+=["J"]

#read    
recs=SeqIO.parse(Path_to_db,format="fasta")
chunks=chunk_gen(recs)

#write IL datbase
print("writing "+Path(Output_path).stem)
with open(Output_path,"w+") as f:

    for ic,c in enumerate(chunks):
        print("chunk "+ str(ic))

        chunk_df=pd.DataFrame([[str(r.seq),r.description] for r in c],columns=["seq","description"])
    
        if Equate_IL:        chunk_df["seq"]=chunk_df["seq"].str.replace("I","L").str.replace("J","L")
        if Remove_ambiguous: chunk_df=chunk_df[~pd.concat([chunk_df["seq"].str.contains(aa,regex=False) for aa in Ambiguous_AAs],axis=1).any(axis=1)]
        if No_Fragments:     chunk_df=chunk_df[~chunk_df["description"].str.contains("(Fragment",regex=False)]
        if GTDB_only:        chunk_df=chunk_df[chunk_df["description"].str.split(Taxid_delimiter).apply(lambda x: x[-1].split(" ")[0]).isin(taxdf.ncbi_taxid)]
        
        if Add_decoy:
            decoy=chunk_df.copy()
            if decoy_method=="scramble": decoy["seq"]=decoy.seq.apply(lambda x: ''.join(random.sample(x,len(x))))
            if decoy_method=="reverse":  decoy["seq"]=decoy.seq.apply(lambda x: x[::-1])
            decoy["description"]=decoy_delimiter+decoy["description"]
            chunk_df=pd.concat([chunk_df,decoy])
            
        f.write("\n"+"\n".join(">"+chunk_df["description"]+"\n"+chunk_df["seq"]))

#Cleanup
#os.remove(Path_to_db)

