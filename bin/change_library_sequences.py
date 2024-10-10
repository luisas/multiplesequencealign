#!/usr/bin/env python3

import pandas as pd
import sys
from Bio import SeqIO
from Bio.Seq import Seq


library = sys.argv[1]
mapping_file = sys.argv[2]
outfile = sys.argv[3]


mapping_df = pd.read_csv(mapping_file, sep = "\t", header = None)
mapping_df.columns = ["id", "seq", "foldseek_seq"]
mapping_df["id"] = mapping_df.id.str.replace("_alphafold.pdb", "", regex = False)
mapping_df.id = mapping_df.id.str.replace(".pdb", "", regex = False)
mapping_df.id = mapping_df.id.str.strip()


file1 = open(library, 'r')
Lines = file1.readlines()

count = 0
header_finished = False
f = open(outfile, "w")

# Strips the newline character
for line in Lines:
    
    if count == 0:
        count += 1
        f.write(line)
        continue
    elif count == 1: 
        nseq = int(line)
        f.write(line)
        count += 1
        continue
    elif count == nseq+2: 
        header_finished = True
    
    if not header_finished: 
        count += 1
        seq_id = line.split()[0]
        seq = line.split()[2]
        # This is to check if it stays the same 
        # print(mapping_df[mapping_df.foldseek_seq == seq])
        oldseq = str(mapping_df[mapping_df.id == seq_id].foldseek_seq.values[0])
        substitute_seq = str(mapping_df[mapping_df.id == seq_id].seq.values[0])
        #print(line)
        #print(oldseq)
        new_line = line.replace(oldseq, substitute_seq)
        f.write(new_line)
    else: 
        f.write(line)
f.close()