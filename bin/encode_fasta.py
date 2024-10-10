#!/usr/bin/env python3

import pandas as pd
import sys
from Bio import SeqIO
from Bio.Seq import Seq


fasta_file = sys.argv[1]
mapping_file = sys.argv[2]
output_file = sys.argv[3]

# Load mapping
mapping = pd.read_csv(mapping_file, sep = "\t", header = None, engine = "python")
mapping.columns = ["id", "seq", "map"]
mapping.id = mapping.id.str.replace(".pdb", "", regex = False)
mapping.id = mapping.id.str.strip()
with open(fasta_file) as original, open(output_file, 'w') as mapped:
    records = SeqIO.parse(fasta_file, 'fasta')
    for record in records:
        fasta_seq = record.seq.strip()
        mapping_seq = mapping[mapping.id == record.id].seq.item().strip()
        if  fasta_seq == mapping_seq :
            record.seq = Seq(mapping[mapping.id == record.id]["map"].item())
        else:
            print("SEQUENCES ARE NOTE MATCHING!!")
        SeqIO.write(record, mapped, 'fasta')