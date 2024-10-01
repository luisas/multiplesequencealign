#!/usr/bin/env python3

import os
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import sys 
mapping_file = sys.argv[1]
output_dir = sys.argv[2]

# Load mapping
mapping = pd.read_csv(mapping_file, sep = "\t", header = None, engine = "python")
mapping.columns = ["id", "seq", "map"]
mapping.id = mapping.id.str.replace("_alphafold.pdb", "", regex = False)
mapping['id'] = mapping['id'].str.split('.pdb').str[0]

# create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# loop through each row of the dataframe
for index, row in mapping.iterrows():
    # create a new SeqRecord object
    seq_record = SeqRecord(
        seq=Seq(row['map']),
        id=row['id'],
        description=''
    )
    # create a new file with the id as the file name in the output directory
    file_path = os.path.join(output_dir, row['id'] )
    with open(file_path, 'w') as f:
        # write the SeqRecord to the file in FASTA format
        SeqIO.write(seq_record, f, 'fasta')