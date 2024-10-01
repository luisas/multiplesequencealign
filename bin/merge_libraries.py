#!/usr/bin/env python3

import pandas as pd
import sys
from Bio import SeqIO
from Bio.Seq import Seq
import sys

def aggsum(x):
    x = x.astype(float)
    x = x.sum()
    return(x.astype(str))

def aggmax(x):
    x = x.astype(float)
    x = x.max()
    return(str(x))

class SequenceEntry: 
    def __init__(self, id, length, sequence):
        self.id = id
        self.length = length
        self.sequence = sequence
        self.seq_3di = ""

    @classmethod
    def parse(cls, line):
        id,length,sequence = line.split()
        length = int(length)
        return cls(id, length, sequence)

    def print(self):
        print(self.id + "\t" + str(self.length) + "\t" + self.sequence)

    def write(self, file): 
        file.write(self.id + "\t" + str(self.length) + "\t" + self.sequence + "\n")

    def print3di(self):
        print(self.id + "\t" + str(self.length) + "\t" + self.seq_3di)

def flip_library_scores(library_scores):
    library_scores_flipped = []
    for library_score in library_scores: 
        library_scores_flipped.append(LibraryScore(library_score.res2, library_score.res1, library_score.score, library_score.extra2, library_score.extra1))
    return(library_scores_flipped)

class LibraryScore: 
    def __init__(self, res1, res2, score, extra1, extra2):
        self.res1 = res1
        self.res2 = res2
        self.score = score
        self.extra1 = extra1
        self.extra2 = extra2

    @classmethod
    def parse(cls, line):
        res1, res2, score, extra1, extra2 = line.split()
        res1 = int(res1)
        res2 = int(res2)
        score = float(score)
        return cls(res1, res2, score, extra1, extra2)

    def print(self):
        print(str(self.res1) + "\t" + str(self.res2) + "\t" + str(self.score) + "\t" + self.extra1 + "\t" + self.extra2)

    def to_df(self):
        df = pd.DataFrame([{'res1': self.res1, 'res2': self.res2, 'score': self.score, 'extra1': self.extra1, 'extra2':self.extra2}])
        return(df)

class LibraryEntry:
    def __init__(self, seq1, seq2, library_scores): 
        self.seq1 = seq1
        self.seq2 = seq2
        self.library_scores = library_scores

    def print(self): 
        print("#" + str(self.seq1) + "\t" + str(self.seq2))
        for score in self.library_scores:
            print(str(score.res1)+ "\t" + str(score.res2) + "\t" + str(score.score) + "\t" + str(score.extra1) + "\t" + str(score.extra2))

    def write(self, file):
        file.write("#" + str(self.seq1) + "\t" + str(self.seq2) + "\n")
        for score in self.library_scores:
            file.write(str(score.res1)+ "\t" + str(score.res2) + "\t" + str(score.score) + "\t" + str(score.extra1) + "\t" + str(score.extra2) + "\n")
    
    def add(self, library_scores, aggfunc):
        entries = pd.DataFrame()
        
        for score in self.library_scores:
            entries = pd.concat([entries, score.to_df()])
        
        for score in library_scores:
            entries = pd.concat([entries, score.to_df()])
        
        if aggfunc == "sum": 
            df = entries.groupby(["res1","res2"]).agg({"score": aggsum, "extra1": 'first', "extra2": 'first'}).reset_index()
        elif aggfunc == "max": 
            df = entries.groupby(["res1","res2"]).agg({"score": aggmax, "extra1": 'first', "extra2": 'first'}).reset_index()

        # Sort 
        df["res1"] = df["res1"].astype(int)
        df["res2"] = df["res2"].astype(int)
        df = df.sort_values(["res1","res2"])
        df_list = df.values.tolist()

        # Update the library scores
        self.library_scores = []
        for entry in df_list: 
            score = LibraryScore(entry[0], entry[1], entry[2], entry[3], entry[4])
            self.library_scores.append(score)
    
    def to_df(self):
        df = pd.DataFrame()
        for score in self.library_scores:
            df_current_score = score.to_df()
            df_current_score["seq1"] = self.seq1
            df_current_score["seq2"] = self.seq2
            df = pd.concat([df, df_current_score])
        return(df)

class Library: 
    def __init__(self, start_line, n_sequences,sequence_entries, library_entries, suffix):
        self.start_line = start_line
        self.n_sequences = n_sequences
        self.sequence_entries = sequence_entries
        self.library_entries = library_entries
        self.suffix = suffix



    @classmethod
    def parse(cls, library_file):
        # INIT
        header_finished = False
        sequence_entries = []
        library_scores = []
        library_entries = []
        suffix = []
        start_line = ""
        n_sequences = 0   

        file1 = open(library_file, 'r')
        Lines = file1.readlines()

        for index, line in enumerate(Lines): 
            # Parse first 2 lines
            if(index == 0):
                start_line = line.strip()
                continue
            elif(index == 1):
                n_sequences = int(line)
                continue
            # Here check if we are done with the header or not
            elif line.startswith("#"):
                if(header_finished):
                    library_entry = LibraryEntry(sequence_1_n, sequence_2_n, library_scores)
                    library_entries.append(library_entry)
                sequence_1_n, sequence_2_n = line.replace("#", "").split()
                sequence_1_n = int(sequence_1_n)
                sequence_2_n = int(sequence_2_n)
                library_scores = [] 
                
                header_finished = True
                continue    
            elif not header_finished: 
                # Here parse the sequence entries 
                seq_entry = SequenceEntry.parse(line)
                sequence_entries.append(seq_entry)
                continue  
            elif line.startswith("!") and header_finished == True:
                if suffix == []:
                    library_entry = LibraryEntry(sequence_1_n, sequence_2_n, library_scores)
                    library_entries.append(library_entry)
                suffix.append(line.strip())
                continue
            else:
                if not line.startswith("#") and not len(line) == 1:
                    library_score = LibraryScore.parse(line)
                    library_scores.append(library_score)
                    continue
        return cls(start_line, n_sequences, sequence_entries, library_entries, suffix)

    # Modification operations on library values
    def shift_residue_indeces(self, seq1, seq2, seq_2_shift, shift):
        for library_entry in self.library_entries:
            if (library_entry.seq1 == seq1 and library_entry.seq2 == seq2) or (library_entry.seq1 == seq2 and library_entry.seq2 == seq1):
                for library_score in library_entry.library_scores:
                    if library_entry.seq1 == seq_2_shift:
                        library_score.res1 = library_score.res1 + shift
                    elif library_entry.seq2 == seq_2_shift: 
                        library_score.res2 = library_score.res2 + shift
                       
                    
    def scale(self, factor, minval = -float("inf"), maxval = float("inf") ): 
        for library_entry in self.library_entries:
            for library_score in library_entry.library_scores:      
                newval =  library_score.score * float(factor)          
                library_score.score = max( minval, newval )
                library_score.score = min( maxval, library_score.score)

    def unify(self, newval): 
        for library_entry in self.library_entries:
            for library_score in library_entry.library_scores:      
                library_score.score = newval

    def unify_above(self, newval, threshold): 
            for library_entry in self.library_entries:
                for library_score in library_entry.library_scores:
                    if library_score.score > threshold:
                        library_score.score = newval

    def get_n_sequence_entry(self, sequence_id):
        i = 1 
        for sequence_entry in self.sequence_entries:
            if sequence_entry.id == sequence_id:
                return i
            i += 1

    def get_sequence_entry(self, n):
        i = 1 
        for sequence_entry in self.sequence_entries:
            if i == n:
                return sequence_entry
            i += 1


    def get_sequence_id(self, n):
       return( self.get_sequence_entry(n).id)     

    def to_df(self):
        entries = pd.DataFrame()
        for entry in self.library_entries:
            entries = pd.concat([entries, entry.to_df()])
        entries = entries.reset_index(drop=True)
        # reorder columns
        entries = entries[["seq1", "seq2", "res1", "res2", "score"]]
        # convert to int
        entries["res1"] = entries["res1"].astype(int)
        entries["res2"] = entries["res2"].astype(int)
        entries["score"] = entries["score"].astype(float)
        return entries

    def plot_library(self, n = 1 ):
        cmap = "flare"
        size_fig = 0.7
        df = self.to_df()
        df = df[df.seq1 == (n)]
        df["seq2"] = df["seq2"].astype(str)
        f, ax = plt.subplots(figsize=(17*size_fig,5*size_fig ))
        ax = sns.scatterplot(data = df, x = "res1", y = "seq2", s = 20, hue = "score", palette = cmap) 
        ax.legend(bbox_to_anchor=(1.02, 1.05))
        ax.set(xlabel = "residue of seq" + str(n), ylabel = "residue of other sequences" )
        return(ax)

    # Modification operations on sequences 
    def change_sequence_entries(self, fasta_file): 
        records = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        for sequence_entry in self.sequence_entries: 
            sequence_entry.sequence = str(records[sequence_entry.id].seq)
            sequence_entry.length = len(sequence_entry.sequence)

    def add_3d(self,fasta): 
        records = SeqIO.to_dict(SeqIO.parse(fasta, "fasta"))
        for sequence_entry in self.sequence_entries: 
            sequence_entry.seq_3di = str(records[sequence_entry.id].seq)

    # Combination operations on libraries
    def combine(self, library, aggfunc = "max"):
        # Itrate over all library entries in library to be added
        for library_entry in library.library_entries:
            seq1_n = library_entry.seq1
            seq2_n = library_entry.seq2
            # Get sequence ids for the library entry
            seq1 = library.get_sequence_id(seq1_n)
            seq2 = library.get_sequence_id(seq2_n)
            library_scores = library_entry.library_scores

            self.add_library_entry(seq1, seq2, library_scores, aggfunc)

    def add_library_entry(self, seq1, seq2, library_scores, aggfunc):
        # Identify if the library entry already exists in the current library
        for library_entry in self.library_entries:
            seq1_id = self.get_sequence_id(library_entry.seq1)
            seq2_id = self.get_sequence_id(library_entry.seq2)
            print(seq1_id, seq2_id)
            if seq1_id == seq1 and seq2_id == seq2:
                library_entry.add(library_scores, aggfunc)
                return
            elif seq1_id == seq2 and seq2_id == seq1:
                library_scores_flipped = flip_library_scores(library_scores)
                library_entry.add(library_scores_flipped, aggfunc)
                return

        library_entry = LibraryEntry(seq1, seq2, library_scores)
        self.library_entries.append(library_entry)


    # PRINTING
    def print_sequence_entries(self):
        for seq_entry in self.sequence_entries:
            seq_entry.print()

    def print_str_entries(self):
        for seq_entry in self.sequence_entries:
            seq_entry.print3di()

    def write_sequence_entries(self, file):
        for seq_entry in self.sequence_entries:
            seq_entry.write(file)

    def print_library_entries(self):
        for library_entry in self.library_entries:
            library_entry.print()

    def write_library_entries(self, file):
        for library_entry in self.library_entries:
            library_entry.write(file)
            
    def print(self):
        print(self.start_line)
        print(self.n_sequences)
        self.print_sequence_entries()
        self.print_library_entries()
        for suffix_line in self.suffix:
            print(suffix_line)

    def write(self, outfile):
        with open(outfile, 'w') as file:
            file.write(self.start_line)
            file.write('\n')
            file.write(str(self.n_sequences))
            file.write('\n')
            self.write_sequence_entries(file)
            self.write_library_entries(file)
            for suffix_line in self.suffix:
                file.write(suffix_line)
                file.write('\n')


library1 = sys.argv[1]
library2 = sys.argv[2]
aggfunc = sys.argv[3]
outname = sys.argv[4]

library1 = Library.parse(library1)
library2 = Library.parse(library2)

library1.combine(library2, aggfunc)

library1.write(outname)