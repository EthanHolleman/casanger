from Bio import SeqIO
import pandas as pd
from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord

def make_kmers(string, size):
    for i in range(len(string)-size):
        yield string[i:i+size], i, i+size


def is_palindrome(candidate):
    if candidate == candidate[::-1]:
        return True
    else:
        return False


def find_all_palindromes(string, sizes):
    palindromes = []
    for each_size in sizes:
        for each_kmer, start, end in make_kmers(string, each_size):
            if is_palindrome(each_kmer):
                palindromes.append(
                    {
                        'sequence': each_kmer,
                        'start': start,
                        'end': end,
                        'size': each_size
                     }
                )
    return pd.DataFrame(palindromes)



def palindrome_to_feature(pal, color=None):
    return GraphicFeature(
        start=pal['start'], end=pal['end'], color=color, label=pal['sequence']
    )


def main():
    
    record = SeqIO.read(snakemake.params['genbank'], 'genbank')
    # Target site location
    search_seq = str(record.seq[target_start-extend:target_end+extend])
    pal_df = find_all_palindromes(search_seq, snakemake.params['sizes'])
    
    
    
    
    

