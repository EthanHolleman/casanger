from Bio import SeqIO
import numpy as np
import pandas as pd
from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
from colour import Color


def make_kmers(string, size):
    for i in range(len(string) - size):
        yield string[i : i + size], i, i + size


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
                        "sequence": each_kmer,
                        "start": start,
                        "end": end,
                        "size": each_size,
                    }
                )
    return pd.DataFrame(palindromes)


def modify_features(graphic_record):
    mod_features = []
    for f in graphic_record.features:
        if f.label != None:
            if "source" == f.label:
                continue
            if "|" in f.label:
                f.label = f.label.split("|")[0]
        mod_features.append(f)
    graphic_record.features = mod_features
    return graphic_record


def dot_plot(seq):
    
    def delta(x,y):
        return 0 if x == y else 1
    
    def M(seq1,seq2,i,j,k):
        return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))


    def makeMatrix(seq,k):
        n = len(seq)
        m = len(seq)
        return [[M(seq,seq,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]

    return np.array(makeMatrix(seq,1))



def plot_palindrome_map(palindromes, template_record, target_start, target_end, extended, seq):
    # genbank something here and then add all features
    lengths = set(palindromes["size"])

    def palindrome_to_feature(pal, color="#4287f5"):
        return GraphicFeature(start=pal["start"], end=pal["end"], color=color)

    search_start, search_end = target_start - extended, target_end + extended

    palindromes["start"] += search_start 
    palindromes["end"] += search_start

    start_color = Color("#4287f5")
    end_color = Color("#51f542")

    gradient = [c.hex for c in start_color.range_to(end_color, len(palindromes))]

    palindrome_feats = [
        palindrome_to_feature(p, gradient.pop()) for i, p in palindromes.iterrows()
    ]
    
    # make feaure for target site
    target_feat = GraphicFeature(
        start=target_start, end=target_end, label='Target sequence', 
        color='#c95a55'
        )

    # set different color for each palindrome length
    fig, (ax1, ax2) = plt.subplots(
        2,  # traces plus plasmid map
        1,
        figsize=(15, 15),
        sharex=False,
        gridspec_kw={"height_ratios": [1, 1.25]}
    )
    graphic_record = BiopythonTranslator().translate_record(template_record)
    print(palindrome_feats)
    graphic_record.features.append(target_feat)
    graphic_record.features += palindrome_feats
    graphic_record = modify_features(graphic_record)
    print(graphic_record.features)
    print(search_start, search_end)
    crop = graphic_record.crop((search_start, search_end))
    
    dotplot_matrix = dot_plot(seq)
    ax2.imshow(dotplot_matrix)
    
    crop.plot(ax=ax1, with_ruler=True, strand_in_label_threshold=4,
               plot_sequence=True, annotate_inline=True,
               elevate_outline_annotations=False)
    
    
    print(ax1, ax2)
    
    return fig



def main():

    record = SeqIO.read(snakemake.params["genbank"], "genbank")
    # Target site location
    target = pd.read_csv(snakemake.input["target"], sep="\t").iloc[0]
    target_start, target_end = target["start"], target["end"]
    extend = 50

    search_seq = str(record.seq[target_start-extend:target_end+extend])
    pal_df = find_all_palindromes(search_seq, snakemake.params["sizes"])

    fig = plot_palindrome_map(pal_df, record, target_start, target_end, extend, search_seq)
    fig.savefig(snakemake.output["png"])


if __name__ == "__main__":
    main()

