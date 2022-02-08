import pandas as pd
import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord
from Bio import SeqIO
from collections import defaultdict
import numpy as np

CHANNELS = ["DATA9", "DATA10", "DATA11", "DATA12"]
COLORS = ["black", "black", "black", "black"]


def blast_row_to_feature(row, color="#ffd700"):
    row_dict = row.to_dict(orient="records")
    assert len(row_dict) == 1
    row_dict = row_dict[0]

    start, end = row_dict["sstart"], row_dict["ssend"]
    # determine strand
    if start < end:
        strand = 1
    else:
        strand = -1
        end, start = row_dict["sstart"], row_dict["ssend"]
    return GraphicFeature(
        start=start, end=end, strand=strand, color=color, label=row_dict["read_type"]
    )


# add features to read in genbank file


# plot traces


def prepare_trace_dict(row):
    # read abi file, extract traces at every fifth value to get
    # basepair resolution of reads then slice to locations
    # of start and end of each alignment on the query and finally
    # pad the sequence based on the start of the alignment to the
    # the template
    row_dict = row.to_dict(orient="records")
    assert len(row_dict) == 1
    row_dict = row_dict[0]

    template_record = SeqIO.parse(row["template_gb"], "gb")
    template_len = len(template_record)

    def extract_traces(abi_path):
        # https://biopython.org/wiki/ABI_traces
        record = SeqIO.read(abi_path, "abi")
        trace = defaultdict(list)
        for c in CHANNELS:
            trace[c] = record.annotations["abif_raw"][c][::5]
        return trace

    def slice_trace(trace):
        alignment_read_span = [row_dict["qstart"], row_dict["qend"]]
        for each_channel in trace:
            trace[each_channel] = trace[each_channel][
                min(alignment_read_span) : max(alignment_read_span)
            ]
        return trace

    def pad_trace(trace):
        sub_span = [row_dict["qstart"], row_dict["qend"]]
        left = [0] * (template_len - min(sub_span))
        right = [0] * (template_len - max(sub_span))
        for each_channel in trace:
            trace[each_channel] = left + trace[each_channel] + right
        return trace


def big_plot(reads):

    fig, subs = plt.subplots(
        len(reads) + 1,
        1,
        figsize=(12, 3),
        sharex=True,
        gridspec_kw={"height_ratios": [4, 1]},
    )
    # plot traces

    def plot_trace(trace_dict, axis):
        for i, each_channel in enumerate(trace):
            axis.plot(trace[each_channel], color=COLORS[i])
        axis.axis("off")
        return axis

    for i, each_trace in enumerate(traces):
        subs[i] = plot_trace(each_trace)  # need way to set title (read_type)

    # finally plot the plasmid map with alignment locations
    record = SeqIO.read("example_sequence.gb", "genbank")
    graphic_record = BiopythonTranslator().translate_record(record)
    graphic_record.plot(ax=subs[-1], with_ruler=False, strand_in_label_threshold=4)

    return fig, subs


def main():

    # read the concatenated BLAST file
    blast_tab = pd.read_csv(snakemake.input["blast"])
    alignment_feats = [blast_row_to_feature(row) for i, row in blast_tab.iterrows()]
    traces = [prepare_trace_dict(row) for i, row in blast.iterrows()]
    big_plot(blast_tab, alignment_feats, traces)


if __name__ == "__main__":
    main()
