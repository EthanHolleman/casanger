import pandas as pd
import matplotlib.pyplot as plt
from dna_features_viewer import BiopythonTranslator, GraphicFeature, GraphicRecord
from Bio import SeqIO
from collections import defaultdict
import numpy as np
from colour import Color

CHANNELS = ["DATA9", "DATA10", "DATA11", "DATA12"]  # nucleotide channels (traces)
COLORS = ["black", "red", "green", "blue"]  # channel colors 


# blast tab output column names
COL_NAMES = (
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "evalue",
    "bitscore",
)


def blast_row_to_feature(row, label, color="#ffd700"):
    """Convert a row (passed as a pandas df) of blast tabular output
    to a GraphicRecord instance.

    Args:
        row (DataFrame): Row of BLAST tabular output as pd.DataFrame
        label (str): Column name containing value to use as the feature label.
        color (str, optional): Color of feature. Defaults to "#ffd700".

    Returns:
        GraphicRecord: Alignment as a GraphicRecord instance.
    """

    row_dict = dict(row)

    start, end = row_dict["sstart"], row_dict["send"]
    # determine strand
    if start < end:
        strand = 1
    else:
        strand = -1
        end, start = row_dict["sstart"], row_dict["send"]
    return GraphicFeature(
        start=start, end=end, strand=strand, color=color, label=row_dict[label]
    )


def local_alignment_to_feature(row, color="#ffd700"):
    """Analogous to `blast_row_to_feature` but converts local alignments
    produced by `shortLocal.py` output to GraphicRecords.   

    Args:
        row (DataFrame): Row of BLAST tabular output as pd.DataFrame.
        color (str, optional): Color of feature. Defaults to "#ffd700".

    Returns:
        GraphicRecord: Alignment as a GraphicRecord instance.
    """
    row_dict = dict(row)
    return GraphicFeature(
        start=row["start"], end=row["end"], color=color, label=row["name"],
        strand=row["strand"]
    )


def prepare_trace_dict(row, strand):
    """Read filepath in `abi` column of `row` argument and use BioPython to
    parse the trace (abi) file and convert to a dictionary containing one
    entry per channel (4 channels 1 for each nucleotide species). Reverse trace
    values if strand is negative (-1). abi files will contain 10 values per
    base call, currently this method simplifies plotting by grabbing the
    fifth value from each channel.

    Args:
        row (DataFrame): Row of BLAST tabular output as pd.DataFrame.
        strand (int): 1 for + strand, -1 for negative strand read.

    Returns:
        dict: Dictionary with one entry per channel.
    """

    row_dict = dict(row)

    template_record = SeqIO.read(row["template_gb"], "gb")
    template_len = len(template_record.seq)

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
            if strand == -1:  # reverse trace is strand is negative 
                trace[each_channel] = trace[each_channel][::-1]
        return trace

    def pad_trace(trace):
        sub_span = [row_dict["sstart"], row_dict["send"]]
        left = [0] * min(sub_span)
        right = [0] * (template_len - max(sub_span))
        for each_channel in trace:
            trace[each_channel] = tuple(left + list(trace[each_channel]) + right)
            # assert len(trace[each_channel]) == template_len, f'padded {len(trace[each_channel])}, template: {template_len}'
        return trace

    trace = extract_traces(row_dict["abi"])
    trace = slice_trace(trace)
    trace = pad_trace(trace)

    return trace


def big_plot(features, traces, template_gb, trace_titles, read_strand, sgRNA_name,
             template_name, template_mass, cas9_species, cas9_concentration,
             sgRNA_concentration):

    fig, subs = plt.subplots(
        len(traces) + 1,  # traces plus plasmid map
        1,
        figsize=(12, 12),
        sharex=True
        # gridspec_kw={"height_ratios": [4, 1]},
    )
    # plot traces

    def plot_trace(trace_dict, axis, title="", axis_off=True, alpha=1, fill=False):
        for i, each_channel in enumerate(trace_dict):
            axis.plot(trace_dict[each_channel], color=COLORS[i], alpha=alpha)
            axis.set_title(title, loc='left')
            axis.title.set_size(20)
        if axis_off:
            axis.axis("off")
        return axis

    def get_alignment_bounds(pad=100):
        start, end = float("inf"), 0
        for each_trace in traces:
            for each_channel in each_trace:
                non_zero = np.nonzero(each_trace[each_channel])
                if non_zero[0][0] < start:
                    start = non_zero[0][0]
                if non_zero[0][-1] > end:
                    end = non_zero[0][-1]

        start, end = start - pad, end + pad
        if start < 0:
            start = 0
        return start, end

    
    def determine_zoom_coordinates(traces, read_strand):
        # want to zoom in on the end of the shortest read
        
        zoom_start, zoom_end = 0, 0
        # identify the shortest trace
        min_read_length = float('inf')
        min_read_end = None
        for each_trace in traces:
            # only look at first channel not interested in actual nucleotides
            channel = each_trace[CHANNELS[0]]
            non_zero = np.nonzero(channel)[0]
            left, right = min(non_zero), max(non_zero)
            read_len = abs(left - right)
            if read_len < min_read_length:
                min_read_length = read_len
                if read_strand == 1:
                    min_read_end = right
                else:
                    min_read_end = left
        
        return (min_read_end - 20, min_read_end + 20)
                
            
                
    def add_zoom(axis, trace_dict, target_region_start, target_region_end):
        # https://matplotlib.org/devdocs/gallery/subplots_axes_and_figures/zoom_inset_axes.html
        axins = axis.inset_axes([0.2, 0.6, 0.4, 0.4])
        axins = plot_trace(trace_dict, axins, axis_off=False, fill=True)
        
        # region to zoom in on
        x1, y1, x2, y2 = target_region_start, 0, target_region_end, axis.get_ylim()[-1]
        axins.set_xlim(x1, x2)
        axins.set_ylim(y1, y2)
        axins.set_xticklabels([])
        axins.set_yticklabels([])
        rec, vis = axis.indicate_inset_zoom(axins, edgecolor="black")
        vis = (True, True, False, False)
        
        return axis

    zoom_start, zoom_end = determine_zoom_coordinates(traces, read_strand)

    for i, each_trace in enumerate(traces):
        subs[i] = plot_trace(
            each_trace, subs[i], trace_titles[i], alpha=0.3
        )  # need way to set title (read_type)
        
        subs[i] = add_zoom(subs[i], each_trace, zoom_start, zoom_end)
    
    
    def modify_features(graphic_record):
        mod_features = []
        for f in graphic_record.features:
            if 'source' == f.label:
                continue
            if '|' in f.label:
                f.label = f.label.split('|')[0]
            mod_features.append(f)
        graphic_record.features = mod_features
        return graphic_record

    
    view_start, view_end = get_alignment_bounds()

    # finally plot the plasmid map with alignment locations
    record = SeqIO.read(template_gb, "genbank")
    graphic_record = BiopythonTranslator().translate_record(record)
    graphic_record.features += features
    graphic_record = modify_features(graphic_record)
    
    graphic_record = graphic_record.crop((int(view_start), int(view_end)))
    graphic_record.plot(ax=subs[-1], with_ruler=True, strand_in_label_threshold=4)
    
    text = f"Verification of Cas9 nickase activity with {sgRNA_name} via Sanger \
            sequencing. {template_mass} of {template_name} was treated with or without \
            {cas9_concentration} of {cas9_species} and with or without {sgRNA_concentration} \
            of sgRNA. Samples were then subjected to Sanger sequencing utilizing \
            a primer binding downstream of the {sgRNA_name} target site. Read traces \
            were converted to fasta files and then aligned back to the {template_name} \
            sequence using BLASTn."
    text = text.replace('             ', '\n')
    
    fig.text(
        0.15, 0, text, ha='left'
        )
    

    return fig, subs


def main():

    # read the concatenated BLAST file
    blast_tab = pd.read_csv(snakemake.input["blast"], sep="\t")

    target_tab = pd.read_csv(snakemake.input["target"], sep="\t")
    primer_tab = pd.read_csv(snakemake.input["primer"], sep="\t")

    assert len(target_tab) == 1
    assert len(primer_tab) == 1
    
    red = Color("red")
    colors = list(red.range_to(Color("green"), 2 + len(blast_tab) + 4))
    colors_hex = [c.hex for c in colors]
    
    target_feat = local_alignment_to_feature(target_tab.iloc[0], color=colors_hex.pop())
    primer_feat = local_alignment_to_feature(primer_tab.iloc[0], color=colors_hex.pop())

    alignment_feats = [
        blast_row_to_feature(row, label="read_type", color=colors_hex.pop()) for i, row in blast_tab.iterrows()
    ]
    alignment_feats += [target_feat, primer_feat]
    
    read_strand = primer_tab.iloc[0]['strand']

    traces = [prepare_trace_dict(row, read_strand) for i, row in blast_tab.iterrows()]

    template = set(blast_tab["template_gb"])
    assert len(template) == 1  # all should be same template
    template = list(template)[0]

    plt, subs = big_plot(
        alignment_feats, traces, template, list(blast_tab["read_type"]), 
        read_strand, sgRNA_name=snakemake.params['sgRNA'],
        template_name=snakemake.params['template_name'],
        template_mass=snakemake.params['template_mass'],
        cas9_species=snakemake.params['cas9_species'],
        cas9_concentration=snakemake.params['cas9_concentration'],
        sgRNA_concentration=snakemake.params['sgRNA_concentration']
    )

    plt.savefig(snakemake.output["png"])
    plt.savefig(snakemake.output["pdf"])


if __name__ == "__main__":
    main()
