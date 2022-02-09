from Bio import pairwise2
from Bio import SeqIO


def main():

    subject = SeqIO.read(snakemake.input["subject"], snakemake.params["subject_frmt"])
    query = SeqIO.read(snakemake.input["query"], snakemake.params["query_frmt"])

    query_rc = query.reverse_complement()

    align_fwd = pairwise2.align.localxx(str(subject.seq), str(query.seq))
    align_rev = pairwise2.align.localxx(str(subject.seq), str(query_rc.seq))

    # pick the best alignment out of all
    # return as some tabular file that is ref to the subject file
