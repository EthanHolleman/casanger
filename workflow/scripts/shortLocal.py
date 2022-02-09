from Bio import pairwise2
from Bio import SeqIO
import pandas as pd


def main():

    subject = SeqIO.read(snakemake.input["subject"], snakemake.params["subject_frmt"])
    query = SeqIO.read(snakemake.input["query"], snakemake.params["query_frmt"])

    query_rc = query.reverse_complement()

    # https://biopython.org/docs/1.76/api/Bio.pairwise2.html
    align_fwd = pairwise2.align.localms(
        str(subject.seq), str(query.seq), 2, -1, -0.5, -0.1
    )
    align_rev = pairwise2.align.localms(
        str(subject.seq), str(query_rc.seq), 2, -1, -0.5, -0.1
    )

    max_fwd = max(align_fwd, key=lambda a: a.score)
    max_rev = max(align_rev, key=lambda a: a.score)

    if max_fwd.score >= max_rev.score:
        best_alignment = max_fwd
        strand = 1
    else:
        best_alignment = max_rev
        strand = -1

    df = pd.DataFrame(
        [
            {
                "name": snakemake.params["name"],
                "strand": strand,
                "start": best_alignment.start,
                "end": best_alignment.end,
            }
        ]
    )

    df.to_csv(str(snakemake.output), sep="\t", index=False)


if __name__ == "__main__":
    main()

    # pick the best alignment out of all
    # return as some tabular file that is ref to the subject file
