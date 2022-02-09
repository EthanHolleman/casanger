from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO


def main():
    print(snakemake.params["name"])
    record = SeqRecord(
        Seq(snakemake.params["seq"]),
        id=snakemake.params["name"],
        name="",
        description="",
    )
    SeqIO.write([record], str(snakemake.output), "fasta")


if __name__ == "__main__":
    main()
