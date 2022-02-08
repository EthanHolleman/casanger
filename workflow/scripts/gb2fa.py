from Bio import SeqIO


def main():

    gb = str(snakemake.params["gb_path"])
    fa = str(snakemake.output)

    record = SeqIO.read(gb, "gb")
    SeqIO.write(record, fa, "fasta")


if __name__ == "__main__":
    main()
