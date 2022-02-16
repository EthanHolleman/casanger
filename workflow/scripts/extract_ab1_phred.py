from Bio import SeqIO
import pandas as pd


def extract_phred(ab1_path, sgRNA, read_type):
    record = SeqIO.read(ab1_path, "abi")
    return pd.DataFrame(
        [
            {"phred": phred, "position": i + 1, "sgRNA": sgRNA, "read_type": read_type}
            for i, phred in enumerate(record.letter_annotations["phred_quality"])
        ]
    )


def main():

    abi_file = str(snakemake.input)
    phred_table = extract_phred(
        abi_file, snakemake.params['sgRNA'], snakemake.params['read_type']
        )
    phred_table.to_csv(str(snakemake.output), sep="\t", index=False)


if __name__ == "__main__":
    main()
