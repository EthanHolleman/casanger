import pandas as pd


def main():
    input_tables = [
        pd.read_csv(t, sep=snakemake.params["delim"]) for t in list(snakemake.input)
    ]
    concat = pd.concat(input_tables)
    concat.to_csv(str(snakemake.output), sep=snakemake.params["delim"], index=False)


if __name__ == "__main__":
    main()
