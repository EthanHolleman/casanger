import pandas as pd

# BLAST output format 6 default headers
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

def main():
    
    blast_results = pd.read_csv(
        str(snakemake.input), sep='\t',
        header=None
    )
    blast_results.columns = COL_NAMES
    blast_results['read_date'] = snakemake.params['read_date']
    blast_results['read_type'] = snakemake.params['read_type']
    blast_results['template_gb'] = snakemake.params['template_gb']
    blast_results['template_name'] = snakemake.params['template_name']
    blast_results['sgRNA'] = snakemake.params['sgRNA']
    blast_results['abi'] = snakemake.params['abi']
    
    
    blast_results.to_csv(
        str(snakemake.output), sep='\t', index=False
    )

if __name__ == '__main__':
    main()