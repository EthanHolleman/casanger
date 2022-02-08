from Bio import SeqIO

def main():
    
    abi = str(snakemake.input)
    fasta = str(snakemake.output)
    
    record = SeqIO.read(abi, 'abi')
    SeqIO.write(record, fasta, 'fasta')

if __name__ == '__main__':
    main()