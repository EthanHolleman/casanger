
rule abi_to_fastq:
    conda:
        '../envs/Py.yml'
    input:
        lambda wildcards: samples.loc[samples[config['sgRNA_name_column']] == wildcards.sgRNA][wildcards.read_type].values[0]
    output:
        'output/{run_name}/fastq/{sgRNA}/{sgRNA}_{read_type}_{read_date}.fa'
    script:'../scripts/abI2fa.py'


rule template_gb_to_fasta:
    conda:
        '../envs/Py.yml'
    output:
        'output/{run_name}/templates/{template}.fa'
    params:
        gb_path=lambda wildcards: samples.loc[samples['template_name'] == wildcards.template]['template_gb'].values[0]
    script:'../scripts/gb2fa.py'


rule make_blast_db:
    conda:
        '../envs/blast.yml'
    input:
        'output/{run_name}/templates/{template}.fa'
    output:
        expand(
            'output/{run_name}/blast-dbs/{template}/{template}.{blast_suffi}', 
            allow_missing=True, blast_suffi=BLAST_SUFFI
            )
    params:
        db_dir=lambda wildcards: f'output/{wildcards.run_name}/blast-dbs/{wildcards.template}',
        db_name=lambda wildcards: f'output/{wildcards.run_name}/blast-dbs/{wildcards.template}/{wildcards.template}',
        seq_file_name = lambda wildcards: wildcards.template

    shell:'''
    mkdir -p {params.db_dir}
    makeblastdb -in {input} -parse_seqids -title "{params.seq_file_name}" -dbtype nucl -out {params.db_name}
    '''


rule run_blast:
    conda:
        '../envs/blast.yml'
    input:
        db=expand(
            'output/{run_name}/blast-dbs/{template}/{template}.{blast_suffi}', 
            allow_missing=True, blast_suffi=BLAST_SUFFI
            ),
        read='output/{run_name}/fastq/{sgRNA}/{sgRNA}_{read_type}_{read_date}.fa'
    params:
        db_path=lambda wildcards: f'output/{wildcards.run_name}/blast-dbs/{wildcards.template}/{wildcards.template}',
        output_dir='output/blast-results'
    output:
        'output/{run_name}/blast-results/{sgRNA}/{sgRNA}_{read_type}_{read_date}.blast.{template}.tsv'
    shell:'''
    mkdir -p {params.output_dir}
    cat {input.read} | blastn -db {params.db_path} -perc_identity 0 -outfmt 6 > {output}
    '''


rule clean_blast_results:
    conda:
        '../envs/Py.yml'
    input:
        'output/{run_name}/blast-results/{sgRNA}/{sgRNA}_{read_type}_{read_date}.blast.{template}.tsv'
    output:
        'output/{run_name}/blast-results-clean/{sgRNA}/{sgRNA}_{read_type}_{read_date}.blast.{template}.clean.tsv'
    params:
        read_date=lambda wildcards: wildcards.read_date,
        read_type=lambda wildcards: wildcards.read_type,
        template_name=lambda wildcards: wildcards.template,
        template_gb=lambda wildcards: samples.loc[samples['template_name'] == wildcards.template]['template_gb'].values[0],
        sgRNA=lambda wildcards: wildcards.sgRNA,
        abi=lambda wildcards: samples.loc[samples['sgRNA'] == wildcards.sgRNA][wildcards.read_type].values[0]
    script:'../scripts/cleanBlast.py'