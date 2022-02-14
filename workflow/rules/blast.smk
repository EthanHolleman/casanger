
rule abi_to_fastq:
    conda:
        '../envs/Py.yml'
    input:
        lambda wildcards: samples.loc[samples[config['sgRNA_name_column']] == wildcards.sgRNA][wildcards.read_type].values[0]
    output:
        'output/{run_name}/fastq/{sgRNA}/{sgRNA}_{read_type}.fa'
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
        read='output/{run_name}/fastq/{sgRNA}/{sgRNA}_{read_type}.fa'
    params:
        db_path=lambda wildcards: f'output/{wildcards.run_name}/blast-dbs/{wildcards.template}/{wildcards.template}',
        output_dir='output/blast-results'
    output:
        'output/{run_name}/blast-results/{sgRNA}/{sgRNA}_{read_type}.blast.{template}.tsv'
    shell:'''
    mkdir -p {params.output_dir}
    cat {input.read} | blastn -db {params.db_path} -perc_identity 0 -outfmt 6 > {output}
    '''


rule clean_blast_results:
    conda:
        '../envs/Py.yml'
    input:
        'output/{run_name}/blast-results/{sgRNA}/{sgRNA}_{read_type}.blast.{template}.tsv'
    output:
        'output/{run_name}/blast-results-clean/{sgRNA}/{sgRNA}_{read_type}.blast.{template}.clean.tsv'
    params:
        read_date=lambda wildcards: samples.loc[samples[config['sgRNA_name_column']] == wildcards.sgRNA]['read_date'].values[0],
        read_type=lambda wildcards: wildcards.read_type,
        template_name=lambda wildcards: wildcards.template,
        template_gb=lambda wildcards: samples.loc[samples['template_name'] == wildcards.template]['template_gb'].values[0],
        sgRNA=lambda wildcards: wildcards.sgRNA,
        abi=lambda wildcards: samples.loc[samples[config['sgRNA_name_column']] == wildcards.sgRNA][wildcards.read_type].values[0]
    script:'../scripts/cleanBlast.py'


rule concat_blast_results:
    conda:
        '../envs/Py.yml'
    input:
        expand(
            'output/{run_name}/blast-results-clean/{sgRNA}/{sgRNA}_{read_type}.blast.{template}.clean.tsv',
            read_type=READ_TYPES, allow_missing=True
        )
    output:
        'output/{run_name}/concat-blast-results-clean/{sgRNA}/{sgRNA}.blast.{template}.clean.concat.tsv'
    params:
        delim='\t'
    script:'../scripts/concatDelim.py'


# rule make_primer_fasta_file:
#     conda:
#         '../envs/Py.yml'
#     output:
#         'output/{run_name}/fasta/primers/{primer_name}.fa'
#     params:
#         name=lambda wildcards: wildcards.primer_name,
#         seq=lambda wildcards: samples.loc[samples['primer_name'] == wildcards.primer_name]['primer_seq'].values[0]
#     script:'../scripts/makeFasta.py'


# rule blast_primer:
#     conda:
#         '../envs/blast.yml'
#     input:
#         read='output/{run_name}/fasta/primers/{primer_name}.fa',
#         db=expand(
#             'output/{run_name}/blast-dbs/{template}/{template}.{blast_suffi}', 
#             allow_missing=True, blast_suffi=BLAST_SUFFI
#             )
#     output:
#         'output/{run_name}/blast-results/primers/{primer_name}.{template}.blast.tsv'

#     params:
#         db_path=lambda wildcards: f'output/{wildcards.run_name}/blast-dbs/{wildcards.template}/{wildcards.template}',
#         output_dir='output/blast-results'
#     shell:'''
#     mkdir -p {params.output_dir}
#     cat {input.read} | blastn -db {params.db_path} -perc_identity 0 -outfmt 6 > {output}
#     '''
    
    
# rule make_target_site_fasta_file:
#     conda:
#         '../envs/Py.yml'
#     output:
#         'output/{run_name}/fasta/targets/{sgRNA}.target.fa'
#     params:
#         name=lambda wildcards: f'{wildcards.sgRNA}-target-with-PAM',
#         seq=lambda wildcards: samples.loc[samples[config['sgRNA_name_column']] == wildcards.sgRNA]['target_PAM'].values[0]
#     script:
#         '../scripts/makeFasta.py'



# rule blast_cas9_target_sites:
#     conda:
#         '../envs/blast.yml'
#     input:
#         read='output/{run_name}/fasta/targets/{sgRNA}.target.fa',
#         db=expand(
#             'output/{run_name}/blast-dbs/{template}/{template}.{blast_suffi}', 
#             allow_missing=True, blast_suffi=BLAST_SUFFI
#             )
#     output:
#         'output/{run_name}/blast-results/targets/{sgRNA}.target.{template}.blast.tsv'

#     params:
#         db_path=lambda wildcards: f'output/{wildcards.run_name}/blast-dbs/{wildcards.template}/{wildcards.template}',
#         output_dir='output/blast-results'
#     shell:'''
#     mkdir -p {params.output_dir}
#     cat {input.read} | blastn -db {params.db_path} -perc_identity 0 -dust no -expect 10 -outfmt 6 > {output}
#     '''

