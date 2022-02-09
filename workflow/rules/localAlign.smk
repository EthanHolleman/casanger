rule make_primer_fasta_file:
    conda:
        '../envs/Py.yml'
    output:
        'output/{run_name}/fasta/primers/{primer_name}.fa'
    params:
        name=lambda wildcards: wildcards.primer_name,
        seq=lambda wildcards: samples.loc[samples['primer_name'] == wildcards.primer_name]['primer_seq'].values[0]
    script:'../scripts/makeFasta.py'

    
rule make_target_site_fasta_file:
    conda:
        '../envs/Py.yml'
    output:
        'output/{run_name}/fasta/targets/{sgRNA}.target.fa'
    params:
        name=lambda wildcards: f'{wildcards.sgRNA}-target-with-PAM',
        seq=lambda wildcards: samples.loc[samples[config['sgRNA_name_column']] == wildcards.sgRNA]['target_PAM'].values[0]
    script:
        '../scripts/makeFasta.py'


rule align_local:
    conda:
        '../envs/Py.yml'
    input:
        query='output/{run_name}/fasta/primers/{primer_name}.fa',
        subject='output/{run_name}/templates/{template}.fa'
    output:
        'output/{run_name}/local-align/primers/{primer_name}.{template}.align.tsv'
    params:
        subject_frmt='fasta',
        query_frmt='fasta',
        name=lambda wildcards: f'Primer {wildcards.primer_name}'
    script:'../scripts/shortLocal.py'



rule align_target:
    conda:
        '../envs/Py.yml'
    input:
        query='output/{run_name}/fasta/targets/{sgRNA}.target.fa',
        subject='output/{run_name}/templates/{template}.fa'
    output:
        'output/{run_name}/local-align/targets/{sgRNA}.{template}.align.tsv'
    params:
        subject_frmt='fasta',
        query_frmt='fasta',
        name=lambda wildcards: f'{wildcards.sgRNA} target site'
    script:'../scripts/shortLocal.py'
    

