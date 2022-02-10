
rule draw_plots:
    conda:
        '../envs/bioDraw.yml'
    input:
        blast='output/{run_name}/concat-blast-results-clean/{sgRNA}/{sgRNA}.blast.{template}.clean.concat.tsv',
        target='output/{run_name}/local-align/targets/{sgRNA}.{template}.align.tsv',
        primer='output/{run_name}/local-align/primers/{primer_name}.{template}.align.tsv'
    output:
        png='output/{run_name}/plots/{sgRNA}/{sgRNA}.{template}.{primer_name}.sanger.blast.align.png',
        pdf='output/{run_name}/plots/{sgRNA}/{sgRNA}.{template}.{primer_name}.sanger.blast.align.pdf'
    params:
        sgRNA=lambda wildcards: wildcards.sgRNA,
        template_name=lambda wildcards: wildcards.template,
        template_mass=lambda wildcards: samples.loc[samples[config['sgRNA_name_column']] == wildcards.sgRNA]['template_mass'].values[0],
        cas9_species=lambda wildcards: samples.loc[samples[config['sgRNA_name_column']] == wildcards.sgRNA]['cas9_species'].values[0],
        cas9_concentration=lambda wildcards: samples.loc[samples[config['sgRNA_name_column']] == wildcards.sgRNA]['cas9_concentration'].values[0],
        sgRNA_concentration=lambda wildcards: samples.loc[samples[config['sgRNA_name_column']] == wildcards.sgRNA]['sgRNA_concentration'].values[0]
    script:'../scripts/drawAlignments.py'