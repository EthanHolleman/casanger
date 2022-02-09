
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
    script:'../scripts/drawAlignments.py'