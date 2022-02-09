
rule draw_plots:
    conda:
        '../envs/bioDraw.yml'
    input:
        blast='output/{run_name}/concat-blast-results-clean/{sgRNA}/{sgRNA}.blast.{template}.clean.concat.tsv',
        target='output/{run_name}/blast-results/targets/{sgRNA}.target.{template}.blast.tsv',
        primer='output/{run_name}/blast-results/primers/{primer_name}.{template}.blast.tsv'
    output:
        png='output/{run_name}/plots/{sgRNA}/{sgRNA}.{template}.{primer_name}.sanger.blast.align.png',
        pdf='output/{run_name}/plots/{sgRNA}/{sgRNA}.{template}.{primer_name}.sanger.blast.align.pdf'
    script:'../scripts/drawAlignments.py'