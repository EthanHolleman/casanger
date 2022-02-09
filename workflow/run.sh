mkdir -p logs
source ~/anaconda3/etc/profile.d/conda.sh
snakemake -j 10 --cluster-config config/cluster.yml \
--cluster "sbatch -p {cluster.partition} -t {cluster.time} -N {cluster.nodes} \
-n {cluster.cpus} --mem {cluster.mem} -J {cluster.name} -o {cluster.output} \
-e {cluster.output} --mail-type ALL --mail-user {cluster.email}" \
--conda-frontend=mamba --latency-wait 30 \
--use-conda --configfile config/config.yml