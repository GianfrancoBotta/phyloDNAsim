# MosaicSim
benchdir = "/cluster/work/bewi/members/gbotta/phyloDNAsim3/workflow/benchmarks/MosaicSim/wgs/bulk"
resdir = "/cluster/work/bewi/members/gbotta/phyloDNAsim3/results/benchmark/MosaicSim/wgs/bulk"

num_clones_list = snakemake.params['num_clones_list']
num_cores_list = snakemake.params['num_cores_list']
num_cells = snakemake.params['num_cells']

MosaicSim_stats_list = {}
for num_clones in num_clones_list:
    time = get_time(os.path.join(benchdir, f'clone_{num_clones}', 'MosaicSim.txt')) + get_time(os.path.join(benchdir, f'clone_{num_clones}', 'generate_genomes.txt'))
    size = get_folder_size(os.path.join(resdir, f'clone_{num_clones}'))
    num_reads = count_reads(os.path.join(resdir, f'clone_{num_clones}', 'tumor_0/samplenum_0/bulkleft.fasta'), gzipped=False) * 2
    MosaicSim_stats_list[f'MosaicSim', 2*num_clones-1] = (num_reads, time, size)

# phyloDNAsim
phyloDNAsim_stats_list = {}
for num_cores in num_cores_list:
    benchdir = f"/cluster/work/bewi/members/gbotta/phyloDNAsim3/workflow/benchmarks/phyloDNAsim/core_{num_cores}/wgs/bulk"
    resdir = f"/cluster/work/bewi/members/gbotta/phyloDNAsim3/results/benchmark/phyloDNAsim/core_{num_cores}/wgs/bulk"
    for num_clones in num_clones_list:
        num_reads = count_reads(os.path.join(resdir, f'clone_{num_clones}', 'sample_1/bulkleft.fq.gz'), gzipped=True) * 2
        time = get_time(os.path.join(benchdir, f'clone_{num_clones}', 'phyloDNAsim.txt'))
        size = get_folder_size(os.path.join(resdir, f'clone_{num_clones}'))
        phyloDNAsim_stats_list[(f'phyloDNAsim {num_cores} cores', 2*num_clones-1)] = (num_reads, time, size)