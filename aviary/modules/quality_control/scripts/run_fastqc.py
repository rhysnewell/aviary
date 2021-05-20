import subprocess
import multiprocessing as mp


def spawn_fastqc(file):
    subprocess.Popen('fastqc -o www/fastqc/ %s' % (file), shell=True).wait()

if os.path.exists('data/short_reads.fastq.gz'):
    subprocess.Popen("fastqc -o www/fastqc/ data/short_reads.fastq.gz")
elif snakemake.config['short_reads_2'] != 'none': # paired end
    pool = mp.Pool(snakemake.threads)

    mp_results = [pool.apply_async(spawn_fastqc, args=reads)
                  for reads in
                  snakemake.config['short_reads_1'] + snakemake.config['short_reads_2']]

    for result in mp_results:
        result.get()

    pool.close()
    pool.join()
elif snakemake.config['short_reads_1'] != 'none': # interleaved
    pool = mp.Pool(snakemake.threads)

    mp_results = [pool.apply_async(spawn_fastqc, args=reads)
                  for reads in
                  snakemake.config['short_reads_1']]

    for result in mp_results:
        result.get()

    pool.close()
    pool.join()
else:
    subprocess.Popen('echo "no short reads" > www/fastqc/short_reads_fastqc.html')

subprocess.Popen('touch www/fastqc/done')
