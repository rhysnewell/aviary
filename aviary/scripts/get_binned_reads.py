import subprocess
import os
import multiprocessing as mp

def process_short_reads(file, args=""):
    subprocess.Popen('eval $(printf "seqtk seq %s %s | mfqe --fastq-read-name-lists "; '
                     'for file in data/binned_reads/*.short.list; do printf "$file "; done; '
                     'printf " --output-fastq-files "; for file in data/binned_reads/*.short.list; '
                     'do printf "${{file:0:-5}}.1.fastq.gz "; done; printf "\n")' %
                     (args, file), shell=True).wait()

if os.path.exists(snakemake.input.long_reads):
    subprocess.Popen('eval $(printf "zcat %s | mfqe --fastq-read-name-lists "; '
                     'for file in data/binned_reads/*.long.list; do printf "$file "; done; '
                     'printf " --output-fastq-files "; for file in data/binned_reads/*.long.list; '
                     'do printf "${{file:0:-5}}.fastq.gz "; done; '
                     'printf "\n")' % (snakemake.input.long_reads), shell=True).wait()

if os.path.exists('data/short_reads.fastq.gz'):
    pool = mp.Pool(snakemake.threads)
    # extract reads from filtered short reads
    # take end 1 and 2 separately
    mp_results = [pool.apply_async(process_short_reads,
                                   args=('data/short_reads.fastq.gz', n))
                                   for n in ['-1', '-2']]
    for result in mp_results:
        result.get()

    pool.close()
    pool.join()
elif snakemake.config['short_reads_2'] != 'none':
    pool = mp.Pool(snakemake.threads)
    # extract reads from all short read files
    mp_results = [pool.apply_async(process_short_reads, args=(file, ''))
                  for file in snakemake.config['short_reads_1'] + snakemake.config['short_reads_2']]

    for result in mp_results:
        result.get()

    pool.close()
    pool.join()

elif snakemake.config['short_reads_1'] != 'none':
    pool = mp.Pool(snakemake.threads)
    # extract reads from all short read files
    mp_results = [pool.apply_async(process_short_reads, args=(file, ''))
                  for file in snakemake.config['short_reads_1']]

    for result in mp_results:
        result.get()

    pool.close()
    pool.join()

subprocess.Popen('touch %s' % (snakemake.output), shell=True).wait()

