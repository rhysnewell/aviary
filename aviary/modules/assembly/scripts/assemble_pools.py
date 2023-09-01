import os
import subprocess

def assemble_pools(
    input_list: str,
    input_fasta: str,
    output_fasta: str,
    metabat_done: str,
    threads: int,
    log: str):
    """
    Assemble metagenome bins using unicycler
    """
        
    try:
        os.makedirs("data/final_assemblies")
    except FileExistsError:
        pass

    out_assemblies = []
    with open(input_list) as f:
        for line in f:
            if len(line.split()) == 6:
                mb_bin, long_reads, length, bases_nano, short_reads, bases_ill = line.split()
            else:
                mb_bin, long_reads, length, bases_nano = line.split()
                short_reads, bases_ill = 'none', 0
            if os.stat(long_reads).st_size == 0:
                no_long = True
            else:
                no_long = False
            long_reads = long_reads[:-5] + '.fastq.gz'
            long_reads = os.path.abspath(long_reads)
            short_reads_1 = short_reads[:-5] + '.1.fastq.gz'
            short_reads_2 = short_reads[:-5] + '.2.fastq.gz'
            short_reads_1 = os.path.abspath(short_reads_1)
            short_reads_2 = os.path.abspath(short_reads_2)
            length, bases_nano, bases_ill = float(length), float(bases_nano), float(bases_ill)
            out_assemblies.append('data/final_assemblies/%s_unicyc/assembly.fasta' % mb_bin)
            if not os.path.exists('data/final_assemblies/%s_unicyc/assembly.fasta' % mb_bin):
                with open(log, 'a') as logf:
                    if short_reads == 'none':
                        subprocess.run(
                            f"unicycler --verbosity 0 -t {threads} -l {long_reads} -o data/final_assemblies/{mb_bin}_unicyc".split(),
                            stdout=logf, stderr=subprocess.STDOUT
                            )
                    elif no_long:
                        subprocess.run(
                            f"unicycler --verbosity 0 -t {threads} -1 {short_reads_1} -2 {short_reads_2} -o data/final_assemblies/{mb_bin}_unicyc".split(),
                            stdout=logf, stderr=subprocess.STDOUT
                            )
                    else:
                        subprocess.run(
                            f"unicycler --verbosity 0 -t {threads} -1 {short_reads_1} -2 {short_reads_2} -l {long_reads} -o data/final_assemblies/{mb_bin}_unicyc".split(),
                            stdout=logf, stderr=subprocess.STDOUT
                            )

    unbinned_set = set()
    if os.path.exists(metabat_done[:-4] + "binned_contigs.unbinned"):
        with open(metabat_done[:-4] + "binned_contigs.unbinned") as f:
            for line in f:
                unbinned_set.add(line.rstrip())


    with open(output_fasta, 'w') as o:
        count = 0
        getseq = False
        with open(input_fasta) as f:
            for line in f:
                if line.startswith('>') and line.split()[0][1:] in unbinned_set:
                    getseq = True
                    o.write('>unbinned_' + str(count) + '\n')
                    count += 1
                elif line.startswith('>'):
                    getseq = False
                elif getseq:
                    o.write(line)
        for i in out_assemblies:
            if not os.path.exists(i):
                with open(i[:-14] + 'unicycler.log') as f:
                    lastline = f.readlines()[-1]
                    if lastline.startswith("Error: SPAdes failed to produce assemblies. See spades_assembly/assembly/spades.log for more info.") or \
                        lastline.startswith("Error: none of the SPAdes graphs were suitable for scaffolding in Unicycler") or \
                        lastline.startswith("Error: miniasm assembly failed"):
                        continue
            with open(i) as assembly:
                for line in assembly:
                    if line.startswith('>'):
                        count += 1
                        o.write('>unicycler_' + str(count) + '\n')
                    else:
                        o.write(line)

    if not os.path.exists(output_fasta):
        open(output_fasta, 'a').close()



if __name__ == '__main__':
    input_list = snakemake.input.list
    input_fasta = snakemake.input.fasta
    output_fasta = snakemake.output.fasta
    metabat_done = snakemake.input.metabat_done

    threads = snakemake.threads
    log = snakemake.log[0]

    with open(log, 'w') as logf: pass

    assemble_pools(input_list, input_fasta, output_fasta, metabat_done, threads, log)