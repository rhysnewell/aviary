ruleorder: skip_long_assembly > get_reads_list_ref > copy_reads > short_only
ruleorder: filter_illumina_assembly > short_only
ruleorder: filter_illumina_ref > filter_illumina_ref_interleaved > ill_copy_reads > ill_copy_reads_interleaved
ruleorder: fastqc > fastqc_long
ruleorder: combine_assemblies > combine_long_only
# ruleorder: instrain > instrain_long
ruleorder: skip_long_assembly > get_high_cov_contigs > short_only
ruleorder: skip_long_assembly > filter_illumina_assembly

onsuccess:
    print("Assembly finished, no error")

onerror:
    print("An error occurred")

onstart:
    import os
    import sys

    from snakemake.utils import logger, min_version

    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"../../scripts"))

    # minimum required snakemake version
    min_version("6.0")
    long_reads = config["long_reads"]
    fasta = config["fasta"]
    short_reads_1 = config["short_reads_1"]
    short_reads_2 = config["short_reads_2"]
    threads = config["max_threads"]
    ## pplacer deadlocks on too many threads
    pplacer_threads = min(48, int(config["pplacer_threads"]))

    if long_reads == "none" and short_reads_1 == "none":
        sys.exit("Need at least one of long_reads or short_reads_1")
    if long_reads != "none" and not os.path.exists(long_reads[0]):
        sys.exit("long_reads does not point to a file")
    if short_reads_1 != "none" and not os.path.exists(short_reads_1[0]):
        sys.exit("short_reads_1 does not point to a file")
    if short_reads_2 != "none" and not os.path.exists(short_reads_2[0]):
        sys.exit("short_reads_2 does not point to a file")
    
# Filter reads against a reference, i.e. for removing host contamination from the metagenome
rule map_reads_ref:
    input:
        fastq = config["long_reads"],
        reference_filter = config["reference_filter"]
    output:
        "data/raw_mapped_ref.bam"
    conda:
        "../../envs/coverm.yaml"
    threads:
         config["max_threads"]
    shell:
        "minimap2 -ax map-ont -t {threads} {input.reference_filter} {input.fastq} | samtools view -b  > {output}"


# Get a list of reads that don't map to genome you want to filter
rule get_umapped_reads_ref:
    input:
        "data/raw_mapped_ref.bam"
    output:
        "data/unmapped_to_ref.list"
    params:
        "no_full"
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/filter_read_list.py"

# If you don't want to filter the reads using a genome just copy them into the folder
rule copy_reads:
    input:
        fastq = config["long_reads"],
    output:
        "data/long_reads.fastq.gz"
    run:
        if input.fastq[0][-3:] == ".gz" and len(input.fastq) == 1: # Check if only one longread sample
            shell("ln -s {input.fastq} {output}")
        elif len(input.fastq) != 1:
            shell("cat {input.fastq} | gzip > {output}")
        elif len(input.fastq) > 1:
            for long_sample in input.fastq:
                if long_sample[-3:] == ".gz":
                    shell("zcat %s | gzip >> {output}" % long_sample)
                else:
                    shell("cat %s | gzip >> {output}" % long_sample)

# If no reference provided merge the short reads and copy to working directory
rule ill_copy_reads:
    input:
        fastq_1 = config["short_reads_1"],
        fastq_2 = config["short_reads_2"]
    output:
        "data/short_reads.fastq.gz"
    conda:
        "../../envs/seqtk.yaml"
    shell:
        "rename.sh prefix=SLAM in={input.fastq_1} in2={input.fastq_2} out={output} addpairnum=f"


rule ill_copy_reads_interleaved:
    input:
        fastq = config["short_reads_1"]
    output:
        "data/short_reads.fastq.gz"
    conda:
        "../../envs/seqtk.yaml"
    shell:
        "rename.sh prefix=SLAM in={input.fastq_1} out={output} addpairnum=f int=t"


# Create new read file with filtered reads
rule get_reads_list_ref:
    input:
        fastq = config["long_reads"],
        list = "data/unmapped_to_ref.list"
    output:
        "data/long_reads.fastq.gz"
    conda:
        "../../envs/seqtk.yaml"
    shell:
        "seqtk subseq {input.fastq} {input.list} | gzip > {output}"


# Assembly long reads with metaflye
rule flye_assembly:
    input:
        fastq = "data/long_reads.fastq.gz"
    output:
        fasta = "data/flye/assembly.fasta",
        graph = "data/flye/assembly_graph.gfa",
        info = "data/flye/assembly_info.txt"
    params:
        genome_size = config["meta_genome_size"]
    conda:
        "../../envs/flye.yaml"
    threads:
        config["max_threads"]
    shell:
        "flye --nano-raw {input.fastq} --meta -o data/flye -t {threads} -g {params.genome_size}"


# Polish the long reads assembly with Racon
rule polish_metagenome_racon:
    input:
        fastq = "data/long_reads.fastq.gz",
        fasta = "data/flye/assembly.fasta"
    conda:
        "../../envs/racon.yaml"
    threads:
        config["max_threads"]
    params:
        prefix = "racon",
        maxcov = 200,
        rounds = 3,
        illumina = False
    output:
        fasta = "data/assembly.pol.rac.fasta"
    script:
        "../../scripts/racon_polish.py"


### Steps if illumina data exists
rule filter_illumina_ref:
    input:
        fastq_1 = config["short_reads_1"],
        fastq_2 = config["short_reads_2"],
        reference_filter = config["reference_filter"]
    output:
        bam = "data/short_unmapped_ref.bam",
        fastq = "data/short_reads.fastq.gz"
    conda:
        "../../envs/minimap2.yaml"
    threads:
         config["max_threads"]
    shell:
        """
        minimap2 -ax sr -t {threads} {input.reference_filter} {input.fastq_1} {input.fastq_2}  |
        samtools view -b -f 12 > {output.bam} && \
        samtools bam2fq {output.bam} | gzip > {output.fastq}
        """


rule filter_illumina_ref_interleaved:
    input:
        fastq_1 = config["short_reads_1"],
        reference_filter = config["reference_filter"]
    output:
        bam = "data/short_unmapped_ref.bam",
        fastq = "data/short_reads.fastq.gz"
    conda:
        "../../envs/minimap2.yaml"
    threads:
         config["max_threads"]
    shell:
        """
        minimap2 -ax sr -t {threads} {input.reference_filter} {input.fastq_1} |
        samtools view -b -f 12 > {output.bam} && \
        samtools bam2fq {output.bam} | gzip > {output.fastq}
        """


# The racon polished long read assembly is polished again with the short reads using Pilon
rule polish_meta_pilon:
    input:
        reads = "data/short_reads.fastq.gz",
        fasta = "data/assembly.pol.rac.fasta"
    output:
        fasta = "data/assembly.pol.pil.fasta",
        bam = "data/pilon.sort.bam"
    threads:
        config["max_threads"]
    params:
        pilon_memory = config["pilon_memory"]
    conda:
        "../../envs/pilon.yaml"
    shell:
        """
        minimap2 -ax sr -t {threads} {input.fasta} {input.reads} | samtools view -b | 
        samtools sort -o {output.bam} - && \
        samtools index {output.bam} && \
        pilon -Xmx{params.pilon_memory}000m --genome {input.fasta} --frags data/pilon.sort.bam \
        --threads {threads} --output data/assembly.pol.pil --fix bases
        """


# The assembly polished with Racon and Pilon is polished again with the short reads using Racon
rule polish_meta_racon_ill:
    input:
        fastq = "data/short_reads.fastq.gz",
        fasta = "data/assembly.pol.pil.fasta"
    output:
        fasta = "data/assembly.pol.fin.fasta",
        paf = "data/racon_polishing/alignment.racon_ill.0.paf"
    threads:
        config["max_threads"]
    conda:
        "../../envs/racon.yaml"
    params:
        prefix = "racon_ill",
        maxcov = 200,
        rounds = 1,
        illumina = True
    script:
        "../../scripts/racon_polish.py"


# High coverage contigs are identified
rule get_high_cov_contigs:
    input:
        info = "data/flye/assembly_info.txt",
        fasta = "data/assembly.pol.fin.fasta",
        graph = "data/flye/assembly_graph.gfa",
        paf = "data/racon_polishing/alignment.racon_ill.0.paf"
    output:
        fasta = "data/flye_high_cov.fasta"
    params:
        min_cov_long = 20.0,
        min_cov_short = 10.0,
        short_contig_size = 200000,
        long_contig_size = 500000
    run:
        ill_cov_dict = {}
        with open(input.paf) as paf:
            for line in paf:
                query, qlen, qstart, qend, strand, ref, rlen, rstart, rend = line.split()[:9]
                ref = ref[:-6]
                if not ref in ill_cov_dict:
                    ill_cov_dict[ref] = 0.0
                ill_cov_dict[ref] += (int(rend) - int(rstart)) / int(rlen)
        count = 0
        with open(input.info) as f:
            f.readline()
            high_cov_set = set()
            short_edges = {}
            for line in f:
                if int(line.split()[1]) > params.long_contig_size:
                    high_cov_set.add(line.split()[0])
                elif int(line.split()[1]) < params.short_contig_size:
                    se1 = line.split()[6].split(',')[0]
                    if se1.startswith('-'):
                        se1 = ("edge_" + se1[1:], True)
                    else:
                        se1 = ("edge_" + se1, False)
                    se2 = line.split()[6].split(',')[-1]
                    if se2.startswith('-'):
                        se2 = ("edge_" + se2[1:], False)
                    else:
                        se2 = ("edge_" + se2, True)
                    if not se1 in short_edges:
                        short_edges[se1] = []
                    short_edges[se1].append(line.split()[0])
                    if not se2 in short_edges:
                        short_edges[se2] = []
                    short_edges[se2].append(line.split()[0])
                if float(line.split()[2]) >= params.min_cov_long or not line.split()[0] in ill_cov_dict or ill_cov_dict[line.split()[0]] <= params.min_cov_short:
                    high_cov_set.add(line.split()[0])
        with open(input.graph) as f:
            filtered_contigs = set()
            for line in f:
                if line.startswith("L"):
                    if (line.split()[1], line.split()[2] == '+') in short_edges and (line.split()[3], line.split()[4] == '-') in short_edges and not line.split()[1] == line.split()[3]:
                        for i in short_edges[(line.split()[1], line.split()[2] == '+')]:
                            filtered_contigs.add(i)
                        for i in short_edges[(line.split()[3], line.split()[4] == '-')]:
                            filtered_contigs.add(i)
        for i in filtered_contigs:
            try:
                high_cov_set.remove(i)
            except KeyError:
                pass
        with open(input.fasta) as f, open(output.fasta, 'w') as o:
            write_line = False
            for line in f:
                if line.startswith('>') and line.split()[0][1:-6] in high_cov_set:
                    write_line = True
                elif line.startswith('>'):
                    write_line = False
                if write_line:
                    o.write(line)


# Illumina reads are filtered against the nanopore assembly.
# Specifically, short reads that do not map to the high coverage long contigs are collected
rule filter_illumina_assembly:
    input:
        reference = "data/flye_high_cov.fasta",
        fastq = "data/short_reads.fastq.gz"
    output:
        bam = "data/sr_vs_long.sort.bam",
        fastq = "data/short_reads.filt.fastq.gz"
    conda:
        "../../envs/minimap2.yaml"
    threads:
         config["max_threads"]
    shell:
        """
        minimap2 -ax sr -t {threads} {input.reference} {input.fastq} |  samtools view -b |
        samtools sort -o {output.bam} - && \
        samtools index {output.bam} && \
        samtools bam2fq -f 12 {output.bam} | gzip > {output.fastq}
        """

# If unassembled long reads are provided, skip the long read assembly
rule skip_long_assembly:
    input:
        fastq = "data/short_reads.fastq.gz",
        unassembled_long = config["unassembled_long"]
    output:
        fastq = "data/short_reads.filt.fastq.gz",
        fasta = "data/flye_high_cov.fasta",
        long_reads = "data/long_reads.fastq.gz"
    shell:
        """
        ln {input.fastq} {output.fastq} && \
        touch {output.fasta} && \
        ln {input.unassembled_long} {output.long_reads}
        """


# If only short reads are provided
rule short_only:
    input:
        fastq = "data/short_reads.fastq.gz"
    output:
        fastq = "data/short_reads.filt.fastq.gz",
        fasta = "data/flye_high_cov.fasta",
        long_reads = "data/long_reads.fastq.gz"
    shell:
        """
        ln {input.fastq} {output.fastq} && \
        touch {output.fasta} && \
        touch {output.long_reads}
        """

# Short reads that did not map to the long read assembly are hybrid assembled with metaspades
# If no long reads were provided, long_reads.fastq.gz will be empty
rule spades_assembly:
    input:
        fastq = "data/short_reads.filt.fastq.gz",
        long_reads = "data/long_reads.fastq.gz"
    output:
        fasta = "data/spades_assembly.fasta"
    threads:
         config["max_threads"]
    params:
         max_memory = config["max_memory"]
    conda:
        "../../envs/spades.yaml"
    shell:
        """
        minimumsize=500000 && \
        actualsize=$(stat -c%s data/short_reads.filt.fastq.gz) && \
        if [ $actualsize -ge $minimumsize ]
        then
            spades.py --memory {params.max_memory} --meta --nanopore {input.long_reads} --12 {input.fastq} \
            -o data/spades_assembly -t {threads} -k 21,33,55,81,99,127 && \
            ln data/spades_assembly/scaffolds.fasta data/spades_assembly.fasta
        else
            touch {output.fasta}
        fi 
        """



# Short reads that did not map to the long read assembly are hybrid assembled with metaspades
# If no long reads were provided, long_reads.fastq.gz will be empty
rule spades_assembly_short:
    output:
        fasta = "data/final_contigs.fasta"
    threads:
         config["max_threads"]
    params:
         max_memory = config["max_memory"]
    conda:
        "../../envs/spades.yaml"
    script:
        "../../scripts/spades_assembly_short.py"


# Short reads are mapped to the spades assembly and jgi_summarize_bam_contig_depths from metabat
# used to calculate the mean coverage of the contigs. The coverage and assembly are then used to bin with metabat.
rule spades_assembly_coverage:
    input:
         fastq = "data/short_reads.filt.fastq.gz",
         fasta = "data/spades_assembly.fasta"
    output:
         assembly_cov = "data/short_read_assembly.cov",
         short_vs_mega = "data/short_vs_mega.bam"
    conda:
         "../../envs/coverm.yaml"
    threads:
         config["max_threads"]
    shell:
        """
        coverm contig -m metabat -t {threads} -r {input.fasta} -i {input.fastq} --bam-file-cache-directory cached_bams/ > data/short_read_assembly.cov;
        ln -s cached_bams/*.bam ./short_vs_mega.bam
        """

rule metabat_binning_short:
    input:
         assembly_cov = "data/short_read_assembly.cov",
         fasta = "data/spades_assembly.fasta"
    output:
         metabat_done = "data/metabat_bins/done"
    conda:
         "../../envs/metabat2.yaml"
    threads:
         config["max_threads"]
    shell:
         """
         mkdir -p data/metabat_bins && \
         metabat --seed 89 --unbinned -m 1500 -l -i {input.fasta} -a {input.assembly_cov} \
         -o data/metabat_bins/binned_contigs && \
         touch data/metabat_bins/done
         """

# Long reads are mapped to the spades assembly
rule map_long_mega:
    input:
        fastq = "data/long_reads.fastq.gz",
        fasta = "data/spades_assembly.fasta"
    output:
        bam = "data/long_vs_mega.bam"
    threads:
        config["max_threads"]
    conda:
        "../../envs/minimap2.yaml"
    shell:
        """
        minimap2 -t {threads} -ax map-ont -a {input.fasta} {input.fastq} |  samtools view -b |
        samtools sort -o {output.bam} - && \
        samtools index {output.bam}
        """

# Long and short reads that mapped to the spades assembly are pooled (binned) together
rule pool_reads:
    input:
        long_bam = "data/long_vs_mega.bam",
        short_bam = "data/short_vs_mega.bam",
        metabat_done = "data/metabat_bins/done"
    output:
        list = "data/list_of_lists.txt"
    conda:
        "../../envs/pysam.yaml"
    script:
        "../../scripts/pool_reads.py"

# Binned read lists are processed to extract the reads associated with each bin
rule get_read_pools:
    input:
        long_reads = "data/long_reads.fastq.gz",
        short_reads = "data/short_reads.fastq.gz",
        list = "data/list_of_lists.txt"
    output:
        "data/binned_reads/done"
    conda:
         "../../envs/mfqe.yaml"
    shell:
         'eval $(printf "zcat {input.long_reads} | mfqe --fastq-read-name-lists "; for file in data/binned_reads/*.long.list; do printf "$file "; done;'\
         ' printf " --output-fastq-files "; for file in data/binned_reads/*.long.list; do printf "${{file:0:-5}}.fastq.gz "; done; printf "\n") && ' \
         'eval $(printf "seqtk seq -1 {input.short_reads} | mfqe --fastq-read-name-lists "; for file in data/binned_reads/*.short.list; do printf "$file "; done;'\
         ' printf " --output-fastq-files "; for file in data/binned_reads/*.short.list; do printf "${{file:0:-5}}.1.fastq.gz "; done; printf "\n") && ' \
         'eval $(printf "seqtk seq -2 {input.short_reads} | mfqe --fastq-read-name-lists "; for file in data/binned_reads/*.short.list; do printf "$file "; done;'\
         ' printf " --output-fastq-files "; for file in data/binned_reads/*.short.list; do printf "${{file:0:-5}}.2.fastq.gz "; done; printf "\n") && touch {output} || touch {output}'

# Short and long reads for each bin are hybrid assembled with Unicycler
rule assemble_pools:
    input:
        fastq = "data/binned_reads/done",
        list = "data/list_of_lists.txt",
        fasta = "data/spades_assembly.fasta",
        metabat_done = "data/metabat_bins/done"
    threads:
        config["max_threads"]
    output:
        fasta = "data/unicycler_combined.fa"
    conda:
        "../../envs/final_assembly.yaml"
    script:
        "scripts/assemble_pools.py"

# The long read high coverage assembly from flye and hybrid assembly from unicycler are combined.
# Long and short reads are mapped to this combined assembly.
rule combine_assemblies:
    input:
        short_reads = "data/short_reads.fastq.gz",
        long_reads = "data/long_reads.fastq.gz",
        unicyc_fasta = "data/unicycler_combined.fa",
        flye_fasta = "data/flye_high_cov.fasta"
    output:
        short_bam = "data/final_short.sort.bam",
        fasta = "data/final_contigs.fasta",
        long_bam = "data/final_long.sort.bam"
    conda:
        "../../envs/coverm.yaml"
    priority: 1
    threads:
        config["max_threads"]
    shell:
        """
        cat {input.flye_fasta} {input.unicyc_fasta} > {output.fasta} && \
        minimap2 -t {threads} -ax map-ont -a {output.fasta} {input.long_reads} |  samtools view -b | 
        samtools sort -o {output.long_bam} - && samtools index {output.long_bam} && \
        minimap2 -ax sr -t {threads} {output.fasta} {input.short_reads} |  samtools view -b |
        samtools sort -o {output.short_bam} - && \
        samtools index {output.short_bam}
        """

rule combine_long_only:
    input:
        long_reads = "data/long_reads.fastq.gz",
        fasta = "data/assembly.pol.rac.fasta"
    output:
        fasta = "data/final_contigs.fasta",
        bam = "data/final_long.sort.bam"
    priority: 1
    conda:
        "../../envs/coverm.yaml"
    threads:
        config["max_threads"]
    shell:
        """
        minimap2 -t {threads} -ax map-ont -a {input.fasta} {input.long_reads} |  samtools view -b | 
        samtools sort -o {output.bam} - && \
        samtools index {output.bam} && \
        ln {input.fasta} {output.fasta}
        """

rule fastqc:
    input:
        "data/short_reads.fastq.gz"
    output:
        "www/short_reads_fastqc.html"
    conda:
        "../../envs/fastqc.yaml"
    threads:
        config["max_threads"]
    shell:
        "fastqc -o www {input}"

rule fastqc_long:
    output:
        "www/short_reads_fastqc.html"
    shell:
        'echo "no short reads" > {output}'


rule nanoplot:
    input:
        "data/long_reads.fastq.gz"
    output:
        "www/nanoplot/longReadsNanoPlot-report.html"
    conda:
        "../../envs/nanoplot.yaml"
    threads:
        config["max_threads"]
    shell:
        "NanoPlot -o www/nanoplot -p longReads --fastq {input}"