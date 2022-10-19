ruleorder: filter_illumina_assembly > short_only
# ruleorder: fastqc > fastqc_long
ruleorder: skip_unicycler_with_qc > skip_unicycler > combine_assemblies > combine_long_only
ruleorder: skip_long_assembly > get_high_cov_contigs > short_only
ruleorder: skip_long_assembly > filter_illumina_assembly
ruleorder: filter_illumina_ref > no_ref_filter
ruleorder: skip_unicycler_with_qc > skip_unicycler > combine_assemblies > combine_long_only > assemble_short_reads
ruleorder: skip_unicycler_with_qc > skip_unicycler > complete_assembly_with_qc > complete_assembly
ruleorder: skip_unicycler_with_qc > skip_unicycler > combine_assemblies > move_spades_assembly

# onsuccess:
#     print("Assembly finished, no error")
#
# onerror:
#     print("An error occurred")

onstart:
    import os
    import sys

    from snakemake.utils import logger, min_version
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"../../scripts"))
    sys.path.append(os.path.join(os.path.dirname(os.path.abspath(workflow.snakefile)),"scripts"))

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
    group: 'assembly'
    output:
        temp("data/raw_mapped_ref.bam")
    conda:
        "../../envs/coverm.yaml"
    benchmark:
        "benchmarks/map_reads_ref.benchmark.txt"
    threads:
         config["max_threads"]
    shell:
        "minimap2 -ax map-ont --split-prefix=tmp -t {threads} {input.reference_filter} {input.fastq} | samtools view -@ {threads} -b > {output}"


# Get a list of reads that don't map to genome you want to filter
rule get_umapped_reads_ref:
    input:
        "data/raw_mapped_ref.bam"
    group: 'assembly'
    output:
        temp("data/unmapped_to_ref.list")
    params:
        "no_full"
    conda:
        "envs/pysam.yaml"
    benchmark:
        "benchmarks/get_unmapped_reads_ref.benchmark.txt"
    script:
        "scripts/filter_read_list.py"


# Create new read file with filtered reads
rule get_reads_list_ref:
    input:
        fastq = config["long_reads"],
        list = "data/unmapped_to_ref.list"
    group: 'assembly'
    output:
        temp("data/long_reads.fastq.gz")
    threads:
        config['max_threads']
    conda:
        "envs/seqtk.yaml"
    benchmark:
        "benchmarks/get_reads_list_ref.benchmark.txt"
    shell:
        "seqtk subseq {input.fastq} {input.list} | pigz -p {threads} > {output}"

# if no reference filter output this done file just to keep the DAG happy
rule no_ref_filter:
    group: 'assembly'
    output:
        filtered = temp("data/short_filter.done")
    shell:
        "touch {output.filtered}"

# Assembly long reads with metaflye
rule flye_assembly:
    input:
        fastq = "data/long_reads.fastq.gz"
    group: 'assembly'
    output:
        fasta = "data/flye/assembly.fasta",
        graph = "data/flye/assembly_graph.gfa",
        info = "data/flye/assembly_info.txt",
        junk1 = temp(directory("data/flye/00-assembly/")),
        junk2 = temp(directory("data/flye/10-consensus/")),
        junk3 = temp(directory("data/flye/20-repeat/")),
        junk4 = temp(directory("data/flye/30-contigger/")),
        junk5 = temp(directory("data/flye/40-polishing/"))
    params:
        long_read_type = config["long_read_type"]
    resources:
        mem_mb=int(config["max_memory"])*1024
    conda:
        "envs/flye.yaml"
    benchmark:
        "benchmarks/flye_assembly.benchmark.txt"
    threads:
        config["max_threads"]
    script:
        "scripts/run_flye.py"


# Polish the long reads assembly with Racon
rule polish_metagenome_racon:
    input:
        fastq = "data/long_reads.fastq.gz",
        fasta = "data/flye/assembly.fasta",
    conda:
        "envs/racon.yaml"
    threads:
        config["max_threads"]
    params:
        prefix = "racon",
        maxcov = 200,
        rounds = 3,
        illumina = False,
        coassemble = config["coassemble"]
    resources:
        mem_mb=int(config["max_memory"])*1024
    group: 'assembly'
    output:
        fasta = "data/assembly.pol.rac.fasta"
    benchmark:
        "benchmarks/polish_metagenome_racon.benchmark.txt"
    script:
        "scripts/racon_polish.py"


### Filter illumina reads against provided reference
rule filter_illumina_ref:
    input:
        reference_filter = config["reference_filter"]
    group: 'assembly'
    output:
        bam = temp("data/short_unmapped_ref.bam"),
        fastq = temp("data/short_reads.fastq.gz"),
        filtered = temp("data/short_filter.done")
    params:
        coassemble = config["coassemble"]
    conda:
        "../../envs/minimap2.yaml"
    threads:
         config["max_threads"]
    benchmark:
        "benchmarks/filter_illumina_ref.benchmark.txt"
    script:
        "scripts/filter_illumina_reference.py"


# Generate BAM file for pilon, discard unmapped reads
rule generate_pilon_sort:
    input:
        fastq = config['short_reads_1'],
        filtered = "data/short_filter.done",
        fasta = "data/assembly.pol.rac.fasta"
    params:
        coassemble = config["coassemble"]
    group: 'assembly'
    output:
        bam = temp("data/pilon.sort.bam"),
        bai = temp("data/pilon.sort.bam.bai")
    threads:
        config["max_threads"]
    resources:
        mem_mb=int(config["max_memory"])*1024
    conda:
        "envs/pilon.yaml"
    benchmark:
        "benchmarks/generate_pilon_sort.benchmark.txt"
    script:
        "scripts/generate_pilon_sort.py"

# The racon polished long read assembly is polished again with the short reads using Pilon
rule polish_meta_pilon:
    input:
        fasta = "data/assembly.pol.rac.fasta",
        bam = "data/pilon.sort.bam",
        bai = "data/pilon.sort.bam.bai"
    group: 'assembly'
    output:
        fasta = "data/assembly.pol.pil.fasta"
    resources:
        mem_mb=int(config["max_memory"])*512
    # threads: # Threads no longer supported for pilon
    #     config["max_threads"]
    params:
        pilon_memory = int(config["max_memory"])*512
    conda:
        "envs/pilon.yaml"
    benchmark:
        "benchmarks/polish_meta_pilon.benchmark.txt"
    shell:
        """
        pilon -Xmx{params.pilon_memory}m --genome {input.fasta} --frags data/pilon.sort.bam --output data/assembly.pol.pil --fix bases >data/pilon.err
        """


# The assembly polished with Racon and Pilon is polished again with the short reads using Racon
rule polish_meta_racon_ill:
    input:
        fastq = config['short_reads_1'], # check short reads are here
        fasta = "data/assembly.pol.pil.fasta"
    group: 'assembly'
    output:
        fasta = "data/assembly.pol.fin.fasta",
        paf = temp("data/racon_polishing/alignment.racon_ill.0.paf")
    resources:
        mem_mb=int(config["max_memory"])*1024
    threads:
        config["max_threads"]
    conda:
        "envs/racon.yaml"
    params:
        prefix = "racon_ill",
        maxcov = 200,
        rounds = 1,
        illumina = True,
        coassemble = config["coassemble"]
    benchmark:
        "benchmarks/polish_meta_racon_ill.benchmark.txt"
    script:
        "scripts/racon_polish.py"


# High coverage contigs are identified
rule get_high_cov_contigs:
    input:
        info = "data/flye/assembly_info.txt",
        fasta = "data/assembly.pol.fin.fasta",
        graph = "data/flye/assembly_graph.gfa",
        paf = "data/racon_polishing/alignment.racon_ill.0.paf"
    group: 'assembly'
    output:
        fasta = "data/flye_high_cov.fasta"
    benchmark:
        "benchmarks/get_high_cov_contigs.benchmark.txt"
    params:
        min_cov_long = int(config["min_cov_long"]),
        min_cov_short = int(config["min_cov_short"]),
        exclude_contig_cov = int(config["exclude_contig_cov"]),
        exclude_contig_size = int(config["exclude_contig_size"]),
        long_contig_size = int(config["long_contig_size"])
    run:
        ill_cov_dict = {}
        # populate illumina coverage dictionary using PAF
        with open(input.paf) as paf:
            for line in paf:
                query, qlen, qstart, qend, strand, ref, rlen, rstart, rend = line.split()[:9]
                ref = ref[:-6]
                if not ref in ill_cov_dict:
                    ill_cov_dict[ref] = 0.0
                ill_cov_dict[ref] += (int(rend) - int(rstart)) / int(rlen)
        count = 0
        high_cov_set = set()
        short_edges = {}
        with open(input.info) as f:
            # Based on the flye assembly_info.txt, retrieve contigs > long_contig_size
            f.readline()
            for line in f:
                contig_name, contig_length, long_coverage, circular, repeat, multiplicity, alt_group, graph_path = line.split()
                is_circular = circular == "Y"
                if int(contig_length) >= params.long_contig_size or is_circular:
                    high_cov_set.add(contig_name)
                elif int(contig_length) < params.long_contig_size:
                    # Take first and last edge in info file for this contig
                    se1 = graph_path.split(',')[0]
                    # Filter the '-' sign from edge name
                    if se1.startswith('-'):
                        se1 = ("edge_" + se1[1:], True)
                    else:
                        se1 = ("edge_" + se1, False)
                    se2 = graph_path.split(',')[-1]
                    # Filter the '-' sign from edge name
                    if se2.startswith('-'):
                        se2 = ("edge_" + se2[1:], False)
                    else:
                        se2 = ("edge_" + se2, True)

                    # Append contig to associated edges in short_edges dict
                    if not se1 in short_edges:
                        short_edges[se1] = []
                    short_edges[se1].append(contig_name)
                    if not se2 in short_edges:
                        short_edges[se2] = []
                    short_edges[se2].append(contig_name)
                # Filtering:
                # 1:
                # exclude contigs with long read coverage less than `exclude_contig_cov` and
                # shorter than `exclude_contig_size`
                # 2:
                # include if a contig is covered by >= min_cov_long long reads place in high cov set
                # Also place in high coverage set if contig was not covered by illumina reads
                # Also place in high coverage set if illumina coverage was <= min_cov_short
                if not (float(long_coverage) <= params.exclude_contig_cov and int(contig_length) <= params.exclude_contig_size):
                    if float(long_coverage) >= params.min_cov_long or not contig_name in ill_cov_dict \
                            or ill_cov_dict[contig_name] <= params.min_cov_short:
                        high_cov_set.add(contig_name)

        filtered_contigs = set()
        # Populate filtered contigs with short contigs that were connected to two edges in short_edges dict in the correct
        # orientation: edge_1 + to edge_2 - and edge_1 != edge_2. So the start of a sequences complements links to
        #               the end of a sequences reverse complement.
        with open(input.graph) as f:
            for line in f:
                if line.startswith("L"):
                    if (line.split()[1], line.split()[2] == '+') in short_edges \
                            and (line.split()[3], line.split()[4] == '-') in short_edges \
                            and not line.split()[1] == line.split()[3]: # and not line.split()[5] == '0M':
                        for i in short_edges[(line.split()[1], line.split()[2] == '+')]:
                            filtered_contigs.add(i)
                        for i in short_edges[(line.split()[3], line.split()[4] == '-')]:
                            filtered_contigs.add(i)

        # Remove contigs that were filtered from the high coverage set
        for i in filtered_contigs:
            try:
                high_cov_set.remove(i)
            except KeyError:
                pass

        # Write high coverage contigs
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
        fastq = config['short_reads_1'], # check short reads were supplied
        fasta = "data/flye_high_cov.fasta"
    group: 'assembly'
    output:
        bam = temp("data/sr_vs_long.sort.bam"),
        bai = temp("data/sr_vs_long.sort.bam.bai"),
        fastq = temp("data/short_reads.filt.fastq.gz")
    params:
        coassemble = config["coassemble"]
    conda:
        "../../envs/minimap2.yaml"
    threads:
         config["max_threads"]
    benchmark:
        "benchmarks/filter_illumina_assembly.benchmark.txt"
    script:
        "scripts/filter_illumina_assembly.py"


# If unassembled long reads are provided, skip the long read assembly
rule skip_long_assembly:
    input:
        unassembled_long = config["unassembled_long"]
    group: 'assembly'
    output:
        # fastq = "data/short_reads.filt.fastq.gz",
        fasta = "data/flye_high_cov.fasta",
        long_reads = temp("data/long_reads.fastq.gz")
    shell:
        """
        touch {output.fasta} && \
        ln {input.unassembled_long} {output.long_reads}
        """


# If only short reads are provided
rule short_only:
    input:
        fastq = config["short_reads_1"]
    group: 'assembly'
    output:
        fasta = "data/flye_high_cov.fasta",
        # long_reads = temp("data/long_reads.fastq.gz")
    shell:
        """
        touch {output.fasta}
        """
        # touch {output.long_reads}

# Short reads that did not map to the long read assembly are hybrid assembled with metaspades
# If no long reads were provided, long_reads.fastq.gz will be empty
rule spades_assembly:
    input:
        fastq = "data/short_reads.filt.fastq.gz",
        long_reads = "data/long_reads.fastq.gz"
    group: 'assembly'
    output:
        fasta = "data/spades_assembly.fasta",
        spades_folder = temp(directory("data/spades_assembly/"))
    threads:
        config["max_threads"]
    resources:
        mem_mb=int(config["max_memory"])*1024
    params:
        max_memory = config["max_memory"],
        long_read_type = config["long_read_type"],
        kmer_sizes = " ".join(config["kmer_sizes"]),
        tmpdir = config["tmpdir"]
    conda:
        "envs/spades.yaml"
    benchmark:
        "benchmarks/spades_assembly.benchmark.txt"
    shell:
        """
        rm -rf data/spades_assembly/tmp; 
        minimumsize=500000;
        actualsize=$(stat -c%s data/short_reads.filt.fastq.gz);
        if [ -d "data/spades_assembly/" ]
        then
            spades.py --restart-from last --memory {params.max_memory} -t {threads} -o data/spades_assembly -k {params.kmer_sizes} --tmp-dir {params.tmpdir} && \
            cp data/spades_assembly/scaffolds.fasta data/spades_assembly.fasta
        elif [ $actualsize -ge $minimumsize ]
        then
            if [ {params.long_read_type} = "ont" ] || [ {params.long_read_type} = "ont_hq" ]
            then
                spades.py --checkpoints all --memory {params.max_memory} --meta --nanopore {input.long_reads} --12 {input.fastq} \
                -o data/spades_assembly -t {threads}  -k {params.kmer_sizes} --tmp-dir {params.tmpdir} 2>data/spades.err && \
                cp data/spades_assembly/scaffolds.fasta data/spades_assembly.fasta
            else
                spades.py --checkpoints all --memory {params.max_memory} --meta --pacbio {input.long_reads} --12 {input.fastq} \
                -o data/spades_assembly -t {threads}  -k {params.kmer_sizes} --tmp-dir {params.tmpdir} 2>data/spades.err && \
                cp data/spades_assembly/scaffolds.fasta data/spades_assembly.fasta
            fi
        else
            touch {output.fasta}
        fi 
        """


# Perform shrot read assembly only with no other steps
rule assemble_short_reads:
    input:
        fastq = config["short_reads_1"]
    group: 'assembly'
    output:
        fasta = "data/short_read_assembly/scaffolds.fasta",
        # We cannot mark the output_folder as temp as then it gets deleted,
        # causing the "TypeError: '>' not supported between instances of
        # 'TBDString' and 'int'" error
        output_folder = directory("data/short_read_assembly/")
        # reads1 = temporary("data/short_reads.1.fastq.gz"),
        # reads2 = temporary("data/short_reads.2.fastq.gz")
    threads:
         config["max_threads"]
    params:
         max_memory = config["max_memory"],
         kmer_sizes = config["kmer_sizes"],
         use_megahit = config["use_megahit"],
         coassemble = config["coassemble"],
         tmpdir = config["tmpdir"],
         final_assembly = True
    resources:
        mem_mb=int(config["max_memory"])*1024
    conda:
        "envs/spades.yaml"
    benchmark:
        "benchmarks/short_read_assembly_short.benchmark.txt"
    script:
        "scripts/assemble_short_reads.py"


rule move_spades_assembly:
    input:
        assembly = "data/short_read_assembly/scaffolds.fasta"
    group: 'assembly'
    output:
        out = "data/final_contigs.fasta"
    shell:
        "cp {input.assembly} {output.out} && rm -rf data/short_read_assembly"


# Short reads are mapped to the spades assembly and jgi_summarize_bam_contig_depths from metabat
# used to calculate the mean coverage of the contigs. The coverage and assembly are then used to bin with metabat.
rule spades_assembly_coverage:
    input:
         fastq = "data/short_reads.filt.fastq.gz",
         fasta = "data/spades_assembly.fasta"
    group: 'assembly'
    output:
         assembly_cov = temp("data/short_read_assembly.cov"),
         bam = temp("data/short_vs_mega.bam"),
         bai = temp("data/short_vs_mega.bam.bai")
    resources:
        mem_mb=int(config["max_memory"])*1024
    params:
         tmpdir = config["tmpdir"]
    conda:
         "../../envs/coverm.yaml"
    threads:
         config["max_threads"]
    benchmark:
        "benchmarks/spades_assembly_coverage.benchmark.txt"
    shell:
        """
        TMPDIR={params.tmpdir} coverm contig -m metabat -t {threads} -r {input.fasta} --interleaved {input.fastq} --bam-file-cache-directory data/cached_bams/ > {output.assembly_cov};
        mv data/cached_bams/*.bam {output.bam} && samtools index -@ {threads} {output.bam}
        """

rule metabat_binning_short:
    input:
         assembly_cov = "data/short_read_assembly.cov",
         fasta = "data/spades_assembly.fasta"
    group: 'assembly'
    output:
         metabat_done = "data/metabat_bins/done"
    conda:
         "../binning/envs/metabat2.yaml"
    threads:
         config["max_threads"]
    benchmark:
        "benchmarks/metabat_binning_short.benchmark.txt"
    shell:
         """
         mkdir -p data/metabat_bins && \
         metabat --seed 89 --unbinned -m 1500 -l -i {input.fasta} -t {threads} -a {input.assembly_cov} \
         -o data/metabat_bins/binned_contigs && \
         touch data/metabat_bins/done
         """

# Long reads are mapped to the spades assembly
rule map_long_mega:
    input:
        fastq = "data/long_reads.fastq.gz",
        fasta = "data/spades_assembly.fasta"
    group: 'assembly'
    output:
        bam = temp("data/long_vs_mega.bam"),
        bai = temp("data/long_vs_mega.bam.bai")
    resources:
        mem_mb=int(config["max_memory"])*1024
    threads:
        config["max_threads"]
    conda:
        "../../envs/minimap2.yaml"
    benchmark:
        "benchmarks/map_long_mega.benchmark.txt"
    shell:
        """
        minimap2 -t {threads} --split-prefix=tmp -ax map-ont -a {input.fasta} {input.fastq} |  samtools view -@ {threads} -b |
        samtools sort -@ {threads} -o {output.bam} - && \
        samtools index -@ {threads} {output.bam}
        """

# Long and short reads that mapped to the spades assembly are pooled (binned) together
rule pool_reads:
    input:
        long_bam = "data/long_vs_mega.bam",
        long_bai = "data/long_vs_mega.bam.bai",
        short_bam = "data/short_vs_mega.bam",
        short_bai = "data/short_vs_mega.bam.bai",
        metabat_done = "data/metabat_bins/done",
    group: 'assembly'
    output:
        list = "data/list_of_lists.txt"
    conda:
        "envs/pysam.yaml"
    benchmark:
        "benchmarks/pool_reads.benchmark.txt"
    script:
        "scripts/pool_reads.py"

# Binned read lists are processed to extract the reads associated with each bin
rule get_read_pools:
    input:
        list = "data/list_of_lists.txt"
    group: 'assembly'
    output:
        "data/binned_reads/done"
    conda:
         "envs/mfqe.yaml"
    threads:
         config['max_threads']
    benchmark:
        "benchmarks/get_read_pools.benchmark.txt"
    script:
         'scripts/get_binned_reads.py'

# Short and long reads for each bin are hybrid assembled with Unicycler
rule assemble_pools:
    input:
        fastq = "data/binned_reads/done",
        list = "data/list_of_lists.txt",
        fasta = "data/spades_assembly.fasta",
        metabat_done = "data/metabat_bins/done"
    threads:
        config["max_threads"]
    group: 'assembly'
    resources:
        mem_mb=int(config["max_memory"])*1024
    output:
        fasta = "data/unicycler_combined.fa"
    conda:
        "envs/final_assembly.yaml"
    benchmark:
        "benchmarks/assemble_pools.benchmark.txt"
    script:
        "scripts/assemble_pools.py"

# The long read high coverage assembly from flye and hybrid assembly from unicycler are combined.
# Long and short reads are mapped to this combined assembly.
rule combine_assemblies:
    input:
        short_fasta = "data/unicycler_combined.fa",
        flye_fasta = "data/flye_high_cov.fasta"
    group: 'assembly'
    output:
        fasta = "data/final_contigs.fasta",
    priority: 1
    threads:
        config["max_threads"]
    script:
        "scripts/combine_assemblies.py"



rule combine_long_only:
    input:
        long_reads = "data/long_reads.fastq.gz",
        fasta = "data/assembly.pol.rac.fasta"
    group: 'assembly'
    output:
        fasta = "data/final_contigs.fasta",
        # long_bam = "data/final_long.sort.bam"
    priority: 1
    conda:
        "../../envs/minimap2.yaml"
    threads:
        config["max_threads"]
    script:
        "scripts/combine_assemblies.py"


rule skip_unicycler:
    input:
        short_fasta = "data/spades_assembly.fasta",
        flye_fasta = "data/flye_high_cov.fasta"
    group: 'assembly'
    output:
        fasta = "data/final_contigs.fasta",
        final_link = 'assembly/final_contigs.fasta',
        sizes = "www/assembly_stats.txt",
        unicycler_skipped = temp("data/unicycler_skipped")
    priority: 1
    threads:
        config["max_threads"]
    shell:
        'cat {input.short_fasta} {input.flye_fasta} > {output.fasta}; '
        'touch data/unicycler_skipped; '
        'mkdir -p www/; '
        'stats.sh {output.fasta} > {output.sizes}; '
        'mkdir -p assembly; '
        'cd assembly; '
        'ln -s ../data/final_contigs.fasta ./; '

rule skip_unicycler_with_qc:
    input:
        short_fasta = "data/spades_assembly.fasta",
        flye_fasta = "data/flye_high_cov.fasta",
        qc_done = 'data/qc_done'
    group: 'assembly'
    output:
        fasta = "data/final_contigs.fasta",
        final_link = 'assembly/final_contigs.fasta',
        sizes = "www/assembly_stats.txt",
        unicycler_skipped = temp("data/unicycler_skipped")
    priority: 1
    threads:
        config["max_threads"]
    shell:
        'cat {input.short_fasta} {input.flye_fasta} > {output.fasta}; '
        'touch data/unicycler_skipped; '
        'mkdir -p www/; '
        'stats.sh {output.fasta} > {output.sizes}; '
        'mkdir -p assembly; '
        'cd assembly; '
        'ln -s ../data/final_contigs.fasta ./; '

rule complete_assembly:
    input:
        fasta = 'data/final_contigs.fasta'
    group: 'assembly'
    output:
        final_link = 'assembly/final_contigs.fasta',
        sizes = "www/assembly_stats.txt"
    shell:
        'mkdir -p www/; '
        'stats.sh {input.fasta} > {output.sizes}; '
        'mkdir -p assembly; '
        'cd assembly; '
        'ln -s ../data/final_contigs.fasta ./; '
        'rm -rf data/racon_polishing; '

rule complete_assembly_with_qc:
    input:
        fasta = 'data/final_contigs.fasta',
        qc_done = 'data/qc_done'
    group: 'assembly'
    output:
        final_link = 'assembly/final_contigs.fasta',
        sizes = "www/assembly_stats.txt"
    shell:
        'mkdir -p www/; '
        'stats.sh {input.fasta} > {output.sizes}; '
        'mkdir -p assembly; '
        'cd assembly; '
        'ln -s ../data/final_contigs.fasta ./; '
        'rm -rf data/racon_polishing; '

rule reset_to_spades_assembly:
    output:
         temp('data/reset_spades')
    shell:
         'rm -rf data/spades*; '
         'rm -rf assembly/; '
         'rm -rf data/final_contigs.fasta; '
         'rm -rf data/flye_high_cov.fasta; '
         'rm -rf data/list_of*; '
         'rm -rf data/binned_reads; '
         'rm -rf data/final_assemblies; '
         'rm -rf data/sr_vs*; '
         'rm -rf data/unicycler*; '
         'rm -rf data/metabat_bins; '
         'rm -rf data/cached_bams; '
         'touch data/reset_spades'

rule remove_final_contigs:
    output:
          temp('data/rewind_time')
    shell:
         'rm -rf assembly/; '
         'rm -rf data/final_contigs.fasta; '
         'touch data/rewind_time; '