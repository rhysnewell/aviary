onsuccess:
    print("Strain analysis finished, no error")

onerror:
    print("An error occurred")


rule generate_bams:
    input:
        bins_directory = 'data/galah_bins/',
    output:
        finished_mapping = 'data/binned_bams/{sample}.bam'
    threads:
        config['max_threads']
    conda:
        '../../envs/coverm.yaml'
    shell:
        'coverm genome '


rule lorikeet:
    input:
         finished_binning = 'data/galah_bins/done',
         bins_directory = 'data/galah_bins/',
         short_reads = config[]