---
title: Performance and resources
---

# Performance and resources

Aviary workloads vary with read volume, assembly complexity, sample count,
enabled binners and reference databases. The repository does not provide a
portable benchmark dataset, so resource choices should be based on a dry run,
local benchmark records and representative pilot samples.

## CPU capacity

`--max-threads` limits the threads assigned to an individual process. The
default is 8. `--n-cores` controls the total CPU capacity Snakemake may schedule
and defaults to 16. Aviary raises `--n-cores` when it is lower than
`--max-threads`, because a scheduled process must fit within total capacity.

Not every upstream program scales to the requested thread count. Increasing
`--n-cores` can improve throughput when independent rules are ready to run, but
only when memory and storage bandwidth can support that concurrency.

## Memory

`--max-memory` is expressed in gigabytes and defaults to 250. It caps workflow
resource requests; it does not reserve that memory or predict consumption.
Reference-heavy stages, large assemblies and concurrent jobs can dominate RAM.
On a scheduler, ensure the coordinator and submitted jobs use limits consistent
with the Aviary arguments.

## Storage and temporary files

Assembly, mapping and binning can create substantial intermediate FASTQ, BAM
and index data. `--tmpdir` selects temporary storage; when omitted, Aviary uses
`TMPDIR`. Prefer node-local scratch when it is large enough and retained for the
duration of the rule. The default `--clean` behaviour removes declared
temporary files after successful consumption.

## Multiple samples

Coverage calculation can be split across jobs. With the default
`--coverage-job-strategy`, Aviary splits coverage work when more than 10 samples
are present; `--coverage-samples-per-job` defaults to 5. `always` and `never`
override that selection. More jobs can improve scheduler utilisation but add
scheduling and file-system overhead.

## GPU execution

GPU-enabled binners require compatible Aviary GPU environments and hardware.
`--request-gpu` requests a GPU only for cluster execution, and GPU-specific
binner options affect only their corresponding stages. Assembly, mapping and
many annotation steps remain CPU workloads.

## Measure the current system

Snakemake writes per-rule measurements beneath `benchmarks/`. Use several
representative runs to set scheduler resources, and retain input sizes,
software versions and database releases with any reported measurements.
