#!/usr/bin/env python3

import pandas as pd
import re

def create_conversion_dict(ids, labels):
    map_ids_metabuli = {}
    taxstring = []
    tax_regex = re.compile("^[pcofgs]__")
    for k, v in zip(ids, labels):
        tax = v.strip()
        if tax in ["unclassified", "root"]:
            continue
        if tax.startswith("d__"):
            taxstring = [tax]
        else:
            tax_tag = re.match(tax_regex, tax)
            if tax_tag:
                replace_pos = [t.startswith(tax_tag[0]) for t in taxstring]
                if sum(replace_pos) > 0:
                    taxstring = taxstring[:replace_pos.index(True)]
                taxstring.append(tax)
            else:
                # For assignments to specific genome ids
                pass

        map_ids_metabuli[k] = ",".join(taxstring)

    return map_ids_metabuli

def read_classifications(classifications):
    """Read Metabuli's tax_classifications.tsv into a DataFrame with integer
    column labels 0..7.

    usecols=range(8) pins the schema to Metabuli's 8 real columns so that both
    classified (is_classified=1) and unclassified (is_classified=0) rows parse
    identically. Unclassified rows carry a stray trailing tab -> a phantom
    empty 9th field; without usecols the C parser sees inconsistent field
    counts (8 vs 9) and raises "Expected 8 fields, saw 9". Selecting the 8 real
    columns absorbs the phantom field (it is always empty) without discarding
    any real data. The leading '#is_classified' header line is read as row 0
    (header=None); it carries a non-numeric taxID and is dropped downstream by
    the inner merge in process_classifications, same as before.
    """
    return pd.read_csv(classifications, delimiter='\t', header=None, usecols=range(8))


def process_classifications(df_clas, map_ids_metabuli, coverm_out):
    df_clas['contigs'] = df_clas[1]
    df_clas['predictions'] = df_clas[2].map(map_ids_metabuli)

    out = pd.merge(df_clas, coverm_out, on="contigs", how="inner")

    return out[['contigs', 'predictions']]


if __name__ == '__main__':
    report = snakemake.params.report
    classifications = snakemake.params.classifications
    filt_cov = snakemake.input.filt_cov
    output_file = snakemake.output[0]
    threads = snakemake.threads
    # Prefer resources.log_path (new pattern); fallback to log[0] if present
    log = getattr(getattr(snakemake, 'resources', object()), 'log_path', None)
    if not log:
        log = snakemake.log[0]

    with open(log, "w") as logf:
        logf.write("Starting conversion of metabuli predictions\n")

    df_report = pd.read_csv(report, delimiter='\t', header=None)

    map_ids_metabuli = create_conversion_dict(df_report[4], df_report[5])
    with open(log, 'a') as logf:
        logf.write(f"Conversion dictionary created\n")

    # df_clas: classified_flag, read_id, taxonomy_identifier, ...
    df_clas = read_classifications(classifications)
    with open(log, 'a') as logf:
        logf.write(f"{len(df_clas)} classifications read\n")

    coverm_out = pd.read_csv(filt_cov, sep='\t')[['contigName']]
    coverm_out.rename(columns={"contigName": "contigs"}, inplace=True)

    df_clas = process_classifications(df_clas, map_ids_metabuli, coverm_out)
    with open(log, 'a') as logf:
        logf.write(f"Conversion complete\n")
        logf.write(f"{len(df_clas)} contigs passed filtering\n")

    df_clas.to_csv(output_file, sep='\t', index=False)
