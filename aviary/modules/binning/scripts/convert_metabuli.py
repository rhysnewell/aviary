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
    log = snakemake.log[0]

    with open(log, "w") as logf:
        logf.write("Starting conversion of metabuli predictions\n")

    df_report = pd.read_csv(report, delimiter='\t', header=None)

    map_ids_metabuli = create_conversion_dict(df_report[4], df_report[5])
    with open(log, 'a') as logf:
        logf.write(f"Conversion dictionary created\n")

    # df_clas: classified_flag, read_id, taxonomy_identifier, ...
    df_clas = pd.read_csv(classifications, delimiter='\t', header=None)
    with open(log, 'a') as logf:
        logf.write(f"{len(df_clas)} classifications read\n")

    coverm_out = pd.read_csv(filt_cov, sep='\t')[['contigName']]
    coverm_out.rename(columns={"contigName": "contigs"}, inplace=True)

    df_clas = process_classifications(df_clas, map_ids_metabuli, coverm_out)
    with open(log, 'a') as logf:
        logf.write(f"Conversion complete\n")
        logf.write(f"{len(df_clas)} contigs passed filtering\n")

    df_clas.to_csv(output_file, sep='\t', index=False)
