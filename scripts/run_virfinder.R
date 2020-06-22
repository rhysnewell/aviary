#!/usr/bin/env Rscript

# Title     : Wrapper for virfinder
# Objective : Run virfinder
# Created by: uqrnewe1
# Created on: 18/06/2020

library(VirFinder)

fasta <- commandArgs(TRUE)[1]
out <- commandArgs(TRUE)[2]

predResult <- VF.pred(fasta)
# predResult$qvalue <- VF.qvalue(predResult$pvalue)
write.table(predResult, out, sep='\t', row.names=F, quote=F)