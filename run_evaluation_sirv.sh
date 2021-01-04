#!/bin/bash

CONF="config.yml"
#rm -fr evaluation/pipeline-nanopore-ref-isoforms
snakemake --use-conda -j 3 all --configfile  $CONF
