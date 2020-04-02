#!/bin/bash

CONF="config.yml"
rm -fr evaluation/pipeline-nanopore-ref-isoforms
snakemake --use-conda -j 20 all --configfile  $CONF
