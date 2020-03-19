#!/bin/bash

CONF="config.yml"
rm -fr evaluation/pipeline-evaluation-sirv
snakemake -j 20 all --configfile  $CONF
