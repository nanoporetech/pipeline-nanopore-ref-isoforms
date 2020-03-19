#!/bin/bash

FASTQ_URL="http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/005/ERR3588905/ERR3588905_1.fastq.gz"
REF_URL="http://ftp.ensembl.org/pub/release-99/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz"
CONF="evaluation/config_dmel.yml"
CORES=30

OUT_DIR="evaluation/pipeline-evaluation-dmel/results"
RES_DIR="$OUT_DIR/results"
DATA_DIR="$OUT_DIR/data"
FASTQ="$DATA_DIR/ERR3588905_1.fastq"
REF="$DATA_DIR/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz"

rm -fr $OUT_DIR/results
mkdir -p $OUT_DIR/data

if [ ! -f $REF ];
then (cd $DATA_DIR; wget $REF_URL)
fi

if [ ! -f $FASTQ ];
then (cd $DATA_DIR; wget $FASTQ_URL; gzip -d ${FASTQ}.gz)
fi

if [ ! -f $FASTQ_FL ];
then cdna_classifier.py -t $CORES ${FASTQ} $FASTQ_FL
fi

snakemake -j $CORES all --configfile $CONF
