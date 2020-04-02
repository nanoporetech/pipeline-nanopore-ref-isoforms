#!/bin/bash

FASTQ_URL="http://ftp.sra.ebi.ac.uk/vol1/fastq/ERR358/005/ERR3588905/ERR3588905_1.fastq.gz"
REF_URL="http://ftp.ensembl.org/pub/release-99/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa.gz"
GFF_URL="http://ftp.ensembl.org/pub/release-99/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.28.99.gff3.gz"
CONF="evaluation/config_dmel.yml"
CORES=30

OUT_DIR="evaluation/pipeline-evaluation-dmel"
RES_DIR="$OUT_DIR/results"
DATA_DIR="$OUT_DIR/data"
FASTQ="$DATA_DIR/ERR3588905_1.fastq"
REF="$DATA_DIR/Drosophila_melanogaster.BDGP6.28.dna.toplevel.fa"
GFF="$DATA_DIR/Drosophila_melanogaster.BDGP6.28.99.gff3"

rm -fr $OUT_DIR/results
mkdir -p $OUT_DIR/data

if [ ! -f $REF ];
then (cd $DATA_DIR; curl -L -C - -O $REF_URL); gzip -d ${REF}.gz
fi

if [ ! -f $GFF ]
then
	(cd $DATA_DIR; curl -L -C - -O $GFF_URL); gzip -d ${GFF}.gz
fi

if [ ! -f $FASTQ ];
then (cd $DATA_DIR; curl -L -C - -O $FASTQ_URL); gzip -d ${FASTQ}.gz
fi

snakemake --use-conda -j $CORES all --configfile $CONF
