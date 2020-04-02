![ONT_logo](/ONT_logo.png)
-----------------------------

Pipeline for annotating genomes using long read transcriptomics data with stringtie and other tools
===================================================================================================

Getting Started
===============

## Input

- The input reads must be in fastq format. 
- The input genome must be in fasta format.

## Output

The pipeline produces the following output:

- TODO

## Depedencies

- [miniconda](https://conda.io/miniconda.html)
- [snakemake](http://snakemake.readthedocs.io/en/latest/) - easily installed via conda
- The rest of the dependencies are installed via conda.

## Installation

Clone the pipeline and the pinfish toolset by issuing:

```bash
git clone --recursive https://github.com/nanoporetech/pipeline-nanopore-ref-isoforms.git
```

## Usage

Edit `config.yml` to set the input genome, input fastq and parameters, then issue:

```bash
snakemake --use-conda -j <num_cores> all
```

Results
=======

## Performance on D. melanogaster PCS109 data

Help
=====

## Licence and Copyright

(c) 2020 Oxford Nanopore Technologies Ltd.

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this
file, You can obtain one at http://mozilla.org/MPL/2.0/.

## FAQs and tips

- The [GFF2](https://www.ensembl.org/info/website/upload/gff.html) files can be visualised using [IGV](http://software.broadinstitute.org/software/igv).

## References and Supporting Information

See the post announcing the tool at the Oxford Nanopore Technologies community [here](https://community.nanoporetech.com/posts/new-transcriptomics-analys).
