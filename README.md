![ONT_logo](/ONT_logo.png)
-----------------------------

Pipeline for annotating genomes using long read transcriptomics data with stringtie and other tools
===================================================================================================

This `snakemake` pipeline creates GFF annotation from Oxford Nanopore cDNA or direct RNA reads.
The cDNA reads are optionally processed by [pychopper](https://github.com/nanoporetech/pychopper) for trimming and orientation. The prcoessed reads are mapped to the reference genome using [minimap2](https://github.com/lh3/minimap2), and then processed by [stringtie](http://ccb.jhu.edu/software/stringtie) in long read mode (with or withouth using a guide annotation) to generate the GFF annotation.

If an existing annotation is available, the inferred annotation is compared to it using [gffcompare](http://ccb.jhu.edu/software/stringtie/gffcompare.shtml) and a report is genereated.

Getting Started
===============

## Input

- The input reads must be in fastq format. 
- The input genome must be in fasta format.
- The existing annotations must be in GFF3 or GFF2 format.

## Output
- The output is of a GFF ('query') file format - this is the annotated genome prodcuced by the pipeline.

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

## Performance on a D. melanogaster PCS109 data

### Performance with guide annotation

Run `./run_evaluation_dmel.sh` to reproduce these results.

![dmel_guide](https://raw.githubusercontent.com/nanoporetech/pipeline-nanopore-ref-isoforms/master/evaluation/results/dmel/guide/str_gffcmp_report_f.png)

See [this folder](https://github.com/nanoporetech/pipeline-nanopore-ref-isoforms/tree/master/evaluation/results/dmel/guide) for raw gffcompare output and the full PDF report.

### Performance without guide annotation

![dmel_no_guide](https://raw.githubusercontent.com/nanoporetech/pipeline-nanopore-ref-isoforms/master/evaluation/results/dmel/no_guide/str_gffcmp_report_f.png)

See [this folder](https://github.com/nanoporetech/pipeline-nanopore-ref-isoforms/tree/master/evaluation/results/dmel/no_guide) for raw gffcompare output and the full PDF report.
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

### Research Release

Research releases are provided as technology demonstrators to provide early access to features or stimulate Community development of tools. Support for this software will be minimal and is only provided directly by the developers. Feature requests, improvements, and discussions are welcome and can be implemented by forking and pull requests. However much as we would like to rectify every issue and piece of feedback users may have, the developers may have limited resource for support of this software. Research releases may be unstable and subject to rapid iteration by Oxford Nanopore Technologies.
