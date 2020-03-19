
import os
from os import path

configfile: "config.yml"
workdir: path.join(config["workdir_top"], config["pipeline"])

WORKDIR = path.join(config["workdir_top"], config["pipeline"])
SNAKEDIR = path.dirname(workflow.snakefile)

include: "snakelib/utils.snake"

in_fastq = config["reads_fastq"]
if not path.isabs(in_fastq):
    in_fastq = path.join(SNAKEDIR, in_fastq)

in_genome = config["genome_fasta"]
if not path.isabs(in_genome):
    in_genome = path.join(SNAKEDIR, in_genome)

in_annotation = config["existing_annotation"]
if not path.isabs(in_annotation) and in_annotation != "":
    in_annotation = path.join(SNAKEDIR, in_annotation)

rule dump_versions:
    output:
        ver = "versions.txt"
    conda: "env.yml"
    shell:"""
    command -v conda > /dev/null && conda list > {output.ver}
    """

rule preprocess_reads:
    input:
        fq = in_fastq
    output:
        pq = "processed_reads/input_reads.fq",
        flq = "processed_reads/full_length_reads.fq",
    params:
        pc = config["run_pychopper"],
        read_stats = config["read_qc"],
        pc_opts = config["pychopper_opts"],
        concat = config["concatenate"],
    threads: config["threads"]
    run:
        shell("mkdir -p processed_reads")
        if params.concat:
            shell("find %s  -regextype posix-extended -regex '.*\.(fastq|fq)$' -exec cat {{}} \\; > processed_reads/input_reads.fq" % str(input.fq))
        else:
            shell("ln -s `realpath %s` processed_reads/input_reads.fq" % str(input.fq))
        if params.pc:
           shell("cd processed_reads; cdna_classifier.py -t %d %s input_reads.fq full_length_reads.fq" % (threads, params.pc_opts))
        else:
            shell("ln -s `realpath processed_reads/input_reads.fq` processed_reads/full_length_reads.fq")

rule build_minimap_index:
    input:
        genome = in_genome
    output:
        index = "index/genome_index.mmi"
    params:
        opts = config["minimap_index_opts"]
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
        minimap2 -t {threads} {params.opts} -I 1000G -d {output.index} {input.genome}
    """

rule map_reads:
    input:
       index = rules.build_minimap_index.output.index,
       fastq = rules.preprocess_reads.output.flq,
    output:
       bam = "alignments/reads_aln_sorted.bam"
    params:
        opts = config["minimap2_opts"],
        min_mq = config["minimum_mapping_quality"],
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} -ax splice {params.opts} {input.index} {input.fastq}\
    | samtools view -q {params.min_mq} -F 2304 -Sb | samtools sort -@ {threads} - -o {output.bam};
    samtools index {output.bam}
    """

skip_bundles = "NO"
if config["bundle_min_reads"] is False:
    skip_bundles = "YES"

checkpoint split_bam:
    input:
        bam = rules.map_reads.output.bam,
    output:
        directory("bam_bundles")
    conda: "env.yml"
    params:
        min_reads = config["bundle_min_reads"]
    threads: config["threads"]
    shell:"""
    if [ {skip_bundles} != "YES" ]
    then
        seqkit bam -j {threads} -N {params.min_reads} {input.bam} -o {output}
    else
        mkdir -p {output}
        ln -s {input.bam} {output}/000000000_ALL:0:1_bundle.bam
    fi
    """

rule run_stringtie:
    input:
        bundle = "bam_bundles/{bundle}.bam"
    output:
        gff = "gff_bundles/{bundle}.gff",
    params:
        opts = config["stringtie_opts"],
        guide = config["use_guide_annotation"],
        ann = in_annotation
    conda: "env.yml"
    threads: config["threads"]
    run:
        if params.guide and params.ann != "":
            print("Using guide annotation from: ", in_annotation)
            params.opts += " -G {} ".format(in_annotation)
        shell("""
            stringtie -L -v -p {} {} -o {} {}
        """.format(threads, params.opts, output.gff, input.bundle))

def gff_bundle_list(wildcards):
    checkpoint_output = checkpoint_output = checkpoints.split_bam.get(**wildcards).output[0]
    bundle = glob_wildcards(os.path.join(checkpoint_output, '{bundle}.bam')).bundle
    bundle = sorted(bundle, key= lambda x: int(x.split("_")[0]))
    gffs = expand('gff_bundles/{bundle}.gff', bundle=bundle)
    return gffs

rule merge_gff_bundles:
    input:
        gffs = gff_bundle_list
    output:
        merged = "results/annotation/str_merged.gff"
    shell: """
    mkdir -p results/annotation
    echo '#gff-version 2' >> {output.merged}
    echo '#pipeline-nanopore-isoforms: stringtie' >> {output.merged}
    for i in {input.gffs}
    do
        grep -v '#' $i >> {output.merged}
    done
    """

rule run_gffcompare:
    input:
        gff = rules.merge_gff_bundles.output.merged,
    output:
        cmp_dir = directory("results/gffcompare")
    params:
        ann = in_annotation,
        opts = config["gffcompare_opts"],
    run:
        shell("mkdir -p {}".format(output.cmp_dir))
        gp = os.path.basename(input.gff)
        ing = os.path.abspath(input.gff)
        if in_annotation != "":
            shell("cd {}; ln -s {} {};gffcompare -o str_merged -r {} {} {}".format(output.cmp_dir, ing, gp, params.ann, params.opts, gp))
            shell("cd {}; {} -r {} str_merged.stats".format(output.cmp_dir, os.path.join(SNAKEDIR,"scripts","plot_gffcmp_stats.py"), "str_gffcmp_report.pdf"))

rule run_gffread:
    input:
        gff = rules.merge_gff_bundles.output.merged,
        genome = in_genome,
        gff_cmp = rules.run_gffcompare.output.cmp_dir,
    output:
        fas = "results/str_transcriptome.fas",
        merged_fas = "results/merged_transcriptome.fas",
    shell: """
        gffread -g {input.genome} -w {output.fas} {input.gff}
        if [ -f {input.gff_cmp}/str_merged.annotated.gtf ]
        then
            gffread -F -g {input.genome} -w {output.merged_fas} {input.gff_cmp}/str_merged.annotated.gtf
        else
            touch {output.merged_fas}
        fi
    """


rule all: ## run the whole pipeline
    input:
        ver = rules.dump_versions.output.ver,
        index = rules.build_minimap_index.output.index,
        aligned_reads = rules.map_reads.output.bam,
        bam_bundles = rules.split_bam.output,
        str_gff = rules.merge_gff_bundles.output,
        cmp_dir = rules.run_gffcompare.output.cmp_dir,
        trs = rules.run_gffread.output.fas,
