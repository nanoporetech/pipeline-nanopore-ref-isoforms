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
        pc = "True" if config["run_pychopper"] else "False",
        pc_opts = config["pychopper_opts"],
        concat = "True" if config["concatenate"] else "False",
    threads: config["threads"]
    conda: "env.yml"
    shell:
        """
        mkdir -p processed_reads;

        if [[ {params.concat} == "True" ]];
        then
            find {input.fq}  -regextype posix-extended -regex '.*\.(fastq|fq)$' -exec cat {{}} \\; > processed_reads/input_reads.fq
        else
            ln -s `realpath {input.fq}` processed_reads/input_reads.fq
        fi

        if [[ {params.pc} == "True" ]];
        then
            cd processed_reads; cdna_classifier.py -t {threads} {params.pc_opts} input_reads.fq full_length_reads.fq
        else
            ln -s `realpath processed_reads/input_reads.fq` processed_reads/full_length_reads.fq
        fi
        """

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

ContextFilter = """AlnContext: { Ref: "%s", LeftShift: -%d, RightShift: %d, RegexEnd: "[Aa]{%d,}", Stranded: True, Invert: True, Tsv: "alignments/internal_priming_fail.tsv"} """ % (in_genome,
config["poly_context"], config["poly_context"], config["max_poly_run"])

ContextDumper = """AlnContext: { Ref: "%s", LeftShift: -%d, RightShift: %d, RegexEnd: ".+", Stranded: False, Invert: True, Tsv: "alignments/context.tsv"} """ % (in_genome, config["poly_context"], config["poly_context"])

rule map_reads:
    input:
       index = rules.build_minimap_index.output.index,
       fastq = rules.preprocess_reads.output.flq,
    output:
       bam = "alignments/reads_aln_sorted.bam",
       stats = "alignments/read_aln_stats.tsv",
    params:
        opts = config["minimap2_opts"],
        min_mq = config["minimum_mapping_quality"],
        flt = lambda x : ContextFilter,
        dump = lambda x : ContextDumper,
        context_plt = config["context_plt_alns"]
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    minimap2 -t {threads} -ax splice {params.opts} {input.index} {input.fastq}\
    | samtools view -q {params.min_mq} -F 2304 -Sb -\
    | seqkit bam -j {threads} -x -T '{params.flt}' -\
    | samtools sort -@ {threads} -o {output.bam} -;
    samtools index {output.bam}

    ((seqkit bam -s -j {threads} {output.bam} 2>&1)  | tee {output.stats} ) || true

    if [[ -s alignments/internal_priming_fail.tsv ]];
    then
        tail -n +2 alignments/internal_priming_fail.tsv | awk '{{print ">" $1 "\\n" $4 }}' - > alignments/context_internal_priming_fail_start.fasta
        tail -n +2 alignments/internal_priming_fail.tsv | awk '{{print ">" $1 "\\n" $6 }}' - > alignments/context_internal_priming_fail_end.fasta
    fi

    if [[ {params.context_plt} -gt 0 ]];
    then
        (seqkit bam -j {threads} -T '{params.dump}' {output.bam} >/dev/null) || true
        tail -n +2 alignments/context.tsv | shuf - > alignments/context_shuff.tsv
        csvtk -t -H filter2 -f '$3 == 1' alignments/context_shuff.tsv > alignments/context_shuff_plus.tsv
        LINES_PLUS={params.context_plt}
        TOTAL_PLUS=`wc -l alignments/context_shuff_plus.tsv | cut -d$' ' -f 1`
        if [[ $LINES_PLUS -gt $TOTAL_PLUS ]];
        then
            LINES_PLUS=$TOTAL_PLUS
        fi
        head -n $LINES_PLUS alignments/context_shuff_plus.tsv | awk '{{print ">" $1 "\\n" $4 }}' - > alignments/context_shuff_plus_start.fasta
        (seqkit sana -i fasta alignments/context_shuff_plus_start.fasta.tmp 2>/dev/null) > alignments/context_shuff_plus_start.fasta || true
        head -n $LINES_PLUS alignments/context_shuff_plus.tsv | awk '{{print ">" $1 "\\n" $6 }}' - > alignments/context_shuff_plus_end.fasta.tmp
        (seqkit sana -i fasta alignments/context_shuff_plus_end.fasta.tmp 2>/dev/null) > alignments/context_shuff_plus_end.fasta || true
        csvtk -t -H filter2 -f '$3 == "-1"' alignments/context_shuff.tsv > alignments/context_shuff_minus.tsv
        LINES_MINUS={params.context_plt}
        TOTAL_MINUS=`wc -l alignments/context_shuff_minus.tsv| cut -d$' ' -f 1`
        if [[ $LINES_MINUS -gt $TOTAL_MINUS ]];
        then
            LINES_MINUS=$TOTAL_MINUS
        fi
        head -n $LINES_MINUS alignments/context_shuff_minus.tsv | awk '{{print ">" $1 "\\n" $4 }}' - > alignments/context_shuff_minus_start.fasta.tmp
        (seqkit sana -i fasta alignments/context_shuff_minus_start.fasta.tmp 2>/dev/null) > alignments/context_shuff_minus_start.fasta || true
        head -n $LINES_MINUS alignments/context_shuff_minus.tsv | awk '{{print ">" $1 "\\n" $6 }}' - > alignments/context_shuff_minus_end.fasta.tmp
        (seqkit sana -i fasta alignments/context_shuff_minus_end.fasta.tmp 2>/dev/null) > alignments/context_shuff_minus_end.fasta || true
        rm alignments/context*.tsv
        rm alignments/*.tmp
    fi

    for fas in alignments/context*.fasta;
    do
        if [[ -s $fas ]];
        then
            NAME="${{fas%.*}}"
            spoa -r 1 -l 1 $fas | sed '1d' | awk '{{print ">" NR "\\n" $1}}' > ${{fas}}_aln
            hmmbuild -n $NAME --dna ${{NAME}}.hmm ${{fas}}_aln > /dev/null
            hmmlogo ${{NAME}}.hmm > ${{NAME}}.logo
            {SNAKEDIR}/scripts/skylign.py -r ${{NAME}}.png `realpath ${{NAME}}.hmm`
        fi
    done
    rm alignments/context*fasta*
    """


rule plot_aln_stats:
    input:
        stats = rules.map_reads.output.stats
    output:
        pdf = "alignments/read_aln_stats.pdf"
    conda: "env.yml"
    threads: config["threads"]
    shell:"""
    {SNAKEDIR}/scripts/plot_aln_stats.py {input.stats} -r {output.pdf}
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

use_guide = "NO"
if config["use_guide_annotation"] is True:
    use_guide = "YES"

str_threads = 10
if config["threads"] < str_threads:
    str_threads = config["threads"]

rule run_stringtie:
    input:
        bundle = "bam_bundles/{bundle}.bam"
    output:
        gff = "gff_bundles/{bundle}.gff",
    params:
        opts = config["stringtie_opts"],
        guide = use_guide,
        ann = in_annotation,
        trp = lambda x: "STR.{}.".format(int(str(x.bundle).split("_")[0]))
    conda: "env.yml"
    threads: str_threads
    shell:
        """
        G_FLAG=""
        if [[ {params.guide} == "YES" ]];
        then
            G_FLAG="-G {params.ann}"
        fi
        stringtie --rf $G_FLAG -l {params.trp} -L -v -p {threads} {params.opts} -o {output.gff} {input.bundle} 2>/dev/null
        """

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
    echo '#pipeline-nanopore-ref-isoforms: stringtie' >> {output.merged}
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
        gp = lambda x: os.path.basename(rules.merge_gff_bundles.output.merged),
        ing = lambda x: os.path.abspath(rules.merge_gff_bundles.output.merged),
    conda: "env.yml"
    shell:
        """
        mkdir -p {output.cmp_dir};
        cd {output.cmp_dir};
        ln -s {params.ing} {params.gp};
        if [[ -f {params.ann} ]];
        then
            gffcompare -o str_merged -r {params.ann} {params.opts} {params.gp}
            # {SNAKEDIR}/scripts/plot_gffcmp_stats.py -r str_gffcmp_report.pdf -t str_merged.tracking str_merged.stats;
        fi
        """

rule run_gffread:
    input:
        gff = rules.merge_gff_bundles.output.merged,
        genome = in_genome,
        gff_cmp = rules.run_gffcompare.output.cmp_dir,
    output:
        fas = "results/str_transcriptome.fas",
        merged_fas = "results/merged_transcriptome.fas",
    conda: "env.yml"
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
        stats = rules.plot_aln_stats.output.pdf,
        bam_bundles = rules.split_bam.output,
        str_gff = rules.merge_gff_bundles.output,
        cmp_dir = rules.run_gffcompare.output.cmp_dir,
        trs = rules.run_gffread.output.fas,
