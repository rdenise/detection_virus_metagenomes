##########################################################################
##########################################################################


rule build_index:
    input:
        ref=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "viral_contigs",
            "viral_contigs_over_3kb.fna",
        ),
    output:
        multiext(
            os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "bowtie2",
                "index",
                "viral_contigs_over_3kb",
            ),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    params:
        extra="",  # optional parameters    
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bowtie2_build",
            "viral_contig.bowtie2_build.log",
        ),
    resources:
        cpus=1,
    threads: 8
    wrapper:
        "v1.24.0/bio/bowtie2/build"

##########################################################################
##########################################################################


rule map_reads:
    input:
        idx=multiext(
            os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "bowtie2",
                "index",
                "viral_contigs_over_3kb",
            ),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
        sample=[
            os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "reads_trimmed",
                "paired",
                "{sample}_R1.fastq.gz",
            ),
            os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "reads_trimmed",
                "paired",
                "{sample}_R2.fastq.gz",
            ),
        ]
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "samtools",
            "{sample}",                
            "{sample}.bam",
        )
    params:
        extra="",  # optional parameters      
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bowtie2",
            "{sample}.bowtie2.log",
        ),
    threads: 10  # Use at least two threads
    wrapper:
        "v1.24.0/bio/bowtie2/align" # It does bowtie and samtool view


##########################################################################
##########################################################################

rule samtools_sort:
    input:
        os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "samtools",
            "{sample}",
            "{sample}.bam",
        ),
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.bam",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "{sample}.samtools.sort.log",
        ),
    params:
        extra="-m 4G",
    threads: 10
    wrapper:
        "v1.24.0/bio/samtools/sort"


##########################################################################
##########################################################################

rule samtools_index:
    input:
        os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.bam",
        ),
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.bam.bai",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "{sample}.samtools.index.log",
        ),
    params:
        extra="-b",  # optional params string
    threads: 10  # This value - 1 will be sent to -@
    wrapper:
        "v1.24.0/bio/samtools/index"

##########################################################################
##########################################################################

rule samtools_idxstats:
    input:
        bam=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.bam",
        ),
        idx=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.bam.bai",
        ),
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.idxstats",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "{sample}.samtools.idxstats.log",
        ),
    params:
        extra="",  # optional params string
    wrapper:
        "v1.24.0/bio/samtools/idxstats"


##########################################################################
##########################################################################

rule samtools_flagstat:
    input:
        os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.bam",
        ),
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.flags",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "{sample}.samtools.flagstat.log",
        ),
    params:
        extra="",  # optional params string
    wrapper:
        "v1.24.0/bio/samtools/flagstat"


##########################################################################
##########################################################################

rule samtools_pileup:
    input:
        sortedbam=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.bam",
        ),
    output:
        pileup=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.pileup",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "samtools",
            "{sample}.samtools.mpileup.log",
        ),
    threads:
        10
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools mpileup {input.sortedbam:q} > {output.pileup:q} 2> {log:q}
        """

##########################################################################
##########################################################################