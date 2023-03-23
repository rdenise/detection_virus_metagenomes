##########################################################################
##########################################################################

rule fastqc:
    input:
        os.path.join(
            FASTQ_FOLDER,
            f"{{sample}}{FASTQ_EXT}",
        ),
    output:
        html=os.path.join(
            OUTPUT_FOLDER,
            "qc",
            "fastqc",
            "untrimmed",
            "{sample}.html",
        ),
        zip=os.path.join(
            OUTPUT_FOLDER,
            "qc",
            "fastqc",
            "untrimmed",
            "{sample}_fastqc.zip",
        )
    params: 
        "--quiet"
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "fastqc",
            "untrimmed",
            "{sample}.log",
        ),
    threads: 1
    wrapper:
        "v1.24.0/bio/fastqc"


##########################################################################
##########################################################################

rule multiqc:
    input: 
        expand(
            os.path.join(
                OUTPUT_FOLDER,
                "qc",
                "fastqc",
                "untrimmed",
                f"{{sample}}{FASTQ_SEP}{{sens}}.html",
            ), 
            sample=FASTQ_SAMPLE,
            sens=["1","2"]
        ),
        expand(
            os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "reads_trimmed",
                "report",
                "json",
                "{sample}.json",
            ), 
            sample=FASTQ_SAMPLE,
        )
    output: 
        os.path.join(
            OUTPUT_FOLDER,
            "qc",
            "multiqc_report.untrimmed.html",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "multiqc",
            "multiqc_report.untrimmed.log",
        ),
    threads: 1
    wrapper: "v1.24.0/bio/multiqc"

##########################################################################
##########################################################################

rule fastp_pe:
    input:
        sample=[
            os.path.join(
                FASTQ_FOLDER,
                f"{{sample}}{FASTQ_SEP}1{FASTQ_EXT}"
            ),
            os.path.join(
                FASTQ_FOLDER,
                f"{{sample}}{FASTQ_SEP}2{FASTQ_EXT}"
            ),
        ]
    output:
        trimmed=[
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
        ],
        # Unpaired reads separately
        unpaired1=os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "reads_trimmed",
                "unpaired",
                "{sample}.u1.fastq.gz",
            ),
        unpaired2=os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "reads_trimmed",
                "unpaired",
                "{sample}.u2.fastq.gz",
            ),
        # or in a single file
        # unpaired="reads_trimmed/pe/{sample}.singletons.fastq",
        merged=os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "reads_trimmed",
                "merged",
                "{sample}.merged.fastq.gz",
            ),
        failed=os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "reads_trimmed",
                "failed",
                "{sample}.failed.fastq.gz",
            ),
        html=os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "reads_trimmed",
                "report",
                "html",
                "{sample}.html",
            ),
        json=os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "reads_trimmed",
                "report",
                "json",
                "{sample}.json",
            ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "fastp",
            "{sample}.log",
        ),
    params:
        adapters="--detect_adapter_for_pe",
        extra="--merge --length_required 60 --cut_front --cut_tail --cut_window_size 4 --cut_mean_quality 20"
    threads: 5
    wrapper:
        "v1.24.0/bio/fastp"

##########################################################################
##########################################################################

rule fastqc_trimmed:
    input:
        os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "reads_trimmed",
            "paired",
            "{sample}.fastq.gz",
        ),
    output:
        html=os.path.join(
            OUTPUT_FOLDER,
            "qc",
            "fastqc",
            "trimmed",
            "{sample}.html",
        ),
        zip=os.path.join(
            OUTPUT_FOLDER,
            "qc",
            "fastqc",
            "trimmed",
            "{sample}_fastqc.zip",
        ),
    params: 
        "--quiet"
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "fastqc",
            "trimmed",
            "{sample}.log",
        ),
    threads: 1
    wrapper:
        "v1.24.0/bio/fastqc"


##########################################################################
##########################################################################

rule multiqc_trimmed:
    input: 
        expand(
            os.path.join(
                OUTPUT_FOLDER,
                "qc",
                "fastqc",
                "trimmed",
                "{sample}_R{sens}.html",
            ), 
            sample=FASTQ_SAMPLE,
            sens=["1","2"]
        ),
        expand(
            os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "reads_trimmed",
                "report",
                "json",
                "{sample}.json",
            ), 
            sample=FASTQ_SAMPLE,
        ),
        expand(
            os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "samtools",
                "{sample}",
                "{sample}_sorted.idxstats",
            ),
            sample=FASTQ_SAMPLE,
        ),
        expand(
            os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "samtools",
                "{sample}",
                "{sample}_sorted.flags",
            ),
            sample=FASTQ_SAMPLE,
        )
    output: 
        os.path.join(
            OUTPUT_FOLDER,
            "qc",
            "multiqc_report.trimmed.html",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "multiqc",
            "multiqc_report.trimmed.log",
        ),
    threads: 1
    wrapper: "v1.24.0/bio/multiqc"

##########################################################################
##########################################################################
