##########################################################################
##########################################################################

rule run_metaspades:
    input:
        reads=[
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
    output:
        contigs=os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "assemblies",
                "{sample}",
                "metaspades",
                "contigs.fasta",
            ),
        scaffolds=os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "assemblies",
                "{sample}",
                "metaspades",
                "scaffolds.fasta",
            ),
        dir=directory(
            os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "assemblies",
                "{sample}",
                "metaspades",
            )
        ),
    params:
        # all parameters are optional
        k="auto",
        extra="",
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "metaspades",
            "{sample}.log",
        ),
    threads: 8
    resources:
        mem_mem=250000,
        time=60 * 24,
    wrapper:
        "v1.24.0/bio/spades/metaspades"

##########################################################################
##########################################################################


rule assemble_combined_sample:
    input:
        r1=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "reads_trimmed",
            "paired",
            "{sample}_R1.fastq.gz",
        ),
        r2=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "reads_trimmed",
            "paired",
            "{sample}_R2.fastq.gz",
        ),
    output:
        contigs=os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "assemblies",
                "{sample}",
                "megahit",
                "final.contigs.fa",
            ),
        dir=directory(
            os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "assemblies",
                "{sample}",
                "megahit",
            )
        ),
    conda:
        "../envs/megahit.yaml"
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "metaspades",
            "{sample}.log",
        ),
    threads:
        10
    shell:
        """
        megahit -1 {input.r1:q} -2 {input.r2:q} -t {threads} -o {output.dir:q} -f &> {log:q}
        """

##########################################################################
##########################################################################

rule append_sample_names_to_subject_individual_assemblies:
    input:
        contigs=os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "assemblies",
                "{sample}",
                "{software}",
                "{file}",
        ),
    output:
        contigs=temp(
            os.path.join(
                    OUTPUT_FOLDER,
                    "processed_files",
                    "assemblies",
                    "{sample}",
                    "{software}",
                    "renamed.{file}",
            )
        ),
    params: 
        sample="{sample}"
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "sed",
            "{sample}.{software}.{file}.log",
        ),
    threads:
        1
    shell:
        """
        sed 's/^>/>{params}~/g' {input.contigs:q} > {output.contigs:q} 2> {log:q}
        """

##########################################################################
##########################################################################


rule combine_all_named_files:
    input:
        metaspades=expand(
            os.path.join(
                    OUTPUT_FOLDER,
                    "processed_files",
                    "assemblies",
                    "{sample}",
                    "metaspades",
                    "renamed.scaffolds.fasta",
            ), 
            sample=FASTQ_SAMPLE,
        ), 
        megahit=expand(
            os.path.join(
                    OUTPUT_FOLDER,
                    "processed_files",
                    "assemblies",
                    "{sample}",
                    "megahit",
                    "renamed.final.contigs.fa",
            ), 
            sample=FASTQ_SAMPLE,
        ), 
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "assemblies",
            "all_contigs.fasta",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "cat",
            "all_contigs.log",
        )
    shell:
        """
        cat {input.metaspades:q} {input.megahit:q} > {output:q} 2> {log:q}
        """

##########################################################################
##########################################################################

rule merge_existing_contigs:
    input:
        contigs=expand(
            os.path.join(
                CONTIGS_FOLDER,
                "{contigs_files}" + CONTIGS_EXT
            ),
            contigs_files=CONTIGS_NAMES,
        )
    output:
        merged=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "assemblies",
            "merge_contigs.fasta",
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "cat",
            "merge_existing_contigs.log",
        )
    shell:
        """
        cat {input.contigs:q} > {output.merge:q} 2> {log:q}
        """

##########################################################################
##########################################################################

rule remove_contigs_less_than_3kb:
    input:
        contigs=CONTIGS_FILES
    output:
        reduced=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "assemblies",
            "all_contigs.over3kb.fasta",
        )
    conda:
        "../envs/bbmap.yaml"
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bbmap",
            "reformat.over3kb.log",
        )
    shell:
        """
        reformat.sh in={input.contigs:q} out={output.reduced:q} minlength=3000 &> {log:q} 
        """

##########################################################################
##########################################################################
