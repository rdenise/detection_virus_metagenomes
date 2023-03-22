##########################################################################
##########################################################################


rule checkv_setup:
    output:
        directory(
            os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "checkv_db",
                "checkv-db-v1.5",
            )
        ),
    log:
        os.path.join(OUTPUT_FOLDER, "logs", "checkv", "checkv_setup.log"),
    params:
        checkv_db=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "checkv_db",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/checkv.yaml"
    threads: 1
    shell:
        """
        checkv download_database {params.checkv_db:q} &> {log:q}
        """


##########################################################################
##########################################################################


rule checkv_run:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "viral_contigs",
            "all_contigs.nr.selected.fasta",
        ),
        database=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "checkv_db",
            "checkv-db-v1.5",
        ),
    output:
        viruses=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "checkv",
            "viruses.fna",
        ),
        proviruses=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "checkv",
            "proviruses.fna",
        ),
        contamination=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "checkv",
            "contamination.tsv",
        ),
        quality_summary=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "checkv",
            "quality_summary.tsv",
        ),
    params:
        output_dir=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "checkv",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "checkv",
            "checkv_run.log"
        ),
    resources:
        cpus=5,
    conda:
        "../envs/checkv.yaml"
    threads: 10
    shell:
        """
        checkv end_to_end {input.fasta:q} {params.output_dir:q} \
        -t {threads} -d {input.database:q} &> {log:q}

        for myfile in {output}
        do 
            touch $myfile
        done
        """


##########################################################################
##########################################################################


rule combine_checkv_contigs:
    input:
        viruses=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "checkv",
            "viruses.fna",
        ),
        proviruses=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "checkv",
            "proviruses.fna",
        ),
    output:
        combine=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "checkv",
            "combined.fna",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "checkv",
            "checkv_cat.log"
        ),
    threads: 1
    shell:
        "cat {input.proviruses:q} {input.viruses:q} > {output.combine:q} 2> {log:q}"


##########################################################################
##########################################################################