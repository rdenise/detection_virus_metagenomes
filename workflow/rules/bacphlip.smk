##########################################################################
##########################################################################


rule predict_viral_lifestyles:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "viral_contigs",
            "viral_contigs_over_3kb.fna",
        ),
    output:
        bacphlip=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "bacphlip_out",
            "viral_contigs_over_3kb.fna.bacphlip",
        ),
        hmmsearch=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "bacphlip_out",
            "viral_contigs_over_3kb.fna.hmmsearch.tsv",
        ),
    params:
        bacphlip=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "viral_contigs",
            "viral_contigs_over_3kb.fna.bacphlip",
        ),
        hmmsearch=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "viral_contigs",
            "viral_contigs_over_3kb.fna.hmmsearch.tsv",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bacphlip",
            "bacphlip.log",
        ),
    conda:
        "../envs/bacphlip.yaml"
    shell:
        """
        bacphlip -i {input.fasta:q} --multi_fasta -f &> {log:q}
        mv {params.bacphlip:q} {output.bacphlip:q}
        mv {params.hmmsearch:q} {output.hmmsearch:q}
        """


##########################################################################
##########################################################################


rule parse_bacphlip:
    input:
        bacphlip=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "bacphlip_out",
            "viral_contigs_over_3kb.fna.bacphlip",
        ),
    output:
        modify=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "bacphlip_out",
            "viral_contigs_over_3kb.fna.bacphlip.tsv",
        ),
    params:
        cutoff=0.95,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bacphlip",
            "parse_bacphlip.log",
        ),
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/parse_bacphlip_results.py"


##########################################################################
##########################################################################
