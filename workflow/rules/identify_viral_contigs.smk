##########################################################################
##########################################################################


rule run_genomad:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "contigs",
            "human_filtered",
            "all_contigs.over3kb.nr.filtered.fasta",
        ),
    output:
        virus_summary=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "genomad",
            "human_filtered",
            "all_contigs.over3kb.nr.filtered_summary",
            "all_contigs.over3kb.nr.filtered_virus_summary.tsv",
        ),
        provirus_classification=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "genomad",
            "human_filtered",
            "all_contigs.over3kb.nr.filtered_aggregated_classification",
            "all_contigs.over3kb.nr.filtered_provirus_aggregated_classification.tsv",
        ),
        aggregate_classification=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "genomad",
            "human_filtered",
            "all_contigs.over3kb.nr.filtered_aggregated_classification",
            "all_contigs.over3kb.nr.filtered_aggregated_classification.tsv",
        ),
    params:
        db=config["genomad"]["path"],
        composition=config["genomad"]["composition"],
        virus_summary=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "genomad",
            "human_filtered",
        ),
    conda:
        "../envs/genomad.yaml"
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "genomad",
            "genomad.human_filtered.log",
        ),
    threads: 20
    shell:
        """
        genomad end-to-end --splits 12 --cleanup --threads {threads} --composition {params.composition} \
        --enable-score-calibration {input.contigs:q} {params.genomad_dir:q} {params.db:q} &> {log:q}
        """


##########################################################################
##########################################################################


rule reduce_dataset:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "contigs",
            "human_filtered",
            "all_contigs.over3kb.nr.filtered.fasta",
        ),
        virus_summary=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "genomad",
            "human_filtered",
            "all_contigs.over3kb.nr.filtered_summary",
            "all_contigs.over3kb.nr.filtered_virus_summary.tsv",
        ),
    output:
        fasta=temp(
            os.path.join(
                OUTPUT_FOLDER,
                "databases",
                "contigs",
                "human_filtered",
                "all_contigs.over3kb.nr.filtered.blast_ready.fasta",
            )
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "postprocess_detection",
            "all_contigs.combine_all.log",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script:
        "../scripts/remove_genomad_hits.py"


##########################################################################
##########################################################################


rule combine_all:
    input:
        tsv_blast=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "virus",
            f"merge.eval_{blast_evalue:.0e}.cov_{blast_coverage}.pident_{blast_pident}.annotation.blasn.tsv",
        ),
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "contigs",
            "human_filtered",
            "all_contigs.over3kb.nr.filtered.fasta",
        ),
        virus_summary=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "genomad",
            "human_filtered",
            "all_contigs.over3kb.nr.filtered_summary",
            "all_contigs.over3kb.nr.filtered_virus_summary.tsv",
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "viral_contigs",
            "all_contigs.nr.selected.fasta",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "postprocess_detection",
            "all_contigs.combine_all.log",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script:
        "../scripts/combine_all.py"


##########################################################################
##########################################################################


rule removing_provirus_contigs_below_3kb:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "checkv",
            "combined.fna",
        ),
    output:
        reduced=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "viral_contigs",
            "viral_contigs_over_3kb.fna",
        ),
    conda:
        "../envs/bbmap.yaml"
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "bbmap",
            "reformat.virus.over3kb.log",
        ),
    shell:
        """
        reformat.sh in={input.contigs:q} out={output.reduced:q} minlength=3000 &> {log:q} 
        """


##########################################################################
##########################################################################


rule run_genomad_annotate:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "viral_contigs",
            "viral_contigs_over_3kb.fna",
        ),
    output:
        annotate=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "genomad",
            "viral_contigs",
            "all_contigs.over3kb.nr_annotate",
            "all_contigs.over3kb.nr_taxonomy.tsv",
        ),
    params:
        db=config["genomad"]["path"],
        virus_summary=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "genomad",
            "viral_contigs_annotation",
        ),
    conda:
        "../envs/genomad.yaml"
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "genomad",
            "genomad.log",
        ),
    threads: 20
    shell:
        """
        genomad annotate --splits 12 --cleanup --threads {threads} \
        {input.contigs:q} {params.genomad_dir:q} {params.db:q} &> {log:q}
        """


##########################################################################
##########################################################################


rule convert_blast_to_ani:
    input:
        blastout=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "otu",
            "viral_contigs_over_3kb_all_VS_all.blastout",
        ),
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "otu",
            "viral_contigs_over_3kb_all_VS_all.ani",
        ),
    params:
        anicalc="../scripts/anicalc.py",
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "otu",
            "anicalc.log",
        ),
    conda:
        "../envs/biopython.yaml"
    shell:
        """
        python {params.anicalc:q} -i {input.blastout:q} -o {output:q} &> {log:q}
        """


##########################################################################
##########################################################################


rule cluster_viruses_into_vOTUs:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "viral_contigs",
            "viral_contigs_over_3kb.fna",
        ),
        ani=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "otu",
            "viral_contigs_over_3kb_all_VS_all.ani",
        ),
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "otu",
            "viral_contigs_over_3kb_otu.tsv",
        ),
    params:
        aniclust="../scripts/aniclust.py",
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "otu",
            "aniclust.log",
        ),
    conda:
        "../envs/biopython.yaml"
    shell:
        """
        python {params.aniclust:q} --fna {input.fasta:q} --ani {input.ani:q} --out {output:q} &> {log:q}
        """


##########################################################################
##########################################################################


rule get_list_of_representative_vOTUs:
    input:
        os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "otu",
            "viral_contigs_over_3kb_otu.tsv",
        ),
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "otu",
            "final_vOTUs_representative_contigs.txt",
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "otu",
            "get_list_of_representative_vOTUs.log",
        ),
    shell:
        """
        cut -f 1 {input:q} > {output:q} 2> {log:q}
        """


##########################################################################
##########################################################################


rule get_representative_vOTU_contigs:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "viral_contigs",
            "viral_contigs_over_3kb.fna",
        ),
        vOTU_list="data/processed/virome/viral_contig_identification/final_vOTUs_representative_contigs.txt",
    output:
        os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "otu",
            "final_vOTUs_representative_contigs.txt",
        ),
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        filterbyname.sh in={input.fasta:q} out={output:q} names={input.vOTU_list:q} include=t &> {log:q}
        """


##########################################################################
##########################################################################
