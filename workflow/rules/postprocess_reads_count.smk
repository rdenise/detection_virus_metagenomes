##########################################################################
##########################################################################

rule get_contig_coverage_info:
    input:
        pileup=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.pileup",
        ),
        idxstats=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "samtools",
            "{sample}",
            "{sample}_sorted.idxstats",
        ),
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "viral_contigs",
            "viral_contigs_over_3kb.fna",
        ),
    output:
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "mapping",
            "{sample}",
            "{sample}_mpileup_cov.tsv"
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "mapping_postprocess",
            "{sample}.get_contig_coverage_info.log",
        ),
    resources:
        cpus=1,
    threads: 1
    conda:
        "../envs/biopython.yaml"    
    script:
        "../scripts/process_mpileup.py"

##########################################################################
##########################################################################

rule get_total_reads_mapping_counts:
    input:
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
            "processed_files",
            "mapping",
            "total_reads_per_sample.tsv"
        )   
    resources:
        cpus=1,
    threads: 1
    shell:
        """
        echo "sample\ttotal_reads" > {output:q}
        for tmp_file in {input}
        do 
            echo -e "$(basename $tmp_file _sorted.flags)\t$(sed -n 1p $tmp_file | awk '{{print $1}}')" >> {output:q}
        done 
        """

##########################################################################
##########################################################################

rule make_single_count_table:
    input:
        total_reads=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "mapping",
            "total_reads_per_sample.tsv"
        ), 
        mpileup_cov=expand(
            os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "mapping",
                "{sample}",
                "{sample}_mpileup_cov.tsv"
            ), 
            sample=FASTQ_SAMPLE,
        )
    output:
        long_format=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "mapping",
            "total_sum_normalised_long.tsv"
        ), 
        otu=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "mapping",
            "otu_total_sum_normalised.tsv"
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "mapping_postprocess",
            "make_single_count_table.log",
        ),
    resources:
        cpus=1,
    threads: 1
    conda:
        "../envs/pandas_plots.yaml"
    script:
        "../scripts/get_count_table.py"


##########################################################################
##########################################################################
