# Module containing all the ncbi-blast related rules

##########################################################################
##########################################################################


rule blastn:
    input:
        contig=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "contigs",
            "human_filtered",
            "all_contigs.over3kb.nr.filtered.blast_ready.fasta",
        ),
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "virus",
            "all_contigs.nr.evalue_{evalue}.{database}.blastn.outfmt6.txt",
        ),
    params:
        database=lambda wildcards: DB_DICT[wildcards.database]["path"],
        outfmt="6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle",
        evalue="{evalue}",
        options_blast="-num_alignments 25000",
        max_len=2000000,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "virus",
            "nr.evalue_{evalue}.{database}.blastn.outfmt6.log",
        ),
    conda:
        "../envs/blast.yaml"
    threads: 20
    shell:
        """
        if [[ $(stat -c "%s" {input.contig:q}) -gt $(stat -c "%s" {params.database:q}) ]]
        then
            blastn -task megablast -query {input.contig:q} -out {output.blast_out:q} \
                -db {params.database:q} -evalue {params.evalue} -outfmt {params.outfmt:q} \
                -num_threads {threads} -mt_mode 1 {params.options_blast} &> {log:q}
        else
            blastn -task megablast -query {input.contig:q} -out {output.blast_out:q} \
                -db {params.database:q} -evalue {params.evalue} -outfmt {params.outfmt:q} \
                -num_threads {threads} -mt_mode 0 {params.options_blast} &> {log:q}
        fi            
        """


##########################################################################
##########################################################################


rule merge_blastn:
    input:
        all_out=expand(
            os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "blast",
                "virus",
                "all_contigs.nr.evalue_{evalue:.0e}.{database}.blastn.outfmt6.txt",
            ),
            database=DB_DICT.keys(),
            evalue=[blast_evalue],
        ),
    output:
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "virus",
            "merge.eval_{evalue}.cov_{coverage}.pident_{pident}.annotation.blasn.tsv",
        ),
    params:
        minimum_length=config["default_blast_option"]["length_min"],
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "virus",
            "merge_blastn.eval_{evalue}.cov_{coverage}.pid_{pident}.log",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script:
        "../scripts/merge_blastn.py"


##########################################################################
##########################################################################


rule blastn_human:
    input:
        contig=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "assemblies",
            "all_contigs.over3kb.nr.fasta",
        ),
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "human",
            "all_contigs.over3kb.nr.human_blastn.tsv",
        ),
    params:
        database=blast_database,
        outfmt="6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle",
        evalue="0.0001",
        max_len=2000000,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "human",
            "nt.human.blastn.outfmt6.log",
        ),
    conda:
        "../envs/blast.yaml"
    threads: 20
    shell:
        """
        blastn -task megablast -query {input.contig:q} -out {output.blast_out:q} \
               -db {params.database:q} -evalue {params.evalue} -outfmt {params.outfmt:q} \
               -word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -taxids 9606 \
               -min_raw_gapped_score 100 -perc_identity 90 -soft_masking true -max_target_seqs 10 \
               -num_threads {threads} -mt_mode 1 &> {log:q}
        """


##########################################################################
##########################################################################


rule merge_blastn_human:
    input:
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "human",
            "all_contigs.over3kb.nr.human_blastn.tsv",
        ),
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "assemblies",
            "all_contigs.over3kb.nr.fasta",
        ),
    output:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "databases",
            "contigs",
            "human_filtered",
            "all_contigs.over3kb.nr.filtered.fasta",
        ),   
    params:
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "human",
            "all_contigs.over3kb.nr.annotation.blastn.tsv",
        ),
        minimum_length=config["default_blast_option"]["length_min"],
        coverage=0.6,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "human",
            "all_contigs_blastn.filtered.human.log",
        ),
    resources:
        cpus=1,
    conda:
        "../envs/biopython.yaml"
    threads: 1
    script:
        "../scripts/remove_human.py"


##########################################################################
##########################################################################

rule blast_makedatabase_nucleotide_dereplication:
    input:
        fasta=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "assemblies",
            "all_contigs.over3kb.fasta",
        )
    output:
        multiext(
            os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "assemblies",
                "all_contigs.over3kb.fasta",
            ),
            ".ndb",
            ".not",
            ".ntf",
            ".nto",
            ".njs"
        )
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "dereplicated",
            "blast.makeblastdb.log",
        ),
    params:
        "-input_type fasta -blastdb_version 5 -parse_seqids"
    conda:
        "../envs/blast.yaml"
    shell:
        """
        makeblastdb -in {input.fasta:q} -dbtype nucl {params} -logfile {log:q} -out {input.fasta:q}
        """

##########################################################################
##########################################################################

rule blast_all_vs_all:
    input:
        contig=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "assemblies",
            "all_contigs.over3kb.fasta",
        ),
        blastdb=multiext(
            os.path.join(
                OUTPUT_FOLDER,
                "processed_files",
                "assemblies",
                "all_contigs.over3kb.fasta",
            ),
            ".ndb",
            ".not",
            ".ntf",
            ".nto",
            ".njs"
        ),
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "dereplicated",
            "all_contigs_VS_all_contigs.blastout",
        )
    params:
        database=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "assemblies",
            "all_contigs.over3kb.fasta",
        ),
        outfmt="6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle",
        evalue="1e-10",
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "dereplicated",
            "blast.log",
        ),
    conda:
        "../envs/blast.yaml"
    threads:
        20
    shell:
        """
        blastn -query {input.contig:q} -db {params.database:q} -out {output.blast_out:q} -outfmt {params.outfmt:q} -evalue {params.evalue} -num_threads {threads} > {log:q}
        """


##########################################################################
##########################################################################


rule keep_only_nr_contigs:
    input:
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "assemblies",
            "all_contigs.over3kb.fasta",
        ), 
        blast_file=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "dereplicated",
            "all_contigs_VS_all_contigs.blastout",
        ),
    output:
        tsv=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "dereplicated",
            "all_contigs_VS_all_contigs.blastout.tsv",
        ),
        contigs=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "assemblies",
            "all_contigs.over3kb.nr.fasta",
        ),
    params:
        evalue = 1e-20,
        coverage = 0.85,
        pident = 0.95,
        minimum_length=config["default_blast_option"]["length_min"],
    conda:
        "../envs/biopython.yaml"
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "dereplicated",
            "blast_dereplication.log",
        ),
    threads:
        1
    script:
        "../scripts/dereplicate_contigs_blastn.py"

##########################################################################
##########################################################################

rule blast_otu:
    input:
        contig=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "viral_contigs",
            "viral_contigs_over_3kb.fna",
        )
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "otu",
            "viral_contigs_over_3kb_all_VS_all.blastout",
        ) 
    params:
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend  sstart send evalue bitscore qlen slen",
        evalue="1e-5",
        max_len=2000000,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "otu",
            "blastn.log",
        ),
    conda:
        "../envs/blast.yaml"
    threads: 20
    shell:
        """
        makeblastdb -in {params.contig:q} -dbtype nucl &> {log:q}

        blastn -task megablast -query {input.contig:q} -out {output.blast_out:q} \
                -db {input.contig:q} -evalue {params.evalue} -outfmt {params.outfmt:q} \
                -num_threads {threads} -max_target_seqs 20000 &> {log:q}
        """


##########################################################################
##########################################################################


rule blastn_taxonomy:
    input:
        contig=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "viral_contigs",
            "viral_contigs_over_3kb.fna",
        ),
    output:
        blast_out=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "taxonomy",
            "viral_contigs_over_3kb.taxonomy.evalue_1e-10.{database}.blastn.outfmt6.txt",
        ),
    params:
        database=lambda wildcards: DB_DICT[wildcards.database]["path"],
        tmp_output=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "taxonomy",
            "tmp_{database}",
        ),
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore staxids stitle",
        evalue="1e-10",
        options_blast="-num_alignments 25000",
        max_len=2000000,
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "taxonomy",
            "viral_contigs_over_3kb.evalue_1e-10.{database}.blastn.outfmt6.log",
        ),
    conda:
        "../envs/blast.yaml"
    threads: 20
    shell:
        """
        if [[ $(stat -c "%s" {input.contig:q}) -gt $(stat -c "%s" {params.database:q}) ]]
        then
            blastn -task megablast -query {input.contig:q} -out {output.blast_out:q} \
                -db {params.database:q} -evalue {params.evalue} -outfmt {params.outfmt:q} \
                -num_threads {threads} -mt_mode 1 &> {log:q}
        else
            blastn -task megablast -query {input.contig:q} -out {output.blast_out:q} \
                -db {params.database:q} -evalue {params.evalue} -outfmt {params.outfmt:q} \
                -num_threads {threads} -mt_mode 0 &> {log:q}
        fi            
        """


##########################################################################
##########################################################################

rule extract_species_from_blasts:
    input:
        ictv_blasts=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "taxonomy",
            "viral_contigs_over_3kb.taxonomy.evalue_1e-10.ICTV.blastn.outfmt6.txt",
        ),
        refseq_blasts=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "taxonomy",
            "viral_contigs_over_3kb.taxonomy.evalue_1e-10.refseq_viral.blastn.outfmt6.txt",
        ),
        crass_blasts=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "taxonomy",
            "viral_contigs_over_3kb.taxonomy.evalue_1e-10.crassphage.blastn.outfmt6.txt",
        ),
        img_blasts=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "blast",
            "taxonomy",
            "viral_contigs_over_3kb.taxonomy.evalue_1e-10.IMG_VR.blastn.outfmt6.txt",
        ),
    output:
        outfile=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "taxonomy",
            "blast_table_species.tsv",
        )
    params:
        ictv_db=DB_DICT["ICTV"]["metadata"],
        refseq_db=DB_DICT["refseq_viral"]["metadata"],
        img_db=DB_DICT["IMG_VR"]["metadata"],
        anicalc=workflow.source_path("../scripts/anicalc.py"),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "taxonomy",
            "update_taxonomies.log",
        ),
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/update_taxonomies.py"

##########################################################################
##########################################################################
