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
        tmp_output=os.path.join(
            OUTPUT_FOLDER, "processed_files", "blast", "virus", "tmp"
        ),
        outfmt="6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle",
        evalue="{evalue}",
        options_blast="-num_alignments 25000",
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "virus",
            "all_contigs.nr.evalue_{evalue}.{database}.blastn.outfmt6.log",
        ),   
    resources:
        cpus=5,
    conda:
        "../envs/blast.yaml"
    threads: 10
    script:
        "../scripts/blastn_wrapper.py"


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
        tmp_output=os.path.join(
            OUTPUT_FOLDER, "processed_files", "blast", "human", "tmp"
        ),
        outfmt="6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle",
        evalue="0.0001",
        options_blast=(
                    "-word_size 28 -best_hit_overhang 0.1 -best_hit_score_edge 0.1 -dust yes -taxids 9606 "
                    "-min_raw_gapped_score 100 -perc_identity 90 -soft_masking true -max_target_seqs 10 "
        ),
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "human",
            "all_contigs.nt.human.blastn.outfmt6.log",
        ),
    resources:
        cpus=5,
    conda:
        "../envs/blast.yaml"
    threads: 20
    script:
        "../scripts/blastn_wrapper.py"

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

rule blast_all_vs_all:
    input:
        contig=os.path.join(
            OUTPUT_FOLDER,
            "processed_files",
            "assemblies",
            "all_contigs.over3kb.fasta",
        )
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
        tmp_output=os.path.join(
            OUTPUT_FOLDER, "processed_files", "blast", "dereplicated", "tmp"
        ),
        outfmt="6 qseqid sseqid pident length qlen slen evalue qstart qend sstart send stitle",
        evalue="1e-20",
        options_blast="",
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
    script:
        "../scripts/blastn_wrapper.py"


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
        blast_out=os.path.join(
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
    conda:
        "../envs/biopython.yaml"
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "dereplicated",
            "blast._dereplicationlog",
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
        database=os.path.join(
            OUTPUT_FOLDER,
            "results",
            "viral_contigs",
            "viral_contigs_over_3kb.fna",
        ),
        tmp_output=os.path.join(
            OUTPUT_FOLDER, "processed_files", "blast", "otu", "tmp"
        ),
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend  sstart send evalue bitscore qlen slen",
        evalue="1e-5",
        options_blast="-task megablast -max_target_seqs 20000",
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "otu",
            "blast.log",
        ),
    conda:
        "../envs/blast.yaml"
    threads:
        20
    script:
        "../scripts/blastn_wrapper.py"


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
            OUTPUT_FOLDER, "processed_files", "blast", "taxonomy", "tmp"
        ),
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore staxids stitle",
        evalue="1e-10",
        options_blast="-num_alignments 25000",
    log:
        os.path.join(
            OUTPUT_FOLDER,
            "logs",
            "blast",
            "taxonomy",
            "viral_contigs_over_3kb.evalue_1e-10.{database}.blastn.outfmt6.log",
        ),   
    resources:
        cpus=5,
    conda:
        "../envs/blast.yaml"
    threads: 10
    script:
        "../scripts/blastn_wrapper.py"


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
        anicalc="../scripts/anicalc.py",
    conda:
        "../envs/biopython.yaml"
    script:
        "../scripts/update_taxonomies.py"

##########################################################################
##########################################################################
