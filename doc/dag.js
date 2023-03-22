digraph pipeline {
    // Input
    input [label="Multiple FASTQ files"];

    // Quality control step
    multiqc [label="MultiQC"];

    // Read trimming step
    fastp [label="Fastp"];

    // Assembly programs
    metaspades [label="MetaSPAdes"];
    megahit [label="MEGAHIT"];

    // Dereplication and merge results
    derep [label="Dereplication & Merge"];

    // Virus detection
    genomad [label="Genomad"];

    // Sequences not annotated as virus
    non_viral_seqs [label="Non-viral sequences"];

    // BLAST against different databases
    ictv_blast [label="ICTV BLAST"];
    refseq_viral_blast [label="RefSeq Viral BLAST"];
    imgvr_blast [label="IMG/VR BLAST"];

    // Merge positive hits with GenomAD
    merge_hits [label="Merge positive hits"];

    // CheckV
    checkv [label="CheckV"];

    // EDGES
    input -> multiqc;
    multiqc -> fastp;
    fastp -> metaspades;
    fastp -> megahit;
    metaspades -> derep;
    megahit -> derep;
    derep -> genomad;
    genomad -> non_viral_seqs;
    non_viral_seqs -> ictv_blast;
    non_viral_seqs -> refseq_viral_blast;
    non_viral_seqs -> imgvr_blast;
    ictv_blast -> merge_hits;
    refseq_viral_blast -> merge_hits;
    imgvr_blast -> merge_hits;
    genomad -> merge_hits;
    merge_hits -> checkv;
}