graph snakemake_dag {
    graph[bgcolor=white, margin=0];
    node[shape=box, style=rounded, fontname=sans, fontsize=10, penwidth=2];
    edge[penwidth=2, color=grey];
	1[label = "MultiQC", color = "0.65 0.6 0.85", style="rounded"];
	2[label = "FastQC", color = "0.42 0.6 0.85", style="rounded"];
	8[label = "Fastp PE", color = "0.08 0.6 0.85", style="rounded"];
	13[label = "iPHoP: Predict host", color = "0.17 0.6 0.85", style="rounded"];
	16[label = "CheckV", color = "orangered", style="rounded"];
	17[label = "Merge positive hits", color = "crimson", style="rounded"];
	19[label = "BLASTn\ndatabase: ICTV", color = "0.31 0.6 0.85", style="rounded"];
	20[label = "Non-viral sequences", color = "0.57 0.6 0.85", style="rounded"];
	22[label = "BLASTn against human genome\nOnly keep non-human contigs", color = "0.62 0.6 0.85", style="rounded"];
	23[label = "Dereplication & Merge", color = "0.33 0.6 0.85", style="rounded"];
	24[label = "Only keep contigs bigger than 3kb", color = "0.49 0.6 0.85", style="rounded"];
	27[label = "metaSPAdes", color = "0.11 0.6 0.85", style="rounded"];
	33[label = "MEGAHIT", color = "0.43 0.6 0.85", style="rounded"];
	35[label = "geNomad", color = "0.12 0.6 0.85", style="rounded"];
	36[label = "BLASTn\ndatabase: RefSeq Viral", color = "0.31 0.6 0.85", style="rounded"];
	37[label = "BLASTn\ndatabase: Crassphage DB", color = "0.31 0.6 0.85", style="rounded"];
	38[label = "BLASTn\ndatabase: IMG/VR", color = "0.31 0.6 0.85", style="rounded"];
	39[label = "BLASTn\ndatabase: MGV", color = "0.31 0.6 0.85", style="dashed"];
	40[label = "BLASTn\ndatabase: GPD", color = "0.31 0.6 0.85", style="dashed"];
	41[label = "BLASTn\ndatabase: GVD", color = "0.31 0.6 0.85", style="dashed"];
	43[label = "BACPHLIP: Predict viral lifestyles", color = "0.32 0.6 0.85", style="rounded"];
	44[label = "Get the taxonomy based on\n1) ICTV\n2)RefSeq Viral\n3)IMG/VR\n4)Crassphage DB", color = "palevioletred", style="rounded"];
	45[label = "BLASTn for taxonomy\ndatabase: ICTV", color = "0.56 0.6 0.85", style="rounded"];
	46[label = "BLASTn for taxonomy\ndatabase: RefSeq Viral", color = "0.56 0.6 0.85", style="rounded"];
	47[label = "BLASTn for taxonomy\ndatabase: Crassphage DB", color = "0.56 0.6 0.85", style="rounded"];
	48[label = "BLASTn for taxonomy\ndatabase: IMG/VR", color = "0.56 0.6 0.85", style="rounded"];
	49[label = "geNomad annotate taxonomy", color = "0.51 0.6 0.85", style="rounded"];
	2 -- 1
	1 -- 8
	8 -- 27
	8 -- 33
	17 -- 16
	22 -- 35
	35 -- 17
	19 -- 17
	36 -- 17
	37 -- 17
	38 -- 17
	39 -- 17
	40 -- 17
	41 -- 17
	20 -- 19
	35 -- 20
	23 -- 22
	24 -- 23
	27 -- 24
	33 -- 24
	20 -- 36
	20 -- 37
	20 -- 38
	20 -- 39
	20 -- 40
	20 -- 41
	45 -- 44
	46 -- 44
	47 -- 44
	48 -- 44
	16 -- 45
	16 -- 46
	16 -- 47
	16 -- 48
	16 -- 49
	16 -- 13
	16 -- 43
}            
