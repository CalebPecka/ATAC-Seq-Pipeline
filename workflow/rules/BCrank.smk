# Load configuration file with sample data.
configfile: "../config/config.yaml"

rule BCrank:
	input:
		bam=expand("{sample}", sample=config["samples"]),
		seq="../results/{input.bam}.results/upstreamPeak_Sequences.fasta"
	output:
		directory("../results/{input.bam}.results/BCrankOutput/")
	params:
		seedSize=config["seedSize"],
		restartSize=config["restartSize"]
	conda:
		"../envs/R.yaml"
	log:
		"logs/{input.bam}.BCrank.log"
	shell:
		"Rscript scripts/BCrank.R {params.seedSize} {params.restartSize} {input.seq} {output}"