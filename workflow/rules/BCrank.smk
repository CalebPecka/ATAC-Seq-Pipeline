# Load configuration file with sample data.
configfile: "../config/config.yaml"

perfectConfig=config["perfectCondaReplication"]
if perfectConfig == "TRUE":
	condaLocation="../envs/R.yaml"
else:
	condaLocation="../envs/Rconfig.yaml"

rule BCrank:
	input:
		seq="../results/{sample}.results/upstreamPeak_Sequences.fasta"
	output:
		directory("../results/{sample}.results/BCrankOutput/")
	params:
		seedSize=config["seedSize"],
		restartSize=config["restartSize"]
	conda:
		condaLocation
	log:
		"logs/{sample}.BCrank.log"
	shell:
		"mkdir -p {output} ;"
		"Rscript scripts/BCrank.R {params.seedSize} {params.restartSize} {input.seq} ../results/{wildcards.sample}.results/BCrankOutput/"
