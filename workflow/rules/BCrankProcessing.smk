# Load configuration file with sample data.
configfile: "../config/config.yaml"

perfectConfig=config["perfectCondaReplication"]
if perfectConfig == "TRUE":
	condaLocation="../envs/R.yaml"
else:
	condaLocation="../envs/Rconfig.yaml"

rule BCrankProcessing:
	input:
		directory="../results/{sample}.results/BCrankOutput/",
		seq="../results/{sample}.results/upstreamPeak_Sequences.fasta",
		NApeak="../results/{sample}.results/NA_peaks.narrowPeak"
	output:
		"../results/{sample}.results/matchingSiteTable.csv"
	params:
		seedSize=config["seedSize"]
	conda:
		condaLocation
	log:
		"logs/{sample}.BCrankProcessing.log"
	shell:
		"Rscript scripts/BCrankProcessing.R {params.seedSize} {input.directory} {input.seq} {input.NApeak} ../results/{wildcards.sample}.results/matchingSiteTable.csv"
