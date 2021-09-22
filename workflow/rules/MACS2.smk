# Load configuration file with sample data.
configfile: "../config/config.yaml"

perfectConfig=config["perfectCondaReplication"]
if perfectConfig == "TRUE":
	condaLocation="../envs/macs.yaml"
else:
	condaLocation="../envs/macsConfig.yaml"

rule MACS2:
	input:
		bam=expand("{sample}", sample=config["samples"])
	output:
		NApeak="../results/{sample}.results/NA_peaks.narrowPeak"
	params:
		dirName="../results/{sample}.results/",
		genomeType=config["genomeType"],
		peakFileType=config["peakFileType"],
		peakPVal=config["peakPVal"],
		seed="0"
	conda:
		"../envs/macs.yaml"
	log:
		"logs/{sample}.MACS2.log"
	shell:
		"macs2 callpeak -t /{sample} -g {params.genomeType} -f {params.peakFileType} -p {params.peakPVal} --seed {params.seed} --outdir {params.dirName}"
