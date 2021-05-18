# Load configuration file with sample data.
configfile: "../config/config.yaml"

rule MACS2:
	input:
		bam=expand("{sample}", sample=config["samples"])
	output:
		NApeak="../results/{input.bam}.results/NA_peaks.narrowPeak"
	params:
		dirName="../results/{input.bam}.results/",
		genomeType=config["genomeType"],
		peakFileType=config["peakFileType"],
		peakPVal=config["peakPVal"],
		seed="0"
	conda:
		"../envs/macs.yaml"
	log:
		"logs/{input.bam}.MACS2.log"
	shell:
		"macs2 callpeak -t {input.bam} -g {params.genomeType} -f {params.peakFileType} -p {params.peakPVal} --seed {params.seed} --outdir {params.dirName}"
