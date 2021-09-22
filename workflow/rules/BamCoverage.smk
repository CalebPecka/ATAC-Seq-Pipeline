# Load configuration file with sample data.
configfile: "../config/config.yaml"

rule BamCoverage:
	input:
		bam=expand("{sample}", sample=config["samples"])
	output:
		"../results/{sample}.results/bigWigCoverage.bw"
	conda:
		"../envs/pyGenomeTracks.yaml"
	log:
		"logs/{sample}.BamCoverage.log"
	shell:
		"bamCoverage -b {wildcards.sample} -o {output}"