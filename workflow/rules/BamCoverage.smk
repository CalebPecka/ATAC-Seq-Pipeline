rule BamCoverage:
	input:
		bam=expand("{sample}", sample=config["samples"])
	output:
		"../results/{input.bam}.results/bigWigCoverage.bw"
	conda:
		"../envs/pyGenomeTracks.yaml"
	log:
		"logs/{input.bam}.BamCoverage.log"
	shell:
		"bamCoverage -b {input.bam} -o {output}"