# Load configuration file with sample data.
configfile: "../config/config.yaml"

rule CreateRefGenome:
	output:
		"../config/HG38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
	log:
		"logs/CreateRefGenome.log"
	shell:
		"curl -LJO ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > {output}"

rule SummitFASTA:
	input:
		bam=expand("{sample}", sample=config["samples"]),
		refGenome={rules.CreateRefGenome.output},
		upstreamPeaks="../results/{input.bam}.results/upstreamPeaks.tsv"
	output:
		"../results/{input.bam}.results/upstreamPeak_Sequences.fasta"
	params:
		motifSize=config["motifSize"]
	conda:
		"../envs/R.yaml"
	log:
		"logs/{input.bam}.SummitFASTA.log"
	shell:
		"Rscript scripts/SummitFASTA.R {input.refGenome} {input.upstreamPeaks} {params.motifSize} {output}"