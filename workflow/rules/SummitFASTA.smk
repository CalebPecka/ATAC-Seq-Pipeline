# Load configuration file with sample data.
configfile: "../config/config.yaml"

perfectConfig=config["perfectCondaReplication"]
if perfectConfig == "TRUE":
	condaLocation="../envs/R.yaml"
else:
	condaLocation="../envs/Rconfig.yaml"

rule CreateRefGenome:
	output:
		"../config/HG38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
	log:
		"logs/CreateRefGenome.log"
	shell:
		"curl -LJO ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > {output}"

rule SummitFASTA:
	input:
		refGenome={rules.CreateRefGenome.output},
		upstreamPeaks="../results/{sample}.results/upstreamPeaks.tsv"
	output:
		"../results/{sample}.results/upstreamPeak_Sequences.fasta"
	params:
		motifSize=config["motifSize"]
	conda:
		condaLocation
	log:
		"logs/{sample}.SummitFASTA.log"
	shell:
		"Rscript scripts/SummitFASTA.R {input.refGenome} {input.upstreamPeaks} {params.motifSize} ../results/{wildcards.sample}.results/upstreamPeak_Sequences.fasta"
