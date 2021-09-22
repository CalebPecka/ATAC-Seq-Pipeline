# Load configuration file with sample data.
configfile: "../config/config.yaml"

perfectConfig=config["perfectCondaReplication"]
if perfectConfig == "TRUE":
	condaLocation="../envs/memeSuite.yaml"
else:
	condaLocation="../envs/memeSuiteConfig.yaml"

rule MotifComparison:
	input:
		newDirectory="../results/{sample}.results/BCrankOutput/"
	output:
		directory("../results/{sample}.results/FimoOutput/")
	params:
		validatedList=config["TFBSreference"]
	conda:
		condaLocation
	log:
		"logs/{sample}.MotifComparison.log"
	shell:
		"mkdir -p {output} ;"
		"awk 'FNR==1 && NR!=1{{next;}}{{print}}' ../results/{wildcards.sample}.results/BCrankOutput/*.meme > {output}/allBCrank.meme ;"
		"tomtom -norc -oc {output} -verbosity 1 -min-overlap 5 {output}/allBCrank.meme {params.validatedList}"