# ATAC-Seq-Pipeline


# Project Overview
This pipeline is designed to take a set of input BAM files and perform ATAC-seq/chip-seq analysis on the data. The pipeline will determine accessible regions of chromatin, motif sequences within those regions, indicative of transcription factor binding sites (TFBS), compare those TFBS to known TFBS from the JASPAR database, align those sequences to the genome, create tracks to dynamically visualize the data, and perform statistical analysis to verify our results.

**The entire pipeline is run using the Python tool "Snakemake".** Snakemake automatically determines missing input files, conda environment requirements/dependencies, and more. Below is a diagram of the files in this directory and there importance. **Most notably,** "_example_files_" contains a subset of the formated data you should expect from the inputs and outputs of this pipeline. "Config" contains a series of parameters and the species reference genome, allowing the user to quickly test how different parameters modify the final results. "Results" contains all of the final outputs in folders separated for each BAM file input into the pipeline. "Workflow" contains the scripts and functions used to run the pipeline in a streamline manner. We do not advise modifying the "Workflow" folder. Instead, rely on the configuration files and Snakemake parameters to modify the pipeline to your needs, if possible.

![ImageOfProjectDirectory](https://github.com/CalebPecka/ATAC-Seq-Pipeline/blob/master/__graphics__/Elongated_Pipeline_Workflow.png)

# Requirements

**Snakemake is the only required installation to run this pipeline.** Configuration files within this project will automatically be installed using Snakemake based on the dependencies that existed during the creation of this project. In short, we automatically install all of the requirements for you, exactly when each script calls on a dependency. The following commands will install Snakemake to a separate environment for you:
  - conda install -n base -c conda-forge mamba
  - conda activate base
  - mamba create -c conda-forge -c bioconda -n snakemake snakemake
  - conda activate snakemake
  - snakemake --help

If you have any issues with installation, see more info at: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

# Workflow Overview

![ImageOfWorkflow](https://github.com/CalebPecka/ATAC-Seq-Pipeline/blob/master/__graphics__/ATACseq_Block_Diagram.png)

This repository assumes you have already processed your raw fastq files, and your data is formatted as a series of .bam files, one per sample. Afterwards we use MACS2 to perform _Peak Calling_, a technique that identifies accessible regions of chromatin in DNA. From there we perform _Motif Identification_ to determine putative transcription factor binding sites with these peaks. JASPAR is used as a database of known transcription factors for _Motif Comparison_. _Site Matching_ is performed to align these sequences against our FASTA sequences to determine which transcription factors bind upstream to which genes. Finally, all of these files can be incorperated together into a _Track Based Visualization_ like the example seen below:

![Track](https://github.com/CalebPecka/ATAC-Seq-Pipeline/blob/master/__graphics__/SLCO4A1.3c39.genomeTrack.png)

The entire pipeline can be run using Snakemake. If you encounter difficulties, documentation for each command can be seen in the sections below.

# How to Use

**To run the entire pipeline, use the following commands:**
 - cd workflow
 - sh scripts/createCondaConfigurations.sh
 - snakemake "../results/{A,B,C,D,E}.results/matchingSiteTable.csv" --use-conda --cores N
 
Where A, B, C, D, and E are any number of BAM files you wish to analyze, delimited by commas. _--cores N_ is used to indicate the number of cores you wish to use. Snakemake will automatically determine which steps of the pipeline can be performed based on the files that have been created, and start a new processes using new cores to optimize the speed of the workflow.

The pipeline must be run from within the "workflow" directory (where Snakefile is located). The shell scripts (scripts/createCondaConfigurations) will download all necessary dependencies in a series of conda environments for you. We've precreated several conda configuration files that perfectly replicate the original conditions of the experiment on macOS. If you would rather try those dependencies, go into "config/config.yaml" and change the "perfectCondaReplication" parameter to "TRUE". Additional parameters for the pipeline can be modified in the "config/config.yaml" file, including file types, p-value thresholds, gtf files, and more. **The location of sample BAM files {A,B,C,D,E}, should also be modified in this location before running the pipeline.**

Once finished, your *tracks.ini* file can be used to visualize open chromatin regions side-by-side with bp gene locations. When you want to plot a region of the genome, use: *pyGenomeTracks --tracks tracks.ini --region chr15:2,500,000-2,800,000 - bigwig.png*, with variables adjusted to visualize the correct chromosome region. I recommend using a base pair range of ~300,000 so that the gene names are easily readible.

The *tracks.ini* file has many visual customization options. Here is some example code to give you an idea of what objects you might want to include: https://pygenometracks.readthedocs.io/en/latest/content/examples.html

![ImageOfWorkflow](https://github.com/CalebPecka/ATAC-Seq-Pipeline/blob/master/__graphics__/Color_Coded_Visualization.png)

Each of the programs, inputs, and outputs for this program can be visualized in the diagrams seen above. Each step of the pipeline has been color coded based on the goal we are trying to achieve (See matching color diagram in "Workflow Overview"). If you choose to perform the pipeline manually, there is no set order when executing each script. Instead, it is easier to follow the arrows on the diagram, and ensure each required input exists in your directory before executing the next script (_Note that Snakemake does this for you!_). In the section below, you will find two sections for **Scripts** and **Input/Output Files**, each containing explanations and requirements for each object in the diagram. The diagram below is a copy of the diagram above, except Input/output files have been colored yellow, scripts are colored blue, and databases are colored red.

# Scripts

![ImageOfWorkflow](https://github.com/CalebPecka/ATAC-Seq-Pipeline/blob/master/__graphics__/Programatic_Visualization.png)

## MACS2
Command: **macs2 callpeak -t {sample} -g {params.genomeType} -f {params.peakFileType} -p {params.peakPVal} --seed {params.seed} --outdir {params.dirName}**

*For the command above, {sample} refers to the location of the BAM files you're processing.*

MACS2 is our peakcalling method. Open chromatin regions can be identified by regions where reads have piled up more than the background coverage.

  --'t' is the parameter used for our input file.
  
  --'g hs' identifies that we're using a human sequence.
  
  --'f' identifies our input file type, in this case bam.
  
  --'p' is our pvalue for identifying if the piled reads are significantly more accessible than the background coverage. These results are tabulated in a narrowPeaks file for future usage.
  
  --'seed' is a randomly generated started value which can be used to create reproducible results.
  
  --'bdg' indicates that we'd like to also create a bedgraph file, which we later convert to a bigwig file, a common file type for visualizing open chromatin regions. The existing pipeline uses a slightly different approach, and we will instead create a bigwig file using the BamCoverage command from DeepTools. 
  
  --'outdir' is our ouput directory. By including a dot, '.', we are telling the computer to send all output files to the current directory.
  
The final output is a narrowPeak file listing a series of peaks and some information about the significance/relative fold change of each peak. An example of the output can be found here: https://github.com/CalebPecka/ATAC-Seq-Pipeline/blob/master/__example_files__/NA_peaks%5BEXAMPLE%5D.narrowPeak

## NarrowPeakSummitTracker.R
Command: **Rscript scripts/NarrowPeakSummitTracker.R {input.gtf} {input.NApeak} {output}**

The R scripts in this pipeline function off of an ordered series of arguments that is passed to the script.

  --'{input.gtf}' is a precreated file in the GitHub repository that lists the locations of genes relative to the Hg38 _Homo sapiens_ reference genome.
  
  --'{input.NApeak}' identifies that we're using a human sequence.
  
  --'{output}' is the location where the output file will be saved..

The goal of this function is to locate every instance of an accessible chromatin region that is within a 1000 base upstream region of genes in the human genome (likely locations to include a transcription factor binding site). The *GFFgenes.bed* file is used to find start locations of genes, and find summits from the MACS2 output *(NA_peaks.narrowPeak)* that are within the desired range of the gene. The output file *upstreamPeaks.tsv* stores information about the accessible chromatin regions location, including chromosome number, and start/stop base-pair (bp) numbers. *upstreamPeaks.tsv* also stored valuable information such as the name of genes that were discovered downstream, and the q-score of the peak region where it was identified. The q-score helps us estimate how confident the program was that it identified a chromatin accessible peak region.

An example of the output can be found here: https://github.com/CalebPecka/ATAC-Seq-Pipeline/blob/master/__example_files__/upstreamPeaks%5BEXAMPLE%5D.tsv

## SummitFASTA.R
Command: **curl -LJO ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz > {output}**
  **Rscript scripts/SummitFASTA.R {input.refGenome} {input.upstreamPeaks} {params.motifSize} {output}**

This command requires 2 steps. First, we check if you currently have the Hg38 reference genome build installed. If not, it is downloaded to a separate directory in the Config folder of this repository. Only then do we call the R script for the next step.

  --'{input.refGenome}' is the location of the reference genome we just installed.
  
  --'{input.upstreamPeaks}' is the output of NarrowPeakSummitTracker, detailing the locations of peaks that are likely to regulate downstream genes.
  
  --'{params.motifSize}' is the size of the FASTA sequences surrounding the summit location of the peaks from the input.upstreamPeaks file. For example, if we were to pick a value of 100, the script would output a FASTA sequence 100 base-pairs long surrounding the highest point of chromatin accessible from every row in the input.upstreamPeaks file.
  
  --'{output}' is the location where the output file will be saved..
  
Once we've collected the LOCATION of potential transcription factor binding sites (TFBS) from *NarrowPeakSummitTracker.R*, we need to determine the sequence of the TFBS. Problematically, the real TFBS can be located anywhere within the accessible chromatin region, which can sometimes be thousands of bases long. At the peak location of each accessible chromatin region (defined by the highest density of sequence fragment overlap), the surrounding nucleotide sequence is stored in a new FASTA file called *upstreamPeak_sequences.tsv*. The size of the nucleotide search space can be modified by the user. By default, the script searches for a 100 bp region surrounding each summit.
A value of 100 is consistent with the default requirements of other motif-searching tools like the tools used by MEME Suite. This value can be modified by changing the MotifSize variable within the *SummitFASTA.R* script.

An example of the output can be found here: https://github.com/CalebPecka/ATAC-Seq-Pipeline/blob/master/__example_files__/upstreamPeak_Sequences%5BEXAMPLE%5D.fasta

## BCrank.R
Command: **Rscript scripts/BCrank.R {params.seedSize} {params.restartSize} {input.seq} {output}**

  --'{params.seedSize}' is used to determine the number of motifs generated for each sample. The total number of motifs equals {params.seedSize} multiplied by {params.restartSize}. **This value has a maximum of 20**.
  
  --'{params.restartSize}' is used to determine the number of motifs generated for each sample. The total number of motifs equals {params.seedSize} multiplied by {params.restartSize}. For example, a seedSize of 3 and a restartSize of 5 will output __3*5=15 seeds__.
  
  --'{input.seq}' is a FASTA file output from SummitFASTA.R.
  
  --'{output}' is the location where the output file will be saved..

BCRANK is tool that, when provided a list of FASTA sequences, determines frequently repeated motif sequences. In our case, repeat motif sequences are likely to be important structural components to a TFBS. BCRank requires the input to be ordered according to confidence level. Our script automatically sorts the data according to q-score, the confidence level provided by MACS2. By default, the BCrank program is repeated 12 times with different randomly generated seeds, each producing 100 motifs. The motifs are stored in a new directory *BCrankOutput/*.

An example of the output can be found here: https://github.com/CalebPecka/ATAC-Seq-Pipeline/tree/master/__example_files__/BCrankOutput%5BEXAMPLE%5D

## BCrankProcessing.R
Command: **Rscript scripts/BCrankProcessing.R {params.seedSize} {input.directory} {input.seq} {input.NApeak} {output}**

  --'{params.seedSize}' is used to determine the number of motifs generated for each sample. The total number of motifs equals {params.seedSize} multiplied by {params.restartSize}. This value should be the same as used in BCrank.R.
  
  --'{input.directory}' is the directory of BCrank motifs output from BCrank.R.
  
  --'{input.seq}' is a FASTA file output from SummitFASTA.R.
  
  --'{input.seq}' is the NarrowPeak data output from MACS2.
  
  --'{output}' is the location where the output file will be saved..

The results from the BCrankOutput/ directory are useless by themselves. The BCrankProcessing.R script concatenates and subsets the files within the directory to obtain a list of all 1200 motifs that *BCrank.R* previously identified. If you modify the BCrankOutput/ directory in any way, it is possible you can interfere with the programming in BCrankProcessing.R. Once the motifs have been extracted, the script searches for instances of every motif within the *upstreamPeak_sequences.tsv* file. *upstreamPeak_sequences.tsv* is a necessary file to run BCrankProcessing.R, but the file has been removed from the GitHub repository diagram to decrease clutter. BCrankProcessing.R also verifies that each matching site is location within the region of a narrow peak. This calculation is performed by comparing the genomic location data with the previously created file, *NA_peaks.narrowPeak*. The output file *matchingSites.csv* stores any information necessary for remaining scripts.

An example of the outpu can be found here: https://github.com/CalebPecka/ATAC-Seq-Pipeline/blob/master/__example_files__/matchingSiteTable%5BEXAMPLE%5D.csv

## Analysis.R
Command: **...**

Finally, the results can be analyzed. Analysis should be performed uniquely according to the design of each experiment. In our case, researchers were interested in verifying the discovery of TFBS by comparing the results to known differential gene expression patterns. Example files have been included in this project to illustrate the creation of a chi-squared analysis, as well as a normal curve distribution for the occupancy of TFBS at each gene. 

## BamCoverage
Command: **bamCoverage -b {wildcards.sample} -o {output}**

*For the command above, {wildcards.sample} refers to the location of the BAM files you're processing.*

BamCoverage is a method of converting a bedgraph file to a bigwig file. Our visualization tool, PyGenomeTracks, requires an input file BigWig type. If you don't wish to visualize your results, this step can be ignored. 

  --'b' is the parameter used for our input file.
  
  --'o' is your output directory. In the example command above, we are naming our file to "bigWig_coverage.bw" using the "-o" output parameter.
  
## Make Tracks
Command: **make_tracks_file --trackFiles bigWig_coverage.bw NA_peaks.narrowPeak GFFgenes.bed -o tracks.ini**

"Make Tracks" layers a series of rows together for the PyGenomeTracks visualization. Using chromosomal positioning, the "tracks.ini" output stores the position of relevant information relative to the genome. The resulting "tracks.ini" file is easily editable, see documentation here: https://pygenometracks.readthedocs.io/en/latest/content/examples.html. Only one parameter is used, 'trackFiles'. Each file you include AFTER the parameter is layered on top of the data from the previous file. In this case, our visualization includes the bigWig coverage (peak visualization), narrowPeaks (bp location of peak visualized as a box plot), and a bed file for the start/stop regions of human genes. The command can be modified to include other results throughout the pipeline, and only requires 1 or more file to be included as a parameter.

_We've already created a version of the tracks.ini output for you that's compatible with all samples from the pipeline. This command is not rerun during the pipeline._

## PyGenomeTracks
Command: **pyGenomeTracks --tracks tracks.ini --region chr1:1000000-4000000 -o image.png**

PyGenomeTracks is the visualization command for the pipeline. It will plot a layered track of the data from the input "track.ini".

  --'tracks' is the parameter to indicate your input file. Must be type .ini.
  
  --'region' indicates the chromosomal location you want to visualize. It is formated where the abbreviated chromosome name is included first, followed by the base-pair region, separating start and stop position by a hyphen. The chromosomal range shown above starts at 1 million and ends at 4 million. 
  
  --'o' is your output directory/file name.
  
## TomTom (Part 1)
Command: **meme upstreamPeak_Sequences.fasta -dna -oc . -nostatus -time 18000 -mod zoops -nmotifs 100 -minw 6 -maxw 25 -objfun classic -revcomp -markov_order 0**

MEME is another powerful tool for identifying motif sequences. The tool is hosted under the website: https://meme-suite.org/meme/tools/meme. MEME can also be ran using a command line, see software versions here: https://meme-suite.org/meme/doc/download.html. After running the command above, you will receive a series MEME.txt file which can be converted using the MEMEtxtToPWMconverter.R file. Simply run the Rscript after modifying the input variable within the script. This will return a series of position weighted matrices of predicted motifs, similar to BCrank. In BCrank, we extracted the consensus sequences from the position weighted matrices. 

## TomTom (Part 2)
Command: **tomtom -o dir BCRANK_motifs.meme jasperMemeFormatted.meme**

TomTom is tool that compares a list of motifs against a database of known motifs. The tool is hosted under the website: https://meme-suite.org/meme/doc/tomtom.html?man_type=web, and can be installed from the same software versions in "TomTom (Part 1)". The steps from TomTom are not automatically ran by the pipeline, and are instead suggested as a useful tool for understand how predicted motifs compare to known transcription factor motifs.
  
  
# Input/Output Files

An example of each input/output file has been included in the *Example Inputs and Outputs* directory. By removing the "[EXAMPLE]" text from each file name and moving the files into the main directory, you can follow along and run each script to see how they function.

## Narrow Peaks
Real File Name: **NA_peaks.narrowPeak**

A file describing the genomic location of chromatin accessible regions. Chromatin accessible regions are discovered by using DNA fragment-based pile-ups, producing "peaks". You can read more about fragment pile-ups and peak-calling here: https://hbctraining.github.io/Intro-to-ChIPseq/lessons/05_peak_calling_macs.html. The MACS2 narrowPeak file type includes 10 columns. Each row designates a separate chromatin region. Column 1 describes the chromosomal location. Column 2 describes the start location of the peak. Column 3 describes the stop location of the peak. Column 10 is the summit location in terms of base pairs from the start location. In other words, the value of column 10 plus the value of column 2 is the summit location for the peak.

## GTF Genes
Real File Name: **GTFgenes.bed**

A file describing the genomic location of genes. Genes are identified according to each transcript signature (ENST), one ENST for each row. The file includes a header, and **this header must have the exact column names or the pipeline will fail.** Column 1 describes the chromosomal location. Column 2 describes the start location of the ENST. Column 3 describes the stop location of the ENST. Column 4 is the name of each ENST. Later steps of the program will remove any data contained AFTER the period in each value of the GeneName column, (i.e. *ENST00000450305.2* will be converted to *ENST00000450305*).

## Upstream Peaks
Real File Name: **upstreamPeaks.tsv**

A file containing the identity of peak summits located within a 1000 base upstream region of genes. Regions 1000 bases upstream of genes are likely to contain transcription factor binding sites (TFBS), locations that regulate the nearby gene. We expect that the accessible chromatin regions defined by each row of the upstreamPeaks.tsv file contain TFBS. Column 1 describes the chromosomal location. Column 2 describes the ENST that is being regulated (found 1000 bases downstream of peak). Column 3 describes the base pair position of the summit. Column 4 describes the Q score of the summit, a representation of the quality of peak. The Q score is a product of MACS2.

## Hg38 Chromosome Sequences
Real File Name: **HG38/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna**

A file containing the human genome in FASTA format. This is downloaded in the *initialization.sh* script so we can collect FASTA sequences from regions that are chromatin accessible.

## Upstream Peak Sequences
Real File Name: **upstreamPeak_sequences.tsv**

A FASTA file of the DNA region surrounding the peaks identified from *upstreamPeaks.tsv*. The header for each sequence is colon delimited, separating the genomic location from the name of each ENST.

## BCrankOutput/
Real File Name: **BCrankOutput/**

A directory containing motif sequences (repeat sequences of DNA that are speculated to be TFBS). The default pipeline runs BCrank with 12 different seeds, identifying 100 motifs in each seed. The outputs are processed and sent to the BcrankOutput. The number preceding each file name in the BCrankOutput directory is used to distinguish each seed trial. The directory includes two types of files, .txt and .meme. Each file type includes the same set of 100 motifs, but formatted for different programs. TomTom only accepts MEME formatted files for position weighted matrices, and many FASTA search tools only accept consensus sequences (.txt).

## Matching Site Table
Real File Name: **matchingSiteTable.csv**

A file containing instances of each motif in the *upstreamPeak_sequences.tsv* file. Column 1 describes the consensus sequence. Column 2 describes the chromosomal location of the sequence. Column 3 describes describes the ENST that is found 1000 bases downstream of the peak. Column 4 describes the start location of the sequence. Column 5 describes the end location of the sequence. Column 6 describes whether or not the sequence is contained within the peak region of the narrowPeak results from MACS2. This is an error check, and no matching sites are included if they do not have a peak region associated with their chromosomal location.

