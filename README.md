# ATAC-Seq-Pipeline


# Requirements

  **R version >= 4.0.0 is required for this project. The correct R version is automatically installed in a new conda environment (ATAC) using the *initialization.sh* script. The *initialization.sh* script includes automated installation that prevents conflicts with package dependencies in R. If you would prefer to set up your own environment, you will have to modify *initialization.sh* and install packages for "Biostrings", "BCRANK", and "seqinr".**
  
  **This GitHub Pipeline includes an optional branch for track-based visualization of chromatin accessibility. To use this commands, you will have to install the following dependencies in the "ATAC" environment created by *initialization.sh*:**
  - Python >= 3.6
  - numpy >= 1.8.0 
  - scipy >= 0.17.0
  - py2bit >= 0.1.0
  - pyBigWig >= 0.3.4
  - pysam >= 0.8
  - matplotlib >= 3.1.1
  - intervaltree >= 2.1.0
  - hicmatrix >= 0.14
  - gffutils >= 0.9

After the depencies have been installed, the following commands will install the PyGenomeTracks visualization tool to your environment.

  - conda install -c bioconda -c conda-forge pygenometracks python=3.7
  - conda install -c bioconda deeptools

# Workflow Overview

![ImageOfWorkflow](https://github.com/CalebPecka/ATAC-Seq-Pipeline/blob/master/ATACseq_Visualization.png)

This repository assumes you have already processed your raw fastq files, and your data is formatted as a series of .bam files, one per sample. The entire pipeline can be run with the *macro.sh* script after performing initialization. If you encounter difficulties, documentation for each command can be seen in the sections below.

**Using the Macro Scripts**

In your main directory, run the *initialization.sh* script. This creates a conda environment for you which installs package depencies for R version 4.0 and downloads the HG38 human reference genome.

If you wish to visualize genomic tracks in future steps, use the following commands to install to your environment:
  - conda install -c bioconda -c conda-forge pygenometracks python=3.7
  - conda install -c bioconda deeptools
  
For each sample, run the *macro.sh* script with 3 arguments. The first being your input .bam file, the second being an output directory.
For example:

  -*sh macro.sh /data/fooFileLocation/3c39.bam /MACS2_out/*
  
  Make sure to include a "/" at the end of your directory reference, seen in the example above.
  
This command generally takes an hour or two to complete.

Once finished, your *tracks.ini* file can be used to visualize open chromatin regions side-by-side with bp gene locations. When you want to plot a region of the genome, use: *pyGenomeTracks --tracks tracks.ini --region chr15:2,500,000-2,800,000 - bigwig.png*, with variables adjusted to visualize the correct chromosome region. I recommend using a base pair range of ~300,000 so that the gene names are easily readible.

The *tracks.ini* file has many visual customization options. Here is some example code to give you an idea of what objects you might want to include: https://pygenometracks.readthedocs.io/en/latest/content/examples.html

# How to Use

Each of the programs, inputs, and outputs for this program can be visualized in the diagram at the top of this documentation. Input/output files are colored yellow, and scripts are colored blue. Some files are required inputs for multiple steps. If you choose to perform the pipeline manually, there is no set order when executing each script. Instead, it is easier to follow the arrows on the diagram, and ensure each required input exists in your directory before executing the next script. Below you will find two sections for **Scripts** and **Input/Output Files**, each containing descriptions and requirements for each object in the diagram. To run the pipeline, follow each blue script box, and execute their commands included below.

# Scripts

## MACS2
Command: **macs2 callpeak -t $1 -g hs -f BAM -p 0.05 --seed 0 --bdg --outdir .**
*For the command above, $1 refers to the .bam file you're processing.*

MACS2 is our peakcalling method. Open chromatin regions can be identified by regions where reads have piled up more than the background coverage.

  --'t' is the parameter used for our input file.
  --'g hs' identifies that we're using a human sequence.
  --'f' identifies our input file type, in this case bam.
  --'p' is our pvalue for identifying if the piled reads are significantly more accessible than the background coverage. These results are tabulated in a narrowPeaks file for future usage.
  --'seed' is a randomly generated started value which can be used to create reproducible results.
  --'bdg' indicates that we'd like to also create a bedgraph file, which we later convert to a bigwig file, a common file type for visualizing open chromatin regions.
  --'outdir' is our ouput directory. By including a dot, '.', we are telling the computer to send all output files to the current directory.

## BdgToBigWig
**bamCoverage -b $1 -o bigWig_coverage.bw**
*For the command above, $1 refers to the .bam file you're processing.*

BamCoverage is a method of converting a bedgraph file to a bigwig file. Our visualization tool, PyGenomeTracks, requires an input file BigWig type. If you don't wish to visualize your results, this step can be ignored. 

  --'b' is the parameter used for our input file.
  --'o' is your output directory. In the example command above, we are naming our file to "bigWig_coverage.bw" using the "-o" output parameter.
  
## Make Tracks
**make_tracks_file --trackFiles bigWig_coverage.bw NA_peaks.narrowPeak GFFgenes.bed -o tracks.ini**

"Make Tracks" layers a series of rows together for the PyGenomeTracks visualization. Using chromosomal positioning, the "tracks.ini" output stores the position of relevant information relative to the genome. The resulting "tracks.ini" file is easily editable, see documentation here: https://pygenometracks.readthedocs.io/en/latest/content/examples.html. Only one parameter is used, 'trackFiles'. Each file you include AFTER the parameter is layered on top of the data from the previous file. In this case, our visualization includes the bigWig coverage (peak visualization), narrowPeaks (bp location of peak visualized as a box plot), and a bed file for the start/stop regions of human genes. The command can be modified to include other results throughout the pipeline, and only requires 1 or more file to be included as a parameter.

## PyGenomeTracks
**pyGenomeTracks --tracks tracks.ini --region chr1:1000000-4000000 -o image.png**

PyGenomeTracks is the visualization command for the pipeline. It will plot a layered track of the data from the input "track.ini".

  --'tracks' is the parameter to indicate your input file. Must be type .ini.
  --'region' indicates the chromosomal location you want to visualize. It is formated where the abbreviated chromosome name is included first, followed by the base-pair region, separating start and stop position by a hyphen. The chromosomal range shown above starts at 1 million and ends at 4 million. 
  --'o' is your output directory/file name.

## NarrowPeakSummitTracker.R
**Rscript NarrowPeakSummitTracker.R**

The goal of this function is to locate every instance of an accessible chromatin region that is within a 1000 base upstream region of genes in the human genome (likely locations to include a transcription factor binding site). The *GFFgenes.bed* file is used to find start locations of genes, and find summits from the MACS2 output *(NA_peaks.narrowPeak)* that are within the desired range of the gene. The output file *upstreamPeaks.tsv* stores information about the accessible chromatin regions location, including chromosome number, and start/stop base-pair (bp) numbers. *upstreamPeaks.tsv* also stored valuable information such as the name of genes that were discovered downstream, and the q-score of the peak region where it was identified. The q-score helps us estimate how confident the program was that it identified a chromatin accessible peak region.

## SummitFASTA.R
**Rscript SummitFASTA.R**

Once we've collected the LOCATION of potential transcription factor binding sites (TFBS) from *NarrowPeakSummitTracker.R*, we need to determine the sequence of the TFBS. Problematically, the real TFBS can be located anywhere within the accessible chromatin region, which can sometimes be thousands of bases long. At the peak location of each accessible chromatin region (defined by the highest density of sequence fragment overlap), the surrounding nucleotide sequence is stored in a new FASTA file called *upstreamPeak_sequences.tsv*. The size of the nucleotide search space can be modified by the user. By default, the script searches for a 100 bp region surrounding each summit.
A value of 100 is consistent with the default requirements of other motif-searching tools like the tools used by MEME Suite. This value can be modified by changing the MotifSize variable within the *SummitFASTA.R* script.

## BCrank.R
**Rscript BCrank.R**

BCRANK is tool that, when provided a list of FASTA sequences, determines frequently repeated motif sequences. In our case, repeat motif sequences are likely to be important structural components to a TFBS. BCRank requires the input to be ordered according to confidence level. Our script automatically sorts the data according to q-score, the confidence level provided by MACS2. By default, the BCrank program is repeated 12 times with different randomly generated seeds, each producing 100 motifs. The motifs are stored in a new directory *BCrankOutput/*.

## TomTom
