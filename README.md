# ATAC-Seq-Pipeline


# Requirements

  **- R version >= 4.0.0 is required for this project, however it is automatically installed in a new conda environment using the *initialization.sh* script. This is done to prevent conflicts with package dependencies. If you would prefer to set up your own environment, you will have to modify the *initialization.sh* and install packages for "Biostrings", "BCRANK", and "seqinr".**
  
  **In the same environment, install:**
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

**Conda Environment**

PyGenomeTracks is a useful tool for visualizing open chromatin regions and gene regions concurrently. The following commands will add these features to your conda environment.

  - conda install -c bioconda -c conda-forge pygenometracks python=3.7
  - conda install -c bioconda deeptools

# Workflow Overview

![ImageOfWorkflow](https://github.com/CalebPecka/ATAC-Seq-Pipeline/blob/master/WorkflowDiagram.png)

This repository assumes you have already processed your raw fastq files, and your data is formatted as a series of .bam files, one per sample. The entire process can be run with the *initialization.sh* and *macro.sh* scripts. If you encounter difficulties, documentation for each command can be seen in the section below.

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

# How to Use / In depth Command Guide

*For each command, $1 refers to the .bam file you're processing, and $2 refers to your output directory.*

**macs2 callpeak -t $1 -g hs -f BAM -p 0.05 --seed 0 --bdg --outdir $2**

MACS2 is our peakcalling method. Open chromatin regions can be identified by regions where reads have piled up more than the background coverage.

  --'t' is the parameter used for our input file.
  
  --'g hs' identifies that we're using a human sequence.
  
  --'f' identifies our input file type, in this case bam.
  
  --'p' is our pvalue for identifying if the piled reads are significantly more accessible than the background coverage. These results are tabulated in a narrowPeaks file for future usage.
  
  --'seed' is a randomly generated started value which can be used to create reproducible results.
  
  --'bdg' indicates that we'd like to also create a bedgraph file, which we later convert to a bigwig file, a common file type for visualizing open chromatin regions.
  
  --'outdir' is our ouput directory.


**bamCoverage -b $1 -o $2/bigWig_coverage.bw**

BamCoverage is a method of converting a bedgraph file to a bigwig file. It is a component of the pygenometracks installation. If you don't wish to visualize your pygenometracks, this step can be ignored.

  --'b' is the parameter used for our input file.
  
  --'o' is our output directory.
  
  
**make_tracks_file --trackFiles $2/bigWig_coverage.bw $2/NA_peaks.narrowPeak GFFgenes.bed -o $2/tracks.ini**

This layers a series of tracks together for the pygenometracks visualization. The resulting .ini file is easily editable, see here: https://pygenometracks.readthedocs.io/en/latest/content/examples.html. Only one parameter is used, 'trackFiles'. Each file you include layer on top of each other. In this case, our visualization includes the bigWig coverage (peak visualization), narrowPeaks (bp location of peak visualized as a box plot), and a bed file for the start/stop regions of human genes.

**Rscript NarrowPeakSummitTracker.R $2/NA_peaks.narrowPeak GFFgenes.bed $2**

This R script will take the location of each peak summit (highest point of read coverage overlap) which was found 1000 bps upstream of human genes. We can often identify these peaks as promoter regions for their respective downstream genes. The resulting 90 bp sequence around each summit location is tabulated in a fasta file. Each fasta header including the start/stop position of the sequence, as well as the gene they theoretically promote.
