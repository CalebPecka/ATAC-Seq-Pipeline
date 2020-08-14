# ATAC-Seq-Pipeline


**Requirements**

  - Python >= 3.6
  - R >= 4.0
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

  - conda create -n pygenometracks -c bioconda -c conda-forge pygenometracks python=3.7
  - conda install -c bioconda deeptools
  
**How to Use**

  - In your main directory, run the *initialization.sh* script. This will download hg38 chromosomal data which we need for future steps.
  - For each sample, run the *macro.sh* script with 3 arguments. The first being an input bam file, the second being an output directory.
  - For example:
    *sh macro.sh /data/fakeFileLocation/3c39.bam /MACS2_out*
  
  
  - Then chill for like an hour.
  - Once complete, your *tracks.ini* has a lot of visual customization options. Here is some example code to give you an idea of what objects you might want to include: https://pygenometracks.readthedocs.io/en/latest/content/examples.html
  - When you want to plot a region of the genome, use: *pyGenomeTracks --tracks bigwig_track.ini --region chr15:2,500,000-2,800,000 -o bigwig.png*, with variables adjusted as you desire. I recommend using a base pair range of 300,000 for the best possible graph readability.
