# ATAC-Seq-Pipeline


*Requirements*

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

*Conda Environment*

  - conda create -n pygenometracks -c bioconda -c conda-forge pygenometracks python=3.7
  - conda install -c bioconda deeptools
  
*How to Use*

  - Run the macro.sh script with 3 arguments. The first being an input bam file, the second being an output directory.
  - For example:
    sh macro.sh /data/fakeFileLocation/3c39.bam /MACS2_out
  
  
  - Then chill for like an hour.
