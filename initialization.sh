#This shell script is used to intialize the project. The R package install install system fails if a single package dependency is not upgraded, as it is not capable of automatically selecting the latest compatible version. Each 'problematic' package dependency is manually installed from the source package below:
conda env create ATAC
conda activate ATAC
conda install -c conda-forge r=4.0.0
r -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/BiocGenerics_0.34.0.tar.gz", repos = NULL, type = "source")' 
r -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/S4Vectors_0.26.1.tar.gz", repos = NULL, type = "source")'
r -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/IRanges_2.22.2.tar.gz", repos = NULL, type = "source")'
r -e 'install.packages("https://www.bioconductor.org/packages/release/bioc/src/contrib/zlibbioc_1.34.0.tar.gz", repos = NULL, type = "source")'
r -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/XVector_0.28.0.tar.gz", repos = NULL, type = "source")'
r -e 'install.packages("https://cran.r-project.org/src/contrib/crayon_1.3.4.tar.gz", repos = NULL, type = "source")'
r -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/Biostrings_2.56.0.tar.gz", repos = NULL, type = "source")'
r -e 'install.packages("https://bioconductor.org/packages/release/bioc/src/contrib/BCRANK_1.50.0.tar.gz", repos = NULL, type = "source")'
r -e 'install.packages("https://cran.r-project.org/src/contrib/pixmap_0.4-11.tar.gz", repos = NULL, type = "source")'
r -e 'install.packages("https://cran.r-project.org/src/contrib/sp_1.4-2.tar.gz", repos = NULL, type = "source")'
r -e 'install.packages("https://cran.r-project.org/src/contrib/ade4_1.7-15.tar.gz", repos = NULL, type = "source")'
r -e 'install.packages("https://cran.r-project.org/src/contrib/segmented_1.2-0.tar.gz", repos = NULL, type = "source")'
r -e 'install.packages("https://cran.r-project.org/src/contrib/seqinr_3.6-1.tar.gz", repos = NULL, type = "source")'

#We also want to create a separate directory for the gene bank we download from. The link below can be modified for any reference genome you'd rather use for your project.
mkdir HG38
cd HG38
curl -LJO "ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"
#If you're using a different reference genome than the link provided above, make sure to replace the GFFgenes.bed file in the initial repository. The GFFgenes.bed file is a list of the start/stop positions for every gene in the provided genome sequence. If they are not provided from the same source, your gene sequences will have an incorrect output!
