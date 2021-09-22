# Create a configuration environment for MACS2
conda activate base
conda create -n macsConfig
conda activate macsConfig
conda install -c bioconda macs2
conda env export > ../envs/macsConfig.yaml

# Create a configuration environment for R
conda activate base
conda create -n Rconfig
conda activate Rconfig
conda install -c r r
conda install -c bioconda bioconductor-bcrank
conda install -c bioconda bioconductor-biostrings
conda env export > ../envs/Rconfig.yaml

# Create a configuration environment for pyGenomeTracks
conda activate base
conda create -n pygenometracksConfig -c bioconda -c conda-forge pygenometracks python=3.7
conda activate pygenometracksConfig
conda env export > ../envs/pyGenomeTracksConfig.yaml

# Create a configuration environment for memeSuite
conda activate base
conda create -n memeSuiteConfig
conda activate memeSuiteConfig
conda install -c bioconda meme
conda env export > ../envs/memeSuiteConfig.yaml


# Return to the snakemake environment
conda activate snakemake