# atlas
## install
conda config --add channels defaults  
conda config --add channels bioconda  
conda config --add channels conda-forge  
mamba create -y -n atlasenv metagenome-atlas=2.15.0  



atlas init --db-dir databases path/to/fastq/files  
atlas run all --keep-going  --report   
