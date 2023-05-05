# atlas
## install
mamba env create -n atlas.2.15.2 --file atlasenv.yml  
conda activate atlas.2.15.2   
pip install --editable .  

atlas init --db-dir databases path/to/fastq/files    
atlas run all --keep-going  --report     

####################################################################


conda config --add channels defaults  
conda config --add channels bioconda  
conda config --add channels conda-forge  
mamba create -y -n atlasenv metagenome-atlas=2.15.0  



atlas init --db-dir databases path/to/fastq/files  
atlas run all --keep-going  --report   
