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


# Note
## checkm2 error  
checkm2 database --setdblocation /data/zhiyu/Database/atlas/atlas2.15.0/databases/CheckM2/CheckM2_database/uniref100.KO.1.dmnd   
checkm2 testrun --database_path /data/zhiyu/Database/atlas/atlas2.15.0/databases/CheckM2/CheckM2_database/uniref100.KO.1.dmnd 

## wrapper-prefix
如果一些github不能下载，采用如下方式：  
atlas run all --keep-going --wrapper-prefix "/data/zhiyu/data/software/atlas-2.15.2/snakemake-wrappers/" --omit-from run_decontamination   
atlas run all --wrapper-prefix "/data/zhiyu/data/software/atlas-2.15.2/snakemake-wrappers/"  

## dag
流程图,sample.tsv最好放一个样本，否则图很大  
atlas run all --dryrun --dag | dot -Tpdf > dag.pdf  
