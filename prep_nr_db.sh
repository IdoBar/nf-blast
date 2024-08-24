echo "source ~/.bashrc
set -Eeo pipefail 
DB_PATH=/project/aerc/ncbi/v5
conda activate aDNA
blastdbcmd -entry all -db \${DB_PATH}/nr -out \${DB_PATH}/nr.fasta 
printf \"Finished extracting sequences from nr database as 'nr.fasta'\n\"
diamond makedb -p 6 --in \${DB_PATH}/nr.fasta -d \${DB_PATH}/nr_diamond --taxonmap \${DB_PATH}/prot.accession2taxid.FULL.gz --taxonnodes \${DB_PATH}/nodes.dmp --taxonnames \${DB_PATH}/names.dmp
printf \"Finished preparing nr Diamond database as 'nr_diamond'\n\"
rm \${DB_PATH}/nr.fasta
printf \"Removed 'nr.fasta', script completed...\n\"" | qsub -l select=1:ncpus=8:mem=32GB,walltime=15:00:00 -N prep_nr_db


echo "source ~/.bashrc
set -Eeo pipefail 
DB_PATH=/project/aerc/ncbi/v5
conda activate aDNA
gzip -v -t \${DB_PATH}/prot.accession2taxid.FULL.gz
printf \"\$(date) - Finished validating 'prot.accession2taxid.FULL.gz'.\n\"
diamond makedb -p 6 --in \${DB_PATH}/nr.fasta -d \${DB_PATH}/nr --taxonmap \${DB_PATH}/prot.accession2taxid.FULL.gz --taxonnodes \${DB_PATH}/nodes.dmp --taxonnames \${DB_PATH}/names.dmp
printf \"$(date) - Finished preparing Diamond database 'nr.dmnd'. \n\"
rm \${DB_PATH}/nr.fasta
printf \"$(date) - Removed 'nr.fasta', script completed. \n\"" | qsub -l select=1:ncpus=8:mem=32GB,walltime=10:00:00 -N prep_nr_db