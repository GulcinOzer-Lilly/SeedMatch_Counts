These scripts run at Brainiac

LISA transcriptome files for each genome:
/lrlhps/users/biorels_infra/LISA_PRD/PRD_FILES/LISA2.0/*/*_transcriptome.fasta


Jeremy prepared below tables to subset UTR and coding regions: 
/lrlhps/users/biorels_infra/LISA_PRD/PRD_FILES/LISA2.0/TR_LIST.csv
/lrlhps/users/biorels_infra/LISA_PRD/PRD_FILES/LISA2.0/TR_POS_TYPE.csv

Input: genomes.txt

#Generate all possible seeds and seeds with wobbles, save under data_prelim
./01_generate_seeds.sh

#Read each transcriptome and region definitions, save full transcripts and 3'UTR as separate Rdata under data_prelim
./02_run_prep_transcriptome_data.sh

#Generate seed match counts 
./03_run_SeedMatch_counts.sh
