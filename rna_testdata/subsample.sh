
curl -L https://osf.io/p4fy5/download -o nema_subset_0Hour.zip
curl -L https://osf.io/ewyv5/download -o nema_subset_6Hour.zip
unzip nema_subset_0Hour.zip
unzip nema_subset_6Hour.zip
mkdir -p orig
mv *gz ./orig

lines=5000

gunzip -c orig/0Hour_ATCACG_L002_R1_001.fastq.gz | head -$lines | gzip -9 > 0Hour_ATCACG_L002_R1_001.fastq.gz
gunzip -c orig/0Hour_ATCACG_L002_R1_002.fastq.gz | head -$lines | gzip -9 > 0Hour_ATCACG_L002_R1_002.fastq.gz
gunzip -c orig/0Hour_ATCACG_L002_R1_003.fastq.gz | head -$lines | gzip -9 > 0Hour_ATCACG_L002_R1_003.fastq.gz
gunzip -c orig/0Hour_ATCACG_L002_R1_004.fastq.gz | head -$lines | gzip -9 > 0Hour_ATCACG_L002_R1_004.fastq.gz
gunzip -c orig/0Hour_ATCACG_L002_R1_005.fastq.gz | head -$lines | gzip -9 > 0Hour_ATCACG_L002_R1_005.fastq.gz
gunzip -c orig/0Hour_ATCACG_L002_R2_001.fastq.gz | head -$lines | gzip -9 > 0Hour_ATCACG_L002_R2_001.fastq.gz
gunzip -c orig/0Hour_ATCACG_L002_R2_002.fastq.gz | head -$lines | gzip -9 > 0Hour_ATCACG_L002_R2_002.fastq.gz
gunzip -c orig/0Hour_ATCACG_L002_R2_003.fastq.gz | head -$lines | gzip -9 > 0Hour_ATCACG_L002_R2_003.fastq.gz
gunzip -c orig/0Hour_ATCACG_L002_R2_004.fastq.gz | head -$lines | gzip -9 > 0Hour_ATCACG_L002_R2_004.fastq.gz
gunzip -c orig/0Hour_ATCACG_L002_R2_005.fastq.gz | head -$lines | gzip -9 > 0Hour_ATCACG_L002_R2_005.fastq.gz 
gunzip -c orig/6Hour_CGATGT_L002_R1_001.fastq.gz | head -$lines | gzip -9 > 6Hour_CGATGT_L002_R1_001.fastq.gz
gunzip -c orig/6Hour_CGATGT_L002_R1_002.fastq.gz | head -$lines | gzip -9 > 6Hour_CGATGT_L002_R1_002.fastq.gz
gunzip -c orig/6Hour_CGATGT_L002_R1_003.fastq.gz | head -$lines | gzip -9 > 6Hour_CGATGT_L002_R1_003.fastq.gz
gunzip -c orig/6Hour_CGATGT_L002_R1_004.fastq.gz | head -$lines | gzip -9 > 6Hour_CGATGT_L002_R1_004.fastq.gz
gunzip -c orig/6Hour_CGATGT_L002_R1_005.fastq.gz | head -$lines | gzip -9 > 6Hour_CGATGT_L002_R1_005.fastq.gz
gunzip -c orig/6Hour_CGATGT_L002_R2_001.fastq.gz | head -$lines | gzip -9 > 6Hour_CGATGT_L002_R2_001.fastq.gz
gunzip -c orig/6Hour_CGATGT_L002_R2_002.fastq.gz | head -$lines | gzip -9 > 6Hour_CGATGT_L002_R2_002.fastq.gz
gunzip -c orig/6Hour_CGATGT_L002_R2_003.fastq.gz | head -$lines | gzip -9 > 6Hour_CGATGT_L002_R2_003.fastq.gz
gunzip -c orig/6Hour_CGATGT_L002_R2_004.fastq.gz | head -$lines | gzip -9 > 6Hour_CGATGT_L002_R2_004.fastq.gz
gunzip -c orig/6Hour_CGATGT_L002_R2_005.fastq.gz | head -$lines | gzip -9 > 6Hour_CGATGT_L002_R2_005.fastq.gz
