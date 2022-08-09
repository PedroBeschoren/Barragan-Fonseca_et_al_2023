#use dada2 to trim low-quality ends, remove chimeric sequences, remove singletons, join denoised paired-end reads, and then dereplicate into ASVs. trim size based on fastqc plot. migrate to a slliding window aproach when possible
qiime dada2 denoise-paired \
--i-demultiplexed-seqs preDada2_TrimmedPrimers.qza \
--p-trunc-len-f 180 \
--p-trunc-len-r 90 \
--p-max-ee-f 1 \
--p-max-ee-r 1 \
--p-n-threads 0 \
--o-table feature-table_trunc.qza \
--o-representative-sequences rep-seqs_trunc.qza \
--o-denoising-stats denoising-stats_trunc.qza
   
# visualize artefact, feature table  
qiime feature-table summarize \
--i-table feature-table_trunc.qza \
--o-visualization feature-table_trunc.qzv \
--m-sample-metadata-file Mapping_file.txt
  
# visualize artefact, representative sequences table  
qiime feature-table tabulate-seqs \
--i-data rep-seqs_trunc.qza \
--o-visualization rep-seqs_trunc.qzv
  
# visualize artefact, denoising stats 
qiime metadata tabulate \
--m-input-file denoising-stats_trunc.qza \
--o-visualization denoising-stats_trunc.qzv
