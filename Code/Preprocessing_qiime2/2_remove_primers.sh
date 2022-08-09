# run cutadapt plugin of qiime, removing primer sequences and discarding reads without primers
qiime cutadapt trim-paired \
--p-cores 2 \
--i-demultiplexed-sequences predada2_untrimmed_withprimers.qza \
--p-front-f GTGCCAGCMGCCGCGGTAA \
--p-front-r GGACTACHVGGGTWTCTAAT \
--p-error-rate 0 \
--p-no-indels \
--p-discard-untrimmed \
--o-trimmed-sequences preDada2_Untrimmed_noPrimers.qza \
--verbose

#create visualization of artefact
qiime demux summarize \
--i-data preDada2_Untrimmed_noPrimers.qza \
--o-visualization preDada2_Untrimmed_noPrimers.qzv
