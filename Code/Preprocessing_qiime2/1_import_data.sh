#imports untrimmed sequences
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path manifest_untrimmed.txt \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path predada2_untrimmed_withprimers.qza

 #create visualization of demultiplxed files
 qiime demux summarize \
  --i-data predada2_untrimmed_withprimers.qza \
  --o-visualization predada2_untrimmed_withprimers.qzv


#imports trimmed sequences
#qiime tools import \
#--type 'SampleData[PairedEndSequencesWithQuality]' \
#--input-path manifest_trimmed.txt \
#--input-format PairedEndFastqManifestPhred33V2 \
#--output-path predada2_trimmed_withprimers.qza

 #create visualization of demultiplxed files
# qiime demux summarize \
#  --i-data predada2_trimmed_withprimers.qza \
#  --o-visualization predada2_trimmed_withprimers.qzv
