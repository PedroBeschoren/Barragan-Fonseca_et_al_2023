#assing ttaxonomy ranks on your representative sequences, based on the trained SILVA reference
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-ssu-nr99-515f-806r-classifier.qza \
  --i-reads rep-seqs_trunc.qza \
  --o-classification taxonomy_trunc.qza

# create a vizualization of the taxonomy
qiime metadata tabulate \
  --m-input-file taxonomy_trunc.qza \
  --o-visualization taxonomy_trunc.qzv

