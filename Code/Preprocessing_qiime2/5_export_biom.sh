#export artefact into a feature table in the BIOM format, creating a new folder in the process
qiime tools export \
--input-path feature-table_trunc.qza \
--output-path 16S_phyloseq_input_FeatureTable

#Export the taxonomy from the .qza taxonomy assingment trained with the correct priemr set, it also crates a new folder in the process 
qiime tools export \
--input-path taxonomy_trunc.qza \
--output-path 16S_phyloseq_input_taxonomy

# export the representative sequences into a FASTA format...it should be able to be loaded into a BIOM file, but Pedro has not figured out how to do that yet. For now load it with phyloseq::import_biom()
qiime tools export \
  --input-path rep-seqs_trunc.qza \
  --output-path exported-rep-seqs_16S

#renames file of representative sequences
mv exported-rep-seqs_16S/dna-sequences.fasta exported-rep-seqs_16S/16S_dna-sequences.fasta

# adjust headers of taxonomy and mapping files
sed -i 's/Confidence/confidence/' 16S_phyloseq_input_taxonomy/taxonomy.tsv
sed -i 's/Taxon/taxonomy/' 16S_phyloseq_input_taxonomy/taxonomy.tsv
sed -i 's/Feature ID/#ASVID/' 16S_phyloseq_input_taxonomy/taxonomy.tsv
sed -i 's/sample-id/#sample-id/' Mapping_file.txt

# add metadata to the biom file
biom add-metadata \
-i 16S_phyloseq_input_FeatureTable/feature-table.biom \
-o 16S_phyloseq_input_FeatureTable_metadata_taxonomy.biom \
-m Mapping_file.txt \
--observation-metadata-fp 16S_phyloseq_input_taxonomy/taxonomy.tsv

# restore the original header of the mapping file to prevent conflicts with other mapping file ascessions
sed -i 's/#sample-id/sample-id/' Mapping_file.txt


