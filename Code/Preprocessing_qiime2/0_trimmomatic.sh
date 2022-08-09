# run trimmomatic to remove low-quality ends of reads
cd raw_sequences

# Run trimmomatic with a forward loop, on all root sample files in a folder
for infile in *R1.fastq.gz
do
base=$(basename ${infile} R1.fastq.gz)
java -jar /home/qiime2/kathe/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE ${infile} ${base}R2.fastq.gz \
             ${base}R1.trim.fastq.gz ${base}R1.untrim.fastq.gz \
             ${base}R2.trim.fastq.gz ${base}R2.untrim.fastq.gz \
             SLIDINGWINDOW:5:30 MINLEN:25
done




# takes ~ 10min on pedro's PC
