# this manifest file indicates where forward are reverse reads are in your filesystem  for each one of your samples
# you will need this to import your data into qiime2. your folder with raw sequences should only have the raw sequences files
# these are the commands performed bellow:

# set your header fieldsnano
# fills "sample-id" with the 5th of separation of your filemanes, with element separated by  ".", and then removes "_R" # note: (filenames are collected by ls, and then filtere$
# fills "forward-absolute-filepath" with forward reads of 16S files
# fills "reverse-absolute-filepath" with reveresr reads of 16S files
# saves output, appending new files names on the document

printf "%s\t%s\t%s\n%s\n" "sample-id" "forward-absolute-filepath" "reverse-absolute-filepath" \
"$(paste <(awk -F '_R' '{print $1}' <(awk -F '.' '{print $5}' <(fgrep "R1" <(ls raw_sequences)))) \
<(paste <(sed -e 's/^/$PWD\/raw_sequences\//' <(fgrep "R1" <(ls raw_sequences))) \
<(sed -e 's/^/$PWD\/raw_sequences\//' <(fgrep "R2" <(ls raw_sequences)))))" \
>> manifest_untrimmed.txt
