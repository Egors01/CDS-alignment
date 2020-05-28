
# input from raw gb batch download
file_gff3=./swine_demo/swine_annot.gff3
genbank_fasta=./swine_demo/swine_seq.fasta
# output
out_result=swine_demo/final_alignment.fa

# intermediate files for cleaning the raw
filebed_out=./swine_demo/swine_annot_bed
annot_file_name=./swine_demo/annot_cds.tsv
cleaned_fasta=./swine_demo/swine_input.fasta
saved_fastaheaders=./swine_demo/swine_names_index.tsv


# clean fasta headers
python ./gb_fastaheader2id.py -i $saved_fastaheaders -o $cleaned_fasta  -f $genbank_fasta
#create CDS annotation table
gff2bed --do-not-sort < $file_gff3 > $filebed_out
python ./create_cds_annotation_table_from_gff.py -b $filebed_out -o $annot_file_name

# stop and show input here
echo ./cds_alignment.py -f $cleaned_fasta -a $annot_file_name --out $out_result --align
python ./cds_alignment.py -f $cleaned_fasta -a $annot_file_name --out $out_result --align

# restore after external align
#python ./cds_alignment.py -f $cleaned_fasta -a $annot_file_name --out $out_result --no-align --only-restore