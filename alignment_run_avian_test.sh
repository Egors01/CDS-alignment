
file_gff3=/home/egors/Projects/CDS-alignment/avian_test/sequence\(4\).gff3
input_fasta_raw=/home/egors/Projects/CDS-alignment/avian_test/sequence\(7\).fasta
filebed_out=/home/egors/Projects/CDS-alignment/avian_test/annot.bed
annot_file_name=/home/egors/Projects/CDS-alignment/avian_test/annot_cds.tsv
cleaned_fasta=avian_test/avian_input.fasta
out_result=avian_test/final_avian_alignment.fa
# clean fasta headers
python ./gb_fastaheader2id.py -i ./index_names -o $cleaned_fasta  -f $input_fasta_raw

#create CDS annotation table
gff2bed --do-not-sort < $file_gff3 > $filebed_out


python ./create_cds_annotation_table_from_gff.py -b $filebed_out -o $annot_file_name
python ./cds_alignment.py -f $cleaned_fasta -a $annot_file_name --out $out_result