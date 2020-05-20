
file_gff3=/home/egors/Projects/CDS-alignment/gb_swine_NA/swine_annot.gff3
filebed_out=/home/egors/Projects/CDS-alignment/gb_swine_NA/swine_annot.bed
annot_file_name=/home/egors/Projects/CDS-alignment/gb_swine_NA/swine_cds_annot_table.tsv
cleaned_fasta=./gb_swine_NA/swine_demo.fasta
# clean fasta headers
python ./gb_fastaheader2id.py -i ./gb_swine_NA/names_index.tsv -o $cleaned_fasta  -f ./gb_swine_NA/swine_seq.fasta

#create CDS annotation table
gff2bed --do-not-sort < $file_gff3 > $filebed_out


python ./create_cds_annotation_table_from_gff.py -b $filebed_out -o $annot_file_name

python ./cds_alignment.py -f $cleaned_fasta -a $annot_file_name --out ./gb_swine_NA/ALN.fasta --align