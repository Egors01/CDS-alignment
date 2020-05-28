input_fasta=./sabine_test/test_input.fasta
annot_file_name=./sabine_test/annot_cds.tsv
# output
out_result=./sabine_test/final_alignment_sabine_test.fa


# stop and show input here
echo ./cds_alignment.py -f $input_fasta -a $annot_file_name --out $out_result --align
python ./cds_alignment.py -f $input_fasta -a $annot_file_name --out $out_result --align

# restore after external align
#python ./cds_alignment.py -f $input_fasta -a $annot_file_name --out $out_result --no-align --only-restore