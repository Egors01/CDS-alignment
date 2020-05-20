import argparse
import pandas as pd
from copy import deepcopy
#input_bed_file = 'genbank_sample_data/complex_data_segment7/sample_annotations.bed'
#output_file = 'genbank_sample_data/complex_data_segment7/annotation_table.tsv'



def create_cds_annotation(input_bed_file,output_file):
    ref_data = pd.read_csv(input_bed_file, sep='\t', header=None)
    ref_data = ref_data[[0, 1, 2, 3, 7, 9]]
    ref_data.columns = ['alignment_id', 'start', 'end', 'feature_id', 'feature_type', 'expand_annot']
    # bypassing the hole for influenza
    try:
        ref_data['segment'] = ref_data['expand_annot'].str.extract('segment=(\d)', expand=True)
    except:
        ref_data['segment'] = '_'
    cds_df = pd.DataFrame()
    for accid in ref_data['alignment_id'].unique():
        segment_data = deepcopy(ref_data.loc[ref_data['alignment_id'] == accid])
        # bypassing the hole for influenza segments
        try:
            segment_number = str(segment_data.dropna().segment.unique()[0])
            segment_data['segment'] = segment_data['segment'].fillna(segment_number)
        except:
            segment_data['segment'] = '_'
        segment_data = segment_data.loc[segment_data['feature_type'] == 'CDS']
        cds_df = cds_df.append(segment_data)
    result_cds_df = cds_df.sort_values(by=['alignment_id','start','feature_id'])[
        ['alignment_id', 'start', 'end', 'feature_id', 'feature_type','segment', 'expand_annot']]
    print('INFO:\tCreated a CDS annotation table {}'.format(output_file))
    result_cds_df.to_csv(output_file, sep='\t', index=False,header=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bed', type=str, help='input  bed', required=True, action='store')
    parser.add_argument('-o', '--out', type=str, help='output annotation file', required=True, action='store')
    args = parser.parse_args()

    input_file = args.bed
    output_file = args.out
    # input_file = '/home/egors/Projects/CDS-alignment/genbank_sample_data/sequences.fasta'
    # output_file = '/home/egors/Projects/CDS-alignment/genbank_sample_data/to_align.fasta'
    # index_file = '/home/egors/Projects/CDS-alignment/genbank_sample_data/fasta_header_to_id.tsv'
    create_cds_annotation(input_bed_file=input_file ,output_file=output_file)
