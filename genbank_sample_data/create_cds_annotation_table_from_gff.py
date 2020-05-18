import pandas as pd
from copy import deepcopy
input_bed_file = '/home/egors/Projects/CDS-alignment/genbank_sample_data/sample_annotations.bed'
output_file ='./annotation_table.tsv'
ref_data = pd.read_csv(input_bed_file, sep='\t', header=None)
ref_data = ref_data[[0, 1, 2, 3, 7, 9]]
ref_data.columns = ['gene_id', 'start', 'end', 'feature_id', 'feature_type', 'expand_annot']
# bypassing the hole for influenza
try:
    ref_data['segment'] = ref_data['expand_annot'].str.extract('segment=(\d)', expand=True)
except:
    ref_data['segment'] = '_'
cds_df = pd.DataFrame()
for accid in ref_data['gene_id'].unique():
    segment_data = deepcopy(ref_data.loc[ref_data['gene_id'] == accid])
    # bypassing the hole for influenza segments
    try:
        segment_number = str(segment_data.dropna().segment.unique()[0])
        segment_data['segment'] = segment_data['segment'].fillna(segment_number)
    except:
        segment_data['segment'] = '_'
    segment_data = segment_data.loc[segment_data['feature_type'] == 'CDS']
    cds_df = cds_df.append(segment_data)
result_cds_df = cds_df.sort_values(by=['gene_id','start','feature_id'])[
    ['gene_id', 'start', 'end', 'segment', 'feature_id', 'feature_type', 'expand_annot']]
print('INFO:\tCreated a CDS annotation table {}'.format(output_file))
result_cds_df.to_csv(output_file, sep='\t', index=False)