import argparse
def rename_headers(input_file ,output_file,index_file):
    from Bio import SeqIO
    import pandas as pd

    fasta_list = list(SeqIO.parse(input_file,format='fasta'))
    new_records = []
    new_id_to_old = {}
    for record in fasta_list:
        new_id_to_old[record.id] = record.description
        record.description = ''
        record.name = ''
        new_records.append(record)
    df = pd.DataFrame.from_dict(new_id_to_old,orient='index')
    df.to_csv(index_file,header=None,sep='\t')
    with open(output_file,'w') as handle:
        SeqIO.write(new_records,handle,format='fasta')
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', type=str, help='input  fasta', required=True, action='store')
    parser.add_argument('-ind', '--index', type=str, help='output table with <new_header><old header> records', required=True, action='store')
    parser.add_argument('-o', '--out', type=str, help='output fasta file with new headers', required=True, action='store')
    args = parser.parse_args()

    input_file = args.fasta
    output_file = args.out
    index_file = args.index
    rename_headers(input_file ,output_file,index_file)
