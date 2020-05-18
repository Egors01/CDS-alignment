from Bio import SeqIO
import pandas as pd
class CdsAlignment():
    def __init__(self,fasta_file,annotation_table):

        self.fasta_list = list(SeqIO.parse(fasta_file,format='fasta'))
        self.annotation = pd.read_csv(annotation_table,sep='\t')