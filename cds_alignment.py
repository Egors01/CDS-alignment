import argparse
import logging
import operator
import os
import subprocess
from collections import Counter
import warnings
import pandas as pd
from Bio import SeqIO

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


class CdsAlignment():
    def __init__(self, fasta_file, annotation_table, output_file, force_three_partitions=True, do_align=True,
                 only_restore=False):

        self.fasta_file = fasta_file
        self.init_folders()
        self.init_empty_variables()

        self.final_alignment_file = output_file

        # Read fasta
        self.fasta_list = list(SeqIO.parse(fasta_file, format='fasta'))

        # save its record order and full descriptions
        self.fasta_list_records_id_ordered = [record.id for record in self.fasta_list]
        self.fasta_list_full_descriptions_from_id = {record.id: record.description for record in self.fasta_list}
        self.n_sequences = len(self.fasta_list_records_id_ordered)

        # read annotation
        self.annotation = pd.read_csv(annotation_table, sep='\t', header=None, usecols=[0, 1, 2],
                                      names=['alignment_id', 'start', 'end'])
        # check if annotation match fasta
        self.avaliable_ids = dict(
            zip(self.annotation['alignment_id'], ['' for i in range(0, self.annotation.shape[0])]))

        for record in self.fasta_list:
            if record.id not in self.avaliable_ids:
                raise ValueError(
                    'There is no record in CDS annotation that match the >{} fasta header'.format(record.id))

        # Start partitioning

        self.force_three_partitions = force_three_partitions

        # EXECUTION
        # print('INFO:\tUsing 5-UTR CDS 3-UTR partitioning')
        self.run_partitioning()
        self.apply_split_to_three_partitions()
        self.save_files_to_align()
        if do_align and not only_restore:
            # go the full procedure if we do alignment
            logger.info('Running alignment...')
            # print('INFO:\tRunning alignment...')
            self.align()
            self.restore_translated()
            self.restore_aligned_partitions()
        elif only_restore:
            logger.info('Looking for pre-existing alignment <names>.fasta_al')
            self.restore_translated()
            self.restore_aligned_partitions()
        elif do_align and only_restore:
            raise LookupError('Confusing argument. Allowed: One of only-restore or do-align should be False')

    def init_folders(self):
        fasta_file_short_name = os.path.basename(self.fasta_file).replace('.fa', '').replace('.fasta', '')

        self.working_dir = os.path.dirname(self.fasta_file)
        self.tmp_dir = os.path.join(self.working_dir, 'align.{}.parts'.format(fasta_file_short_name))

        if not os.path.exists(self.tmp_dir):
            os.makedirs(self.tmp_dir)

        self.split_analysis_file = os.path.join(self.tmp_dir, 'applied_orf_splitting.tsv')
        self.alignnments_input_output_file = os.path.join(self.tmp_dir, 'files_list.tsv')
        self.empty_ends_file_five = os.path.join(self.tmp_dir, '5_UTR_emtpy_ends')
        self.empty_ends_file_three = os.path.join(self.tmp_dir, '3_UTR_emtpy_ends')
        self.restored_cds_fasta_file = os.path.join(self.tmp_dir, 'cds_nucleotde_restored.fasta_al')

        self.split_filenames = ['five_utr.fasta', 'cds_translated.fasta', 'three_utr.fasta', 'cds_nucleotide.fasta']

    def init_empty_variables(self):
        self.transcript_partitioning = {}
        self.five_primes_list = []
        self.three_primes_list = []
        self.cds_aa_list = []
        self.cds_nuc_list = []
        self.five_primes_empty = []
        self.three_primes_empty = []
        self.splittted_files_types = {}
        self.n_partitions = 0
        self.restored_cds_aligned = []
        self.alignnments_input_output_df = pd.DataFrame()

    def run_partitioning(self):

        splits_per_record = {}
        for record in self.fasta_list:
            self.transcript_partitioning[record.id] = []
            transcript_df = self.annotation.loc[self.annotation['alignment_id'] == record.id]

            partition_start = transcript_df.iloc[0]['start']
            partition_end = transcript_df.iloc[0]['end']

            for index, row in transcript_df.iterrows():

                cds_start = row['start']
                cds_end = row['end']
                # cds_name = row['feature_id']
                if self.force_three_partitions:
                    self.transcript_partitioning[record.id].append([partition_start, partition_end])
                    splits_per_record[record.id] = len(self.transcript_partitioning[record.id]) * 2 + 1
                    break
                else:
                    # FInding optimal partitioning for complex ORF structure
                    # if we already captured cds in partition
                    if (cds_start >= partition_start) and (cds_end <= partition_end):
                        pass
                    else:
                        # extend right
                        if cds_start >= partition_start and (cds_end >= partition_end) and cds_start < partition_end:
                            partition_end = cds_end
                        if cds_start < partition_start and (cds_end <= partition_end):
                            partition_start = cds_start

                        # if the start and end of new cds is
                        # downstream the end of current partition we make new partition for this cds
                        if (cds_start > partition_end) and (cds_end > partition_end):
                            self.transcript_partitioning[record.id].append([partition_start, partition_end])
                            partition_start = cds_start
                            partition_end = cds_end
                            # print('suggesting new partition for cds {} record id {}'.format(cds_name,record.id))
                self.transcript_partitioning[record.id].append([partition_start, partition_end])
                splits_per_record[record.id] = len(self.transcript_partitioning[record.id]) * 2 + 1
        # complete partitioning with five and three UTR
        for record in self.fasta_list:
            seq_len = len(record.seq)
            first_cds_start = self.transcript_partitioning[record.id][0][0]
            last_cds_end = self.transcript_partitioning[record.id][-1][-1]
            self.transcript_partitioning[record.id].insert(0, [0, first_cds_start])
            self.transcript_partitioning[record.id].append([last_cds_end, seq_len])

        suggested_splits = dict(Counter([x for x in splits_per_record.values()]))
        # Save record with the applied partitioning
        f = open(self.split_analysis_file, 'w')
        for key, value in splits_per_record.items():
            f.write('{}\t{}\t{}\n'.format(key, value, self.transcript_partitioning[key]))

        # Set n_partitions
        if self.force_three_partitions:
            self.n_partitions = 3

        else:
            # apply partitioning based on cds mapping
            self.n_partitions = max(suggested_splits.items(), key=operator.itemgetter(1))[0]
            logger.info(
                'Suggested splits based on povided annotation <n_splits>:<count>\n\tids and count are in{}'.format(
                    self.split_analysis_file))
            for key, value in suggested_splits.items():
                logger.info('{}:{}'.format(key, value + 2))

        logger.info('Partitioning complete')

    def apply_split_to_three_partitions(self):

        # check if partitioning is okay
        for record in self.fasta_list:
            if record.id not in self.transcript_partitioning:
                logger.error('There is no header >{} in partitioning record'.format(record))
                raise ValueError()
        if self.n_partitions == 3:
            pass
        else:
            logger.error('Three (3) partitions for the split selected. However there is {} in the split.'.format(
                self.n_partitions))
            raise ValueError()

        # spliting into three parts
        ids_to_drop_from_input = []
        for record in self.fasta_list:
            [start, end] = self.transcript_partitioning[record.id][1]
            five_seq = SeqRecord(seq=record.seq[:start], id=record.id, name='', description='')
            orf_aa_seq = SeqRecord(seq=record.seq[start:end].translate(), id=record.id, name='', description='')
            three_seq = SeqRecord(seq=record.seq[end:], id=record.id, name='', description='')
            orf_nuc_seq = SeqRecord(seq=record.seq[start:end], id=record.id, name='', description='')
            orf_mrna = str(orf_nuc_seq.seq).replace('T', 'U')

            # check if only one stop is present and if CDS ends with stop
            if orf_aa_seq.seq[-1] == '*' and str(orf_aa_seq.seq).count('*') == 1:
                pass

            else:
                logger.warning(
                    'CDS does not end with stop codon or multiple stop codon encountered.\n Sequence id {}\n Sequence will be dropped'.format(
                        record.id))
                ids_to_drop_from_input.append(record.id)

            # check start codon                                                                                      orf_aa_seq.seq),RuntimeWarning,stacklevel=2)
            if orf_aa_seq.seq[0] == 'M' and orf_aa_seq.seq[-1] == '*':
                pass
            elif orf_mrna[0:3] == 'GUG' or \
                    orf_mrna[0:3] == 'UUG' or \
                    orf_mrna[0:3] == 'AUU' or \
                    orf_mrna[0:3] == 'AUA':
                logger.warning(
                    'The start codon for sequence id {} is not AUG. Encountered start codon {} '.format(record.id,
                                                                                                        orf_mrna[0:3]))
            else:
                text = 'Incorrect protein translation for {}. Expected it to be M->* AA sequence or with alternative GUG | UUG | AUU | AUA start codons'.format(
                    record.id)
                logger.warning(text)
                ids_to_drop_from_input.append(record.id)

                # go to the next record
                continue

            self.cds_aa_list.append(orf_aa_seq)
            self.cds_nuc_list.append(orf_nuc_seq)

            if not len(five_seq.seq) == 0:
                self.five_primes_list.append(five_seq)
            else:
                logger.debug('Record {} does not have 5\'-UTR'.format(record.id))
                self.five_primes_empty.append(record.id)

            if not len(three_seq.seq) == 0:
                self.three_primes_list.append(three_seq)
            else:
                logger.debug('Record {} does not have 3\'-UTR'.format(record.id))
                self.three_primes_empty.append(record.id)

        # drop the ignored sequences from the input info lists (they are already not appended to alignments)

        self.fasta_list = [record for record in self.fasta_list if
                           record.id not in set(ids_to_drop_from_input)]

        self.fasta_list_records_id_ordered = [record_id for record_id in self.fasta_list_records_id_ordered if
                                              record_id not in set(ids_to_drop_from_input)]
        self.n_sequences = len(self.fasta_list_records_id_ordered)
        logger.info('Ingnored sequences {}'.format(set(ids_to_drop_from_input)))

        return

    def save_files_to_align(self):

        # saving alignments types and assigning filenames
        if self.n_partitions == 3:
            alignment_types = ['nuc', 'aa', 'nuc', 'nuc_cds']
            files = [os.path.join(self.tmp_dir, file) for file in self.split_filenames]
            lists = [self.five_primes_list, self.cds_aa_list, self.three_primes_list, self.cds_nuc_list]
        else:
            # catch unwritten processing of multiple cds split
            raise RuntimeError('Should not be here')
            pass

        n_sequences_in_partition = {}

        # saving
        for file, fasta_list, aln_type in zip(files, lists, alignment_types):
            self.splittted_files_types[file] = aln_type
            n_sequences_in_partition[file] = len(fasta_list)
            with open(file, 'w') as handle:
                SeqIO.write(fasta_list, handle, 'fasta')

        # save empty ends
        f = open(self.empty_ends_file_five, 'w')
        for x in self.five_primes_empty:
            f.write('{}\n'.format(x))
        f.close()
        f = open(self.empty_ends_file_three, 'w')
        for x in self.three_primes_empty:
            f.write('{}\n'.format(x))
        f.close()

        # save input-output files names for alignment outside
        df_in = pd.DataFrame.from_dict(self.splittted_files_types, orient='index').reset_index()
        df_in.columns = ['input_files', 'alignment_type']
        df_in['output_files'] = df_in['input_files'].str.replace('.fasta', '').str.replace('.fa', '') + '.fasta_al'
        df_in['n_sequences'] = pd.Series(
            [n_sequences_in_partition[inpf] for inpf in df_in['input_files'].values.tolist()])

        self.alignnments_input_output_df = df_in[['input_files', 'output_files', 'alignment_type', 'n_sequences']]
        self.alignnments_input_output_df.to_csv(self.alignnments_input_output_file, sep='\t', index=False)

        # create emtpy output files (to simplify processing content-empty partitions)
        for empty_output_file in df_in['output_files'].values.tolist():
            f = open(empty_output_file, 'w')
            f.close()

    def align(self):

        for i, row in self.alignnments_input_output_df.iterrows():
            input = row['input_files']
            output = row['output_files']
            aln_type = row['alignment_type']
            cline = ''
            if logger.level == logging.DEBUG:
                v = '-v'
            else:
                v = ''
            require_alignment = bool(row['n_sequences'] > 0)
            if aln_type == 'aa' and require_alignment:
                cline = 'clustalo -i {input} -o {out}  --threads 2 {v} --force --seqtype protein '.format(input=input,
                                                                                                          out=output,
                                                                                                          v=v)
            elif (aln_type == 'nuc' or aln_type == 'nuc_cds') and require_alignment:
                cline = 'clustalo -i {input} -o {out}  --threads 2 {v} --force --seqtype DNA '.format(input=input,
                                                                                                      out=output, v=v)
            logger.info(cline)
            return_code = subprocess.run(cline, shell=True).returncode
            if return_code != 0:
                raise RuntimeError('Alignment program had non-zero return code for run\n\t {}'.format(cline))

        return

    def restore_translated(self):
        fasta_nuc = self.cds_nuc_list

        algfilename = self.alignnments_input_output_df.loc[self.alignnments_input_output_df['alignment_type'] == 'aa'][
            'output_files'].to_list()[0]
        fasta_aa = list(SeqIO.parse(algfilename, format='fasta'))

        for record_aa, record_nuc in zip(fasta_aa, fasta_nuc):
            restored_sequence = ''
            start = 0
            stop = 3
            if str(record_nuc.id) != str(record_aa.id):
                raise LookupError('Record order mismatch CDS NUC->{} CDS AA -> {}'.format(record_nuc.id, record_aa.id))
            for aa in record_aa.seq:
                # take codon from original fasta
                if aa != '-':
                    codon = record_nuc.seq[start:stop]
                    restored_sequence += codon
                    start += 3
                    stop += 3
                    if codon.translate() == aa:
                        pass
                    else:
                        raise ValueError('wrong translation {}'.format(record_nuc.id))
                # append three gaps
                elif aa == '-':
                    restored_sequence += '---'
                else:
                    LookupError('Should be either AA or gap in AA cds alignment. But found {}'.format(aa))
            self.restored_cds_aligned.append(
                SeqRecord(seq=restored_sequence, id=record_nuc.id, name='', description=''))
        with open(self.restored_cds_fasta_file, 'w') as handle:
            SeqIO.write(self.restored_cds_aligned, handle, 'fasta')

        logger.info('Restored nucleotide seqeunce from protein alignment for CDS ')

        return

    def restore_aligned_partitions(self):
        algfilenames = \
            self.alignnments_input_output_df.loc[self.alignnments_input_output_df['alignment_type'] == 'nuc'][
                'output_files'].to_list()
        five_utr_file = algfilenames[0]
        three_utr_file = algfilenames[1]

        fasta_fives = dict(SeqIO.index(five_utr_file, format='fasta'))
        fasta_threes = dict(SeqIO.index(three_utr_file, format='fasta'))
        fasta_cds_aln = {record.id: record for record in self.restored_cds_aligned}

        # create mock deletion sequnces taht will stand for missing 5 and 3 UTR
        if not fasta_threes.__len__() == 0:
            mock_three_seq = '-' * len(fasta_threes.get(list(fasta_threes)[0]).seq)
        else:
            mock_three_seq = ''
        if not fasta_fives.__len__() == 0:
            mock_five_seq = '-' * len(fasta_fives[list(fasta_fives.keys())[0]].seq)
        else:
            mock_five_seq = ''

        concatenated_fasta_list = []

        # adding mock sequences to 5 and 3 ends
        for missing_id in self.fasta_list_records_id_ordered:
            if missing_id not in fasta_threes.keys():  # extra check
                fasta_threes[missing_id] = SeqRecord(seq=mock_three_seq, id=missing_id)
            if missing_id not in fasta_fives.keys():
                fasta_fives[missing_id] = SeqRecord(seq=mock_five_seq, id=missing_id)

        # creating new list with concatenated aligned parts
        for fasta_id in self.fasta_list_records_id_ordered:
            new_seq = fasta_fives[fasta_id].seq + fasta_cds_aln[fasta_id].seq + fasta_threes[fasta_id].seq
            description = self.fasta_list_full_descriptions_from_id[fasta_id]
            concatenated_fasta_list.append(SeqRecord(seq=new_seq, id='', description=description, name=''))

        # check if the original equals ungapped final
        for i, fasta_id in enumerate(self.fasta_list_records_id_ordered):
            if self.fasta_list[i].seq != concatenated_fasta_list[i].seq.ungap('-'):
                raise LookupError(
                    'Restoring process failed. record {} with index {}: input and ungapped output doesnot match'.format(
                        fasta_id, i))

        # save concatenated result
        # write formatted fasta output
        handle = open(self.final_alignment_file, 'w')
        for fasta_rec in concatenated_fasta_list:
            lines = []
            for i in range(0, len(fasta_rec.seq), 60):
                lines.append(str(fasta_rec.seq[i: i + 60] + "\n"))
            name = '>' + fasta_rec.description.strip() + '\n'
            handle.write(name)
            handle.write("".join(lines))
            handle.write('\n')
        handle.close()


# input_file = '/home/egors/Projects/CDS-alignment/genbank_sample_data/segment2.fasta'
# annot = '/home/egors/Projects/CDS-alignment/genbank_sample_data/annotation_segment2_test.tsv'
# align=True
# only_restore=False
# output_file = './out.fa'
# fasta = '/home/egors/Projects/CDS-alignment/gb_swine_NA/swine_sample.fasta'
# annot = '/home/egors/Projects/CDS-alignment/gb_swine_NA/swine_cds_annot_table.tsv'

# a = CdsAlignment(fasta, annot,do_align=False)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--fasta', type=str, help='input fasta', required=True, action='store')
    parser.add_argument('-o', '--out', type=str, help='output aligned fasta', required=True, action='store')
    parser.add_argument('-a', '--annot', type=str, help='cds annotation table', required=True, action='store')

    parser.add_argument('--align', dest='align', action='store_true', help='does split to UTRs and CDS and alignment')
    parser.add_argument('--no-align', dest='align', action='store_false', help=
    'skip alignment and assume we already have all files in <fasta_name>.part directory')
    parser.add_argument('--only-restore', dest='only_restore', action='store_true',
                        help='Looks for alignment in <fasta_name>.part dir and put aligned sequences together ')
    parser.add_argument("-v", "--verbose", help="stdout verbosity",
                        action="store_true")
    parser.set_defaults(align=True, only_restore=False, verbose=True)
    args = parser.parse_args()

    input_file = args.fasta
    output_file = args.out
    annot = args.annot
    align = args.align
    only_restore = args.only_restore

    args = parser.parse_args()
    if args.verbose:
        logger = logging.getLogger()
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            '%(asctime)s %(levelname)-8s %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.INFO)

    CdsAlignment(fasta_file=input_file, annotation_table=annot, output_file=output_file,
                 do_align=align, only_restore=only_restore)
