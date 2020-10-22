#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Purpose
-------



Expected input
--------------


Code documentation
------------------
"""


import os
import sys
import csv
import argparse
import subprocess
import datetime as dt

from Bio import SeqIO
from Bio.Seq import Seq


def reverse_complement(dna_sequence):
    """ Determines the reverse complement of a given DNA strand.

        Parameters
        ----------
        dna_sequence : str
            String representing a DNA sequence.

        Returns
        -------
        reverse_complement_strand : str
            The reverse complement of the DNA sequence (without
            lowercase letters).
    """

    base_complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A',
                       'a': 'T', 'c': 'G', 'g': 'C', 't': 'A',
                       'n': 'N', 'N': 'N'}

    # convert string into list with each character as a separate element
    bases = list(dna_sequence)

    # determine complement strand
    # default to 'N' if nucleotide is not in base_complement dictionary
    bases = [base_complement.get(base, 'N') for base in bases]

    complement_strand = ''.join(bases)

    # reverse strand
    reverse_complement_strand = reverse_str(complement_strand)

    return reverse_complement_strand


def reverse_str(string):
    """ Reverses character order in input string.

        Parameters
        ----------
        string : str
            String to be reversed.

        Returns
        -------
        revstr :str
            Reverse of input string.
    """

    revstr = string[::-1]

    return revstr


def translate_sequence(dna_str, table_id):
    """ Translate a DNA sequence using the BioPython package.

        Parameters
        ----------
        dna_str :str
            DNA sequence as string type.
        table_id : int
            Translation table identifier.

        Returns
        -------
        protseq : str
            Protein sequence created by translating
            the input DNA sequence.
    """

    myseq_obj = Seq(dna_str)
    protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

    return protseq


def translate_dna_aux(dna_sequence, method, table_id):
    """ Tries to translate an input DNA sequence in specified orientation
        and stores exceptions when the input sequence cannot be translated.

        Parameters
        ----------
        dna_sequence : str
            String representing a DNA sequence.
        method : str
            A string specifying the way the sequence will
            be oriented to attempt translation.
        table_id : int
            Translation table identifier.

        Returns
        -------
        List with following elements if translation is successful:
            protseq : str
                String representing the translated DNA sequence.
            myseq : str
                String representing the DNA sequence in the
                orientation used to translate it.
        Otherwise, returns string derived from captured exception.
    """

    myseq = dna_sequence
    # try to translate original sequence
    if method == 'original':
        try:
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse complement
    elif method == 'revcomp':
        try:
            myseq = reverse_complement(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse
    elif method == 'rev':
        try:
            myseq = reverse_str(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh
    # try to translate the reverse reverse complement
    elif method == 'revrevcomp':
        try:
            myseq = reverse_str(myseq)
            myseq = reverse_complement(myseq)
            protseq = translate_sequence(myseq, table_id)
        except Exception as argh:
            return argh

    return [protseq, myseq]


def check_str_alphabet(string, alphabet):
    """ Determine if a string only has characters from specified
        alphabet.

        Parameters
        ----------
        string : str
            Input string.
        alphabet : str
            String that has all characters from desired
            alphabet.

        Returns
        -------
        "True" if sequence only has characters from specified
        alphabet and string "ambiguous or invalid characters" if
        it any of its characters is not in the alphabet.
    """

    valid_chars = alphabet
    if all(n in valid_chars for n in string) is True:
        return True
    else:
        return 'ambiguous or invalid characters'


def check_str_multiple(string, number):
    """ Determine if length of input string is multiple of
        a specified number.

        Parameters
        ----------
        string : str
            Input string.
        number : int
            Integer that will be used to check if sequence
            length is multiple of.

        Returns
        -------
        "True" if the length of the sequence is a multiple of the
        specified number and "sequence length is not a multiple of number"
        if condition is not satisfied.
    """

    if len(string) % number == 0:
        return True
    else:
        return 'sequence length is not a multiple of {0}'.format(number)


def translate_dna(dna_sequence, table_id, min_len):
    """ Checks if sequence is valid and attempts to translate it,
        calling several functions to ensure that the sequence only has
        'ACTG', is multiple of 3 and that it can be translated in any of 4
        different orientations. Stores exceptions so that it is possible to
        understand the sequence could not be translated.

        Parameters
        ----------
        dna_sequence : str
            String representing a DNA sequence.
        table_id : int
            Translation table identifier.
        min_len : int
            Minimum accepted length. Sequences with length below
            this value will not be translated.

        Returns
        -------
        If the sequence can be translated:
            sequence : list
                A list with two elemets, the protein sequence
                and the DNA sequence in the correct orientation.
            coding_strand : str
                The strand orientation that had could be
                translated.
        Otherwise:
            exception_str : str
                A string containing the exceptions that
                determined that the sequence could not be translated.
    """

    original_seq = dna_sequence.upper()
    exception_collector = []
    strands = ['sense', 'antisense', 'revsense', 'revantisense']
    translating_methods = ['original', 'revcomp', 'rev', 'revrevcomp']

    # check if the string is DNA, without ambiguous bases
    valid_dna = check_str_alphabet(original_seq, 'ACTG')
    if valid_dna is not True:
        return valid_dna

    # check if sequence is multiple of three
    valid_length = check_str_multiple(original_seq, 3)
    if valid_length is not True:
        return valid_length

    # check if sequence is not shorter than the accepted minimum length
    if len(original_seq) < min_len:
        return 'sequence shorter than {0} nucleotides'.format(min_len)

    # try to translate in 4 different orientations
    # or reach the conclusion that the sequence cannot be translated
    i = 0
    translated = False
    while translated is False:
        sequence, exception_collector = retranslate(original_seq,
                                                    translating_methods[i],
                                                    table_id, strands[i],
                                                    exception_collector)

        i += 1
        if i == len(strands) or isinstance(sequence, list) is True:
            translated = True

    coding_strand = strands[i-1]

    # if the sequence could be translated, return list with protein and DNA
    # sequence in correct orientation
    if isinstance(sequence, list):
        return [sequence, coding_strand]
    # if it could not be translated, return the string with all exception
    # that were collected
    else:
        exception_str = ','.join(exception_collector)
        return exception_str


def retranslate(sequence, method, table_id, strands, exception_collector):
    """ Sends sequence for translation and collects exceptions when
        the sequence cannot be translated.

        Parameters
        ----------
        sequence : str
            String representing a DNA sequence.
        method : str
            A string specifying the sequence orientation
            that should be used to attempt translation.
        table_id : int
            Translation table identifier.
        strands : list
            List with 4 different orientations that can
            be checked.
        exception_collector : list
            List used to store all exceptions arising from
            translation attempts.

        Returns
        -------
        A list with following elements:
            translated_seq : list
                A list with the protein sequence and with the DNA
                sequence in the orientation used for translation.
            exception_collector : list
                A list with the exceptions that are captured when
                the sequence could not be translated.
        Otherwise:
            translated_seq : str
                A string with the exception/reason why the sequence
                could not be translated.
            exception_collector : list
                List with all exception that have been captured during
                translation attempts of the current sequence.
    """

    translated_seq = translate_dna_aux(sequence, method, table_id)
    if not isinstance(translated_seq, list):
        exception_collector.append('{0}({1})'.format(strands,
                                                     translated_seq.args[0]))

    return [translated_seq, exception_collector]


def execute_augustus(input_file, species, output_file):
    """ Executes AUGUSTUS to predict genes in the input
        genome.

        Parameters
        ----------
        input_file : str
            Path to the input FASTA file with genome contigs.
        species : str
            Identifier of the species passed to choose
            prediction model used by AUGUSTUS.
        output_file : str
            Path to the output file that will store stdout
            from AUGUSTUS.
    """

    out_handle = open(output_file, 'w')
    proc = subprocess.call(['augustus', '--species={0}'.format(species),
                            '--gff3=on', '--stopCodonExcludedFromCDS=false',
                            '--strand=both', '--genemodel=complete',
                            '--noInFrameStop=true', input_file],
                           stdout=out_handle,
                           stderr=subprocess.PIPE)
    out_handle.close()

    # return exit code
    return proc


def make_blast_db(input_fasta, output_path, db_type):
    """ Creates a BLAST database.

        Parameters
        ----------
        input_fasta : str
            Path to the input file with sequences.
        output_path : str
            Path to the output database.
        db_type : str
            Type of the database, nucleotide (nuc) or
            protein (prot).

        Returns
        -------
        Creates a BLAST database with the input sequences.
    """

    blastdb_cmd = ['makeblastdb', '-in', input_fasta, '-out', output_path,
                   '-parse_seqids', '-dbtype', db_type]

    makedb_cmd = subprocess.Popen(blastdb_cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)
    # stdout = makedb_cmd.stdout.readlines()
    # stderr = makedb_cmd.stderr.readlines()
    # print(stdout, stderr)

    makedb_cmd.wait()


def run_blast(blastp_path, blast_db, fasta_file, blast_output,
              max_hsps=1, threads=1, ids_file=None):
    """ Runs BLASTp.

        Parameters
        ----------
        blastp_path : str
            Path to BLASTp executables.
        blast_db : str
            Path to the BLASTdb to use.
        fasta_file : str
            Path to the FASTA file that will be passed as query.
        blast_output : str
            Path to the file that will store the results.
        max_hsp : int
            Maximum number o High Scoring Pairs to determine when
            comparing each query/subject.
        threads : int
            Number of threads to use to run BLASTp.
        ids_file : str
            Path to a file with a set of identifiers for the
            sequences that should be used as subjects.

        Returns
        -------
        stderr : str
            String with information about any errors that ocurred.
    """

    blast_args = [blastp_path, '-db', blast_db, '-query', fasta_file,
                  '-out', blast_output, '-outfmt', '6 qseqid sseqid score',
                  '-max_hsps', str(max_hsps), '-num_threads', str(threads),
                  '-evalue', '0.001']

    if ids_file is not None:
        blast_args.extend(['-seqidlist', ids_file])

    blast_proc = subprocess.Popen(blast_args,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    stderr = blast_proc.stderr.readlines()

    return stderr


def read_blast_tabular(blast_tabular_file):
    """ Read a file with BLAST results in tabular format.

        Parameters
        ----------
        blast_tabular_file : str
            Path to output file of BLAST.

        Returns
        -------
        blasting_results : list
            A list with a sublist per line in the input file.
    """

    with open(blast_tabular_file, 'r') as blastout:
        reader = csv.reader(blastout, delimiter='\t')
        blasting_results = [row for row in reader]

    return blasting_results


def apply_bsr(inputs):
    """ Computes the BSR value for each alignment in a list of
        alignments and selects or excludes sequence identifers
        based on BSR value and sequence length.

        Parameters
        ----------
        inputs : list
            A list with the following elements:
                - Path to the file with BLASTp results.
                - Path to the FASTA file with sequences that were
                  aligned.
                - BSR value used as threshold.

        Returns
        -------
        excluded_alleles : list
            A list with the identifiers of the sequences that were
            excluded based on the BSR threshold and length values.
    """

    blast_file = inputs[0]
    blast_results = read_blast_tabular(blast_file)
    self_scores = {r[0]: r[2] for r in blast_results if r[0] == r[1]}
    # do not include self-scores lines, no need to evaluate those hits
    blast_results = [r for r in blast_results if r[0] != r[1]]

    fasta_file = inputs[1]
    lengths = {}
    for k in self_scores:
        record = fasta_file.get(k)
        sequence = str(record.seq)
        lengths[k] = len(sequence)

    bsr = inputs[2]

    excluded_alleles = []
    for res in blast_results:

        query = res[0]
        hit = res[1]
        score = res[2]

        if query not in excluded_alleles:
            # try to apply BSR strategy
            try:
                self_blast_score = self_scores[query]

                query_length = lengths[query]
                hit_length = lengths[hit]
                blast_score_ratio = float(score) / float(self_blast_score)

                # BSR has to be greater than threshold, just as in the original function
                if blast_score_ratio >= bsr and hit not in excluded_alleles:

                    if hit_length > query_length and query not in excluded_alleles:
                        excluded_alleles.append(query)

                    elif hit_length <= query_length:
                        excluded_alleles.append(hit)
            # it might not work because there is no self score for
            # some sequences due to low complexity regions...
            except Exception:
                excluded_alleles.append(query)

    return excluded_alleles


def main(input_file, output_path, species, bsr, threads):

    print('\n{0}\n{1}\n{0}'.format('-'*9, ' euTyper'))

    start_date = dt.datetime.now()
    start_date_str = dt.datetime.strftime(start_date, '%Y-%m-%dT%H:%M:%S')
    print('\nStarted at: {0}\n'.format(start_date_str))

    print('Number of genomes/assemblies: {0}'.format(1))
    print('Species: {0}'.format(species))
    print('BLAST Score Ratio: {0}'.format(bsr))
    print('Number of threads: {0}'.format(threads))

    augustus_outfile = os.path.join(output_path, 'augustus_results.gff')
    print('\nRunning AUGUSTUS...')
    exit_code = execute_augustus(input_file, species, augustus_outfile)
    if exit_code != 0:
        sys.exit('AUGUSTUS returned exit code != 0. Exited.')
    print('Finished running AUGUSTUS.')

    # import input genomes
    chrs = {rec.id: str(rec.seq).upper()
            for rec in SeqIO.parse(input_file, 'fasta')}

    print('Extracting CDSs from AUGUSTUS results and creating '
          'full gene sequences...')
    # import AUGUSTUS results
    with open(augustus_outfile, 'r') as af:
        augustus_lines = af.readlines()

    # start getting genes
    genes = {}
    gene_id = 0
    for l in augustus_lines:
        if '# start gene' in l:
            gene_id += 1
            genes[gene_id] = []
        elif 'CDS' in l and '#' not in l:
            cds_line = l.split('\t')
            genes[gene_id].append((cds_line[0], cds_line[2],
                                   cds_line[3], cds_line[4],
                                   cds_line[6], cds_line[7]))

    # the last value is the codon phase information this information
    # does not mean that we should discard 0, 1 or 2 bases from the
    # start of the sequence. It means that 0, 1 or 2 bases of the sequence
    # will complete the last codon from the previous CDS and the first complete
    # codon of the new CDS only starts after those initial bases.
    # So, a CDS might not have length that is a multiple of 3 and will not
    # end in a complete codon. We need to concatenate all gene CDSs to obtain
    # the full sequence that can be translated.

    # get cds sequences for each gene
    gene_cds = {}
    for gid, ginfo in genes.items():
        cds_sequences = []
        strands = []
        for c in ginfo:
            start = (int(c[2])-1)
            end = int(c[3])
            cds = chrs[c[0]][start:end]
            if c[4] != '+':
                cds = reverse_complement(cds)

            cds_sequences.append(cds)
            strands.append(c[4])

        gene_cds[gid] = (cds_sequences, list(set(strands)))

    # join CDS sequences for each gene and translate
    # must concatenate sequences in the right order and orientation
    gene_dna_prot = {}
    for gid, gcds in gene_cds.items():
        strand = gcds[1]
        if len(strand) == 1:
            if strand[0] == '+':
                dna_seq = ''.join(gcds[0])
            elif strand[0] == '-':
                dna_seq = ''.join(gcds[0][::-1])

            prot_seq = translate_dna(dna_seq, 1, 0)
            gene_dna_prot[gid] = [dna_seq, str(prot_seq[0][0])]
        elif len(strand) > 1:
            print(gid)

    # save DNA and Protein sequences to file
    dna_records = ['>{0}\n{1}'.format(k, v[0])
                   for k, v in gene_dna_prot.items()]
    protein_records = ['>{0}\n{1}'.format(k, v[1])
                       for k, v in gene_dna_prot.items()]

    proteins_file = os.path.join(output_path, 'proteins.fasta')
    with open(proteins_file, 'w') as pf:
        pf.write('\n'.join(protein_records))

    dna_file = os.path.join(output_path, 'dna.fasta')
    with open(dna_file, 'w') as df:
        df.write('\n'.join(dna_records))

    print('Aligning translated genes with BLASTp...')
    # create BLASTdb
    blastdb = os.path.join(output_path, 'prots')
    make_blast_db(proteins_file, blastdb, 'prot')

    # align proteins with BLASTp
    blast_output = os.path.join(output_path, 'blastout.tsv')
    blasterr = run_blast('blastp', blastdb, proteins_file, blast_output,
                         max_hsps=1, threads=threads, ids_file=None)

    # index FASTA file with DNA sequences
    indexed_fasta = SeqIO.index(dna_file, 'fasta')

    # apply BSR
    excluded_ids = apply_bsr([blast_output, indexed_fasta, bsr])
    print('Excluded {0} sequences highly similar to larger '
          'sequences.'.format(len(excluded_ids)))

    # determine distinct set of gene ids
    distinct_ids = [i for i in gene_dna_prot.keys()
                    if str(i) not in excluded_ids]
    
    print('Constructing final schema structure...')
    schema_dir = os.path.join(output_path, 'schema_seed')
    schema_short_dir = os.path.join(schema_dir, 'short')
    os.makedirs(schema_short_dir)

    # save files
    for seqid in distinct_ids:
        locus_file = 'gene{0}.fasta'.format(seqid)
        schema_file = os.path.join(schema_dir, locus_file)
        locus_short_file = 'gene{0}_short.fasta'.format(seqid)
        schema_short_file = os.path.join(schema_short_dir, locus_short_file)

        header = '>gene{0}_1'.format(seqid)
        seq = str(indexed_fasta.get(str(seqid)).seq)

        record = '\n'.join([header, seq]) + '\n'

        with open(schema_file, 'w') as sf:
            sf.write(record)

        with open(schema_short_file, 'w') as ssf:
            ssf.write(record)

    end_date = dt.datetime.now()
    end_date_str = dt.datetime.strftime(end_date, '%Y-%m-%dT%H:%M:%S')

    delta = end_date - start_date
    minutes, seconds = divmod(delta.total_seconds(), 60)

    print('\nFinished at: {0}'.format(end_date_str))
    print('Created schema with {0} genes based on {1} genome in'
          '{2: .0f}m{3: .0f}s.'.format(len(distinct_ids), 1,
                                       minutes, seconds))


def parse_arguments():

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-i', type=str, required=True,
                        dest='input_files',
                        help='')

    parser.add_argument('-o', type=str, required=True,
                        dest='output_path',
                        help='')

    parser.add_argument('-s', type=str, required=True,
                        dest='species',
                        help='')

    parser.add_argument('--bsr', type=float, required=False,
                        default=0.6, dest='bsr',
                        help='')
    
    parser.add_argument('--t', type=int, required=False,
                        default=1, dest='threads',
                        help='')

    args = parser.parse_args()

    return [args.input_files, args.output_path, args.species,
            args.bsr, args.threads]


if __name__ == '__main__':

    args = parse_arguments()

    main(args[0], args[1], args[2], args[3], args[4])
