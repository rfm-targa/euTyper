#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
AUTHOR

    Rafael Mamede
    github: @rfm-targa

DESCRIPTION



"""


import os
import sys
import time
import argparse
import subprocess

from Bio import SeqIO
from Bio.Seq import Seq


def reverse_complement(dna_sequence):
    """ Determines the reverse complement of given DNA strand.

        Args:
            strDNA (str): string representing a DNA sequence.

        Returns:
            revC_dna (str): the reverse complement of the DNA sequence, without
            lowercase letters.

        Example:
            >>> reverse_complement('ATCGgcaNn')
            'NNTGCCGAT'
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
    """ Reverse character order in input string.

        Args:
            string (str): string to be reversed.

        Returns:
            revstr (str): reverse of input string.
    """

    revstr = string[::-1]

    return revstr


def translate_sequence(dna_str, table_id):
    """ Translate a DNA sequence using the BioPython package.

        Args:
            dna_str (str): DNA sequence as string type.
            table_id (int): translation table identifier.

        Returns:
            protseq (str): protein sequence created by translating
            the input DNA sequence.
    """

    myseq_obj = Seq(dna_str)
    protseq = Seq.translate(myseq_obj, table=table_id, cds=True)

    return protseq


def translate_dna_aux(dna_sequence, method, table_id):
    """ Tries to translate an input DNA sequence in specified orientation
        and stores exceptions when the input sequence cannot be translated.

        Args:
            dna_sequence (str): string representing a DNA sequence.
            method (str): a string specifying the way the sequence will
            be oriented to attempt translation.
            table_id (int): translation table identifier.

        Returns:
            List with following elements if translation is successful:
                protseq (str): string representing the translated DNA sequence.
                myseq (str): string representing the DNA sequence in the
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

        Args:
            string (str): input string.
            alphabet (str): string that has all characters from desired
            alphabet.

        Returns:
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

        Args:
            string (str): input string.
            number (int): integer that will be used to check if sequence
            length is multiple of.

        Returns:
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

        Args:
            dna_sequence (str):
            table_id (int):

        Returns:
            If the sequence can be translated,
            a list with following elements:
                sequence (list): a list with two elemets, the protein sequence
                and the DNA sequence in the correct orientation.
                coding_strand (str): the strand orientation that had could be
                translated.
            Otherwise:
                exception_str (str): a string containing the exceptions that
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

        Args:
            sequence (str): string representing a DNA sequence.
            method (str): a string specifying the sequence orientation
            that should be used to attempt translation.
            table_id (int): translation table identifier.
            strands (list): list with 4 different orientations that can
            be checked.
            exception_collector (list): list used to store all exceptions
            arising from translation attempts.

        Returns:
            A list with following elements, if the sequence can be translated:
                translated_seq (list): a list with the protein sequence and
                with the DNA sequence in the orientation used for translation.
                exception_collector (list): a list with the exceptions that are
                captured when the sequence could not be translated.
            Otherwise:
                translated_seq (str): a string with the exception/reason why
                the sequence could not be translated.
                exception_collector (list): list with all exception that have
                been captured during translation attempts of the current
                sequence.
    """

    translated_seq = translate_dna_aux(sequence, method, table_id)
    if not isinstance(translated_seq, list):
        exception_collector.append('{0}({1})'.format(strands,
                                                     translated_seq.args[0]))

    return [translated_seq, exception_collector]


def execute_augustus(input_file, species, output_file):
    """
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


def main(input_file, output_path, species):

    augustus_outfile = os.path.join(output_path, 'augustus_results.gff')
    exit_code = execute_augustus(input_file, species, augustus_outfile)
    if exit_code != 0:
        sys.exit('AUGUSTUS returned exit code != 0. Exited.')

    # import input genomes
    chrs = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(input_file, 'fasta')}

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

    # the last value is the codon phase information
    # this information does not mean that we should discard 0, 1 or 2 bases
    # from the start of the sequence. It means that 0, 1 or 2 bases of the sequence
    # will complete the last codon from the previous CDS and the first complete
    # codon of the new CDS only starts after those initial bases.
    # So, a CDS might not have length that is a multiple of 3 and will not end in a complete
    # codon. We need to concatenate all gene CDSs to obtain the full sequence that can be
    # translated.

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
    # not working for sequences in antisense strand
    # I might be concatenating CDS in the antisense as if they were
    # in the sense strand...noob me :'|
    # I have to concatenate sequences in the right order and orientation!
    # Seems to work fine now!!! :)
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

    schema_dir = os.path.join(output_path, 'schema_seed')
    schema_short_dir = os.path.join(schema_dir, 'short')
    os.makedirs(schema_short_dir)

    # save files
    for locus, seqs in gene_dna_prot.items():
        locus_file = 'gene{0}.fasta'.format(locus)
        schema_file = os.path.join(schema_dir, locus_file)
        locus_short_file = 'gene{0}_short.fasta'.format(locus)
        schema_short_file = os.path.join(schema_short_dir, locus_short_file)

        header = '>gene{0}_1'.format(locus)
        seq = seqs[0]

        record = '\n'.join([header, seq]) + '\n'

        with open(schema_file, 'w') as sf:
            sf.write(record)

        with open(schema_short_file, 'w') as ssf:
            ssf.write(record)

    print('Created schema seed!')


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

    args = parser.parse_args()

    return [args.input_files, args.output_path, args.species]


if __name__ == '__main__':

    args = parse_arguments()

    main(args[0], args[1], args[2])
