#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#######################
# euTyper main script #
#######################

#################
# import modules#
#################


import datetime as dt
from Bio.Seq import Seq
import os
import sys
import argparse
import itertools
from multiprocessing import Pool
import time
import traceback
from collections import Counter
import platform
import subprocess
from Bio import SeqIO
import csv
from Bio.SeqIO import FastaIO
import shutil


#################
# Pre Processing#
#################

def retrieve_sample_data(sample_metadata_path):
    """ Retrieve file names from sample metadata.

        Parameters
        ----------
        - path to sample meta data

        Returns: list of sample names (genomes / assemblies)
    """
    sample_information = {}

    with open(sample_metadata_path) as file:
        for line in file:
            line = line.rstrip().split('\t')
            sample_information[line[0]] = (line[1],line[2])

    return(sample_information)

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

    coding_strand = strands[i - 1]

    # if the sequence could be translated, return list with protein and DNA
    # sequence in correct orientation
    if isinstance(sequence, list):
        return [sequence, coding_strand]
    # if it could not be translated, return the string with all exception
    # that were collected
    else:
        exception_str = ','.join(exception_collector)
        return exception_str

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

def execute_augustus(sample_information, threads, output_path, pt_to_input_folder):
    """ Executes AUGUSTUS to predict genes in the input
        genome.

        Parameters
        ----------
        input_file : str
            Path to sample metadata containing required data.
        Threads : int
            amount of CPU the augustus parralelisation is allowed to use.
        output_file : str
            Path to the output file that will store stdout
            from AUGUSTUS.
        path_to_input_folder : str
            path to directory containing input folders

        Returns
        ----------
        temporary directory paths
    """
    #######################################################
    # preparing for parralisation to reduce computing time#
    #######################################################

    temp_dir = os.path.join(output_path,'temp')
    os.makedirs(temp_dir, exist_ok=True)
    
    temp_dir_split = os.path.join(temp_dir,'1_splitted_fasta')
    os.makedirs(temp_dir_split, exist_ok=True)


    print('Preparing fasta files for parralelisation')
    # split all fasta fies:
    for sample in sample_information:
        sample_path = os.path.join(pt_to_input_folder, sample)
        split_Mfasta(sample_path, temp_dir_split, sample)

    files = os.listdir(temp_dir_split)


    # create Augustus job list
    temp_dir_2, file_to_delete = create_augustus_job_lists(files, temp_dir_split, temp_dir, sample_information)

    print('Start performing Augustus parralel')
    os.system("parallel -j"+ str(threads) +" --bar --no-notice 'nice ./{}' < job.lst")

    #delete files
    for file in file_to_delete:
        os.remove(file)

    return(temp_dir_2, temp_dir)

def create_augustus_job_lists(files, temp_dir_split, temp_dir, sample_information):
    '''
    create augustus job lists to enable paralelisation

    Parameters
    ----------
    files: list
        splitted fasta files

    temp_dir_split: str
        path of the splitted fasta files

    temp_dir: str
        path to general tmp directory

    sample_information: dict
        dic containing for each sample the species used by augustus

   return:
    ---------------
    tmp_dir_augustus: str
        path to tmp dir with augustus files

    files: list
    list with files who are created
    '''
    files_created = []
    count = 0
    files = sorted(files)
    temp_dir_augustus_output = os.path.join(temp_dir,'2_augustus_output')
    print(sample_information)
    os.makedirs(temp_dir_augustus_output, exist_ok=True)
    print(files)
    for file in list(files):
            try:
                count +=1
                species = sample_information[(file.split('_split_'))[0]][0]
                augcall = 'augustus --species=' + species + ' --gff3=on  --stopCodonExcludedFromCDS=false --strand=both --genemodel=complete --noInFrameStop=true'
                job_name = (('Augustus_job_' +str(count) + '_' + file[:-3])).replace('.','_')
                input_file = os.path.join(temp_dir_split,file)
                output_file = os.path.join(temp_dir_augustus_output,(file +'.gff'))
                errfile = os.path.join(temp_dir_augustus_output, (file + '.err'))

                with open(job_name,'w+') as f:
                    f.write('#\n')
                    f.write((augcall + ' ' + input_file + ' --outfile=' + output_file +' --errfile=' + errfile))

                    files_created.append(job_name)

            except:
                print('##################################')
                print('sample {0} could not be found in the sample meta data \n'
                      'this is probably due contig/scaffold name change in fasta file'.format(file))
                print('####################################')

    for file in files_created:
        os.system("chmod +x {0}".format(file))


    with open('job.lst', 'w') as f:
        for file in files_created:
            f.write((file) +'\n')

    files_created.append('job.lst')

    return (temp_dir_augustus_output, files_created)

def split_Mfasta(sample_path,temp_dir_split, sample):
    '''split the fasta files into their contigs
    which allows for optimized parralelisation of Augustus

    sample path: str
        path to directory containing the fasta files

    tmp dir split: str
        path to temporary directory for the splitted fasta files to be stored

    sample: str
        unique sample name

    :return
    ----------
    NA


    '''

    count = 0
    file = False
    with open(sample_path) as fasta_in:
        for line in fasta_in:
            if '>' in line:
                if file:
                    fasta_out.close()
                count +=1
                out = os.path.join(temp_dir_split,(sample + "_split_" + str(count) +'.fa'))
                fasta_out = open(out,'w')
                fasta_out.write(line)
                file = True
            else:
                fasta_out.write(line)

def retrieve_CDS_AA_from_augustus(temp_dir_2, sample_information ,temp_dir, output_path, prob):
    '''
    retrieves the CDS and AA sequences generated from the augustus files.
    It will be stored to a dict with each gene identifier containing the CDS and AA sequence

    :parameter
    ----------------
    :param temp_dir_2: str
        path to splitted fasta files
    :param sample_information: dict
        dict containg each sample, the species, and codon table nr used to calculate
        the protein sequence (AA)

    :param temp_dir: str
        path to general temporary directory

    :param output_path: str
        path to temporary directory where all splitted fasta will be stored

    :param prob: int
        cutoff for the probabillity of augustus used as cutoff (default = 0.8)

    :return:
    -----------------
    path where splitted fasta files will be stored in tmp directory
    '''
    os.makedirs(os.path.join(output_path, 'CDS_AA'), exist_ok=True)
    #############################################
    # linking augustus output to splitted genome#
    #############################################
    # key = splitted fa files
    # value is linked augustus file
    files = os.listdir(temp_dir_2)
    combo = {}
    names = sample_information.keys()
    for name in names:
        combo[name] = []
        for file in files:
            if name != file and name in file and file[-3:] != 'err':
                combo[name].append(file)

    print('Extracting CDSs from AUGUSTUS results and creating '
          'full gene sequences...')

    # sort list to prevent reassigning segID
    for key in combo.keys():
        list_files = combo[key]
        combo[key] = sorted(list_files)
    ################################################
    # retrieving CDS
    ################################################

    for key in combo.keys():
        Mgid = 0
        gene_id = 0
        files = combo[key]
        print("#####################")
        print(combo)
        for value in files:
            genes = {}
            # import input genomes
            chrs = {rec.id: str(rec.seq).upper()
                    for rec in SeqIO.parse((os.path.join(temp_dir, '1_splitted_fasta', value))[:-4], 'fasta')}

            # import AUGUSTUS results
            with open(os.path.join(temp_dir_2, value), 'r') as af:
                augustus_lines = af.readlines()

                # start getting genes
                # Filter genes on overall probabillity check.
                # at start gene, proabillity = false. when overall gene probabillity is above cutoff it will be True
                # CDS will be only noted down when probabillity is above cutoff
                for l in augustus_lines:
                    if '# start gene' in l:
                        probabillity_check = False

                    # add proabillity score to filter later by
                    elif 'gene' in l and '#' not in l and float((l.split('\t'))[5]) >= prob:
                        gene_id += 1
                        genes[gene_id] = []
                        probabillity_check = True
                    elif 'CDS' in l and '#' not in l and probabillity_check:
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
                    start = (int(c[2]) - 1)
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
                Mgid += 1
                strand = gcds[1]
                if len(strand) == 1:
                    if strand[0] == '+':
                        dna_seq = ''.join(gcds[0])
                    elif strand[0] == '-':
                        dna_seq = ''.join(gcds[0][::-1])


                    prot_seq = translate_dna(dna_seq, (sample_information[key])[1] , 0)
                    gene_dna_prot[Mgid] = [dna_seq, str(prot_seq[0][0])]
                elif len(strand) > 1:
                    print(gid)

            # save DNA and Protein sequences to file
            dna_records = ['>{0}\n{1}'.format(k, v[0])
                           for k, v in gene_dna_prot.items()]

            protein_records = ['>{0}\n{1}'.format(k, v[1])
                               for k, v in gene_dna_prot.items()]

            os.makedirs(os.path.join(os.path.join(temp_dir,'3_splitted_fastas')), exist_ok=True)

            proteins_file = os.path.join(os.path.join(temp_dir,'3_splitted_fastas'), (value + '_proteins.fasta'))
            with open(proteins_file, 'w') as pf:
                pf.write('\n'.join(protein_records))

            dna_file = os.path.join(os.path.join(temp_dir,'3_splitted_fastas'), (value + '_dna.fasta'))
            with open(dna_file, 'w') as df:
                df.write('\n'.join(dna_records))

    return(os.path.join(temp_dir,'3_splitted_fastas'))

def combine_fastas(path_to_tmp3, sample_information, output_path):
    '''

    combine all splitted fasta files beloning to the same species

    parameters
    ---------
    :param path_to_tmp3: str
        path to tmp directory where splitted fasta will be stored

    :param sample_information: dict
        dictornary containing the sample name, species and codon table nr.

    :param output_path:
        path to CDS_AA where all combined fasta file per species will be stored


    :return:
        path to CDS_AA directory
    '''
    print('Combining created fasta files')
    files = os.listdir(path_to_tmp3)

    ###########################
    # find corrosponding files#
    ###########################
    for key in sample_information.keys():
        proteins = []
        CDS = []
        for file in files:
            if key in file and 'proteins' in file:
                proteins.append(file)
            elif key in file and 'dna' in file:
                CDS.append(file)

        # sorting algorithm do
        proteins = sorted(proteins)
        CDS = sorted(CDS)
    ##########################
    # combine the fasta files#
    ##########################

        output = os.path.join(output_path, 'CDS_AA', (key +'_protein'))
        with open(output, 'w') as f_out:
            for file in proteins:
                file = os.path.join(path_to_tmp3,file)
                with open(file) as f:
                    for line in f:
                        if "\n" in line:
                            f_out.write((line))
                        else:
                            f_out.write((line + "\n"))

        output = os.path.join(output_path, 'CDS_AA', (key +'_CDS'))
        with open(output, 'w') as f_out:
            for file in CDS:
                file = os.path.join(path_to_tmp3,file)
                with open(file) as f:
                    for line in f:
                        if "\n" in line:
                            f_out.write((line))
                        else:
                            f_out.write((line + "\n"))

    return(os.path.join(output_path, 'CDS_AA'))




############
# universal#
############

def replace_multiple_characters(input_string, replacements):
    """ Replaces multiple characters in a string.
        Parameters
        ----------
        input_string : str
            String with characters to be replaced.
        Returns
        -------
        replaced : str
            Input string without replaced characters.
    """

    for r in replacements:
        if r[0] in input_string:
            input_string = input_string.replace(*r)

    return input_string

def split_fasta(fasta_path, output_path, num_seqs, filenames):
    """ Splits a FASTA file.
        Parameters
        ----------
        fasta_path : str
            Path to a FASTA file.
        output_path : str
            Path to the output directory where new FASTA
            files will be created.
        num_seqs : int
            Split FASTA file into files with this number
            of sequences.
        filenames : gen
            Generator with names to attribute to new files.
        Returns
        -------
        splitted_files : list
            List with paths to the new files that were
            created by splitting the input FASTA file.
    """
    CHAR_REPLACEMENTS = [("|", "_"), ("_", "-"), ("(", ""),
                         (")", ""), ("'", ""), ("\"", ""), (":", "")]

    splitted_files = []
    current_recs = []
    records = [rec for rec in SeqIO.parse(fasta_path, 'fasta')]
    for record in records:
        current_recs.append(record)
        if len(current_recs) == num_seqs or record.id == records[-1].id:
            file_name = filenames.__next__()
            file_name = replace_multiple_characters(file_name, CHAR_REPLACEMENTS)

            new_file = join_paths(output_path,
                                  ['{0}{1}'.format(file_name, '.fasta')])

            splitted_files.append(new_file)

            write_records(current_recs, new_file)

            current_recs = []

    return splitted_files

def replace_multiple_characters(input_string, replacements):
    """ Replaces multiple characters in a string.
        Parameters
        ----------
        input_string : str
            String with characters to be replaced.
        Returns
        -------
        replaced : str
            Input string without replaced characters.
    """

    for r in replacements:
        if r[0] in input_string:
            input_string = input_string.replace(*r)

    return input_string

def write_records(records, output_file):
    """ Writes FASTA records (BioPython SeqRecord) to a file.
        Parameters
        ----------
        records : list
            List with BioPython SeqRecord objects.
        output_file : str
            Path to the output file.
    """

    with open(output_file, 'w') as output_handle:
        fasta_out = FastaIO.FastaWriter(output_handle, wrap=None)
        fasta_out.write_file(records)

def create_directory(directory_path):
    """ Creates a diretory if it does not exist."""

    if not os.path.exists(directory_path):
        os.makedirs(directory_path)

def read_tabular(input_file, delimiter='\t'):
    """ Read tabular file.
        Parameters
        ----------
        input_file : str
            Path to a tabular file.
        delimiter : str
            Delimiter used to separate file fields.
        Returns
        -------
        lines : list
            A list with a sublist per line in the input file.
            Each sublist has the fields that were separated by
            the defined delimiter.
    """

    with open(input_file, 'r') as infile:
        reader = csv.reader(infile, delimiter=delimiter)
        lines = [line for line in reader]

    return lines

def retrieve_CDS_AA(paths):
    """ Creates a dictionary with the CDS coupled to their AA sequence

         Parameters
         ----------
        paths: dictionary of genomes:[ CDS file, protein file

        :return
        ---------
        CDS_AA dict {gene ID (species_name + gene number) : CDS, AA}
        :
     """
    CDS_AA = {}
    counterA = 0
    for sample in paths:
        counter = 0
        with open(paths[sample][0]) as CDS:
            for line in CDS:
                if ">" in line:
                    tmp_id = ((line[1:]).rstrip() + "_" + sample)
                    CDS_AA.setdefault(tmp_id, [])
                else:
                    CDS_AA[tmp_id].append(line.rstrip())

        with open(paths[sample][1]) as AA:
            for line in AA:
                if ">" in line:
                    tmp_id = ((line[1:]).rstrip() + "_" + sample)

                else:
                    CDS_AA[tmp_id].append(line.rstrip())

        #small check if CDS_AA is correct
        for i in CDS_AA.values():
            counter += 1
            counterA += 1

            if len(i) != 2:
                print('WARNING SOME PROTEINS COULD NOT BE COUPLED TO THEIR CORROSPONDING SEQ')

        print('sample {0} contains {1} sequentions \n'.format(sample, counter))
    print (' in total {0} sequences has been found'.format(counterA))

    return (CDS_AA)

def remove_duplicates(CDS_AA):
    '''
    removes duplications in the CDS_AA DIC to reduce computing time

    parameters
    :param CDS_AA: dict
            CDS_AA dict {gene ID (species_name + gene number) : CDS, AA}

    :return:
            CDS_AA dict {gene ID (species_name + gene number) : CDS, AA}

    '''
    len1 = len(CDS_AA)
    # deletes duplicating DNA seq and AA seq out of item do reduce amount of computing time
    toPOP = []
    CDS = {}
    AA = {}
    # Duplicate items will not be set in the new dic
    for ID, value in CDS_AA.items():
        if value[0] not in CDS.values():
            CDS[ID] = value[0]
        else:
            toPOP.append(ID)

    for item in toPOP:
        CDS_AA.pop(item)

    toPOP = []

    # removing duplicate AA
    for ID, value in CDS_AA.items():
        if value[1] not in AA.values():
            AA[ID] = value[1]
        else:
            toPOP.append(ID)

    for item in toPOP:
        CDS_AA.pop(item)

    len2 = len(CDS_AA)
    print(' finisht removing of ' + str(len1 - len2) + ' duplications')

    return (CDS_AA)

def remove_small(CDS_AA, MinSeqlen):
    '''
    removes CDS sequences who are smaller then 201 nucleotides

    :param CDS_AA: dict
                CDS_AA dict {gene ID (species_name + gene number) : CDS, AA}

    :param MinSeqlen: int
        minimal sequence length to be not filterd out

    :return:
    ----------
         CDS_AA dict {gene ID (species_name + gene number) : CDS, AA}

    '''
    toPOP = []

    # goes through every DNA sequence. if smaller it ID will be added to pop list
    for ID, value in CDS_AA.items():
        if len(value[0]) < MinSeqlen:
            toPOP.append(ID)

    for item in toPOP:
        CDS_AA.pop(item)

    return (CDS_AA)
    print(' finisht removing of ' + str(len(toPOP)) + ' sequences')

def join_list(lst, link):
    """ Joins all elements in a list into a single string.
        Parameters
        ----------
        lst : list
            List with elements to be joined.
        link : str
            Character used to join list elements.
        Returns
        -------
        joined_list : str
            A single string with all elements in the input
            list joined by the character chosen as link.
    """

    joined_list = link.join(lst)

    return joined_list

def write_to_file(text, output_file, write_mode, end_char):
    """ Writes a single string to a file.
        Parameters
        ----------
        text : str
            A single string to write to the output file.
        output_file : str
            Path to the output file.
        write_mode : str
            Write mode can be 'w', writes text and overwrites
            any text in file, or 'a', appends text to text
            already in file.
        end_char : str
            Character added to the end of the file.
    """

    with open(output_file, write_mode) as out:
        out.write(text + end_char)

def map_async_parallelizer(inputs, function, cpu, callback='extend',
                           chunksize=1, show_progress=True):
    """ Parallelizes function calls by creating several processes
        and distributing inputs.
        Parameters
        ----------
        inputs : list
            List with inputs to process.
        function
            Function to be parallelized.
        cpu : int
            Number of processes to create (based on the
            number of cores).
        callback : str
            Results can be appended, 'append', to the
            list that stores results or the list of results
            can be extended, 'extend'.
        chunksize : int
            Size of input chunks that will be passed to
            each process. The function will create groups
            of inputs with this number of elements.
        show_progress: bool
            True to show a progress bar with the percentage
            of inputs that have been processed, False
            otherwise.
        Returns
        -------
        results : list
            List with the results returned for each function
            call.
    """

    results = []
    pool = Pool(cpu)
    if callback == 'extend':
        rawr = pool.map_async(function, inputs,
                              callback=results.extend, chunksize=chunksize)
    elif callback == 'append':
        rawr = pool.map_async(function, inputs,
                              callback=results.append, chunksize=chunksize)

    if show_progress is True:
        completed = False
        while completed is False:
            completed = progress_bar(rawr, len(inputs))

    rawr.wait()

    return results
def progress_bar(process, total, tickval=5, ticknum=20, completed=False):
    """ Creates and prints progress bar to stdout.
        Parameters
        ----------
        process : multiprocessing.pool.MapResult
            Multiprocessing object.
        total : int
            Total number of inputs that have to be processed.
        tickval : int
            Progress completion percentage value for each
            tick.
        ticknum : int
            Total number of ticks in progress bar.
        completed : bool
            Boolean indicating if process has completed.
        Returns
        -------
        completed : bool
            Boolean indicating if process has completed.
    """

    # check if process has finished
    if (process.ready()):
        # print full progress bar and satisfy stopping condition
        progress_bar = '[{0}] 100%'.format('=' * ticknum)
        completed = True

    # check how many inputs have been processed
    remaining = process._number_left
    if remaining == total:
        # print empty progress bar
        progress_bar = '[{0}] 0%'.format(' ' * ticknum)
    else:
        # print progress bar, incremented by 5%
        progress = int(100 - (remaining / total) * 100)
        progress_tick = progress // tickval
        progress_bar = '[{0}{1}] {2}%'.format('=' * progress_tick,
                                              ' ' * (ticknum - progress_tick),
                                              progress)

    print('\r', progress_bar, end='')
    time.sleep(0.5)

    return completed

def flatten_list(list_to_flatten):
    """ Flattens one level of a nested list.
        Parameters
        ----------
        list_to_flatten : list
            List with nested lists.
        Returns
        -------
        flattened_list : str
            Input list flattened by one level.
    """

    flattened_list = list(itertools.chain(*list_to_flatten))

    return flattened_list

def function_helper(input_args):
    """ Runs function by passing set of provided inputs and
        captures exceptions raised during function execution.
        Parameters
        ----------
        input_args : list
            List with function inputs and function object to call
            in the last index.
        Returns
        -------
        results : list
            List with the results returned by the function.
            If an exception is raised it returns a list with
            the name of the function and the exception traceback.
    """

    try:
        results = input_args[-1](*input_args[0:-1])
    except Exception as e:
        func_name = (input_args[-1]).__name__
        traceback_lines = traceback.format_exception(etype=type(e), value=e,
                                                     tb=e.__traceback__)
        traceback_text = ''.join(traceback_lines)
        print('Error on {0}:\n{1}\n'.format(func_name, traceback_text))
        results = [func_name, traceback_text]

    return results

def write_lines(lines, output_file, joiner='\n'):
    """ Writes a list of strings to a file. The strings
        are joined with newlines before being written to
        file.
        Parameters
        ----------
        lines : list
            List with the lines/strings to write to the
            output file.
        output_file : str
            Path to the output file.
    """

    joined_lines = join_list(lines, joiner)

    write_to_file(joined_lines, output_file, 'a', '\n')

def join_paths(parent_path, child_paths):
    """ Creates a new path by joining a parent directory
        and a list with child paths."""

    joined_paths = os.path.join(parent_path, *child_paths)

    return joined_paths

def create_short(schema_files, schema_dir):
    """ Creates the 'short' directory for a schema.
        Creates the directory and copies schema files
        to the directory (should be used when the schema
        only has 1 sequence per gene/locus).
        Parameters
        ----------
        schema_files : list
            List with paths to all FASTA files in the schema.
        schema_dir : str
            Path to the schema's directory.
        Returns
        -------
        True on completion.
    """

    short_path = join_paths(schema_dir, ['short'])
    create_directory(short_path)

    for file in schema_files:
        short_file = join_paths(short_path, [file_basename(file)])
        short_file = short_file.replace('.fasta', '_short.fasta')
        shutil.copy(file, short_file)

    return True

def sort_data(data, sort_key=None, reverse=False):
    """ Sorts an iterable.
        Parameters
        ----------
        data : iter
            Iterable to sort.
        sort_key
            If provided, data will be sorted based
            on this function.
        reverse : bool
            If sorting order should be inverted.
        Returns
        -------
        sorted_data
            List with sorted elements.
    """

    if sort_key is None:
        sorted_data = sorted(data, reverse=reverse)
    elif sort_key is not None:
        sorted_data = sorted(data, key=sort_key, reverse=reverse)

    return sorted_data

def file_basename(file_path, suffix=True):
    """ Extract file basename from path.
        Parameters
        ----------
        file_path : str
            Path to the file.
        suffix : bool
            Specify if the basename should include the file
            extension.
        Returns
        -------
        basename : str
            File basename extracted from input path.
    """

    basename = os.path.basename(file_path)

    if suffix is False:
        basename = basename.split('.')[0]

    return basename


##############
# clustering##
##############

def cluster_sequences(sequences, word_size, window_size, clustering_sim,
                      representatives, grow_clusters, kmer_offset,
                      seq_num_cluster, temp_directory, cpu_cores,
                      file_prefix, divide, position):
    """ Clusters sequences based on the proportion of shared minimizers.
        Parameters
        ----------
        sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values.
        word_size : int
            Value k for the k-mer size.
        window_size : int
            Value for the window size/number of consecutive
            k-mers per window.
        clustering_sim : float
            Similarity threshold to add a sequence to
            a cluster.
        representatives : dict
            Dictionary with k-mers as keys and a list with
            identifiers of sequences that contain that k-mer
            as values.
        grow_clusters : bool
            If it is allowed to create new clusters.
        kmer_offset : int
            Value to indicate offset of consecutive kmers.
        seq_num_cluster : int
            Maximum number of clusters that a sequence can be
            added to.
        temp_directory : str
            Path to the directory where the clustering results
            will be saved to.
        cpu_cores : int
            Number of clustering processes to run in parallel.
        file_prefix : str
            A prefix to include in the names of created files.
        divide : bool
            If input sequences should be divided into smaller
            groups that can be processed in parallel.
        position : bool
            True if the start position for each k-mer should be saved.
            False otherwise.
        Returns
        -------
        clusters : dict
            Dictionary with representative sequence identifiers
            as keys and lists of tuples as values. Each tuple
            contains a sequence identifier of a clustered
            sequence, the proportion of shared kmers with the
            representative and the length of the clustered
            sequence.
    """

    # sort sequences by length
    sorted_seqs = {k: v for k, v in sort_data(sequences.items(),
                                              sort_key=lambda x: len(x[1]),
                                              reverse=True)}

    if divide is True:
        # divide sequences into sublists
        # do not divide based on number of available cores as it may
        # lead to different results with different number of cores
        cluster_inputs = split_iterable(sorted_seqs,
                                        int(len(sorted_seqs) / 40 + 10))
    else:
        cluster_inputs = [sorted_seqs]

    common_args = [word_size, window_size, clustering_sim,
                   representatives, grow_clusters, kmer_offset,
                   position, seq_num_cluster,
                   clusterer]
    cluster_inputs = [[c, *common_args] for c in cluster_inputs]

    # cluster proteins in parallel
    print('\nClustering sequences based on the proportion '
          'of shared distinct minimizers...')
    clustering_results = map_async_parallelizer(cluster_inputs,
                                                function_helper,
                                                cpu_cores)

    # merge clusters
    clusters = [d[0] for d in clustering_results]
    clusters = merge_dictionaries(clusters)
    rep_sequences = [d[1] for d in clustering_results]
    rep_sequences = merge_dictionaries(rep_sequences)

    # perform clustering with representatives
    if len(cluster_inputs) > 1:
        # cluster representatives
        rep_clusters = clusterer(rep_sequences, word_size,
                                 window_size, clustering_sim,
                                 representatives, grow_clusters,
                                 kmer_offset, position,
                                 seq_num_cluster)

        merged_clusters = {}
        for k, v in rep_clusters[0].items():
            # merge clusters whose representatives are similar
            for n in v:
                # representatives from other clusters are added with
                # similarity score against new representative
                # clustered sequences from other clusters are added
                # with similarity score against their representative
                add_seqids = [n] + [s for s in clusters[n[0]] if s[0] != n[0]]
                merged_clusters.setdefault(k, []).extend(add_seqids)

        clusters = merged_clusters

    print(' \n Clustered {0} sequences into {1} '
          'clusters.'.format(len(sorted_seqs), len(clusters)))

    # sort clusters
    clusters = {k: v for k, v in sort_data(clusters.items())}

    # write file with clustering results
    clusters_out = os.path.join(temp_directory,
                                '{0}.txt'.format(file_prefix))
    write_clusters(clusters, clusters_out)

    return clusters

def write_clusters(clusters, outfile):
    """ Writes information about clusters to file.
        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the decimal proportion of shared
            distinct kmers/minimizers and the length of the
            clustered sequence.
        outfile : str
            Path to the file that will be created to save
            information about clusters.
    """

    cluster_lines = []
    for rep, seqids in clusters.items():
        current_cluster = []
        current_cluster.append('>{0}'.format(rep))
        clustered = [', '.join(['{}'] * len(s)).format(*s)
                     for s in seqids]
        current_cluster.extend(clustered)
        cluster_lines.append(current_cluster)

    # sort by number of lines to get clusters with more sequences first
    cluster_lines = sort_data(cluster_lines,
                              sort_key=lambda x: len(x), reverse=True)
    cluster_lines = flatten_list(cluster_lines)
    cluster_text = join_list(cluster_lines, '\n')

    write_to_file(cluster_text, outfile, 'w', '\n')

def clusterer(sorted_sequences, word_size, window_size,
              clustering_sim, representatives, grow,
              offset, position, seq_num_cluster):
    """ Cluster sequences based on the decimal proportion of
        shared distinct minimizers.
        Parameters
        ----------
        sorted_sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values. Sorted by decreasing sequence
            length.
        word_size : int
            Value k for the kmer size.
        window_size : int
            Window size used to determine minimizers.
        clustering_sim : float
            Similarity threshold to cluster a sequence into
            a cluster.
        representatives : dict
            Dictionary with kmers as keys and a list with
            identifiers of sequences that contain that kmer
            as values.
        grow : bool
            If it is allowed to create new clusters.
        offset : int
            Value to indicate offset of consecutive kmers.
        position : bool
            If minimizer sequence position should be stored.
        seq_num_cluster : int
            Maximum number of clusters that a sequence can be
            added to.
        Returns
        -------
        A list with the following elements:
            clusters : dict
                Dictionary with the identifiers of sequences
                that are cluster representatives as keys and
                a list with tuples as values. Each tuple has
                the identifier of a sequence that was added to
                the cluster, the decimal proportion of shared
                distinct minimizers and the length of the clustered
                sequence.
            reps_sequences : dict
                Dictionary with the identifiers of sequences
                that are cluster representatives as keys and
                their sequences as values.
    """

    clusters = {}
    reps_sequences = {}
    if representatives is None:
        reps_groups = {}
    else:
        reps_groups = representatives

        clusters_ids = set(list(itertools.chain(*(list(reps_groups.values())))))
        clusters = {rep: [] for rep in clusters_ids}

    cluster_results = minimizer_clustering(sorted_sequences, word_size,
                                           window_size, position,
                                           offset, clusters,
                                           reps_sequences, reps_groups,
                                           seq_num_cluster, clustering_sim)
    return cluster_results[0:2]

def select_representatives(kmers, reps_groups, clustering_sim):
    """ Determines the set of clusters that a sequence
        can be added to based on the decimal proportion
        of shared distinct kmers.
        Parameters
        ----------
        kmers : list or set
            Set of kmers determined by decomposing a single
            sequence.
        reps_groups : dict
            Dictionary with kmers as keys and sequence
            identifiers of sequences that contain that
            kmer as values.
        clustering_sim : float
            Sequences are added to clusters if they
            share a minimum decimal proportion of
            distinct kmers with a cluster representative.
        Returns
        -------
        selected_reps : list
            List with a tuple per cluster/representative
            that the sequence can be added to. Each tuple
            has the identifier of the cluster representative
            and the decimal proportion of shared distinct
            kmers.
    """

    current_reps = [reps_groups[k] for k in kmers if k in reps_groups]
    current_reps = flatten_list(current_reps)

    # count number of kmer hits per representative
    counts = Counter(current_reps)
    selected_reps = [(k, v / len(kmers))
                     for k, v in counts.items()
                     if v / len(kmers) >= clustering_sim]

    # sort by identifier and then by similarity to always get same order
    selected_reps = sorted(selected_reps, key=lambda x: x[0])
    selected_reps = sorted(selected_reps, key=lambda x: x[1], reverse=True)

    return selected_reps

def minimizer_clustering(sorted_sequences, word_size, window_size, position,
                         offset, clusters, reps_sequences, reps_groups,
                         seq_num_cluster, clustering_sim):
    """ Cluster sequences based on the decimal proportion of
        shared distinct minimizers.
        Parameters
        ----------
        sorted_sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values. Sorted by decreasing sequence
            length.
        word_size : int
            Value k for the kmer size.
        window_size : int
            Window size used to determine minimizers.
        position : bool
            If minimizer sequence position should be stored.
        offset : int
            Value to indicate offset of consecutive kmers.
        clusters : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the decimal proportion of shared
            distinct minimizers and the length of the clustered
            sequence. This dictionary should be empty at the
            start of the clustering step during the CreateSchema
            process. For the AlleleCall process, the dictionary
            should contain the identifiers of the loci
            representatives.
        reps_sequences : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            their sequences as values. This dictionary should
            be empty at the start of the clustering step
            during the CreateSchema process. For the AlleleCall
            process, the dictionary should contain the
            identifiers of the loci representatives as keys
            and their sequences as values.
        reps_groups : dict
            Dictionary with kmers as keys and a list with
            identifiers of sequences that contain that kmer
            as values. This dictionary should be empty at the
            start of the clustering step during the CreateSchema
            process. For the AlleleCall process, the dictionary
            should contain all kmers contained in the set of
            the schema's representatives sequences.
        seq_num_cluster : int
            Maximum number of clusters that a sequence can be
            added to.
        clustering_sim : float
            Similarity threshold to cluster a sequence into
            a cluster.
        Returns
        -------
        A list with the following elements:
            clusters : dict
                Dictionary with the identifiers of sequences
                that are clusters representatives as keys and
                a list with tuples as values. Each tuple has
                the identifier of a sequence that was added to
                the cluster, the decimal proportion of shared
                distinct minimizers and the length of the clustered
                sequence.
            reps_sequences : dict
                Dictionary with the identifiers of sequences
                that are cluster representatives as keys and
                their sequences as values.
            reps_groups : dict
                Dictionary with kmers as keys and a list with
                identifiers of sequences that contain that kmer
                as values.
    """

    for protid, protein in sorted_sequences.items():

        minimizers = determine_minimizers(protein, window_size,
                                          word_size, offset=offset,
                                          position=position)

        distinct_minimizers = set(minimizers)

        selected_reps = select_representatives(distinct_minimizers,
                                               reps_groups,
                                               clustering_sim)

        top = (len(selected_reps)
               if len(selected_reps) < seq_num_cluster
               else seq_num_cluster)

        # sort to get most similar at index 0
        if len(selected_reps) > 0:
            for i in range(0, top):
                clusters[selected_reps[i][0]].append((protid,
                                                      selected_reps[i][1],
                                                      len(protein),
                                                      len(minimizers),
                                                      len(distinct_minimizers)))
        else:
            for k in distinct_minimizers:
                reps_groups.setdefault(k, []).append(protid)

            clusters[protid] = [(protid, 1.0, len(protein),
                                 len(minimizers), len(distinct_minimizers))]
            reps_sequences[protid] = protein

    return [clusters, reps_sequences, reps_groups]

def sequence_kmerizer(sequence, k_value, offset=1, position=False):
    """ Decomposes a sequence into kmers.
        Parameters
        ----------
        sequence : str
            Sequence to divide into kmers.
        k_value : int
            Value for the size of kmers.
        offset : int
            Value to indicate offset of consecutive kmers.
        position : bool
            If the start position of the kmers in the sequence
            should be stored.
        Returns
        -------
        kmers : list
            List with the kmers determined for the input
            sequence. The list will contain strings if
            it is not specified that positions should be
            stored and tuples of kmer and start position
            if the position is stored.
    """

    if position is False:
        kmers = [sequence[i:i + k_value]
                 for i in range(0, len(sequence) - k_value + 1, offset)]


    elif position is True:
        kmers = [(sequence[i:i + k_value], i)
                 for i in range(0, len(sequence) - k_value + 1, offset)]

    return kmers

def determine_minimizers(sequence, adjacent_kmers, k_value, offset=1,
                         position=False):
    """ Determines minimizers for a sequence based on
        lexicographical order. Skips windows that
        cannot have a minimizer based on the minimizer
        computed in the previous iteration.
        Parameters
        ----------
        sequence : str
            String representing the sequence.
        adjacent_kmers : int
            Window size value. Number of adjacent kmers per group.
        k_value : int
            Value of k for the kmer size.
        offset : int
            Value to indicate offset of consecutive kmers.
        position : bool
            If the start position of the kmers in the sequence
            should be stored.
        Returns
        -------
        minimizers : list
            A list with the set of minimizers determined
            for the input sequence.
    """

    # break sequence into kmers
    kmers = sequence_kmerizer(sequence, k_value,
                              offset=offset, position=position)

    i = 0
    previous = None
    sell = False
    minimizers = []
    # determine total number of windows
    last_window = (len(kmers) - adjacent_kmers)
    while i <= last_window:
        # get kmers in current window
        window = kmers[i:i + adjacent_kmers]
        # pick smallest kmer as minimizer
        minimizer = [min(window)]
        # get position in window of smallest minimizer
        minimizer_idx = window.index(minimizer[0])
        # sliding window that does not included last minimizer
        if previous is None:
            # simply store smallest minimizer
            minimizers.extend(minimizer)

        # sliding window includes last minimizer because we
        # skipped some sliding windows
        else:
            # check if minimizer is the same as the one picked
            # in the last window
            # Do not store minimizer if it is the same
            if minimizer[0] != previous:
                # get kmers smaller than last minimizer
                skipped = window[1:minimizer_idx]
                # determine if any of the smaller kmers is
                # the minimizer of a skipped window
                minimal = previous
                for m in skipped:
                    if m < minimal:
                        minimizer.append(m)
                        minimal = m
                minimizers.extend(minimizer)

        # slide by 1 if minimizer has index 0 in window
        if minimizer_idx == 0:
            i += 1
            previous = None
        # skip sliding windows based on minimizer position
        else:
            i += minimizer_idx
            # if adding minimizer index surpasses last window value we
            # might miss one last minimizer because it will fail the condition
            # find a better way to control this condition!
            if i > last_window and sell is False:
                i = last_window
                sell = True
            previous = minimizer[0]

    return minimizers

def merge_dictionaries(dictionaries_list):
    """ Merges several dictionaries into a single dictionary.
        Parameters
        ----------
        dictionaries_list : list
            A list with the dictionaries to merge.
        Returns
        -------
        merged_dicts : dict
            A dictionary resulting from merging
            all input dictionaries.
    """

    merged_dicts = {}
    for d in dictionaries_list:
        merged_dicts = {**merged_dicts, **d}

    return merged_dicts

def cluster_representative_filter(clusters, representative_filter,
                                  output_directory, file_prefix):
    """ Excludes sequences from clusters based on the proportion
        of shared kmers with the representative. After removing
        highly similar sequences, excludes clusters that are
        singletons (only contain the representative).
        Parameters
        ----------
        clusters : dict
            Dictionary with representative sequence identifiers
            as keys and lists of tuples as values. Each tuple
            contains a sequence identifier of a clustered
            sequence, the proportion of shared kmers with the
            representative and the length of the clustered
            sequence.
        representative_filter : float
            Similarity threshold value. Sequences with
            equal or greater similarity value with the cluster's
            representative are excluded from clusters.
        output_directory : str
            Path to the directory where the clustering results
            will be saved to.
        file_prefix : str
            A prefix to include in the names of created files.
        Returns
        -------
        A list with the following elements:
            pruned_clusters : dict
                Clusters without the sequences that were highly
                similar to the cluster's representative and without
                the clusters that were singletons (only contained
                the representative).
            excluded_seqids : list
                List with the sequence identifiers of the sequences
                that were excluded from the clusters.
    """

    # remove sequences that are very similar to representatives
    pruning_results = representative_pruner(clusters,
                                            representative_filter)

    pruned_clusters, excluded_seqids = pruning_results

    # get identifiers of excluded sequences
    # determine set because same seqids could be in several clusters
    excluded_seqids = set([e[0] for e in excluded_seqids])

    print('Removed {0} sequences based on high similarity with '
          'the cluster representative.'.format(len(excluded_seqids)))

    # remove excluded seqids from clusters without high representative
    # similarity
    pruned_clusters = {k: [e for e in v if e[0] not in excluded_seqids]
                       for k, v in pruned_clusters.items()}

    # write file with pruning results
    pruned_out = os.path.join(output_directory,
                              '{0}_clusters.txt'.format(file_prefix))
    write_clusters(pruned_clusters, pruned_out)

    # identify singletons and exclude those clusters
    singletons = select_clusters(pruned_clusters, 0)
    print('Identified and removed {0} singletons.'.format(len(singletons)))

    pruned_clusters = remove_entries(pruned_clusters, singletons)

    # determine number of sequences that still need to be evaluated
    # +1 to include representative
    clustered_sequences = sum([len(v) + 1 for k, v in pruned_clusters.items()])
    print('Remaining sequences after representative and singleton '
          'pruning: {0}'.format(clustered_sequences))

    return [pruned_clusters, excluded_seqids]

def split_iterable(iterable, size):
    """ Splits a dictionary.
        Parameters
        ----------
        iterable : dict
            Dictionary to split.
        size : int
            Size of dictionaries created from the input
            dictionary.
        Returns
        -------
        chunks : list
            List with dictionaries of defined size
            resulting from splitting the input dictionary.
    """

    chunks = []
    it = iter(iterable)
    for i in range(0, len(iterable), size):
        chunks.append({k: iterable[k] for k in itertools.islice(it, size)})

    return chunks

def representative_pruner(clusters, sim_cutoff):
    """ Removes sequences from clusters based on a similarity
        threshold.
        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the percentage of shared
            kmers/minimizers and the length of the clustered
            sequence.
        sim_cutoff : float
            Similarity threshold value. Sequences with
            equal or greater similarity value are excluded
            from clusters.
        Returns
        -------
        A list with the following elements:
            pruned_clusters : dict
                Input dictionary without values for the
                sequences that had similarity values
                equal or greater that defined threshold.
            excluded : list
                List with a list per cluster. Each list
                has the sequences that were excluded from
                a cluster.
    """

    excluded = []
    pruned_clusters = {}
    for rep, seqids in clusters.items():
        pruned_clusters[rep] = [seqid
                                for seqid in seqids
                                if seqid[1] < sim_cutoff]
        excluded.extend([seqid
                         for seqid in seqids
                         if seqid[1] >= sim_cutoff and seqid[0] != rep])

    return [pruned_clusters, excluded]

def remove_entries(dictionary, keys):
    """ Creates new dictionary without entries with
        specified keys.
        Parameters
        ----------
        dictionary : dict
            Input dictionary.
        keys : list
            List of keys for the entries that should
            not be included in the new dictionary.
        Returns
        -------
        new_dict : dict
            Dictionary without entries with keys in
            the input list.
    """

    new_dict = {k: v for k, v in dictionary.items() if k not in keys}

    return new_dict

def select_clusters(clusters, cluster_size):
    """ Determines clusters that contain a specified number
        of sequences.
        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the percentage of shared
            kmers/minimizers and the length of the clustered
            sequence.
        cluster_size : int
            Number of sequences in the clusters that will
            be selected.
        Returns
        -------
        clusters_ids : list
            List with cluster identifiers for clusters that
            contain a number of sequences equal to specified
            cluster size.
    """

    clusters_ids = [k for k, v in clusters.items() if len(v) == cluster_size]

    return clusters_ids

def cluster_intra_filter(clusters, sequences, word_size,
                         intra_filter, output_directory,
                         file_prefix):
    """ Determines similarity between clustered sequences and
        excludes sequences that are highly similar to other clustered
        sequences.
        Parameters
        ----------
        clusters : dict
            Dictionary with representative sequence identifiers
            as keys and lists of tuples as values. Each tuple
            contains a sequence identifier of a clustered
            sequence, the proportion of shared kmers with the
            representative and the length of the clustered
            sequence.
        sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values.
        word_size : int
            Value k for the k-mer size.
        intra_filter : float
            Similarity threshold value. Sequences with
            equal or greater similarity value with other
            clustered sequences are excluded from clusters.
        output_directory : str
            Path to the directory where the clustering results
            will be saved to.
        file_prefix : str
            A prefix to include in the names of created files.
        Returns
        -------
        A list with the following elements:
            clusters : dict
                Clusters without the sequences that were highly
                similar to other clustered sequences.
            intra_excluded : list
                List with the identifiers of the sequences that
                were excluded.
    """

    # identify clusters with more than 1 sequence
    intra_clusters = {k: v for k, v in clusters.items() if len(v) > 1}

    excluded_seqids, excluded_sims = intra_cluster_sim(intra_clusters,
                                                       sequences,
                                                       word_size,
                                                       intra_filter)

    intra_excluded = [v for k, v in excluded_seqids.items()]
    intra_excluded = flatten_list(intra_excluded)
    # get identifiers of excluded sequences
    # determine set because same seqids could be in several clusters
    intra_excluded = set(intra_excluded)
    print('Removed {0} sequences based on high similarity with '
          'other clustered sequences.'.format(len(intra_excluded)))

    # remove excluded seqids from clusters without high intra-similarity
    pruned_clusters = {k: [e for e in v if e[0] not in intra_excluded]
                       for k, v in clusters.items()}

    # write excluded to file
    intrasim_out = os.path.join(output_directory,
                                '{0}_excluded.txt'.format(file_prefix))
    write_clusters(excluded_sims, intrasim_out)
    # write clusters to file
    intrasim_out = os.path.join(output_directory,
                                '{0}_clusters.txt'.format(file_prefix))
    write_clusters(pruned_clusters, intrasim_out)

    # add key because it is representative identifier
    clustered_sequences = sum([len(v) + 1 for k, v in pruned_clusters.items()])
    print('Remaining sequences after intra-cluster pruning: '
          '{0}'.format(clustered_sequences))

    return [pruned_clusters, intra_excluded]

def intra_cluster_sim(clusters, sequences, word_size, intra_filter):
    """ Determines the percentage of shared kmers/minimizers
        between sequences in the same cluster and excludes
        sequences that are similar to other sequences in the
        cluster.
        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            a list with tuples as values. Each tuple has
            the identifier of a sequence that was added to
            the cluster, the percentage of shared
            kmers/minimizers and the length of the clustered
            sequence.
        sequences : dict
            Dictionary with sequence identifiers as keys
            and sequences as values.
        word_size : int
            Value k for the kmer size.
        intra_filter : float
            Similarity threshold value. If two sequences in
            the same cluster have a similarity value equal
            or greater to this value, the shorter sequence
            will be excluded from the cluster.
        Returns
        -------
        excluded_seqids : dict
            Dictionary with the identifiers of sequences
            that are cluster representatives as keys and
            a list with the sequence identifiers of sequences
            that were excluded from the cluster as values.
        excluded_sims : dict
            A dictionary with the identifiers of sequences
            that are cluster representatives as keys
            and a list with tuples as values. Each tuple
            contains the sequence identifier and the
            similarity value for a match that led to
            an exclusion.
    """

    excluded_seqids = {}
    excluded_sims = {}
    for representative, clustered in clusters.items():
        # get identifiers of sequences in the cluster
        clustered_ids = [c[0] for c in clustered]
        clustered_seqs = {seqid: sequences[seqid] for seqid in clustered_ids}

        # create kmer index
        kmers_mapping, cluster_kmers = kmer_index(clustered_seqs, word_size)

        excluded = []
        similarities = []
        # for each sequence in the cluster
        for seqid, kmers in cluster_kmers.items():
            if seqid not in excluded:
                query_kmers = kmers

                # select sequences with same kmers
                sims = select_representatives(query_kmers,
                                              kmers_mapping,
                                              intra_filter)

                ###################################################################
                if len(sims) > 1:
                    # exclude current sequence
                    candidates = [s for s in sims if s[0] != seqid]
                    for c in candidates:
                        picked = pick_excluded(seqid, c,
                                               clustered_seqs[seqid],
                                               clustered_seqs[c[0]])
                        similarities.append(picked)
                        excluded.append(picked[0])

        # convert into set first to remove possible duplicates
        excluded_seqids[representative] = list(set(excluded))
        excluded_sims[representative] = similarities

    return [excluded_seqids, excluded_sims]

def pick_excluded(query, match, query_seq, match_seq):
    """ Determines which sequence should be excluded
        based on sequence length.
        Parameters
        ----------
        query : str
            Query sequence identifier.
        match : tup
            Tuple with a sequence identifier and
            the decimal proportion of shared
            distinct kmers with the query sequence.
        query_seq : str
            Query sequence.
        match_seq : str
            Match sequence.
        Returns
        -------
        picked : tup
            Tuple with the sequence identifier of the
            sequence that should be excluded, the sequence
            identifier of the sequence that matched with
            the excluded sequence and the decimal proportion
            of shared distinct kmers between both sequences.
    """

    # if query sequence is longer, keep it and exclude candidate
    if len(query_seq) >= len(match_seq):
        picked = (match[0], query, match[1])
    # otherwise, exclude query sequence
    elif len(match_seq) > len(query_seq):
        picked = (query, match[0], match[1])

    return picked

def kmer_index(sequences, word_size):
    """ Creates a kmer index based on a set
        of sequences.
        Parameters
        ----------
        sequences : dict
            Dictionary with sequence identifiers
            as keys and sequences as values.
        word_size : int
            Value k for the kmer size.
        Returns
        -------
        kmers_mapping : dict
            Dictionary with kmers as keys and the
            list of sequence identifiers of the
            sequences that contain the kmers as
            values.
        seqs_kmers : dict
            Dictionary with sequence identifiers
            as keys and the set of distinct kmers
            for each sequence as values.
    """

    kmers_mapping = {}
    seqs_kmers = {}
    for seqid, seq in sequences.items():
        minimizers = determine_minimizers(seq, word_size,
                                          word_size, position=False)
        kmers = set(minimizers)

        # dict with sequence indentifiers and kmers
        seqs_kmers[seqid] = kmers

        # create dict with kmers as keys and list
        # of sequences with given kmers as values
        for kmer in kmers:
            kmers_mapping.setdefault(kmer, []).append(seqid)

    return [kmers_mapping, seqs_kmers]

############
# blasting #
############

def blast_clusters(clusters, sequences, output_directory,
                   blastp_path, makeblastdb_path, cpu_cores,
                   file_prefix):
    """ Uses BLAST to align sequences in the same clusters.
        Parameters
        ----------
        clusters : dict
            Dictionary with representative sequence identifiers
            as keys and lists of tuples as values. Each tuple
            contains a sequence identifier of a clustered
            sequence, the proportion of shared kmers with the
            representative and the length of the clustered
            sequence.
        sequences : dict
            Dictionary with sequence identifiers as keys and
            sequences as values.
        output_directory : str
            Path to the directory where the clustering results
            will be saved to.
        blastp_path : str
            Path to the `BLASTp` executable.
        makeblastdb_path : str
            Path to the `makeblastdb` executable.
        cpu_cores : int
            Number of BLASTp processes to run in parallel.
        file_prefix : str
            A prefix to include in the names of created files.
        Returns
        -------
        A list with the following elements:
            blast_results : list
                List with paths to the files with BLASTp results
                (one file per cluster).
            ids_dict : dict
                Dictionary that maps sequence identifiers to
                shorter and unique integer identifiers used
                to avoid errors during BLAST execution related with
                sequence headers/identifiers that exceed the
                length limit allowed by BLAST.
    """

    print('Clusters to BLAST: {0}'.format(len(clusters)))
    # create FASTA file with sequences in clusters
    clustered_seqs_file = join_paths(output_directory,
                                     ['{0}_clustered_proteins.fasta'.format(file_prefix)])
    clustered_sequences = [[k] + [e[0] for e in v] for k, v in clusters.items()]
    clustered_sequences = flatten_list(clustered_sequences)
    # do not include duplicate identifiers
    clustered_sequences = list(set(clustered_sequences))
    get_sequences_by_id(sequences, clustered_sequences, clustered_seqs_file)

    # create FASTA file with replaced headers to avoid header
    # length limitation in BLAST
    integer_clusters = join_paths(output_directory,
                                  ['{0}_clustered_proteins_int.fasta'.format(file_prefix)])
    ids_dict = integer_headers(clustered_seqs_file, integer_clusters)

    # create BLAST DB
    blast_db = join_paths(output_directory,
                          ['{0}_clustered_proteins_int'.format(file_prefix)])
    db_stderr = make_blast_db(makeblastdb_path, integer_clusters,
                              blast_db, 'prot')
    if len(db_stderr) > 0:
        sys.exit(db_stderr)

    blast_results_dir = os.path.join(output_directory,
                                     '{0}_results'.format(file_prefix))
    os.mkdir(blast_results_dir)

    ######################################################################################################
    # create files with replaced sequence identifiers per cluster
    seqids_to_blast = blast_inputs(clusters, blast_results_dir, ids_dict)

    # distribute clusters per available cores
    process_num = 20 if cpu_cores <= 20 else cpu_cores
    splitted_seqids = split_genes_by_core(seqids_to_blast,
                                          process_num,
                                          'seqcount')

    common_args = [integer_clusters, blast_results_dir, blastp_path,
                   blast_db, cluster_blaster]

    splitted_seqids = [[s, *common_args] for s in splitted_seqids]

    # create the FASTA files with the protein sequences before BLAST?
    print('BLASTing protein sequences in each cluster...\n')

    # BLAST each sequences in a cluster against every sequence in that cluster
    blast_results = map_async_parallelizer(splitted_seqids,
                                           function_helper,
                                           cpu_cores,
                                           show_progress=True)

    return [blast_results, ids_dict]

def cluster_fasta(temp_directory, clusters, AA):
    ''' Creates a fasta file of the clusters in the TMP directory

    parameters
    ----------
    temp_directory : str
        path to the tmp directory
    clusters : dict
        Clusters without the sequences that were highly
        similar to other clustered sequences.
    AA: dict
        seqID coupled to their AA

    Return:
    path to created fasta file:
    '''

    # create output file
    fasta_path = os.path.join(temp_directory, 'cluster_fasta')

    # retrieve identifier and sequences and write down in file:
    with open(fasta_path, 'w') as f:
        for k, v in AA.items():
            if k in clusters.keys():
                f.write((('>') + k + "\n"))
                f.write((v + '\n'))

    return (fasta_path)

def get_sequences_by_id(sequences, seqids, out_file, limit=5000):
    """ Retrieves sequences from an indexed FASTA file.
        Parameters
        ----------
        sequences : dict or Bio.File._IndexedSeqFileDict
            Dictionary with seqids as keys and sequences
            as values or a Fasta file index created with
            BioPython.
        seqids : list
            List with the identifiers of the sequences
            that should be retrieved.
        out_file : str
            Path to the FASTA file to which selected
            sequences will be saved.
        limit : int
            Maximum number of sequences that will be
            kept in memory at a time (to avoid keeping
            huge datasets in memory).
        Returns
        -------
        Creates a file with the sequences that have the
        identifiers in the input list.
    """

    if type(sequences) == dict:
        seqs = [(seqid, sequences[seqid]) for seqid in seqids]
    else:
        seqs = [(seqid, str(sequences[seqid].seq)) for seqid in seqids]

    records = []
    for seq in seqs:
        record = fasta_str_record(seq[0], seq[1])
        records.append(record)

        if len(records) == limit or seq[0] == seqids[-1]:
            lines = join_list(records, '\n')
            write_to_file(lines, out_file, 'a', '\n')
            records = []

def integer_headers(input_fasta, output_fasta, start=1, limit=5000):
    """ Switches FASTA records headers in a file by integer
        values.
        Parameters
        ----------
        input_fasta : str
            Path to the a FASTA file.
        output_fasta : str
            Path to the output file with modified headers.
        start : int
            Integer value of first identifier.
        limit : int
            Maximum number of FASTA records to keep in
            memory.
        Returns
        -------
        ids_map : dict
            Dictionary with mapping between integer and original
            headers.
    """

    seqs = []
    ids_map = {}
    exausted = False
    seq_generator = SeqIO.parse(input_fasta, 'fasta')
    while exausted is False:
        record = next(seq_generator, None)
        if record is not None:
            new_id = 'seq_{0}'.format(start)
            ids_map[new_id] = record.id
            sequence = str(record.seq)
            new_rec = '>{0}\n{1}'.format(new_id, sequence)
            seqs.append(new_rec)
            start += 1
        elif record is None:
            exausted = True

        if len(seqs) == limit or exausted is True:
            write_lines(seqs, output_fasta)
            seqs = []

    return ids_map

def make_blast_db(makeblastdb_path, input_fasta, output_path, db_type,
                  ignore=None):
    """ Creates a BLAST database.
        Parameters
        ----------
        makeblastdb_path : str
            Path to the 'maskeblastdb' executable.
        input_fasta : str
            Path to the FASTA file that contains the sequences
            that should be added to the BLAST database.
        output_path : str
            Path to the directory where the database files
            will be created. Database files will have names
            with the path's basemane.
        db_type : str
            Type of the database, nucleotide (nuc) or
            protein (prot).
        ignore : list of None
            List with BLAST warnings that should be ignored.
        Returns
        -------
        Creates a BLAST database with the input sequences.
    """

    blastdb_cmd = [makeblastdb_path, '-in', input_fasta, '-out', output_path,
                   '-parse_seqids', '-dbtype', db_type]

    makedb_cmd = subprocess.Popen(blastdb_cmd,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    stderr = makedb_cmd.stderr.readlines()

    if len(stderr) > 0:
        stderr = decode_str(stderr, 'utf8')
        if ignore is not None:
            stderr = filter_list(stderr, ignore)

    return stderr

def decode_str(str_list, encoding):
    """ Decodes bytes objects in the input list and
        strips decoded strings from whitespaces and
        newlines.
        Parameters
        ----------
        str_list
            List with string or bytes objects to decode
            and strip of whitespaces and newlines.
        encoding : str
            Encoding codec to use.
        Returns
        -------
        decoded : list
            List with strings without whitespaces or
            newlines.
    """

    decoded = [m.decode(encoding).strip()
               if type(m) == bytes
               else m.strip()
               for m in str_list]

    return decoded

def filter_list(lst, remove):
    """ Removes elements from a list.
        Parameters
        ----------
        lst : list
            Input list.
        remove : list
            List of elements to remove from input list.
        Returns
        -------
        filtered_list : list
            List without the removed elements.
    """

    filtered_list = list(set(lst) - set(remove))

    return filtered_list

def cluster_blaster(seqids, sequences, output_directory,
                    blast_path, blastdb_path):
    """ Aligns sequences in the same cluster with BLAST.
        Parameters
        ----------
        seqids : list
            List with cluster identifiers.
        sequences : str
            Path to the FASTA file with the protein sequences
            in all clusters.
        output_directory : str
            Path to the directory where FASTA files with
            the sequences in each cluster and files with
            BLAST results will be written to.
        blast_path : str
            Path to BLAST executables
        blastdb_path : str
            Path to a BLAST database.
        Returns
        -------
        out_files : list
            List with the paths to the files with the BLAST
            results for each cluster.
    """

    indexed_fasta = SeqIO.index(sequences, 'fasta')

    out_files = []
    for cluster in seqids:

        cluster_id = cluster
        ids_file = os.path.join(output_directory,
                                '{0}_ids.txt'.format(cluster_id))

        with open(ids_file, 'r') as clstr:
            cluster_ids = [l.strip() for l in clstr.readlines()]

        fasta_file = os.path.join(output_directory,
                                  '{0}_protein.fasta'.format(cluster_id))
        # create file with protein sequences
        get_sequences_by_id(indexed_fasta, cluster_ids, fasta_file)

        blast_output = os.path.join(output_directory,
                                    '{0}_blast_out.tsv'.format(cluster_id))

        # BLAST warnings to be ignored
        IGNORE_RAISED = ['Warning: [blastp] To obtain better run time performance, please run '
                         'blastdb_aliastool -seqid_file_in <INPUT_FILE_NAME> -seqid_file_out '
                         '<OUT_FILE_NAME> and use <OUT_FILE_NAME> as the argument to -seqidlist']

        # Use subprocess to capture errors and warnings
        stderr = run_blast(blast_path, blastdb_path, fasta_file,
                           blast_output, 1, 1, ids_file,
                           ignore=IGNORE_RAISED)

        if len(stderr) > 0:
            raise ValueError('\n'.join(stderr))

        out_files.append(blast_output)

    return out_files

def fasta_str_record(seqid, sequence):
    """ Creates the string representation of a FASTA record.
        Parameters
        ----------
        seqid : str
            Sequence identifier to include in the header.
        sequence : str
            String representing DNA or Protein sequence.
        Returns
        -------
        record : str
            String representation of the FASTA record.
    """

    record = '>{0}\n{1}'.format(seqid, sequence)

    return record

def run_blast(blast_path, blast_db, fasta_file, blast_output,
              max_hsps=1, threads=1, ids_file=None, blast_task=None,
              max_targets=None, ignore=None):
    """ Execute BLAST to align sequences in a FASTA file
        against a BLAST database.
        Parameters
        ----------
        blast_path : str
            Path to BLAST executables.
        blast_db : str
            Path to the BLAST database.
        fasta_file : str
            Path to the FASTA file with sequences to
            align against the BLAST database.
        blast_output : str
            Path to the file that will be created to
            store BLAST results.
        max_hsps : int
            Maximum number of High Scoring Pairs per
            pair of aligned sequences.
        threads : int
            Number of threads/cores used to run BLAST.
        ids_file : str
            Path to a file with sequence identifiers,
            one per line. Sequences will only be aligned
            to the sequences in the BLAST database that
            have any of the identifiers in this file.
        blast_task : str
            Type of BLAST task.
        max_targets : int
            Maximum number of target/subject sequences
            to align against.
        ignore : list or None
            List with BLAST warnings that should be ignored.
        Returns
        -------
        stderr : str
            String with errors raised during BLAST execution.
    """

    blast_args = [blast_path, '-db', blast_db, '-query', fasta_file,
                  '-out', blast_output, '-outfmt', '6 qseqid sseqid score',
                  '-max_hsps', str(max_hsps), '-num_threads', str(threads),
                  '-evalue', '0.001']

    if ids_file is not None:
        blast_args.extend(['-seqidlist', ids_file])
    if blast_task is not None:
        blast_args.extend(['-task', blast_task])
    if max_targets is not None:
        blast_args.extend(['-max_target_seqs', str(max_targets)])

    blast_proc = subprocess.Popen(blast_args,
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE)

    stderr = blast_proc.stderr.readlines()

    if len(stderr) > 0:
        stderr = decode_str(stderr, 'utf8')
        if ignore is not None:
            stderr = filter_list(stderr, ignore)

    return stderr

def join_list(lst, link):
    """ Joins all elements in a list into a single string.
        Parameters
        ----------
        lst : list
            List with elements to be joined.
        link : str
            Character used to join list elements.
        Returns
        -------
        joined_list : str
            A single string with all elements in the input
            list joined by the character chosen as link.
    """

    joined_list = link.join(lst)

    return joined_list

def blast_inputs(clusters, output_directory, ids_dict):
    """ Creates files with the identifiers of the sequences
        in each cluster.
        Parameters
        ----------
        clusters : dict
            Dictionary with the identifiers of cluster
            representatives as keys and a list with tuples
            as values (each tuple has the identifier of a
            sequence that is in the cluster, the decimal
            proportion of shared minimizers and the length
            of that sequence).
        output_directory : str
            Path to the directory where files with identifiers
            will be created.
        ids_dict : dict
            Dictionary that maps sequence identifiers to
            shorter and unique identifiers that will be
            saved in the files and used as sequence
            identifiers during BLAST to avoid errors
            related with sequence headers/identifiers
            that exceed length limit allowed by BLAST.
        Returns
        -------
        ids_to_blast : list
            List with the identifiers of all clusters.
    """

    rev_ids = {v: k for k, v in ids_dict.items()}

    ids_to_blast = []
    for rep in clusters:
        cluster_file = os.path.join(output_directory,
                                    '{0}_ids.txt'.format(rev_ids[rep]))
        cluster_ids = [rev_ids[rep]] + [rev_ids[seqid[0]] for seqid in clusters[rep]]
        cluster_lines = join_list(cluster_ids, '\n')
        write_to_file(cluster_lines, cluster_file, 'w', '')
        ids_to_blast.append((rev_ids[rep], len(cluster_ids)))

    return ids_to_blast

def split_genes_by_core(inputs, cores, method):
    """ Creates balanced lists of loci to distribute per number
        of available cores. Loci lists can be created based
        on the number of sequence per locus (seqcount), the mean
        length of the sequences in each locus or the product of
        both values.
        Parameters
        ----------
        inputs : list
            List with one sublist per locus. Each sublist has
            a locus identifier, the total number of sequences
            and sequence mean legth for that locus.
        cores : int
            The number of loci groups that should be created.
            Based on the number of CPU cores that will be
            used to process the inputs.
        method : str
            "seqcount" to create loci lists based on the total
            number of sequences, "length" to split based
            on mean length of sequences and "seqcount+length" to
            split based on both criteria.
        Returns
        -------
        splitted_ids : list
            List with sublists that contain loci identifiers.
            Sublists are balanced based on the chosen method.
    """

    # initialize list with sublists to store inputs
    splitted_ids = [[] for cpu in range(cores)]
    # initialize list with chosen criterion values
    # for each sublist of inputs
    splitted_values = [0 for cpu in range(cores)]
    i = 0
    for locus in inputs:
        if method == 'seqcount':
            splitted_values[i] += locus[1]
        elif method == 'length':
            splitted_values[i] += locus[2]
        elif method == 'seqcount+length':
            splitted_values[i] += locus[1] * locus[2]
        splitted_ids[i].append(locus[0])
        # at the end of each iteration, choose the sublist
        # with lowest criterion value
        i = splitted_values.index(min(splitted_values))

    return splitted_ids

def get_sequences_by_id(sequences, seqids, out_file, limit=5000):
    """ Retrieves sequences from an indexed FASTA file.
        Parameters
        ----------
        sequences : dict or Bio.File._IndexedSeqFileDict
            Dictionary with seqids as keys and sequences
            as values or a Fasta file index created with
            BioPython.
        seqids : list
            List with the identifiers of the sequences
            that should be retrieved.
        out_file : str
            Path to the FASTA file to which selected
            sequences will be saved.
        limit : int
            Maximum number of sequences that will be
            kept in memory at a time (to avoid keeping
            huge datasets in memory).
        Returns
        -------
        Creates a file with the sequences that have the
        identifiers in the input list.
    """

    if type(sequences) == dict:
        seqs = [(seqid, sequences[seqid]) for seqid in seqids]
    else:
        seqs = [(seqid, str(sequences[seqid].seq)) for seqid in seqids]

    records = []
    for seq in seqs:
        record = fasta_str_record(seq[0], seq[1])
        records.append(record)

        if len(records) == limit or seq[0] == seqids[-1]:
            lines = join_list(records, '\n')
            write_to_file(lines, out_file, 'a', '\n')
            records = []

def apply_bsr(blast_results, CDS, bsr, ids_dict):
    """ Computes the BLAST Score Ratio value for BLAST
        alignments and returns the identifiers of the
        sequences that are similar to sequences of
        equal or greater size.
        Parameters
        ----------
        blast_results : list
            List with the path to a file with BLAST
            results in tabular format.
        fasta_file : str
            Path to a FASTA file that contains the
            sequences that were aligned.
        bsr : float
            The BSR value to use as threshold
        ids_dict : dict
            Dictionary with the mapping between
            sequence identifiers used for BLAST and
            the original sequence identifiers.
        Returns
        -------
        excluded_alleles : list
            List with the identifiers of the sequences
            that were highly similar to other sequences.
    """

    self_scores = {r[0]: r[2] for r in blast_results if r[0] == r[1]}
    # do not include self-scores lines, no need to evaluate those hits
    blast_results = [r for r in blast_results if r[0] != r[1]]

    lengths = {}

    for k in self_scores:
        sequence = str(CDS[(ids_dict[k])])
        lengths[k] = len(sequence)

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



############################################
# Add loci  & add gene & prep_recluster    #
############################################
def add_loci_to_seed(add_loci, output_path, force, codon_nr):
    '''
    adds a loci to the seed. if Force is True. genes will be directly added to cluster direcotry
    if not then to CDS_AA directory to be evaluated by clustering

    parameters
    ----------
    :param add_loci: str
        path to fasta file (1 gene per >header)

    :param output_path: str
        path to where the sequence will be stored

    :param force: Bool
        bool , if True then CDS will be directtly stored in Clustering DIR.
        False will the sequence be stored at CDS_AA and the AA be calculated

    :param codon_nr: int
        codont table number be used on the CDS

    :return:
    -------
    NA
    '''
    # retrieves the genes from a file in fasta format, where each gene is separated by >
    genes_to_add = retrieve_genes_from_add_loci_file(add_loci)

    # when force parameter is active, Genes will be directly added to Loci pool but with Force-key so that it can be deleted later
    if force:
        # retrieves the genes from already processed augusts file, to prefent unnecacery duplication
        genes = retrieve_genes_of_loci(output_path, 'clustering')
        print('Current set of loci contains {0} genes'.format(len(genes)))

        # work out genes who are already in the CDS_AA folder
        genes_to_add = remove_dup_from_dic(genes_to_add, genes)
        # add loci directly to clustering dir with Forced pre_fix

        add_loci_to_directory(genes_to_add, output_path, 'clustering', 'forced', '.fasta', codon_nr)

        print('finisht adding loci to scheme seed')

    else:
        # retrieves the genes from already processed augusts file, to prefent unnecacery duplication
        genes = retrieve_genes_of_loci(output_path, 'CDS_AA')
        print('Current set of loci contains {0} genes'.format(len(genes)))

        # work out genes who are already in the CDS_AA folder
        genes_to_add = remove_dup_from_dic(genes_to_add, genes)

        # add loci directly to clustering dir with Forced pre_fix
        add_loci_to_directory(genes_to_add, output_path,'CDS_AA','custom_loci','fna_CDS', codon_nr)

        print('finisht adding loci to genes folder\n'
              'it is advised to Rerun the clustering step for optimal results\n'
              'This can be done by using the command ###### INSERT COMMAND HERE############')

def retrieve_genes_of_loci(output_path, dir):
    '''
    retrieve genes from CDS_AA dir or clustering dir to prevent duplication

    parameters
    ---------
    :param output_path: str
        path to euTyper directory
    :param dir: str
        Directory of CDS_AA or clustering

    :return:
    list of all genes in either CDS_AA or Clustering
    '''
    #check weither output path already exist
    path = os.path.join(output_path, dir)
    if not os.path.isdir(path):
        print('path {0} does not exist, try running euTyper normal first withouth -add_loci or correct path'.format(path))
        exit()
    # Retrieve all genomes to prevent duplication of loci's
    genomes = []
    files = os.listdir(path)
    for file in files:
        try:
            with open(os.path.join(path,file),'r') as f:
                seq = ''
                for line in f:
                    if ">" not in line:
                        seq = seq + line.rstrip()
            genomes.append(seq)
        except:
            # small check that it isnt a directory
            if not os.path.isdir(os.path.join(path,file)):
                print('could not open file "{0}"'.format(file))

    return(genomes)

def retrieve_genes_from_add_loci_file(add_loci):
    '''
    retrieve the genes from a fasta file and checks if it is a DNA sequence (to prevent bugs)

    parameters
    ----------
    :param add_loci: str
        path to fasta file to be read in
    :return:
        list of genes to be added
    '''
    DNA_seq =['A','a',"C",'c',"G","g","T",'t']
    # Retrieve gene and genes from add loci file in fasta format
    genes_to_add = {}
    with open(add_loci,'r') as f:
        for line in f:
            if ">" in line:
                gene_name = line[1:].rstrip()
                genes_to_add[gene_name] = ''
            elif line != '\n':
                genes_to_add[gene_name] = genes_to_add[gene_name] + line.rstrip()

    # go through genes to check if they are DNA sequences
    # if it is detected that it isnt a DNA sequence, it will be deleted

    genes_to_remove = []
    for key in genes_to_add.keys():
        for char in genes_to_add[key]:
            if char not in DNA_seq:
                print('found character {0} in DNA sequence {1} \nThis is from given gene name {2}'.format(char,genes_to_add[key],key))
                print('\ngene will be deleted from genes to add to loci')
                genes_to_remove.append(key)
                break

    #remove non DNA genes
    for key in genes_to_remove:
        genes_to_add.pop(key)

    return(genes_to_add)

def add_loci_to_directory(genes_to_add, output_path, dir1, pre_fix, suffix1, codon_nr):
    '''
    adds loci to a certain directory (CDS_AA or clustering)
    if it is clustering the prefix of forced will be used so it can be identified as possible faulty

    parameters
    ---------
    :param genes_to_add: list
        list of genes to add to the directory
    :param output_path: str
        path to the general euTyper directory

    :param dir1: str
        is either CDS_AA or clustering

    :param pre_fix: str
        the prefix of a gene to identify that that genes has be incorperated into euTyper with using this function

    :param suffix1: str
        the suffix of a gene to be added
    :param codon_nr: int
        incase the gene will not be directly added to CLustering dir, but to CDS_AA. it is required to have his AA known.
        so that incan be evaluated by clustering and blasting

    :return:
    NA
    '''
    # generate path structure to Clustering dir and create fasta files containing the sequences
    path = os.path.join(output_path, dir1)
    # debugging incase euTyper has not run yet
    if not os.path.isdir(path):
        print('clustering directory does not exist \nRun Eutyper normal first')

    # loop through all sequences and generate a fasta file in the required foled
    for key, value in genes_to_add.items():
        print(key)
        #an if statement for if Forced is not enforced.
        # it breaks the loop when the protein sequence doesnt make any sense
        if dir1 == 'CDS_AA':
            #calculate protein sequence from CDS
            protein = translate_dna(value, codon_nr, 0)

            if (protein[-1]) != 'sense':
                print('{0} doesnt produce a protein which makes sense, will be discarded'.format(key))
                # continue skips to next iteration
                continue

        # open file to write down the CDS
        file_name = os.path.join(path,(pre_fix+'_' + key + suffix1))
        with open(file_name,'w') as f:
            f.write(('>forced_'+key.rstrip()+'\n'))
            f.write((value.rstrip()))

        # if the folder is clustering etc for forced add also sequence in short foleder
        if dir1 == 'clustering':
            file_name = os.path.join(path,'short',('forced_' + key +'_short' + suffix1))
            with open(file_name,'w') as f:
                f.write(('>forced_'+key.rstrip()+'\n'))
                f.write((value.rstrip()))

        elif dir1 == 'CDS_AA':
            # open file to write down the CDS
            file_name = os.path.join(path, (pre_fix + '_' + key + 'fna_protein'))
            with open(file_name, 'w') as f:
                f.write(('>forced_' + key.rstrip() + '\n'))
                f.write((str(protein[0][0]).rstrip()))



    print('{0} fasta files added to clustering directory with forced pre_fixed'.format(len(genes_to_add)))

def remove_dup_from_dic(genes_to_add, genes):
    '''
    removes duplication if the gene to be added is already in either CDS_AA or clustering

    parameters
    ----------
    :param genes_to_add: list
        list of genes to add
    :param genes: list
        list of genes already in clustering or CDS_AA

    :return:
    ---------
    list of genes to be added
    '''

    to_pop = []
    for key, value in genes_to_add.items():
            if value in genes:
                print('{0} is already in CDS_AA'.format(key))
                to_pop.append(key)

    for item in to_pop:
            genes_to_add.pop(item)

    return (genes_to_add)

def prepare_temp_dir(output_path, TMP_directory):


    '''
    renames the old temporary directory to prevent bugs or gaining weird outcomes.
    instead of deleting a used can check files or manually rever to the old one

    parameters
    :param output_path: str
        path to general euTyper directory

    :param TMP_directory: str
        name of directory to be renamed

    :return:
    -------
    NA
    '''
    print('start renaming TMP files to reduce buggs')
    # retrieve number and dir name to generate new TMP dir name
    new_temp_dir = TMP_directory.split('_')
    #gen prefix
    count = 1
    #prevents running in error
    if os.path.exists(os.path.join(output_path,'temp',TMP_directory)):
        #loops thorugh every number aslong as the path already exists
        while os.path.exists(os.path.join(output_path, 'temp', (new_temp_dir[0] + '.' + str(count) + '_' + new_temp_dir[1]))):
            print(os.path.join(output_path, 'temp', (new_temp_dir[0] + '.' + str(count) + '_' + new_temp_dir[1])))
            count += 1
        new_file_name = (os.path.join(output_path, 'temp', (new_temp_dir[0] + '.' + str(count) + '_' + new_temp_dir[1])))

        print('old temp directory {0} will be renamed as {1}'.format(TMP_directory, new_file_name))
        #renaming directory to prevent bugs
        os.rename(os.path.join(output_path,'temp',TMP_directory),new_file_name)

def prepare_recluster(output_path):
    '''
    prepared reclustering by renaming old directory. 
    this is so an user can still use the old directory incase of unwanted results
    
    :parameter
    -------------
    :param output_path: str
        str to temporary directory of euTyper
        
    :return: 
    NA
    '''
    if os.path.exists(os.path.join(output_path,'clustering')):

        count = 1
        while os.path.exists(os.path.join(output_path, ('clustering_' + str(count)))):
            count += 1
        new_file_name = (os.path.join(output_path, ('clustering_' + str(count))))
        print('old clustering directory {0} will be renamed as {1}'.format(os.path.join(output_path,'clustering'), new_file_name))

        os.rename(os.path.join(output_path,'clustering'),new_file_name)
##################
# Main functions #
##################


def parse_arguments():
    '''
    Parse all arguments 
    :return: 
    list of arguments
    '''

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)


    ###########################
    # --add_loci method method#
    ###########################

    # -add_loci input files path
    parser.add_argument('--add_loci', type=str, required=False,
                        dest='add_loci', default='false',
                        help='path to input file (fasta format > per gene)')
    # --force, parameter to decide weither added loci should be taken into consideration with clustering or directly adding to set of Loci
    parser.add_argument('--force', type=bool, required=False,
                        dest='force', default=False,
                        help='')

    # --codon, codon table parameter required when forced is not true.
    parser.add_argument('--codon', type=int, required=False,
                        dest='codon', default=0,
                        help='')

    ############################
    # --add_genome method      #
    ############################

    # -add_loci input files path
    parser.add_argument('--add_gene', type=str, required=False,
                        dest='add_gene', default='false',
                        help='path to directory containing fasta files of genomes')

    # -recluster input files path
    parser.add_argument('--recluster', type=bool, required=False,
                        dest='recluster', default=False,
                        help='lofical argument wether reclustering should happen')
    # Normal run required #
    #######################


    # -i input files path
    parser.add_argument('-i', type=str, required=False,
                        dest='input_files', default= 'false',
                        help='')
    # -o output path
    parser.add_argument('-o', type=str, required=True,
                        dest='output_path',
                        help='outputh path but with --add_loci parameter als path used to find working Dir')

    # -s sample metadata path
    parser.add_argument('-s', type=str, required=False, default= 'false',
                        dest='sample_metadata_input',
                        help='')

    #######################################################################
    # optional ( increase --t is highly adivised) must make Thread checker#
    #######################################################################

    # --bsr default =  0.6
    parser.add_argument('--bsr', type=float, required=False,
                        default=0.6, dest='bsr',
                        help='')

    # --t Threads = 1
    parser.add_argument('--t', type=int, required=False,
                        default=1, dest='threads',
                        help='')

    # --p probabillity of augustus = 0.8
    parser.add_argument('--p', type=float, required=False,
                        default=0.8, dest='prob',
                        help='')

    # --m min Minimal sequention length required = 201
    parser.add_argument('--m', type=int, required=False,
                        default=201, dest='min_seq_length',
                        help='')

    # --word_size = 5
    parser.add_argument('--word_size', type=int, required=False,
                        default=5, dest='word_size',
                        help='')

    # --window_size = 5
    parser.add_argument('--window', type=int, required=False,
                        default=5, dest='window_size',
                        help='')

    # --clustering_sim = 0.2
    parser.add_argument('--clustering', type=float, required=False,
                        default=0.2, dest='clustering_sim',
                        help='')


    # --Kmer_offset = 1
    parser.add_argument('--Kmer_offset', type=int, required=False,
                        default=1, dest='Kmer_offset',
                        help='')


    # --seq_num_cluster = 1
    parser.add_argument('--seq_num_cluster', type=int, required=False,
                        default=1, dest='seq_num_cluster',
                        help='')

    # --file_prefix = 'clustering'
    parser.add_argument('--file_prefix', type=str, required=False,
                        default= 'clustering', dest='file_prefix',
                        help='')

    # -- representative_filter = 0.9
    parser.add_argument('--representative_filter', type=float, required=False,
                        default=0.9, dest='representative_filter',
                        help='')

    # -- intra_filter = 0.9
    parser.add_argument('--intra_filter', type=float, required=False,
                        default=0.9, dest='intra_filter',
                        help='')


    args = parser.parse_args()



    #############################
    # edit input data of command#
    #############################
    if os.path.isdir(args.input_files[-1]):
        #edit input files incase / is forgotten
        if (args.input_files[-1]) != '/':
            args.input_files = args.input_files +'/'


    return [args.input_files, args.output_path, args.sample_metadata_input,
            args.bsr, args.threads, args.prob, args.min_seq_length, args.word_size
            , args.window_size, args.clustering_sim,
            args.Kmer_offset, args.seq_num_cluster, args.file_prefix,
            args.representative_filter, args.intra_filter, args.add_loci, args.force,
            args.codon, args.add_gene, args.recluster]

def check_arguments(args):
    """ Checks the given parameters so far if it is in the right format.
        or if something is required will print it and exit euTyper to prevent missuse
        Parameters
        ----------
        inputs : list of arguments
        Returns: premature exiting programm
        
            
    """
    input_data_files_sample_metadata = []

    availeble_codon_table =[1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25,26,27,28,29,30,30,31,33]



    #######################
    # Checking  requirements#
    #######################

    #############################
    # check the --add_loci modus#
    #############################
    if args[15] != 'false':

        if not (os.path.isfile(args[15])):
                print('add_loci path doesnt lead to a file')
                print('current add_loci path is "{0}" '.format(args[15]))
                exit()

            # this check only should happen when force is False
        if not args[16]:
                # check weither a codon number is given
                if args[17] == 0:
                    print('add a codon number with --codon')
                    exit()
                # check wether codon table number is in list of available numbers
                elif args[17] not in availeble_codon_table:
                    print('given codon table number {0} is not in current list of available codon table numbers\n'
                          'please choose a different number'.format(args[17]))
                    exit()

        #check wether working directory existis

        if not os.path.isdir(args[1]):
            print('path {0} does not exists, run normal euTyper first or correct path to already an existing working dir'.format(args[1]))
            exit()
    ###########################
    # check the add_gene modus#
    ###########################

    elif args[18] != 'false':
        #check wether output path exists
        if  os.path.isdir(args[18]):
            files_in_input = os.listdir(args[18])
        else:
            print("{0} is not a path to directory \nit is required to give a path to a directory containing fasta"
                  "files of genomes to be processed by Augustus".format(args[18]))
            exit()


        #check sample meta data
        if args[2] != 'false':
                print(args[2])
                with open(args[2]) as f:
                    for line in f:
                        line = line.strip().split('\t')

                        # checks if there are 3 columns
                        if len(line) != 3:
                            print("sample meta data doesnt contain 3 column ( tab delimated) ")
                            exit()

                        # test weither samplemetada codon table number contains int
                        try:
                            codon_table_nr = int(line[2])
                            if codon_table_nr not in availeble_codon_table:
                                print( ' codon table nr ' + (line[2]) + ' is not available ')
                                print (' exiting euTyper')
                                exit()
                        except:
                            print(' codon table doesnt all contain numbers')
                            print (' input ' + (line[2]) + ' is not a valid choise')
                            print(' exiting euTyper')
                            exit()

                    # checks weither the input data names are correct in sample metadata
                        input_data_files_sample_metadata.append(line[0])

                    for element in input_data_files_sample_metadata:
                        if element not in files_in_input:
                            print( (element) + ' is in sample metadata is not found in input data')
                            print(' exiting euTyper')
                            exit()
           # except:
           #     print('sample metadata file (path) {0} could not be openend please check path'.format(args[2]))
           #     exit()
        else:
            print('-s parameter is empty , fill in path to sample_metadata \n'
                  'sample metadata is a TSV file containing samplenames species(augustus) and codon table number')
            exit()

        # check output path

        if not os.path.isdir(args[1]):
            print(' -o path {0} does not exists, run normal euTyper first or correct path to already an existing working dir'.format(args[1]))
            exit()


    ##########################
    # check recluster#
    ########################

    elif args[19]:
        # check output path

        if not os.path.isdir(args[1]):
            print(' -o path {0} does not exists, run normal euTyper first or correct path to already an existing working dir'.format(args[1]))
            exit()
    ############################
    # check normal intended run#
    ############################
    else:

        if args[0] != 'false':
            print(args[0])
            if os.path.isdir(args[0]):
                files_in_input = os.listdir(args[0])
            else:
                print(' input data doesnt lead to directory')
                print('exiting euTyper')
                exit()
        else:
            print('argument -i is empty, please fill in an input directory containing genomes of fasta files \n')
            exit()

            # check sample meta data
        if args[2] != 'false':
            try:
                with open(args[2]) as f:
                    for line in f:
                        line = line.strip().split('\t')

                        # checks if there are 3 columns
                        if len(line) != 3:
                            print("sample meta data doesnt contain 3 column ( tab delimated) ")
                            exit()

                        # test weither samplemetada codon table number contains int
                        try:
                            codon_table_nr = int(line[2])
                            if codon_table_nr not in availeble_codon_table:
                                print(' codon table nr ' + (line[2]) + ' is not available ')
                                print(' exiting euTyper')
                                exit()
                        except:
                            print(' codon table doesnt all contain numbers')
                            print(' input ' + (line[2]) + ' is not a valid choise')
                            print(' exiting euTyper')
                            exit()

                        # checks weither the input data names are correct in sample metadata
                        input_data_files_sample_metadata.append(line[0])

                    for element in input_data_files_sample_metadata:
                        if element not in files_in_input:
                            print((element) + ' is in sample metadata is not found in input data')
                            print(' exiting euTyper')
                            exit()
            except:
                print('sample metadata file (path) {0} could not be openend please check path'.format(args[2]))
                exit()

        else:
            print('-s parameter is empty , fill in path to sample_metadata \n'
                  'sample metadata is a TSV file containing samplenames species(augustus) and codon table number')
            exit()
        #check of probabillity cutoff.
        try:
            if float(args[5]) > 1.0 or float(args[5]) < 0.0:
                print('probabillity cutoff is outside range 0.0 - 1.0 \n input {0} was not accepted'.format(args[5]))
                exit()
        except:
            print('probabillity cutoff could not be converted to an float \n input {0} was not accepted'.format(str(args[5])))
            exit()

        # checks wether output directory already existis
        if os.path.isdir(args[1]):
            print('path {0} does already exists, to prevent overwriting please select a new output folder \n'
                  'euTyper will create the new folder'.format(args[1]))


        #checks if the output folder is linked to a file. if not existign or to an other directory that is okay
        if os.path.isfile(args[1]):
            print(' output directory path is linked to a file, please change : ' + (args[1]))
            exit()

        if os.path.isdir(args[1]):
            print(' output directory path is linked to a directory, please change : ' + (args[1]))
            print(' to prevent overwriting other files')
            exit()

def pre_processing(sample_information, output_path, threads, prob, pt_to_input_folder):
        """ The preprocessing runs Augustus to gain the CDS location of the gnomes and assemblies
            From the CDS the AA will be calculated by the given codon table

            Parameters
            ----------
            - path to sample meta data
            - output path
            - sample name
            - species parameter
            - codon table number
            Returns:
            - Augustus file
            - CDS sequences
            - AA sequences
        """

        print('\n Running AUGUSTUS on ' + str(len(sample_information)) + ' sub samples')
        temp_dir_2, temp_dir = execute_augustus(sample_information, threads, output_path, pt_to_input_folder)

        # retrieve the CDS location and merging the parralel augustus output
        tmp_output = retrieve_CDS_AA_from_augustus(temp_dir_2, sample_information, temp_dir, output_path, prob)

        # combining created files
        fasta_location = combine_fastas(tmp_output, sample_information, output_path)

        return (fasta_location)

def seed_creation(output_directory, MinSeqlen, word_size, window_size, clustering_sim,
                      cpu_cores, kmer_offset, seq_num_cluster, representative_filter, intra_filter, BSR, schema_name):
    '''
    Gereral flow to create the Scheme seed of all the genes and AA found in the CDS_AA folder
    :param output_directory: 
    :param MinSeqlen: 
    :param word_size: 
    :param window_size: 
    :param clustering_sim: 
    :param cpu_cores: 
    :param kmer_offset: 
    :param seq_num_cluster: 
    :param representative_filter: 
    :param intra_filter: 
    :param BSR: 
    :param schema_name: 
    :return: 
    '''
    #################################################
    # Retrieve and check protein files from augustus#
    # if they match in the given sample_metadata    #
    #################################################

    # define directory for temporary files
    temp_directory = join_paths(output_directory, ['temp'])
    os.makedirs(temp_directory, exist_ok=True)
    # create scheme seed
    files = []

    for file in os.listdir(os.path.join(output_directory,'CDS_AA')):
        if "_CDS" in file:
            files.append((file.split('_CDS'))[0])
        elif '_protein' in file:
            files.append((file.split('_protein'))[0])

    #remove duplication
    files = list(set(files))
    # find protein and CDA file of samples and couple to dictionary#
    paths = {}
    # retrieve file names of the protein output path
    path = os.path.join(output_directory, 'CDS_AA')
    files1 = sorted(os.listdir(path))

    for file in files:
        TMP_path = []
        for file1 in files1:
            if file in file1:
                TMP_path.append(os.path.join(path, file1))
        if len(TMP_path) == 2:
            paths[file] = TMP_path

    # generate dictionary with the CDS coupled to their AA
    CDS_AA = retrieve_CDS_AA(paths)

    print(('\nRemoving smaller then ' + str(MinSeqlen) + ' nucleotides '), end='')
    CDS_AA = remove_small(CDS_AA, MinSeqlen)

    # remove duplicates
    print('\nRemoving duplicated sequences...', end='')
    CDS_AA = remove_duplicates(CDS_AA)


    #####################################
    # clustering found protein sequences#
    #####################################

    # create directory to store clustering data
    clustering_dir = join_paths(temp_directory, ['4_clustering'])
    os.makedirs(clustering_dir, exist_ok=True)

    AA = {}
    # create dictionary with only AA
    for k, v in CDS_AA.items():
        AA[k] = v[1]

    CDS = {}
    for k, v in CDS_AA.items():
        CDS[k] = v[0]

    print(('amount of sequences at prior clustering ') + str(len(AA)))

    # first clustering step
    cs_results = cluster_sequences(AA, word_size, window_size,
                                   clustering_sim, None, True,
                                   kmer_offset, seq_num_cluster, clustering_dir, cpu_cores,
                                   'clusters', True, False)

    # Refining the clusters
    cp_results = cluster_representative_filter(cs_results,
                                               representative_filter,
                                               clustering_dir,
                                               'repfilter')

    clusters, excluded_seqids = cp_results

    # remove excluded seqids
    schema_seqids = list(set(AA) - excluded_seqids)

    print(('ID after excluding ' + str(len(schema_seqids))))

    # intra cluster pruner step
    cip_results = cluster_intra_filter(clusters, AA,
                                       word_size, intra_filter,
                                       clustering_dir,
                                       'intrafilter')

    clusters, intra_excluded = cip_results

    # remove excluded seqids - we get set of sequences from clusters
    # plus singletons
    schema_seqids = list(set(schema_seqids) - set(intra_excluded))

    print(('ID after intra filter ' + str(len(schema_seqids))))
    print((('amount of cluster sequences ') + str(len(clusters))))

    #########################
    # Blast P found clusters#
    #########################

    BLASTP_ALIAS = 'blastp.exe' if platform.system() == 'Windows' else 'blastp'
    MAKEBLASTDB_ALIAS = 'makeblastdb.exe' if platform.system() == 'Windows' else 'makeblastdb'
    blast_path = ''

    blastp_path = os.path.join(blast_path, BLASTP_ALIAS)
    makeblastdb_path = os.path.join(blast_path, MAKEBLASTDB_ALIAS)

    if len(clusters) > 0:
        blasting_dir = join_paths(clustering_dir, ['cluster_blaster'])
        create_directory(blasting_dir)

        blast_results, ids_dict = blast_clusters(clusters, AA,
                                                 blasting_dir, blastp_path,
                                                 makeblastdb_path, cpu_cores,
                                                 'blast')

        blast_files = flatten_list(blast_results)
        ########################
        # combine all DNA fasta#
        ########################



        # compute and exclude based on BSR
        blast_excluded_alleles = [apply_bsr(read_tabular(file),
                                            CDS,
                                            BSR,
                                            ids_dict)
                                  for file in blast_files]

        # merge bsr results
        blast_excluded_alleles = flatten_list(blast_excluded_alleles)

        blast_excluded_alleles = [ids_dict[seqid] for seqid in blast_excluded_alleles]
        schema_seqids = list(set(schema_seqids) - set(blast_excluded_alleles))
        print('\n\nRemoved {0} sequences based on high BSR value with '
              'other sequences.'.format(len(set(blast_excluded_alleles))))

        # perform final BLAST to identify similar sequences that do not
        # share many/any kmers
        print('Total of {0} sequences to compare in final BLAST.'.format(len(schema_seqids)))

        # sort seqids before final BLASTp to ensure consistent results
        schema_seqids = sort_data(schema_seqids, sort_key=lambda x: x.lower())

        # create directory for final BLASTp
        final_blast_dir = join_paths(temp_directory, ['5_final_blast'])
        create_directory(final_blast_dir)

        beta_file = os.path.join(final_blast_dir, 'pre_schema_seed.fasta')
        get_sequences_by_id(AA, schema_seqids, beta_file)

        integer_seqids = os.path.join(final_blast_dir, 'pre_schema_seed_int.fasta')
        ids_dict2 = integer_headers(beta_file, integer_seqids)

        blast_db = join_paths(final_blast_dir, ['pre_schema_seed_int'])
        db_stderr = make_blast_db(makeblastdb_path, integer_seqids, blast_db, 'prot')

        if len(db_stderr) > 0:
            sys.exit(db_stderr)

        print('Performing final BLASTp...', end='')
        blast_output = '{0}/{1}_blast_out.tsv'.format(final_blast_dir,
                                                      'pre_schema_seed')
        blast_stderr = run_blast(blastp_path, blast_db, integer_seqids,
                                 blast_output, 1, cpu_cores)

        if len(blast_stderr) > 0:
            sys.exit(blast_stderr)

        final_excluded = apply_bsr(read_tabular(blast_output),
                                   CDS,
                                   BSR,
                                   ids_dict2)
        final_excluded = [ids_dict2[seqid] for seqid in final_excluded]

        schema_seqids = list(set(schema_seqids) - set(final_excluded))

        print('removed {0} sequences that were highly similar '
              'to other sequences.'.format(len(final_excluded)))

        output_schema = os.path.join(final_blast_dir, 'schema_seed.fasta')

        # create file with the schema representative sequences
        get_sequences_by_id(CDS, schema_seqids, output_schema)

        schema_files = create_schema_structure(output_schema, output_directory,
                                               final_blast_dir, schema_name)

        return [schema_files, temp_directory]

def create_schema_structure(schema_seed_fasta, output_directory,
                            temp_directory, schema_name):
    """ Creates the schema seed directory with one FASTA file per
        distinct locus and the `short` directory with the FASTA files
        used to save the representative sequences.
        Parameters
        ----------
        schema_seed_fasta : str
            Path to the FASTA file that contains the sequences that
            constitute the schema seed. Each FASTA record in the file
            is a representative sequence chosen for a locus.
        output_directory : str
            Path to the main output directory of the process.
        temp_directory : str
            Path to the directory where the FASTA file with the
            records for the schema seed will be saved to.
        schema_name : str
            Name for the schema's directory.
        Returns
        -------
        schema_files : list
            List with the paths to the FASTA files in the schema seed.
    """
    CHAR_REPLACEMENTS = [("|", "_"), ("_", "-"), ("(", ""),
                         (")", ""), ("'", ""), ("\"", ""), (":", "")]
    # add allele identifier to all sequences
    schema_records = ['>{0}\n{1}'.format(replace_multiple_characters(rec.id, CHAR_REPLACEMENTS) + '_1', str(rec.seq))
                      for rec in SeqIO.parse(schema_seed_fasta, 'fasta')]

    final_records = os.path.join(temp_directory, 'schema_seed_final.fasta')
    write_lines(schema_records, final_records)

    schema_dir = join_paths(output_directory, [schema_name])
    create_directory(schema_dir)

    # create directory and schema files
    filenames = (record.id[:-2] for record in SeqIO.parse(final_records, 'fasta'))
    schema_files = split_fasta(final_records, schema_dir, 1, filenames)
    create_short(schema_files, schema_dir)

    return schema_files



def main(args):

    ###########################################################
    # giving arguments different variables to increase clarity#
    ###########################################################

    #sample_information, output_path,threads, prob ,path_to_input_folder
    path_to_input_folder = args[0]
    output_path = args[1]
    sample_metadata_path = args[2]
    BSR_cutoff = args[3]
    threads = args[4]
    prob = args[5]
    min_seq_length = args[6]
    word_size = args[7]
    window_size = args[8]
    clustering_sim = args[9]
    kmer_offset = args[10]
    seq_num_cluster = args[11]
    file_prefix = args[12]
    representative_filter = args[13]
    intra_filter = args[14]
    add_loci_path = args[15]
    force = args[16]
    codon = args[17]
    add_gene_path = args[18]
    recluster = args[19]




    ###########
    # add_loci#
    ###########
    if add_loci_path != 'false':

        start_date = dt.datetime.now()
        start_date_str = dt.datetime.strftime(start_date, '%Y-%m-%dT%H:%M:%S')
        print('\nStarted at: {0}\n'.format(start_date_str))
        output_path1 = os.path.join(output_path,'clustering')
        print('started adding genes from {0} to scheme seed in working dir {1}'.format(add_loci_path, output_path1))
        add_loci_to_seed(add_loci_path, output_path, force, codon)

        if not force:
            if recluster:
                # prepare for reclustering

                prepare_temp_dir(output_path, '4_clustering')
                prepare_temp_dir(output_path, '5_final_blast')

                # recluster again
                results = seed_creation(output_path, min_seq_length, word_size, window_size,
                                        clustering_sim,
                                        threads, kmer_offset, seq_num_cluster, representative_filter, intra_filter,
                                        BSR_cutoff, file_prefix)


    ############
    # add genes#
    ############
    elif add_gene_path != 'false':


        start_date = dt.datetime.now()
        start_date_str = dt.datetime.strftime(start_date, '%Y-%m-%dT%H:%M:%S')
        print('\nStarted at: {0}\n'.format(start_date_str))
        output_path1 = os.path.join(output_path, 'CDS_AA')
        # prepare adding genes by editing the TMP directories

        prepare_temp_dir(output_path, '1_splitted_fasta')
        prepare_temp_dir(output_path, '2_augustus_output')
        prepare_temp_dir(output_path, '3_splitted_fastas')

        print('started adding genomes from {0} to working dir {1}'.format(add_gene_path, output_path1))

        sample_information = retrieve_sample_data(sample_metadata_path)

        fasta_location = pre_processing(sample_information, output_path,threads, prob ,add_gene_path)
        print(' finisht adding genomes to "{0}"'.format(fasta_location))


        ######## NOTE IMPLEMENT AGAIN CLUSTERING############ if statement is true.

        if recluster:
            prepare_temp_dir(output_path, '4_clustering')
            prepare_temp_dir(output_path, '5_final_blast')
            prepare_recluster(output_path)
            results = seed_creation(output_path, min_seq_length, word_size, window_size,
                                    clustering_sim,
                                    threads, kmer_offset, seq_num_cluster, representative_filter, intra_filter,
                                    BSR_cutoff, file_prefix)




    ############
    # recluster#
    ############

    elif recluster:
        prepare_temp_dir(output_path, '4_clustering')
        prepare_temp_dir(output_path, '5_final_blast')
        prepare_recluster(output_path)
        results = seed_creation(output_path, min_seq_length, word_size, window_size,
                                clustering_sim,
                                threads, kmer_offset, seq_num_cluster, representative_filter, intra_filter,
                                BSR_cutoff, file_prefix)
    #########
    # Total #
    #########
    else:


        # creating a nice statistic overiew of the euTyper run
        print(' sample metadata and input files passed the check')
        start_date = dt.datetime.now()
        start_date_str = dt.datetime.strftime(start_date, '%Y-%m-%dT%H:%M:%S')
        print('\nStarted at: {0}\n'.format(start_date_str))

        print('Number of genomes/assemblies: {0}'.format( len(os.listdir(path_to_input_folder))))
        print('data will be stored at: {0}'.format(args[1]))

        print('BLAST Score Ratio: {0}'.format(args[3]))
        print('Number of threads: {0}'.format(threads))

        # generating output folder if required
        if not os.path.isdir(output_path):
            os.makedirs(output_path)

        #retrieve the sample names from sample metadata
        sample_information = retrieve_sample_data(sample_metadata_path)

        print(sample_information)

        ######################################
        # start preprocessing all input files#
        ######################################
        fasta_location = pre_processing(sample_information, output_path,threads, prob ,path_to_input_folder)

        print(fasta_location)
        ##############################################
        # CREATE scheme seed by clustering & blasting#
        ##############################################
        results = seed_creation(output_path, min_seq_length, word_size, window_size, clustering_sim,
                      threads, kmer_offset, seq_num_cluster, representative_filter, intra_filter, BSR_cutoff, file_prefix)


        # print message about schema that was created
        print('Created schema seed with {0} loci.'.format(len(results[0])))

        print(results[1])

if __name__ == '__main__':

    print('\n{0}\n{1}\n{0}'.format('-' * 9, ' euTyper'))

    # python3 eutyper.py -i /media/DATA/mmib/Student_DATA/mmibstudent1/Donny_data/input_data  -o /media/DATA/mmib/Student_DATA/mmibstudent1/Donny_data/euTyper_output -s /media/DATA/mmib/Student_DATA/mmibstudent1/Donny_data/sample_metadata
    # parse arguments or add the default for programm
    args = parse_arguments()
    # check args to definine if programm should run#
    check_arguments(args)
    main(args)
