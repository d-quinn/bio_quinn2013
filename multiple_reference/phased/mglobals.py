import os
from os.path import join
import logging
import glob
import csv
from collections import namedtuple
import sys
import pprint
import pickle
import multiprocessing

cpu_count = multiprocessing.cpu_count()


# Project working directory, where samples will be linked to
samples_path = '/path/to/working_directory'
# Path to supporting files
supporting_files = '/path/to/Quinn2013_supporting_files'
original_path = join(samples_path, 'original')
alternate_path = join(samples_path, 'alternate')
coverage_cutoffs = join(supporting_files, 'cov_cutoffs.csv')
dros_gtf = join(supporting_files, 'iGenomes/Drosophila_melanogaster/Ensembl/'
                'BDGP5.25/Annotation/Genes/genes.gtf')
dros_gtf_index = join(supporting_files, 'transcriptome_data')
dros_gtf_cds = join(supporting_files, 'CDSgtf/CDS.gtf')
all_snps = join(supporting_files, 'freeze2_sorted_nohead_small.vcf')
# Changed this to vcf so I could use it as the reference vcf in snps_combine.
# Used convert_Frank_vcfs_to_bed.py to convert.
dgrp_5_lines_sup = join(supporting_files, 'freeze2_Filter2homoState.bins__5.vcf')
dgrp_both_11 = join(supporting_files, 'freeze2_362_765_noINDEL_1-1.vcf')

current_snp_file = dgrp_both_11
current_genes_file = dros_gtf

# Software paths
varscan_path = join(supporting_files, 'VarScan.v2.3.5.jar')
gatk_path = join(supporting_files, 'GenomeAnalysisTK-2.4-9/GenomeAnalysisTK.jar')
trimmomatic_path = join(supporting_files, 'Trimmomatic-0.30/trimmomatic-0.30.jar')

# Create logger
log_filename = join(samples_path, 'pipeline_log.log')
format = '%(levelname)s [%(asctime)s] %(message)s'
datefmt = '%m/%d/%Y %H:%M:%S'
log = logging.getLogger('pipeline')
log.setLevel(logging.DEBUG)
fh = logging.FileHandler(log_filename)
fh.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
formatter = logging.Formatter(fmt=format, datefmt=datefmt)
fh.setFormatter(formatter)
ch.setFormatter(formatter)
log.addHandler(fh)
log.addHandler(ch)

log.info('\n\n')
log.info('{:-^50}'.format('START'))

log.info('Creating directories')
if not os.path.exists(original_path):
    os.mkdir(original_path)
if not os.path.exists(alternate_path):
    os.mkdir(alternate_path)

# Link in sample
os.chdir(samples_path)
try:
    os.symlink(join(supporting_files, 'massouras_FASTA', '16_A12_pUn_down.fastq'),
               '16_A12_pUn_down.fastq')
except OSError:
    pass

# Load in the samples list if it exists, or else build it
try:
    with open(join(samples_path, 'samples_list.p'), 'rb') as pf:
        samples_list = pickle.load(pf)
        num_samples = list({f[:-6] for f in os.listdir(samples_path)
                            if f.endswith('.fastq') and
                            not (f.endswith('trim_R1.fastq') or f.endswith('trim_R2.fastq'))})
        assert len(samples_list) == len(num_samples), 'unpickling error!'
    log.info('unpickled samples_list')
except IOError:
    samples_list = list({f[:-6] for f in os.listdir(samples_path)
                        if f.endswith('.fastq') and not
                        (f.endswith('trim_R1.fastq') or f.endswith('trim_R2.fastq'))})
    with open(join(samples_path, 'samples_list.p'), 'wb') as pf:
        pickle.dump(samples_list, pf)
    log.info('built and pickled samples_list')

pretty_list = pprint.pformat(samples_list)
log.debug('Samples list:\n{}'.format(pretty_list))

if samples_list:
    for sample in samples_list:
        if not os.path.exists(join(original_path, sample)):
            os.mkdir(join(original_path, sample))
        if not os.path.exists(join(alternate_path, sample)):
            os.mkdir(join(alternate_path, sample))
else:
    raise AssertionError('Add samples to samples_path and restart')

original_fastas = ['16_A12_pUn_down_freeze2__orig__.fa']

alternate_fastas = ['16_A12_pUn_down_freeze2__new__.fa']


def original_fasta_link():
    '''
    Link to original reference fastas - these are needed for alignment and
    for building the alternate fastas later.
    '''
    sequence_path = join(supporting_files, 'iGenomes/Drosophila_melanogaster/Ensembl/'
                        'BDGP5.25/Sequence')
    bowtie_files = glob.glob(join(sequence_path, 'Bowtie2Index', 'genome.*'))

    for sample in samples_list:
        os.chdir(join(original_path, sample))
        # Link to bowtie index
        for index_file in bowtie_files:
            try:
                os.symlink(index_file, os.path.basename(index_file))
            except OSError:
                pass
        # Link to fasta index and dictionary (convenience to prevent mpileup from building)
        try:
            os.symlink(join(sequence_path, 'WholeGenomeFasta', 'genome.dict'),
                       'genome.dict')
        except OSError:
            pass
        try:
            os.symlink(join(sequence_path, 'WholeGenomeFasta', 'genome.fa.fai'),
                       'genome.fa.fai')
        except OSError:
            pass

original_fasta_link()

log.info('Setup complete!')
