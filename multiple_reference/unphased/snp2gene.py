#!/usr/bin/python

import pysam
import os
import sys
import csv
import operator
import logging
import collections
from itertools import chain

log = logging.getLogger('pipeline')


def converter(input_f, output_f):
    '''
    Converts a vcf that has been annotated using intersectBed with a
    genes.gtf file into a more easily manipulated csv.
    '''
    with open(input_f, 'rb') as f:
        readfields = ('CHROM', 'POS', 'ID', 'REF', 'ALT',
                      'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE')
        initial_reader = csv.DictReader(f, delimiter='\t', fieldnames=readfields, restkey='GTFINFO')
        posSorted = sorted(initial_reader, key=lambda row: int(row['POS']))
        sortedreader = sorted(posSorted, key=lambda row: row['CHROM'])

        row_fields_ofinterest = ('CHROM', 'POS', 'REF', 'ALT')
        sample_fields_ofinterest = ('RD', 'AD')
        gtf_fields_ofinterest = ('gene_id', 'exon_number', 'gene_name')
        fields_of_interest = row_fields_ofinterest + sample_fields_ofinterest + gtf_fields_ofinterest

        def non_duplicated_positions(r):
            position_buffer = None
            for row in r:
                # Check to make sure the vcf has the correct number of columns
                assert len(row['GTFINFO']) == 9, 'VCF contains too many columns'
                position = (row['CHROM'], row['POS'])
                if position != position_buffer:
                    yield row
                    position_buffer = position

        def reorganized_rows(r):
            for row in r:
                format_list = row['FORMAT'].split(':')
                sample_list = row['SAMPLE'].split(':')
                sample_dict = dict(zip(format_list, sample_list))

                gtf_stuff = row['GTFINFO'][8].split(';')[:-1]
                gtf_list = [item.replace('"', '').strip() for item in gtf_stuff]
                gtf_dict = {title: geneinfo for title, geneinfo in (item.split(' ')
                            for item in gtf_list)}
                newrow = {}
                for field in row_fields_ofinterest:
                    newrow[field] = row[field]
                for field in sample_fields_ofinterest:
                    newrow[field] = sample_dict[field]
                for field in gtf_fields_ofinterest:
                    newrow[field] = gtf_dict[field]
                yield newrow

        non_dups = non_duplicated_positions(sortedreader)
        reorganized_rows = reorganized_rows(non_dups)

        with open(output_f, 'wb') as fout:
            writer = csv.DictWriter(fout, fieldnames=fields_of_interest, lineterminator='\n',
                                    delimiter='\t')
            writer.writerows(reorganized_rows)


class Variant(object):

    def __init__(self, original, new, refdict):
        self.original = original
        self.new = new
        self.chrom = self.original['CHROM']
        self.position = int(self.original['POS'])
        self.coordinates = (self.chrom, self.position)
        self.real_reference_base = refdict[self.coordinates].reference
        self.real_alternate_base = refdict[self.coordinates].alternate
        self.gene_id = self.original['gene_id']
        self.exon_number = self.original['exon_number']
        self.gene_name = self.original['gene_name']
        self.cutoffs = cos  # cos is a global variable defined in the combine_snps function

    def __repr__(self):
        return 'Variant(original={}, new={})'.format(self.original, self.new)

    def get_snp_data(self, orig_bam, new_bam):

        # If this is the second pass, remove SNPs that had greater than two conflicts
        # in the first pass.
        if second_pass:
            if self.coordinates in snps_to_be_removed:
                return

        start_pos = self.position - 1
        end_pos = self.position

        orig_refList = []
        orig_altList = []
        # These will correspond to the *real* reference and alternate bases,
        # not the reversed ones in the new bam
        new_refList = []
        new_altList = []

        orig_sam = pysam.Samfile(orig_bam, 'rb')
        new_sam = pysam.Samfile(new_bam, 'rb')

        for pileupcolumn in orig_sam.pileup(reference=self.chrom, start=start_pos, end=end_pos):
            if pileupcolumn.pos == start_pos:
                bases = set()
                for count, read in enumerate(pileupcolumn.pileups, 1):
                    quality = ord(read.alignment.qqual[read.qpos]) - 33
                    bases.add(read.alignment.seq[read.qpos])
                    if read.alignment.seq[read.qpos] == self.real_reference_base and quality >= 20:
                        orig_refList.append(read.alignment.qname)
                    elif read.alignment.seq[read.qpos] == self.real_alternate_base and quality >= 20:
                        orig_altList.append(read.alignment.qname)
                # Check for three state SNPs
                if len(bases) > 2:
                    return
                # Check coverage
                cutoff = self.cutoffs.get(count, 24)
                if count >= 15 and len(orig_refList) >= cutoff and len(orig_altList) >= cutoff:
                    pass
                else:
                    return

        for pileupcolumn in new_sam.pileup(reference=self.chrom, start=start_pos, end=end_pos):
            if pileupcolumn.pos == start_pos:
                bases = set()
                for count, read in enumerate(pileupcolumn.pileups, 1):
                    quality = ord(read.alignment.qqual[read.qpos]) - 33
                    bases.add(read.alignment.seq[read.qpos])
                    if (read.alignment.seq[read.qpos] == self.real_reference_base and
                        quality >= 20):
                        new_refList.append(read.alignment.qname)
                    elif (read.alignment.seq[read.qpos] == self.real_alternate_base and
                          quality >= 20):
                        new_altList.append(read.alignment.qname)
                # Check for three state SNPs
                if len(bases) > 2:
                    return
                # Check coverage
                cutoff = self.cutoffs.get(count, 24)
                if count >= 15 and len(new_refList) >= cutoff and len(new_altList) >= cutoff:
                    pass
                else:
                    return

        orig_sam.close()
        new_sam.close()

        # If this is the second pass through the script, remove reads that had 1 or 2 conflicts
        if second_pass:
            orig_refList = [x for x in orig_refList if x not in orig_reads_to_be_removed]
            orig_altList = [x for x in orig_altList if x not in orig_reads_to_be_removed]
            new_refList = [x for x in new_refList if x not in new_reads_to_be_removed]
            new_altList = [x for x in new_altList if x not in new_reads_to_be_removed]

        fields = ('id', 'orig_refReads', 'orig_altReads', 'new_refReads', 'new_altReads')
        snp = collections.namedtuple('snp', fields)

        return snp(self.coordinates,
                   orig_refList,
                   orig_altList,
                   new_refList,
                   new_altList)


class Gene(object):

    def __init__(self, current_gene=None):
        self.current_gene = current_gene
        self.ids = []
        self.orig_refReads = {}
        self.orig_altReads = {}
        self.new_refReads = {}
        self.new_altReads = {}
        self.orig_snp_stats = {}
        self.new_snp_stats = {}
        self._ref_set = None
        self._alt_set = None

    def __repr__(self):
        return 'Gene("{}")'.format(self.current_gene)

    def add_snp(self, snp):

        self.orig_refReads[snp.id] = snp.orig_refReads
        self.orig_altReads[snp.id] = snp.orig_altReads
        self.new_refReads[snp.id] = snp.new_refReads
        self.new_altReads[snp.id] = snp.new_altReads

        self.ids.append(snp.id)

        # Create a new key in our snp stats dictionaries
        self.orig_snp_stats[snp.id] = dict(validated=0,
                                           unvalidated=(len(snp.orig_refReads) + len(snp.orig_altReads)),
                                           conflicts=0)
        self.new_snp_stats[snp.id] = dict(validated=0,
                                          unvalidated=(len(snp.new_refReads) + len(snp.new_altReads)),
                                          conflicts=0)

    @property
    def ref_set(self):
        if not self._ref_set:
            read_lists = (lst for lst in (chain(self.orig_refReads.itervalues(),
                                                self.new_refReads.itervalues())))
            reads = chain(*read_lists)
            self._ref_set = {read for read in reads}
        return self._ref_set

    @property
    def alt_set(self):
        if not self._alt_set:
            read_lists = (lst for lst in (chain(self.orig_altReads.itervalues(),
                                                self.new_altReads.itervalues())))
            reads = chain(*read_lists)
            self._alt_set = {read for read in reads}
        return self._alt_set

    def write_gene(self, writer, gene_info):
        # Only write if a snp was successfully added (i.e. not all SNPs returned None)
        if len(self.ids) > 0:
            newrow = {}
            newrow['chrom'] = gene_info.chrom
            newrow['R_Depth'] = len(self.ref_set)
            newrow['A_Depth'] = len(self.alt_set)
            newrow['gene_id'] = gene_info.gene_id
            newrow['exon_number'] = gene_info.exon_number
            newrow['gene_name'] = gene_info.gene_name
            writer.writerow(newrow)

    def write_snp_stats(self, writer):

        orig_ref_list = list(chain(*self.orig_refReads.values()))
        orig_alt_list = list(chain(*self.orig_altReads.values()))
        new_ref_list = list(chain(*self.new_refReads.values()))
        new_alt_list = list(chain(*self.new_altReads.values()))

        orig_refCount = collections.Counter(orig_ref_list)
        orig_altCount = collections.Counter(orig_alt_list)
        new_refCount = collections.Counter(new_ref_list)
        new_altCount = collections.Counter(new_alt_list)

        def validator(counter, reads_dict, stats_dict):
            for read in counter:
                # pick out reads that appear more than once (according to our counter)
                # and are thus validated
                if counter[read] >= 2:
                    # Add this read to our summary validated variable
                    # identify the SNPs each read belongs to.
                    for idd, reads in reads_dict.iteritems():
                        # Add 1 to validation for that SNP if the read is in it.
                        # Remember that our stats dictionaries already have the id
                        # of all the SNPs.
                        # At the same time, subtract 1 from the unvalidated count for that SNP,
                        # since one is now validated
                        if read in reads:
                            stats_dict[idd]['validated'] += 1
                            stats_dict[idd]['unvalidated'] -= 1

        validator(orig_refCount, self.orig_refReads, self.orig_snp_stats)
        validator(orig_altCount, self.orig_altReads, self.orig_snp_stats)
        validator(new_refCount, self.new_refReads, self.new_snp_stats)
        validator(new_altCount, self.new_altReads, self.new_snp_stats)

        orig_refSet = set(orig_ref_list)
        orig_altSet = set(orig_alt_list)
        orig_conflicts = (orig_refSet.intersection(orig_altSet))
        for read in orig_conflicts:
            # Store the number of times this conflicted read appears in a SNP.
            conflicts_for_read = 0
            for idd, reads in chain(self.orig_refReads.iteritems(), self.orig_altReads.iteritems()):
                if read in reads:
                    self.orig_snp_stats[idd]['conflicts'] += 1
                    conflicts_for_read += 1
            # if the read conflicted 1 or 2 times, store this read to be removed on the
            # second pass.
            if 0 < conflicts_for_read <= 2:
                orig_reads_to_be_removed.add(read)

        new_refSet = set(new_ref_list)
        new_altSet = set(new_alt_list)
        new_conflicts = (new_refSet.intersection(new_altSet))
        for read in new_conflicts:
            # Store the number of times this conflicted read appears in a SNP.
            conflicts_for_read = 0
            for idd, reads in chain(self.new_refReads.iteritems(), self.new_altReads.iteritems()):
                if read in reads:
                    self.new_snp_stats[idd]['conflicts'] += 1
                    conflicts_for_read += 1
            # if the read conflicted 1 or 2 times, store this read to be removed on the
            # second pass.
            if 0 < conflicts_for_read <= 2:
                new_reads_to_be_removed.add(read)

        for idd, snp in self.orig_snp_stats.iteritems():
            newrow = {}
            newrow['chrom'] = idd[0]
            newrow['pos'] = idd[1]
            newrow['gene'] = self.current_gene
            newrow['fastq'] = 'original'
            newrow['validated'] = snp['validated']
            newrow['unvalidated'] = snp['unvalidated']
            newrow['conflicts'] = snp['conflicts']
            writer.writerow(newrow)
            # If a SNP has more than 2 conflicts, store it to be removed on the second pass.
            if snp['conflicts'] > 2:
                snps_to_be_removed.add(idd)

        for idd, snp in self.new_snp_stats.iteritems():
            newrow = {}
            newrow['chrom'] = idd[0]
            newrow['pos'] = idd[1]
            newrow['gene'] = self.current_gene
            newrow['fastq'] = 'new'
            newrow['validated'] = snp['validated']
            newrow['unvalidated'] = snp['unvalidated']
            newrow['conflicts'] = snp['conflicts']
            writer.writerow(newrow)
            if snp['conflicts'] > 2:
                snps_to_be_removed.add(idd)


def combine_snps(input_orig, input_new, output_f, orig_bam, new_bam, ref_vcf, snp_stats_f):

    with open(ref_vcf, 'rb') as ref_vcf:
        ref_reader = csv.reader(ref_vcf, delimiter='\t')

        def get_single_base_positions(reader):
            for row in reader:
                if len(row[3]) == 1 and len(row[4]) == 1:
                    yield [row[0], int(row[1]), row[3], row[4]]  # chrom, pos, ref, alt
                else:
                    assert len(row[3]) == len(row[4])
                    position = int(row[1])
                    for refbase, altbase in zip(row[3], row[4]):
                        yield [row[0], position, refbase, altbase]
                        position += 1

        Ref_tup = collections.namedtuple('Ref_tup', ['reference', 'alternate'])
        # Dictionary containing the coordinates and the ref and alt bases from the reference vcf,
        # the known SNPs file for filtering.
        ref_dict = {(row[0], row[1]): Ref_tup(row[2], row[3]) for row in
                    get_single_base_positions(ref_reader)}

    with open(input_orig, 'rb') as f_orig, open(input_new, 'rb') as f_new:

        fields = ('CHROM', 'POS', 'REF', 'ALT', 'RD', 'AD', 'gene_id', 'exon_number', 'gene_name')
        reader_orig = csv.DictReader(f_orig, delimiter='\t', fieldnames=fields)
        reader_new = csv.DictReader(f_new, delimiter='\t', fieldnames=fields)

        orig_row_holder = {(row['CHROM'], row['POS']): row for row in reader_orig}
        new_row_holder = {(row['CHROM'], row['POS']): row for row in reader_new}

        variants = []
        for coord in orig_row_holder:
            if coord in new_row_holder:
                v = Variant(original=orig_row_holder[coord],
                            new=new_row_holder[coord],
                            refdict=ref_dict)
                variants.append(v)

        variants.sort(key=lambda variant: variant.gene_id)

        with open(output_f, 'wb') as fout, open(snp_stats_f, 'wb') as snp_out:

            write_fields = ('chrom', 'R_Depth', 'A_Depth', 'gene_id', 'exon_number', 'gene_name')
            gene_writer = csv.DictWriter(fout, delimiter='\t', fieldnames=write_fields,
                                         lineterminator='\n')
            header = {field: field for field in write_fields}
            gene_writer.writerow(header)

            snp_fields = ('chrom', 'pos', 'gene', 'fastq', 'validated', 'unvalidated', 'conflicts')
            snp_writer = csv.DictWriter(snp_out, delimiter='\t', fieldnames=snp_fields,
                                        lineterminator='\n')
            header = {field: field for field in snp_fields}
            snp_writer.writerow(header)

            genebuffer = None

            for rowcount, variant in enumerate(variants, 1):
                if rowcount % 10000 == 0:
                    print 'rows examined: {}'.format(rowcount)

                if variant.gene_id == genebuffer:
                    # examining a gene we have seen before
                    snp = variant.get_snp_data(orig_bam, new_bam)
                    if snp:
                        gene.add_snp(snp)

                else:
                    # examining a new gene
                    if genebuffer:
                        gene.write_gene(gene_writer, gene_info)
                        gene.write_snp_stats(snp_writer)

                    # Add the new gene's information
                    genebuffer = variant.gene_id
                    Gene_info = collections.namedtuple('Gene_info', ('chrom', 'gene_id', 'exon_number',
                                                                     'gene_name'))
                    gene_info = Gene_info(variant.chrom, variant.gene_id, variant.exon_number,
                                          variant.gene_name)
                    gene = Gene(current_gene=variant.gene_id)
                    snp = variant.get_snp_data(orig_bam, new_bam)
                    if snp:
                        gene.add_snp(snp)

            # Once we reach the end of the csv we have one row left to write,
            # as this will not be triggered by the 'else' statement above
            gene.write_gene(gene_writer, gene_info)
            gene.write_snp_stats(snp_writer)


def snp2gene(input_orig, input_new, output_f, orig_bam, new_bam, ref_vcf, snp_stats_f, cutoff_table):

    def get_coverage_cutoffs(cutoff_table):
        with open(cutoff_table, 'rb') as cf:
            cutoff_reader = csv.reader(cf, delimiter=',')
            next(cutoff_reader)  # remove header
            cutoffs = {int(row[0]): int(row[1]) for row in cutoff_reader}
        return cutoffs

    # Here we store our variable cutoff dictionary,
    # Our reads to be removed on the second pass through the script,
    # And our SNPs to be removed on the second pass through the script.
    # The removal of reads and SNPs is based on the number of conflicts.
    global cos
    cos = get_coverage_cutoffs(cutoff_table)
    global orig_reads_to_be_removed
    global new_reads_to_be_removed
    orig_reads_to_be_removed = set()
    new_reads_to_be_removed = set()
    global snps_to_be_removed
    snps_to_be_removed = set()

    global second_pass
    second_pass = False
    print('Combining snps to genes, first pass')
    combine_snps(input_orig, input_new, output_f, orig_bam, new_bam, ref_vcf, snp_stats_f)

    second_pass = True
    print('Combining snps to genes, second pass')
    combine_snps(input_orig, input_new, output_f, orig_bam, new_bam, ref_vcf,
                 (snp_stats_f + '2'))


# Go to directory '/scratch/Drew/bias_correc_work/testing_refandAlt_fasta/altGenome_b2N1'
# input_orig = 'original/16_A12_pUn_down/16_A12_pUn_down_INTER_py.csv'
# input_new = 'alternate/16_A12_pUn_down/16_A12_pUn_down_INTER_py.csv'
# output_f = 'testing_refactoring'
# original_bam_path = 'original/16_A12_pUn_down/16_A12_pUn_down_thout/filter.bam'
# new_bam_path = 'alternate/16_A12_pUn_down/16_A12_pUn_down_thout/filter.bam'
# coverage_cutoff = 15

# combine_snps(input_orig, input_new, output_f, original_bam_path,
#                   new_bam_path, 'testing_refactoring_snpstats')
