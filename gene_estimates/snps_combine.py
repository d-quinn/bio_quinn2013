import pysam
import csv
import sys
import collections


class Variant(object):

    def __init__(self, original, new, refdict):
        self.original = original
        self.new = new
        self.chrom = self.original['CHROM']
        self.position = int(self.original['POS'])
        self.coordinates = (self.chrom, self.position)
        self.real_reference_base = refdict[self.coordinates].reference
        self.real_alternate_base = refdict[self.coordinates].alternate
        self.cutoffs = cos  # cos is a global variable defined in the combine_snps function

    def __repr__(self):
        return 'Variant(original={}, new={})'.format(self.original, self.new)

    def get_coverage_cutoffs(csv):
        with open(csv, 'rb') as cf:
            cutoff_reader = csv.reader(cf, delimiter=',')
            next(cutoff_reader)  # remove header
            cutoffs = {int(row[0]): int(row[1]) for row in cutoff_reader}
        return cutoffs

    def get_new_row(self, orig_bam, new_bam):

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
                    if (read.alignment.seq[read.qpos] == self.real_reference_base and
                            quality >= 20):
                        orig_refList.append(read.alignment.qname)
                    elif (read.alignment.seq[read.qpos] == self.real_alternate_base and
                          quality >= 20):
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

        reference_depth = len(set(orig_refList + new_refList))
        alternate_depth = len(set(orig_altList + new_altList))

        newrow = {'CHROM': self.chrom,
                  'POS': self.position,
                  'REF': self.real_reference_base,
                  'ALT': self.real_alternate_base,
                  'R_Depth': reference_depth,
                  'A_Depth': alternate_depth}

        # print 'orig_refList:', orig_refList
        # print 'orig_altList:', orig_altList
        # print 'new_refList:', new_refList
        # print 'new_altList:', new_altList
        return newrow


def combine_SNPs(orig_f, new_f, orig_bam, new_bam, ref_vcf, output_f, cutoff_table):

    def get_coverage_cutoffs(cutoff_table):
        with open(cutoff_table, 'rb') as cf:
            cutoff_reader = csv.reader(cf, delimiter=',')
            next(cutoff_reader)  # remove header
            cutoffs = {int(row[0]): int(row[1]) for row in cutoff_reader}
        return cutoffs

    global cos
    cos = get_coverage_cutoffs(cutoff_table)

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


    with open(orig_f, 'rb') as of, open(new_f, 'rb') as nf:
        fields = ('CHROM', 'POS', 'REF', 'ALT', 'RD', 'AD',
                  'gene_id', 'exon_number', 'gene_name')
        oreader = csv.DictReader(of, fields, delimiter='\t')
        nreader = csv.DictReader(nf, fields, delimiter='\t')

        orig_row_holder = {(row['CHROM'], row['POS']): row for row in oreader}
        new_row_holder = {(row['CHROM'], row['POS']): row for row in nreader}

        variants = []
        for coord in orig_row_holder:
            if coord in new_row_holder:
                v = Variant(original=orig_row_holder[coord],
                            new=new_row_holder[coord],
                            refdict=ref_dict)
                variants.append(v)

    with open(output_f, 'wb') as fout:
        fields = ('CHROM', 'POS', 'REF', 'ALT', 'R_Depth', 'A_Depth')
        writer = csv.DictWriter(fout, fields, delimiter='\t', lineterminator='\n')
        writer.writerow({field: field for field in fields})
        for count, var in enumerate(variants, 1):
            if count % 10000 == 0:
                print 'rows examined:', count
            newrow = var.get_new_row(orig_bam, new_bam)
            if newrow:
                writer.writerow(newrow)

# orig_f = '/scratch/Drew/testdir/original/16_A12_pUn_down/16_A12_pUn_down_INTER_py.csv'
# new_f = '/scratch/Drew/testdir/alternate/16_A12_pUn_down/16_A12_pUn_down_INTER_py.csv'
# output_f = 'snpstest.csv'
# orig_bam = '/scratch/Drew/testdir/original/16_A12_pUn_down/16_A12_pUn_down_thout/filter.bam'
# alt_bam = '/scratch/Drew/testdir/alternate/16_A12_pUn_down/16_A12_pUn_down_thout/filter.bam'

# combine_SNPs(orig_f, new_f, orig_bam, alt_bam, output_f)


def quick_mean_propR(input_f):
    with open(input_f, 'rb') as f:
        reader = csv.DictReader(f, delimiter='\t')

        propRs = []
        for i, row in enumerate(reader, 1):
            ref_depth = float(row['R_Depth'])
            alt_depth = float(row['A_Depth'])
            if ref_depth != 0 and alt_depth != 0:
                propR = ref_depth/(ref_depth + alt_depth)
                propRs.append(propR)
        mean_propR = sum(propRs)/len(propRs)
        return (mean_propR, i)

# quick_mean_propR('snps.vcf')

    # i = 0
    # printcount = 0
    # while True and printcount < 2:
    #     row = variants[i].get_new_row(orig_bam, new_bam)
    #     if row:
    #         print variants[i]
    #         print row
    #         print i
    #         printcount += 1
    #     i += 1

    # sys.exit(0)
