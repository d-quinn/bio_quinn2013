import pysam
import csv


def get_new_row(orig_bam, new_bam, row):

    assert '1/1' in row['DGRP-362'] or '1/1' in row['DGRP-765']

    # Check that we actually have a SNP here - i.e. the fastas will have
    # different bases at this position
    if '1/1' in row['DGRP-362'] and '1/1' in row['DGRP-765']:
        return

    ref_base = row['REF']
    alt_base = row['ALT']

    chrom = row['CHROM']
    startpos = int(row['POS']) - 1
    endpos = startpos + 1

    orig_ref_reads = []
    orig_alt_reads = []
    new_ref_reads = []
    new_alt_reads = []

    orig_sam = pysam.Samfile(orig_bam, 'rb')
    new_sam = pysam.Samfile(new_bam, 'rb')

    for pileupcolumn in orig_sam.pileup(reference=chrom, start=startpos, end=endpos):
        if pileupcolumn.pos == startpos:
            bases = set()
            for count, read in enumerate(pileupcolumn.pileups, 1):
                quality = ord(read.alignment.qqual[read.qpos]) - 33
                bases.add(read.alignment.seq[read.qpos])
                if read.alignment.seq[read.qpos] == ref_base and quality >= 20:
                    orig_ref_reads.append(read.alignment.qname)
                elif read.alignment.seq[read.qpos] == alt_base and quality >= 20:
                    orig_alt_reads.append(read.alignment.qname)
            # Check for three state SNPs
            if len(bases) > 2:
                return
            # Check coverage
            cutoff = cos.get(count, 24)
            if count >= 15 and len(orig_ref_reads) >= cutoff and len(orig_alt_reads) >= cutoff:
                pass
            else:
                return

    for pileupcolumn in new_sam.pileup(reference=chrom, start=startpos, end=endpos):
        if pileupcolumn.pos == startpos:
            bases = set()
            for count, read in enumerate(pileupcolumn.pileups, 1):
                quality = ord(read.alignment.qqual[read.qpos]) - 33
                bases.add(read.alignment.seq[read.qpos])
                if read.alignment.seq[read.qpos] == ref_base and quality >= 20:
                    new_ref_reads.append(read.alignment.qname)
                elif read.alignment.seq[read.qpos] == alt_base and quality >= 20:
                    new_alt_reads.append(read.alignment.qname)
            # Check for three state SNPs
            if len(bases) > 2:
                return
            # Check coverage
            cutoff = cos.get(count, 24)
            if count >= 15 and len(new_ref_reads) >= cutoff and len(new_alt_reads) >= cutoff:
                pass
            else:
                return

    orig_sam.close()
    new_sam.close()

    reference_depth = len(set(orig_ref_reads + new_ref_reads))
    alternate_depth = len(set(orig_alt_reads + new_alt_reads))

    newrow = {field: row[field] for field in ('CHROM', 'POS', 'REF', 'ALT')}
    newrow['R_Depth'] = reference_depth
    newrow['A_Depth'] = alternate_depth
    return newrow


def combine_SNPs(input_vcf, orig_bam, new_bam, output_f, cutoff_table):

    def get_coverage_cutoffs(cutoff_table):
        with open(cutoff_table, 'rb') as cf:
            cutoff_reader = csv.reader(cf, delimiter=',')
            next(cutoff_reader)  # remove header
            cutoffs = {int(row[0]): int(row[1]) for row in cutoff_reader}
        return cutoffs

    global cos
    cos = get_coverage_cutoffs(cutoff_table)

    with open(input_vcf, 'rb') as vcf:
        fields = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                  'INFO', 'FORMAT', 'DGRP-362', 'DGRP-765')
        vcfreader = csv.DictReader(vcf, fields, delimiter='\t')

        with open(output_f, 'wb') as fout:
            fields = ('CHROM', 'POS', 'REF', 'ALT', 'R_Depth', 'A_Depth')
            writer = csv.DictWriter(fout, fields, delimiter='\t', lineterminator='\n')
            writer.writerow({field: field for field in fields})
            for count, row in enumerate(vcfreader, 1):
                if count % 10000 == 0:
                    print 'rows examined:', count
                newrow = get_new_row(orig_bam, new_bam, row)
                if newrow:
                    writer.writerow(newrow)


## For generating the parental alignment counts


# input_vcf = '/scratch/Drew/DGRP_Work_FREEZE2/2_genomes_vcfs/freeze2_362_765_noINDEL_1-1.vcf'
# orig_bam = '/scratch/Drew/2_genome_fastas_retry/alternate/16_A12_pUn_down/16_A12_pUn_down_thout/filter.bam'
# new_bam = '/scratch/Drew/2_genome_fastas_retry/original/16_A12_pUn_down/16_A12_pUn_down_thout/filter.bam'

# combine_SNPs(input_vcf, orig_bam, new_bam, 'snps.vcf')

