import os
import glob
import subprocess
import csv
from Bio import SeqIO
import sys
from os.path import join as pjoin

supporting_files = '/path/to/Quinn2013_supporting_files'

def vcf_fix(template_f, input_f, output_f):
    '''
    Removes duplicate alleles from a varScan generated vcf. Also transfers the header
    from a specified template file.
    '''
    with open(template_f, 'rb') as t, open(input_f, 'rb') as f, open(output_f, 'wb') as fout:
        for line in t:
            if line.startswith('##'):
                fout.write(line)
            elif line.startswith('#') and not line.startswith('##'):
                fout.write(line)
                break
        reader = csv.reader(f, delimiter='\t')
        # Sort the vcfs here as well
        sortedreader = sorted(reader, key=lambda row: int(row[1]))
        sortedreader.sort(key=lambda row: row[0])
        writer = csv.writer(fout, delimiter='\t', lineterminator='\n')
        for row in sortedreader:
            writer.writerow(row)

            # if row[3] != row[4]:
            #     newrow = row
            #     alt_alleles = row[4]
            #     newalt_alleles = ','.join(list(set(row[4].split(','))))
            #     newrow[4] = newalt_alleles
            #     writer.writerow(newrow)


def fasta_fix(input_f, output_f):

    def replaceheaders():
        chromList = ['2L', '2R', '3L', '3R', '4', 'M', 'X']
        i = 0
        for rec in SeqIO.parse(input_f, "fasta"):
            rec.id = chromList[i]
            rec.description = ''
            yield rec
            i += 1

    newrecs = (rec for rec in replaceheaders())
    SeqIO.write(newrecs, output_f, "fasta")

# Here's the original freeze2 vcf
freeze2_vcf = pjoin(supporting_files, 'freeze2.vcf')

vcfs = ['/path/to/16_A12_pUn_down_freeze2__orig__.vcf',
        '/path/to/16_A12_pUn_down_freeze2__new__.vcf']

vcf_noextension = [vcf.split('.')[0] for vcf in vcfs]
zipped_vcfs = zip(vcfs, vcf_noextension)

print vcfs
print vcf_noextension
print zipped_vcfs


# Fix the vcfs and add header back
for vcf, vcf_noe in zipped_vcfs:
    vcf_fix(template_f=freeze2_vcf, input_f=vcf, output_f=vcf_noe + '_fixed.vcf')


gatk_path = pjoin(supporting_files, 'GenomeAnalysisTK-2.4-9/GenomeAnalysisTK.jar')

for vcf, vcf_noe in zipped_vcfs:
    fixed_vcf = vcf_noe + '_fixed.vcf'
    new_fasta = vcf_noe + '._unfixed_fa'
    print 'Calling alternatefastamaker with vcf', vcf, 'and vcf_noe', vcf_noe
    subprocess.check_call(['java', '-Xmx2g', '-jar',
                           gatk_path,
                           '-R', 'genome.fa',
                           '-T', 'FastaAlternateReferenceMaker',
                           '-o', new_fasta,
                           '--variant', fixed_vcf])

    final_fasta = vcf_noe + '.fa'
    fasta_fix(input_f=new_fasta, output_f=final_fasta)
    # Delete the unfixed version
    os.remove(new_fasta)
    subprocess.check_call(['bowtie2-build',
                           '-f', final_fasta,
                           vcf_noe])

