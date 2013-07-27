import csv
import pprint


def read_hapcut_output(input_hap):
    with open(input_hap) as f:
        record_lines = []
        for line in f:
            if not line.startswith('********'):
                record_lines.append(line)
            else:
                yield record_lines
                record_lines = []


class Phased_snp:
    def __init__(self, snp):
        self.snp = snp
        self.id = (snp[3], snp[4])

    def get_phase(self):
        if int(self.snp[1]) == 0 and int(self.snp[2]) == 1:
            return 0
        elif int(self.snp[2]) == 0 and int(self.snp[1]) == 1:
            return 1
        else:
            raise AssertionError('Error with phasing, row:', self.snp)


def phased_snps(input_hap):
    phased_dict = {}
    for record in read_hapcut_output(input_hap):
        snps = csv.reader(record[1:], delimiter='\t')
        for snp in snps:
            phased_snp = Phased_snp(snp)
            phased_dict[phased_snp.id] = phased_snp
    return phased_dict


def make_phased_vcfs(input_vcf, input_hap):
    assert input_vcf[-4:] == '.vcf'
    out_orig = input_vcf[:-4] + '__orig__.vcf'
    out_new = input_vcf[:-4] + '__new__.vcf'

    with open(input_vcf, newline='') as i_vcf:
        fields = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                  'INFO', 'FORMAT', 'SAMPLE')
        vcf_reader = csv.DictReader(i_vcf, delimiter='\t', fieldnames=fields)

        with open(out_orig, 'w', newline='') as o_orig, open(out_new, 'w', newline='') as o_new:

            orig_writer = csv.DictWriter(o_orig, fields, delimiter='\t', lineterminator='\n')
            new_writer = csv.DictWriter(o_new, fields, delimiter='\t', lineterminator='\n')

            phased_dict = phased_snps(input_hap)
            phase_count = 0

            for row in vcf_reader:
                # print(row)
                row_id = (row['CHROM'], row['POS'])
                # print(row_id)
                format_list = row['FORMAT'].split(':')
                sample_list = row['SAMPLE'].split(':')
                sample_dict = dict(zip(format_list, sample_list))
                sample_dict['RD'], sample_dict['AD'] = int(sample_dict['RD']), int(sample_dict['AD'])

                # If the SNP is homozygous alternate (i.e. there are only alternate reads),
                # Add the SNP to both vcfs
                if sample_dict['AD'] > 0 and sample_dict['RD'] == 0:
                    orig_writer.writerow(row)
                    new_writer.writerow(row)

                # If the SNP is heterozygous and not phased,
                # Add the SNP to the new vcf
                elif (sample_dict['AD'] > 0 and sample_dict['RD'] > 0 and
                      row_id not in phased_dict):
                    new_writer.writerow(row)

                # If the SNP is heterozygous and phased,
                # Add the SNP to orig_vcf if it is phased as 0, add to new_vcf if phased as 1
                elif (sample_dict['AD'] > 0 and sample_dict['RD'] > 0 and
                      row_id in phased_dict):
                    phase_count += 1
                    if phased_dict[row_id].get_phase() == 0:
                        orig_writer.writerow(row)
                    elif phased_dict[row_id].get_phase() == 1:
                        new_writer.writerow(row)
            print('Number of SNPs from vcf that were phased:', phase_count)



make_phased_vcfs('/path/to/16_A12_pUn_down_freeze2.vcf',
                 '/path/to/output_haplotype_file')
