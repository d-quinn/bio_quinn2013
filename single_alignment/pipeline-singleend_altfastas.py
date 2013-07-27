#!/usr/bin/python

## Global modules
import os
from os.path import join
import logging
import shutil
import pysam

## Local modules
import mglobals
import helpers
import snp2gene

log = logging.getLogger('pipeline')


@helpers.log_func
def my_tophat():

    if mglobals.original:
        log.info('Aligning reads to the original reference fasta')
        fastas = mglobals.original_fastas
    else:
        log.info('Aligning reads to alternate reference fastas')
        fastas = mglobals.alternate_fastas

    @helpers.multiprocess(zip(mglobals.samples_list, fastas))
    def tophat_call(sample, ref_fasta):

        if mglobals.original:
            os.chdir(join(mglobals.original_path, sample))
        else:
            os.chdir(join(mglobals.alternate_path, sample))

        ref_fasta_base = ref_fasta.split('.')[0]

        mismatches = '2'
        number_of_samples = len(mglobals.samples_list)
        threads_per_sample = mglobals.cpu_count//number_of_samples

        threads = str(threads_per_sample)
        log.info('threads per sample ' + threads)

        log.info('tophat: aligning sample {} with ref fasta {}'.format(sample, ref_fasta))
        tophat_params = ['nice', '-n', '5',
                         'tophat',
                         '-p', threads,
                         '-G', mglobals.dros_gtf,
                         '--transcriptome-index=../transcriptome_data/known',
                         '-N', mismatches,
                         '--b2-L', '20',
                         '--b2-N', '1',
                         '--read-edit-dist', mismatches,
                         '-o', (sample + '_thout'),
                         '--no-novel-juncs',
                         ref_fasta_base,
                         join(mglobals.samples_path, (sample + '.fastq'))]
        helpers.sub_call(tophat_params)

        log.info('tophat: finished analyzing sample: {} with ref fasta: {}'.format(sample, ref_fasta))

    # Copy transcriptome index to original_path and alternate_path
    for path in [mglobals.original_path, mglobals.alternate_path]:
        if not os.path.exists(join(path, 'transcriptome_data')):
            log.info('Linking transcriptome data to ' + path)
            os.symlink(mglobals.dros_gtf_index, join(path, 'transcriptome_data'))

    tophat_call()


@helpers.log_func
def my_alignment_filter():

    @helpers.multiprocess(mglobals.samples_list)
    def filter_call(sample):

        if mglobals.original:
            os.chdir(join(mglobals.original_path, sample, (sample + '_thout')))
        else:
            os.chdir(join(mglobals.alternate_path, sample, (sample + '_thout')))

        log.info('Filtering aligned reads for: ' + sample)

        # Index each bamfile
        if not os.path.exists('accepted_hits.bam.bai'):
            pysam.index('accepted_hits.bam')

        # Sort by NH flag
        raw_reads = pysam.Samfile('accepted_hits.bam', 'rb')
        filter_reads = pysam.Samfile('filter.bam', 'wb', template=raw_reads)
        for read in raw_reads.fetch():
            if ('NH', 1) in read.tags:
                filter_reads.write(read)
        raw_reads.close()
        filter_reads.close()

        pysam.index('filter.bam')

    filter_call()


@helpers.log_func
def my_pileup(out_file_extension='.mpileup'):

    if mglobals.original:
        fastas = mglobals.original_fastas
    else:
        fastas = mglobals.alternate_fastas

    @helpers.multiprocess(zip(mglobals.samples_list, fastas))
    def pileup_call(sample, ref_fasta):
        if mglobals.original:
            os.chdir(join(mglobals.original_path, sample))
        else:
            os.chdir(join(mglobals.alternate_path, sample))

        log.info('mpileup: creating .mpileup file for {} with ref fasta: {}'.format(sample, ref_fasta))
        pileup_command = ['nice', '-n', '5',
                          'samtools', 'mpileup',
                          '-B',
                          '-d10000000',
                          '-f', ref_fasta,
                          join((sample + '_thout'), 'filter.bam')]

        output_file = sample + out_file_extension
        with open(output_file, 'w') as output_file:
            helpers.sub_call(pileup_command, stdout=output_file)
        log.info('mpileup: finished for {} with ref fasta: {}'.format(sample, ref_fasta))

    pileup_call()


@helpers.log_func
def my_variant_calls(in_file_extension='.mpileup', out_file_extension='.vcf'):
    '''
    Note, build alternate fastas depends on the out_file_extension being '.vcf'.
    '''
    @helpers.multiprocess(mglobals.samples_list)
    def variant_calls_call(sample):
        if mglobals.original:
            os.chdir(join(mglobals.original_path, sample))
        else:
            os.chdir(join(mglobals.alternate_path, sample))

        log.info('Varscan: creating csv for: ' + sample)

        varscan_command = ['nice', '-n', '5',
                           'java', '-jar', mglobals.varscan_path,
                           'mpileup2snp',
                           (sample + in_file_extension),
                           '--min-coverage', '2',
                           '--min-avg-qual', '20',
                           '--strand-filter', '0',
                           '--p-value', '1',
                           '--min-var-freq', '1e-10',
                           '--output-vcf', '1',
                           ]

        output_file = sample + out_file_extension
        with open(output_file, 'w') as out:
            helpers.sub_call(varscan_command, stdout=out)
        log.info('varscan finished for: ' + sample)

    variant_calls_call()


@helpers.log_func
def cov_and_dgrp_filter(in_file_extension='.vcf', out_file_extension='_freeze2.vcf'):

    @helpers.multiprocess(mglobals.samples_list)
    def filter_call(sample):
        if mglobals.original:
            os.chdir(join(mglobals.original_path, sample))
        else:
            os.chdir(join(mglobals.alternate_path, sample))

        log.info('Filtering {0} by coverage'.format(sample))
        helpers.filter_vcf_by_coverage_cutoffs(vcf=(sample + in_file_extension),
                                               cutoff_table=mglobals.coverage_cutoffs)

        log.info('Filtering {0} according to SNP file: {1}'.format(sample, mglobals.current_snp_file))
        dgrp_intersect_command = ['nice', '-n', '5',
                                  'intersectBed',
                                  '-a', (sample + '_covfil.vcf'),  # the output of the helper
                                                                   # function above.
                                  '-b', mglobals.current_snp_file,
                                  '-wa'
                                  ]
        sample_dgrp_intersect = sample + out_file_extension
        with open(sample_dgrp_intersect, 'w') as out:
            helpers.sub_call(dgrp_intersect_command, stdout=out)

    filter_call()


@helpers.log_func
def annotate_vcf(in_file_extension='_freeze2.vcf', out_file_extension='_geneannot.vcf'):

    @helpers.multiprocess(mglobals.samples_list)
    def annotate_call(sample):
        if mglobals.original:
            os.chdir(join(mglobals.original_path, sample))
        else:
            os.chdir(join(mglobals.alternate_path, sample))

        log.info('Annotating ' + sample + in_file_extension + ' with ' + mglobals.current_genes_file)
        gtf_intersect_command = ['nice', '-n', '5',
                                 'intersectBed',
                                 '-a', (sample + in_file_extension),
                                 '-b', mglobals.current_genes_file,
                                 '-wa',
                                 '-wb'
                                 ]
        sample_gtf_intersect = sample + out_file_extension
        with open(sample_gtf_intersect, 'w') as out:
            helpers.sub_call(gtf_intersect_command, stdout=out)

    annotate_call()


@helpers.log_func
def build_alternate_fastas(in_file_extension='_geneannot.vcf'):
    # If we are doing the original alignment, we can now build the alternate
    # reference fastas for each sample

    @helpers.multiprocess(mglobals.samples_list)
    def build_fastas_call(sample):
        os.chdir(join(mglobals.original_path, sample))
        log.info('Beginning to build alternate fasta for: ' + sample)
        fixed_vcf = sample + '_fix.vcf'
        log.info('Removing duplicated annotations (per transcript annotations)')
        helpers.remove_dups(input_f=(sample + in_file_extension),
                            output_f=(sample + '.temp'))
        log.info('Removing duplicate alleles and adding header')
        # The fact that the original vcf was named sample.vcf is hardcoded
        # here. Be careful.
        helpers.vcf_fix(template_f=(sample + '.vcf'),
                        input_f=(sample + '.temp'),
                        output_f=fixed_vcf)
        # Delete temporary file
        os.remove(sample + '.temp')
        log.info('Creating alternate fasta')
        new_fasta = sample + '_unfixed.fa'
        helpers.sub_call(['nice', '-n', '5',
                          'java', '-Xmx2g', '-jar',
                          mglobals.gatk_path,
                          '-R', 'genome.fa',
                          '-T', 'FastaAlternateReferenceMaker',
                          '-o', new_fasta,
                          '--variant', fixed_vcf])
        # Fix the fasta
        log.info('Fixing gatk fasta')
        # If you change this name, you need to change the alternate fastas list as well.
        final_fasta = sample + '.fa'
        helpers.fasta_fix(input_f=new_fasta, output_f=final_fasta)
        # Delete the unfixed version
        os.remove(new_fasta)
        log.info('Moving new fasta to: ' + join(mglobals.alternate_path, sample))
        shutil.move(final_fasta, join(mglobals.alternate_path, sample))
        log.info('Indexing new fasta')
        os.chdir(join(mglobals.alternate_path, sample))
        helpers.sub_call(['bowtie2-build',
                          '-f', final_fasta,
                          sample])

    build_fastas_call()


@helpers.log_func
def vcf_to_csv(in_file_extension='_geneannot.vcf',
               out_file_extension='_INTER_py.csv'):

    @helpers.multiprocess(mglobals.samples_list)
    def vcf_to_csv_call(sample):
        if mglobals.original:
            os.chdir(join(mglobals.original_path, sample))
        else:
            os.chdir(join(mglobals.alternate_path, sample))

        log.info('Converting vcf to csv for: ' + sample)
        snp2gene.converter(input_f=(sample + in_file_extension),
                           output_f=(sample + out_file_extension))

    vcf_to_csv_call()


@helpers.log_func
def combine_snps(in_file_extension='_INTER_py.csv'):

    os.chdir(mglobals.samples_path)

    @helpers.multiprocess(mglobals.samples_list)
    def combine_snps_call(sample):
        log.info('Combining SNPs for: ' + sample)

        snps_combine.combine_SNPs(orig_f=join(mglobals.original_path, sample,
                                             (sample + in_file_extension)),
                                  new_f=join(mglobals.alternate_path, sample,
                                            (sample + in_file_extension)),
                                  orig_bam=join(mglobals.original_path, sample,
                                               (sample + '_thout'), 'filter.bam'),
                                  new_bam=join(mglobals.alternate_path, sample,
                                              (sample + '_thout'), 'filter.bam'),
                                  ref_vcf=join(mglobals.current_snp_file),
                                  output_f=(sample + '_snps.csv'),
                                  cutoff_table=mglobals.coverage_cutoffs)

        mean_propR, snp_count = snps_combine.quick_mean_propR(sample + '_snps.csv')
        log.info('Mean proportion reference for {} = {}'.format(sample, mean_propR))
        log.info('\tNumber of snps = {}'.format(snp_count))

    combine_snps_call()


@helpers.log_func
def csv_recalibrate(in_file_extension='_INTER_py.csv', out_file_extension='_genes.csv'):

    os.chdir(mglobals.samples_path)

    @helpers.multiprocess(mglobals.samples_list)
    def csv_recalibrate_call(sample):
        log.info('Combining genes for: ' + sample)

        snp2gene.snp2gene(input_orig=join(mglobals.original_path, sample,
                                          (sample + in_file_extension)),
                          input_new=join(mglobals.alternate_path, sample,
                                         (sample + in_file_extension)),
                          output_f=(sample + out_file_extension),
                          orig_bam=join(mglobals.original_path, sample,
                                        (sample + '_thout'), 'filter.bam'),
                          new_bam=join(mglobals.alternate_path, sample,
                                       (sample + '_thout'), 'filter.bam'),
                          ref_vcf=mglobals.current_snp_file,
                          snp_stats_f=(sample + '_snp_stats.csv'),
                          cutoff_table=mglobals.coverage_cutoffs)

    csv_recalibrate_call()


def main():

    mglobals.original = True
    my_tophat()
    my_alignment_filter()
    my_pileup()
    my_variant_calls()
    cov_and_dgrp_filter()
    annotate_vcf()
    # build_alternate_fastas()
    vcf_to_csv()

    mglobals.original = False
    # my_tophat()
    # my_alignment_filter()
    # my_pileup()
    # my_variant_calls()
    # cov_and_dgrp_filter()
    # annotate_vcf()
    # vcf_to_csv()

    # combine_snps()
    # csv_recalibrate()

    log.info('Pipeline completed successfully!')

if __name__ == '__main__':
    main()
