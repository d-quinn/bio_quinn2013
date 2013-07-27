import csv
from sys import exit
from Bio import SeqIO
import subprocess
import multiprocessing
import functools
import logging
import pprint
import mglobals

log = logging.getLogger('pipeline')


def sub_call(command, stdout=None, shell=False):
    '''
    Subprocess manager: will log commands and raise an error
    if return code does not equal 0.
    '''
    pretty_command = pprint.pformat(command, indent=4)
    log.debug('running command:\n{}'.format(pretty_command))
    if stdout:
        log.debug('writing to: ' + stdout.name)
    subprocess.check_call(command, stdout=stdout, shell=shell)


def log_func(func):
    '''
    Decorator which will wrap functions so that the beginning
    and end is logged.
    '''
    @functools.wraps(func)
    def log_wrapper(*args, **kwargs):
        if mglobals.original:
            log.info('{:-^50}'.format(' Original: beginning ' + func.__name__ + ' '))
        else:
            log.info('{:-^50}'.format(' Alternate: beginning ' + func.__name__ + ' '))
        output = func(*args, **kwargs)
        if mglobals.original:
            log.info(' Original: finished ' + func.__name__ + ' ')
        else:
            log.info(' Alternate: finished ' + func.__name__ + ' ')
        return output
    return log_wrapper


def multiprocess(iterator):
    '''
    Returns the multiprocess decorator. This function's only purpose
    is to pass its iterator argument along.
    '''

    def multiprocess_decorator(func):
        '''
        Returns the multiprocess wrapper
        '''

        @functools.wraps(func)
        def multiprocess_wrapper():
            '''
            Multiprocesses a given function with arguments specified
            by the iterator.  This will create a separate process for each
            item in the iterator.
            '''
            log.info('#' * 45 + ' -> start')
            log.info('Calling multi_processor with target: ' + func.__name__)
            jobs = []
            for it in iterator:
                if isinstance(it, tuple):
                    j = multiprocessing.Process(target=func,
                                                name=(func.__name__ + ' ' + str(it)),
                                                args=it)
                else:
                    j = multiprocessing.Process(target=func,
                                                name=(func.__name__ + ' ' + it),
                                                args=(it,))
                jobs.append(j)
                j.start()

            for j in jobs:
                j.join()
            for j in jobs:
                if j.exitcode != 0:
                    raise AssertionError('multi process returned with non-0 '
                                         'exit code')
            log.info('#' * 45 + ' end <-')
        return multiprocess_wrapper

    return multiprocess_decorator


def remove_dups(input_f, output_f):
    with open(input_f, 'rb') as f, open(output_f, 'wb') as fout:
        reader = csv.reader(f, delimiter='\t')

        # sort by chromosome and then position
        possortedreader = sorted(reader, key=lambda row: int(row[1]))
        sortedreader = sorted(possortedreader, key=lambda row: row[0])

        writer = csv.writer(fout, delimiter='\t', lineterminator='\n')
        posBuffer = None
        for row in sortedreader:
            if row[:2] == posBuffer:
                continue
            posBuffer = row[:2]
            writer.writerow(row)


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
            if row[3] != row[4]:
                newrow = row
                alt_alleles = row[4]
                newalt_alleles = ','.join(list(set(row[4].split(','))))
                newrow[4] = newalt_alleles
                writer.writerow(newrow)


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


def filter_vcf_by_coverage_cutoffs(vcf, cutoff_table):
    with open(cutoff_table, 'rb') as cf:
        cutoff_reader = csv.reader(cf, delimiter=',')
        next(cutoff_reader)  # remove header
        cutoffs = {int(row[0]): int(row[1]) for row in cutoff_reader}

    assert vcf[-4:] == '.vcf'
    with open(vcf, 'rb') as vf, open((vcf[:-4] + '_covfil.vcf'), 'wb') as fout:
        # Remove commented lines
        for line in vf:
            if line.startswith('##'):
                continue
            else:
                break
        fields = ('CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER',
                  'INFO', 'FORMAT', 'Sample1')
        vcf_reader = csv.DictReader(vf, delimiter='\t', fieldnames=fields)
        writer = csv.DictWriter(fout, delimiter='\t', fieldnames=fields, lineterminator='\n')
        for row in vcf_reader:
            format_list = row['FORMAT'].split(':')
            sample_list = row['Sample1'].split(':')
            sample_dict = dict(zip(format_list, sample_list))
            num_ref = int(sample_dict['RD'])
            num_alt = int(sample_dict['AD'])
            coverage = num_ref + num_alt
            cutoff = cutoffs.get(coverage, 24)  # Will return 24 if key not in dict
            if coverage >= 10 and num_ref >= cutoff and num_alt >= cutoff:
                writer.writerow(row)

