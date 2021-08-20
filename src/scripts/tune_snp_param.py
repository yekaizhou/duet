#!/usr/bin/env python3
# coding=utf-8

import os
import logging
import argparse
import time
import sys
from duet.utils import set_logging, check_envs
from duet.snp_phasing import snp_phasing
from duet.sv_calling import sv_calling
from duet.sv_phasing import sv_phasing

def parse_args(argv):
    parser = argparse.ArgumentParser(
    #    description = 'SNP-Assisted Phased SV Detection from Low-depth Long-reads'
    )
    parser.add_argument('--skip_full_aln', action = 'store_true',
        help = 'whether skipping full alignment [%(default)s]')
    parser.add_argument('-m', '--min_allele_frequency', type = float, default = 0.25,
        help = 'minimum allele frequency required to call a candidate SNP [%(default)s]')
    parser.add_argument('--skip_indel', action = 'store_true',
        help = 'whether skipping indel calling [%(default)s]')
    parser.add_argument('-t', '--thread', type = int, default = 40,
        help = 'number of threads to use [%(default)s]')
    parser.add_argument('-c', '--cluster_max_distance', type = float, default = 0.9,
        help = 'maximum span-position distance between SV marks in a cluster to call a SV candidates [%(default)s]')
    parser.add_argument('-s', '--sv_min_size', type = int, default = 50,
        help = 'minimum SV size to be reported [%(default)s]')
    parser.add_argument('-r', '--min_support_read', type = int, default = 2,
        help = 'minimum number of reads that support a SV to be reported [%(default)s]')
    parser.add_argument('BAM', type = str,
            help = 'sorted alignment file in .bam format (along with .bai file in the same directory)')
    parser.add_argument('REFERENCE', type = str,
            help = 'indexed reference genome in .fasta format (along with .fai file in the same directory)')
    parser.add_argument('OUTPUT', type = str,
        help = 'working and output directory (existing files in the directory will be overwritten)')
    args = parser.parse_args()
    return args

def snp_calling(home, ref_path, aln_path, maf, thread, pileup, noindel):
    lines = '*************************'
    logging.info(lines + ' SNP CALLING STARTED ' + lines)
    starttime = time.time()
    snp_calling_home = home + '/snp_calling/'
    os.system('mkdir ' + snp_calling_home)
    command = 'run_clair3.sh -b ' + aln_path + ' -f ' + ref_path + ' -m "${CONDA_PREFIX}/bin/models/ont" -t ' + \
              str(thread) + ' -p ont -o ' + snp_calling_home + ' --snp_min_af=' + str(maf)
    if pileup:
        command += ' --pileup_only'
    if noindel:
        command += ' --call_snp_only'
    os.system(command)
    logging.info(lines + ' SNP CALLING COMPLETED IN ' + str(round(time.time() - starttime, 3)) + 's ' + lines)

def main(argv):
    args = parse_args(argv)
    check_envs(args.REFERENCE, args.BAM)
    if not os.path.exists(args.OUTPUT):
        os.makedirs(args.OUTPUT)
    set_logging(args.OUTPUT)
    lines = '*************************'
    starttime = time.time()
    logging.info(lines + ' DUET STARTED ' + lines)
    snp_calling(args.OUTPUT, args.REFERENCE, args.BAM, args.min_allele_frequency, args.thread, args.skip_full_aln, args.skip_indel)
    sv_calling(args.OUTPUT, args.REFERENCE, args.BAM, args.cluster_max_distance, args.sv_min_size)
    snp_phasing(args.OUTPUT, args.REFERENCE, args.BAM, args.thread)
    sv_phasing(args.OUTPUT, args.sv_min_size, args.min_support_read, args.thread)
    logging.info(lines + ' DUET FINISHED IN ' + str(round(time.time() - starttime, 3)) + 's ' + lines)
    logging.info('OUTPUT .VCF FILE AT ' + args.OUTPUT + '/phased_sv.vcf')

if __name__ == '__main__':
    main(sys.argv[1:])